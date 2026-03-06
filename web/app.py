#!/usr/bin/env python3
"""Flask web app: PI co-author network, affiliation network, and PI (author) network by research area."""
from __future__ import annotations

import json
import sys
import threading
from pathlib import Path
from queue import Empty, Queue

# Allow importing from project root and from web/
ROOT = Path(__file__).resolve().parent.parent
WEB_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(WEB_DIR))

from flask import Flask, render_template, request, jsonify, Response

from pi_coauthor_network import (
    parse_scholar_url_or_id,
    get_profile_by_scholar_id,
    get_profile_by_scholar_id_pubmed_fast,
    build_coauthor_network_from_papers,
    search_pi_pmids,
    fetch_articles,
    build_author_regions,
    _region_color,
)
from researchgate_fetcher import get_researchgate_papers, normalize_researchgate_url
from journal_if import filter_records_by_min_if

from fetch_pubmed import (
    search_pmids_by_query,
    fetch_articles as fetch_pubmed_articles,
    build_institution_network,
    build_corresponding_author_network,
    build_author_regions as build_author_regions_pubmed,
    institution_region,
    filter_institution_network_to_qs_top500,
)

app = Flask(__name__, static_folder="static", template_folder="templates")
MAX_AUTHORS = 50
MAX_NODES_DISPLAY = 1000  # When network has more nodes, keep top N by degree


def _user_facing_error(exc: Exception) -> str:
    """Turn a fetch exception into a short, actionable message."""
    msg = str(exc).strip()
    if not msg:
        return "Google Scholar could not be reached. Try again in a few minutes."
    if "cannot fetch" in msg.lower() or "max retries" in msg.lower() or "blocked" in msg.lower():
        return (
            "Google Scholar is blocking or rate-limiting requests. "
            "Try again later, use a different network (e.g. VPN), or run the CLI script locally: "
            "python pi_coauthor_network.py \"<profile_url>\""
        )
    return msg


def graph_to_vis(co_edges, authors, author_region):
    """Convert (co_edges, authors, author_region) to vis-network nodes and edges."""
    nodes = []
    for a in authors:
        region = author_region.get(a, "Others")
        nodes.append({
            "id": a,
            "label": a if len(a) <= 40 else a[:37] + "...",
            "title": a,
            "color": _region_color(region),
        })
    edges = []
    for (a, b), w in co_edges.items():
        edges.append({"from": a, "to": b, "value": w})
    return {"nodes": nodes, "edges": edges}


def _top_k_nodes_by_degree(edges: dict, nodes: set, k: int):
    """Keep only the top k nodes by degree (number of incident edges); return (filtered_edges, kept_nodes)."""
    if len(nodes) <= k:
        return edges, nodes
    degree = {}
    for n in nodes:
        degree[n] = 0
    for (a, b) in edges:
        degree[a] = degree.get(a, 0) + 1
        degree[b] = degree.get(b, 0) + 1
    sorted_nodes = sorted(degree.keys(), key=lambda x: -degree[x])
    keep = set(sorted_nodes[:k])
    filtered_edges = {(a, b): w for (a, b), w in edges.items() if a in keep and b in keep}
    return filtered_edges, keep


def _compute_network_analysis(edges: dict, nodes: set, top_n: int = 20):
    """
    Compute network analysis: top nodes by degree, communities, and basic stats.
    edges: dict (a, b) -> weight; nodes: set of node ids.
    Returns dict with top_degree, communities, n_components, density, avg_clustering.
    """
    import networkx as nx
    from networkx.algorithms.community import label_propagation_communities

    G = nx.Graph()
    G.add_nodes_from(nodes)
    for (a, b), w in edges.items():
        G.add_edge(a, b, weight=w)

    # Top nodes by degree (number of neighbors)
    degree = dict(G.degree())
    sorted_by_degree = sorted(degree.items(), key=lambda x: -x[1])
    top_degree = [{"node": n, "degree": d} for n, d in sorted_by_degree[:top_n]]

    # Communities (label propagation)
    try:
        communities_iter = label_propagation_communities(G)
        communities = []
        for i, comm in enumerate(communities_iter):
            comm_list = list(comm)
            if len(comm_list) >= 2:  # skip singletons for clarity, or include all
                communities.append({"id": i, "size": len(comm_list), "nodes": sorted(comm_list)})
            else:
                communities.append({"id": i, "size": len(comm_list), "nodes": comm_list})
        communities.sort(key=lambda c: -c["size"])
    except Exception:
        communities = []

    # Stats
    n_components = nx.number_connected_components(G)
    density = round(nx.density(G), 4) if G.number_of_nodes() > 0 else 0
    try:
        avg_clustering = round(nx.average_clustering(G), 4)
    except Exception:
        avg_clustering = 0

    return {
        "top_degree": top_degree,
        "communities": communities,
        "n_components": n_components,
        "density": density,
        "avg_clustering": avg_clustering,
    }


def institution_network_to_vis(inst_edges, institutions):
    """Convert (inst_edges, institutions) to vis-network nodes and edges with region colors."""
    nodes = []
    for inst in institutions:
        region = institution_region(inst)
        label = inst if len(inst) <= 40 else inst[:37] + "..."
        nodes.append({
            "id": inst,
            "label": label,
            "title": inst,
            "color": _region_color(region),
        })
    edges = []
    for (a, b), w in inst_edges.items():
        edges.append({"from": a, "to": b, "value": w})
    return {"nodes": nodes, "edges": edges}


# Research area presets for topic-based networks (PubMed)
RESEARCH_AREAS = [
    ("metagenomics", "Metagenomics"),
    ("CRISPR", "CRISPR"),
    ("single cell RNA", "Single cell RNA-seq"),
    ("cancer genomics", "Cancer genomics"),
    ("machine learning biology", "ML in biology"),
    ("microbiome", "Microbiome"),
]


@app.route("/")
def index():
    """Landing page with nav to PI co-author, Affiliation network, PI network."""
    return render_template("home.html")


@app.route("/pi-coauthor")
def pi_coauthor_page():
    return render_template("index.html", app="pi_coauthor")


@app.route("/affiliation")
def affiliation_page():
    return render_template("affiliation.html", research_areas=RESEARCH_AREAS)


@app.route("/pi-network")
def pi_network_page():
    return render_template("pi_network.html", research_areas=RESEARCH_AREAS)


def _run_build_scholar(scholar_id: str, queue: Queue, min_if: float = 10.0) -> None:
    """Run the build pipeline from Google Scholar: titles from Scholar, fetch from PubMed (100% title match) for authors + affiliations (region coloring)."""
    try:
        def progress(msg: str) -> None:
            queue.put({"type": "progress", "message": msg})

        profile_url, pi_name, papers, author_region = get_profile_by_scholar_id_pubmed_fast(scholar_id, progress_callback=progress, min_if=min_if)
        if not papers:
            queue.put({"type": "error", "error": "No papers matched in PubMed (100% title match). Try another source or check the profile."})
            return
        _finish_build(queue, profile_url, pi_name, papers, author_region_override=author_region)
    except Exception as e:
        queue.put({"type": "error", "error": _user_facing_error(e)})
    finally:
        queue.put(None)


def _run_build_researchgate(profile_url: str, queue: Queue) -> None:
    """Run the build pipeline from ResearchGate."""
    try:
        def progress(msg: str) -> None:
            queue.put({"type": "progress", "message": msg})

        profile_url, pi_name, papers = get_researchgate_papers(profile_url, progress_callback=progress)
        _finish_build(queue, profile_url, pi_name, papers)
    except Exception as e:
        queue.put({"type": "error", "error": str(e)})
    finally:
        queue.put(None)


def _run_build_pubmed(author_name: str, affiliation: str, queue: Queue, years: int = 5, max_results: int = 500, min_if: float = 10.0) -> None:
    """Run the build pipeline from PubMed (author search). No blocking; uses NCBI API."""
    try:
        def progress(msg: str) -> None:
            queue.put({"type": "progress", "message": msg})

        progress("Searching PubMed for author...")
        pmids = search_pi_pmids(author_name, affiliation or None, years=years)
        if max_results and len(pmids) > max_results:
            pmids = pmids[:max_results]
        progress(f"Found {len(pmids)} PMIDs. Fetching details...")
        if not pmids:
            queue.put({"type": "error", "error": f"No papers found for this author in the last {years} years. Try different name or affiliation."})
            return
        records = fetch_articles(pmids)
        records = filter_records_by_min_if(records, min_if=min_if, exclude_no_if=True)
        if not records:
            queue.put({"type": "error", "error": f"No papers left after filtering by journal impact factor (IF >= {min_if}). Papers with no IF data are excluded."})
            return
        papers = [{"title": r.get("title", ""), "authors": r.get("authors", [])} for r in records]
        profile_url = ""
        author_region = build_author_regions(records)
        _finish_build(queue, profile_url, author_name.strip(), papers, author_region_override=author_region)
    except Exception as e:
        queue.put({"type": "error", "error": str(e)})
    finally:
        queue.put(None)


def _parse_pasted_papers(text: str) -> tuple[str, list[dict]]:
    """Parse pasted text into (name_or_label, papers). Lines: 'Title | Author1, Author2' or 'Title \\t Author1; Author2'."""
    lines = [ln.strip() for ln in (text or "").strip().splitlines() if ln.strip()]
    papers = []
    for line in lines:
        if "|" in line:
            parts = line.split("|", 1)
            title = parts[0].strip()
            authors_str = parts[1].strip() if len(parts) > 1 else ""
        elif "\t" in line:
            parts = line.split("\t", 1)
            title = parts[0].strip()
            authors_str = parts[1].strip() if len(parts) > 1 else ""
        else:
            continue
        authors = []
        for sep in (";", ","):
            if sep in authors_str:
                authors = [a.strip() for a in authors_str.split(sep) if a.strip()]
                break
        if not authors and authors_str:
            authors = [a.strip() for a in authors_str.split() if a.strip()]
        if title or authors:
            papers.append({"title": title or "(no title)", "authors": authors})
    return "Pasted list", papers


def _run_build_paste(pasted_text: str, queue: Queue) -> None:
    """Run the build pipeline from pasted paper list."""
    try:
        queue.put({"type": "progress", "message": "Parsing pasted list..."})
        pi_name, papers = _parse_pasted_papers(pasted_text)
        if not papers:
            queue.put({"type": "error", "error": "No valid lines found. Use format: Title | Author1, Author2 (one paper per line)."})
            return
        _finish_build(queue, "", pi_name, papers)
    except Exception as e:
        queue.put({"type": "error", "error": str(e)})
    finally:
        queue.put(None)


def _finish_build(
    queue: Queue,
    profile_url: str,
    pi_name: str,
    papers: list,
    author_region_override: dict | None = None,
) -> None:
    queue.put({"type": "progress", "message": "Building co-author network..."})
    if not papers:
        queue.put({"type": "error", "error": "No publications on profile"})
        return
    papers = [p for p in papers if len(p.get("authors", [])) <= MAX_AUTHORS]
    co_edges, authors = build_coauthor_network_from_papers(papers)
    if len(authors) > MAX_NODES_DISPLAY:
        co_edges, authors = _top_k_nodes_by_degree(co_edges, authors, MAX_NODES_DISPLAY)
    if author_region_override is not None:
        author_region = {n: author_region_override.get(n, "Others") for n in authors}
    else:
        author_region = {n: "Others" for n in authors}
    vis = graph_to_vis(co_edges, authors, author_region)
    analysis = _compute_network_analysis(co_edges, authors)
    queue.put({
        "type": "result",
        "data": {
            "profile_url": profile_url,
            "pi_name": pi_name,
            "n_papers": len(papers),
            "n_authors": len(authors),
            "n_edges": len(co_edges),
            "nodes": vis["nodes"],
            "edges": vis["edges"],
            "papers": [{"title": p.get("title", ""), "authors": p.get("authors", [])} for p in papers],
            "analysis": analysis,
        },
    })


@app.route("/build", methods=["POST"])
def build():
    body = request.get_json() or {}
    raw = (body.get("profile_url") or body.get("scholar_url") or request.form.get("profile_url") or request.form.get("scholar_url") or "").strip()
    source = (body.get("source") or request.form.get("source") or "google_scholar").strip().lower()
    if source not in ("google_scholar", "researchgate", "pubmed", "paste"):
        source = "google_scholar"
    stream = body.get("stream", False)
    author_name = (body.get("author_name") or request.form.get("author_name") or "").strip()
    affiliation = (body.get("affiliation") or request.form.get("affiliation") or "").strip()
    pasted_text = (body.get("pasted_text") or request.form.get("pasted_text") or "").strip()
    pubmed_years = max(1, min(15, int(body.get("years", body.get("pubmed_years", 5)))))
    pubmed_max_results = max(100, min(50000, int(body.get("max_results", body.get("pubmed_max_results", 500)))))
    min_if = float(body.get("min_if", 10))
    min_if = max(0.1, min(100.0, min_if))

    if source == "pubmed":
        if not author_name:
            return jsonify({"error": "Author name required for PubMed search"}), 400
    elif source == "paste":
        if not pasted_text:
            return jsonify({"error": "Paste at least one line: Title | Author1, Author2"}), 400
    elif not raw:
        return jsonify({"error": "Profile URL required"}), 400

    if source == "researchgate":
        if not normalize_researchgate_url(raw):
            return jsonify({"error": "Invalid ResearchGate profile URL. Example: https://www.researchgate.net/profile/Firstname-Lastname"}), 400
    elif source == "google_scholar":
        scholar_id = parse_scholar_url_or_id(raw)
        if not scholar_id:
            return jsonify({"error": "Could not parse Google Scholar user ID. Use full URL or user ID (e.g. JhtCUvIAAAAJ)."}), 400

    if stream:
        def generate():
            q = Queue()
            if source == "researchgate":
                t = threading.Thread(target=_run_build_researchgate, args=(raw, q))
            elif source == "pubmed":
                t = threading.Thread(target=_run_build_pubmed, args=(author_name, affiliation, q, pubmed_years, pubmed_max_results, min_if))
            elif source == "paste":
                t = threading.Thread(target=_run_build_paste, args=(pasted_text, q))
            else:
                t = threading.Thread(target=_run_build_scholar, args=(scholar_id, q, min_if))
            t.start()
            while True:
                try:
                    item = q.get(timeout=60)
                except Empty:
                    yield json.dumps({"type": "progress", "message": "Still working…"}) + "\n"
                    continue
                if item is None:
                    break
                if item.get("type") == "error":
                    yield json.dumps({"type": "error", "error": item["error"]}) + "\n"
                    break
                yield json.dumps(item) + "\n"

        return Response(
            generate(),
            mimetype="application/x-ndjson",
            headers={"Cache-Control": "no-cache", "X-Accel-Buffering": "no"},
        )

    try:
        author_region_override = None
        if source == "researchgate":
            profile_url, pi_name, papers = get_researchgate_papers(raw)
        elif source == "pubmed":
            pmids = search_pi_pmids(author_name, affiliation or None, years=pubmed_years)
            if len(pmids) > pubmed_max_results:
                pmids = pmids[:pubmed_max_results]
            if not pmids:
                return jsonify({"error": f"No papers found for this author (last {pubmed_years} years). Try different name or affiliation."}), 404
            records = fetch_articles(pmids)
            records = filter_records_by_min_if(records, min_if=min_if, exclude_no_if=True)
            if not records:
                return jsonify({"error": f"No papers left after filtering by journal impact factor (IF >= {min_if}). Papers with no IF data are excluded."}), 404
            profile_url = ""
            pi_name = author_name
            papers = [{"title": r.get("title", ""), "authors": r.get("authors", [])} for r in records]
            author_region_override = build_author_regions(records)
        elif source == "paste":
            pi_name, papers = _parse_pasted_papers(pasted_text)
            profile_url = ""
            if not papers:
                return jsonify({"error": "No valid lines. Use: Title | Author1, Author2"}), 400
        else:
            profile_url, pi_name, papers, author_region_override = get_profile_by_scholar_id_pubmed_fast(scholar_id, min_if=min_if)
            if not papers:
                return jsonify({"error": "No papers matched in PubMed (100% title match). Try another source."}), 404
    except Exception as e:
        return jsonify({"error": _user_facing_error(e) if source == "google_scholar" else str(e)}), 502
    if not papers:
        return jsonify({"error": "No publications on profile"}), 404
    papers = [p for p in papers if len(p.get("authors", [])) <= MAX_AUTHORS]
    co_edges, authors = build_coauthor_network_from_papers(papers)
    if len(authors) > MAX_NODES_DISPLAY:
        co_edges, authors = _top_k_nodes_by_degree(co_edges, authors, MAX_NODES_DISPLAY)
    if author_region_override is not None:
        author_region = {n: author_region_override.get(n, "Others") for n in authors}
    else:
        author_region = {n: "Others" for n in authors}
    vis = graph_to_vis(co_edges, authors, author_region)
    analysis = _compute_network_analysis(co_edges, authors)
    return jsonify({
        "profile_url": profile_url,
        "pi_name": pi_name,
        "n_papers": len(papers),
        "n_authors": len(authors),
        "n_edges": len(co_edges),
        "nodes": vis["nodes"],
        "edges": vis["edges"],
        "papers": [{"title": p.get("title", ""), "authors": p.get("authors", [])} for p in papers],
        "analysis": analysis,
    })


def _run_build_affiliation(research_area: str, years: int, max_results: int, queue: Queue, min_if: float = 10.0) -> None:
    """Build institution co-affiliation network from PubMed by research area."""
    try:
        def progress(msg: str) -> None:
            queue.put({"type": "progress", "message": msg})
        progress("Searching PubMed...")
        end_year = 2025
        mindate = f"{end_year - years}/01/01"
        maxdate = f"{end_year}/12/31"
        pmids = search_pmids_by_query(research_area, mindate=mindate, maxdate=maxdate, max_results=max_results)
        if not pmids:
            queue.put({"type": "error", "error": "No papers found for this research area. Try a different term or date range."})
            return
        progress(f"Found {len(pmids)} papers. Fetching details...")
        records = fetch_pubmed_articles(pmids)
        if not records:
            queue.put({"type": "error", "error": "Could not fetch article details."})
            return
        records = filter_records_by_min_if(records, min_if=min_if, exclude_no_if=True)
        if not records:
            queue.put({"type": "error", "error": f"No papers left after filtering by journal impact factor (IF >= {min_if}). Papers with no IF data are excluded."})
            return
        progress("Building institution network (QS top 500 only)...")
        inst_edges_full, institutions_full = build_institution_network(records)
        inst_edges_full, institutions_full = filter_institution_network_to_qs_top500(inst_edges_full, institutions_full)
        analysis = _compute_network_analysis(inst_edges_full, institutions_full)
        if len(institutions_full) > MAX_NODES_DISPLAY:
            inst_edges, institutions = _top_k_nodes_by_degree(inst_edges_full, institutions_full, MAX_NODES_DISPLAY)
        else:
            inst_edges, institutions = inst_edges_full, institutions_full
        vis = institution_network_to_vis(inst_edges, institutions)
        queue.put({
            "type": "result",
            "data": {
                "profile_url": "",
                "pi_name": f"Affiliation network: {research_area}",
                "n_papers": len(records),
                "n_authors": len(institutions),
                "n_edges": len(inst_edges),
                "nodes": vis["nodes"],
                "edges": vis["edges"],
                "papers": [{"title": r.get("title", ""), "authors": [c[0] for c in r.get("author_affiliations", [])]} for r in records[:200]],
                "analysis": analysis,
            },
        })
    except Exception as e:
        queue.put({"type": "error", "error": str(e)})
    finally:
        queue.put(None)


def _run_build_pi_network(research_area: str, years: int, max_results: int, queue: Queue, min_if: float = 10.0) -> None:
    """Build co-authorship (PI/author) network from PubMed by research area."""
    try:
        def progress(msg: str) -> None:
            queue.put({"type": "progress", "message": msg})
        progress("Searching PubMed...")
        end_year = 2025
        mindate = f"{end_year - years}/01/01"
        maxdate = f"{end_year}/12/31"
        pmids = search_pmids_by_query(research_area, mindate=mindate, maxdate=maxdate, max_results=max_results)
        if not pmids:
            queue.put({"type": "error", "error": "No papers found for this research area. Try a different term or date range."})
            return
        progress(f"Found {len(pmids)} papers. Fetching details...")
        records = fetch_pubmed_articles(pmids)
        if not records:
            queue.put({"type": "error", "error": "Could not fetch article details."})
            return
        records = filter_records_by_min_if(records, min_if=min_if, exclude_no_if=True)
        if not records:
            queue.put({"type": "error", "error": f"No papers left after filtering by journal impact factor (IF >= {min_if}). Papers with no IF data are excluded."})
            return
        progress("Building PI (corresponding-author) network...")
        co_edges_full, authors_full = build_corresponding_author_network(records)
        analysis = _compute_network_analysis(co_edges_full, authors_full)
        if len(authors_full) > MAX_NODES_DISPLAY:
            co_edges, authors = _top_k_nodes_by_degree(co_edges_full, authors_full, MAX_NODES_DISPLAY)
        else:
            co_edges, authors = co_edges_full, authors_full
        author_region = build_author_regions_pubmed(records)
        author_region = {n: author_region.get(n, "Others") for n in authors}
        vis = graph_to_vis(co_edges, authors, author_region)
        queue.put({
            "type": "result",
            "data": {
                "profile_url": "",
                "pi_name": f"Author network: {research_area}",
                "n_papers": len(records),
                "n_authors": len(authors),
                "n_edges": len(co_edges),
                "nodes": vis["nodes"],
                "edges": vis["edges"],
                "papers": [{"title": r.get("title", ""), "authors": [c[0] for c in r.get("author_affiliations", [])]} for r in records[:200]],
                "analysis": analysis,
            },
        })
    except Exception as e:
        queue.put({"type": "error", "error": str(e)})
    finally:
        queue.put(None)


@app.route("/build-affiliation", methods=["POST"])
def build_affiliation():
    """Build institution co-affiliation network by research area (PubMed)."""
    body = request.get_json() or {}
    research_area = (body.get("research_area") or body.get("custom_query") or "metagenomics").strip()
    if not research_area:
        return jsonify({"error": "Research area or custom query required"}), 400
    years = int(body.get("years", 5))
    years = max(1, min(10, years))
    max_results = int(body.get("max_results", 500))
    max_results = max(100, min(2000, max_results))
    min_if = float(body.get("min_if", 10))
    min_if = max(0.1, min(100.0, min_if))
    stream = body.get("stream", True)

    if stream:
        def generate():
            q = Queue()
            t = threading.Thread(target=_run_build_affiliation, args=(research_area, years, max_results, q, min_if))
            t.start()
            while True:
                try:
                    item = q.get(timeout=120)
                except Empty:
                    yield json.dumps({"type": "progress", "message": "Still working…"}) + "\n"
                    continue
                if item is None:
                    break
                if item.get("type") == "error":
                    yield json.dumps({"type": "error", "error": item["error"]}) + "\n"
                    break
                yield json.dumps(item) + "\n"
        return Response(
            generate(),
            mimetype="application/x-ndjson",
            headers={"Cache-Control": "no-cache", "X-Accel-Buffering": "no"},
        )
    # Non-stream: run synchronously
    end_year = 2025
    mindate = f"{end_year - years}/01/01"
    maxdate = f"{end_year}/12/31"
    pmids = search_pmids_by_query(research_area, mindate=mindate, maxdate=maxdate, max_results=max_results)
    if not pmids:
        return jsonify({"error": "No papers found for this research area."}), 404
    records = fetch_pubmed_articles(pmids)
    if not records:
        return jsonify({"error": "Could not fetch article details."}), 502
    records = filter_records_by_min_if(records, min_if=min_if, exclude_no_if=True)
    if not records:
        return jsonify({"error": f"No papers left after filtering by journal impact factor (IF >= {min_if}). Papers with no IF data are excluded."}), 404
    inst_edges_full, institutions_full = build_institution_network(records)
    inst_edges_full, institutions_full = filter_institution_network_to_qs_top500(inst_edges_full, institutions_full)
    analysis = _compute_network_analysis(inst_edges_full, institutions_full)
    if len(institutions_full) > MAX_NODES_DISPLAY:
        inst_edges, institutions = _top_k_nodes_by_degree(inst_edges_full, institutions_full, MAX_NODES_DISPLAY)
    else:
        inst_edges, institutions = inst_edges_full, institutions_full
    vis = institution_network_to_vis(inst_edges, institutions)
    return jsonify({
        "profile_url": "",
        "pi_name": f"Affiliation network: {research_area}",
        "n_papers": len(records),
        "n_authors": len(institutions),
        "n_edges": len(inst_edges),
        "nodes": vis["nodes"],
        "edges": vis["edges"],
        "papers": [{"title": r.get("title", ""), "authors": [c[0] for c in r.get("author_affiliations", [])]} for r in records[:200]],
        "analysis": analysis,
    })


@app.route("/build-pi-network", methods=["POST"])
def build_pi_network():
    """Build co-authorship (author) network by research area (PubMed)."""
    body = request.get_json() or {}
    research_area = (body.get("research_area") or body.get("custom_query") or "metagenomics").strip()
    if not research_area:
        return jsonify({"error": "Research area or custom query required"}), 400
    years = int(body.get("years", 5))
    years = max(1, min(10, years))
    max_results = int(body.get("max_results", 500))
    max_results = max(100, min(50000, max_results))
    min_if = float(body.get("min_if", 10))
    min_if = max(0.1, min(100.0, min_if))
    stream = body.get("stream", True)

    if stream:
        def generate():
            q = Queue()
            t = threading.Thread(target=_run_build_pi_network, args=(research_area, years, max_results, q, min_if))
            t.start()
            while True:
                try:
                    item = q.get(timeout=120)
                except Empty:
                    yield json.dumps({"type": "progress", "message": "Still working…"}) + "\n"
                    continue
                if item is None:
                    break
                if item.get("type") == "error":
                    yield json.dumps({"type": "error", "error": item["error"]}) + "\n"
                    break
                yield json.dumps(item) + "\n"
        return Response(
            generate(),
            mimetype="application/x-ndjson",
            headers={"Cache-Control": "no-cache", "X-Accel-Buffering": "no"},
        )
    end_year = 2025
    mindate = f"{end_year - years}/01/01"
    maxdate = f"{end_year}/12/31"
    pmids = search_pmids_by_query(research_area, mindate=mindate, maxdate=maxdate, max_results=max_results)
    if not pmids:
        return jsonify({"error": "No papers found for this research area."}), 404
    records = fetch_pubmed_articles(pmids)
    if not records:
        return jsonify({"error": "Could not fetch article details."}), 502
    records = filter_records_by_min_if(records, min_if=min_if, exclude_no_if=True)
    if not records:
        return jsonify({"error": f"No papers left after filtering by journal impact factor (IF >= {min_if}). Papers with no IF data are excluded."}), 404
    co_edges_full, authors_full = build_corresponding_author_network(records)
    analysis = _compute_network_analysis(co_edges_full, authors_full)
    if len(authors_full) > MAX_NODES_DISPLAY:
        co_edges, authors = _top_k_nodes_by_degree(co_edges_full, authors_full, MAX_NODES_DISPLAY)
    else:
        co_edges, authors = co_edges_full, authors_full
    author_region = build_author_regions_pubmed(records)
    author_region = {n: author_region.get(n, "Others") for n in authors}
    vis = graph_to_vis(co_edges, authors, author_region)
    return jsonify({
        "profile_url": "",
        "pi_name": f"Author network: {research_area}",
        "n_papers": len(records),
        "n_authors": len(authors),
        "n_edges": len(co_edges),
        "nodes": vis["nodes"],
        "edges": vis["edges"],
        "papers": [{"title": r.get("title", ""), "authors": [c[0] for c in r.get("author_affiliations", [])]} for r in records[:200]],
        "analysis": analysis,
    })


if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5050, debug=True)
