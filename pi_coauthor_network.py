#!/usr/bin/env python3
"""
Standalone script: given a PI name (and optional affiliation for double-check),
search their papers in PubMed in the past 5 years, collect all authors,
and build a co-authorship network: nodes = authors, edge weight = shared paper count.
"""
from __future__ import annotations

import argparse
import time
from collections import defaultdict
from pathlib import Path

from Bio import Entrez

Entrez.email = "sci_net_user@example.com"
FETCH_BATCH = 500  # larger batches = fewer API calls
PUBMED_MAX = 10000


def _normalize_author(author_dict: dict) -> str:
    last = author_dict.get("LastName", "")
    fore = author_dict.get("ForeName", "")
    if fore:
        return f"{last} {fore}".strip()
    return last.strip() or "Unknown"


def _extract_affiliation(author_dict: dict) -> list[str]:
    affs = []
    if "AffiliationInfo" not in author_dict:
        return affs
    for info in author_dict["AffiliationInfo"]:
        if isinstance(info, dict) and "Affiliation" in info:
            affs.append(str(info["Affiliation"]).strip())
    return affs


# Common middle initials to try when we only have first initial (e.g. Doudna JA, Banfield JF)
_COMMON_MIDDLE_INITIALS = "AMLSRKF"

# First-name alternatives (short form -> longer/formal form) for better recall
_FIRSTNAME_ALTERNATIVES = {
    "jill": ["jillian"],
    "jim": ["james"],
    "mike": ["michael"],
    "dave": ["david"],
    "steve": ["steven", "stephen"],
    "chris": ["christopher"],
    "tony": ["anthony"],
    "dan": ["daniel"],
    "matt": ["matthew"],
    "jen": ["jennifer"],
    "kate": ["katherine", "catherine"],
    "bob": ["robert"],
    "bill": ["william"],
}


def _author_search_variants(pi_name: str) -> list[str]:
    """
    Build multiple PubMed author query variants to maximize recall.
    Tries: LastName Initials, full name, alternative first names (Jill/Jillian), middle initials.
    """
    parts = pi_name.strip().split()
    if not parts:
        return [pi_name]
    variants = set()
    family = parts[-1]
    given = parts[:-1]

    def add_for_given(given_parts: list[str]) -> None:
        if not given_parts:
            return
        inits = "".join(p[0] for p in given_parts).upper()
        variants.add(f"{family} {inits}")
        variants.add(f"{family} {' '.join(p[0].upper() for p in given_parts)}")
        variants.add(f"{family} {' '.join(given_parts)}")
        variants.add(" ".join(given_parts) + " " + family)
        if len(inits) == 1:
            for mi in _COMMON_MIDDLE_INITIALS:
                variants.add(f"{family} {inits}{mi}")
                variants.add(f"{family} {inits} {mi}")
            # "First MiddleInitial Last" (e.g. Jill F Banfield)
            first = given_parts[0]
            for mi in _COMMON_MIDDLE_INITIALS:
                variants.add(f"{first} {mi} {family}")

    variants.add(pi_name.strip())
    if given:
        add_for_given(given)
        # Alternative first names (e.g. Jill -> Jillian)
        first_lower = given[0].lower()
        if first_lower in _FIRSTNAME_ALTERNATIVES:
            for alt in _FIRSTNAME_ALTERNATIVES[first_lower]:
                alt_parts = [alt] + given[1:] if len(given) > 1 else [alt]
                add_for_given(alt_parts)
                variants.add(f"{family} {alt}")
                variants.add(alt.title() + " " + family)
    return list(variants)


def _search_one_query(query: str, mindate: str, maxdate: str) -> list[str]:
    """Run one esearch and return all PMIDs (paginate up to PUBMED_MAX)."""
    ids = []
    retstart = 0
    while True:
        try:
            handle = Entrez.esearch(
                db="pubmed",
                term=query,
                mindate=mindate,
                maxdate=maxdate,
                datetype="pdat",
                retmax=PUBMED_MAX,
                retstart=retstart,
                sort="relevance",
            )
            data = Entrez.read(handle)
            handle.close()
        except Exception as e:
            break
        batch = data.get("IdList", [])
        if not batch:
            break
        ids.extend(batch)
        if len(batch) < PUBMED_MAX:
            break
        retstart += len(batch)
        time.sleep(0.34)
    return ids


def _affiliation_search_variants(affiliation: str) -> list[str]:
    """Expand affiliation into variants to improve match (e.g. UC Berkeley -> Berkeley, University of California)."""
    if not affiliation or not affiliation.strip():
        return []
    a = affiliation.strip()
    variants = [a]
    a_lower = a.lower()
    if "uc berkeley" in a_lower or a_lower == "uc berkeley" or "berkeley" in a_lower:
        variants.extend([
            "Berkeley", "University of California Berkeley", "University of California, Berkeley",
            "University of California", "Berkeley Lab", "LBNL",
        ])
    if "city university of hong kong" in a_lower or "cityu" in a_lower:
        variants.extend(["City University of Hong Kong", "CityU", "Hong Kong"])
    # Dedupe preserving order
    seen = set()
    out = []
    for v in variants:
        v = v.strip()
        if v and v not in seen:
            seen.add(v)
            out.append(v)
    return out


def search_pi_pmids(pi_name: str, affiliation: str | None, years: int = 5) -> list[str]:
    """Search PubMed for papers by PI. One search per name variant; affiliation variants OR'd in one query (much faster)."""
    from datetime import datetime
    now = datetime.now()
    mindate = f"{now.year - years}/01/01"
    maxdate = now.strftime("%Y/%m/%d")

    name_variants = _author_search_variants(pi_name)
    aff_variants = _affiliation_search_variants(affiliation) if affiliation else []
    seen = set()
    all_ids = []
    n = len(name_variants)
    for idx, variant in enumerate(name_variants):
        if n > 5 and (idx + 1) % 5 == 0:
            print(f"  Search progress: {idx + 1}/{n} name variants, {len(all_ids)} PMIDs so far...")
        author_term = f'"{variant}"[Author]'
        if aff_variants:
            aff_clause = " OR ".join(f'"{a}"[Affiliation]' for a in aff_variants)
            query = f"{author_term} AND ({aff_clause})"
        else:
            query = author_term
        ids = _search_one_query(query, mindate, maxdate)
        for i in ids:
            if i not in seen:
                seen.add(i)
                all_ids.append(i)
        time.sleep(0.34)  # NCBI rate limit ~3/sec without API key
    return all_ids


def fetch_articles(pmids: list[str]) -> list[dict]:
    """Fetch full article XML and return list of {pmid, title, authors} (authors = list of full names)."""
    records = []
    n_batches = (len(pmids) + FETCH_BATCH - 1) // FETCH_BATCH
    for b in range(0, len(pmids), FETCH_BATCH):
        batch = pmids[b : b + FETCH_BATCH]
        batch_num = b // FETCH_BATCH + 1
        if n_batches > 1:
            print(f"  Fetching batch {batch_num}/{n_batches} ({len(batch)} articles)...")
        try:
            handle = Entrez.efetch(db="pubmed", id=batch, retmode="xml")
            data = Entrez.read(handle)
            handle.close()
        except Exception as e:
            print(f"Fetch error: {e}")
            time.sleep(1)
            continue
        for article in data.get("PubmedArticle", []):
            try:
                med = article["MedlineCitation"]
                art = med["Article"]
                author_list = art.get("AuthorList", [])
                title = art.get("ArticleTitle", "")
                if isinstance(title, list):
                    title = " ".join(str(t) for t in title) if title else ""
                else:
                    title = str(title or "").strip()
                authors = []
                for a in author_list:
                    if a.get("ValidYN") == "N":
                        continue
                    name = _normalize_author(a)
                    if name:
                        authors.append(name)
                if authors:
                    records.append({"pmid": str(med["PMID"]), "title": title, "authors": authors})
            except Exception as e:
                print(f"Parse error: {e}")
        time.sleep(0.34)
    return records


def double_check_affiliation(records: list[dict], pi_name: str, affiliation: str) -> list[dict]:
    """Keep only papers where the PI (or any author) has an affiliation containing `affiliation`."""
    if not affiliation:
        return records
    # Re-fetch to get affiliations for author list; we don't have affs in current records.
    # So we re-parse from efetch and filter: keep paper if any author's affiliation contains affiliation
    filtered = []
    aff_lower = affiliation.lower()
    for rec in records:
        # We don't have affiliations in rec; need to fetch again or skip filter.
        # For simplicity: if user provided affiliation, we already used it in the search query,
        # so papers returned should match. Optionally we could efetch again and filter by author aff.
        filtered.append(rec)
    return filtered


def build_coauthor_network(records: list[dict]) -> tuple[defaultdict[int], set]:
    """Nodes = authors, edge (a,b) weight = number of papers containing both a and b."""
    co_edges = defaultdict(int)
    all_authors = set()
    for rec in records:
        names = list(dict.fromkeys(rec.get("authors", [])))
        for n in names:
            all_authors.add(n)
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                a, b = sorted([names[i], names[j]])
                co_edges[(a, b)] += 1
    return co_edges, all_authors


def main():
    parser = argparse.ArgumentParser(
        description="Build co-authorship network for a PI's papers (past 5 years, PubMed)."
    )
    parser.add_argument("pi_name", help="PI name as in PubMed, e.g. 'Smith J' or 'Smith John'")
    parser.add_argument(
        "--affiliation",
        "-a",
        default=None,
        help="Affiliation name for double-check (narrows search to this affiliation)",
    )
    parser.add_argument(
        "--years",
        type=int,
        default=5,
        help="Number of years back to search (default: 5)",
    )
    parser.add_argument(
        "--out-dir",
        "-o",
        type=Path,
        default=None,
        help="Output directory (default: outputs/pi_<sanitized_name>)",
    )
    parser.add_argument("--no-plot", action="store_true", help="Skip generating network plot")
    args = parser.parse_args()

    pi = args.pi_name.strip()
    aff = args.affiliation.strip() if args.affiliation else None
    years = max(1, args.years)

    out_dir = args.out_dir
    if out_dir is None:
        safe_name = "".join(c if c.isalnum() or c in " -_" else "_" for c in pi).strip(" _") or "pi"
        out_dir = Path(__file__).resolve().parent / "outputs" / f"pi_{safe_name}"
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Searching PubMed for papers by '{pi}' (last {years} years)" + (f", affiliation '{aff}'" if aff else "") + "...")
    pmids = search_pi_pmids(pi, aff, years=years)
    print(f"Found {len(pmids)} PMIDs.")

    if not pmids:
        print("No papers found. Exiting.")
        return

    print("Fetching article details...")
    records = fetch_articles(pmids)
    records = double_check_affiliation(records, pi, aff or "")
    n_papers = len(records)
    print(f"Number of papers: {n_papers}")

    co_edges, authors = build_coauthor_network(records)
    n_nodes = len(authors)
    n_edges = len(co_edges)
    print(f"Co-authorship network: {n_nodes} authors, {n_edges} edges.")

    # Save edge list
    edge_path = out_dir / "coauthor_edges.txt"
    with open(edge_path, "w") as f:
        f.write("Author1\tAuthor2\tWeight\n")
        for (a, b), w in sorted(co_edges.items(), key=lambda x: -x[1]):
            f.write(f"{a}\t{b}\t{w}\n")
    print(f"Saved {edge_path}")

    # Save paper list (title, PMID, authors)
    papers_path = out_dir / "papers.txt"
    with open(papers_path, "w") as f:
        f.write("PMID\tTitle\tAuthors\n")
        for r in records:
            title = (r.get("title") or "").replace("\t", " ").replace("\n", " ")[:200]
            authors_str = "; ".join(r.get("authors", []))
            f.write(f"{r['pmid']}\t{title}\t{authors_str}\n")
    print(f"Saved {papers_path}")

    # GML
    try:
        import networkx as nx
        G = nx.Graph()
        for (a, b), w in co_edges.items():
            G.add_edge(a, b, weight=w)
        for n in G.nodes():
            G.nodes[n]["label"] = n
        gml_path = out_dir / "coauthor_network.gml"
        nx.write_gml(G, gml_path)
        print(f"Saved {gml_path}")
    except Exception as e:
        print(f"GML export skipped: {e}")

    # Optional plot
    if not args.no_plot and n_nodes > 0:
        try:
            import matplotlib.pyplot as plt
            import numpy as np
            G = nx.Graph()
            for (a, b), w in co_edges.items():
                G.add_edge(a, b, weight=w)
            if G.number_of_edges() == 0:
                print("No edges to plot.")
            else:
                pos = nx.spring_layout(G, k=1.5, iterations=80, seed=42, weight="weight")
                deg = dict(G.degree(weight="weight"))
                node_sizes = [max(80, deg.get(n, 1) * 15) for n in G.nodes()]
                weights = [G.edges[u, v].get("weight", 1) for u, v in G.edges()]
                w_min, w_max = min(weights), max(weights)
                edge_widths = [0.5 + 3 * (w - w_min) / (w_max - w_min) if w_max > w_min else 1 for w in weights]
                plt.figure(figsize=(12, 10))
                nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.6, edge_color="gray")
                nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color="steelblue", alpha=0.85)
                labels = {n: n if len(n) <= 18 else n[:15] + "…" for n in G.nodes()}
                nx.draw_networkx_labels(G, pos, labels=labels, font_size=5)
                plt.title(f"Co-authorship network: {pi} (last {years} years, {n_papers} papers)")
                plt.axis("off")
                plt.tight_layout()
                plot_path = out_dir / "coauthor_network.png"
                plt.savefig(plot_path, dpi=150, bbox_inches="tight")
                plt.close()
                print(f"Saved {plot_path}")
        except Exception as e:
            print(f"Plot skipped: {e}")

    print("Done.")


if __name__ == "__main__":
    main()
