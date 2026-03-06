#!/usr/bin/env python3
"""
Standalone script: given a PI name (and optional affiliation for double-check),
search their papers in PubMed in the past 5 years, collect all authors,
and build a co-authorship network: nodes = authors, edge weight = shared paper count.
"""
from __future__ import annotations

import argparse
import re
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


# Region keywords for coloring nodes by affiliation (same as main project)
_REGION_KEYWORDS = {
    "China": ["China", "Chinese", "Beijing", "Shanghai", "Hong Kong", "Fudan", "Tsinghua", "Zhejiang", "Wuhan", "Nanjing", "CAS", "Peking", "Huazhong", "Nankai"],
    "US": ["USA", "U.S.A", "United States", "America", "Harvard", "Stanford", "MIT", "Yale", "Columbia", "California", "Berkeley", "Michigan", "Cornell", "Duke", "Texas", "Washington", "Chicago", "Boston", "U.S."],
    "UK": ["UK", "U.K", "United Kingdom", "England", "Oxford", "Cambridge", "London", "Edinburgh", "Imperial", "UCL", "Bristol", "Leeds", "Glasgow"],
    "Europe": ["Germany", "French", "France", "Netherlands", "Switzerland", "Sweden", "Denmark", "Italy", "Spain", "Munich", "Berlin", "Paris", "Amsterdam", "Zurich", "ETH", "Max Planck", "CNRS", "INSERM", "Karolinska", "Vienna"],
    "Japan": ["Japan", "Japanese", "Tokyo", "Kyoto", "Osaka", "Waseda", "Tohoku", "Riken"],
    "Australia": ["Australia", "Australian", "Sydney", "Melbourne", "Brisbane", "Perth", "Adelaide", "Canberra", "Queensland", "New South Wales", "Monash", "UNSW", "University of Sydney", "University of Melbourne"],
    "Canada": ["Canada", "Canadian", "Toronto", "Vancouver", "Montreal", "McGill", "British Columbia", "UBC", "Alberta", "Ontario", "Quebec", "Calgary", "Waterloo", "University of Toronto"],
}


def _affiliation_to_region(aff_text: str) -> str:
    """Infer region from affiliation string. Returns China, US, UK, Europe, Japan, Australia, Canada, or Others."""
    if not aff_text:
        return "Others"
    lower = aff_text.lower()
    for region, keywords in _REGION_KEYWORDS.items():
        for kw in keywords:
            if kw.lower() in lower:
                return region
    return "Others"


_REGION_ORDER = ("China", "US", "UK", "Europe", "Japan", "Australia", "Canada", "Others")


def _region_color(region: str) -> str:
    """Distinct colors per region for node coloring."""
    palette = {
        "China": "#e74c3c",
        "US": "#3498db",
        "UK": "#2ecc71",
        "Europe": "#9b59b6",
        "Japan": "#f39c12",
        "Australia": "#1abc9c",
        "Canada": "#e67e22",
        "Others": "#95a5a6",
    }
    return palette.get(region, palette["Others"])


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


def _split_camel_given(token: str) -> list[str] | None:
    """If token looks like two words in one (e.g. PengFan, Pengfan), return ['Peng', 'Fan']. Else None."""
    if len(token) < 3:
        return None
    # Split on internal uppercase: PengFan -> Peng, Fan
    parts = []
    current = []
    for i, c in enumerate(token):
        if c.isupper() and i > 0 and current:
            parts.append("".join(current))
            current = [c]
        else:
            current.append(c)
    if current:
        parts.append("".join(current))
    if len(parts) >= 2 and all(len(p) >= 2 for p in parts):
        return [p[0].upper() + p[1:].lower() for p in parts]
    return None


def _author_search_variants(pi_name: str) -> list[str]:
    """
    Build multiple PubMed author query variants to maximize recall.
    Tries: LastName Initials, full name, alternative first names (Jill/Jillian), middle initials.
    Ensures compound given names (e.g. PengFan) get both one-word and two-word forms (Pengfan Zhang, Peng Fan Zhang, Zhang Pengfan).
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
        # One-word given name that might be written as two words (e.g. PengFan -> Peng Fan)
        if len(given) == 1:
            token = given[0]
            split = _split_camel_given(token)
            if split:
                add_for_given(split)
            # Title-case form so we have "Pengfan Zhang" and "Zhang Pengfan" (no missing word)
            title_given = token[0].upper() + token[1:].lower() if len(token) > 1 else token
            if title_given != token:
                add_for_given([title_given])
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


def _pubmed_pmid_for_title(title: str) -> str | None:
    """Search PubMed by title; return first PMID if found, else None."""
    if not title or not title.strip():
        return None
    # Clean title for search: remove extra spaces, truncate very long
    t = " ".join(title.strip().split())[:200]
    if not t:
        return None
    try:
        handle = Entrez.esearch(db="pubmed", term=f'"{t}"[Title]', retmax=1)
        data = Entrez.read(handle)
        handle.close()
        ids = data.get("IdList", [])
        return ids[0] if ids else None
    except Exception:
        return None


def parse_scholar_url_or_id(s: str) -> str | None:
    """
    Extract Google Scholar user ID from a profile URL or return the string if it's already an ID.
    E.g. 'https://scholar.google.com/citations?user=JhtCUvIAAAAJ&hl=en' -> 'JhtCUvIAAAAJ'
    """
    s = (s or "").strip()
    if not s:
        return None
    if "scholar.google" in s and "user=" in s:
        try:
            from urllib.parse import urlparse, parse_qs
            parsed = urlparse(s)
            q = parse_qs(parsed.query)
            ids = q.get("user", [])
            if ids:
                return ids[0].strip()
        except Exception:
            pass
    # Treat as raw user ID (alphanumeric and maybe underscores)
    if s and all(c.isalnum() or c in "_-" for c in s):
        return s
    return None


def get_profile_by_scholar_id(
    scholar_id: str,
    progress_callback: None = None,
) -> tuple[str, str, list[dict]]:
    """
    Fetch Google Scholar profile by user ID. Returns (profile_url, author_name, list of {title, authors}).
    Co-authors are taken directly from each publication on Scholar (no PubMed).
    If progress_callback is provided, call it with status messages (e.g. "Loaded 5/13 publications...").
    """
    from scholarly import scholarly

    def report(msg: str) -> None:
        if progress_callback:
            progress_callback(msg)
        else:
            print(msg)

    report("Fetching Google Scholar profile...")
    author = scholarly.search_author_id(scholar_id.strip())
    author = scholarly.fill(author, sections=["basics", "publications"])
    sid = (author.get("scholar_id") if isinstance(author, dict) else getattr(author, "scholar_id", None)) or ""
    profile_url = f"https://scholar.google.com/citations?user={sid}" if sid else ""
    name = (author.get("name") if isinstance(author, dict) else getattr(author, "name", None)) or "Unknown"
    pubs = author.get("publications", []) if isinstance(author, dict) else getattr(author, "publications", [])
    n_pubs = len(pubs or [])
    report(f"Found {n_pubs} publications. Loading details...")
    papers = []
    for i, pub in enumerate(pubs or []):
        try:
            filled = scholarly.fill(pub)
            bib = filled.get("bib", {}) if isinstance(filled, dict) else getattr(filled, "bib", {})
            if not isinstance(bib, dict):
                continue
            title = (bib.get("title") or "").strip()
            author_list = bib.get("author") or []
            if isinstance(author_list, str):
                author_list = [a.strip() for a in author_list.split(" and ") if a.strip()]
            else:
                author_list = [str(a).strip() for a in author_list if a]
            authors = author_list
            if title or authors:
                papers.append({"title": title or "(no title)", "authors": authors})
            report(f"Loaded {len(papers)}/{n_pubs} publications...")
            time.sleep(0.5)
        except Exception as e:
            continue
    report(f"Done. {len(papers)} publications with author lists.")
    return profile_url, name, papers


def _normalize_title_for_match(title: str) -> str:
    """Normalize title for 100% match: strip, collapse whitespace, lowercase."""
    if not title:
        return ""
    s = re.sub(r"\s+", " ", str(title).strip()).lower()
    return s


def get_scholar_titles_only(
    scholar_id: str,
    progress_callback: None = None,
) -> tuple[str, str, list[str]]:
    """
    Fetch Google Scholar profile and return only publication titles (no per-pub fill for authors).
    Returns (profile_url, author_name, list of title strings).
    """
    from scholarly import scholarly

    def report(msg: str) -> None:
        if progress_callback:
            progress_callback(msg)
        else:
            print(msg)

    report("Fetching Google Scholar profile (titles only)...")
    author = scholarly.search_author_id(scholar_id.strip())
    author = scholarly.fill(author, sections=["basics", "publications"])
    sid = (author.get("scholar_id") if isinstance(author, dict) else getattr(author, "scholar_id", None)) or ""
    profile_url = f"https://scholar.google.com/citations?user={sid}" if sid else ""
    name = (author.get("name") if isinstance(author, dict) else getattr(author, "name", None)) or "Unknown"
    pubs = author.get("publications", []) if isinstance(author, dict) else getattr(author, "publications", [])
    titles = []
    for pub in pubs or []:
        bib = pub.get("bib", {}) if isinstance(pub, dict) else getattr(pub, "bib", None)
        if isinstance(bib, dict):
            t = (bib.get("title") or "").strip()
        else:
            t = (getattr(bib, "title", None) or "").strip() if bib else ""
        if not t and pub:
            try:
                filled = scholarly.fill(pub)
                bib = filled.get("bib", {}) if isinstance(filled, dict) else getattr(filled, "bib", {})
                if isinstance(bib, dict):
                    t = (bib.get("title") or "").strip()
                if t:
                    report(f"Got title for paper {len(titles) + 1}...")
                time.sleep(0.35)
            except Exception:
                pass
        if t:
            titles.append(t)
    report(f"Found {len(titles)} publication titles from Scholar.")
    return profile_url, name, titles


def _search_pmids_by_title(title: str) -> list[str]:
    """Search PubMed by exact title phrase; return list of PMIDs (may be empty)."""
    if not title or not title.strip():
        return []
    query = f'"{title.strip()}"[Title]'
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=20)
        data = Entrez.read(handle)
        handle.close()
        return list(data.get("IdList", []))
    except Exception:
        return []


def fetch_and_match_records_by_titles(
    titles: list[str],
    progress_callback: None = None,
    min_if: float = 10.0,
) -> list[dict]:
    """
    For each title, search PubMed and fetch articles; keep only records with 100% title match.
    Filters by journal IF >= min_if; excludes papers with no IF data.
    Returns list of records with pmid, title, authors, author_affiliations (for region coloring).
    """
    def report(msg: str) -> None:
        if progress_callback:
            progress_callback(msg)
        else:
            print(msg)

    want_normalized = {_normalize_title_for_match(t) for t in titles if t}
    if not want_normalized:
        return []

    report("Searching PubMed for each title...")
    pmids_to_fetch = []
    for i, title in enumerate(titles):
        if not title or not title.strip():
            continue
        pmids = _search_pmids_by_title(title)
        for pid in pmids:
            if pid not in pmids_to_fetch:
                pmids_to_fetch.append(pid)
        time.sleep(0.34)
        if progress_callback and (i + 1) % 10 == 0:
            report(f"Searched {i + 1}/{len(titles)} titles...")

    if not pmids_to_fetch:
        report("No PMIDs found for these titles.")
        return []

    report(f"Fetching {len(pmids_to_fetch)} articles from PubMed...")
    records = fetch_articles(pmids_to_fetch)
    from journal_if import filter_records_by_min_if

    records = filter_records_by_min_if(records, min_if=min_if, exclude_no_if=True)
    if records:
        report(f"After filtering by journal IF >= {min_if}: {len(records)} papers.")

    # Keep only 100% title match
    matched = []
    for rec in records:
        t = (rec.get("title") or "").strip()
        norm = _normalize_title_for_match(t)
        if norm in want_normalized:
            matched.append(rec)
    report(f"Matched {len(matched)}/{len(records)} papers by exact title.")
    return matched


def get_profile_by_scholar_id_pubmed_fast(
    scholar_id: str,
    progress_callback: None = None,
    min_if: float = 10.0,
) -> tuple[str, str, list[dict], dict[str, str]]:
    """
    Fast path: get publication titles from Google Scholar (minimal calls), then fetch
    each paper from PubMed with 100% title match. Returns (profile_url, pi_name, papers, author_region).
    papers = list of {title, authors}; author_region from PubMed affiliations for region coloring.
    Filters by journal IF >= min_if; excludes papers with no IF data.
    """
    profile_url, name, titles = get_scholar_titles_only(scholar_id, progress_callback=progress_callback)
    if not titles:
        return profile_url, name, [], {}

    records = fetch_and_match_records_by_titles(titles, progress_callback=progress_callback, min_if=min_if)
    if not records:
        return profile_url, name, [], {}

    papers = [
        {
            "title": r.get("title", ""),
            "authors": [a for a, _ in r.get("author_affiliations", [])],
        }
        for r in records
    ]
    author_region = build_author_regions(records)
    return profile_url, name, papers, author_region


def build_coauthor_network_from_papers(papers: list[dict]) -> tuple[defaultdict, set]:
    """Build (co_edges, all_authors) from list of {title, authors} (e.g. from Google Scholar)."""
    co_edges = defaultdict(int)
    all_authors = set()
    for rec in papers:
        names = list(dict.fromkeys(rec.get("authors", [])))
        for n in names:
            all_authors.add(n)
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                a, b = sorted([names[i], names[j]])
                co_edges[(a, b)] += 1
    return co_edges, all_authors


def _single_affiliation_variants(a: str) -> list[str]:
    """Expand a single affiliation string into variants for PubMed search."""
    a = a.strip()
    if not a:
        return []
    variants = [a]
    a_lower = a.lower()
    if "uc berkeley" in a_lower or a_lower == "uc berkeley" or "berkeley" in a_lower:
        variants.extend([
            "Berkeley", "University of California Berkeley", "University of California, Berkeley",
            "University of California", "Berkeley Lab", "LBNL",
        ])
    if "city university of hong kong" in a_lower or "cityu" in a_lower:
        variants.extend(["City University of Hong Kong", "CityU", "Hong Kong"])
    if "ubc" in a_lower or "university of british columbia" in a_lower:
        variants.extend(["UBC", "University of British Columbia", "British Columbia"])
    if "oxford" in a_lower:
        variants.extend(["Oxford", "University of Oxford", "Oxford University"])
    if "westlake" in a_lower:
        variants.extend(["Westlake University", "Westlake", "Westlake U"])
    if "zhejiang" in a_lower:
        variants.extend(["Zhejiang University", "Zhejiang", "Hangzhou"])
    if "innovative genomics" in a_lower or "igi" in a_lower:
        variants.extend(["Innovative Genomics Institute", "IGI", "Innovative Genomics"])
    if "osu" in a_lower or "ohio state" in a_lower:
        variants.extend(["Ohio State", "Ohio State University", "OSU", "The Ohio State University"])
    if "oregon state" in a_lower:
        variants.extend(["Oregon State", "Oregon State University", "OSU"])
    return variants


def _affiliation_search_variants(affiliation: str) -> list[str]:
    """Expand affiliation into variants. Supports multiple institutions: 'UBC, Oxford' or 'UBC and Oxford'."""
    if not affiliation or not affiliation.strip():
        return []
    raw = affiliation.strip()
    # Split on " and " or "," to allow multiple institutions
    parts = [p.strip() for p in raw.replace(" and ", ",").split(",") if p.strip()]
    seen = set()
    out = []
    for part in parts:
        for v in _single_affiliation_variants(part):
            v = v.strip()
            if v and v not in seen:
                seen.add(v)
                out.append(v)
    return out if out else [raw]


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
                journal = ""
                j = art.get("Journal", {}) or {}
                if isinstance(j, dict):
                    journal = (j.get("Title") or j.get("ISOAbbreviation") or "")
                if isinstance(journal, list):
                    journal = " ".join(str(x) for x in journal) if journal else ""
                journal = str(journal or "").strip()
                authors = []
                author_affs = []  # (name, [affiliation strings])
                for a in author_list:
                    if a.get("ValidYN") == "N":
                        continue
                    name = _normalize_author(a)
                    if name:
                        authors.append(name)
                        author_affs.append((name, _extract_affiliation(a)))
                if authors:
                    records.append({"pmid": str(med["PMID"]), "title": title, "journal": journal, "authors": authors, "author_affiliations": author_affs})
            except Exception as e:
                print(f"Parse error: {e}")
        time.sleep(0.34)
    return records


def _pi_author_match_set(pi_name: str) -> set[str]:
    """Set of author-name strings that we consider the PI (for affiliation check)."""
    return set(_author_search_variants(pi_name))


def _affiliation_match_keywords(affiliation: str) -> set[str]:
    """All affiliation variant strings (lowercase) to match against author affiliations."""
    keywords = set()
    for part in (affiliation or "").split(","):
        part = part.strip()
        if not part:
            continue
        for v in _single_affiliation_variants(part):
            if v:
                keywords.add(v.lower())
    return keywords


def double_check_affiliation(records: list[dict], pi_name: str, affiliation: str) -> list[dict]:
    """Keep only papers where the PI appears as author AND has an affiliation matching the given one."""
    if not affiliation or not affiliation.strip():
        return records
    pi_names = _pi_author_match_set(pi_name)
    aff_keywords = _affiliation_match_keywords(affiliation)
    filtered = []
    for rec in records:
        author_affs = rec.get("author_affiliations", [])
        if not author_affs:
            filtered.append(rec)
            continue
        pi_has_match = False
        for name, affs in author_affs:
            if name not in pi_names:
                continue
            aff_text = " ".join(str(a) for a in affs).lower()
            if any(kw in aff_text for kw in aff_keywords):
                pi_has_match = True
                break
        if pi_has_match:
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


def build_author_regions(records: list[dict]) -> dict[str, str]:
    """From author_affiliations in records, infer region per author (majority vote)."""
    from collections import Counter
    author_affs = defaultdict(list)
    for rec in records:
        for name, affs in rec.get("author_affiliations", []):
            if name:
                author_affs[name].extend(affs)
    author_region = {}
    for author, affs in author_affs.items():
        regions = [_affiliation_to_region(a) for a in affs if a]
        regions = [r for r in regions if r != "Others"]
        if not regions:
            author_region[author] = "Others"
        else:
            author_region[author] = Counter(regions).most_common(1)[0][0]
    return author_region


def plot_from_gml(out_dir: Path, title: str = "Co-authorship network") -> None:
    """Load GML from out_dir and save coauthor_network.png. Use when main run hit plot errors."""
    import networkx as nx
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    gml_path = out_dir / "coauthor_network.gml"
    if not gml_path.exists():
        raise FileNotFoundError(f"No {gml_path}")
    G = nx.read_gml(str(gml_path))
    if G.number_of_edges() == 0:
        print("No edges to plot.")
        return
    pos = nx.spring_layout(G, k=1.5, iterations=80, seed=42, weight="weight")
    deg = dict(G.degree(weight="weight"))
    node_sizes = [max(120, deg.get(n, 1) * 22) for n in G.nodes()]
    weights = [G.edges[u, v].get("weight", 1) for u, v in G.edges()]
    w_min, w_max = min(weights), max(weights)
    edge_widths = [0.8 + 4 * (w - w_min) / (w_max - w_min) if w_max > w_min else 1.2 for w in weights]
    node_colors = [_region_color(G.nodes[n].get("region", "Others")) for n in G.nodes()]
    plt.figure(figsize=(24, 18))
    nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.6, edge_color="gray")
    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors, alpha=0.9, edgecolors="white", linewidths=1)
    labels = {n: n if len(n) <= 25 else n[:22] + "…" for n in G.nodes()}
    nx.draw_networkx_labels(G, pos, labels=labels, font_size=9)
    from matplotlib.patches import Patch
    legend_handles = [Patch(facecolor=_region_color(r), label=r) for r in _REGION_ORDER]
    plt.legend(handles=legend_handles, loc="upper left", fontsize=10)
    plt.title(title, fontsize=14)
    plt.axis("off")
    plt.tight_layout()
    plot_path = out_dir / "coauthor_network.png"
    plt.savefig(plot_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Saved {plot_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Build co-authorship network from papers on a Google Scholar profile. Input: profile URL or user ID."
    )
    parser.add_argument(
        "scholar_url_or_id",
        nargs="?",
        help="Google Scholar profile URL (e.g. https://scholar.google.com/citations?user=JhtCUvIAAAAJ) or user ID (e.g. JhtCUvIAAAAJ)",
    )
    parser.add_argument(
        "--out-dir",
        "-o",
        type=Path,
        default=None,
        help="Output directory (default: outputs/pi_<name_from_profile>)",
    )
    parser.add_argument("--no-plot", action="store_true", help="Skip generating network plot")
    parser.add_argument(
        "--max-authors",
        type=int,
        default=50,
        help="Exclude papers with more than this many authors (default: 50)",
    )
    parser.add_argument(
        "--plot-only",
        type=Path,
        default=None,
        metavar="OUT_DIR",
        help="Only generate figure from existing GML in OUT_DIR",
    )
    args = parser.parse_args()

    if args.plot_only is not None:
        out_dir = Path(args.plot_only)
        if not out_dir.is_dir():
            parser.error(f"--plot-only: not a directory: {out_dir}")
        try:
            plot_from_gml(out_dir, title=f"Co-authorship network: {out_dir.name}")
        except Exception as e:
            print(f"Plot failed: {e}")
            raise SystemExit(1)
        print("Done.")
        return

    raw = (args.scholar_url_or_id or "").strip()
    if not raw:
        parser.error("Google Scholar profile URL or user ID required. Example: https://scholar.google.com/citations?user=JhtCUvIAAAAJ")
    scholar_id = parse_scholar_url_or_id(raw)
    if not scholar_id:
        parser.error("Could not parse Google Scholar user ID from input. Use full URL or the user=... part (e.g. JhtCUvIAAAAJ).")

    print(f"Fetching Google Scholar profile (user={scholar_id})...")
    try:
        scholar_profile_url, pi_name, papers = get_profile_by_scholar_id(scholar_id)
    except Exception as e:
        print(f"Failed to fetch profile: {e}")
        raise SystemExit(1)
    print(f"Profile: {pi_name}. Found {len(papers)} publications (co-authors from Scholar only, no PubMed).")

    if not papers:
        print("No publications on profile. Exiting.")
        return

    out_dir = args.out_dir
    if out_dir is None:
        safe_name = "".join(c if c.isalnum() or c in " -_" else "_" for c in pi_name).strip(" _") or "pi"
        out_dir = Path(__file__).resolve().parent / "outputs" / f"pi_{safe_name}"
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    max_authors = max(1, args.max_authors)
    n_before = len(papers)
    papers = [p for p in papers if len(p.get("authors", [])) <= max_authors]
    n_removed = n_before - len(papers)
    if n_removed:
        print(f"Excluded {n_removed} papers with >{max_authors} authors.")
    n_papers = len(papers)
    print(f"Number of papers: {n_papers}")

    co_edges, authors = build_coauthor_network_from_papers(papers)
    author_region = {n: "Others" for n in authors}
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

    # Save paper list (title, authors from Google Scholar)
    papers_path = out_dir / "papers.txt"
    with open(papers_path, "w") as f:
        f.write("Title\tAuthors\n")
        for p in papers:
            title = (p.get("title") or "").replace("\t", " ").replace("\n", " ")[:200]
            authors_str = "; ".join(p.get("authors", []))
            f.write(f"{title}\t{authors_str}\n")
    print(f"Saved {papers_path}")

    # GML
    try:
        import networkx as nx
        G = nx.Graph()
        for (a, b), w in co_edges.items():
            G.add_edge(a, b, weight=w)
        for n in G.nodes():
            G.nodes[n]["label"] = n
            G.nodes[n]["region"] = author_region.get(n, "Others")
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
                node_sizes = [max(120, deg.get(n, 1) * 22) for n in G.nodes()]
                weights = [G.edges[u, v].get("weight", 1) for u, v in G.edges()]
                w_min, w_max = min(weights), max(weights)
                edge_widths = [0.8 + 4 * (w - w_min) / (w_max - w_min) if w_max > w_min else 1.2 for w in weights]
                node_colors = [_region_color(author_region.get(n, "Others")) for n in G.nodes()]
                plt.figure(figsize=(24, 18))
                nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.6, edge_color="gray")
                nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors, alpha=0.9, edgecolors="white", linewidths=1)
                labels = {n: n if len(n) <= 25 else n[:22] + "…" for n in G.nodes()}
                nx.draw_networkx_labels(G, pos, labels=labels, font_size=9)
                from matplotlib.patches import Patch
                legend_handles = [Patch(facecolor=_region_color(r), label=r) for r in _REGION_ORDER]
                plt.legend(handles=legend_handles, loc="upper left", fontsize=10)
                plt.title(f"Co-authorship network: {pi_name} ({n_papers} papers)", fontsize=14)
                plt.axis("off")
                plt.tight_layout()
                plot_path = out_dir / "coauthor_network.png"
                plt.savefig(plot_path, dpi=200, bbox_inches="tight")
                plt.close()
                print(f"Saved {plot_path}")
        except Exception as e:
            print(f"Plot skipped: {e}")

    if scholar_profile_url:
        print(f"Google Scholar profile (double-check): {scholar_profile_url}")
    print("Done.")


if __name__ == "__main__":
    main()
