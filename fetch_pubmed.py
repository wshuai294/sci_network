"""
Fetch metagenomics papers from PubMed and extract
corresponding authors and affiliations for building science networks.
Date range and result limit are configurable (default: last year, no limit).
"""
from __future__ import annotations

import re
import time
from collections import defaultdict
from pathlib import Path

from Bio import Entrez

# NCBI requires an email for API usage
Entrez.email = "sci_net_user@example.com"
Entrez.api_key = None  # set to your API key for higher rate limit

# Batch sizes for API calls (stay under rate limits)
SEARCH_BATCH = 10000
FETCH_BATCH = 200

# Canonical institution name: (canonical, [aliases]) — first matching alias wins (list longer first).
# UC campuses listed explicitly so they stay distinct (not merged into one "University of California").
INSTITUTION_CANONICAL = [
    # University of California campuses (specific first)
    ("UC Berkeley", ["University of California Berkeley", "UC Berkeley", "UCB", "Berkeley"]),
    ("UCLA", ["University of California Los Angeles", "UCLA", "UC Los Angeles"]),
    ("UC San Diego", ["University of California San Diego", "UCSD", "UC San Diego"]),
    ("UC Davis", ["University of California Davis", "UC Davis", "UCD"]),
    ("UC San Francisco", ["University of California San Francisco", "UCSF", "UC San Francisco"]),
    ("UC Irvine", ["University of California Irvine", "UC Irvine", "UCI"]),
    ("UC Santa Barbara", ["University of California Santa Barbara", "UC Santa Barbara", "UCSB"]),
    ("UC Santa Cruz", ["University of California Santa Cruz", "UC Santa Cruz", "UCSC"]),
    ("UC Riverside", ["University of California Riverside", "UC Riverside", "UCR"]),
    ("UC Merced", ["University of California Merced", "UC Merced"]),
    # Other universities (merge common aliases)
    ("Harvard University", ["Harvard Medical School", "Harvard University", "Harvard"]),
    ("Stanford University", ["Stanford University", "Stanford"]),
    ("ETH Zurich", ["ETH Zurich", "ETH Zürich", "Swiss Federal Institute of Technology"]),
]


def _normalize_institution(name: str) -> str:
    """
    Merge equivalent affiliations: strip 'The ', normalize abbrevs (Univ. -> University),
    and map known aliases to a canonical name (e.g. DOE JGI, JGI -> JGI).
    """
    if not name or not name.strip():
        return ""
    s = re.sub(r"\s+", " ", name.strip())
    if s.lower().startswith("the "):
        s = s[4:].strip()
    s_lower = s.lower()
    for canonical, aliases in INSTITUTION_CANONICAL:
        for alias in aliases:
            al = alias.lower()
            if al in s_lower or s_lower in al or s_lower == al:
                return canonical
    s = re.sub(r"\bUniv\.\b", "University", s, flags=re.IGNORECASE)
    s = re.sub(r"\bInst\.\b", "Institute", s, flags=re.IGNORECASE)
    s = re.sub(r"\bDept\.\b", "Department", s, flags=re.IGNORECASE)
    s = re.sub(r"\bColl\.\b", "College", s, flags=re.IGNORECASE)
    s = re.sub(r"\bLab\.\b", "Laboratory", s, flags=re.IGNORECASE)
    return s


# Abbreviations for brief institution names (stored in network)
_BRIEF_ABBREVS = [
    (r"\bUniversity\b", "U"),
    (r"\bUniversities\b", "Us"),
    (r"\bDepartment\b", "Dept"),
    (r"\bCollege\b", "Coll"),
    (r"\bSchool\b", "Sch"),
    (r"\bNational\b", "Natl"),
    (r"\bInternational\b", "Int"),
    (r"\bScience\b", "Sci"),
    (r"\bSciences\b", "Sci"),
    (r"\bTechnology\b", "Tech"),
    (r"\bMedical\b", "Med"),
    (r"\bMedicine\b", "Med"),
    (r"\bResearch\b", "Res"),
    (r"\bCenter\b", "Ctr"),
    (r"\bCentre\b", "Ctr"),
    (r"\bAcademy\b", "Acad"),
    (r"\bFaculty\b", "Fac"),
    (r"\bDivision\b", "Div"),
]


def _brief_institution(name: str, max_len: int = 35) -> str:
    """Shorten institution name for display/storage (University -> U, etc.)."""
    if not name or not name.strip():
        return ""
    s = name.strip()
    for pattern, repl in _BRIEF_ABBREVS:
        s = re.sub(pattern, repl, s, flags=re.IGNORECASE)
    s = re.sub(r"\s+", " ", s).strip()
    if len(s) > max_len:
        s = s[: max_len - 1].rstrip() + "…"
    return s


# PubMed ESearch returns at most 10,000 UIDs per query. To get all papers we search by year and merge.
PUBMED_ESEARCH_MAX = 10000


def search_pmids_by_query(
    query: str,
    mindate: str = "2020/01/01",
    maxdate: str = "2025/12/31",
    max_results: int | None = 500,
) -> list[str]:
    """
    Search PubMed by free-text query (e.g. research area). Returns PMID list.
    If max_results is None, fetches all by splitting date range into year chunks.
    """
    if not query or not query.strip():
        return []
    query = query.strip()
    all_ids = []

    if max_results is None:
        import datetime
        start = datetime.datetime.strptime(mindate[:10].replace("/", "-"), "%Y-%m-%d")
        end = datetime.datetime.strptime(maxdate[:10].replace("/", "-"), "%Y-%m-%d")
        seen = set()
        y = start.year
        while y <= end.year:
            y_mindate = mindate if y == start.year else f"{y}/01/01"
            y_maxdate = maxdate if y == end.year else f"{y}/12/31"
            retstart = 0
            while True:
                try:
                    handle = Entrez.esearch(
                        db="pubmed",
                        term=query,
                        mindate=y_mindate,
                        maxdate=y_maxdate,
                        datetype="pdat",
                        retmax=min(PUBMED_ESEARCH_MAX, PUBMED_ESEARCH_MAX - retstart),
                        retstart=retstart,
                        sort="relevance",
                    )
                    data = Entrez.read(handle)
                    handle.close()
                except Exception as e:
                    print(f"Search error ({y}): {e}")
                    break
                ids = data.get("IdList", [])
                if not ids:
                    break
                for i in ids:
                    if i not in seen:
                        seen.add(i)
                        all_ids.append(i)
                if len(ids) < PUBMED_ESEARCH_MAX:
                    break
                retstart += len(ids)
                time.sleep(0.34)
            y += 1
            time.sleep(0.34)
        return all_ids

    retstart = 0
    chunk = min(SEARCH_BATCH, max_results or 500)
    limit = max_results if max_results is not None else 500
    while retstart < limit:
        try:
            retmax = min(chunk, limit - retstart)
            handle = Entrez.esearch(
                db="pubmed",
                term=query,
                mindate=mindate,
                maxdate=maxdate,
                datetype="pdat",
                retmax=retmax,
                retstart=retstart,
                sort="relevance",
            )
            data = Entrez.read(handle)
            handle.close()
        except Exception as e:
            print(f"Search error: {e}")
            break
        ids = data.get("IdList", [])
        if not ids:
            break
        all_ids.extend(ids)
        retstart += len(ids)
        if len(ids) < chunk:
            break
        time.sleep(0.34)
    return all_ids[:limit]


def search_metagenomics_pmids(
    mindate: str = "2020/01/01",
    maxdate: str = "2025/12/31",
    max_results: int | None = None,
) -> list[str]:
    """
    Search PubMed for metagenomics papers. If max_results is None, fetch ALL by splitting
    the date range into year chunks (PubMed returns at most 10,000 per query).
    """
    return search_pmids_by_query("metagenomics", mindate=mindate, maxdate=maxdate, max_results=max_results)


def _normalize_author(author_dict: dict) -> str:
    """Build full name from Biopython author dict."""
    last = author_dict.get("LastName", "")
    fore = author_dict.get("ForeName", "")
    if fore:
        return f"{last} {fore}".strip()
    return last.strip() or "Unknown"


def _extract_affiliation(author_dict: dict) -> list[str]:
    """Extract affiliation strings from an author."""
    affs = []
    if "AffiliationInfo" not in author_dict:
        return affs
    for info in author_dict["AffiliationInfo"]:
        if isinstance(info, dict) and "Affiliation" in info:
            affs.append(str(info["Affiliation"]).strip())
    return affs


def _is_corresponding(author_dict: dict, index: int, total: int) -> bool:
    """
    Heuristic: corresponding authors are often the last and second-to-last (co-corresponding).
    PubMed XML does not always tag 'corresp'; we use last author and co-last when present.
    """
    if index == total - 1:
        return True
    if total >= 2 and index == total - 2:
        return True
    return False


def _affiliation_to_institutions(aff_text: str) -> list[str]:
    """
    Extract only university names from affiliation string (ignore institutes;
    many institutes belong to a university and would double-count).
    Patterns: "University of X", "X University". Fallback: segment containing "University".
    """
    if not aff_text or not aff_text.strip():
        return []
    text = aff_text.strip()
    institutions = []
    parts = re.split(r"[,;]\s*", text)

    for part in parts:
        part = part.strip()
        # "University of Something"
        m = re.search(
            r"\b(University\s+of\s+[^,;]+?)(?:\s*,\s*[A-Z]{2}\s*\d|\.|$)",
            part,
            re.IGNORECASE,
        )
        if m:
            institutions.append(m.group(1).strip().rstrip(".,"))
        # "Something University"
        m2 = re.search(
            r"\b([A-Za-z\s\-]+University)(?:\s*,\s*|\s*\.|$)",
            part,
        )
        if m2:
            name = m2.group(1).strip().rstrip(".,")
            if len(name) > 5 and "University" in name:
                institutions.append(name)

    # Fallback: only use segment containing "University" (ignore institute-only segments)
    if not institutions and len(text) > 10:
        for seg in parts:
            if "niversity" in seg.lower():
                clean = re.sub(r"\s*\d{5,}.*$", "", seg).strip().rstrip(".,")
                if len(clean) > 8:
                    institutions.append(clean)
                    break

    return list(dict.fromkeys(institutions))  # dedupe order-preserving


def _affiliation_to_institutions_canonical(aff_text: str) -> list[str]:
    """Extract university names only, merge to canonical, then brief for network building."""
    raw = _affiliation_to_institutions(aff_text)
    out = []
    seen = set()
    for inst in raw:
        canonical = _normalize_institution(inst)
        if not canonical:
            continue
        brief = _brief_institution(canonical)
        if brief and brief not in seen:
            seen.add(brief)
            out.append(brief)
    return out


def fetch_articles(pmids: list[str]) -> list[dict]:
    """Fetch full article records from PubMed and return list of parsed records."""
    records = []
    for i in range(0, len(pmids), FETCH_BATCH):
        batch = pmids[i : i + FETCH_BATCH]
        try:
            handle = Entrez.efetch(
                db="pubmed",
                id=batch,
                retmode="xml",
            )
            data = Entrez.read(handle)
            handle.close()
        except Exception as e:
            print(f"Fetch error for batch: {e}")
            time.sleep(1)
            continue

        for article in data.get("PubmedArticle", []):
            try:
                rec = _parse_article(article)
                if rec:
                    records.append(rec)
            except Exception as e:
                print(f"Parse error: {e}")
                continue

        time.sleep(0.34)

    return records


def _parse_article(article: dict) -> dict | None:
    """Parse one PubmedArticle into corresponding authors and affiliations."""
    try:
        med = article["MedlineCitation"]
        art = med["Article"]
        raw_authors = art.get("AuthorList", [])
    except KeyError:
        return None

    # Normalize: Entrez can return a single Author dict when there is one author (not a list).
    if not raw_authors:
        return None
    if isinstance(raw_authors, dict) and "LastName" in raw_authors:
        author_list = [raw_authors]
    elif isinstance(raw_authors, list):
        author_list = [a for a in raw_authors if isinstance(a, dict) and a.get("LastName")]
    else:
        return None

    if not author_list:
        return None

    title = art.get("ArticleTitle", "")
    if isinstance(title, list):
        title = " ".join(str(t) for t in title) if title else ""
    else:
        title = str(title or "").strip()

    journal = ""
    j = art.get("Journal", {}) or {}
    if isinstance(j, dict):
        journal = j.get("Title") or j.get("ISOAbbreviation") or ""
    if isinstance(journal, list):
        journal = " ".join(str(x) for x in journal) if journal else ""
    journal = str(journal or "").strip()

    total = len(author_list)
    corresponding = []
    # All authors and their affiliations for university network
    author_affiliations = []

    for idx, author in enumerate(author_list):
        if author.get("ValidYN") == "N":
            continue
        name = _normalize_author(author)
        affs = _extract_affiliation(author)
        if _is_corresponding(author, idx, total):
            corresponding.append((name, affs))
        author_affiliations.append((name, affs))

    if not corresponding:
        # Fallback: treat last author as corresponding
        last_author = author_list[-1]
        if last_author.get("ValidYN") != "N":
            name = _normalize_author(last_author)
            affs = _extract_affiliation(last_author)
            corresponding.append((name, affs))

    pmid = str(med["PMID"])
    return {
        "pmid": pmid,
        "title": title,
        "journal": journal,
        "corresponding": corresponding,
        "author_affiliations": author_affiliations,
    }


# Region keywords for inferring author/institution region (match pi_coauthor_network and web legend)
_REGION_KEYWORDS = {
    "China": ["China", "Chinese", "Beijing", "Shanghai", "Hong Kong", "Fudan", "Tsinghua", "Zhejiang", "Wuhan", "Nanjing", "Sun Yat-sen", "Chinese Acad", "CAS", "Peking", "Huazhong", "Nankai"],
    "US": ["USA", "U.S.A", "United States", "America", "Harvard", "Stanford", "MIT", "Yale", "Columbia", "California", "Berkeley", "Michigan", "Cornell", "Duke", "Johns Hopkins", "Texas", "Washington", "Chicago", "Boston", "U.S."],
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


def institution_region(inst_name: str) -> str:
    """Infer region from institution name (for network node coloring)."""
    return _affiliation_to_region(inst_name or "")


def build_author_regions(records: list[dict]) -> dict[str, str]:
    """Build author -> region from all author affiliations (majority vote per author)."""
    author_affs: dict[str, list[str]] = defaultdict(list)
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
            from collections import Counter
            author_region[author] = Counter(regions).most_common(1)[0][0]
    return author_region


def build_author_network(records: list[dict]) -> tuple[defaultdict, set]:
    """
    Build co-authorship network (all authors per paper).
    Node = author; Edge weight = number of papers where they co-authored.
    Using all authors avoids underestimating weights (co-corresponding only used last 2).
    """
    co_edges = defaultdict(int)
    all_authors = set()

    for rec in records:
        # Use all authors on the paper (co-authorship), not just corresponding
        names = list(dict.fromkeys(c[0] for c in rec.get("author_affiliations", []) if c[0]))
        for n in names:
            all_authors.add(n)
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                a, b = sorted([names[i], names[j]])
                co_edges[(a, b)] += 1

    return co_edges, all_authors


def build_corresponding_author_network(records: list[dict]) -> tuple[defaultdict, set]:
    """
    Build co-authorship network using only corresponding authors per paper (PI network).
    A paper can have more than one corresponding author (e.g. last + co-last).
    Node = corresponding author; edge weight = number of papers where they were co-corresponding.
    Papers with 0 or 1 corresponding author contribute no edges.
    """
    co_edges = defaultdict(int)
    all_authors = set()

    for rec in records:
        corresponding = rec.get("corresponding") or []
        names = list(dict.fromkeys(c[0] for c in corresponding if c[0]))
        for n in names:
            all_authors.add(n)
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                a, b = sorted([names[i], names[j]])
                co_edges[(a, b)] += 1

    return co_edges, all_authors


# Path to QS top 500 list (one institution name per line). If missing, embedded list is used.
_DATA_DIR = Path(__file__).resolve().parent / "data"
_QS_TOP500_PATH = _DATA_DIR / "qs_top500.txt"

# Embedded QS-style names (top ~500 or so) when data/qs_top500.txt is not present.
_QS_TOP500_EMBEDDED = (
    "Massachusetts Institute of Technology", "University of Cambridge", "University of Oxford",
    "Harvard University", "Stanford University", "Imperial College London",
    "ETH Zurich", "National University of Singapore", "UCL", "UC Berkeley",
    "University of Chicago", "University of Pennsylvania", "Cornell University",
    "University of Melbourne", "California Institute of Technology", "Yale University",
    "Peking University", "Tsinghua University", "Princeton University", "University of New South Wales",
    "University of Sydney", "University of Toronto", "University of Edinburgh", "Columbia University",
    "Université PSL", "Nanyang Technological University", "Johns Hopkins University", "University of Hong Kong",
    "University of Tokyo", "McGill University", "University of Manchester", "University of Michigan",
    "Australian National University", "Northwestern University", "Fudan University", "University of Bristol",
    "KU Leuven", "Technical University of Munich", "London School of Economics", "Delft University of Technology",
    "Monash University", "University of Texas at Austin", "University of Illinois Urbana-Champaign",
    "University of Amsterdam", "Hong Kong University of Science and Technology", "King's College London",
    "Kyoto University", "Seoul National University", "KAIST", "Shanghai Jiao Tong University",
    "University of Wisconsin-Madison", "University of Washington", "Lund University", "KTH Royal Institute of Technology",
    "University of British Columbia", "Institut Polytechnique de Paris", "New York University",
    "Chinese University of Hong Kong", "Zhejiang University", "University of Queensland",
    "Duke University", "University of North Carolina Chapel Hill", "Osaka University",
    "Boston University", "University of Zurich", "Sorbonne University", "Brown University",
    "Penn State University", "University of Leeds", "University of Birmingham", "University of Sheffield",
    "Rice University", "Ohio State University", "University of Western Australia", "University of Southampton",
    "University of Groningen", "University of Oslo", "University of Helsinki", "University of Geneva",
    "University of Warwick", "University of St Andrews", "University of Nottingham", "University of York",
    "University of Glasgow", "University of Birmingham", "University of Liverpool", "Durham University",
    "University of Copenhagen", "University of Barcelona", "Universidad de Buenos Aires",
    "University of California Los Angeles", "University of California San Diego", "University of California Davis",
    "University of California Santa Barbara", "University of Maryland", "Purdue University",
    "University of Pittsburgh", "University of Minnesota", "University of Florida", "Texas A&M University",
    "University of Colorado Boulder", "University of Arizona", "University of Utah", "Michigan State University",
    "Georgia Institute of Technology", "University of Virginia", "University of Iowa", "Indiana University",
    "Vanderbilt University", "Emory University", "Case Western Reserve University", "University of Rochester",
    "Carnegie Mellon University", "University of Southern California", "University of California Irvine",
    "University of California San Francisco", "University of California Santa Cruz", "University of California Riverside",
    "Arizona State University", "Florida State University", "North Carolina State University", "University of Massachusetts",
    "University of Waterloo", "McMaster University", "University of Alberta", "Western University",
    "Simon Fraser University", "University of Calgary", "Queen's University", "Université de Montréal",
    "University of Geneva", "EPFL", "University of Lausanne", "University of Bern", "University of Basel",
    "Heidelberg University", "Ludwig Maximilian University of Munich", "Humboldt University of Berlin",
    "Free University of Berlin", "University of Freiburg", "University of Bonn", "University of Göttingen",
    "University of Hamburg", "University of Cologne", "University of Frankfurt", "University of Stuttgart",
    "University of Tübingen", "University of Würzburg", "University of Erlangen-Nuremberg", "RWTH Aachen",
    "University of Paris", "University of Lyon", "Aix-Marseille University", "University of Strasbourg",
    "University of Montpellier", "University of Lille", "University of Bordeaux", "University of Toulouse",
    "University of Amsterdam", "Utrecht University", "Erasmus University Rotterdam", "Leiden University",
    "Vrije Universiteit Amsterdam", "Wageningen University", "University of Groningen", "Eindhoven University of Technology",
    "KU Leuven", "Ghent University", "University of Brussels", "University of Leuven",
    "University of Barcelona", "Autonomous University of Barcelona", "Complutense University of Madrid",
    "University of Valencia", "University of Granada", "University of Seville", "University of Zaragoza",
    "University of Salamanca", "University of the Basque Country", "Polytechnic University of Catalonia",
    "Sapienza University of Rome", "University of Bologna", "University of Milan", "University of Padua",
    "University of Pisa", "University of Naples Federico II", "University of Turin", "Politecnico di Milano",
    "Stockholm University", "Uppsala University", "University of Gothenburg", "Chalmers University of Technology",
    "University of Oslo", "Norwegian University of Science and Technology", "University of Copenhagen",
    "Aarhus University", "University of Helsinki", "Aalto University", "University of Turku",
    "University of Vienna", "Vienna University of Technology", "University of Zurich", "ETH Zurich",
    "University of Bern", "University of Geneva", "University of Lausanne", "University of Basel",
    "University of Warsaw", "Jagiellonian University", "Warsaw University of Technology", "Charles University Prague",
    "Masaryk University", "Eötvös Loránd University", "University of Belgrade", "University of Athens",
    "National Technical University of Athens", "Aristotle University of Thessaloniki", "University of Lisbon",
    "University of Porto", "University of Coimbra", "Trinity College Dublin", "University College Dublin",
    "University of Edinburgh", "University of Manchester", "University of Bristol", "University of Glasgow",
    "University of Birmingham", "University of Southampton", "University of Leeds", "University of Sheffield",
    "University of Nottingham", "University of St Andrews", "University of Warwick", "University of York",
    "University of Liverpool", "Newcastle University", "Queen Mary University of London", "University of Exeter",
    "University of Sussex", "University of Leicester", "University of East Anglia", "University of Bath",
    "University of Surrey", "Lancaster University", "University of Reading", "University of Aberdeen",
    "University of Dublin", "University of Helsinki", "University of Tokyo", "Kyoto University",
    "Tokyo Institute of Technology", "Osaka University", "Tohoku University", "Nagoya University",
    "Hokkaido University", "Kyushu University", "Waseda University", "Keio University", "University of Tsukuba",
    "Seoul National University", "KAIST", "Yonsei University", "Korea University", "Pohang University of Science and Technology",
    "Sungkyunkwan University", "Hanyang University", "Ewha Womans University", "Sogang University",
    "National Taiwan University", "National Tsing Hua University", "National Yang Ming Chiao Tung University",
    "National Sun Yat-Sen University", "Hong Kong University", "Chinese University of Hong Kong",
    "Hong Kong University of Science and Technology", "Hong Kong Polytechnic University", "City University of Hong Kong",
    "Peking University", "Tsinghua University", "Fudan University", "Zhejiang University", "Shanghai Jiao Tong University",
    "University of Science and Technology of China", "Nanjing University", "Wuhan University", "Harbin Institute of Technology",
    "Sun Yat-sen University", "Sichuan University", "Tongji University", "Beijing Normal University",
    "Nankai University", "Tianjin University", "Xiamen University", "Jilin University", "East China Normal University",
    "Shanghai University", "Huazhong University of Science and Technology", "Dalian University of Technology",
    "National University of Singapore", "Nanyang Technological University", "Singapore Management University",
    "University of Malaya", "Universiti Kebangsaan Malaysia", "Universiti Putra Malaysia", "Chulalongkorn University",
    "Mahidol University", "Bandung Institute of Technology", "University of Indonesia", "Gadjah Mada University",
    "University of the Philippines", "Indian Institute of Technology Bombay", "Indian Institute of Technology Delhi",
    "Indian Institute of Science", "Indian Institute of Technology Madras", "University of Delhi",
    "Jawaharlal Nehru University", "University of Calcutta", "University of Hyderabad", "Banaras Hindu University",
    "University of Sydney", "University of Melbourne", "University of Queensland", "Monash University",
    "University of New South Wales", "University of Western Australia", "University of Adelaide",
    "Australian National University", "University of Auckland", "University of Otago", "Victoria University of Wellington",
    "University of Cape Town", "University of Witwatersrand", "Stellenbosch University", "University of Pretoria",
    "Cairo University", "American University in Cairo", "Tel Aviv University", "Hebrew University of Jerusalem",
    "Technion", "Bar-Ilan University", "Ben-Gurion University", "King Saud University", "King Abdullah University of Science and Technology",
    "American University of Beirut", "Bogazici University", "Bilkent University", "Sabanci University",
    "Koç University", "Middle East Technical University", "University of Tehran", "Sharif University of Technology",
    "Pontifical Catholic University of Chile", "University of Chile", "University of São Paulo", "State University of Campinas",
    "Federal University of Rio de Janeiro", "Universidad de los Andes", "National University of Colombia",
    "University of the Andes", "Tecnológico de Monterrey", "National Autonomous University of Mexico",
    "University of Georgia", "University of Delaware", "University of Oregon", "Oregon State University",
    "University of Iowa", "Iowa State University", "University of Connecticut", "University of Maryland Baltimore",
    "Rutgers University", "Stony Brook University", "University at Buffalo", "University of Illinois Chicago",
    "University of Kansas", "University of Kentucky", "University of Nebraska Lincoln", "University of Oklahoma",
    "University of South Carolina", "University of Tennessee", "Virginia Tech", "Washington State University",
    "University of Alberta", "Dalhousie University", "Université Laval", "University of Ottawa",
    "University of Victoria", "York University", "Concordia University", "Université du Québec",
    "University of Leeds", "University of Ulster", "University of Portsmouth", "University of Hull",
    "University of Northumbria", "University of Bayreuth", "University of Bremen", "University of Kiel",
    "University of Ulm", "University of Genoa", "University of Crete", "University of Costa Rica",
    "University of Madras", "ITMO University", "Peter the Great St Petersburg Polytechnic University",
    "Moscow State Institute of International Relations", "Lahore University of Management Sciences",
    "Edith Cowan University", "University of Palermo", "University of Alcalá", "University of Utara Malaysia",
    "National University de La Plata", "City University of New York", "Macau University of Science and Technology",
    "Catholic University of the Sacred Heart", "Kyungpook National University", "University of Texas at Dallas",
    "China Agricultural University", "Universidad de Palermo", "Pontifical Catholic University Argentina",
    "University of Technology Brunei", "School of Oriental and African Studies",
)


def _get_qs_top500_brief_set() -> set:
    """Return set of brief institution names that are in QS top 500 (for filtering affiliation network)."""
    brief_set = set()
    if _QS_TOP500_PATH.exists():
        try:
            with open(_QS_TOP500_PATH, encoding="utf-8") as f:
                for line in f:
                    name = line.strip()
                    if name and not name.startswith("#"):
                        b = _brief_institution(_normalize_institution(name))
                        if b:
                            brief_set.add(b)
        except Exception:
            pass
    if not brief_set:
        for name in _QS_TOP500_EMBEDDED:
            b = _brief_institution(_normalize_institution(name))
            if b:
                brief_set.add(b)
    return brief_set


def filter_institution_network_to_qs_top500(
    inst_edges: defaultdict,
    institutions: set,
) -> tuple[defaultdict, set]:
    """
    Keep only institutions that match QS World University Rankings top 500.
    Returns (filtered_inst_edges, filtered_institutions).
    """
    allowed = _get_qs_top500_brief_set()
    keep = institutions & allowed
    if not keep:
        return inst_edges, institutions  # no match: return as-is to avoid empty network
    filtered = defaultdict(int)
    for (a, b), w in inst_edges.items():
        if a in keep and b in keep:
            filtered[(a, b)] = w
    return filtered, keep


def build_institution_network(records: list[dict]) -> tuple[defaultdict, set]:
    """
    Build institution co-affiliation network.
    Node = institution (university etc.); Edge weight = number of papers where both appear.
    """
    co_edges = defaultdict(int)
    all_inst = set()

    for rec in records:
        # Per paper: collect all institutions (canonical/merged names only)
        paper_inst = set()
        for _name, affs in rec.get("author_affiliations", []):
            for aff_text in affs:
                for inst in _affiliation_to_institutions_canonical(aff_text):
                    if inst:
                        paper_inst.add(inst)
                        all_inst.add(inst)
        paper_list = sorted(paper_inst)
        for i in range(len(paper_list)):
            for j in range(i + 1, len(paper_list)):
                a, b = paper_list[i], paper_list[j]
                co_edges[(a, b)] += 1

    return co_edges, all_inst


def _write_papers_multiple_corresponding_or_affiliations(records: list[dict], out_dir: Path) -> None:
    """
    Write to file all papers that have >1 corresponding author OR >1 distinct university
    (institutes ignored in counting). Output: paper title, PMID, corresponding author names,
    and university affiliation names (brief, canonical).
    """
    out_path = out_dir / "papers_multiple_corresponding_or_affiliations.tsv"
    lines = []
    lines.append("PMID\tTitle\tCorrespondingAuthors\tAffiliationNames")
    count_included = 0
    for rec in records:
        corresponding = rec.get("corresponding", [])
        n_corr = len(corresponding)
        # Distinct universities on this paper (ignore institutes; use same canonical/brief as network)
        paper_unis = set()
        for _name, affs in rec.get("author_affiliations", []):
            for aff_text in affs:
                for inst in _affiliation_to_institutions_canonical(aff_text):
                    if inst:
                        paper_unis.add(inst)
        n_aff = len(paper_unis)
        if n_corr <= 1 and n_aff <= 1:
            continue
        count_included += 1
        title = (rec.get("title") or "").replace("\t", " ").replace("\n", " ")
        corr_names = "; ".join(c[0] for c in corresponding if c[0])
        aff_names = "; ".join(sorted(paper_unis))
        lines.append(f"{rec['pmid']}\t{title}\t{corr_names}\t{aff_names}")
    with open(out_path, "w") as f:
        f.write("\n".join(lines))
    print(f"Papers with >1 corresponding author or >1 university: {count_included} -> {out_path}")


def main(
    mindate: str = "2020/01/01",
    maxdate: str = "2025/12/31",
    max_results: int | None = None,
):
    import json

    out_dir = Path(__file__).resolve().parent / "data"
    out_dir.mkdir(exist_ok=True)

    limit_str = "all (by year chunks to exceed 10k limit)" if max_results is None else f"max {max_results}"
    print(f"Searching PubMed for metagenomics papers ({mindate} to {maxdate}, {limit_str})...")
    pmids = search_metagenomics_pmids(mindate=mindate, maxdate=maxdate, max_results=max_results)
    print(f"Found {len(pmids)} PMIDs. Fetching abstracts...")

    records = fetch_articles(pmids)
    n_papers = len(records)
    print(f"Number of papers: {n_papers}")
    print(f"Parsed {n_papers} articles.")

    with open(out_dir / "records.json", "w") as f:
        # Serialize for records (convert defaultdicts if any in nested structure)
        out = []
        for r in records:
            out.append(
                {
                    "pmid": r["pmid"],
                    "title": r.get("title", ""),
                    "corresponding": r["corresponding"],
                    "author_affiliations": r["author_affiliations"],
                }
            )
        json.dump(out, f, indent=0)

    co_edges, authors = build_author_network(records)
    inst_edges, insts = build_institution_network(records)
    author_region = build_author_regions(records)

    # Save network edge lists
    with open(out_dir / "author_edges.txt", "w") as f:
        f.write("Author1\tAuthor2\tWeight\n")
        for (a, b), w in sorted(co_edges.items(), key=lambda x: -x[1]):
            f.write(f"{a}\t{b}\t{w}\n")

    with open(out_dir / "author_regions.txt", "w") as f:
        f.write("Author\tRegion\n")
        for a in sorted(author_region.keys()):
            f.write(f"{a}\t{author_region[a]}\n")

    with open(out_dir / "institution_edges.txt", "w") as f:
        f.write("Institution1\tInstitution2\tWeight\n")
        for (a, b), w in sorted(inst_edges.items(), key=lambda x: -x[1]):
            f.write(f"{a}\t{b}\t{w}\n")

    print(f"Author network: {len(authors)} nodes, {len(co_edges)} edges")
    print(f"Institution network: {len(insts)} nodes, {len(inst_edges)} edges")

    # Output papers with >1 corresponding author or >1 university (ignore institutes)
    _write_papers_multiple_corresponding_or_affiliations(records, out_dir)

    print("Data saved under data/")
    return records, co_edges, authors, inst_edges, insts


def rebuild_networks_from_records(data_dir: Path | None = None) -> None:
    """Rebuild author/institution edges and author_regions from data/records.json (no API calls)."""
    import json

    data_dir = data_dir or Path(__file__).resolve().parent / "data"
    records_path = data_dir / "records.json"
    if not records_path.exists():
        print("No data/records.json found. Run fetch first (no --rebuild-only).")
        return
    with open(records_path) as f:
        records = json.load(f)
    print(f"Loaded {len(records)} records from data/records.json. Rebuilding networks with merged institutions...")

    co_edges, authors = build_author_network(records)
    inst_edges, insts = build_institution_network(records)
    author_region = build_author_regions(records)

    with open(data_dir / "author_edges.txt", "w") as f:
        f.write("Author1\tAuthor2\tWeight\n")
        for (a, b), w in sorted(co_edges.items(), key=lambda x: -x[1]):
            f.write(f"{a}\t{b}\t{w}\n")
    with open(data_dir / "author_regions.txt", "w") as f:
        f.write("Author\tRegion\n")
        for a in sorted(author_region.keys()):
            f.write(f"{a}\t{author_region[a]}\n")
    with open(data_dir / "institution_edges.txt", "w") as f:
        f.write("Institution1\tInstitution2\tWeight\n")
        for (a, b), w in sorted(inst_edges.items(), key=lambda x: -x[1]):
            f.write(f"{a}\t{b}\t{w}\n")

    _write_papers_multiple_corresponding_or_affiliations(records, data_dir)

    print(f"Author network: {len(authors)} nodes, {len(co_edges)} edges")
    print(f"Institution network: {len(insts)} nodes, {len(inst_edges)} edges (merged/ canonical names)")
    print("Data saved under data/")


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] in ("--rebuild-only", "-r"):
        rebuild_networks_from_records()
    else:
        mindate = sys.argv[1] if len(sys.argv) > 1 else "2020/01/01"
        maxdate = sys.argv[2] if len(sys.argv) > 2 else "2025/12/31"
        max_results = int(sys.argv[3]) if len(sys.argv) > 3 else None
        main(mindate=mindate, maxdate=maxdate, max_results=max_results)
