# Science network: Metagenomics (PubMed, last 5 years)

This project builds two networks from PubMed metagenomics papers (2020–2025):

1. **Co-corresponding author network**  
   - **Nodes:** Corresponding authors (heuristic: last and second-to-last author per paper).  
   - **Edges:** Two authors are connected if they were co-corresponding on the same paper; edge weight = number of such papers.

2. **University/institution co-affiliation network**  
   - **Nodes:** Universities/institutions parsed from author affiliations.  
   - **Edges:** Two institutions are connected if they co-appeared on the same paper; edge weight = number of such papers.

## Setup

```bash
cd /Users/shuai-wang/methy/sci_net
python -m venv .venv
source .venv/bin/activate   # or .venv\Scripts\activate on Windows
pip install -r requirements.txt
```

## Usage

1. **Fetch data and build networks** (PubMed API; may take several minutes):

   ```bash
   python fetch_pubmed.py
   ```

   This will:
   - Search PubMed for "metagenomics" with publication date 2020–2025.
   - Fetch article XML and parse corresponding authors and affiliations.
   - Write `data/records.json`, `data/author_edges.txt`, and `data/institution_edges.txt`.

2. **Visualize both networks**:

   ```bash
   python visualize_networks.py
   ```

   If `data/` is missing, this script runs `fetch_pubmed.py` first. Images are saved under `outputs/`:
   - `outputs/author_network.png` — co-corresponding author network.
   - `outputs/institution_network.png` — university co-affiliation network.

## Outputs

- `data/records.json` — Per-article corresponding authors and author–affiliation list.
- `data/author_edges.txt` — Author1, Author2, weight (tab-separated).
- `data/institution_edges.txt` — Institution1, Institution2, weight (tab-separated).
- `outputs/author_network.png` — Plot of the author network.
- `outputs/institution_network.png` — Plot of the institution network.

## Notes

- **Corresponding author:** PubMed XML does not always tag corresponding authors. We use the last author (and second-to-last as co-corresponding) as a proxy.
- **Rate limits:** Without an NCBI API key, use at most ~3 requests per second. The script uses small batches and sleeps between calls.
- **Institutions:** Parsed from free-text affiliations with simple patterns (e.g. "University of X", "X University", "X Institute"). You can extend `_affiliation_to_institutions` in `fetch_pubmed.py` for better extraction.
