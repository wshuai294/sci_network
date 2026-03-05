#!/usr/bin/env python3
"""Run fetch then visualization to produce both network figures."""
from fetch_pubmed import main as fetch_main
from visualize_networks import main as viz_main

if __name__ == "__main__":
    print("Step 1: Fetch PubMed metagenomics data and build networks...")
    fetch_main()
    print("Step 2: Visualize networks...")
    viz_main()
    print("Done. See outputs/ for author_network.png and institution_network.png.")
