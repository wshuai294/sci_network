#!/usr/bin/env python3
"""Standalone: read coauthor_network.gml from an output dir and save coauthor_network.png."""
from pathlib import Path
import argparse

_REGION_ORDER = ("China", "US", "UK", "Europe", "Japan", "Australia", "Canada", "Others")
_REGION_COLORS = {
    "China": "#e74c3c",
    "US": "#3498db",
    "UK": "#2ecc71",
    "Europe": "#9b59b6",
    "Japan": "#f39c12",
    "Australia": "#1abc9c",
    "Canada": "#e67e22",
    "Others": "#95a5a6",
}


def main():
    parser = argparse.ArgumentParser(description="Plot co-author network from existing GML.")
    parser.add_argument("out_dir", type=Path, help="Output directory containing coauthor_network.gml")
    args = parser.parse_args()
    out_dir = args.out_dir.resolve()
    gml_path = out_dir / "coauthor_network.gml"
    if not gml_path.exists():
        raise SystemExit(f"Not found: {gml_path}")

    import networkx as nx
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

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
    node_colors = [_REGION_COLORS.get(G.nodes[n].get("region", "Others"), _REGION_COLORS["Others"]) for n in G.nodes()]
    plt.figure(figsize=(24, 18))
    nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.6, edge_color="gray")
    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors, alpha=0.9, edgecolors="white", linewidths=1)
    labels = {n: n if len(n) <= 25 else n[:22] + "…" for n in G.nodes()}
    nx.draw_networkx_labels(G, pos, labels=labels, font_size=9)
    from matplotlib.patches import Patch
    legend_handles = [Patch(facecolor=_REGION_COLORS[r], label=r) for r in _REGION_ORDER]
    plt.legend(handles=legend_handles, loc="upper left", fontsize=10)
    plt.title(f"Co-authorship network: {out_dir.name}", fontsize=14)
    plt.axis("off")
    plt.tight_layout()
    plot_path = out_dir / "coauthor_network.png"
    plt.savefig(plot_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Saved {plot_path}")


if __name__ == "__main__":
    main()
