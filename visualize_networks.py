"""
Visualize co-corresponding author network and institution co-affiliation network.
Shows only top 200 nodes by degree per network and writes corresponding edge files.
Institution network: abbreviated labels, region colors, overlap-resistant labels.
"""
import re
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

TOP_NODES = 200

# Abbreviations for institution labels (word-boundary replacement)
ABBREVS = [
    (r"\bUniversity\b", "U"),
    (r"\bUniversities\b", "Us"),
    (r"\bInstitute\b", "Inst"),
    (r"\bInstitutes\b", "Insts"),
    (r"\bDepartment\b", "Dept"),
    (r"\bCollege\b", "Coll"),
    (r"\bSchool\b", "Sch"),
    (r"\bLaboratory\b", "Lab"),
    (r"\bLaboratories\b", "Labs"),
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
    (r"\bKey Laboratory\b", "Key Lab"),
]

# Region keywords (order matters: more specific first)
REGION_KEYWORDS = {
    "China": [
        "China", "Chinese", "Beijing", "Shanghai", "Hong Kong", "Fudan", "Tsinghua",
        "Zhejiang", "Wuhan", "Nanjing", "Sun Yat-sen", "Harbin", "Xiamen", "Sichuan",
        "Tianjin", "Shandong", "Hunan", "Zhongshan", "Guangzhou", "Shenzhen", "Hainan",
        "Chinese Acad", "CAS", "Chinese Academy", "Peking", "Huazhong", "Nankai",
    ],
    "US": [
        "USA", "U.S.A", "United States", "America", "Harvard", "Stanford", "MIT",
        "Yale", "Columbia", "California", "Berkeley", "Michigan", "Cornell", "Duke",
        "Johns Hopkins", "Pennsylvania", "North Carolina", "Texas", "Washington",
        "Chicago", "Boston", "San Francisco", "USA.", "U.S.",
    ],
    "UK": [
        "UK", "U.K", "United Kingdom", "England", "Scotland", "Wales", "Oxford",
        "Cambridge", "London", "Edinburgh", "Manchester", "Birmingham", "Imperial",
        "King's College London", "UCL", "Queen Mary", "Bristol", "Leeds", "Glasgow",
    ],
    "Japan": [
        "Japan", "Japanese", "Tokyo", "Kyoto", "Osaka", "Waseda", "Tohoku", "Nagoya",
        "Hokkaido", "Kyushu", "Keio", "Riken",
    ],
    "Europe": [
        "Germany", "German", "France", "French", "Netherlands", "Dutch", "Switzerland",
        "Swedish", "Sweden", "Denmark", "Danish", "Italy", "Italian", "Spain", "Spanish",
        "Austria", "Belgium", "Norway", "Finland", "Ireland", "Portugal", "Poland",
        "Munich", "Berlin", "Paris", "Amsterdam", "Zurich", "Copenhagen", "Madrid",
        "Leiden", "Utrecht", "Groningen", "ETH", "EMBL", "Max Planck", "CNRS", "INSERM",
        "Pasteur", "Karolinska", "Vienna", "Helsinki",
    ],
}


def abbreviate_institution(name: str, max_len: int = 28) -> str:
    """Shorten institution name with abbreviations (University→U, etc.) and optional truncation."""
    s = name
    for pattern, repl in ABBREVS:
        s = re.sub(pattern, repl, s, flags=re.IGNORECASE)
    s = re.sub(r"\s+", " ", s).strip()
    if len(s) > max_len:
        s = s[: max_len - 1].rstrip() + "…"
    return s


def get_region(inst_name: str) -> str:
    """Infer region from institution name (China, US, UK, Europe, Japan, Others)."""
    name_lower = inst_name.lower()
    for region, keywords in REGION_KEYWORDS.items():
        for kw in keywords:
            if kw.lower() in name_lower:
                return region
    return "Others"


def get_region_color(region: str) -> str:
    """Distinct colors per region for institution nodes."""
    palette = {
        "China": "#e74c3c",
        "US": "#3498db",
        "UK": "#2ecc71",
        "Europe": "#9b59b6",
        "Japan": "#f39c12",
        "Others": "#95a5a6",
    }
    return palette.get(region, palette["Others"])


def load_author_network(data_dir: Path):
    edges_path = data_dir / "author_edges.txt"
    if not edges_path.exists():
        return None
    G = nx.Graph()
    with open(edges_path) as f:
        next(f)  # header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                a, b, w = parts[0], parts[1], int(parts[2])
                G.add_edge(a, b, weight=w)
    return G


def load_institution_network(data_dir: Path):
    edges_path = data_dir / "institution_edges.txt"
    if not edges_path.exists():
        return None
    G = nx.Graph()
    with open(edges_path) as f:
        next(f)
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                a, b, w = parts[0], parts[1], int(parts[2])
                G.add_edge(a, b, weight=w)
    return G


def top_k_subgraph(G: nx.Graph, k: int = TOP_NODES) -> nx.Graph:
    """Return subgraph induced by the top k nodes by weighted degree."""
    if G is None or G.number_of_nodes() == 0:
        return G
    deg = G.degree(weight="weight")
    top_nodes = sorted(deg, key=lambda x: x[1], reverse=True)[:k]
    top_names = [n for n, _ in top_nodes]
    return G.subgraph(top_names).copy()


def _write_gml(G: nx.Graph, path: Path, node_region: dict, node_label=None):
    """Export graph to GML with node attributes region and label. node_label(n) optional."""
    if G is None:
        return
    H = G.copy()
    for n in H.nodes():
        H.nodes[n]["region"] = node_region.get(n, "Others")
        H.nodes[n]["label"] = (node_label(n) if node_label else n)
    for u, v in H.edges():
        if "weight" not in H.edges[u, v]:
            H.edges[u, v]["weight"] = 1
    try:
        nx.write_gml(H, path)
        print(f"Saved {path}")
    except Exception as e:
        # GML requires string/integer attributes; escape if needed
        for n in H.nodes():
            for k in list(H.nodes[n].keys()):
                val = H.nodes[n][k]
                if not isinstance(val, (str, int, float)):
                    H.nodes[n][k] = str(val)
        nx.write_gml(H, path)
        print(f"Saved {path}")


def write_edge_file(G: nx.Graph, path: Path, header: str = "Node1\tNode2\tWeight\n"):
    """Write edge list (Node1, Node2, Weight) for graph G."""
    if G is None or G.number_of_edges() == 0:
        with open(path, "w") as f:
            f.write(header)
        return
    with open(path, "w") as f:
        f.write(header)
        for u, v, d in sorted(G.edges(data=True), key=lambda e: -e[2].get("weight", 1)):
            w = d.get("weight", 1)
            f.write(f"{u}\t{v}\t{w}\n")
    print(f"Saved {path}")


def _shorten_label(label: str, max_len: int = 22) -> str:
    if len(label) <= max_len:
        return label
    return label[: max_len - 3].rstrip() + "…"


def _compute_layout(H):
    """Compute a clean layout with good spacing."""
    try:
        pos = nx.kamada_kawai_layout(H, scale=2.0, weight="weight")
    except Exception:
        pos = nx.spring_layout(H, k=2.0, iterations=100, seed=42, weight="weight")
    # Center and scale to [-1, 1] for consistent aspect
    pos = np.asarray([pos[n] for n in H.nodes()])
    if len(pos) > 0:
        pos = pos - pos.mean(axis=0)
        s = np.abs(pos).max()
        if s > 0:
            pos = pos / s
    return dict(zip(H.nodes(), pos))


def _compute_layout_multicomponent(H):
    """Layout for disconnected graph: layout each component then place with offset so they don't overlap."""
    comps = list(nx.connected_components(H))
    if len(comps) <= 1:
        return _compute_layout(H)
    pos = {}
    cols = max(1, int(np.ceil(np.sqrt(len(comps)))))
    gap = 2.4
    for idx, comp in enumerate(comps):
        sub = H.subgraph(comp).copy()
        try:
            p = nx.kamada_kawai_layout(sub, scale=1.0, weight="weight")
        except Exception:
            p = nx.spring_layout(sub, k=1.2, iterations=80, seed=42 + idx, weight="weight")
        arr = np.asarray([p[n] for n in sub.nodes()])
        if len(arr) > 0:
            arr = arr - arr.mean(axis=0)
            s = np.abs(arr).max()
            if s > 0:
                arr = arr / (s + 1e-6) * 0.9
        row, col = idx // cols, idx % cols
        offset = np.array([col * gap, row * gap])
        for i, n in enumerate(sub.nodes()):
            pos[n] = arr[i] + offset
    # Center overall
    if pos:
        all_pos = np.array(list(pos.values()))
        center = all_pos.mean(axis=0)
        pos = {n: pos[n] - center for n in pos}
        s = np.abs(np.array(list(pos.values()))).max()
        if s > 0:
            scale = 1.0 / s
            pos = {n: pos[n] * scale for n in pos}
    return pos


def load_author_regions(data_dir: Path) -> dict[str, str]:
    """Load author -> region from data/author_regions.txt."""
    path = data_dir / "author_regions.txt"
    if not path.exists():
        return {}
    out = {}
    with open(path) as f:
        next(f)  # header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                out[parts[0]] = parts[1]
    return out


def draw_author_network(
    G,
    title: str,
    out_path: Path,
    author_region: dict[str, str],
    *,
    min_edge_weight: int = 1,
    node_size_scale: float = 70,
    figsize=(14, 11),
    max_nodes_for_labels: int = 120,
):
    """Draw author network with all components (no largest-CC collapse), colored by affiliation region."""
    if G is None or G.number_of_nodes() == 0:
        fig, ax = plt.subplots(figsize=figsize, facecolor="#fafafa")
        ax.set_facecolor("#fafafa")
        ax.text(0.5, 0.5, "No network data. Run fetch_pubmed.py first.", ha="center", va="center", fontsize=12)
        ax.axis("off")
        fig.savefig(out_path, dpi=200, bbox_inches="tight", facecolor="#fafafa")
        plt.close()
        return

    H = nx.Graph()
    for u, v, d in G.edges(data=True):
        w = d.get("weight", 1)
        if w >= min_edge_weight:
            H.add_edge(u, v, weight=w)
    for u, v in list(H.edges()):
        H.add_node(u)
        H.add_node(v)
    if H.number_of_nodes() == 0:
        H = G.copy()

    pos = _compute_layout_multicomponent(H)
    weights = np.array([H.edges[u, v].get("weight", 1) for u, v in H.edges()])
    deg = dict(H.degree(weight="weight"))
    node_sizes = np.array([max(80, deg.get(n, 1) * node_size_scale) for n in H.nodes()])
    nodes = list(H.nodes())
    regions = [author_region.get(n, "Others") for n in nodes]
    node_colors = [get_region_color(r) for r in regions]

    w_min = weights.min() if len(weights) else 0
    w_max = weights.max() if len(weights) else 0
    if w_max > w_min:
        edge_widths = 0.8 + 3.5 * (weights - w_min) / (w_max - w_min)
    else:
        edge_widths = np.ones(len(weights)) * 1.5 if len(weights) else []
    edge_widths = np.clip(edge_widths, 0.6, 4.5) if len(edge_widths) else []
    if w_max > w_min:
        edge_alphas = 0.35 + 0.5 * (weights - w_min) / (w_max - w_min)
    else:
        edge_alphas = np.ones(len(weights)) * 0.5 if len(weights) else []
    edge_color = "#2d3748"

    fig, ax = plt.subplots(figsize=figsize, facecolor="#fafafa")
    ax.set_facecolor("#fafafa")

    for (u, v), w, alpha in zip(H.edges(), edge_widths, edge_alphas):
        ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]], color=edge_color, linewidth=float(w), alpha=float(alpha), zorder=0)

    for i, n in enumerate(nodes):
        ax.scatter(pos[n][0], pos[n][1], s=node_sizes[i], c=[node_colors[i]], edgecolors="white", linewidths=1.2, zorder=1, alpha=0.92)

    n_show = min(max_nodes_for_labels, H.number_of_nodes())
    sorted_by_deg = sorted(nodes, key=lambda x: deg.get(x, 0), reverse=True)[:n_show]
    labels = {n: _shorten_label(n, 20) for n in sorted_by_deg}
    texts = []
    for n, lbl in labels.items():
        t = ax.annotate(lbl, (pos[n][0], pos[n][1]), fontsize=5.5, ha="center", va="center", fontweight="500", color="#1a202c", zorder=2,
                        bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="none", alpha=0.9))
        texts.append(t)
    try:
        from adjustText import adjust_text
        xs = np.array([pos[n][0] for n in labels])
        ys = np.array([pos[n][1] for n in labels])
        adjust_text(texts, x=xs, y=ys, ax=ax, expand=(1.1, 1.1), arrowprops=None, force_text=(0.2, 0.4), force_points=(0.1, 0.2))
    except Exception:
        pass

    all_pos = list(pos.values())
    if all_pos:
        arr = np.array(all_pos)
        lim = max(np.abs(arr).max() * 1.15, 0.01)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, fontsize=15, fontweight="600", color="#1a202c", pad=14)

    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    legend_elements = [
        Line2D([0], [0], color=edge_color, linewidth=2, alpha=0.7, label="Co-occurrence (thicker = more)"),
        Patch(facecolor=get_region_color("China"), label="China"),
        Patch(facecolor=get_region_color("US"), label="US"),
        Patch(facecolor=get_region_color("UK"), label="UK"),
        Patch(facecolor=get_region_color("Europe"), label="Europe"),
        Patch(facecolor=get_region_color("Japan"), label="Japan"),
        Patch(facecolor=get_region_color("Others"), label="Others"),
    ]
    ax.legend(handles=legend_elements, loc="upper left", frameon=True, fancybox=True, framealpha=0.95, fontsize=7, ncol=2)
    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight", facecolor="#fafafa")
    plt.close()
    print(f"Saved {out_path}")


def draw_network(
    G,
    title: str,
    out_path: Path,
    *,
    min_edge_weight: int = 1,
    node_size_scale: float = 50,
    figsize=(14, 11),
    shorten_labels: bool = False,
    max_nodes_for_labels: int = 80,
    use_largest_cc: bool = True,
):
    if G is None or G.number_of_nodes() == 0:
        fig, ax = plt.subplots(figsize=figsize, facecolor="#fafafa")
        ax.set_facecolor("#fafafa")
        ax.text(0.5, 0.5, "No network data. Run fetch_pubmed.py first.", ha="center", va="center", fontsize=12)
        ax.axis("off")
        fig.savefig(out_path, dpi=200, bbox_inches="tight", facecolor="#fafafa")
        plt.close()
        return

    # Build subgraph with edges above threshold
    H = nx.Graph()
    for u, v, d in G.edges(data=True):
        w = d.get("weight", 1)
        if w >= min_edge_weight:
            H.add_edge(u, v, weight=w)
    for u, v in list(H.edges()):
        H.add_node(u)
        H.add_node(v)
    if H.number_of_nodes() == 0:
        H = G.copy()
    elif use_largest_cc:
        H = H.subgraph(max(nx.connected_components(H), key=len)).copy()

    pos = _compute_layout(H)
    weights = np.array([H.edges[u, v].get("weight", 1) for u, v in H.edges()])
    deg = dict(H.degree(weight="weight"))
    node_sizes = np.array([max(80, deg.get(n, 1) * node_size_scale) for n in H.nodes()])
    # Normalize for color (log scale for degree)
    deg_vals = np.array([deg.get(n, 1) for n in H.nodes()])
    deg_log = np.log1p(deg_vals)
    if deg_log.max() > deg_log.min():
        node_colors = (deg_log - deg_log.min()) / (deg_log.max() - deg_log.min())
    else:
        node_colors = np.ones(H.number_of_nodes()) * 0.5

    # Edge width: scale weight with cap
    w_min, w_max = weights.min(), weights.max()
    if w_max > w_min:
        edge_widths = 0.8 + 3.5 * (weights - w_min) / (w_max - w_min)
    else:
        edge_widths = np.ones(len(weights)) * 1.5
    edge_widths = np.clip(edge_widths, 0.6, 4.5)
    # Edge color: light to dark by weight
    if w_max > w_min:
        edge_alphas = 0.35 + 0.5 * (weights - w_min) / (w_max - w_min)
    else:
        edge_alphas = np.ones(len(weights)) * 0.5

    # Colormap: nodes (blue–teal), edges (neutral gray)
    cmap_nodes = plt.cm.viridis
    edge_color = "#2d3748"

    fig, ax = plt.subplots(figsize=figsize, facecolor="#fafafa")
    ax.set_facecolor("#fafafa")

    # Draw edges first (so they sit under nodes)
    for (u, v), w, alpha in zip(H.edges(), edge_widths, edge_alphas):
        x = [pos[u][0], pos[v][0]]
        y = [pos[u][1], pos[v][1]]
        ax.plot(x, y, color=edge_color, linewidth=float(w), alpha=float(alpha), zorder=0)

    # Draw nodes with white border for clarity
    nodes = list(H.nodes())
    for i, n in enumerate(nodes):
        ax.scatter(
            pos[n][0],
            pos[n][1],
            s=node_sizes[i],
            c=[cmap_nodes(node_colors[i])],
            edgecolors="white",
            linewidths=1.2,
            zorder=1,
            alpha=0.92,
        )

    # Labels: show for top-degree nodes to avoid clutter
    if H.number_of_nodes() <= max_nodes_for_labels:
        label_nodes = set(nodes)
    else:
        sorted_by_deg = sorted(nodes, key=lambda x: deg.get(x, 0), reverse=True)
        label_nodes = set(sorted_by_deg[: max_nodes_for_labels])
    labels = {}
    for n in label_nodes:
        lbl = _shorten_label(n, 20) if shorten_labels else n
        labels[n] = lbl

    for n, lbl in labels.items():
        x, y = pos[n]
        ax.annotate(
            lbl,
            (x, y),
            fontsize=6,
            ha="center",
            va="center",
            fontweight="500",
            color="#1a202c",
            zorder=2,
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="none", alpha=0.85),
        )

    ax.set_xlim(-1.15, 1.15)
    ax.set_ylim(-1.15, 1.15)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, fontsize=15, fontweight="600", color="#1a202c", pad=14)

    # Legend: node size and edge weight
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    legend_elements = [
        Line2D([0], [0], color=edge_color, linewidth=2, alpha=0.7, label="Co-occurrence (thicker = more)"),
        Patch(facecolor=cmap_nodes(0.2), edgecolor="white", linewidth=1, label="Node size ∝ activity"),
    ]
    ax.legend(handles=legend_elements, loc="upper left", frameon=True, fancybox=True, framealpha=0.95, fontsize=8)

    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight", facecolor="#fafafa")
    plt.close()
    print(f"Saved {out_path}")


def draw_institution_network(
    G,
    title: str,
    out_path: Path,
    *,
    min_edge_weight: int = 1,
    node_size_scale: float = 28,
    figsize=(14, 11),
    max_nodes_for_labels: int = 200,
):
    """Draw institution network with region colors, abbreviated labels, and overlap-resistant labels."""
    if G is None or G.number_of_nodes() == 0:
        fig, ax = plt.subplots(figsize=figsize, facecolor="#fafafa")
        ax.set_facecolor("#fafafa")
        ax.text(0.5, 0.5, "No network data. Run fetch_pubmed.py first.", ha="center", va="center", fontsize=12)
        ax.axis("off")
        fig.savefig(out_path, dpi=200, bbox_inches="tight", facecolor="#fafafa")
        plt.close()
        return

    H = nx.Graph()
    for u, v, d in G.edges(data=True):
        w = d.get("weight", 1)
        if w >= min_edge_weight:
            H.add_edge(u, v, weight=w)
    for u, v in list(H.edges()):
        H.add_node(u)
        H.add_node(v)
    if H.number_of_nodes() == 0:
        H = G.copy()
    else:
        H = H.subgraph(max(nx.connected_components(H), key=len)).copy()

    pos = _compute_layout(H)
    weights = np.array([H.edges[u, v].get("weight", 1) for u, v in H.edges()])
    deg = dict(H.degree(weight="weight"))
    node_sizes = np.array([max(60, deg.get(n, 1) * node_size_scale) for n in H.nodes()])
    nodes = list(H.nodes())

    # Node colors by region
    regions = [get_region(n) for n in nodes]
    node_colors = [get_region_color(r) for r in regions]

    w_min, w_max = weights.min(), weights.max()
    if w_max > w_min:
        edge_widths = 0.8 + 3.5 * (weights - w_min) / (w_max - w_min)
    else:
        edge_widths = np.ones(len(weights)) * 1.5
    edge_widths = np.clip(edge_widths, 0.6, 4.5)
    if w_max > w_min:
        edge_alphas = 0.35 + 0.5 * (weights - w_min) / (w_max - w_min)
    else:
        edge_alphas = np.ones(len(weights)) * 0.5
    edge_color = "#2d3748"

    fig, ax = plt.subplots(figsize=figsize, facecolor="#fafafa")
    ax.set_facecolor("#fafafa")

    for (u, v), w, alpha in zip(H.edges(), edge_widths, edge_alphas):
        ax.plot(
            [pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
            color=edge_color, linewidth=float(w), alpha=float(alpha), zorder=0,
        )

    for i, n in enumerate(nodes):
        ax.scatter(
            pos[n][0], pos[n][1],
            s=node_sizes[i],
            c=[node_colors[i]],
            edgecolors="white",
            linewidths=1.2,
            zorder=1,
            alpha=0.92,
        )

    # Labels: abbreviated, then adjust to avoid overlap
    n_show = min(max_nodes_for_labels, H.number_of_nodes())
    sorted_by_deg = sorted(nodes, key=lambda x: deg.get(x, 0), reverse=True)[:n_show]
    texts = []
    for n in sorted_by_deg:
        lbl = abbreviate_institution(n, max_len=26)
        t = ax.annotate(
            lbl,
            (pos[n][0], pos[n][1]),
            fontsize=5.5,
            ha="center",
            va="center",
            fontweight="500",
            color="#1a202c",
            zorder=2,
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="none", alpha=0.9),
        )
        texts.append(t)

    try:
        from adjustText import adjust_text
        xs = np.array([pos[n][0] for n in sorted_by_deg])
        ys = np.array([pos[n][1] for n in sorted_by_deg])
        adjust_text(
            texts,
            x=xs,
            y=ys,
            ax=ax,
            expand=(1.15, 1.15),
            arrowprops=None,  # no arrows, labels only (reduces overlap)
            force_text=(0.3, 0.5),
            force_points=(0.15, 0.25),
        )
    except Exception:
        pass  # fallback: labels stay at node positions

    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.2, 1.2)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, fontsize=15, fontweight="600", color="#1a202c", pad=14)

    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    legend_elements = [
        Line2D([0], [0], color=edge_color, linewidth=2, alpha=0.7, label="Co-occurrence (thicker = more)"),
        Patch(facecolor=get_region_color("China"), label="China"),
        Patch(facecolor=get_region_color("US"), label="US"),
        Patch(facecolor=get_region_color("UK"), label="UK"),
        Patch(facecolor=get_region_color("Europe"), label="Europe"),
        Patch(facecolor=get_region_color("Japan"), label="Japan"),
        Patch(facecolor=get_region_color("Others"), label="Others"),
    ]
    ax.legend(handles=legend_elements, loc="upper left", frameon=True, fancybox=True, framealpha=0.95, fontsize=7, ncol=2)

    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight", facecolor="#fafafa")
    plt.close()
    print(f"Saved {out_path}")


def main():
    data_dir = Path(__file__).resolve().parent / "data"
    out_dir = Path(__file__).resolve().parent / "outputs"
    out_dir.mkdir(exist_ok=True)

    if not (data_dir / "author_edges.txt").exists():
        print("No data found. Running fetch_pubmed to fetch PubMed data...")
        from fetch_pubmed import main as fetch_main

        fetch_main()

    G_author = load_author_network(data_dir)
    G_inst = load_institution_network(data_dir)

    # Restrict to top 200 nodes by degree and write edge files
    H_author = top_k_subgraph(G_author, TOP_NODES)
    H_inst = top_k_subgraph(G_inst, TOP_NODES)

    write_edge_file(H_author, out_dir / "author_edges_top200.txt", header="Author1\tAuthor2\tWeight\n")
    write_edge_file(H_inst, out_dir / "institution_edges_top200.txt", header="Institution1\tInstitution2\tWeight\n")

    author_region = load_author_regions(data_dir)
    if not author_region and (data_dir / "records.json").exists():
        # Build from records if author_regions.txt was not generated yet
        import json
        from collections import Counter, defaultdict
        with open(data_dir / "records.json") as f:
            records = json.load(f)
        author_affs = defaultdict(list)
        for rec in records:
            for name, affs in rec.get("author_affiliations", []):
                if name:
                    author_affs[name].extend(affs)
        for author, affs in author_affs.items():
            regions = [get_region(a) for a in affs if a]
            regions = [r for r in regions if r != "Others"]
            author_region[author] = Counter(regions).most_common(1)[0][0] if regions else "Others"

    draw_author_network(
        H_author,
        f"Metagenomics: Co-authorship network (top {TOP_NODES} by degree)",
        out_dir / "author_network.png",
        author_region,
        min_edge_weight=1,
        node_size_scale=70,
        figsize=(14, 11),
        max_nodes_for_labels=120,
    )
    draw_institution_network(
        H_inst,
        f"Metagenomics: University co-affiliation network (top {TOP_NODES} by degree)",
        out_dir / "institution_network.png",
        min_edge_weight=1,
        node_size_scale=28,
        figsize=(14, 11),
        max_nodes_for_labels=200,
    )
    # Export GML for each network (with node attributes: region, label)
    _write_gml(H_author, out_dir / "author_network.gml", author_region, node_label=lambda n: n)
    _write_gml(H_inst, out_dir / "institution_network.gml", {n: get_region(n) for n in H_inst.nodes()}, node_label=lambda n: abbreviate_institution(n))
    print("Done. See outputs/ for PNGs, edge files (*_top200.txt), and GML files.")


if __name__ == "__main__":
    main()
