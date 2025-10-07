"""Interactive plotting utilities based on Plotly without runtime dependency."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable, List, Sequence

from .analysis import DendrogramNode

PLOTLY_CDN = "https://cdn.plot.ly/plotly-2.27.0.min.js"


def _write_plotly_html(fig: dict, output: str | Path) -> None:
    target = Path(output)
    target.parent.mkdir(parents=True, exist_ok=True)
    figure_json = json.dumps(fig)
    html = f"""
<!DOCTYPE html>
<html lang=\"en\">
<head>
    <meta charset=\"utf-8\" />
    <script src=\"{PLOTLY_CDN}\"></script>
    <title>Interactive plot</title>
    <style>
        html, body {{ height: 100%; margin: 0; }}
        #plot {{ width: 100%; height: 100%; }}
    </style>
</head>
<body>
    <div id=\"plot\"></div>
    <script>
        const figure = {figure_json};
        Plotly.newPlot('plot', figure.data, figure.layout, {{responsive: true}});
    </script>
</body>
</html>
"""
    target.write_text(html, encoding="utf-8")


def pca_scatter(scores: Sequence[Sequence[float]], individuals: Sequence[str], explained: Sequence[float], output: str | Path) -> None:
    if not scores:
        raise ValueError("scores must contain at least one entry")
    x = [row[0] for row in scores]
    y = [row[1] if len(row) > 1 else 0.0 for row in scores]
    hover = [f"{ind}<br>PC1: {row[0]:.3f}<br>PC2: {row[1]:.3f}" for ind, row in zip(individuals, scores)]
    var1 = explained[0] * 100 if explained else 0.0
    var2 = explained[1] * 100 if len(explained) > 1 else 0.0
    fig = {
        "data": [
            {
                "type": "scatter",
                "mode": "markers",
                "x": x,
                "y": y,
                "text": hover,
                "hovertemplate": "%{text}<extra></extra>",
            }
        ],
        "layout": {
            "title": "PCA Scatter",
            "xaxis": {"title": f"PC1 ({var1:.2f}% var)"},
            "yaxis": {"title": f"PC2 ({var2:.2f}% var)"},
            "template": "plotly_white",
        },
    }
    _write_plotly_html(fig, output)


def gwas_manhattan(results: Sequence[dict], output: str | Path) -> None:
    if not results:
        raise ValueError("results must not be empty")
    positions = list(range(len(results)))
    neg_log_p = [_neg_log10(entry.get("pvalue", 1.0)) for entry in results]
    markers = [entry.get("marker", f"M{i}") for i, entry in enumerate(results)]
    fig = {
        "data": [
            {
                "type": "scatter",
                "mode": "markers",
                "x": positions,
                "y": neg_log_p,
                "text": markers,
                "hovertemplate": "Marker: %{text}<br>-log10(p): %{y:.3f}<extra></extra>",
            }
        ],
        "layout": {
            "title": "GWAS Manhattan Plot",
            "xaxis": {"title": "Variant Index"},
            "yaxis": {"title": "-log10(pvalue)"},
            "template": "plotly_white",
        },
    }
    _write_plotly_html(fig, output)


def _assign_leaf_positions(node: DendrogramNode, order: List[str], positions: dict) -> None:
    if node.is_leaf():
        positions[node.name] = len(order)
        order.append(node.name)
        return
    if node.left is not None:
        _assign_leaf_positions(node.left, order, positions)
    if node.right is not None:
        _assign_leaf_positions(node.right, order, positions)
def dendrogram_plot(root: DendrogramNode, output: str | Path) -> None:
    order: List[str] = []
    positions: dict = {}
    _assign_leaf_positions(root, order, positions)
    leaves = order
    if not leaves:
        raise ValueError("Dendrogram must contain at least one leaf")
    max_distance = max((node.distance for node in _iter_nodes(root)), default=0.0)
    if max_distance == 0:
        max_distance = 1.0

    x_coords: List[float] = []
    y_coords: List[float] = []
    for name in leaves:
        x_coords.append(positions[name])
        y_coords.append(0.0)

    lines_x: List[List[float]] = []
    lines_y: List[List[float]] = []

    def _traverse(node: DendrogramNode) -> None:
        if node.is_leaf():
            return
        if node.left is not None:
            _traverse(node.left)
        if node.right is not None:
            _traverse(node.right)
        if node.left is None or node.right is None:
            return
        left_names = node.left.members()
        right_names = node.right.members()
        left_x = [positions[name] for name in left_names]
        right_x = [positions[name] for name in right_names]
        left_y = node.left.distance
        right_y = node.right.distance
        current_y = node.distance
        lines_x.append([min(left_x), max(left_x)])
        lines_y.append([left_y, left_y])
        lines_x.append([min(right_x), max(right_x)])
        lines_y.append([right_y, right_y])
        lines_x.append([sum(left_x) / len(left_x), sum(right_x) / len(right_x)])
        lines_y.append([current_y, current_y])
        lines_x.append([sum(left_x) / len(left_x), sum(left_x) / len(left_x)])
        lines_y.append([left_y, current_y])
        lines_x.append([sum(right_x) / len(right_x), sum(right_x) / len(right_x)])
        lines_y.append([right_y, current_y])

    _traverse(root)

    traces = [
        {
            "type": "scatter",
            "mode": "markers",
            "x": x_coords,
            "y": y_coords,
            "text": leaves,
            "hovertemplate": "%{text}<extra></extra>",
        }
    ]

    for x_vals, y_vals in zip(lines_x, lines_y):
        traces.append(
            {
                "type": "scatter",
                "mode": "lines",
                "x": x_vals,
                "y": y_vals,
                "line": {"color": "#636efa"},
                "hoverinfo": "skip",
            }
        )

    fig = {
        "data": traces,
        "layout": {
            "title": "Dendrogram",
            "xaxis": {"tickvals": list(range(len(leaves))), "ticktext": leaves, "tickangle": 45},
            "yaxis": {"title": "Distance", "range": [0, max_distance * 1.05]},
            "template": "plotly_white",
        },
    }
    _write_plotly_html(fig, output)


def _iter_nodes(node: DendrogramNode) -> Iterable[DendrogramNode]:
    stack = [node]
    while stack:
        current = stack.pop()
        yield current
        if current.left is not None:
            stack.append(current.left)
        if current.right is not None:
            stack.append(current.right)


def admixture_plot(q_matrix: Sequence[Sequence[float]], individuals: Sequence[str], output: str | Path) -> None:
    if not q_matrix:
        raise ValueError("q_matrix must not be empty")
    n_components = len(q_matrix[0])
    if n_components == 0:
        raise ValueError("q_matrix must contain at least one component")
    traces = []
    for comp in range(n_components):
        values = [row[comp] for row in q_matrix]
        traces.append(
            {
                "type": "bar",
                "x": individuals,
                "y": values,
                "name": f"Cluster {comp + 1}",
            }
        )
    layout = {
        "barmode": "stack",
        "title": "Admixture Proportions",
        "xaxis": {"title": "Individuals", "tickangle": 45},
        "yaxis": {"title": "Proportion", "range": [0, 1]},
        "template": "plotly_white",
    }
    fig = {"data": traces, "layout": layout}
    _write_plotly_html(fig, output)


def _neg_log10(value: float) -> float:
    if value <= 0:
        return 0.0
    import math

    return -math.log10(value)
