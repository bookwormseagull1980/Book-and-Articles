"""Create manuscript figures from a three-dimensional field-generation run."""

from __future__ import annotations

import argparse
import csv
from html import escape
import math
from pathlib import Path
import shutil
import statistics
import sys
from typing import Iterable, Mapping, Sequence

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from cosmic_cg.intermediate import write_output_inventory
from cosmic_cg.io.tables import read_rows

try:
    from reportlab.lib import colors
    from reportlab.pdfgen import canvas as pdf_canvas
except ModuleNotFoundError:  # pragma: no cover
    colors = None
    pdf_canvas = None


WIDTH = 7.2 * 72.0
HEIGHT = 5.2 * 72.0
PANEL_GAP = 46.0
COLORS = {
    "ink": "#1a1a1a",
    "muted": "#666666",
    "grid": "#d7d7d7",
    "blue": "#2369a7",
    "red": "#a94438",
    "green": "#2f7d4f",
    "purple": "#6d5aa8",
    "orange": "#b36b22",
    "gray": "#9a9a9a",
    "lightblue": "#d9e8f5",
    "lightred": "#efd9d5",
    "lightgreen": "#d8eadf",
    "lightpurple": "#e6e0f0",
}


class Figure:
    def __init__(self, path: Path, title: str):
        self.path = path
        self.title = title
        self.svg: list[str] = [
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH:.0f}" height="{HEIGHT:.0f}" viewBox="0 0 {WIDTH:.0f} {HEIGHT:.0f}">',
            '<rect x="0" y="0" width="100%" height="100%" fill="white"/>',
        ]
        self.pdf = None
        if pdf_canvas is not None:
            self.pdf = pdf_canvas.Canvas(str(path.with_suffix(".pdf")), pagesize=(WIDTH, HEIGHT))
            self.pdf.setTitle(title)
            self.pdf.setFillColor(_pdf_color(COLORS["ink"]))
            self.pdf.setStrokeColor(_pdf_color(COLORS["ink"]))

    def save(self) -> None:
        self.svg.append("</svg>")
        self.path.with_suffix(".svg").write_text("\n".join(self.svg) + "\n", encoding="utf-8")
        if self.pdf is not None:
            self.pdf.showPage()
            self.pdf.save()

    def text(self, x: float, y: float, text: str, size: float = 9.0, color: str = "ink", anchor: str = "start", vertical: bool = False) -> None:
        anchor_map = {"start": "start", "middle": "middle", "end": "end"}
        transform = f' transform="rotate(-90 {x:.2f} {y:.2f})"' if vertical else ""
        self.svg.append(
            f'<text x="{x:.2f}" y="{y:.2f}" font-family="Helvetica,Arial,sans-serif" '
            f'font-size="{size:.2f}" fill="{COLORS[color]}" text-anchor="{anchor_map[anchor]}"{transform}>{escape(text)}</text>'
        )
        if self.pdf is not None:
            self.pdf.saveState()
            self.pdf.setFillColor(_pdf_color(COLORS[color]))
            self.pdf.setFont("Helvetica", size)
            px, py = x, HEIGHT - y
            if vertical:
                self.pdf.translate(px, py)
                self.pdf.rotate(90)
                px, py = 0, 0
            if anchor == "middle":
                self.pdf.drawCentredString(px, py, text)
            elif anchor == "end":
                self.pdf.drawRightString(px, py, text)
            else:
                self.pdf.drawString(px, py, text)
            self.pdf.restoreState()

    def line(self, x1: float, y1: float, x2: float, y2: float, color: str = "ink", width: float = 0.8) -> None:
        self.svg.append(
            f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" '
            f'stroke="{COLORS[color]}" stroke-width="{width:.2f}" stroke-linecap="round"/>'
        )
        if self.pdf is not None:
            self.pdf.setStrokeColor(_pdf_color(COLORS[color]))
            self.pdf.setLineWidth(width)
            self.pdf.line(x1, HEIGHT - y1, x2, HEIGHT - y2)

    def rect(self, x: float, y: float, w: float, h: float, fill: str = "gray", stroke: str | None = None, width: float = 0.7) -> None:
        stroke_value = "none" if stroke is None else COLORS[stroke]
        self.svg.append(
            f'<rect x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" '
            f'fill="{COLORS[fill]}" stroke="{stroke_value}" stroke-width="{width:.2f}"/>'
        )
        if self.pdf is not None:
            self.pdf.setFillColor(_pdf_color(COLORS[fill]))
            self.pdf.setStrokeColor(_pdf_color(COLORS[stroke] if stroke else COLORS[fill]))
            self.pdf.setLineWidth(width)
            self.pdf.rect(x, HEIGHT - y - h, w, h, stroke=1 if stroke else 0, fill=1)

    def circle(self, x: float, y: float, r: float, fill: str = "blue", stroke: str | None = None, width: float = 0.6) -> None:
        stroke_value = "none" if stroke is None else COLORS[stroke]
        self.svg.append(
            f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{r:.2f}" '
            f'fill="{COLORS[fill]}" stroke="{stroke_value}" stroke-width="{width:.2f}"/>'
        )
        if self.pdf is not None:
            self.pdf.setFillColor(_pdf_color(COLORS[fill]))
            self.pdf.setStrokeColor(_pdf_color(COLORS[stroke] if stroke else COLORS[fill]))
            self.pdf.setLineWidth(width)
            self.pdf.circle(x, HEIGHT - y, r, stroke=1 if stroke else 0, fill=1)

    def polyline(self, points: Sequence[tuple[float, float]], color: str = "blue", width: float = 1.5) -> None:
        if not points:
            return
        coords = " ".join(f"{x:.2f},{y:.2f}" for x, y in points)
        self.svg.append(
            f'<polyline points="{coords}" fill="none" stroke="{COLORS[color]}" '
            f'stroke-width="{width:.2f}" stroke-linejoin="round" stroke-linecap="round"/>'
        )
        if self.pdf is not None and len(points) > 1:
            self.pdf.setStrokeColor(_pdf_color(COLORS[color]))
            self.pdf.setLineWidth(width)
            path = self.pdf.beginPath()
            path.moveTo(points[0][0], HEIGHT - points[0][1])
            for x, y in points[1:]:
                path.lineTo(x, HEIGHT - y)
            self.pdf.drawPath(path, stroke=1, fill=0)

    def polygon(self, points: Sequence[tuple[float, float]], fill: str = "lightblue", stroke: str | None = None) -> None:
        if not points:
            return
        coords = " ".join(f"{x:.2f},{y:.2f}" for x, y in points)
        stroke_value = "none" if stroke is None else COLORS[stroke]
        self.svg.append(f'<polygon points="{coords}" fill="{COLORS[fill]}" stroke="{stroke_value}" stroke-width="0.6"/>')
        if self.pdf is not None:
            self.pdf.setFillColor(_pdf_color(COLORS[fill]))
            self.pdf.setStrokeColor(_pdf_color(COLORS[stroke] if stroke else COLORS[fill]))
            path = self.pdf.beginPath()
            path.moveTo(points[0][0], HEIGHT - points[0][1])
            for x, y in points[1:]:
                path.lineTo(x, HEIGHT - y)
            path.close()
            self.pdf.drawPath(path, stroke=1 if stroke else 0, fill=1)


def main() -> int:
    args = parse_args()
    run_dir = Path(args.run_dir)
    paper_dir = Path(args.paper_figures_dir) if args.paper_figures_dir else None
    output_dir = Path(args.output_dir) if args.output_dir else paper_dir if paper_dir is not None else run_dir / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    if paper_dir is not None and paper_dir.resolve() != output_dir.resolve():
        paper_dir.mkdir(parents=True, exist_ok=True)

    tables = _load_tables(run_dir)
    renderers = {1: _figure1, 2: _figure2, 4: _figure4, 5: _figure5, 6: _figure6, 8: _figure8}
    written = [renderers[number](output_dir, tables) for number in dict.fromkeys(args.figures)]
    if paper_dir is not None and paper_dir.resolve() != output_dir.resolve():
        for path in written:
            for suffix in (".svg", ".pdf"):
                src = path.with_suffix(suffix)
                if src.is_file():
                    shutil.copy2(src, paper_dir / src.name)
    if output_dir.resolve().is_relative_to(run_dir.resolve()):
        write_output_inventory(run_dir)
    for path in written:
        print(path.with_suffix(".svg"))
        if path.with_suffix(".pdf").is_file():
            print(path.with_suffix(".pdf"))
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--output-dir", default="", help="Override the rendering destination.")
    parser.add_argument("--paper-figures-dir", default="", help="Render here unless --output-dir is set; otherwise copy here. Defaults to RUN_DIR/figures.")
    parser.add_argument("--figures", type=int, nargs="+", choices=(1, 2, 4, 5, 6, 8), default=(1, 2, 4, 5, 6, 8))
    return parser.parse_args()


def _load_tables(root: Path) -> dict[str, list[dict[str, str]]]:
    return {
        "window": read_rows(root / "tables" / "field_window_table.csv"),
        "mode": read_rows(root / "tables" / "field_mode_table.csv"),
        "object": read_rows(root / "tables" / "field_object_table.csv"),
        "transfer": read_rows(root / "tables" / "field_transfer_table.csv"),
        "summary": read_rows(root / "tables" / "field_relation_summary_table.csv"),
        "epoch": read_rows(root / "tables" / "field_epoch_relation_table.csv"),
        "closure": read_rows(root / "tables" / "field_closure_table.csv"),
        "cost": read_rows(root / "tables" / "field_cost_table.csv"),
    }


def _figure1(output_dir: Path, tables: Mapping[str, Sequence[Mapping[str, str]]]) -> Path:
    path = output_dir / "fig01_state_architecture"
    fig = Figure(path, "Finite-window description")
    fig.text(34, 34, "Finite-window description", 13)
    fig.text(34, 54, "Retained modes, residual power and peak statistics", 8, "muted")
    nodes = [
        ("Active field", "k_NL, M_NL, S_AW", 70, 112, "lightblue"),
        ("Residual tail", "F_T, shell sectors", 282, 112, "lightred"),
        ("Objects", "M_h, r_peak, B_frac", 70, 238, "lightgreen"),
        ("Peak correlations", "A_fld, G_fld, spin proxy", 282, 238, "lightpurple"),
    ]
    for title, subtitle, x, y, color in nodes:
        fig.rect(x, y, 168, 72, color, "grid")
        fig.text(x + 12, y + 28, title, 10)
        fig.text(x + 12, y + 49, subtitle, 7, "muted")
    fig.rect(176, 192, 168, 34, "gray", "grid")
    fig.text(260, 214, "State(a)", 14, "ink", "middle")
    for x1, y1, x2, y2 in ((154, 184, 218, 192), (366, 184, 302, 192), (154, 238, 218, 226), (366, 238, 302, 226)):
        fig.line(x1, y1, x2, y2, "muted", 1.0)
    fig.save()
    return path


def _figure2(output_dir: Path, tables: Mapping[str, Sequence[Mapping[str, str]]]) -> Path:
    path = output_dir / "fig02_active_window_history"
    fig = Figure(path, "Active-window, scale and mass history")
    fig.text(34, 34, "Active-window, scale and mass history", 13)
    rows = _sort_by_a(tables["window"])
    a = _values(rows, "a")
    k = _values(rows, "k_nl_hmpc")
    panels = _grid_panels(2, 2)
    series = [
        ("A", "k_NL [h Mpc^-1]", k, True, "blue"),
        ("B", "lambda_NL [Mpc h^-1]", [2.0 * math.pi / value for value in k], True, "blue"),
        ("C", "M_NL [Msun h^-1]", _values(rows, "m_nl_msunh"), True, "green"),
        ("D", "span", _values(rows, "s_aw_global"), False, "purple"),
    ]
    sigma = _values(rows, "measured_sigma_delta")
    for panel, (label, y_label, y_values, log_y, color) in zip(panels, series):
        _draw_line_panel(fig, panel, a, y_values, y_label, color, log_x=True, log_y=log_y, panel_label=label)
        if label == "D":
            scaled = _rescale(sigma, y_values)
            points = _panel_points(panel, a, scaled, log_x=True, log_y=False, y_domain=_domain(y_values, False))
            fig.polyline(points, "orange", 1.1)
            _legend(fig, panel, [("S_AW", "purple"), ("sigma_delta affine", "orange")])
    z0 = _nearest(rows, "z", 0.0)
    fig.text(494, 51, f"z=0: lambda_NL={2.0 * math.pi / _float(z0, 'k_nl_hmpc'):.2f} Mpc h^-1", 8, "muted", "end")
    fig.text(284, 370, f"{len(rows)} epochs; a > 1 follows the prescribed continuation.", 7, "muted", "middle")
    fig.save()
    return path


def _figure4(output_dir: Path, tables: Mapping[str, Sequence[Mapping[str, str]]]) -> Path:
    path = output_dir / "fig04_residual_tail_audit"
    fig = Figure(path, "Residual spectral power")
    fig.text(34, 34, "Residual spectral power", 13)
    rows = _sort_by_a(tables["window"])
    panels = _stacked_panels(2)
    panel = panels[0]
    a = _values(rows, "a")
    active = _values(rows, "active_power_fraction")
    tail = _values(rows, "tail_power_fraction")
    _draw_axes(fig, panel, a, [0.0, 1.0], "log10 a", "power fraction", log_x=True, log_y=False, y_limits=(0.0, 1.0))
    points_tail = _panel_points(panel, a, tail, log_x=True, log_y=False, y_domain=(0.0, 1.0))
    bottom = panel[1] + panel[3]
    polygon = [(points_tail[0][0], bottom)] + points_tail + [(points_tail[-1][0], bottom)]
    fig.polygon(polygon, "lightred")
    fig.polyline(points_tail, "red", 1.6)
    fig.polyline(_panel_points(panel, a, active, log_x=True, log_y=False, y_domain=(0.0, 1.0)), "blue", 1.6)
    _legend(fig, panel, [("tail", "red"), ("active", "blue")])
    stage = _first(rows, lambda row: _float(row, "tail_power_fraction") >= _float(row, "active_power_fraction"))
    if stage:
        x = _scale_x(panel, _float(stage, "a"), _domain(a, True), True)
        fig.line(x, panel[1], x, bottom, "muted", 0.7)
        fig.text(x + 5, panel[1] + 13, f"z={_float(stage, 'z'):.0f}", 7, "muted")
    closure_rows = _sort_by_a(tables["closure"])
    if closure_rows:
        closure_a = _values(closure_rows, "a")
        closure_y = [math.log10(1.0 + max(0.0, _float(row, "closure_error"))) for row in closure_rows]
        _draw_line_panel(
            fig,
            panels[1],
            closure_a,
            closure_y,
            "log10(1+epsilon_P)",
            "purple",
            log_x=True,
            log_y=False,
            panel_label="B",
        )
    fig.save()
    return path


def _figure5(output_dir: Path, tables: Mapping[str, Sequence[Mapping[str, str]]]) -> Path:
    path = output_dir / "fig05_object_fate"
    fig = Figure(path, "Peak growth classification")
    fig.text(34, 34, "Peak growth classification", 13)
    fig.text(34, 54, "Scale-band and rank matching", 8, "muted")
    rows = _aggregate_summary_by_a(tables["summary"])
    panel = (62.0, 78.0, 420.0, 236.0)
    a = [row["a"] for row in rows]
    forming = [row["forming_fraction"] for row in rows]
    balanced = [row["balanced_fraction"] for row in rows]
    stripping = [row["stripping_fraction"] for row in rows]
    _draw_axes(fig, panel, a, [0.0, 1.0], "log10 a", "object fraction", log_x=True, log_y=False, y_limits=(0.0, 1.0))
    _draw_stacked_fraction(fig, panel, a, [forming, balanced, stripping], ["lightred", "lightgreen", "lightpurple"], ["red", "green", "purple"])
    _legend(fig, panel, [("forming", "red"), ("balanced", "green"), ("stripping", "purple")])
    stage = _first(rows, lambda row: row["balanced_fraction"] >= 0.5)
    if stage:
        x = _scale_x(panel, stage["a"], _domain(a, True), True)
        fig.line(x, panel[1], x, panel[1] + panel[3], "muted", 0.7)
        fig.text(x - 5, panel[1] + 13, f"balanced z={stage['z']:.1f}", 7, "muted", "end")
    fig.save()
    return path


def _figure6(output_dir: Path, tables: Mapping[str, Sequence[Mapping[str, str]]]) -> Path:
    path = output_dir / "fig06_transfer_channels"
    fig = Figure(path, "Transfer channels")
    fig.text(34, 34, "Transfer channels", 13)
    rows = _aggregate_transfer_by_a(tables["transfer"])
    panels = _stacked_panels(2)
    a = [row["a"] for row in rows]
    _draw_line_panel(fig, panels[0], a, [row["a_pq_median"] for row in rows], "median A_fld", "blue", log_x=True, log_y=True, panel_label="A")
    _draw_line_panel(fig, panels[1], a, [row["abs_g_pq_median"] for row in rows], "median |G_fld|", "purple", log_x=True, log_y=False, panel_label="B")
    fig.save()
    return path


def _figure8(output_dir: Path, tables: Mapping[str, Sequence[Mapping[str, str]]]) -> Path:
    path = output_dir / "fig08_field_relations"
    fig = Figure(path, "Three-dimensional field relations")
    fig.text(34, 34, "Three-dimensional field relations", 13)
    z0_modes = [row for row in tables["mode"] if abs(_float(row, "z")) < 1.0e-8]
    z0_objects = [row for row in tables["object"] if abs(_float(row, "z")) < 1.0e-8]
    z0_transfers = [row for row in tables["transfer"] if abs(_float(row, "z")) < 1.0e-8]
    panels = _grid_panels(2, 2)
    panels[2] = (panels[2][0], panels[2][1], panels[1][0] + panels[1][2] - panels[2][0], panels[2][3])
    fig.text(494, 51, f"z=0: {len(z0_modes)} shells, {len(z0_objects)} peaks, {len(z0_transfers)} pairs", 8, "muted", "end")
    _draw_line_panel(
        fig,
        panels[0],
        _values(z0_modes, "k_mid_hmpc"),
        _values(z0_modes, "power_fraction"),
        "shell power fraction",
        "blue",
        log_x=True,
        log_y=True,
        panel_label="A",
        x_label="k [h Mpc^-1]",
    )
    _draw_scatter_panel(
        fig,
        panels[1],
        [math.log10(_float(row, "mass_msunh")) for row in z0_objects],
        [_float(row, "b_h") for row in z0_objects],
        "log10 M_h [Msun h^-1]",
        "B_frac",
        "green",
        "B",
    )
    _draw_hist_panel(
        fig,
        panels[2],
        [abs(_float(row, "g_pq")) for row in z0_transfers],
        "abs G_fld",
        "pairs",
        "purple",
        "C",
        bins=16,
    )
    fig.save()
    return path


def _draw_line_panel(
    fig: Figure,
    panel: tuple[float, float, float, float],
    x_values: Sequence[float],
    y_values: Sequence[float],
    y_label: str,
    color: str,
    *,
    log_x: bool = False,
    log_y: bool = False,
    panel_label: str = "",
    x_label: str = "log10 a",
) -> None:
    _draw_axes(fig, panel, x_values, y_values, x_label, y_label, log_x=log_x, log_y=log_y)
    points = _panel_points(panel, x_values, y_values, log_x=log_x, log_y=log_y)
    fig.polyline(points, color, 1.4)
    for x, y in points:
        fig.circle(x, y, 2.2, color)
    if panel_label:
        fig.text(panel[0], panel[1] - 10, panel_label, 10)


def _draw_scatter_panel(
    fig: Figure,
    panel: tuple[float, float, float, float],
    x_values: Sequence[float],
    y_values: Sequence[float],
    x_label: str,
    y_label: str,
    color: str,
    panel_label: str,
) -> None:
    x_domain = _domain(x_values, False)
    y_domain = _domain(y_values + [0.0, 1.0], False)
    _draw_axes(fig, panel, x_values, y_values, x_label, y_label, y_limits=y_domain)
    for x_value, y_value in zip(x_values, y_values):
        fig.circle(_scale_x(panel, x_value, x_domain, False), _scale_y(panel, y_value, y_domain, False), 1.8, color)
    fig.line(panel[0], _scale_y(panel, 1.0, y_domain, False), panel[0] + panel[2], _scale_y(panel, 1.0, y_domain, False), "grid", 0.7)
    fig.text(panel[0], panel[1] - 10, panel_label, 10)


def _draw_hist_panel(
    fig: Figure,
    panel: tuple[float, float, float, float],
    values: Sequence[float],
    x_label: str,
    y_label: str,
    color: str,
    panel_label: str,
    *,
    bins: int,
) -> None:
    lo, hi = 0.0, 1.0
    counts = [0] * bins
    for value in values:
        if math.isfinite(value):
            index = min(bins - 1, max(0, int((value - lo) / (hi - lo) * bins)))
            counts[index] += 1
    y_domain = (0.0, max(1.0, 1.1 * max(counts)))
    _draw_axes(fig, panel, [lo, hi], counts, x_label, y_label,
               x_limits=(lo, hi), y_limits=y_domain)
    for i, count in enumerate(counts):
        x0 = panel[0] + i * panel[2] / bins
        width = panel[2] / bins * 0.82
        y = _scale_y(panel, count, y_domain, False)
        fig.rect(x0 + 2, y, width, panel[1] + panel[3] - y, color)
    fig.text(panel[0], panel[1] - 10, panel_label, 10)


def _draw_stacked_fraction(
    fig: Figure,
    panel: tuple[float, float, float, float],
    x_values: Sequence[float],
    layers: Sequence[Sequence[float]],
    fills: Sequence[str],
    strokes: Sequence[str],
) -> None:
    x_domain = _domain(x_values, True)
    baseline = [0.0 for _ in x_values]
    for layer, fill, stroke in zip(layers, fills, strokes):
        top = [base + value for base, value in zip(baseline, layer)]
        top_points = [(_scale_x(panel, x, x_domain, True), _scale_y(panel, y, (0.0, 1.0), False)) for x, y in zip(x_values, top)]
        bottom_points = [(_scale_x(panel, x, x_domain, True), _scale_y(panel, y, (0.0, 1.0), False)) for x, y in zip(reversed(x_values), reversed(baseline))]
        fig.polygon(top_points + bottom_points, fill)
        fig.polyline(top_points, stroke, 1.0)
        baseline = top


def _draw_axes(
    fig: Figure,
    panel: tuple[float, float, float, float],
    x_values: Sequence[float],
    y_values: Sequence[float],
    x_label: str,
    y_label: str,
    *,
    log_x: bool = False,
    log_y: bool = False,
    categorical: bool = False,
    x_limits: tuple[float, float] | None = None,
    y_limits: tuple[float, float] | None = None,
) -> None:
    x, y, w, h = panel
    x_domain = x_limits if x_limits is not None else _domain(x_values, log_x)
    y_domain = y_limits if y_limits is not None else _domain(y_values, log_y)
    fig.line(x, y + h, x + w, y + h, "ink", 0.8)
    fig.line(x, y, x, y + h, "ink", 0.8)
    for tick in _ticks(y_domain, log_y):
        py = _scale_y(panel, tick, y_domain, log_y)
        fig.line(x, py, x + w, py, "grid", 0.35)
        fig.text(x - 6, py + 3, _tick_label(tick), 8, "muted", "end")
    if not categorical:
        for tick in _ticks(x_domain, log_x):
            px = _scale_x(panel, tick, x_domain, log_x)
            fig.line(px, y + h, px, y + h + 4, "ink", 0.5)
            display_tick = math.log10(tick) if log_x and "log10" in x_label else tick
            fig.text(px, y + h + 13, _tick_label(display_tick), 8, "muted", "middle")
    fig.text(x + w / 2, y + h + 27, x_label, 8, "ink", "middle")
    fig.text(max(16, x - 56), y + h / 2, y_label, 8, "ink", "middle", vertical=True)


def _panel_points(
    panel: tuple[float, float, float, float],
    x_values: Sequence[float],
    y_values: Sequence[float],
    *,
    log_x: bool,
    log_y: bool,
    y_domain: tuple[float, float] | None = None,
) -> list[tuple[float, float]]:
    x_domain = _domain(x_values, log_x)
    yd = y_domain if y_domain is not None else _domain(y_values, log_y)
    return [
        (_scale_x(panel, x, x_domain, log_x), _scale_y(panel, y, yd, log_y))
        for x, y in zip(x_values, y_values)
        if _valid_for_scale(x, log_x) and _valid_for_scale(y, log_y)
    ]


def _stacked_panels(count: int) -> list[tuple[float, float, float, float]]:
    left, top, width = 76.0, 68.0, 386.0
    height = (260.0 - PANEL_GAP * (count - 1)) / count
    return [(left, top + i * (height + PANEL_GAP), width, height) for i in range(count)]


def _grid_panels(rows: int, cols: int) -> list[tuple[float, float, float, float]]:
    left, top = 64.0, 66.0
    column_gap, row_gap = 64.0, 54.0
    width = (430.0 - column_gap * (cols - 1)) / cols
    height = (266.0 - row_gap * (rows - 1)) / rows
    panels = []
    for r in range(rows):
        for c in range(cols):
            panels.append((left + c * (width + column_gap), top + r * (height + row_gap), width, height))
    return panels


def _scale_x(panel: tuple[float, float, float, float], value: float, domain: tuple[float, float], log_scale: bool) -> float:
    x, _y, w, _h = panel
    v = math.log10(value) if log_scale else value
    lo = math.log10(domain[0]) if log_scale else domain[0]
    hi = math.log10(domain[1]) if log_scale else domain[1]
    return x + (v - lo) / (hi - lo) * w


def _scale_y(panel: tuple[float, float, float, float], value: float, domain: tuple[float, float], log_scale: bool) -> float:
    _x, y, _w, h = panel
    v = math.log10(value) if log_scale else value
    lo = math.log10(domain[0]) if log_scale else domain[0]
    hi = math.log10(domain[1]) if log_scale else domain[1]
    return y + h - (v - lo) / (hi - lo) * h


def _domain(values: Iterable[float], log_scale: bool) -> tuple[float, float]:
    selected = [float(value) for value in values if math.isfinite(float(value)) and (float(value) > 0 if log_scale else True)]
    if not selected:
        return (1.0e-3, 1.0) if log_scale else (0.0, 1.0)
    lo, hi = min(selected), max(selected)
    if lo == hi:
        delta = abs(lo) * 0.2 if lo else 1.0
        lo, hi = lo - delta, hi + delta
    if log_scale:
        log_lo, log_hi = math.log10(lo), math.log10(hi)
        pad = 0.08 * (log_hi - log_lo if log_hi > log_lo else 1.0)
        return 10 ** (log_lo - pad), 10 ** (log_hi + pad)
    pad = 0.08 * (hi - lo)
    return lo - pad, hi + pad


def _ticks(domain: tuple[float, float], log_scale: bool) -> list[float]:
    lo, hi = domain
    if log_scale:
        start = math.ceil(math.log10(lo))
        end = math.floor(math.log10(hi))
        ticks = [10.0**power for power in range(start, end + 1)]
        if len(ticks) < 3:
            ticks = [10 ** (math.log10(lo) + i * (math.log10(hi) - math.log10(lo)) / 3.0) for i in range(4)]
        return ticks[:5]
    return [lo + i * (hi - lo) / 4.0 for i in range(5)]


def _tick_label(value: float) -> str:
    if value == 0:
        return "0"
    return f"{value:.4g}"


def _legend(fig: Figure, panel: tuple[float, float, float, float], items: Sequence[tuple[str, str]]) -> None:
    x = panel[0] + panel[2] - 92
    y = panel[1] + 12
    for index, (label, color) in enumerate(items):
        yy = y + index * 12
        fig.line(x, yy - 3, x + 12, yy - 3, color, 1.4)
        fig.text(x + 16, yy, label, 6.6, "muted")


def _sort_by_a(rows: Sequence[Mapping[str, str]]) -> list[Mapping[str, str]]:
    return sorted(rows, key=lambda row: _float(row, "a"))


def _aggregate_summary_by_a(rows: Sequence[Mapping[str, str]]) -> list[dict[str, float]]:
    grouped: dict[tuple[float, float], list[Mapping[str, str]]] = {}
    for row in rows:
        z = _float(row, "z")
        grouped.setdefault((_scale_factor(row), z), []).append(row)
    output = []
    for (a, z), selected in sorted(grouped.items()):
        total = sum(int(float(row["object_count"])) for row in selected)
        if total <= 0:
            continue
        output.append(
            {
                "a": a,
                "z": z,
                "forming_fraction": sum(int(float(row["object_count"])) * _float(row, "forming_fraction") for row in selected) / total,
                "balanced_fraction": sum(int(float(row["object_count"])) * _float(row, "balanced_fraction") for row in selected) / total,
                "stripping_fraction": sum(int(float(row["object_count"])) * _float(row, "stripping_fraction") for row in selected) / total,
            }
        )
    return output


def _aggregate_transfer_by_a(rows: Sequence[Mapping[str, str]]) -> list[dict[str, float]]:
    grouped: dict[tuple[float, float], list[Mapping[str, str]]] = {}
    for row in rows:
        z = _float(row, "z")
        grouped.setdefault((_scale_factor(row), z), []).append(row)
    output = []
    for (a, z), selected in sorted(grouped.items()):
        output.append(
            {
                "a": a,
                "z": z,
                "a_pq_median": statistics.median(_values(selected, "a_pq")),
                "abs_g_pq_median": statistics.median(abs(_float(row, "g_pq")) for row in selected),
            }
        )
    return output


def _table_sizes(run_dir: Path) -> list[tuple[str, float]]:
    rows = []
    for path in (run_dir / "tables").glob("*.csv"):
        rows.append((path.name, float(path.stat().st_size)))
    return sorted(rows, key=lambda item: item[1], reverse=True)


def _values(rows: Sequence[Mapping[str, str]], key: str) -> list[float]:
    values = []
    for row in rows:
        try:
            value = float(row[key])
        except (KeyError, TypeError, ValueError):
            continue
        if math.isfinite(value):
            values.append(value)
    return values


def _float(row: Mapping[str, str], key: str) -> float:
    return float(row[key])


def _scale_factor(row: Mapping[str, str]) -> float:
    if "a" in row and row["a"] != "":
        return float(row["a"])
    return 1.0 / (1.0 + float(row["z"]))


def _nearest(rows: Sequence[Mapping[str, str]], key: str, value: float) -> Mapping[str, str]:
    return min(rows, key=lambda row: abs(_float(row, key) - value))


def _first(rows: Sequence[Mapping[str, str]] | Sequence[Mapping[str, float]], predicate) -> Mapping[str, str] | Mapping[str, float] | None:
    for row in rows:
        if predicate(row):
            return row
    return None


def _rescale(values: Sequence[float], target: Sequence[float]) -> list[float]:
    v_domain = _domain(values, False)
    t_domain = _domain(target, False)
    if v_domain[0] == v_domain[1]:
        return [0.5 * (t_domain[0] + t_domain[1]) for _ in values]
    return [t_domain[0] + (value - v_domain[0]) / (v_domain[1] - v_domain[0]) * (t_domain[1] - t_domain[0]) for value in values]


def _valid_for_scale(value: float, log_scale: bool) -> bool:
    return math.isfinite(value) and (value > 0.0 if log_scale else True)


def _pdf_color(value: str):
    if colors is None:
        return None
    value = value.lstrip("#")
    red = int(value[0:2], 16) / 255.0
    green = int(value[2:4], 16) / 255.0
    blue = int(value[4:6], 16) / 255.0
    return colors.Color(red, green, blue)


if __name__ == "__main__":
    raise SystemExit(main())
