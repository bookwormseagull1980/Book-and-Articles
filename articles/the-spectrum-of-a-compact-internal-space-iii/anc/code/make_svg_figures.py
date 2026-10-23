"""Create vector figures from saved figure-source tables."""

from __future__ import annotations

import argparse
from html import escape
import math
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from cosmic_cg.intermediate import write_output_inventory
from cosmic_cg.io.tables import read_rows


WIDTH = 1600
HEIGHT = 1100
MARGIN = {"left": 170, "right": 80, "top": 95, "bottom": 150}
COLORS = {
    "ink": "#151515",
    "muted": "#676767",
    "grid": "#d9d9d9",
    "endpoint": "#1b6fb8",
    "tail": "#b4472f",
    "entropy": "#2f7d4f",
    "forming": "#b4472f",
    "balanced": "#2f7d4f",
    "stripping": "#4b5f9e",
}


def main() -> int:
    args = parse_args()
    run_dir = Path(args.run_dir)
    out_dir = run_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)

    written = []
    figure2 = run_dir / "figure_data" / "figure2_moving_endpoint.csv"
    figure3 = run_dir / "figure_data" / "figure3_tail_activity.csv"
    figure4 = run_dir / "figure_data" / "figure4_object_fate.csv"
    if figure2.is_file():
        path = out_dir / "figure2_moving_endpoint.svg"
        path.write_text(_figure2(read_rows(figure2), run_dir.name), encoding="utf-8")
        written.append(path)
    if figure3.is_file():
        path = out_dir / "figure3_tail_activity.svg"
        path.write_text(_figure3(read_rows(figure3), run_dir.name), encoding="utf-8")
        written.append(path)
    if figure4.is_file():
        path = out_dir / "figure4_object_fate.svg"
        path.write_text(_figure4(read_rows(figure4), run_dir.name), encoding="utf-8")
        written.append(path)

    for path in written:
        print(path)
    if written:
        write_output_inventory(run_dir)
    if args.strict and not written:
        return 1
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--strict", action="store_true")
    return parser.parse_args()


def _figure2(rows: list[dict[str, str]], run_id: str) -> str:
    z_values = [_num(row, "z") for row in rows]
    k_values = [_num(row, "k_NL_hmpc") for row in rows]
    s_values = [_num(row, "S_AW") for row in rows]
    x_min, x_max = _range(z_values)
    y_min, y_max = _range(k_values + [max(value, 1.0e-8) for value in k_values])
    y_min = max(y_min * 0.75, 1.0e-8)
    y_max = y_max * 1.25
    canvas = _canvas("Moving Nonlinear Endpoint", run_id, "z", "k_NL [h/Mpc]")
    canvas.extend(_axes(x_min, x_max, y_min, y_max, log_y=True))
    points = [(_x(z, x_min, x_max), _y_log(k, y_min, y_max)) for z, k in zip(z_values, k_values)]
    canvas.append(_polyline(points, COLORS["endpoint"]))
    for (x, y), row in zip(points, rows):
        canvas.append(_circle(x, y, 9, COLORS["endpoint"]))
        canvas.append(_label(x + 18, y - 15, f"S_AW={_num(row, 'S_AW'):.3g}", 26, COLORS["muted"]))
    if s_values:
        canvas.append(_label(WIDTH - 430, 150, f"active window: {min(s_values):.3g}-{max(s_values):.3g}", 28, COLORS["entropy"]))
    return _finish(canvas)


def _figure3(rows: list[dict[str, str]], run_id: str) -> str:
    canvas = _canvas("Tail Activity Audit", run_id, "diagnostic", "value")
    labels = ["R_h", "Omega_tail"]
    values = [
        _num(rows[0], "R_h_median") if rows else 0.0,
        _num(rows[0], "Omega_tail") if rows else 0.0,
    ]
    lows = [_num(rows[0], "R_h_p16") if rows else 0.0, values[1]]
    highs = [_num(rows[0], "R_h_p84") if rows else 0.0, values[1]]
    y_min, y_max = 0.0, max(max(highs) * 1.25, 1.0)
    plot_left = MARGIN["left"]
    plot_right = WIDTH - MARGIN["right"]
    plot_bottom = HEIGHT - MARGIN["bottom"]
    bar_width = 170
    canvas.extend(_axes(0.0, 1.0, y_min, y_max, categorical=True))
    for i, label in enumerate(labels):
        x = plot_left + (i + 1) * (plot_right - plot_left) / 3.0
        y = _y(values[i], y_min, y_max)
        bottom = plot_bottom
        color = COLORS["tail"] if i == 0 else COLORS["endpoint"]
        canvas.append(_rect(x - bar_width / 2, y, bar_width, bottom - y, color))
        canvas.append(_label(x - 45, plot_bottom + 55, label, 30, COLORS["ink"]))
        canvas.append(_label(x - 50, y - 20, f"{values[i]:.3g}", 27, COLORS["ink"]))
        if lows[i] != highs[i]:
            y_low = _y(lows[i], y_min, y_max)
            y_high = _y(highs[i], y_min, y_max)
            canvas.append(_line(x, y_low, x, y_high, COLORS["ink"], 4))
            canvas.append(_line(x - 42, y_low, x + 42, y_low, COLORS["ink"], 4))
            canvas.append(_line(x - 42, y_high, x + 42, y_high, COLORS["ink"], 4))
    return _finish(canvas)


def _figure4(rows: list[dict[str, str]], run_id: str) -> str:
    masses = [math.log10(_num(row, "mass_msun_h")) for row in rows if _positive(row, "mass_msun_h")]
    b_values = [_num(row, "B_h") for row in rows if _finite(row, "B_h")]
    x_min, x_max = _range(masses)
    y_min, y_max = _range(b_values + [0.0, 1.0])
    y_min -= 0.15 * (y_max - y_min)
    y_max += 0.15 * (y_max - y_min)
    canvas = _canvas("Object Fate Diagram", run_id, "log10 M [Msun/h]", "B_h")
    canvas.extend(_axes(x_min, x_max, y_min, y_max))
    canvas.append(_line(MARGIN["left"], _y(0.0, y_min, y_max), WIDTH - MARGIN["right"], _y(0.0, y_min, y_max), COLORS["grid"], 3))
    canvas.append(_line(MARGIN["left"], _y(1.0, y_min, y_max), WIDTH - MARGIN["right"], _y(1.0, y_min, y_max), COLORS["grid"], 3))
    for row in rows:
        if not _positive(row, "mass_msun_h") or not _finite(row, "B_h"):
            continue
        x = _x(math.log10(_num(row, "mass_msun_h")), x_min, x_max)
        y = _y(_num(row, "B_h"), y_min, y_max)
        klass = row.get("B_class", "balanced")
        canvas.append(_circle(x, y, 10, COLORS.get(klass, COLORS["muted"])))
    canvas.extend(_legend([("forming", COLORS["forming"]), ("balanced", COLORS["balanced"]), ("stripping", COLORS["stripping"])]))
    return _finish(canvas)


def _canvas(title: str, run_id: str, x_label: str, y_label: str) -> list[str]:
    return [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{HEIGHT}" viewBox="0 0 {WIDTH} {HEIGHT}">',
        '<rect width="100%" height="100%" fill="white"/>',
        _label(MARGIN["left"], 60, title, 44, COLORS["ink"]),
        _label(WIDTH - MARGIN["right"] - 420, 60, run_id, 24, COLORS["muted"]),
        _label(WIDTH / 2 - 80, HEIGHT - 45, x_label, 32, COLORS["ink"]),
        f'<text x="55" y="{HEIGHT / 2 + 95:.1f}" font-family="Arial" font-size="32" fill="{COLORS["ink"]}" transform="rotate(-90 55,{HEIGHT / 2 + 95:.1f})">{escape(y_label)}</text>',
    ]


def _axes(x_min: float, x_max: float, y_min: float, y_max: float, *, log_y: bool = False, categorical: bool = False) -> list[str]:
    left = MARGIN["left"]
    right = WIDTH - MARGIN["right"]
    top = MARGIN["top"]
    bottom = HEIGHT - MARGIN["bottom"]
    lines = [
        _line(left, bottom, right, bottom, COLORS["ink"], 3),
        _line(left, top, left, bottom, COLORS["ink"], 3),
    ]
    for i in range(6):
        value = y_min + (y_max - y_min) * i / 5.0
        y = _y_log(value, y_min, y_max) if log_y else _y(value, y_min, y_max)
        lines.append(_line(left, y, right, y, COLORS["grid"], 1))
        label = f"{value:.2g}" if log_y else f"{value:.2f}"
        lines.append(_label(left - 95, y + 10, label, 24, COLORS["muted"]))
    if not categorical:
        for i in range(5):
            value = x_min + (x_max - x_min) * i / 4.0
            x = _x(value, x_min, x_max)
            lines.append(_line(x, bottom, x, bottom + 12, COLORS["ink"], 3))
            lines.append(_label(x - 25, bottom + 45, f"{value:.2g}", 24, COLORS["muted"]))
    return lines


def _legend(items: list[tuple[str, str]]) -> list[str]:
    y0 = 145
    lines = []
    for i, (label, color) in enumerate(items):
        y = y0 + i * 42
        lines.append(_circle(WIDTH - 310, y - 8, 9, color))
        lines.append(_label(WIDTH - 285, y, label, 25, COLORS["ink"]))
    return lines


def _range(values: list[float]) -> tuple[float, float]:
    finite = [value for value in values if math.isfinite(value)]
    if not finite:
        return 0.0, 1.0
    lo = min(finite)
    hi = max(finite)
    if lo == hi:
        delta = abs(lo) * 0.2 if lo != 0 else 1.0
        return lo - delta, hi + delta
    return lo, hi


def _x(value: float, x_min: float, x_max: float) -> float:
    left = MARGIN["left"]
    right = WIDTH - MARGIN["right"]
    return left + (value - x_min) / (x_max - x_min) * (right - left)


def _y(value: float, y_min: float, y_max: float) -> float:
    top = MARGIN["top"]
    bottom = HEIGHT - MARGIN["bottom"]
    return bottom - (value - y_min) / (y_max - y_min) * (bottom - top)


def _y_log(value: float, y_min: float, y_max: float) -> float:
    return _y(math.log10(value), math.log10(y_min), math.log10(y_max))


def _label(x: float, y: float, text: str, size: int, color: str) -> str:
    return f'<text x="{x:.1f}" y="{y:.1f}" font-family="Arial" font-size="{size}" fill="{color}">{escape(text)}</text>'


def _line(x1: float, y1: float, x2: float, y2: float, color: str, width: int) -> str:
    return f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{color}" stroke-width="{width}"/>'


def _polyline(points: list[tuple[float, float]], color: str) -> str:
    if not points:
        return ""
    coords = " ".join(f"{x:.1f},{y:.1f}" for x, y in points)
    return f'<polyline points="{coords}" fill="none" stroke="{color}" stroke-width="5" stroke-linejoin="round"/>'


def _circle(x: float, y: float, radius: float, color: str) -> str:
    return f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{radius:.1f}" fill="{color}"/>'


def _rect(x: float, y: float, width: float, height: float, color: str) -> str:
    return f'<rect x="{x:.1f}" y="{y:.1f}" width="{width:.1f}" height="{height:.1f}" fill="{color}"/>'


def _num(row: dict[str, str], key: str) -> float:
    return float(row.get(key, "nan") or "nan")


def _finite(row: dict[str, str], key: str) -> bool:
    try:
        return math.isfinite(_num(row, key))
    except ValueError:
        return False


def _positive(row: dict[str, str], key: str) -> bool:
    return _finite(row, key) and _num(row, key) > 0


def _finish(lines: list[str]) -> str:
    return "\n".join(line for line in lines if line) + "\n</svg>\n"


if __name__ == "__main__":
    raise SystemExit(main())
