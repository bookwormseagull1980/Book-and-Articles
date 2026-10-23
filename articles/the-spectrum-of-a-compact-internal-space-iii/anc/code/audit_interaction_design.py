"""Audit the interaction layer used by generated-universe runs."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import json
from pathlib import Path
import re
import sys
from typing import Iterable, Mapping


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from cosmic_cg.io.audit import audit_run
from cosmic_cg.io.field_audit import audit_field_run
from cosmic_cg.io.pm_audit import audit_pm_run


@dataclass(frozen=True)
class SourceRule:
    pattern: str
    category: str
    reason: str


def _token(*parts: str) -> str:
    return "".join(parts)


SOURCE_RULES = (
    SourceRule(rf"\b{_token('peak', '_z')}\b", "manual_epoch", "interaction channels must be selected from stored states"),
    SourceRule(rf"\b{_token('width', '_ln', '_a')}\b", "manual_epoch", "interaction channels must be selected from stored states"),
    SourceRule(rf"\b{_token('epoch', '_phase')}\b", "manual_epoch", "interaction phases must be derived from stored states"),
    SourceRule(rf"\b{_token('applicable', '_epoch')}\b", "manual_epoch", "interaction strengths determine active channels"),
    SourceRule(rf"\b{_token('exit', '_condition')}\b", "manual_exit", "channels are evaluated as continuous records"),
    SourceRule(rf"\b{_token('handoff', '_from')}\b", "manual_handoff", "handoffs must be inferred from dominance changes"),
    SourceRule(rf"\b{_token('handoff', '_to')}\b", "manual_handoff", "handoffs must be inferred from dominance changes"),
    SourceRule(_token("适用", "时期"), "manual_epoch", "Chinese manual-period marker"),
    SourceRule(_token("退出", "条件"), "manual_exit", "Chinese exit-condition marker"),
    SourceRule(rf"\b{_token('TNG', '_API', '_KEY')}\b", "external_data", "the generator must not read external survey credentials"),
    SourceRule(r"\billustris\b", "external_data", "the generator must not depend on external simulations"),
    SourceRule(r"\btng\b", "external_data", "the generator must not depend on external simulations"),
    SourceRule(r"\babacus\b", "external_data", "the generator must not depend on external simulations"),
    SourceRule(r"\bquijote\b", "external_data", "the generator must not depend on external simulations"),
    SourceRule(r"\brequests\b", "external_data", "the generator must not call remote data services"),
    SourceRule(r"\burllib\b", "external_data", "the generator must not call remote data services"),
)

SOURCE_DIRS = ("cosmic_cg", "scripts")
SKIP_SOURCE_FILES = {Path("scripts") / "audit_interaction_design.py"}


def main() -> int:
    args = parse_args()
    root = Path(args.source_root).resolve()
    run_dir = Path(args.run_dir)
    if not run_dir.is_absolute():
        run_dir = root / run_dir

    source_matches = scan_sources(root)
    run_summary = audit_interaction_run(run_dir)
    issues = []
    issues.extend(source_matches)
    issues.extend(run_summary["issues"])
    status = "pass" if not any(item["severity"] == "error" for item in issues) else "fail"

    output = {
        "status": status,
        "source_root": str(root),
        "run_dir": str(run_dir),
        "source_scan_match_count": len(source_matches),
        "table_summary": run_summary["tables"],
        "coverage": run_summary["coverage"],
        "issues": issues,
    }
    write_reports(output, run_dir / "quality_checks")
    print(json.dumps(output, ensure_ascii=False, indent=2))
    if args.strict and status != "pass":
        return 1
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, help="generated-universe result directory")
    parser.add_argument("--source-root", default=PROJECT_ROOT, help="project root to scan")
    parser.add_argument("--strict", action="store_true", help="return nonzero when an error is found")
    return parser.parse_args()


def scan_sources(root: Path) -> list[dict[str, object]]:
    matches = []
    compiled = [(rule, re.compile(rule.pattern, re.IGNORECASE)) for rule in SOURCE_RULES]
    for directory in SOURCE_DIRS:
        base = root / directory
        if not base.is_dir():
            continue
        for path in sorted(base.rglob("*.py")):
            relative = path.relative_to(root)
            if relative in SKIP_SOURCE_FILES or "__pycache__" in relative.parts:
                continue
            text = path.read_text(encoding="utf-8", errors="ignore")
            for line_no, line in enumerate(text.splitlines(), 1):
                for rule, pattern in compiled:
                    if pattern.search(line):
                        matches.append(
                            {
                                "severity": "error",
                                "category": rule.category,
                                "file": str(relative),
                                "line": line_no,
                                "reason": rule.reason,
                            }
                        )
    return matches


def audit_interaction_run(run_dir: Path) -> dict[str, object]:
    issues: list[dict[str, object]] = []
    builtin = audit_base_run(run_dir)
    if builtin.status != "pass":
        issues.append(
            {
                "severity": "error",
                "category": "run_audit",
                "file": str(run_dir),
                "reason": "generated-universe run audit did not pass",
            }
        )

    tables = {}
    strength_rows = read_rows(run_dir / "tables" / "interaction_strength_table.csv", issues)
    dominance_rows = read_rows(run_dir / "tables" / "interaction_dominance_table.csv", issues)
    budget_rows = read_rows(run_dir / "tables" / "interaction_budget_table.csv", issues)
    crossing_rows = read_rows(run_dir / "tables" / "interaction_crossing_table.csv", issues)
    carrier_rows = read_rows(run_dir / "tables" / "carrier_history_table.csv", issues)
    event_rows = read_rows(run_dir / "tables" / "emergence_event_table.csv", issues)
    table_rows = {
        "interaction_strength_table.csv": strength_rows,
        "interaction_dominance_table.csv": dominance_rows,
        "interaction_budget_table.csv": budget_rows,
        "interaction_crossing_table.csv": crossing_rows,
        "carrier_history_table.csv": carrier_rows,
        "emergence_event_table.csv": event_rows,
    }
    for name, rows in table_rows.items():
        tables[name] = {"rows": len(rows), "columns": list(rows[0].keys()) if rows else []}

    check_strength_table(strength_rows, issues)
    check_dominance_table(strength_rows, dominance_rows, issues)
    check_table_headers(table_rows, issues)

    coverage = coverage_summary(strength_rows, budget_rows, crossing_rows, carrier_rows, event_rows)
    if coverage["channel_count"] < 10:
        issues.append(error("coverage", "interaction_strength_table.csv", "fewer than ten interaction channels"))
    if coverage["band_count"] < 5:
        issues.append(error("coverage", "interaction_strength_table.csv", "fewer than five scale bands"))
    if coverage["z_max"] < 1.0e10:
        issues.append(error("coverage", "interaction_strength_table.csv", "early high-redshift states are missing"))
    if coverage["z_min"] > 1.0e-12:
        issues.append(error("coverage", "interaction_strength_table.csv", "present or future states are missing"))
    if coverage["crossing_count"] == 0:
        issues.append(error("coverage", "interaction_crossing_table.csv", "no dominance changes were recorded"))
    if coverage["event_count"] == 0:
        issues.append(error("coverage", "emergence_event_table.csv", "no emergence events were recorded"))

    return {"issues": issues, "tables": tables, "coverage": coverage}


def audit_base_run(run_dir: Path):
    tables = run_dir / "tables"
    if (tables / "field_window_table.csv").is_file():
        return audit_field_run(run_dir)
    if (tables / "pm_snapshot_table.csv").is_file():
        return audit_pm_run(run_dir)
    return audit_run(run_dir, "generated_universe")


def read_rows(path: Path, issues: list[dict[str, object]]) -> list[dict[str, str]]:
    if not path.is_file():
        issues.append(error("missing_table", str(path), "required interaction table is missing"))
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def check_strength_table(rows: Iterable[Mapping[str, str]], issues: list[dict[str, object]]) -> None:
    grouped: dict[tuple[str, str, str], list[Mapping[str, str]]] = {}
    for row in rows:
        key = (row["z"], row["a"], row["band_id"])
        grouped.setdefault(key, []).append(row)
        for column in ("carrier", "mediator", "source_term", "sink_term", "output_table"):
            if row[column] == "":
                issues.append(error("blank_field", "interaction_strength_table.csv", f"{column} is blank"))
        if float(row["raw_strength"]) < 0.0:
            issues.append(error("strength_range", "interaction_strength_table.csv", "negative raw strength"))
        chi = float(row["chi"])
        if chi < 0.0 or chi > 1.0:
            issues.append(error("strength_range", "interaction_strength_table.csv", "chi lies outside [0, 1]"))
    for key, selected in grouped.items():
        chi_sum = sum(float(row["chi"]) for row in selected)
        if abs(chi_sum - 1.0) > 1.0e-8:
            issues.append(error("normalization", "interaction_strength_table.csv", f"chi sum is {chi_sum} at {key}"))


def check_dominance_table(
    strength_rows: Iterable[Mapping[str, str]],
    dominance_rows: Iterable[Mapping[str, str]],
    issues: list[dict[str, object]],
) -> None:
    grouped: dict[tuple[str, str, str], list[Mapping[str, str]]] = {}
    for row in strength_rows:
        grouped.setdefault((row["z"], row["a"], row["band_id"]), []).append(row)
    dominance = {(row["z"], row["a"], row["band_id"]): row for row in dominance_rows}
    for key, selected in grouped.items():
        if key not in dominance:
            issues.append(error("dominance", "interaction_dominance_table.csv", f"missing dominance row at {key}"))
            continue
        strongest = max(selected, key=lambda row: float(row["chi"]))
        if dominance[key]["dominant_channel"] != strongest["channel_id"]:
            issues.append(error("dominance", "interaction_dominance_table.csv", f"dominant channel mismatch at {key}"))


def check_table_headers(table_rows: Mapping[str, list[Mapping[str, str]]], issues: list[dict[str, object]]) -> None:
    blocked_fragments = (
        _token("peak", "_z"),
        _token("width", "_ln", "_a"),
        _token("epoch", "_phase"),
        _token("applicable", "_epoch"),
        _token("exit", "_condition"),
        _token("handoff", "_from"),
        _token("handoff", "_to"),
    )
    for name, rows in table_rows.items():
        if not rows:
            continue
        for column in rows[0].keys():
            if any(fragment in column for fragment in blocked_fragments):
                issues.append(error("manual_column", name, f"manual scheduling column {column} is present"))


def coverage_summary(
    strength_rows: list[Mapping[str, str]],
    budget_rows: list[Mapping[str, str]],
    crossing_rows: list[Mapping[str, str]],
    carrier_rows: list[Mapping[str, str]],
    event_rows: list[Mapping[str, str]],
) -> dict[str, object]:
    z_values = [float(row["z"]) for row in strength_rows] or [0.0]
    return {
        "state_count": len({(row["z"], row["a"]) for row in strength_rows}),
        "band_count": len({row["band_id"] for row in strength_rows}),
        "channel_count": len({row["channel_id"] for row in strength_rows}),
        "carrier_count": len({row["carrier"] for row in carrier_rows}),
        "budget_row_count": len(budget_rows),
        "crossing_count": len(crossing_rows),
        "event_count": len(event_rows),
        "z_min": min(z_values),
        "z_max": max(z_values),
    }


def error(category: str, file_name: str, reason: str) -> dict[str, object]:
    return {"severity": "error", "category": category, "file": file_name, "reason": reason}


def write_reports(report: Mapping[str, object], qc_dir: Path) -> None:
    qc_dir.mkdir(parents=True, exist_ok=True)
    rendered = json.dumps(report, ensure_ascii=False, indent=2)
    (qc_dir / "interaction_design_audit.json").write_text(rendered + "\n", encoding="utf-8")
    coverage = report["coverage"]
    lines = [
        "# Interaction Design Audit",
        "",
        f"status: {report['status']}",
        f"source scan matches: {report['source_scan_match_count']}",
        f"states: {coverage['state_count']}",
        f"scale bands: {coverage['band_count']}",
        f"channels: {coverage['channel_count']}",
        f"carriers: {coverage['carrier_count']}",
        f"crossings: {coverage['crossing_count']}",
        f"events: {coverage['event_count']}",
        f"z range: {coverage['z_min']} to {coverage['z_max']}",
        "",
        "Issues",
        "",
    ]
    issues = report["issues"]
    if issues:
        lines.extend(f"- {item['severity']} {item['category']} {item['file']} {item['reason']}" for item in issues)
    else:
        lines.append("- none")
    (qc_dir / "interaction_design_audit.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
