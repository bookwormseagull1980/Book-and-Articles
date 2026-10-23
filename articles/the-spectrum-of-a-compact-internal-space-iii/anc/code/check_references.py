"""Check the local literature register for duplicate keys and missing fields."""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REFERENCES = ROOT / "references"
BIB_PATH = REFERENCES / "references.bib"
MATRIX_PATH = REFERENCES / "literature_matrix.csv"


ENTRY_RE = re.compile(r"@(?P<kind>[A-Za-z]+)\s*\{\s*(?P<key>[^,\s]+)\s*,", re.MULTILINE)
FIELD_RE = re.compile(r"^\s*(?P<field>[A-Za-z][A-Za-z0-9_-]*)\s*=", re.MULTILINE)


def _entry_spans(text: str) -> list[tuple[str, str]]:
    matches = list(ENTRY_RE.finditer(text))
    spans: list[tuple[str, str]] = []
    for index, match in enumerate(matches):
        start = match.start()
        end = matches[index + 1].start() if index + 1 < len(matches) else len(text)
        spans.append((match.group("key"), text[start:end]))
    return spans


def _fields(block: str) -> set[str]:
    return {match.group("field").lower() for match in FIELD_RE.finditer(block)}


def read_bib(path: Path = BIB_PATH) -> dict[str, set[str]]:
    text = path.read_text(encoding="utf-8")
    entries: dict[str, set[str]] = {}
    duplicates: set[str] = set()
    for key, block in _entry_spans(text):
        if key in entries:
            duplicates.add(key)
        entries[key] = _fields(block)
    if duplicates:
        entries["__DUPLICATES__"] = duplicates
    return entries


def read_matrix_keys(path: Path = MATRIX_PATH) -> set[str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return {row["key"].strip() for row in reader if row.get("key", "").strip()}


def check() -> dict[str, object]:
    entries = read_bib()
    duplicates = sorted(entries.pop("__DUPLICATES__", set()))
    matrix_keys = read_matrix_keys()
    bib_keys = set(entries)

    required = {"author", "title", "year"}
    missing_required = {
        key: sorted(required - fields)
        for key, fields in entries.items()
        if required - fields
    }
    missing_locator = {
        key: sorted({"doi", "eprint", "url"} - fields)
        for key, fields in entries.items()
        if not ({"doi", "eprint", "url"} & fields)
    }

    report = {
        "status": "ok",
        "bib_entries": len(entries),
        "matrix_entries": len(matrix_keys),
        "duplicate_keys": duplicates,
        "matrix_missing_bib": sorted(matrix_keys - bib_keys),
        "bib_missing_matrix": sorted(bib_keys - matrix_keys),
        "missing_required_fields": missing_required,
        "missing_locator_fields": missing_locator,
    }

    if (
        duplicates
        or report["matrix_missing_bib"]
        or report["bib_missing_matrix"]
        or missing_required
    ):
        report["status"] = "fail"
    elif missing_locator:
        report["status"] = "warn"
    return report


def main() -> int:
    report = check()
    print(json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True))
    return 0 if report["status"] in {"ok", "warn"} else 1


if __name__ == "__main__":
    raise SystemExit(main())

