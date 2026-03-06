"""
Journal Impact Factor lookup and paper filtering.
Loads data/journal_if.csv (format: journal_title,impact_factor) and filters
records to keep only papers from journals with IF >= min_if.
When exclude_no_if=True, papers with no IF data are also excluded.
"""
from __future__ import annotations

import csv
from pathlib import Path

_DATA_DIR = Path(__file__).resolve().parent / "data"
_JOURNAL_IF_PATH = _DATA_DIR / "journal_if.csv"
_MIN_IF_DEFAULT = 10.0


def _normalize_journal(s: str) -> str:
    if not s:
        return ""
    return " ".join(str(s).strip().split()).lower()


def load_journal_if_map() -> dict[str, float]:
    """Load journal -> impact factor from data/journal_if.csv. Keys are normalized journal names."""
    out = {}
    if not _JOURNAL_IF_PATH.exists():
        return out
    try:
        with open(_JOURNAL_IF_PATH, encoding="utf-8") as f:
            for row in csv.reader(f):
                if not row or row[0].strip().startswith("#"):
                    continue
                name = _normalize_journal(row[0])
                if not name:
                    continue
                try:
                    if_val = float(row[1].strip().replace(",", "."))
                    out[name] = if_val
                except (ValueError, IndexError):
                    continue
    except Exception:
        pass
    return out


def filter_records_by_min_if(
    records: list[dict],
    min_if: float = _MIN_IF_DEFAULT,
    journal_if_map: dict[str, float] | None = None,
    exclude_no_if: bool = True,
) -> list[dict]:
    """
    Keep only records where journal is in the IF map and IF >= min_if.
    When exclude_no_if=True (default), also exclude papers with no IF data (journal not in list or empty).
    When exclude_no_if=False, keep papers with unknown IF; only remove when known IF < min_if.
    Each record must have a "journal" key (string, can be empty).
    """
    if journal_if_map is None:
        journal_if_map = load_journal_if_map()
    if not journal_if_map:
        return records
    kept = []
    for r in records:
        journal = _normalize_journal(r.get("journal") or "")
        if not journal:
            if not exclude_no_if:
                kept.append(r)
            continue
        if journal not in journal_if_map:
            if not exclude_no_if:
                kept.append(r)
            continue
        if journal_if_map[journal] >= min_if:
            kept.append(r)
    return kept
