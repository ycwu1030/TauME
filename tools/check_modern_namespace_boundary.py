#!/usr/bin/env python3
"""Reject dependencies from the modern TauOO path onto legacy compatibility code."""

from __future__ import annotations

from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
PRODUCTION_ROOTS = (ROOT / "include/tauamp/tautau", ROOT / "include/tauamp/delphes")
TEST_PATTERNS = (
    "tests/tautau/rho_a1*_test.cpp",
    "tests/tautau/ordered_channel*_test.cpp",
    "tests/tautau/common_channel*_test.cpp",
    "tests/delphes/common_channel*_test.cpp",
)
ALLOWED_INTERNAL_PREFIXES = ("tauamp/tautau/", "tauamp/delphes/")
FORBIDDEN_HEADERS = {
    "tauamp/constants.h",
    "tauamp/decayplane.h",
    "tauamp/eventreader.h",
    "tauamp/matrixelements.h",
    "tauamp/tauamp.h",
    "tauamp/taudecay.h",
    "tauamp/utilities.h",
}
FORBIDDEN_SYMBOLS = (
    "tauamp::TauDecay_",
    "tauamp::ME_CPV_",
    "tauamp::EventReader",
    "dlv_t",
    "clv_t",
    "cd_t",
)
INCLUDE_RE = re.compile(r'^\s*#\s*include\s+"([^"]+)"', re.MULTILINE)


def files_to_check() -> list[Path]:
    files: set[Path] = set()
    for root in PRODUCTION_ROOTS:
        files.update(root.rglob("*.h"))
        files.update(root.rglob("*.cpp"))
    for pattern in TEST_PATTERNS:
        files.update(ROOT.glob(pattern))
    return sorted(files)


def main() -> int:
    errors: list[str] = []
    checked = files_to_check()
    for path in checked:
        text = path.read_text(encoding="utf-8")
        relative = path.relative_to(ROOT)
        for header in INCLUDE_RE.findall(text):
            if header in FORBIDDEN_HEADERS:
                errors.append(f"{relative}: forbidden legacy header {header}")
            if header.startswith("tauamp/") and not header.startswith(ALLOWED_INTERNAL_PREFIXES):
                errors.append(f"{relative}: internal include is outside modern namespace boundary: {header}")
        for symbol in FORBIDDEN_SYMBOLS:
            if symbol in text:
                errors.append(f"{relative}: forbidden legacy symbol {symbol}")
    if errors:
        print("modern TauOO namespace boundary: FAIL", file=sys.stderr)
        for error in errors:
            print(error, file=sys.stderr)
        return 1
    print(f"modern TauOO namespace boundary: PASS ({len(checked)} files checked)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
