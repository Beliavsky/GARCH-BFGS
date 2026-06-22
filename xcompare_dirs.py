#!/usr/bin/env python3
"""Compare two directories by selected file types."""

from __future__ import annotations

import argparse
import filecmp
from collections import OrderedDict
from pathlib import Path


DEFAULT_PATTERNS = ["*.f90", "*.py", "*.c", "*.cpp", "*.r", "*make*"]


def matching_files(root: Path, patterns: list[str]) -> dict[str, set[Path]]:
    """Return pattern -> relative file paths matching that pattern under root."""
    out: dict[str, set[Path]] = OrderedDict()
    for pattern in patterns:
        matches: set[Path] = set()
        for path in root.rglob(pattern):
            if path.is_file():
                matches.add(path.relative_to(root))
        out[pattern] = matches
    return out


def print_group(title: str, rows: list[Path]) -> None:
    """Print one named group of relative paths."""
    print(title)
    print("-" * len(title))
    if rows:
        for row in rows:
            print(row.as_posix())
    else:
        print("(none)")
    print()


def compare_pattern(dir1: Path, dir2: Path, pattern: str, files1: set[Path], files2: set[Path]) -> None:
    """Print files only in one directory and common files with different contents."""
    only1 = sorted(files1 - files2)
    only2 = sorted(files2 - files1)
    common = sorted(files1 & files2)
    distinct = [rel for rel in common if not filecmp.cmp(dir1 / rel, dir2 / rel, shallow=False)]

    print("=" * 80)
    print(f"File type: {pattern}")
    print("=" * 80)
    print_group(f"Only in {dir1}", only1)
    print_group(f"Only in {dir2}", only2)
    print_group("Common but distinct", distinct)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="List selected files only in directory 1, only in directory 2, and common files with different contents."
    )
    parser.add_argument("dir1", type=Path)
    parser.add_argument("dir2", type=Path)
    parser.add_argument(
        "--patterns",
        nargs="+",
        default=DEFAULT_PATTERNS,
        help="glob patterns to compare, default: *.f90 *.py *.c *.cpp *.r *make*",
    )
    args = parser.parse_args()

    dir1 = args.dir1.resolve()
    dir2 = args.dir2.resolve()
    if not dir1.is_dir():
        raise SystemExit(f"not a directory: {dir1}")
    if not dir2.is_dir():
        raise SystemExit(f"not a directory: {dir2}")

    files1 = matching_files(dir1, args.patterns)
    files2 = matching_files(dir2, args.patterns)
    print(f"Directory 1: {dir1}")
    print(f"Directory 2: {dir2}")
    print(f"Patterns: {' '.join(args.patterns)}")
    print()

    for pattern in args.patterns:
        compare_pattern(dir1, dir2, pattern, files1[pattern], files2[pattern])


if __name__ == "__main__":
    main()
