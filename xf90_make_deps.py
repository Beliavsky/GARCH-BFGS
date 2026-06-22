#!/usr/bin/env python3
"""List local Fortran sources by the number of Makefile targets depending on them."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


VAR_REF_RE = re.compile(r"\$\(([^)]+)\)|\${([^}]+)}")


def logical_make_lines(path: Path) -> list[str]:
    """Return Makefile logical lines with trailing backslash continuations joined."""
    lines: list[str] = []
    current = ""
    for raw in path.read_text().splitlines():
        line = raw.rstrip()
        if line.endswith("\\"):
            current += line[:-1] + " "
        else:
            lines.append(current + line)
            current = ""
    if current:
        lines.append(current)
    return lines


def strip_comment(line: str) -> str:
    """Remove comments; this Makefile does not use escaped hash characters."""
    return line.split("#", 1)[0].strip()


def parse_variables(lines: list[str]) -> dict[str, str]:
    """Parse simple Makefile variable assignments."""
    variables: dict[str, str] = {}
    for line in lines:
        text = strip_comment(line)
        if not text or text.startswith("\t"):
            continue
        match = re.match(r"^([A-Za-z_][A-Za-z0-9_]*)\s*(?::=|\+=|=)\s*(.*)$", text)
        if not match:
            continue
        name, value = match.groups()
        if "+=" in text.split(name, 1)[1].split(value, 1)[0]:
            variables[name] = (variables.get(name, "") + " " + value).strip()
        else:
            variables[name] = value.strip()
    return variables


def expand_vars(text: str, variables: dict[str, str], depth: int = 0) -> str:
    """Expand simple $(VAR) and ${VAR} references."""
    if depth > 20:
        return text

    def repl(match: re.Match[str]) -> str:
        name = match.group(1) or match.group(2)
        return expand_vars(variables.get(name, ""), variables, depth + 1)

    expanded = VAR_REF_RE.sub(repl, text)
    if expanded == text:
        return expanded
    return expand_vars(expanded, variables, depth + 1)


def rule_target_dependencies(lines: list[str], variables: dict[str, str]) -> list[tuple[list[str], list[str]]]:
    """Return Makefile rule targets and expanded prerequisite tokens."""
    rules: list[tuple[list[str], list[str]]] = []
    for line in lines:
        text = strip_comment(line)
        if not text or text.startswith("\t") or ":" not in text:
            continue
        if re.match(r"^[A-Za-z_][A-Za-z0-9_]*\s*(?::=|\+=|=)", text):
            continue
        left, right = text.split(":", 1)
        if not left.strip():
            continue
        if "|" in right:
            right = right.split("|", 1)[0]
        targets = expand_vars(left, variables).split()
        deps = expand_vars(right, variables).split()
        if targets:
            rules.append((targets, deps))
    return rules


def source_from_dependency(dep: str, f90_files: set[str]) -> str | None:
    """Map a prerequisite token to a local .f90 source name when possible."""
    name = Path(dep).name
    lower = name.lower()
    if lower.endswith(".f90") and name in f90_files:
        return name
    if lower.endswith(".o"):
        candidate = name[:-2] + ".f90"
        if candidate in f90_files:
            return candidate
    return None


def dependency_counts(makefile: Path, source_dir: Path) -> dict[str, set[str]]:
    """Return local .f90 source -> set of Makefile targets depending on it."""
    f90_files = {path.name for path in source_dir.glob("*.f90")}
    lines = logical_make_lines(makefile)
    variables = parse_variables(lines)
    rules = rule_target_dependencies(lines, variables)
    targets_by_source = {name: set() for name in f90_files}
    for targets, deps in rules:
        sources = {source_from_dependency(dep, f90_files) for dep in deps}
        sources.discard(None)
        for source in sources:
            targets_by_source[source].update(targets)
    return targets_by_source


def executable_targets(makefile: Path) -> tuple[set[str], set[str]]:
    """Return all .exe rule targets and the .exe targets listed by all."""
    lines = logical_make_lines(makefile)
    variables = parse_variables(lines)
    rules = rule_target_dependencies(lines, variables)
    exe_targets = {target for targets, _ in rules for target in targets if target.lower().endswith(".exe")}
    all_exes: set[str] = set()
    for targets, deps in rules:
        if "all" in targets:
            all_exes.update(dep for dep in deps if dep.lower().endswith(".exe"))
    return exe_targets, all_exes


def print_x_program_audit(makefile: Path, source_dir: Path) -> None:
    """Print coverage of x*.f90 source files by .exe targets and make all."""
    exe_targets, all_exes = executable_targets(makefile)
    x_sources = sorted(path.name for path in source_dir.glob("x*.f90"))
    expected_exes = {source[:-4] + ".exe" for source in x_sources}
    missing_rules = sorted(source for source in x_sources if source[:-4] + ".exe" not in exe_targets)
    omitted_from_all = sorted(exe for exe in expected_exes & exe_targets if exe not in all_exes)

    print(f"x*.f90 files:                         {len(x_sources)}")
    print(f"x*.f90 files with matching exe rule:  {len(x_sources) - len(missing_rules)}")
    print(f"x*.f90 files without exe rule:        {len(missing_rules)}")
    print(f"matching x*.exe rules omitted by all: {len(omitted_from_all)}")

    if missing_rules:
        print("\nx*.f90 files without matching .exe rule")
        print("--------------------------------------")
        for source in missing_rules:
            print(source)

    if omitted_from_all:
        print("\nmatching x*.exe rules omitted by all")
        print("------------------------------------")
        for exe in omitted_from_all:
            print(exe)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Count how many Makefile targets depend on each .f90 file in the current directory."
    )
    parser.add_argument("--makefile", default="makefile", help="Makefile path, default: makefile")
    parser.add_argument("--show-targets", action="store_true", help="also print the dependent target names")
    parser.add_argument("--audit-x", action="store_true", help="audit x*.f90 files against .exe targets and make all")
    args = parser.parse_args()

    makefile = Path(args.makefile)
    if args.audit_x:
        print_x_program_audit(makefile, Path("."))
        return

    counts = dependency_counts(makefile, Path("."))
    rows = sorted(counts.items(), key=lambda item: (-len(item[1]), item[0].lower()))

    if args.show_targets:
        print(f"{'file':36s} {'targets':>7s}  target_names")
        print("-" * 100)
        for source, targets in rows:
            print(f"{source:36s} {len(targets):7d}  {' '.join(sorted(targets))}")
    else:
        print(f"{'file':36s} {'targets':>7s}")
        print("-" * 46)
        for source, targets in rows:
            print(f"{source:36s} {len(targets):7d}")


if __name__ == "__main__":
    main()
