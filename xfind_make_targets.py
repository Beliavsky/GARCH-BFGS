#!/usr/bin/env python3
"""Find make targets that build and run a Fortran main program."""

from __future__ import annotations

import argparse
import os
import re
import subprocess
from pathlib import Path


VAR_REF_RE = re.compile(r"\$\(([^)]+)\)")
ASSIGN_RE = re.compile(r"^([A-Za-z_][A-Za-z0-9_]*)\s*(?::=|\?=|\+=|=)\s*(.*)$")
TARGET_RE = re.compile(r"^([^#\s][^:=]*?)\s*:\s*(.*)$")
MODULE_RE = re.compile(r"^\s*module\s+([A-Za-z_][A-Za-z0-9_]*)\b", re.IGNORECASE)
USE_RE = re.compile(r"^\s*use(?:\s*,\s*[^:]+)?\s*(?:::)?\s*([A-Za-z_][A-Za-z0-9_]*)\b", re.IGNORECASE)


def expand_vars(text: str, variables: dict[str, str], limit: int = 20) -> str:
    """Expand simple $(NAME) make variables using assignments found in the file."""

    def repl(match: re.Match[str]) -> str:
        return variables.get(match.group(1), match.group(0))

    old = text
    for _ in range(limit):
        new = VAR_REF_RE.sub(repl, old)
        if new == old:
            return new
        old = new
    return old


def strip_inline_comment(line: str) -> str:
    """Remove simple inline comments. This is enough for this makefile style."""

    return line.split("#", 1)[0].rstrip()


def parse_makefile(path: Path) -> tuple[dict[str, str], list[tuple[int, list[str], str, str]]]:
    variables: dict[str, str] = {}
    targets: list[tuple[int, list[str], str, str]] = []

    logical_lines: list[tuple[int, str]] = []
    current = ""
    start_line = 0
    for lineno, raw in enumerate(path.read_text().splitlines(), start=1):
        line = raw.rstrip()
        if current:
            current += line[:-1].strip() if line.endswith("\\") else line.strip()
        else:
            start_line = lineno
            current = line[:-1].strip() if line.endswith("\\") else line
        if line.endswith("\\"):
            current += " "
            continue
        logical_lines.append((start_line, current))
        current = ""
    if current:
        logical_lines.append((start_line, current))

    for lineno, raw in logical_lines:
        line = strip_inline_comment(raw)
        if not line or line.startswith("\t"):
            continue

        assign = ASSIGN_RE.match(line)
        if assign:
            variables[assign.group(1)] = assign.group(2).strip()
            continue

        target_match = TARGET_RE.match(line)
        if target_match:
            lhs = target_match.group(1).strip()
            rhs = target_match.group(2).strip()
            names = [name for name in lhs.split() if name]
            targets.append((lineno, names, rhs, raw))

    return variables, targets


def command_prefix(makefile: Path) -> str:
    if makefile.name.lower() == "makefile" and makefile.parent in (Path("."), Path("")):
        return "make"
    return f"make -f {makefile}"


def append_unique(items: list, item) -> None:
    if item not in items:
        items.append(item)


def object_source_map(
    targets: list[tuple[int, list[str], str, str]], variables: dict[str, str]
) -> dict[str, str]:
    source_by_object: dict[str, str] = {}
    for _, names, deps, _ in targets:
        expanded_names = [expand_vars(name, variables) for name in names]
        dep_words = expand_vars(deps, variables).split()
        source_words = [word for word in dep_words if word.lower().endswith(".f90")]
        if not source_words:
            continue
        for name in expanded_names:
            if name.lower().endswith(".o"):
                source_by_object[name] = source_words[0]
    return source_by_object


def shell_quote(text: str) -> str:
    if not text:
        return '""'
    if re.search(r'[\s"&|<>^]', text):
        return '"' + text.replace('"', r'\"') + '"'
    return text


def source_modules(source: str) -> tuple[set[str], set[str]]:
    path = Path(source)
    if not path.exists():
        return set(), set()

    provides: set[str] = set()
    uses: set[str] = set()
    for raw in path.read_text(errors="ignore").splitlines():
        line = raw.split("!", 1)[0]
        module_match = MODULE_RE.match(line)
        if module_match and not re.match(r"^\s*module\s+procedure\b", line, re.IGNORECASE):
            provides.add(module_match.group(1).lower())
            continue
        use_match = USE_RE.match(line)
        if use_match:
            uses.add(use_match.group(1).lower())
    return provides, uses


def order_sources_by_modules(sources: list[str]) -> list[str]:
    source_set = set(sources)
    provider: dict[str, str] = {}
    uses_by_source: dict[str, set[str]] = {}

    for source in sources:
        provides, uses = source_modules(source)
        uses_by_source[source] = uses
        for module in provides:
            provider[module] = source

    dependencies: dict[str, set[str]] = {}
    for source in sources:
        dependencies[source] = {
            provider[module]
            for module in uses_by_source[source]
            if module in provider and provider[module] != source
        }

    ordered: list[str] = []
    remaining = list(sources)
    while remaining:
        ready = [source for source in remaining if not (dependencies[source] & source_set)]
        if not ready:
            ordered.extend(remaining)
            break
        for source in ready:
            append_unique(ordered, source)
            remaining.remove(source)
            source_set.remove(source)
    return ordered


def single_compile_command(
    exe: str,
    build_targets: list[tuple[int, str, str]],
    targets: list[tuple[int, list[str], str, str]],
    variables: dict[str, str],
    no_options: bool = False,
) -> tuple[str, list[str]] | None:
    if not build_targets:
        return None

    build_line = build_targets[0][2]
    target_match = TARGET_RE.match(strip_inline_comment(build_line))
    if not target_match:
        return None

    dep_words = expand_vars(target_match.group(2).strip(), variables).split()
    source_by_object = object_source_map(targets, variables)
    sources: list[str] = []
    for dep in dep_words:
        if dep.lower().endswith(".f90"):
            append_unique(sources, dep)
        elif dep.lower().endswith(".o"):
            source = source_by_object.get(dep)
            if source is None:
                candidate = f"{Path(dep).stem}.f90"
                if Path(candidate).exists():
                    source = candidate
            if source is not None:
                append_unique(sources, source)

    if not sources:
        return None
    sources = order_sources_by_modules(sources)

    fc = expand_vars(variables.get("FC", "gfortran"), variables)
    if no_options:
        pieces = [fc, *sources]
    else:
        fflags = expand_vars(variables.get("FFLAGS", ""), variables).split()
        pieces = [fc, *fflags, "-o", exe, *sources]
    pieces = [piece for piece in pieces if piece]
    return " ".join(shell_quote(piece) for piece in pieces), pieces


def compile_single_source_command(
    exe: str, single_source: Path, variables: dict[str, str], no_options: bool = False
) -> tuple[str, list[str]]:
    fc = expand_vars(variables.get("FC", "gfortran"), variables)
    if no_options:
        pieces = [fc, str(single_source)]
    else:
        fflags = expand_vars(variables.get("FFLAGS", ""), variables).split()
        pieces = [fc, *fflags, "-o", exe, str(single_source)]
    pieces = [piece for piece in pieces if piece]
    return " ".join(shell_quote(piece) for piece in pieces), pieces


def source_dependencies(
    build_targets: list[tuple[int, str, str]],
    targets: list[tuple[int, list[str], str, str]],
    variables: dict[str, str],
) -> list[str] | None:
    if not build_targets:
        return None

    build_line = build_targets[0][2]
    target_match = TARGET_RE.match(strip_inline_comment(build_line))
    if not target_match:
        return None

    dep_words = expand_vars(target_match.group(2).strip(), variables).split()
    source_by_object = object_source_map(targets, variables)
    sources: list[str] = []
    for dep in dep_words:
        if dep.lower().endswith(".f90"):
            append_unique(sources, dep)
        elif dep.lower().endswith(".o"):
            source = source_by_object.get(dep)
            if source is None:
                candidate = f"{Path(dep).stem}.f90"
                if Path(candidate).exists():
                    source = candidate
            if source is not None:
                append_unique(sources, source)

    if not sources:
        return None
    return order_sources_by_modules(sources)


def write_single_source(output_file: Path, sources: list[str]) -> None:
    with output_file.open("w", newline="\n") as out:
        out.write("! This file was generated by find_make_targets.py.\n")
        out.write("! Source files are concatenated in module dependency order.\n\n")
        for source in sources:
            out.write(f"! ===== begin {source} =====\n")
            out.write(Path(source).read_text())
            out.write(f"\n! ===== end {source} =====\n\n")


def executable_command(exe: str) -> list[str]:
    path = Path(exe)
    if path.is_absolute() or path.parent != Path("."):
        return [str(path)]
    prefix = ".\\" if os.name == "nt" else "./"
    return [f"{prefix}{exe}"]


def default_compiler_executable() -> str:
    return "a.exe" if os.name == "nt" else "a.out"


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Find make targets that build and run a Fortran main program."
    )
    parser.add_argument("main_program", help="Fortran main program, e.g. xfit_gen_garch_returns.f90")
    parser.add_argument("--f", default="Makefile", help="makefile to scan; default: Makefile")
    parser.add_argument("--verbose", action="store_true", help="print build target, object rule, and related entries")
    parser.add_argument("--single-command", action="store_true", help="print one compiler command using source dependencies")
    parser.add_argument("--no-options", action="store_true", help="omit compiler flags and output option in single-file commands")
    parser.add_argument(
        "--single-source",
        nargs="?",
        const="",
        metavar="OUT.f90",
        help="write one source file with dependencies prepended; default: <main>_single.f90",
    )
    action = parser.add_mutually_exclusive_group()
    action.add_argument("--build", action="store_true", help="build the executable target")
    action.add_argument("--run", action="store_true", help="build and run using the run target")
    args = parser.parse_args()

    makefile = Path(args.f)
    if not makefile.exists():
        raise SystemExit(f"makefile not found: {makefile}")

    main_path = Path(args.main_program)
    stem = main_path.name
    if stem.lower().endswith(".f90"):
        stem = stem[:-4]
    elif stem.lower().endswith(".exe"):
        stem = stem[:-4]
    elif stem.lower().endswith(".o"):
        stem = stem[:-2]
    source = f"{stem}.f90"
    default_obj = f"{stem}.o"
    default_exe = f"{stem}.exe"

    variables, targets = parse_makefile(makefile)
    build_targets: list[tuple[int, str, str]] = []
    run_targets: list[tuple[int, str, str]] = []
    object_targets: list[tuple[int, str, str]] = []
    related_lines: list[tuple[int, str]] = []
    object_names: set[str] = {default_obj}
    executable_names: set[str] = {default_exe}
    object_variables: set[str] = set()

    for name, value in variables.items():
        expanded = expand_vars(value, variables)
        words = expanded.split()
        if source in words or default_obj in words or default_exe in words or expanded == default_exe:
            append_unique(related_lines, (0, f"{name} = {value}"))

    for lineno, names, deps, raw in targets:
        expanded_names = [expand_vars(name, variables) for name in names]
        expanded_deps = expand_vars(deps, variables)
        dep_words = expanded_deps.split()

        if source in dep_words:
            for expanded_name in expanded_names:
                if expanded_name.endswith(".o"):
                    object_names.add(expanded_name)
                    append_unique(object_targets, (lineno, expanded_name, raw))
        if default_obj in expanded_names:
            append_unique(object_targets, (lineno, default_obj, raw))

        if any(token in raw for token in (source, default_obj, default_exe)):
            append_unique(related_lines, (lineno, raw))

    for name, value in variables.items():
        expanded = expand_vars(value, variables)
        words = set(expanded.split())
        if words & object_names:
            object_variables.add(name)
            append_unique(related_lines, (0, f"{name} = {value}"))

    for lineno, names, deps, raw in targets:
        expanded_names = [expand_vars(name, variables) for name in names]
        expanded_deps = expand_vars(deps, variables)
        dep_words = set(expanded_deps.split())
        raw_deps = set(deps.split())

        has_object_dep = bool(dep_words & object_names) or any(f"$({name})" in raw_deps for name in object_variables)
        if has_object_dep:
            for expanded_name in expanded_names:
                if expanded_name.endswith(".exe") and (lineno, expanded_name, raw) not in build_targets:
                    append_unique(build_targets, (lineno, expanded_name, raw))
                    executable_names.add(expanded_name)
                    append_unique(related_lines, (lineno, raw))

    for lineno, names, deps, raw in targets:
        expanded_deps = expand_vars(deps, variables)
        dep_words = set(expanded_deps.split())
        if dep_words & executable_names:
            for name in names:
                expanded_name = expand_vars(name, variables)
                if expanded_name == "run" or expanded_name.startswith("run_"):
                    append_unique(run_targets, (lineno, expanded_name, raw))

    prefix = command_prefix(makefile)
    make_args = ["make"]
    if not (makefile.name.lower() == "makefile" and makefile.parent in (Path("."), Path(""))):
        make_args.extend(["-f", str(makefile)])

    if args.single_source is not None:
        if not build_targets:
            raise SystemExit("build target not found")
        sources = source_dependencies(build_targets, targets, variables)
        if sources is None:
            raise SystemExit("source dependencies could not be inferred")
        output_file = Path(args.single_source) if args.single_source else Path(f"{stem}_single.f90")
        write_single_source(output_file, sources)
        print("Single source file:")
        print(f"  {output_file}")
        if args.verbose:
            print("\nIncluded sources:")
            for source_name in sources:
                print(f"  {source_name}")
        if args.build or args.run:
            command, command_args = compile_single_source_command(
                build_targets[0][1], output_file, variables, no_options=args.no_options
            )
            print("Single-source build command:", flush=True)
            print(f"  {command}", flush=True)
            result = subprocess.run(command_args)
            if result.returncode != 0 or args.build:
                return result.returncode
            run_exe = default_compiler_executable() if args.no_options else build_targets[0][1]
            run_cmd = executable_command(run_exe)
            print("Run command:", flush=True)
            print(f"  {' '.join(shell_quote(piece) for piece in run_cmd)}", flush=True)
            return subprocess.run(run_cmd).returncode
        return 0

    if args.single_command:
        if not build_targets:
            raise SystemExit("build target not found")
        command_info = single_compile_command(
            build_targets[0][1], build_targets, targets, variables, no_options=args.no_options
        )
        if command_info is None:
            raise SystemExit("single compiler command could not be inferred")
        command, command_args = command_info
        print("Single compiler command:", flush=True)
        print(f"  {command}", flush=True)
        if args.build or args.run:
            result = subprocess.run(command_args)
            if result.returncode != 0 or args.build:
                return result.returncode
            run_exe = default_compiler_executable() if args.no_options else build_targets[0][1]
            run_cmd = executable_command(run_exe)
            print("Run command:", flush=True)
            print(f"  {' '.join(shell_quote(piece) for piece in run_cmd)}", flush=True)
            return subprocess.run(run_cmd).returncode
        return 0

    if args.build:
        if not build_targets:
            raise SystemExit("build target not found")
        cmd = make_args + [build_targets[0][1]]
        print("Build command:", flush=True)
        print(f"  {' '.join(cmd)}", flush=True)
        return subprocess.run(cmd).returncode

    if args.run:
        if not run_targets:
            raise SystemExit("build-and-run target not found")
        cmd = make_args + [run_targets[0][1]]
        print("Build-and-run command:", flush=True)
        print(f"  {' '.join(cmd)}", flush=True)
        return subprocess.run(cmd).returncode

    if run_targets:
        if args.verbose:
            print(f"Main program: {source}")
            print(f"Makefile:     {makefile}")

            if build_targets:
                print("\nBuild target:")
                for _, target, _ in build_targets:
                    print(f"  {target}")
                print("\nBuild command:")
                print(f"  {prefix} {build_targets[0][1]}")
            else:
                print("\nBuild target: not found")

            print("\nRun target:")
            for _, target, _ in run_targets:
                print(f"  {target}")
            print()

        print("Build-and-run command:")
        print(f"  {prefix} {run_targets[0][1]}")
    else:
        if args.verbose:
            print(f"Main program: {source}")
            print(f"Makefile:     {makefile}")
            if build_targets:
                print("\nBuild target:")
                for _, target, _ in build_targets:
                    print(f"  {target}")
                print("\nBuild command:")
                print(f"  {prefix} {build_targets[0][1]}")
            else:
                print("\nBuild target: not found")
            print("\nRun target: not found")
        else:
            print("Build-and-run command: not found")

    if args.verbose and object_targets:
        print("\nObject/source rule:")
        for lineno, _, raw in object_targets:
            print(f"  line {lineno}: {raw}")

    if args.verbose and related_lines:
        print("\nRelated makefile entries:")
        for lineno, raw in related_lines:
            label = "variable" if lineno == 0 else f"line {lineno}"
            print(f"  {label}: {raw}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
