#!/usr/bin/env python3
from pathlib import Path
import re
import subprocess
import sys

NAMES = [
    "BBFE_elemmat_fluid_sups_coef_metric_tensor",
    "BBFE_elemmat_fluid_sups_coef_metric_tensor_LSIC",
    "BBFE_elemmat_fluid_sups_mat_NR",
    "BBFE_elemmat_fluid_sups_vec_NR",
]

def fail(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    raise SystemExit(2)

def strip_leading_storage(signature):
    # Keep the actual return type and parameter spelling; only remove a
    # possible leading "static", because these provider symbols must be global.
    return re.sub(r'^\s*static\s+', '', signature, count=1)

def extract_definition_signature(text, name):
    # Locate a function definition by name, then preserve the exact parameter
    # spelling from core_FOM.c.  This matters for C++ mangling (const, pointer
    # rank, array parameter decay, etc.).
    pattern = re.compile(
        r'(?m)^[ \t]*(?:static[ \t]+)?(?:double|void)[ \t]+'
        + re.escape(name)
        + r'[ \t]*\('
    )
    m = pattern.search(text)
    if not m:
        fail(f"definition start not found for {name}")

    open_paren = text.find('(', m.start())
    if open_paren < 0:
        fail(f"opening parenthesis not found for {name}")

    depth = 0
    i = open_paren
    in_string = False
    in_char = False
    in_line_comment = False
    in_block_comment = False
    escape = False

    while i < len(text):
        c = text[i]
        n = text[i + 1] if i + 1 < len(text) else ''

        if in_line_comment:
            if c == '\n':
                in_line_comment = False
        elif in_block_comment:
            if c == '*' and n == '/':
                in_block_comment = False
                i += 1
        elif in_string:
            if escape:
                escape = False
            elif c == '\\':
                escape = True
            elif c == '"':
                in_string = False
        elif in_char:
            if escape:
                escape = False
            elif c == '\\':
                escape = True
            elif c == "'":
                in_char = False
        else:
            if c == '/' and n == '/':
                in_line_comment = True
                i += 1
            elif c == '/' and n == '*':
                in_block_comment = True
                i += 1
            elif c == '"':
                in_string = True
            elif c == "'":
                in_char = True
            elif c == '(':
                depth += 1
            elif c == ')':
                depth -= 1
                if depth == 0:
                    close_paren = i
                    break
        i += 1
    else:
        fail(f"unbalanced parameter list for {name}")

    # Verify this is a definition, not merely a declaration.
    j = close_paren + 1
    while j < len(text) and text[j].isspace():
        j += 1
    if j >= len(text) or text[j] != '{':
        # Search the next occurrence in case a prototype preceded the definition.
        rest = text[m.end():]
        next_match = pattern.search(rest)
        if not next_match:
            fail(f"found declaration but not definition for {name}")
        shifted = m.end() + next_match.start()
        return extract_definition_signature(text[shifted:], name)

    signature = text[m.start():close_paren + 1].strip()
    signature = strip_leading_storage(signature)
    return signature + ";"

def nm_output(obj):
    try:
        return subprocess.check_output(
            ["nm", "-g", str(obj)],
            text=True,
            stderr=subprocess.STDOUT,
        )
    except subprocess.CalledProcessError as exc:
        print(exc.output, file=sys.stderr)
        fail(f"nm failed for {obj}")

def classify_linkage(nm_text, name):
    candidates = [line.strip() for line in nm_text.splitlines() if name in line]
    if not candidates:
        fail(f"provider object does not export a symbol containing {name}")

    # Mach-O C symbol is normally "_NAME". ELF may expose "NAME".
    for line in candidates:
        token = line.split()[-1]
        if token == name or token == "_" + name:
            return "c", line

    # C++ Itanium ABI names contain _Z / __Z and the identifier text.
    for line in candidates:
        token = line.split()[-1]
        if ("_Z" in token or "__Z" in token) and name in token:
            return "cpp", line

    fail(
        f"cannot classify linkage for {name}; candidates:\n  "
        + "\n  ".join(candidates)
    )

def main():
    if len(sys.argv) != 4:
        fail(
            "usage: generate_nr_api_from_core_fom.py "
            "<core_FOM.c> <provider.o> <output.h>"
        )

    src = Path(sys.argv[1])
    obj = Path(sys.argv[2])
    out = Path(sys.argv[3])

    if not src.is_file():
        fail(f"source not found: {src}")
    if not obj.is_file():
        fail(f"provider object not found: {obj}")

    text = src.read_text()
    signatures = {name: extract_definition_signature(text, name) for name in NAMES}

    raw_nm = nm_output(obj)
    modes = {}
    raw_lines = {}
    for name in NAMES:
        mode, line = classify_linkage(raw_nm, name)
        modes[name] = mode
        raw_lines[name] = line

    unique_modes = sorted(set(modes.values()))
    if len(unique_modes) != 1:
        fail(
            "mixed C/C++ linkage across the four provider functions:\n"
            + "\n".join(f"  {n}: {modes[n]}" for n in NAMES)
        )
    linkage = unique_modes[0]

    lines = [
        "#pragma once",
        "",
        "/*",
        " * AUTO-GENERATED from fluid_sups_core/core_FOM.c.",
        " * Do not hand-edit this file.",
        " *",
        " * The exact parameter spelling is copied from the verified provider",
        " * definitions so C++ name mangling cannot drift from the provider.",
        " */",
        "",
    ]

    if linkage == "c":
        lines += [
            "#ifdef __cplusplus",
            'extern "C" {',
            "#endif",
            "",
        ]

    for name in NAMES:
        lines.append(signatures[name])
        lines.append("")

    if linkage == "c":
        lines += [
            "#ifdef __cplusplus",
            "}",
            "#endif",
            "",
        ]

    out.write_text("\n".join(lines))

    print("Generated exact NR API header")
    print(f"  source   : {src}")
    print(f"  provider : {obj}")
    print(f"  output   : {out}")
    print(f"  linkage  : {linkage}")
    print("")
    for name in NAMES:
        print(f"---- {name}")
        print(signatures[name])
        print(f"nm: {raw_lines[name]}")
        print("")

if __name__ == "__main__":
    main()
