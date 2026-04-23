#!/usr/bin/env python3
import os
import argparse
from collections import OrderedDict

CODENAME = "surf_bc_merge >"
FILENAME_D_BC = "D_bc.dat"
DEF_DIRECTORY = "."


def read_bc_file_grouped(filepath: str):
    with open(filepath, "r", encoding="utf-8") as f:
        lines = [line.strip() for line in f if line.strip()]

    if not lines:
        raise ValueError(f'Empty file: "{filepath}"')

    header = lines[0].split()
    if len(header) < 2:
        raise ValueError(f'Invalid header in "{filepath}"')

    nbc = int(header[0])
    b = int(header[1])

    grouped = OrderedDict()
    for i, line in enumerate(lines[1:], start=2):
        parts = line.split()
        if len(parts) < 3:
            raise ValueError(f'Invalid line {i} in "{filepath}": {line}')

        nnum = int(parts[0])
        bnum = int(parts[1])
        val = float(parts[2])

        if nnum not in grouped:
            grouped[nnum] = []
        grouped[nnum].append((nnum, bnum, val))

    actual_nbc = sum(len(v) for v in grouped.values())
    if actual_nbc != nbc:
        print(
            f'{CODENAME} Warning: header count ({nbc}) != actual line count ({actual_nbc}) in "{os.path.basename(filepath)}"'
        )

    return b, grouped


def merge_bc_files(infiles, directory, prefer="last"):
    merged = OrderedDict()
    base_b = None

    for infile in infiles:
        path = os.path.join(directory, infile)
        b, grouped = read_bc_file_grouped(path)

        if base_b is None:
            base_b = b
        elif b != base_b:
            print(
                f'{CODENAME} Warning: header second value differs: "{infile}" has {b}, first file has {base_b}'
            )

        count = sum(len(v) for v in grouped.values())
        print(f'{CODENAME} The number of BCs in "{infile}": {count}')

        for nnum, block in grouped.items():
            if prefer == "first":
                if nnum not in merged:
                    merged[nnum] = block
            else:
                merged[nnum] = block

    flat_data = []
    for block in merged.values():
        flat_data.extend(block)

    return base_b, flat_data


def write_bc_file(outfile, directory, b, data):
    outpath = os.path.join(directory, outfile)
    with open(outpath, "w", encoding="utf-8") as f:
        f.write(f"{len(data)} {b}\n")
        for nnum, bnum, val in data:
            f.write(f"{nnum} {bnum} {val:.15e}\n")


def main():
    parser = argparse.ArgumentParser(
        prog="surf_bc_merge.py",
        description="Merge multiple B.C. files and remove duplicated nodes."
    )
    parser.add_argument(
        "infiles",
        nargs="+",
        help="input B.C. files"
    )
    parser.add_argument(
        "-d", "--directory",
        default=DEF_DIRECTORY,
        help="input & output directory"
    )
    parser.add_argument(
        "-o", "--outfile",
        default=FILENAME_D_BC,
        help="output filename"
    )
    parser.add_argument(
        "--prefer",
        choices=["first", "last"],
        default="last",
        help="which duplicated node entry to keep (default: last)"
    )

    args = parser.parse_args()

    print()
    if args.directory == DEF_DIRECTORY:
        print(f"{CODENAME} Input & output directory: {args.directory} (default)")
    else:
        print(f"{CODENAME} Input & output directory: {args.directory}")

    if args.outfile == FILENAME_D_BC:
        print(f"{CODENAME} Output filename for Dirichlet B.C.: {args.outfile} (default)")
    else:
        print(f"{CODENAME} Output filename for Dirichlet B.C.: {args.outfile}")

    print(f"{CODENAME} Duplicate node policy: keep {args.prefer}")
    print()

    b, merged_data = merge_bc_files(args.infiles, args.directory, prefer=args.prefer)

    print(f'{CODENAME} The number of merged BCs after deduplication: {len(merged_data)}')

    write_bc_file(args.outfile, args.directory, b, merged_data)

    print()


if __name__ == "__main__":
    main()