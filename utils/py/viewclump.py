#!/usr/bin/env python3
"""Sort a RAMSES text table (clump files, ...) on one column,
keeping the header line on top.

Usage: viewclump.py [-a] <file|output_dir|output_number> [column] [prefix]
  file|output_dir  one table file, or an output directory whose
                   <prefix>.* files (one per rank) are merged
  output_number    shortcut for output_NNNNN in the current directory
  column           column number (from 1) or header name
                   (default: mpatch, or column 5 if there is none)
  prefix           file prefix in an output directory (default: clump)
  -a               ascending order (default: descending)
  --halo ID        instead of sorting, list the properties of the clump
                   whose index is ID, followed by the amr2map options
                   -xce -yce -zce -rad (peak position and r200, in units
                   of boxlen, read from info.txt when there is one)

Examples:
  viewclump.py 8
  viewclump.py 8 mass
  viewclump.py output_00008/clump.00001
  viewclump.py output_00008 mass
  viewclump.py -a output_00008 rho_max
  viewclump.py 8 --halo 642425

The result is shown in less -S on a terminal, printed otherwise
(so it can be piped, e.g. into head).
"""
import glob
import os
import subprocess
import sys


def read_boxlen(table_file):
    """boxlen from info.txt next to the table, 1 if there is none."""
    info = os.path.join(os.path.dirname(table_file) or ".", "info.txt")
    try:
        with open(info) as fh:
            for line in fh:
                if line.split("=")[0].strip() == "boxlen":
                    return float(line.split("=")[1])
    except (OSError, ValueError, IndexError):
        pass
    return 1.0


def show_halo(header, rows, halo, table_file):
    """Print the properties of one clump and the amr2map options."""
    names = header.split()
    target = int(halo)
    row = None
    for line in rows:
        first = line.split(None, 1)
        if first and first[0].isdigit() and int(first[0]) == target:
            row = line.split()
            break
    if row is None:
        sys.stderr.write(f"No clump with index {halo}\n")
        return 1

    width = max(len(n) for n in names)
    for n, v in zip(names, row):
        print(f"{n:>{width}s} = {v}")

    need = ["pos_x", "pos_y", "pos_z", "r200"]
    if not all(n in names for n in need):
        sys.stderr.write("No pos_x/pos_y/pos_z/r200 columns: no amr2map options\n")
        return 0
    val = dict(zip(names, row))
    boxlen = read_boxlen(table_file)
    rad = float(val["r200"])
    if rad <= 0 and "rmax" in val and float(val["rmax"]) > 0:
        sys.stderr.write("r200 is 0 (no mass profile): using rmax\n")
        rad = float(val["rmax"])
    opts = (f"-xce {float(val['pos_x']) / boxlen:.9e} "
            f"-yce {float(val['pos_y']) / boxlen:.9e} "
            f"-zce {float(val['pos_z']) / boxlen:.9e}")
    if rad > 0:
        opts += f" -rad {rad / boxlen:.9e}"
    else:
        sys.stderr.write("r200 and rmax are 0: no -rad\n")
    print()
    print(opts)
    return 0


def main():
    args = sys.argv[1:]
    descending = True
    if "-a" in args:
        descending = False
        args.remove("-a")
    halo = None
    if "--halo" in args:
        i = args.index("--halo")
        if i + 1 >= len(args) or not args[i + 1].isdigit():
            sys.stderr.write("--halo needs a clump index\n")
            return 1
        halo = args[i + 1]
        del args[i:i + 2]
    if not args:
        sys.stderr.write(__doc__)
        return 1
    src = args[0]
    col = args[1] if len(args) > 1 else None
    prefix = args[2] if len(args) > 2 else "clump"

    # A bare output number means output_NNNNN in the current directory
    if src.isdigit() and not os.path.exists(src):
        src = f"output_{int(src):05d}"

    if os.path.isdir(src):
        files = sorted(glob.glob(os.path.join(src, prefix + ".*")))
        if not files:
            sys.stderr.write(f"No {prefix}.* files in {src}\n")
            return 1
    elif os.path.isfile(src):
        files = [src]
    else:
        sys.stderr.write(f"No such file: {src}\n")
        return 1

    header = None
    rows = []
    for f in files:
        with open(f) as fh:
            lines = fh.readlines()
        if not lines:
            continue
        if header is None:
            header = lines[0]
        rows.extend(lines[1:])
    if header is None:
        sys.stderr.write("All files are empty\n")
        return 1

    if halo is not None:
        return show_halo(header, rows, halo, files[0])

    # Turn a column name into its number
    names = header.split()
    if col is None:
        col = "mpatch" if "mpatch" in names else "5"
    if col.isdigit():
        icol = int(col)
        if not 1 <= icol <= len(names):
            sys.stderr.write(f"Column {icol} out of range 1-{len(names)}\n")
            return 1
    elif col in names:
        icol = names.index(col) + 1
    else:
        sys.stderr.write(f"No column named {col}. Columns are:\n")
        for i, n in enumerate(names, 1):
            sys.stderr.write(f"{i:3d} {n}\n")
        return 1

    def key(line):
        fields = line.split(None, icol)
        if len(fields) < icol:
            return float("-inf")
        s = fields[icol - 1]
        try:
            return float(s)
        except ValueError:
            try:
                return float(s.replace("D", "E").replace("d", "e"))
            except ValueError:
                return float("-inf")

    rows.sort(key=key, reverse=descending)

    if sys.stdout.isatty():
        pager = subprocess.Popen(["less", "-S"], stdin=subprocess.PIPE, text=True)
        out = pager.stdin
    else:
        pager = None
        out = sys.stdout
    try:
        out.write(header)
        out.writelines(rows)
        out.flush()
    except BrokenPipeError:
        pass
    finally:
        if pager is not None:
            try:
                out.close()
            except BrokenPipeError:
                pass
            pager.wait()
        else:
            # Quiet exit when the reader (e.g. head) closes the pipe early
            try:
                sys.stdout.flush()
            except BrokenPipeError:
                devnull = os.open(os.devnull, os.O_WRONLY)
                os.dup2(devnull, sys.stdout.fileno())
    return 0


if __name__ == "__main__":
    sys.exit(main())
