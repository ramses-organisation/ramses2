#!/usr/bin/env python3
"""
Thermal history of a uniform universe (namelist/thermal_history.nml),
read from the log file of a serial 1D run, where each output is printed
as a table of cells (lev x d u P [xHII xHeII xHeIII]) instead of files.

T/mu = P/d * scale_T2(aexp), with scale_T2 = mH/kB * scale_v**2 and
scale_v = 1e7 * boxlen_ini / aexp [cm/s] (cosmological code units).
With non-equilibrium chemistry, the ionization fractions give
mu = 1 / (X*(1+xHII) + Y/4*(1+xHeII+2*xHeIII)) and the temperature T.
Without them (equilibrium cooling), only T/mu is available, unless --zreion
is given: T is then estimated with mu=1.22 (neutral) before z_reion and
mu=0.59 (ionized) after.

Usage:
    python3 thermal_history.py run1.log [run2.log ...] [-o thermal_history.png]
"""
import argparse
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

mH = 1.6605390e-24  # same constants as cooling/constants.f90
kB = 1.3806490e-16


def read_log(filename, X=0.76, Y=0.24):
    """Return aexp, T/mu, T (None without ionization fractions) and xHII (or None)."""
    with open(filename) as f:
        lines = f.readlines()
    aexp, T2, T, xH = [], [], [], []
    a = None
    boxlen = None
    k = 0
    while k < len(lines):
        line = lines[k]
        if line.startswith(" aexp=") and a is None:
            a = float(re.search(r"aexp=\s*(\S+)", line).group(1))
        elif line.startswith(" box size="):
            boxlen = float(re.search(r"=\s*(\S+)", line).group(1))
        elif line.startswith(" Fine step="):
            a = float(re.search(r" a=\s*(\S+)", line).group(1))
        elif line.startswith(" Output ") and line.strip().endswith("cells"):
            ncell = int(line.split()[1])
            head = lines[k + 1].split()
            rows = np.array([l.split() for l in lines[k + 2:k + 2 + ncell]], dtype=float)
            col = {name: rows[:, i] for i, name in enumerate(head)}
            scale_T2 = mH / kB * (1e7 * boxlen / a) ** 2
            t2 = np.mean(col["P"] / col["d"]) * scale_T2
            aexp.append(a)
            T2.append(t2)
            if "xHII" in col:
                xh = np.mean(col["xHII"])
                xhe1 = np.mean(col.get("xHeII", 0.0))
                xhe2 = np.mean(col.get("xHeIII", 0.0))
                mu = 1.0 / (X * (1 + xh) + Y / 4 * (1 + xhe1 + 2 * xhe2))
                T.append(t2 * mu)
                xH.append(xh)
            k += 2 + ncell
            continue
        k += 1
    has_x = len(T) == len(T2) and len(T) > 0
    return (np.array(aexp), np.array(T2),
            np.array(T) if has_x else None, np.array(xH) if has_x else None)


def main():
    parser = argparse.ArgumentParser(description="Temperature versus aexp from RAMSES log files")
    parser.add_argument("logs", nargs="+", help="log file(s) of thermal_history.nml runs")
    parser.add_argument("-o", "--output", default="thermal_history.png", help="output figure")
    parser.add_argument("--zreion", type=float, default=None,
                        help="for logs without ionization data: T = T/mu * (1.22 before, 0.59 after z_reion)")
    args = parser.parse_args()

    fig, (ax, bx) = plt.subplots(2, 1, figsize=(7.5, 7), sharex=True,
                                 gridspec_kw=dict(height_ratios=[2.2, 1]))
    a = np.logspace(np.log10(5e-4), 0, 200)
    ax.plot(a, 2.726 / a, "k:", lw=1, label="$T_{CMB}$")
    for log in args.logs:
        aexp, T2, T, xH = read_log(log)
        if T is not None:
            ax.plot(aexp, T, "-", label=log + ": T")
            bx.plot(aexp, xH, "-", label=log)
            np.savetxt(log + ".T.txt", np.column_stack((aexp, T2, T, xH)),
                       header="aexp T/mu[K] T[K] xHII")
        elif args.zreion is not None:
            T = T2 * np.where(1 / aexp - 1 < args.zreion, 0.59, 1.22)
            ax.plot(aexp, T, "--", label=log + ": T ($\\mu$=1.22, 0.59 after $z_{reion}$)")
            np.savetxt(log + ".T.txt", np.column_stack((aexp, T2, T)), header="aexp T/mu[K] T[K]")
        else:
            ax.plot(aexp, T2, "--", label=log + ": T/mu (no ionization data)")
            np.savetxt(log + ".T.txt", np.column_stack((aexp, T2)), header="aexp T/mu[K]")
        print(log, ":", len(aexp), "outputs")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel("T [K]  (mean density)")
    sec = ax.secondary_xaxis("top", functions=(lambda x: 1 / x - 1, lambda z: 1 / (1 + z)))
    zt = [1000, 300, 100, 30, 10, 5, 2, 1, 0]
    sec.set_xticks(zt)
    sec.set_xticklabels([str(z) for z in zt])
    sec.set_xlabel("z")
    ax.legend(loc="lower left", fontsize=8)
    bx.set_yscale("log")
    bx.set_ylabel("xHII")
    bx.set_xlabel("aexp")
    fig.tight_layout()
    fig.savefig(args.output, dpi=130)
    print("Saved", args.output)


if __name__ == "__main__":
    main()
