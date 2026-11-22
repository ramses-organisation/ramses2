#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from plot_backend import set_backend
import argparse
import ramses as ram

parser = argparse.ArgumentParser(description="Compare the halo occupation fractions (fraction of halos hosting at least one galaxy) of several runs from the clump catalogs")
parser.add_argument("nout", help="enter output number")
parser.add_argument("paths", nargs="+", help="run directories containing output_*")
parser.add_argument("--labels", nargs="+", help="legend labels, one per run")
parser.add_argument("--nbin", help="number of bins per dex")
parser.add_argument("--npart", help="minimum number of dark matter particles per halo (default 100)")
parser.add_argument("--nstar", help="minimum number of star particles per galaxy (default 1)")
parser.add_argument("--mstar", help="minimum stellar mass per galaxy in Msol (default 0)")
parser.add_argument("--out", help="output a png image")
args = parser.parse_args()
no_display = set_backend()

nout = args.nout
labels = args.labels if args.labels else args.paths
if len(labels) != len(args.paths):
    print("Need one label per run")
    exit()
if args.nbin==None:
    nbin=4
else:
    nbin=int(args.nbin)
npart_min = int(args.npart) if args.npart else 100
nstar_min = int(args.nstar) if args.nstar else 1
mstar_min = float(args.mstar) if args.mstar else 0.0

# Halo masses in Msol and occupation flags for each run
masses = []
occupied = []
for path in args.paths:
    print("Reading output number ",nout," in ",path)
    i=ram.rd_info(nout,path=path+"/")
    c=ram.rd_clump(nout,path=path+"/")
    s=ram.rd_part(nout,path=path+"/",prefix='star',peak=True)
    unit_m = i.unit_d * i.unit_l**3 / 2e33

    # Number of star particles and stellar mass in each peak (peak id = clump index)
    inpeak = s.peak_id > 0
    pid, inv = np.unique(s.peak_id[inpeak], return_inverse=True)
    npeak = np.bincount(inv)
    mpeak = np.bincount(inv, weights=s.mass[inpeak])

    # Match peaks to clumps in the catalog
    order = np.argsort(c.index)
    pos = np.clip(np.searchsorted(c.index, pid, sorter=order), 0, len(c.index)-1)
    found = c.index[order[pos]] == pid
    nstar = np.zeros(len(c.index), dtype=int)
    mstar = np.zeros(len(c.index))
    nstar[order[pos[found]]] = npeak[found]
    mstar[order[pos[found]]] = mpeak[found]

    # A galaxy is a clump (central or satellite) with enough stars
    galaxy = (nstar >= nstar_min) & (mstar * unit_m >= mstar_min) & (nstar > 0)

    # A halo is occupied if any of its clumps is a galaxy
    hpos = order[np.searchsorted(c.index, c.halo, sorter=order)]
    occ = np.zeros(len(c.index), dtype=bool)
    np.logical_or.at(occ, hpos, galaxy)

    # Halos are central clumps (index = parent = halo), mass = mpatch
    halo = (c.index == c.parent) & (c.index == c.halo) & (c.npart >= npart_min) & (c.mpatch > 0)
    print("Halos with npart>=%d: N=%d  with a galaxy: %d" % (npart_min, np.count_nonzero(halo), np.count_nonzero(halo & occ)))
    masses.append(c.mpatch[halo] * unit_m)
    occupied.append(occ[halo])

# Common logarithmic bins
lmin = np.floor(np.log10(min([np.min(m) for m in masses]))*nbin)/nbin
lmax = np.ceil(np.log10(max([np.max(m) for m in masses]))*nbin)/nbin
bins = np.linspace(lmin, lmax, int(round((lmax-lmin)*nbin))+1)
centers = 0.5*(bins[:-1]+bins[1:])

nhalo = [np.histogram(np.log10(m), bins=bins)[0] for m in masses]
nocc = [np.histogram(np.log10(m[o]), bins=bins)[0] for m, o in zip(masses, occupied)]
fracs = []
for nh, no in zip(nhalo, nocc):
    with np.errstate(divide='ignore', invalid='ignore'):
        fracs.append(np.where(nh > 0, no/nh, np.nan))

# Summary and table: occupied / total halos and occupation fraction per bin
print()
for label, m, o in zip(labels, masses, occupied):
    f = np.count_nonzero(o)/len(m) if len(m) > 0 else 0.0
    print("%-12s N=%7d  occupied=%7d  fraction=%6.3f" % (label, len(m), np.count_nonzero(o), f))
print()
print("log10 M    " + "".join(["%18s" % l for l in labels]))
for ib in range(len(centers)):
    line = "%6.2f     " % centers[ib]
    for nh, no, f in zip(nhalo, nocc, fracs):
        if nh[ib] > 0:
            line += "%18s" % ("%d/%d %5.3f" % (no[ib], nh[ib], f[ib]))
        else:
            line += "%18s" % "-"
    print(line)

fig, ax1 = plt.subplots(figsize=(8, 4.8))

for label, f, nh in zip(labels, fracs, nhalo):
    line, = ax1.step(centers, f, where='mid', label=label)
    # Binomial error on the fraction
    with np.errstate(divide='ignore', invalid='ignore'):
        err = np.where(nh > 0, np.sqrt(f*(1-f)/nh), np.nan)
    ax1.errorbar(centers, f, yerr=err, fmt='none', ecolor=line.get_color(), lw=0.8)

ax1.set_ylim(0, 1.05)
ax1.set_ylabel(r'$f_{\rm occ}$ (halos with $\geq$ 1 galaxy)')
ax1.legend(fontsize='small', ncol=1 if len(labels) <= 12 else 2,
           loc='upper left', bbox_to_anchor=(1.02, 1), borderaxespad=0)
ax1.set_xlabel(r'log$_{10}$ M$_{\rm halo}$ [Msol]')
plt.tight_layout()

if args.out:
    plt.savefig(args.out)

if not no_display:
    plt.show()
