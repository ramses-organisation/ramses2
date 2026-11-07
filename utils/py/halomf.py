#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from plot_backend import set_backend
import argparse
import ramses as ram

parser = argparse.ArgumentParser(description="Compare the halo mass functions of several runs from the clump catalogs")
parser.add_argument("nout", help="enter output number")
parser.add_argument("paths", nargs="+", help="run directories containing output_*")
parser.add_argument("--labels", nargs="+", help="legend labels, one per run")
parser.add_argument("--nbin", help="number of bins per dex")
parser.add_argument("--npart", help="minimum number of dark matter particles per halo (default 100)")
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

# Halo masses in Msol and comoving box volume in Mpc^3 for each run
masses = []
vols = []
for path in args.paths:
    print("Reading output number ",nout," in ",path)
    i=ram.rd_info(nout,path=path+"/")
    c=ram.rd_clump(nout,path=path+"/")
    unit_m = i.unit_d * i.unit_l**3 / 2e33
    # Halos are central clumps (index = parent = halo), mass = mpatch
    halo = (c.index == c.parent) & (c.index == c.halo) & (c.npart >= npart_min) & (c.mpatch > 0)
    print("Halos with npart>=%d: N=%d" % (npart_min, np.count_nonzero(halo)))
    masses.append(c.mpatch[halo] * unit_m)
    vols.append((i.boxlen * i.unit_l / i.aexp / 3.08568e24)**3)

# Common logarithmic bins
lmin = np.floor(np.log10(min([np.min(m) for m in masses]))*nbin)/nbin
lmax = np.ceil(np.log10(max([np.max(m) for m in masses]))*nbin)/nbin
bins = np.linspace(lmin, lmax, int(round((lmax-lmin)*nbin))+1)
bin_size = bins[1]-bins[0]
centers = 0.5*(bins[:-1]+bins[1:])

hists = [np.histogram(np.log10(m), bins=bins)[0] for m in masses]

# Summary and table, first run as reference
print()
for label, m in zip(labels, masses):
    print("%-12s N=%7d  Mtot=%10.4e Msol  M50=%10.4e  M90=%10.4e  Mmax=%10.4e" %
          (label, len(m), np.sum(m), np.percentile(m,50), np.percentile(m,90), np.max(m)))
print()
print("log10 M    " + "".join(["%12s" % l for l in labels]) +
      "".join(["%12s" % ("ratio "+str(k+1)) for k in range(1,len(labels))]))
for ib in range(len(centers)):
    line = "%6.2f     " % centers[ib] + "".join(["%12d" % h[ib] for h in hists])
    for h in hists[1:]:
        if hists[0][ib] > 0:
            line += "%12.3f" % (h[ib]/hists[0][ib])
        else:
            line += "%12s" % "-"
    print(line)

fig, ax1 = plt.subplots(figsize=(8, 4.8))

for label, h, vol in zip(labels, hists, vols):
    # Differential mass function per unit comoving volume
    with np.errstate(divide='ignore'):
        logphi = np.where(h > 0, np.log10(h/bin_size/vol), np.nan)
    ax1.step(centers, logphi, where='mid', label=label)

ax1.yaxis.set_major_locator(MultipleLocator(0.6))
ax1.set_ylabel(r'log$_{10}$ $\phi$ [Mpc$^{-3}$ dex$^{-1}$]')
ax1.legend(fontsize='small', ncol=1 if len(labels) <= 12 else 2,
           loc='upper left', bbox_to_anchor=(1.02, 1), borderaxespad=0)
ax1.set_xlabel(r'log$_{10}$ M$_{\rm halo}$ [Msol]')
plt.tight_layout()

if args.out:
    plt.savefig(args.out)

if not no_display:
    plt.show()
