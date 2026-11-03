#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from plot_backend import set_backend
import argparse
import ramses as ram

parser = argparse.ArgumentParser(description="Plot the stellar mass fraction of each clump versus clump mass")
parser.add_argument("nout", help="enter output number")
parser.add_argument("--path", help="run directory containing output_* (default ./)")
parser.add_argument("--npart", help="minimum number of dark matter particles per clump (default 100)")
parser.add_argument("--out", help="output a png image")
args = parser.parse_args()
no_display = set_backend()

nout = args.nout
path = args.path+"/" if args.path else "./"
npart_min = int(args.npart) if args.npart else 100

print("Reading output number ",nout," in ",path)
i=ram.rd_info(nout,path=path)
unit_m = i.unit_d * i.unit_l**3 / 2e33

# Star particles with their peak id from the clump finder
s=ram.rd_part(nout,path=path,prefix='star',peak=True)

# Clump catalog
c=ram.rd_clump(nout,path=path)

# Total stellar mass in each peak (peak id = clump index)
inpeak = s.peak_id > 0
pid, inv = np.unique(s.peak_id[inpeak], return_inverse=True)
mpeak = np.bincount(inv, weights=s.mass[inpeak])

# Match peaks to clumps in the catalog
order = np.argsort(c.index)
pos = np.clip(np.searchsorted(c.index, pid, sorter=order), 0, len(c.index)-1)
found = c.index[order[pos]] == pid
mstar = np.zeros(len(c.index))
mstar[order[pos[found]]] = mpeak[found]

print()
print("Stars: N=%d  Mtot=%10.4e Msol" % (s.npart, np.sum(s.mass)*unit_m))
print("Stars in a peak: N=%d  M=%10.4e Msol" % (np.count_nonzero(inpeak), np.sum(s.mass[inpeak])*unit_m))
print("Peaks with stars not in the catalog: %d  M=%10.4e Msol" % (np.count_nonzero(~found), np.sum(mpeak[~found])*unit_m))
print("Clumps: N=%d  with stars: %d" % (len(c.index), np.count_nonzero(mstar > 0)))

# Maximum stellar mass over all clumps of the same halo, stored on the halo clump
hpos = order[np.searchsorted(c.index, c.halo, sorter=order)]
mstar_max = np.zeros(len(c.index))
np.maximum.at(mstar_max, hpos, mstar)
print("Halos where the max stellar mass is not in the central clump: %d" % np.count_nonzero(mstar_max > mstar))
mstar = mstar_max

# Keep only central clumps (index = parent = halo) with at least npart_min particles
central = (c.index == c.parent) & (c.index == c.halo) & (c.npart >= npart_min)
print("Central clumps with npart>=%d: N=%d  with stars: %d" % (npart_min, np.count_nonzero(central), np.count_nonzero(central & (mstar > 0))))

# Stellar mass fraction versus clump mass
mclump = c.mpatch * unit_m
ok = central & (mstar > 0) & (c.mpatch > 0)
ratio = mstar[ok] / c.mpatch[ok]

fig, ax = plt.subplots(figsize=(6,5))
ax.scatter(mclump[ok], ratio, s=2, alpha=0.5)

# Median in log-spaced mass bins (bins with at least 5 halos)
lm = np.log10(mclump[ok])
bins = np.linspace(lm.min(), lm.max(), 16)
ib = np.digitize(lm, bins[1:-1])
kept = [k for k in range(len(bins)-1) if np.count_nonzero(ib==k) >= 5]
if kept:
    xm = [np.median(lm[ib==k]) for k in kept]
    ym = [np.median(ratio[ib==k]) for k in kept]
    ax.plot(10**np.array(xm), ym, color='k', lw=2, label='median')
ax.axhline(i.omega_b/i.omega_m, color='k', ls='--', lw=1, label=r'$\Omega_b/\Omega_m$')
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r'M$_{\rm clump}$ [Msol]')
ax.set_ylabel(r'M$_{\rm star}$/M$_{\rm clump}$')
ax.set_title("output %s  z=%.2f" % (str(nout).zfill(5), i.redshift))
ax.legend()
plt.tight_layout()

if args.out:
    plt.savefig(args.out)

if not no_display:
    plt.show()
