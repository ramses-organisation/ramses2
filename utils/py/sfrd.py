#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from plot_backend import set_backend
import argparse
import ramses as ram

parser = argparse.ArgumentParser(description="Plot the cosmic star formation rate density versus redshift")
parser.add_argument("nout", help="enter output number")
parser.add_argument("--path", help="specify a path")
parser.add_argument("--out", help="output a png image")
parser.add_argument("--bin", help="specify the bin size in Gyr")
args = parser.parse_args()
no_display = set_backend()

path = args.path+"/" if args.path else "./"
bin_size = float(args.bin) if args.bin else 0.1  # in Gyr

nout = args.nout
print("Reading output number ",nout)

s=ram.rd_part(nout,path=path,prefix='star')
i=ram.rd_info(nout,path=path)

if not bool(np.any(s.birth_date < 0)):
    print("sfrd.py requires a cosmological simulation")
    exit()

# Proper time in units of 1/H0 with t=0 at a=1 (same convention as friedman in init_time.f90)
a_tab = np.logspace(-4, 0, 100001)
dtda = 1.0 / (a_tab * np.sqrt(i.omega_m/a_tab**3 + i.omega_k/a_tab**2 + i.omega_l))
t_tab = np.concatenate(([0.0], np.cumsum(0.5*(dtda[1:]+dtda[:-1])*np.diff(a_tab))))
t_tab = t_tab - t_tab[-1]

# Conversion from 1/H0 to Gyr
unit_t_Gyr = i.unit_t / i.aexp**2 / (365 * 24 * 3600 * 1e9)

tbirth = s.birth_date
tnow = i.texp
span = (tnow - np.min(tbirth)) * unit_t_Gyr
if span < bin_size:
    print("Reduce bin size, bin=", bin_size, " time span=", span)
    exit()

n_bin = int(span / bin_size)
bins = np.linspace(np.min(tbirth), tnow, n_bin+1)
dt_yr = (bins[1] - bins[0]) * unit_t_Gyr * 1e9

# Stellar mass in Msol and comoving box volume in Mpc^3
unit_m = i.unit_d * i.unit_l**3 / 2e33
vol = (i.boxlen * i.unit_l / i.aexp / 3.08568e24)**3

mass, _ = np.histogram(tbirth, weights=s.mass * unit_m, bins=bins)
sfrd = mass / dt_yr / vol

# Bin edges in redshift (increasing order)
z_edges = (1.0 / np.interp(bins, t_tab, a_tab) - 1.0)[::-1]
sfrd = sfrd[::-1]

print("Box size = %10.4e comoving Mpc  Mstar = %10.4e Msol" % (vol**(1./3.), np.sum(s.mass)*unit_m))

fig, ax = plt.subplots()
with np.errstate(divide='ignore'):
    logsfrd = np.where(sfrd > 0, np.log10(sfrd), np.nan)
ax.stairs(logsfrd, z_edges, baseline=None)
ax.set_xlim(0, 8)
ax.set_ylim(-3.2, -1.15)
ax.xaxis.set_major_locator(MultipleLocator(2))
ax.yaxis.set_major_locator(MultipleLocator(0.2))
ax.set_xlabel('Redshift')
ax.set_ylabel(r'log$_{10}$ SFRD [Msol/yr/Mpc$^3$]')

if args.out:
    plt.savefig(args.out)

if not no_display:
    plt.show()
