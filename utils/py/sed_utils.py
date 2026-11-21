#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Utilities for converting spectral energy distribution (SED) tables to
RAMSES-RT readable format, from the formats of Bruzual & Charlot (2003),
Starburst99 and BPASS (v2.2/v2.3). Routines are also included to plot the RAMSES-readable
SEDs and the derived radiation group properties.

The RAMSES format is a directory containing thee files:
    == metallicity_bins.dat ==
ascii file with the number of metallicity bins in the 1st line and then
the bin values (mass fraction), one per line
    == age_bins.dat ==
ascii, # of age bins in 1st line, then stellar pop ages (yr), one per line
    == all_seds.dat ==
f77 unformatted file, with the number of wavelength bins (followed by a
dummy 0, both int32), then wavelengths (Angstrom), then nZ*nAge arrays of
length nLambda with luminosities per bin, in units of L_sun/Angstrom/Msun,
where L_sun=3.828d33 ergs/s (same as L_sun in cooling/constants.f90),
nZ=# of metallicity bins, nAge=# of age bins, and nLambda=# of wavelength
bins.

"Public" routines:
-convertS99toRamses
-convertBC03toRamses
-convertBPASStoRamses
-plotRamsesSEDs
-plotRamsesGroups

"Private" routines, for file reading
-readRamsesSEDs
-readBC03Array

To convert and plot, either use the command line (sed_utils.py -h) or
import this module and call the individual routines. Note that in order
to plot the radiation group properties, RAMSES must first be run on the
(converted) SEDs.
"""

import numpy as np
import os
import sys
from scipy.io import FortranFile

__author__ = "Maxime Trebitsch and Joki Rosdahl"
__credits__ = ["Maxime Trebitsch and Joki Rosdahl"]
__license__ = "BSD"
__version__ = "1.0"
__maintainer__ = "Maxime Trebitsch, Joki Rosdahl"
__email__ = "maxime.trebitsch@ens-lyon.org, karl-joakim.rosdahl@univ-lyon1.fr"
__status__ = "Mint condition"

L_SUN = 3.828e33  # Solar luminosity [erg/s], must match L_sun in RAMSES


def writeLambdaCount(sedFile, nLambda):
    """Write the first record of all_seds.dat, read by RAMSES as
       two default integers: nLs, dum."""
    sedFile.write_record(np.array([nLambda, 0], dtype=np.int32))


##########################################################################
##########################################################################
def convertS99toRamses(files=None, ZBins=None, Mtot=None, outDir=None):
    """Convert Starburst99 outputs to ramses-readable SEDs.

    Parameters (user will be prompted for those if not present):
    ----------------------------------------------------------------------
    files: paths to N files containing the S99 SED ascii tables, with
           column format time (yr), wavelength (A), log total (erg/s/A),
           log stellar (erg/s/A), log nebular (erg/s/A)
    ZBins: N values giving the metallicity (mass fraction) for each S99
           SED file.
    Mtot:  Total stellar mass, in Solar masses, in the S99 model.
    outDir: Directory to put in the new ramses-rt files.

    """

    ### Prompt user for arguments if not provided-------------------------
    if outDir==None:
        print('Please write the path to where to generate new files')
        outDir=input('> ')
        if(outDir==''):
            print('No path given, nothing do to')
            return
    if files==None or ZBins==None:
        print('Please write paths and metallicities (mass fraction'
                  ,'units), separated by returns, and end with return',)
        files=[]
        ZBins=[]
        while True:
            filename=input('filename >')
            if filename == '': break
            metallicity=input('metal mass fraction >')
            if metallicity == '': break
            files.append(filename)
            ZBins.append(float(metallicity))
        print('ok checking now')
        if(files==[]):
            print('No S99 files given, nothing do to')
            return
    if Mtot==None:
        print('Please write the stellar mass used in the S99, in units of'
                  ,'Solar masses and in the format <x>e<y>.\n'
                  ,'The default in S99 is 1e6 Msun.')
        Mtot=float(input('>'))

    ### Check for consistency and initialise -----------------------------
    nFiles = len(files)                    #           Number of SED files
    nZ = len(ZBins)                        #    Number of metallicity bins
    if nZ != nFiles:
        print('The number of Z-bins must equal number of files')
        return
    ageBins = None
    lambdaBins = None

    ### Create output directory and files --------------------------------
    os.makedirs(outDir, exist_ok=True)

    sedFile = FortranFile(outDir+'/all_seds.dat','w')
    ZFile = open(outDir+'/metallicity_bins.dat', 'w')
    ageFile = open(outDir+'/age_bins.dat', 'w')

    ### Read the S99 data and convert to the ramses format----------------
    print('Reading S99 data and converting...')
    for file in files:  # Loop over SED tables for different metallicities
        print('Converting file ',file)
        data=np.loadtxt(file, skiprows=6)                       # Read SED
        if ageBins is None:                          # Some initialisation
            ageBins = sorted(set(data[:,0]))
            nAge = len(ageBins)                  # Read number of age bins
            lambdaBins = np.array(sorted(set(data[:,1])))
            nLambda = len(lambdaBins)          # Number of wavelength bins
            writeLambdaCount(sedFile, nLambda) # Write
            sedFile.write_record(lambdaBins)

        # Read luminosities and write to 'all_seds.dat'
        for iAge in range(0,nAge):
            ind_l = iAge*nLambda                 # upper and lower indexes
            ind_u = (iAge+1)*nLambda             # for one spectrum
            ages    = data[ind_l:ind_u,0]
            lambdas = data[ind_l:ind_u,1]        # read data
            lums    = data[ind_l:ind_u,2]
            # Check age bin consistency
            if any(age != ageBins[iAge] for age in ages):
                print('Age bins are not identical everywhere!!!')
                print('CANCELLING CONVERSION!!!')
                return
            # Check lambda bin consistency
            if not np.array_equal(lambdas,lambdaBins):
                print('Wavelength bins are not identical everywhere!!!')
                print('CANCELLING CONVERSION!!!')
                return

            # Convert and add to SED file
            sedFile.write_record((10.**lums)/L_SUN/Mtot)
            progress=(iAge+1)/nAge
            sys.stdout.write("\rProgress: [{0:50s}] {1:.1f}%".format(
                                '#' * int(progress * 50), progress * 100))
        print(' ')

    sedFile.close()

    ###  Write metallicities to 'metallicity_bins.dat' (units = mass frac)
    ZFile.write("%8d\n" % nZ)
    for Z in ZBins: ZFile.write("%14.6e\n" % Z)
    ZFile.close()

    ### Write stellar population ages to 'age_bins.dat' (in units of yr)
    ageFile.write("%8d\n" % nAge)
    for age in ageBins: ageFile.write("%14.6e\n" % age)
    ageFile.close()


##########################################################################
##########################################################################
def convertBC03toRamses(files=None, outDir=None):
    """Convert BC03 outputs to ramses-readable SEDs.

    Parameters (user will be prompted for those if not present):
    ----------------------------------------------------------------------
    files: list of each BC03 SED ascii file, typically named
           bc2003_xr_mxx_xxxx_ssp.ised_ASCII
    outDir: Directory to put in the new ramses-rt files.
    """
    import re

    # Prompt user for files and outDir if not provided--------------------
    if outDir==None:
        print('Please write the path to where to generate new files')
        outDir=input('> ')
        if(outDir==''):
            print('No path given, nothing do to')
            return
    if files==None:
        print('Please write the model to read',)
        files=[]
        while True:
            filename=input('filename >')
            if filename == '': break
            files.append(filename)
        print('ok checking now')
        if not len(files):
            print('No BC03 files given, nothing do to')
            return

    # Initialise ---------------------------------------------------------
    ageBins = None
    lambdaBins = None
    ZBins=[]

    # Create output directory and files ----------------------------------
    os.makedirs(outDir, exist_ok=True)

    sedFile = FortranFile(outDir+'/all_seds.dat','w')
    ZFile = open(outDir+'/metallicity_bins.dat', 'w')
    ageFile = open(outDir+'/age_bins.dat', 'w')

    # Read the BC03 files and convert to ramses format -------------------
    print('Reading BC03 files and converting...')
    for fileName in files:  # Loop SED tables for different metallicities
        print('Converting file ',fileName)
        file = open(fileName, 'r')
        ages,lastLine=readBC03Array(file) # Read age bins
        nAge=len(ages)
        if ageBins is None: ageBins=ages
        if not np.array_equal(ages,ageBins):  # check for consistency
            print('Age bins are not identical everywhere!!!')
            print('CANCELLING CONVERSION!!!')
            return
        # Read four (?) useless lines
        line=file.readline()
        line=file.readline()
        line=file.readline()
        line=file.readline()
        line=file.readline()
        ### These last three lines are identical and contain the metallicity
        ZZ, = re.search(r'Z=([0-9]+\.?[0-9]*)', line).groups()
        ZBins.append(float(ZZ))
        # Read wavelength bins
        lambdas,lastLine=readBC03Array(file,lastLineFloat=lastLine)
        if lambdaBins is None: # Write wavelengths to sed file
            lambdaBins=lambdas
            writeLambdaCount(sedFile, len(lambdaBins))
            sedFile.write_record(lambdaBins)
        if not np.array_equal(lambdas,lambdaBins):  # check for consistency
            print('Wavelength bins are not identical everywhere!!!')
            print('CANCELLING CONVERSION!!!')
            return
        # Read luminosities and write to 'all_seds.dat'
        for iAge in range(0,nAge):
            lums,lastLine = readBC03Array(file,lastLineFloat=lastLine)
            if len(lums)!=len(lambdaBins):
                print('Inconsistent number of wavelength bins in BC03')
                print('STOPPING!!')
                return
            # Read useless array
            tmp,lastLine = readBC03Array(file,lastLineFloat=lastLine)
            # Write luminosities to ramses format
            sedFile.write_record(lums)
            progress=(iAge+1)/nAge
            sys.stdout.write("\rProgress: [{0:50s}] {1:.1f}%".format(
                                '#' * int(progress * 50), progress * 100))
        file.close()
        print(' ')
        lastLine=None
    sedFile.close()

    ###  Write metallicities to 'metallicity_bins.dat' (units = mass frac)
    ZFile.write("%8d\n" % len(ZBins))
    for Z in ZBins: ZFile.write("%14.6e\n" % Z)
    ZFile.close()

    ### Write stellar population ages to 'age_bins.dat' (in units of yr)
    ageFile.write("%8d\n" % len(ageBins))
    for age in ageBins: ageFile.write("%14.6e\n" % age)
    ageFile.close()


##########################################################################
##########################################################################
def bpassMetallicity(fileName):
    """Metal mass fraction from a BPASS file name, e.g.
       spectra-bin-imf_chab300.z020.dat(.gz) -> 0.020, ...zem4... -> 1e-4.
    """
    import re
    match = re.search(r'\.z(em\d|\d{3})\.', os.path.basename(fileName))
    if match is None:
        raise ValueError('Cannot find metallicity in file name '+fileName)
    zz = match.group(1)
    if zz.startswith('em'): return 10.**(-int(zz[2:]))
    return int(zz)*1e-3


def convertBPASStoRamses(files, outDir, lambdaMax=None):
    """Convert BPASS (v2.2 or v2.3) spectra to ramses-readable SEDs.

    Parameters:
    ----------------------------------------------------------------------
    files: BPASS spectra files, one per metallicity, e.g.
           spectra-bin-imf_chab300.z*.dat (or .dat.gz). The metallicity
           is read from the file name and the bins are sorted in Z.
           Each file has a wavelength column (Angstrom) followed by 51
           columns for log(age/yr)=6.0,6.1,...,11.0, in units of
           L_sun/Angstrom for an instantaneous burst of 1e6 Msun.
    outDir: Directory to put in the new ramses-rt files.
    lambdaMax: optional maximum wavelength (Angstrom) to keep. RAMSES
           reads the full table on every MPI rank at startup, so this
           reduces the size of all_seds.dat (full BPASS: ~40 MB per
           metallicity). Must be redward of the lowest photon group
           energy (group_L0): lambda[A] = 12398.4/E[eV].
    """
    MBURST = 1e6                       # BPASS burst mass in Msun

    files = sorted(files, key=bpassMetallicity)
    ZBins = [bpassMetallicity(f) for f in files]
    if len(set(ZBins)) != len(ZBins):
        print('Duplicate metallicities in the BPASS file list!!!')
        print('CANCELLING CONVERSION!!!')
        return

    os.makedirs(outDir, exist_ok=True)
    sedFile = FortranFile(outDir+'/all_seds.dat','w')

    print('Reading BPASS files and converting...')
    lambdaBins = None
    for fileName in files:
        print('Converting file ',fileName,' Z=',bpassMetallicity(fileName))
        data = np.loadtxt(fileName)            # Handles .gz transparently
        lambdas = data[:,0]
        if lambdaMax is not None:
            keep = lambdas <= lambdaMax
            data = data[keep,:]
            lambdas = lambdas[keep]
        nAge = data.shape[1]-1
        if lambdaBins is None:
            lambdaBins = lambdas
            ageBins = 10.**(6.+0.1*np.arange(nAge))      # yr
            writeLambdaCount(sedFile, len(lambdaBins))
            sedFile.write_record(lambdaBins)
        if not np.array_equal(lambdas,lambdaBins):
            print('Wavelength bins are not identical everywhere!!!')
            print('CANCELLING CONVERSION!!!')
            sedFile.close()
            return
        if nAge != len(ageBins):
            print('Age bins are not identical everywhere!!!')
            print('CANCELLING CONVERSION!!!')
            sedFile.close()
            return
        for iAge in range(nAge):     # L_sun/A per 1e6 Msun -> per Msun
            sedFile.write_record(np.ascontiguousarray(data[:,1+iAge])/MBURST)
    sedFile.close()

    ###  Write metallicities to 'metallicity_bins.dat' (units = mass frac)
    with open(outDir+'/metallicity_bins.dat', 'w') as ZFile:
        ZFile.write("%8d\n" % len(ZBins))
        for Z in ZBins: ZFile.write("%14.6e\n" % Z)

    ### Write stellar population ages to 'age_bins.dat' (in units of yr)
    with open(outDir+'/age_bins.dat', 'w') as ageFile:
        ageFile.write("%8d\n" % len(ageBins))
        for age in ageBins: ageFile.write("%14.6e\n" % age)


##########################################################################
##########################################################################
def plotRamsesSEDs(sedDir,iZ=1,iAge=None, pdf=None):
    """Plot SED (luminosity per solar mass and unit wavelength versus
       wavelength), read from ramses format.

    Parameters:
    ----------------------------------------------------------------------
    sedDir: Directory containing the ramses-readable SED tables
    iZ: metallicity index for which to plot (1 to number of Z bins)
    iAge: age index(es) for which to plot (1 to number of age bins). If
          omitted, curves are plotted for a selection of stellar pop ages.
    """
    import matplotlib.pyplot as p
    import matplotlib
    from matplotlib.backends.backend_pdf import PdfPages

    ### Read the SEDs ----------------------------------------------------
    seds = readRamsesSEDs(sedDir)
    ZBins=seds['ZBins']
    ageBins_myr=np.array(seds['ageBins'])/1e6
    lambdaBins = seds['lambdaBins']
    spectra = seds['spectra']

    ### Set up the stellar population ages and corresponding colors ------
    if iAge==None: # Select some ages
        iAge=set()
        ages_myr=[0, 1, 3, 5, 10, 50, 100, 250, 500, 1000, 10000]
        for age_myr in ages_myr:
            age_ind=np.argmin(np.abs(age_myr-ageBins_myr))
            iAge.add(age_ind+1)
        iAge=sorted(iAge)
    if np.isscalar(iAge): iAge=[iAge]   # if scalar argument was passed in
    cm = matplotlib.colormaps['viridis'] # Curve colors
    NUM_COLORS=len(iAge)
    p.gca().set_prop_cycle('color'
                        ,[cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])

    ### Set up the plot --------------------------------------------------
    matplotlib.rcParams.update({'font.size': 16})
    xrg=[8e1,1.7e6] # x-axis range in Angstrom
    yrg=[1e-7,1e1] # y-axis range in L_sun/M_sun/Angstrom
    p.title('Ramses SED with metal mass fraction of '+str(ZBins[iZ-1]))
    p.xlabel(r'$\lambda$ [$\mathrm{\AA}$]')
    p.xscale('log')
    p.ylabel(r'J$_{\lambda}$ [L$_{\odot}$/M$_{\odot}$/$\mathrm{\AA}$]')
    p.yscale('log')
    p.axis([xrg[0],xrg[1],yrg[0],yrg[1]])

    ### Plot the spectra -------------------------------------------------
    for i in iAge:
        label=str(round(ageBins_myr[i-1]))+' Myr'
        p.plot(lambdaBins,spectra[:,i-1,iZ-1],label=label)

    ### Show legend for the age(s) ---------------------------------------
    leg=p.legend(frameon=False,prop={'size':10})

    if pdf != None:
        fig=p.figure(1)
        pp=PdfPages(pdf)
        pp.savefig(fig)
        pp.close()
        p.close()
    else: p.show()

##########################################################################
##########################################################################
def plotRamsesGroups(groupFiles=['SEDtable1.list'],isEgy=False, pdf=None):
    """Plot photon group properties and stellar luminosities, as generated
       by ramses-rt from SED tables, and found in SEDtableX.list files
       in the ramses run directory.
       The format of each (ascii) file is
       age (Gyr), Z, luminosity, time-integrated luminosity, group energy
       (eV), cross sections (cm2). The luminosity is in units of
       photons/s/Msun by default, but erg/s/Msun if the ramses-run was
       made with SED_isEgy=.true.

    Parameters:
    ----------------------------------------------------------------------
    groupFiles: Files containing group properties (one per group)
    isEgy: If stellar luminosity was run in energy mode (as opposed to the
           default photon number conserving mode -- see ramses namelist
           parameter of the same name).
    """
    import matplotlib.pyplot as p
    import matplotlib
    from matplotlib.backends.backend_pdf import PdfPages
    matplotlib.rcParams['figure.figsize'] = 13, 10

    p.subplots_adjust(left=0.07,right=0.97,hspace=0.04, top=0.95
                      ,bottom=0.06, wspace=0.04)
    xrg=[2e-1,2e4]
    lumrg=[1e39,1e52]
    lumLabel=r'Lum [#/s/M$_{\odot}$]'
    ilumLabel=r'$\Pi$ [#/M$_{\odot}$]'
    ilumrg=[1e53,1e65]
    lumConv=1.
    egyrg=[1e-1,1e3]
    sigrg=[1e-20,1e-17]
    if isEgy:
        lumrg=[1e30,1e37]
        lumLabel=r'Lum [erg/s/M$_{\odot}$]'
        lumConv=1.60218e-12                                   # ergs in eV
        ilumrg=[1e45,1e53]
        ilumLabel=r'$\Pi$ [erg/M$_{\odot}$]'
    ### Read the files ---------------------------------------------------
    nGroups = len(groupFiles)
    for i,fileName in enumerate(groupFiles):
        file = open(fileName, 'r')
        headerLine=file.readline().split()
        nAge=int(headerLine[0])
        nZ=int(headerLine[1])
        data=np.loadtxt(file)
        # Plot the luminosity
        ax1=p.subplot2grid((4,nGroups),(0,i))
        p.setp(ax1.get_xticklabels(), visible=False)
        if i>0: p.setp(ax1.get_yticklabels(), visible=False)
        if i==0: ax1.set_ylabel(lumLabel)
        ax1.set_ylim([lumrg[0],lumrg[1]])
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        p.title('Group '+str(i+1))
        for iZ in range(0,nZ):
            label='Z='+str(data[iZ*nAge,1])             # Curve metallicity
            ax1.plot(data[iZ*nAge:(iZ+1)*nAge,0]*1e3
                        ,data[iZ*nAge:(iZ+1)*nAge,2]*lumConv
                        ,label=label)
        # Plot the integrated luminosity
        ax=p.subplot2grid((4,nGroups),(1,i), sharex=ax1)
        p.setp(ax.get_xticklabels(), visible=False)
        if i>0: p.setp(ax.get_yticklabels(), visible=False)
        if i==0: ax.set_ylabel(ilumLabel)
        ax.set_ylim([ilumrg[0],ilumrg[1]])
        ax.set_yscale('log')
        for iZ in range(0,nZ):
            ax.plot(data[iZ*nAge:(iZ+1)*nAge,0]*1e3
                        ,data[iZ*nAge:(iZ+1)*nAge,3]*lumConv)
        # Plot the photon group energy
        ax=p.subplot2grid((4,nGroups),(2,i), sharex=ax1)
        p.setp(ax.get_xticklabels(), visible=False)
        if i>0: p.setp(ax.get_yticklabels(), visible=False)
        if i==0: ax.set_ylabel(r'egy [eV]')
        ax.set_ylim([egyrg[0],egyrg[1]])
        ax.set_yscale('log')
        for iZ in range(0,nZ):
            label='Z='+str(data[iZ*nAge,1])             # Curve metallicity
            ax.plot(data[iZ*nAge:(iZ+1)*nAge,0]*1e3
                        ,data[iZ*nAge:(iZ+1)*nAge,4],label=label)
        if(i==0): leg=ax.legend(frameon=False,prop={'size':10})
        # Plot cross sections
        ax=p.subplot2grid((4,nGroups),(3,i), sharex=ax1)
        if i>0: p.setp(ax.get_yticklabels(), visible=False)
        ax.set_xlabel(r'Time [Myr]')
        if i==0: ax.set_ylabel(r'$\sigma$ [cm$^2$]')
        ax.set_ylim([sigrg[0],sigrg[1]])
        ax.set_yscale('log')
        for iZ in range(0,nZ):
            ax.plot(data[iZ*nAge:(iZ+1)*nAge,0]*1e3
                        ,data[iZ*nAge:(iZ+1)*nAge,5])
            ax.plot(data[iZ*nAge:(iZ+1)*nAge,0]*1e3
                        ,data[iZ*nAge:(iZ+1)*nAge,7]
                        ,ls='--')
            ax.plot(data[iZ*nAge:(iZ+1)*nAge,0]*1e3
                        ,data[iZ*nAge:(iZ+1)*nAge,9]
                        ,ls=':')
        ax.set_xlim(xrg[0],xrg[1])

    if pdf != None:
        fig=p.figure(1)
        pp=PdfPages(pdf)
        pp.savefig(fig)
        pp.close()
        p.close()
    else: p.show()


##########################################################################
##########################################################################
def readRamsesSEDs(sedDir):
    """Read SED in ramses format and return

    Parameters:
    ----------------------------------------------------------------------
    sedDir: Directory containing the SED tables
    """
    # Read metallicity bins
    ZFile = open(sedDir+'/metallicity_bins.dat', 'r')
    nZ = int(ZFile.readline())
    ZBins = []
    for Z in range(0,nZ): ZBins.append(float(ZFile.readline()))
    ZFile.close()

    # Read age bins
    ageFile = open(sedDir+'/age_bins.dat', 'r')
    nAge = int(ageFile.readline())
    ageBins = []
    for age in range(0,nAge): ageBins.append(float(ageFile.readline()))
    ageFile.close()

    # Read wavelength bins and spectra
    sedFile = FortranFile(sedDir+'/all_seds.dat','r')
    nLambda = sedFile.read_ints()[0]
    lambdaBins = sedFile.read_reals()
    spectra = np.empty([nLambda,nAge,nZ])
    for iZ in range(0,nZ):
        for iAge in range(0,nAge):
            spectrum = sedFile.read_reals()
            spectra[:,iAge,iZ] = spectrum

    return {'ZBins':ZBins, 'ageBins':ageBins, 'lambdaBins':lambdaBins
                ,'spectra':spectra}


##########################################################################
##########################################################################
def readBC03Array(file, lastLineFloat=None):
    """Read a record from bc03 ascii file. The record starts with the
       number of elements N and is followed by N numbers. The record may
       or may not start within a line, i.e. a line need not necessarily
       start with a record.
    Parameters:
    ----------------------------------------------------------------------
    file: handle on open bc03 ascii file
    lastLineFloat: still open line from last line read, in case of a
                   record starting mid-line.
    Returns array, lastLine, where:
    ----------------------------------------------------------------------
    array = The array values read from the file
    lastLine = The remainder of the last line read (in floating format),
               for continued reading of the file
    """
    if lastLineFloat==None or len(lastLineFloat)==0:
        ### Nothing in last line, so read next line
        line=file.readline()
        lineStr = line.split()
        lastLineFloat = [float(x) for x in lineStr]
    ### Read array 'header' (i.e. number of elements)
    arrayCount = int(lastLineFloat[0])          # Length of returned array
    array=np.empty(arrayCount)                  #     Initialise the array
    lastLineFloat=lastLineFloat[1:len(lastLineFloat)]
    iA=0 # Running array index
    while True: # Read numbers until array is full
        for iL in range(0,len(lastLineFloat)):  #     Loop numbers in line
            array[iA]=lastLineFloat[iL]
            iA=iA+1
            if iA >= arrayCount:                #  Array is full so return
                return array,lastLineFloat[iL+1:]
        line=file.readline()   # Went through the line so get the next one
        lineStr = line.split()
        lastLineFloat = [float(x) for x in lineStr]


##########################################################################
##########################################################################
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description='Convert BC03, Starburst99 or BPASS SEDs to RAMSES-RT format'
                    ' and plot them.')
    sub = parser.add_subparsers(dest='cmd', required=True)

    pb = sub.add_parser('bc03', help='convert BC03 *.ised_ASCII files')
    pb.add_argument('outdir', help='output SED directory')
    pb.add_argument('files', nargs='+',
                    help='one BC03 ASCII file per metallicity')

    pB = sub.add_parser('bpass', help='convert BPASS spectra-*.dat files')
    pB.add_argument('outdir', help='output SED directory')
    pB.add_argument('files', nargs='+',
                    help='BPASS spectra files (.dat or .dat.gz), one per'
                         ' metallicity, Z read from the file name')
    pB.add_argument('--lmax', type=float, default=None,
                    help='keep wavelengths <= LMAX Angstrom (default all)')

    ps = sub.add_parser('s99', help='convert Starburst99 spectrum files')
    ps.add_argument('outdir', help='output SED directory')
    ps.add_argument('files', nargs='+',
                    help='one S99 spectrum file per metallicity')
    ps.add_argument('--Z', type=float, nargs='+', required=True,
                    help='metal mass fraction of each file')
    ps.add_argument('--mtot', type=float, default=1e6,
                    help='S99 total stellar mass in Msun (default 1e6)')

    pp = sub.add_parser('plot', help='plot SEDs from a RAMSES SED directory')
    pp.add_argument('seddir')
    pp.add_argument('--iZ', type=int, default=1,
                    help='metallicity index, starting at 1 (default 1)')
    pp.add_argument('--pdf', default=None, help='save to this PDF file')

    pg = sub.add_parser('groups',
                        help='plot SEDtableX.list files written by RAMSES')
    pg.add_argument('files', nargs='+')
    pg.add_argument('--pdf', default=None, help='save to this PDF file')

    args = parser.parse_args()
    if args.cmd == 'bc03':
        convertBC03toRamses(files=args.files, outDir=args.outdir)
    elif args.cmd == 'bpass':
        convertBPASStoRamses(args.files, args.outdir, lambdaMax=args.lmax)
    elif args.cmd == 's99':
        convertS99toRamses(files=args.files, ZBins=args.Z,
                           Mtot=args.mtot, outDir=args.outdir)
    elif args.cmd == 'plot':
        plotRamsesSEDs(args.seddir, iZ=args.iZ, pdf=args.pdf)
    elif args.cmd == 'groups':
        plotRamsesGroups(groupFiles=args.files, pdf=args.pdf)
