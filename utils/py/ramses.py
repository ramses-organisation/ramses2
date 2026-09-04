import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from scipy.io import FortranFile
from astropy.io import ascii
import os
import re
import glob

import time

class Cool:
    """
    This is the class for RAMSES cooling table.
    """
    def __init__(self,n1,n2):
        """
        This function initialize the cooling table. 
        Args:
            n1: number of points for the gas density axis
            n2: number of points for the gas temperature axis
        """
        self.n1 = n1
        self.n2 = n2
        self.nH = np.zeros([n1])
        self.T2 = np.zeros([n2])
        self.cool = np.zeros([n1,n2])
        self.heat = np.zeros([n1,n2])
        self.spec = np.zeros([n1,n2,6])
        self.xion = np.zeros([n1,n2])

def clean(dat,n1,n2):
    dat = np.array(dat)
    dat = dat.reshape(n2, n1)
    return dat

def clean_spec(dat,n1,n2):
    dat = np.array(dat)
    dat = dat.reshape(6, n2, n1)
    return dat

def rd_cool(filename):
    """This function reads a RAMSES cooling table file (unformatted Fortran binary) 
    and store it in a cooling object.

    Args:
        filename: the complete path (including the name) of the cooling table file.

    Returns:
        A cooling table (Cool) object.

    Authors: Romain Teyssier (Princeton University, October 2022)
    """
    with FortranFile(filename, 'r') as f:
        n1, n2 = f.read_ints('i')
        c = Cool(n1,n2)
        nH = f.read_reals('f8')
        T2 = f.read_reals('f8')
        cool = f.read_reals('f8')
        heat = f.read_reals('f8')
        cool_com = f.read_reals('f8')
        heat_com = f.read_reals('f8')
        metal = f.read_reals('f8')
        cool_prime = f.read_reals('f8')
        heat_prime = f.read_reals('f8')
        cool_com_prime = f.read_reals('f8')
        heat_com_prime = f.read_reals('f8')
        metal_prime = f.read_reals('f8')
        mu = f.read_reals('f8')
        n_spec = f.read_reals('f8')
        c.nH = nH
        c.T2 = T2
        c.cool = clean(cool,n1,n2)
        c.heat = clean(heat,n1,n2)
        c.metal = clean(metal,n1,n2)
        c.spec = clean_spec(n_spec,n1,n2)
        c.xion = c.spec[0]
        for i in range(0,n2):
            c.xion[i,:] = c.spec[0,i,:] - c.nH
        return c

def test_cool(filename):
    """This function reads a binary file produced by the non-equilibrium chemistry solver
    of the RAMSES-RT code and make a plot.

    Args:
        filename: the complete path (including the name) of the cooling table file.

    Authors: Romain Teyssier (Princeton University, March 2025)
    """
    np.set_printoptions(linewidth=120)
    n_d = np.fromfile(filename,dtype=np.int32,count=1,offset=0)[0]
    n_t = np.fromfile(filename,dtype=np.int32,count=1,offset=4)[0]
    n_x = np.fromfile(filename,dtype=np.int32,count=1,offset=8)[0]
    n_time = np.fromfile(filename,dtype=np.int32,count=1,offset=12)[0]
    n_vec = np.fromfile(filename,dtype=np.int32,count=1,offset=16)[0]
    skip = 20
    nh = np.fromfile(filename,dtype=np.float64,count=n_d,offset=skip)
    skip = skip + n_d*8
    times = np.fromfile(filename,dtype=np.float64,count=n_time,offset=skip)
    skip = skip + n_time*8
    size = n_d*n_t*n_x*n_time*n_vec
    data = np.fromfile(filename,dtype=np.float64,count=size,offset=skip)
    data = np.reshape(data,(n_d,n_t,n_x,n_time,n_vec),order='F')
    print('shape of data =',n_d,n_t,n_x,n_time,n_vec)
    print('n_H[H/cc] =',nh)
    print('T_ini[K] =',data[0,:,0,0,0])

    figure, axis = plt.subplots(n_d, n_t, figsize=(15, 15))
    plt.subplots_adjust(wspace=0,hspace=0)
    i=0
    ix=n_x-1
    for id in range(0,n_d):
        for it in range(0,n_t):
            axis[id,it].plot([],color='r',label='xHII')
            axis[id,it].plot([],color='g',label='xHeII')
            axis[id,it].plot([],color='b',label='xHeIII')
            for ix in range(0,n_x):
                axis[id,it].plot(np.log10(times),data[id,it,ix,:,1],color='r')
                axis[id,it].plot(np.log10(times),data[id,it,ix,:,2],color='g')
                axis[id,it].plot(np.log10(times),data[id,it,ix,:,3],color='b')

            axis[id,it].set_xlim([-3.2,3.2])
            if it>0:
                axis[id,it].set_yticklabels([])
            else:
                axis[id,it].set_ylabel("fraction")
            if id<n_d-1:
                axis[id,it].set_ylim([0.0001,1.2])
            else:
                axis[id,it].set_ylim([0,1.2])
                axis[id,it].set_xlabel("log time [Myr]")
            if it==n_t-1 and id==n_d-1:
                axis[id,it].legend(loc="lower right")
            axis[id,it].set_title('log nH = ' + str(np.log10(nh[id])) + 
                                  ' log T = ' + str(np.log10(data[id,it,ix,0,0])), y=0.9, va="top")
            i=i+1

class Map:
    """This class defines a map object.
    """
    def __init__(self,nx,ny):
        """This function initalize a map object.
        
        Args:
            nx: number of pixels in the x direction
            ny: number of pixels in the y direction
        """
        self.nx = nx
        self.ny = ny
        self.data = np.zeros([nx,ny])

def rd_map(filename):
    """This function reads a RAMSES map file (unformatted Fortran binary)
    as produced by the RAMSES utilities amr2map or part2map and store it in a map object.

    Args:
        filename: the complete path (including the name) of the map file.

    Returns:
        A map (class Map) object.

    Example:
        import miniramses as ram
        map = ram.rd_map("dens.map")
        plt.imshow(map.data,origin="lower")

    Authors: Romain Teyssier (Princeton University, October 2022)
    """
    with FortranFile(filename, 'r') as f:
        t, dx, dy, dz = f.read_reals('f8')
        nx, ny = f.read_ints('i')
        dat = f.read_reals('f4')
    
    dat = np.array(dat)
    dat = dat.reshape(ny, nx)
    m = Map(nx,ny)
    m.data = dat.T
    m.time = t
    m.nx = nx
    m.ny = ny
    
    return m

class Histo:
    """This class defines a histogram object.
    """
    def __init__(self,nx,ny):
        """This function initalize a histogram object.

        Args:
            nx: number of pixels in the x direction
            ny: number of pixels in the y direction
        """
        self.nx = nx
        self.ny = ny
        self.h = np.zeros([nx,ny])

def rd_histo(filename):
    """This function reads a RAMSES histogram file (unformatted Fortran binary)
    as produced by the RAMSES utilities histo and store it in a Histo object.

    Args:
        filename: the complete path (including the name) of the histo file.

    Returns:
        A histogram (class Histo) object.

    Example:
        import miniramses as ram
        h = ram.rd_histo("histo.dat")
        plt.imshow(h.data,origin="lower")
    
    Authors: Romain Teyssier (Princeton University, October 2022)
    """
    with FortranFile(filename, 'r') as f:
        nx, ny = f.read_ints('i')
        dat = f.read_reals('f4')
        lxmin, lxmax = f.read_reals('f8')
        lymin, lymax = f.read_reals('f8')

    dat = np.array(dat)
    dat = dat.reshape(ny, nx)
    h = Histo(nx,ny)
    h.data = dat
    h.nx = nx
    h.ny = ny
    h.lxmin = lxmin
    h.lxmax = lxmax
    h.lymin = lymin
    h.lymax = lymax

    return h

class Part:
    def __init__(self,nnp,nndim,star=False,sink=False,tree=False,peak=False):
        self.npart = nnp
        self.ndim = nndim
        self.pos = np.zeros([nndim,nnp])
        self.vel = np.zeros([nndim,nnp])
        self.mass = np.zeros([nnp])
        self.level = np.zeros([nnp])
        self.birth_id = np.zeros([nnp])
        if(star):
            self.metallicity = np.zeros([nnp])
            self.birth_date = np.zeros([nnp])
        if(tree):
            self.birth_date = np.zeros([nnp])
            self.merging_date = np.zeros([nnp])
            self.merging_id = np.zeros([nnp])
            self.tracking_id = np.zeros([nnp])
        if(sink):
            self.angmom = np.zeros([nndim,nnp])
            self.accel = np.zeros([nndim,nnp])
            self.birth_date = np.zeros([nnp])
        if(peak):
            self.halo_id = np.zeros([nnp],dtype=np.int32)
            self.peak_id = np.zeros([nnp],dtype=np.int32)
            
def rd_part(nout,**kwargs):
    """This function reads a RAMSES particle file (unformatted Fortran binary) 
    as produced by the RAMSES code in the snapshot directory output_00* 
    and store it in a variable containing all the particle information (Part object).

    Args:
        nout: the RAMSES snapshot number. For example output_000012 corresponds to nout=12.

    Optional args:

        path:   a string containing the relative path of the output folder.
        prefix: a string describing the particle type. It corresponds to the filenames uswd in
                the output folder. prefix='part' by default. It takes the followinf  possible values:
                prefix='part', 'star', 'tree', 'sink'
        center: a numpy array containing the coordinates of the center of the sphere restricting
                the region to read in data.
        radius: the radius of the sphere restricting the region to read in data.
        peak:   a logical True or False to read also the peak information of the clump finder.

    Returns:
        A variable p (class Part) object defined as:
            p.npart: number of particles
            p.ndim: number of space dimensions
            p.pos: coordinates of the particles. p.pos[0] gives the x coordinate as a numpy array.
            p.vel: velocities of the particles. p.vel[0] gives the x-component as a numpy array.
            p.mass: array containing the particle masses
        The number of fields depends on the particle type defined by prefix.

    Example:
        import miniramses as ram
        p = ram.rd_part(12,center=[0.5,0.5,0.5],radius=0.1)
        print(np.max(p.pos[0]))
    
    Authors: Romain Teyssier (Princeton University, October 2022)
    """
    
    prefix = kwargs.get("prefix","part")
    backup = kwargs.get("backup",False)
    center = kwargs.get("center")
    radius = kwargs.get("radius")
    path = kwargs.get("path","./")
    peak = kwargs.get("peak",False)
    silent = kwargs.get("silent",False)

    car1 = str(nout).zfill(5)
    i = rd_info(nout,path=path,backup=backup)
    ncpu = i.ncpu
    ndim = i.ndim
    levelmin = i.levelmin
    nlevelmax = i.nlevelmax
    boxlen = i.boxlen

    #if ( not (center is None)  and not (radius is None) ):
    #    info = rd_info(nout)
    #    cpulist = get_cpu_list(info,**kwargs)
    #    print("Will open only",len(cpulist),"files")
    #else:
    #    cpulist = range(1,ncpu+1)
    cpulist = range(1,ncpu+1)

    star = False
    sink = False
    tree = False
    if(prefix=="star"):
        star=True
    if(prefix=="sink"):
        sink=True
    if(prefix=="tree"):
        tree=True
    prefix2="/"+prefix+"."
    npart = 0
    for icpu in cpulist:
        car2 = str(icpu).zfill(5)
        if(backup):
            filename = path+"/backup_"+car1+prefix2+car2
        else:
            filename = path+"/output_"+car1+prefix2+car2

        npart2 = np.fromfile(filename,dtype=np.int32,count=1,offset=4)[0]
        npart = npart + npart2

    if silent==False:
        txt = "Found "+str(npart)+" particles"
        print(txt)

    p = Part(npart,ndim,star,sink,tree,peak)
    p.npart = npart
    p.ndim = ndim

    ipart = 0
    for icpu in cpulist:
        car2 = str(icpu).zfill(5)
        if(backup):
            filename = path+"/backup_"+car1+prefix2+car2
        else:
            filename = path+"/output_"+car1+prefix2+car2

        npart2 = np.fromfile(filename,dtype=np.int32,count=1,offset=4)[0]

        offset = np.int64(8) # prevent overflow
        
        # read particle positions
        for idim in range(0,ndim):
            if(backup):
                xp = np.fromfile(filename,dtype=np.float64,count=npart2,offset=offset)
                offset = offset + npart2*8
            else:
                xp = np.fromfile(filename,dtype=np.float32,count=npart2,offset=offset)
                offset = offset + npart2*4

            p.pos[idim,ipart:ipart+npart2] = xp

        # read particle velocities
        for idim in range(0,ndim):
            if(backup):
                xp = np.fromfile(filename,dtype=np.float64,count=npart2,offset=offset)
                offset = offset + npart2*8
            else:
                xp = np.fromfile(filename,dtype=np.float32,count=npart2,offset=offset)
                offset = offset + npart2*4

            p.vel[idim,ipart:ipart+npart2] = xp

        # read particle masses
        if(backup):
            xp = np.fromfile(filename,dtype=np.float64,count=npart2,offset=offset)
            offset = offset + npart2*8
        else:
            xp = np.fromfile(filename,dtype=np.float32,count=npart2,offset=offset)
            offset = offset + npart2*4

        p.mass[ipart:ipart+npart2] = xp

        if(star):
            # read particle metallicities
            if(backup):
                xp = np.fromfile(filename,dtype=np.float64,count=npart2,offset=offset)
                offset = offset + npart2*8
            else:
                xp = np.fromfile(filename,dtype=np.float32,count=npart2,offset=offset)
                offset = offset + npart2*4

            p.metallicity[ipart:ipart+npart2] = xp

            # read particle birth times
            if(backup):
                xp = np.fromfile(filename,dtype=np.float64,count=npart2,offset=offset)
                offset = offset + npart2*8
            else:
                xp = np.fromfile(filename,dtype=np.float32,count=npart2,offset=offset)
                offset = offset + npart2*4

            p.birth_date[ipart:ipart+npart2] = xp

        if(sink):
            # read particle accelerations
            for idim in range(0,ndim):
                if(backup):
                    xp = np.fromfile(filename,dtype=np.float64,count=npart2,offset=offset)
                    offset = offset + npart2*8
                else:
                    xp = np.fromfile(filename,dtype=np.float32,count=npart2,offset=offset)
                    offset = offset + npart2*4

                p.accel[idim,ipart:ipart+npart2] = xp

            # read particle angular momentum
            for idim in range(0,ndim):
                if(backup):
                    xp = np.fromfile(filename,dtype=np.float64,count=npart2,offset=offset)
                    offset = offset + npart2*8
                else:
                    xp = np.fromfile(filename,dtype=np.float32,count=npart2,offset=offset)
                    offset = offset + npart2*4

                p.angmom[idim,ipart:ipart+npart2] = xp

            # read particle birth times
            if(backup):
                xp = np.fromfile(filename,dtype=np.float64,count=npart2,offset=offset)
                offset = offset + npart2*8
            else:
                xp = np.fromfile(filename,dtype=np.float32,count=npart2,offset=offset)
                offset = offset + npart2*4

            p.birth_date[ipart:ipart+npart2] = xp

        if(tree):
            # read particle birth times
            if(backup):
                xp = np.fromfile(filename,dtype=np.float64,count=npart2,offset=offset)
                offset = offset + npart2*8
            else:
                xp = np.fromfile(filename,dtype=np.float32,count=npart2,offset=offset)
                offset = offset + npart2*4

            p.birth_date[ipart:ipart+npart2] = xp

            # read particle merging times
            if(backup):
                xp = np.fromfile(filename,dtype=np.float64,count=npart2,offset=offset)
                offset = offset + npart2*8
            else:
                xp = np.fromfile(filename,dtype=np.float32,count=npart2,offset=offset)
                offset = offset + npart2*4

            p.merging_date[ipart:ipart+npart2] = xp

        # read particle level
        xp = np.fromfile(filename,dtype=np.int32,count=npart2,offset=offset)
        offset = offset + npart2*4

        p.level[ipart:ipart+npart2] = xp

        # read particle id
        xp = np.fromfile(filename,dtype=np.int32,count=npart2,offset=offset)
        offset = offset + npart2*4

        p.birth_id[ipart:ipart+npart2] = xp

        # read particle merging id and tracking id
        if(tree):
            xp = np.fromfile(filename,dtype=np.int32,count=npart2,offset=offset)
            offset = offset + npart2*4

            p.merging_id[ipart:ipart+npart2] = xp

            xp = np.fromfile(filename,dtype=np.int32,count=npart2,offset=offset)
            offset = offset + npart2*4

            p.tracking_id[ipart:ipart+npart2] = xp

        ipart = ipart + npart2

    if(peak):
        prefix2="/peak_part."        
        if(star):
            prefix2="/peak_star."
        if(sink):
            prefix2="/peak_sink."
        if(tree):
            prefix2="/peak_tree."
        ipart = 0
        for icpu in cpulist:
            car2 = str(icpu).zfill(5)
            filename = path+"/output_"+car1+prefix2+car2
            npart2 = np.fromfile(filename,dtype=np.int32,count=1,offset=4)[0]
            offset = np.int64(8)

            # read particle halo id
            hid = np.fromfile(filename,dtype=np.int32,count=npart2,offset=offset)
            p.halo_id[ipart:ipart+npart2] = hid
            offset = offset + npart2*4

            # read particle peak id
            pid = np.fromfile(filename,dtype=np.int32,count=npart2,offset=offset)
            p.peak_id[ipart:ipart+npart2] = pid
            offset = offset + npart2*4

            ipart = ipart + npart2

    # Filtering particles within the input sphere
    if ( not (center is None)  and not (radius is None) ):

        # Periodic boundaries
        for idim in range(0,ndim):
            xp = p.pos[idim]-center[idim]
            xp[xp>boxlen/2]=xp[xp>boxlen/2]-boxlen
            xp[xp<-boxlen/2]=xp[xp<-boxlen/2]+boxlen
            p.pos[idim] = xp+center[idim]
        if ndim==1:
            r = np.sqrt((p.pos[0]-center[0])**2)
        if ndim==2:
            r = np.sqrt((p.pos[0]-center[0])**2+(p.pos[1]-center[1])**2)
        if ndim==3:
            r = np.sqrt((p.pos[0]-center[0])**2+(p.pos[1]-center[1])**2+(p.pos[2]-center[2])**2)
        p.npart = np.count_nonzero(r < radius)
        p.mass = p.mass[r < radius]
        p.pos = p.pos[:,r < radius]
        p.vel = p.vel[:,r < radius]
        p.level = p.level[r < radius]
        p.birth_id = p.birth_id[r < radius]
        if(star):
            p.metallicity = p.metallicity[r < radius]
            p.birth_date = p.birth_date[r < radius]
        if(sink):
            p.accel = p.accel[:,r < radius]
            p.birth_date = p.birth_date[r < radius]
            p.angmom = p.angmom[r < radius]
        if(tree):
            p.birth_date = p.birth_date[r < radius]
            p.merging_date = p.merging_date[r < radius]
            p.merging_id = p.merging_id[r < radius]
            p.tracking_id = p.tracking_id[r < radius]
        if(peak):
            p.peak_id = p.peak_id[r < radius]
            p.halo_id = p.halo_id[r < radius]
        if(silent==False):
            txt = "Kept "+str(p.npart)+" particles"
            print(txt)

    return p

class Level:
    def __init__(self,nndim):
        self.level = 0
        self.ngrid = 0
        self.ndim = nndim
        self.xg = np.empty(shape=(nndim,0))
        self.refined = np.empty(shape=(2**nndim,0),dtype=bool)

def rd_amr(nout,**kwargs):

    backup = kwargs.get("backup",False)
    center = kwargs.get("center")
    radius = kwargs.get("radius")
    path = kwargs.get("path","./")

    car1 = str(nout).zfill(5)
    i = rd_info(nout,path=path,backup=backup)
    ncpu = i.ncpu
    ndim = i.ndim
    levelmin = i.levelmin
    nlevelmax = i.nlevelmax

    txt = "ncpu="+str(ncpu)+" ndim="+str(ndim)+" nlevelmax="+str(nlevelmax)
    print(txt)
    print("Time=",i.texp)
    print("Reading grid data...")

    #if ( not (center is None)  and not (radius is None) ):
    #    info = rd_info(nout)
    #    cpulist = get_cpu_list(info,**kwargs)
    #    print("Will open only",len(cpulist),"files")
    #else:
    #    cpulist = range(1,ncpu+1)
    cpulist = range(1,ncpu+1)

    amr=[]
    for ilevel in range(0,nlevelmax):
        amr.append(Level(ndim))
        
    amr[0].boxlen = i.boxlen
    
    numbl = np.zeros([nlevelmax,ncpu],dtype=np.int32)
    
    # Reading and computing total AMR grids count
    for icpu in cpulist:

        car1 = str(nout).zfill(5)
        car2 = str(icpu).zfill(5)

        if(backup):
            filename = path+"/backup_"+car1+"/amr."+car2
        else:
            filename = path+"/output_"+car1+"/amr."+car2

        skip = 12
        for ilevel in range(levelmin-1,nlevelmax):
            offset = skip+4*(ilevel+1-levelmin)
            numbl[ilevel,icpu-1] = np.fromfile(filename,dtype=np.int32,count=1,offset=offset)[0]
            amr[ilevel].ngrid = amr[ilevel].ngrid + numbl[ilevel,icpu-1]

    # Allocating memory
    for ilevel in range(0,nlevelmax):
        amr[ilevel].xg = np.zeros([ndim,amr[ilevel].ngrid],dtype=float)
        amr[ilevel].refined = np.zeros([2**ndim,amr[ilevel].ngrid],dtype=bool)

    iskip = np.zeros(nlevelmax, dtype=int)
    nvar = ndim+1

    # Reading and storing data
    for icpu in cpulist:

        car1 = str(nout).zfill(5)
        car2 = str(icpu).zfill(5)
        if(backup):
            filename = path+"/backup_"+car1+"/amr."+car2
        else:
            filename = path+"/output_"+car1+"/amr."+car2

        offset = 12 + 4*(nlevelmax+1-levelmin)
        for ilevel in range(levelmin-1,nlevelmax):
            ncache = numbl[ilevel,icpu-1]

            transfer = np.fromfile(filename,dtype=np.int32,count=nvar*ncache,offset=offset)
            transfer = np.reshape(transfer,(ncache,nvar))
            transfer = np.transpose(transfer)

            # Store grid Cartesian index
            for idim in range(0,ndim):
                amr[ilevel].xg[idim,iskip[ilevel]:iskip[ilevel]+ncache] = transfer[idim]
                
            # Store cell refinement map
            refined_int = transfer[ndim]
            for ind in range(0,2**ndim):
                amr[ilevel].refined[ind,iskip[ilevel]:iskip[ilevel]+ncache] = (refined_int >> ind) & 1
            
            offset = offset + ncache*nvar*4
            iskip[ilevel] = iskip[ilevel] + ncache
            
    return amr

class Hydro:
    def __init__(self,nndim,nnvar):
        self.level = 0
        self.ngrid = 0
        self.ndim = nndim
        self.nvar = nnvar
        self.u = np.empty(shape=(nnvar,2**nndim,0))

def rd_hydro(nout,**kwargs):

    prefix = kwargs.get("prefix","hydro")
    backup = kwargs.get("backup",False)
    center = kwargs.get("center")
    radius = kwargs.get("radius")
    path = kwargs.get("path","./")

    car1 = str(nout).zfill(5)
    i = rd_info(nout,path=path,backup=backup)
    ncpu = i.ncpu
    ndim = i.ndim
    levelmin = i.levelmin
    nlevelmax = i.nlevelmax

    #if ( not (center is None)  and not (radius is None) ):
    #    info = rd_info(nout)
    #    cpulist = get_cpu_list(info,**kwargs)
    #    print("Will open only",len(cpulist),"files")
    #else:
    #    cpulist = range(1,ncpu+1)
    cpulist = range(1,ncpu+1)

    # Get number of hydro variables
    car1 = str(nout).zfill(5)
    if(backup):
        filename = path+"/backup_"+car1+"/"+prefix+".00001"
    else:
        filename = path+"/output_"+car1+"/"+prefix+".00001"

    nvar = np.fromfile(filename,dtype=np.int32,count=1,offset=4)[0]
    
    txt = "Found nvar="+str(nvar)
    print(txt)
    print("Reading "+prefix+" data...")

    hydro=[]
    for ilevel in range(0,nlevelmax):
        hydro.append(Hydro(ndim,nvar))
        hydro[ilevel].level = ilevel
        
    numbl = np.zeros([nlevelmax,ncpu],dtype=np.int32)
    
    # Reading and computing total AMR grids count
    for icpu in cpulist:

        car2 = str(icpu).zfill(5)
        if(backup):
            filename = path+"/backup_"+car1+"/"+prefix+"."+car2
        else:
            filename = path+"/output_"+car1+"/"+prefix+"."+car2

        skip = 16
        for ilevel in range(levelmin-1,nlevelmax):
            offset = skip+4*(ilevel+1-levelmin)
            numbl[ilevel,icpu-1] = np.fromfile(filename,dtype=np.int32,count=1,offset=offset)[0]
            hydro[ilevel].ngrid = hydro[ilevel].ngrid + numbl[ilevel,icpu-1]

    # Allocating memory
    for ilevel in range(0,nlevelmax):
        hydro[ilevel].u = np.zeros([nvar,2**ndim,hydro[ilevel].ngrid],dtype=float)
        hydro[ilevel].nvar = nvar

    iskip = np.zeros(nlevelmax, dtype=int)
    nvartot = nvar*2**ndim
    
    # Reading and storing data
    for icpu in cpulist:

        car2 = str(icpu).zfill(5)
        if(backup):
            filename = path+"/backup_"+car1+"/"+prefix+"."+car2
        else:
            filename = path+"/output_"+car1+"/"+prefix+"."+car2

        offset = 16 + 4*(nlevelmax+1-levelmin)
        
        for ilevel in range(levelmin-1,nlevelmax):
            ncache = numbl[ilevel,icpu-1]

            if(backup):
                transfer = np.fromfile(filename,dtype=np.float64,count=nvartot*ncache,offset=offset)
            else:
                transfer = np.fromfile(filename,dtype=np.float32,count=nvartot*ncache,offset=offset)

            transfer = np.reshape(transfer,(ncache,nvar,2**ndim))            
            transfer = np.transpose(transfer,(1,2,0))

            # Store cell hydro variables
            for ivar in range(0,nvar):
                for ind in range(0,2**ndim):
                    hydro[ilevel].u[ivar,ind,iskip[ilevel]:iskip[ilevel]+ncache] = transfer[ivar,ind]

            if(backup):
                offset = offset + ncache*nvartot*8
            else:
                offset = offset + ncache*nvartot*4

            iskip[ilevel] = iskip[ilevel] + ncache

    return hydro

def mk_image(x,y,dx,var):
    """
    Function to make image from cell data
    """
    xmin = np.min(x-dx/2)
    xmax = np.max(x+dx/2)
    ymin = np.min(y-dx/2)
    ymax = np.max(y+dx/2)

    dxmin = np.min(dx)
    dxmax = np.max(dx)

    nx = int((xmax-xmin)/dxmax)*int(dxmax/dxmin)
    ny = int((ymax-ymin)/dxmax)*int(dxmax/dxmin)

    nlev = int(np.log(dxmax/dxmin)/np.log(2))+1

    print("Making image of size: ",nx,ny)
    
    image = np.zeros((nx,ny))
    
    for lev in range(0,nlev):

        dxloc = dxmax/2**lev

        # Filter cells on the level
        filt = dx == dxloc

        # Skip levels without cells
        if (filt.sum() < 1):
            continue

        up_samp = int(2**(nlev-lev-1))

        # Setup the bins
        nxloc = int(nx/up_samp)
        nyloc = int(ny/up_samp)

        bins = (nxloc,nyloc)

        # Create the image
        H, _, _ = np.histogram2d(x[filt],
                                 y[filt],
                                 bins=bins,
                                 range=((xmin,xmax),(ymin,ymax)),weights=var[filt])

        if lev < nlev:
            H = H.repeat(up_samp, axis=1).repeat(up_samp, axis=0)

        image += H

    return image.T

def mk_cube(x,y,z,dx,var):
    """
    Function to make Cartesian cube from cell data
    """
    xmin = np.min(x-dx/2)
    xmax = np.max(x+dx/2)
    ymin = np.min(y-dx/2)
    ymax = np.max(y+dx/2)
    zmin = np.min(z-dx/2)
    zmax = np.max(z+dx/2)

    dxmin = np.min(dx)
    dxmax = np.max(dx)

    nx = int((xmax-xmin)/dxmax)*int(dxmax/dxmin)
    ny = int((ymax-ymin)/dxmax)*int(dxmax/dxmin)
    nz = int((zmax-zmin)/dxmax)*int(dxmax/dxmin)

    nlev = int(np.log(dxmax/dxmin)/np.log(2))+1

    print("Making cube of size: ",nx,ny,nz)
    
    cube = np.zeros((nx,ny,nz))
    
    for lev in range(0,nlev):

        dxloc = dxmax/2**lev

        # Filter cells on the level
        filt = dx == dxloc

        # Skip levels without cells
        if (filt.sum() < 1):
            continue

        up_samp = int(2**(nlev-lev-1))

        # Setup the bins
        nxloc = int(nx/up_samp)
        nyloc = int(ny/up_samp)
        nzloc = int(nz/up_samp)

        bins = (nxloc,nyloc,nzloc)

        points = np.column_stack((x[filt],y[filt],z[filt]))

        # Create the image
        C, _ = np.histogramdd(points,
                              bins=bins,
                              range=((xmin,xmax),(ymin,ymax),(zmin,zmax)),
                              weights=var[filt])

        if lev < nlev:
            C = C.repeat(up_samp, axis=2).repeat(up_samp, axis=1).repeat(up_samp, axis=0)

        cube += C

    return cube.T

def rotate_view(c,**kwargs):
    """This function rotate the input cells into a view where the z-axis is aligned with the
    angular momentum vector and the x- and y-axis are in the rotation plane.

    Args:
        c: an object of type cell.

    Optional args:

        center: a numpy array containing the coordinates of the center

        velocity: a numpy array containing the velocity of the center

    Returns:
        x, y, z: the 3-coordinates of the input cells after the rotation.

    Example:
        import miniramses as ram
        c = ram.rd_cell(12,center=[0.5,0.5,0.5],radius=0.1)
        x, y, z = ram.rotate_view(c,center=[0.5,0.5,0.5],velocity=[0,0,0])

    Authors: Carlos Sarkis (Princeton University, July 2025)
    """

    center = kwargs.get("center")
    velocity = kwargs.get("velocity")

    if(center is None):
        xc=np.mean(c.x[0])
        yc=np.mean(c.x[1])
        zc=np.mean(c.x[2])
    else:
        xc=center[0]
        yc=center[1]
        zc=center[2]

    if(velocity is None):
        uc=np.mean(c.u[1])
        vc=np.mean(c.u[2])
        wc=np.mean(c.u[3])
    else:
        uc=velocity[0]
        vc=velocity[1]
        wc=velocity[2]

    x0=c.x[0]-xc
    y0=c.x[1]-yc
    z0=c.x[2]-zc
    u0=c.u[1]-uc
    v0=c.u[2]-vc
    w0=c.u[3]-wc

    m = c.u[0] * c.dx**3  # mass of each cell
    Lx = np.sum(m * (y0 * w0 - z0 * v0))
    Ly = np.sum(m * (z0 * u0 - x0 * w0))
    Lz = np.sum(m * (x0 * v0 - y0 * u0))
    L = np.array([Lx, Ly, Lz])
    L_hat = L / np.linalg.norm(L)  # normalized direction vector

    # New basis: u3 = disk normal, u1 and u2 = in-plane axes
    u3 = L_hat
    u1 = np.cross(u3, [0, 0, 1])
    if np.linalg.norm(u1) == 0:
        u1 = np.cross(u3, [0, 1, 0])
    u1 /= np.linalg.norm(u1)
    u2 = np.cross(u3, u1)

    # Rotation matrix: columns = new basis vectors
    R = np.vstack([u1, u2, u3]).T  # shape (3, 3)

    # Apply rotation to all positions
    coords = np.vstack([x0, y0, z0])        # shape: (3, N)
    rotated = R.T @ coords                  # rotate into new frame

    return rotated[0], rotated[1], rotated[2]

class Cell:
    def __init__(self,nndim,nnvar):
        self.ncell = 0
        self.ndim = nndim
        self.nvar = nnvar
        self.x = np.empty(shape=(nndim,0))
        self.u = np.empty(shape=(nnvar,0))
        self.dx = np.empty(shape=(0))
        self.level = np.empty(shape=(0),dtype=np.int8)

def rd_cell(nout,**kwargs):
    """This function reads RAMSES AMR and hydro files (unformatted Fortran binary) 
    as produced by the RAMSES code in the snapshot directory output_00* 
    and store it in a variable containing all the hydro leaf cells information (Cell object).

    Args:
        nout: the RAMSES snapshot number. For example output_000012 corresponds to nout=12.

    Optional args:

        center: a numpy array containing the coordinates of the center of the sphere restricting the region to read in data.

        radius: the radius of the sphere restricting the region to read in data.

    Returns:
        A variable c (class Cell) object defined as:
            c.ncell: number of AMR cells
            c.ndim: number of space dimensions
            c.nvar: number of hydro variables
            c.x: coordinates of the cells. c.x[0] gives the x coordinate as a numpy array.
            c.u: hydro variables in each cell. For example, c.u[0] gives the gas density as a numpy array.
            c.dx: array containing the individual AMR cell sizes.
            c.level: refinement levels of cells.

    Example:
        import miniramses as ram
        c = ram.rd_cell(12,center=[0.5,0.5,0.5],radius=0.1)
        print(np.max(c.dx))

    Authors: Romain Teyssier (Princeton University, October 2022)
    """
    
    path = kwargs.get("path","./")
    center = kwargs.get("center")
    radius = kwargs.get("radius")
    geom = kwargs.get("geom","circle")

    a = rd_amr(nout,**kwargs)
    h = rd_hydro(nout,**kwargs)

    nlevelmax = len(a)
    ndim = a[0].ndim
    nvar = h[0].nvar
    boxlen = a[0].boxlen
    
    offset = np.zeros([ndim,2**ndim])
    if (ndim == 1):
        offset[0,:]=[-0.5,0.5]
    if (ndim == 2):
        offset[0,:]=[-0.5,0.5,-0.5,0.5]
        offset[1,:]=[-0.5,-0.5,0.5,0.5]
    if (ndim == 3):
        offset[0,:]=[-0.5,0.5,-0.5,0.5,-0.5,0.5,-0.5,0.5]
        offset[1,:]=[-0.5,-0.5,0.5,0.5,-0.5,-0.5,0.5,0.5]
        offset[2,:]=[-0.5,-0.5,-0.5,-0.5,0.5,0.5,0.5,0.5]

    ncell = 0
    for ilev in range(0,nlevelmax):
        ncell = ncell + np.count_nonzero(a[ilev].refined == False)

    print("Found",ncell,"leaf cells")
    print("Extracting leaf cells...")

    c = Cell(ndim,nvar)
    c.ncell = ncell

    c.x = np.zeros([ndim, ncell])
    c.u = np.zeros([nvar, ncell])
    c.dx = np.zeros(ncell)
    c.level = np.zeros(ncell, dtype=np.int8)

    icell = 0
    for ilev in range(0,nlevelmax):
        dx = 0.5*boxlen/2**ilev
        for ind in range(0,2**ndim):
            mask = a[ilev].refined[ind] == False
            nc = np.count_nonzero(mask)
            if (nc > 0):
                for idim in range(0,ndim):
                    c.x[idim, icell:icell+nc] = (2*a[ilev].xg[idim,mask]+1+offset[idim,ind])*dx
                for ivar in range(0,nvar):
                    c.u[ivar, icell:icell+nc] = h[ilev].u[ivar,ind,mask]
                c.dx[icell:icell+nc] = dx
                c.level[icell:icell+nc] = ilev
                icell = icell + nc

    # Filtering cells
    if ( not (center is None)  and not (radius is None) ):

        # Periodic boundaries
        for idim in range(0,ndim):
            xx = c.x[idim]-center[idim]
            xx[xx>boxlen/2]=xx[xx>boxlen/2]-boxlen
            xx[xx<-boxlen/2]=xx[xx<-boxlen/2]+boxlen
            c.x[idim] = xx+center[idim]
        if geom == "circle":
            if ndim==1:
                r = np.sqrt((c.x[0]-center[0])**2) - dx
            elif ndim==2:
                r = np.sqrt((c.x[0]-center[0])**2+(c.x[1]-center[1])**2) - dx
            elif ndim==3:
                r = np.sqrt((c.x[0]-center[0])**2+(c.x[1]-center[1])**2+(c.x[2]-center[2])**2) - dx
        elif geom == "square":
            if ndim==1:
                r = np.abs(c.x[0]-center[0]) - dx
            elif ndim==2:
                r = np.maximum.reduce([np.abs(c.x[0]-center[0]), np.abs(c.x[1]-center[1])]) - dx
            elif ndim==3:
                r = np.maximum.reduce([np.abs(c.x[0]-center[0]), np.abs(c.x[1]-center[1]), np.abs(c.x[2]-center[2])]) - dx
        c.ncell = np.count_nonzero(r < radius)
        c.u  = c.u[:,r < radius]
        c.x  = c.x[:,r < radius]
        c.dx = c.dx[r < radius]
        c.level = c.level[r < radius]

    if(ndim==1):
        c.x = c.x[0]
        ind = np.argsort(c.x)
        c.x = c.x[ind]
        c.dx = c.dx[ind]
        c.level = c.level[ind]
        for ivar in range(0,nvar):
            c.u[ivar]=c.u[ivar,ind]
        
    return c

def save_cell(c,filename):

    with open(filename,'wb') as f:
        np.save(f,c.ncell)
        np.save(f,c.ndim)
        np.save(f,c.nvar)
        np.save(f,c.dx)
        np.save(f,c.x)
        np.save(f,c.u)
        np.save(f,c.level)

def load_cell(filename):

    with open(filename,'rb') as f:
        ncell = np.load(f)
        ndim = np.load(f)
        nvar = np.load(f)
        c = Cell(ndim,nvar)
        c.ncell = ncell
        c.ndim = ndim
        c.nvar = nvar
        c.dx = np.append(c.dx,np.load(f))
        c.x  = np.append(c.x, np.load(f),axis=1)
        c.u  = np.append(c.u, np.load(f),axis=1)
        c.level  = np.append(c.level, np.load(f),axis=1)

    return c

class Snap1d:
    pass

def rd_log(filename,**kwargs):
    """This function reads the standard ouput (aka log file)
    as produced by the RAMSES code for 1D simulations.

    Args:
        filename: the log file name (usually run.log).

    Optional args:

        None

    Returns:
        A variable r (class Run) object defined as:
            c.ncell: number of AMR cells.
            c.lev: level of refinement of the cells.
            c.x: coordinates of the cells.
            c.d: density of the cells.
            c.u: velocity field x-component.
            c.v: velocity field y-component.
            c.w: velocity field z-component.
            c.d: pressure of the cells.
            c.A: magnetic field x-component.
            c.B: magnetic field y-component.
            c.C: magnetic field z-component.

    Example:
        import miniramses as ram
        r = ram.rd_log("run.log")
        plt.plot(r["x"],r["d"]))
    
    Authors: Romain Teyssier (Princeton University, October 2022)
    """
    cmd="grep -n Output "+filename+" > /tmp/out.txt"
    os.system(cmd)
    lines = ascii.read("/tmp/out.txt")
    print("Found "+str(len(lines))+" output(s)")

    r=[]
    for out in range(0,len(lines)):

        i = int(lines["col1"][out][:-1])
        n = int(lines["col3"][out])

        data = ascii.read(filename,header_start=i-3,data_start=i-2,data_end=i+n-2)

        r.append(Snap1d())
        r[out].ncell=n
        r[out].lev=np.array(data["lev"],dtype='int')
        r[out].x=np.array(data["x"])
        r[out].d=np.array(data["d"])
        r[out].p=np.array(data["P"])
        r[out].u=np.array(data["u"])
        if len(data.columns)>5:
            r[out].v=np.array(data["v"])
            r[out].w=np.array(data["w"])
            r[out].A=np.array(data["A"])
            r[out].B=np.array(data["B"])
            r[out].C=np.array(data["C"])

    return r

class Info:
    def __init__(self,nncpu):
        self.bound_key = np.zeros(shape=(nncpu+1),dtype=np.double)
        
def rd_info(nout,**kwargs):
    
    backup = kwargs.get("backup",False)
    path = kwargs.get("path","./")
    units = kwargs.get("units",False)

    car1 = str(nout).zfill(5)
    if(backup):
        filename = path+"/backup_"+car1+"/info.txt"
    else:
        filename = path+"/output_"+car1+"/info.txt"

    info=ascii.read(filename,delimiter="=",format='no_header')

    nfile=int(info[0][1])
    ncpu =int(info[1][1])
    if( not backup ):
        ncpu=nfile

    i = Info(ncpu)

    i.path=path
    i.nout=nout
    i.ncpu=ncpu
    i.ndim=int(info[2][1])
    i.levelmin=int(info[3][1])
    i.nlevelmax=int(info[4][1])
    #
    i.boxlen=info[7][1]
    i.time=info[8][1]
    i.texp=info[9][1]
    i.aexp=info[10][1]
    i.redshift= 0.0 if (1.0/i.aexp)-1.0 < 0 else (1.0/i.aexp)-1.0 
    #
    i.H0=info[11][1]
    i.h= i.H0 / 100.0 
    #
    i.omega_m=info[12][1]
    i.omega_l=info[13][1]
    i.omega_k=info[14][1]
    i.omega_b=info[15][1]
    i.gamma=info[16][1]
    #
    i.unit_l=info[17][1]
    i.unit_d=info[18][1]
    i.unit_t=info[19][1]
    i.unit_v = i.unit_l / i.unit_t
    i.unit_m = i.unit_d * i.unit_l**3

    i.G = 3 / (8 * np.pi) * i.omega_m * i.aexp # gravitational constant in code units
    i.rho_m = 1.0

    if units:
        import astropy.units as u
        import astropy.constants as const

        # convert to physical units
        i.unit_l *= u.cm
        i.unit_d *= u.g / u.cm**3
        i.unit_t *= u.s
        i.unit_v *= u.cm / u.s
        i.unit_m *= u.g 
        
        i.H0     *= u.km / u.s / u.Mpc
        i.boxlen *= i.unit_l
        i.time   *= i.unit_t 
        i.texp   *= i.unit_t

        i.c = const.c.to(u.cm / u.s).value / i.unit_v.value # speed of light in code units


    # Get the temperature conversion
    i.unit_T2 = ((i.unit_l / i.unit_t)**2) * 1.6605390e-24 / 1.3806490e-16

    # Read all particle-type header files (*_header.txt) in the output
    # directory and store the total number of particles of each type as an
    # info attribute named "n"+prefix (e.g. part_header.txt -> i.npart,
    # tree_header.txt -> i.ntree, peak_part_header.txt -> i.npeak_part).
    # Non particle-type headers (e.g. grav_header.txt) are skipped.
    outdir = os.path.dirname(filename)
    i.nheader = {}
    for header_file in sorted(glob.glob(os.path.join(outdir,"*_header.txt"))):
        prefix = os.path.basename(header_file)[:-len("_header.txt")]
        with open(header_file) as f:
            lines = [line.strip() for line in f]
        # only particle-type headers report a total particle count
        if "Total number of particles" not in lines:
            continue
        idx = lines.index("Total number of particles")
        nparticles = int(lines[idx+1])
        setattr(i,"n"+prefix,nparticles)
        i.nheader[prefix] = nparticles

    rd_rt_info = kwargs.get("rt",False)
    if rd_rt_info:
      if(backup):
          rt_filename = path+"/backup_"+car1+"/rt_info.txt"
      else:
          rt_filename = path+"/output_"+car1+"/rt_info.txt"

      rt_info=ascii.read(rt_filename,delimiter="=",format='no_header')

      i.nrtvar = int(rt_info[0][1])
      i.nrtgrp = int(rt_info[1][1])
      i.nion   = int(rt_info[2][1])
      i.iion   = int(rt_info[3][1])
      i.x_h    = rt_info[4][1]
      i.y_he   = rt_info[5][1]
      i.unit_np= rt_info[6][1]
      i.unit_fp= rt_info[7][1]
      i.rt_c_fraction= np.array(rt_info[8][1].split(),dtype=float)
      i.groupL0 = np.array(rt_info[10][1].split(),dtype=float)
      i.groupL1 = np.array(rt_info[11][1].split(),dtype=float)
      i.group_egy = np.zeros(i.nrtgrp)
      i.group_csn = np.zeros((i.nrtgrp, i.nion))
      i.group_cse = np.zeros((i.nrtgrp, i.nion))

      for igrp in range(i.nrtgrp):
        iline = 14 + igrp*4
        print(iline,igrp)
        i.group_egy[igrp] = rt_info[iline][1]
        i.group_csn[igrp,:] = np.array(rt_info[iline+1][1].split(),dtype=float)
        i.group_cse[igrp,:] = np.array(rt_info[iline+2][1].split(),dtype=float)

    return i

def hilbert3d(x,y,z,bit_length):
    
    state_diagram = [ 1, 2, 3, 2, 4, 5, 3, 5,
                      0, 1, 3, 2, 7, 6, 4, 5,
                      2, 6, 0, 7, 8, 8, 0, 7,
                      0, 7, 1, 6, 3, 4, 2, 5,
                      0, 9,10, 9, 1, 1,11,11,
                      0, 3, 7, 4, 1, 2, 6, 5,
                      6, 0, 6,11, 9, 0, 9, 8,
                      2, 3, 1, 0, 5, 4, 6, 7,
                      11,11, 0, 7, 5, 9, 0, 7,
                      4, 3, 5, 2, 7, 0, 6, 1,
                      4, 4, 8, 8, 0, 6,10, 6,
                      6, 5, 1, 2, 7, 4, 0, 3,
                      5, 7, 5, 3, 1, 1,11,11,
                      4, 7, 3, 0, 5, 6, 2, 1,
                      6, 1, 6,10, 9, 4, 9,10,
                      6, 7, 5, 4, 1, 0, 2, 3,
                      10, 3, 1, 1,10, 3, 5, 9,
                      2, 5, 3, 4, 1, 6, 0, 7,
                      4, 4, 8, 8, 2, 7, 2, 3,
                      2, 1, 5, 6, 3, 0, 4, 7,
                      7, 2,11, 2, 7, 5, 8, 5,
                      4, 5, 7, 6, 3, 2, 0, 1,
                      10, 3, 2, 6,10, 3, 4, 4,
                      6, 1, 7, 0, 5, 2, 4, 3]

    state_diagram = np.array(state_diagram)
    state_diagram = state_diagram.reshape((8,2,12),order='F')

    n = len(x)
    order = np.zeros(n,dtype="double")
    x_bit_mask = np.zeros(bit_length  ,dtype="bool")
    y_bit_mask = np.zeros(bit_length  ,dtype="bool")
    z_bit_mask = np.zeros(bit_length  ,dtype="bool")
    i_bit_mask = np.zeros(3*bit_length,dtype=bool)
    
    for ip in  range(0,n):
        
        for i in range(0,bit_length):
            x_bit_mask[i] = x[ip] & (1 << i)
            y_bit_mask[i] = y[ip] & (1 << i)
            z_bit_mask[i] = z[ip] & (1 << i)
            
        for i in range(0,bit_length):
            i_bit_mask[3*i+2] = x_bit_mask[i]
            i_bit_mask[3*i+1] = y_bit_mask[i]
            i_bit_mask[3*i  ] = z_bit_mask[i]
            
        cstate = 0
        for i in range(bit_length-1,-1,-1):
            b2 = 0
            if (i_bit_mask[3*i+2]):
                b2 = 1
            b1 = 0
            if (i_bit_mask[3*i+1]):
                b1 = 1
            b0 = 0
            if (i_bit_mask[3*i  ]):
                b0 = 1
            sdigit = b2*4 + b1*2 + b0
            nstate = state_diagram[sdigit,0,cstate]
            hdigit = state_diagram[sdigit,1,cstate]
            i_bit_mask[3*i+2] = hdigit & (1 << 2)
            i_bit_mask[3*i+1] = hdigit & (1 << 1)
            i_bit_mask[3*i  ] = hdigit & (1 << 0)
            cstate = nstate
            
        order[ip]= 0
        for i in range(0,3*bit_length):
            b0 = 0
            if (i_bit_mask[i]):
                b0 = 1
            order[ip] = order[ip] + float(b0)*2.**i
                
    return order

def hilbert2d(x,y,bit_length):
    
    state_diagram = [ 1, 0, 2, 0, 
                      0, 1, 3, 2, 
                      0, 3, 1, 1, 
                      0, 3, 1, 2, 
                      2, 2, 0, 3, 
                      2, 1, 3, 0, 
                      3, 1, 3, 2, 
                      2, 3, 1, 0 ]
    
    state_diagram = np.array(state_diagram)    
    state_diagram = state_diagram.reshape((4,2,4), order='F')
    
    n = len(x)
    order = np.zeros(n,dtype="double")
    x_bit_mask = np.zeros(bit_length  ,dtype="bool")
    y_bit_mask = np.zeros(bit_length  ,dtype="bool")
    i_bit_mask = np.zeros(2*bit_length,dtype=bool)
    
    for ip in  range(0,n):
        
        for i in range(0,bit_length):
            x_bit_mask[i] = bool(x[ip] & (1 << i))
            y_bit_mask[i] = bool(y[ip] & (1 << i))
            
        for i in range(0,bit_length):
            i_bit_mask[2*i+1] = x_bit_mask[i]
            i_bit_mask[2*i  ] = y_bit_mask[i]
            
        cstate = 0
        for i in range(bit_length-1,-1,-1):
            b1 = 0
            if (i_bit_mask[2*i+1]):
                b1 = 1
            b0 = 0
            if (i_bit_mask[2*i  ]):
                b0 = 1
            sdigit = b1*2 + b0
            nstate = state_diagram[sdigit,0,cstate]
            hdigit = state_diagram[sdigit,1,cstate]
            i_bit_mask[2*i+1] = hdigit & (1 << 1)
            i_bit_mask[2*i  ] = hdigit & (1 << 0)
            cstate = nstate
            
        order[ip]= 0
        for i in range(0,2*bit_length):
            b0 = 0
            if (i_bit_mask[i]):
                b0 = 1
            order[ip] = order[ip] + float(b0)*2.**i
                
    return order

def get_cpu_list(info,**kwargs):

    center = kwargs.get("center")
    radius = kwargs.get("radius")
    center = np.array(center)
    radius = float(radius)
    
    for ilevel in range(0,info.nlevelmax):
        dx = 1/2**ilevel
        if (dx < 2*radius/info.boxlen):
            break

    levelmin = np.max([ilevel,1])
    bit_length = levelmin-1
    nmax = 2**bit_length
    ndim = info.ndim
    ncpu = info.ncpu
    nlevelmax = info.nlevelmax
    dkey = 2**(ndim*(nlevelmax+1-bit_length))
    ibound = [0, 0, 0, 0, 0, 0]
    if(bit_length > 0):
        ibound[0:3] = (center-radius)*nmax/info.boxlen
        ibound[3:6] = (center+radius)*nmax/info.boxlen
        ibound[0:3] = np.array(ibound[0:3]).astype(int)
        ibound[3:6] = np.array(ibound[3:6]).astype(int)
        ndom = 8
        idom = [ibound[0], ibound[3], ibound[0], ibound[3], ibound[0], ibound[3], ibound[0], ibound[3]]
        jdom = [ibound[1], ibound[1], ibound[4], ibound[4], ibound[1], ibound[1], ibound[4], ibound[4]]
        kdom = [ibound[2], ibound[2], ibound[2], ibound[2], ibound[5], ibound[5], ibound[5], ibound[5]]
        order_min = hilbert3d(idom,jdom,kdom,bit_length)
    else:
        ndom = 1
        order_min = np.array([0.])
        
    bounding_min = order_min*dkey
    bounding_max = (order_min+1)*dkey

    cpu_min = np.zeros(ndom, dtype=int)
    cpu_max = np.zeros(ndom, dtype=int)
    for icpu in range(0,ncpu):
        for idom in range(0,ndom):
            if( (info.bound_key[icpu] <= bounding_min[idom]) and (info.bound_key[icpu+1] > bounding_min[idom]) ):
                cpu_min[idom] = icpu+1
            if( (info.bound_key[icpu] < bounding_max[idom]) and (info.bound_key[icpu+1] >= bounding_max[idom]) ):
                cpu_max[idom] = icpu+1


    ncpu_read = 0
    cpu_read = np.zeros(ncpu, dtype=bool)
    cpu_list = []
    for idom in range(0,ndom):
        for icpu in range(cpu_min[idom]-1,cpu_max[idom]):
            if ( not cpu_read[icpu] ):
                cpu_list.append(icpu+1)
                ncpu_read = ncpu_read+1
                cpu_read[icpu] = True

    return cpu_list

def visu(x,y,dx,v,**kwargs):
    '''The simple visualization function visu() make a 2D scatter plot from RAMSES AMR data. 

    Args:

        x: the x-coordinate of the cells to show on the scatter plot.
        y: the y-coordinate of the cells to show on the scatter plot.
        dx: the size of the cells to show on the scatter plot.
        v: the value to show as a color square contained in the cell.

    Optional args:

        vmin: minimum value for the input array v to use in the color range
        vmax: maximum value for the input array v to use in the color range 
        log: when set, use the log of the input array v in the color range
        colorbar: when True, draw a colorbar (default: True)
        log_floor: lower bound applied to |v| before log10 (default 0)
        sort: useful only for 3D data. Plot the square symbola in the scatter plot in increasing order of array sort.

    Returns:

        Output a scatter plot figure of size 1000 pixels aside.

    Example:

        Example for a 2D or 3D RAMSES dataset using variable c from the object Cell. 
        import miniramses as ram
        c=ram.rd_cell(2)
        ram.visu(c.x[0],c.x[1],c.dx,c.u[0],sort=c.u[0],log=1,vmin=-3,vmax=1)

    Authors: Romain Teyssier (Princeton University, October 2022)
    '''

    xmin=np.min(x-dx/2)
    xmax=np.max(x+dx/2)
    ymin=np.min(y-dx/2)
    ymax=np.max(y+dx/2)
    
    log = kwargs.get("log",None)
    vmin = kwargs.get("vmin",None)
    vmax = kwargs.get("vmax",None)
    sort = kwargs.get("sort",None)
    cmap = kwargs.get("cmap",'viridis')
    grid = kwargs.get("grid",None)
    log_floor = kwargs.get("log_floor",0)
    show_colorbar = kwargs.get("colorbar",True)
    
    if( not (log is None)):
        # Standard log scaling: log data; transform limits consistently
        v = np.log10(np.maximum(np.abs(v), float(log_floor)))
        if not (vmin is None):
            vmin = np.log10(float(vmin))
        if not (vmax is None):
            vmax = np.log10(float(vmax))

    print("min=",np.min(v)," max=",np.max(v))

    if( not (sort is None)):
        ind = np.argsort(sort)
    else:
        ind = np.arange(0,v.size)

    olddpi = plt.rcParams['figure.dpi']
    plt.rcParams['figure.dpi'] = 58
    px = 1/plt.rcParams['figure.dpi'] 
    fig, ax = plt.subplots(figsize=(1000*px,1000*px))
    ax.set_xlim([xmin,xmax])
    ax.set_ylim([ymin,ymax])
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
    plt.scatter(x,y,s=0.0001)
    rescale=np.maximum(xmax-xmin,ymax-ymin)        
    ax.set_aspect("equal")
    edgec = None
    linew = None
    if( not (grid is None)):
        edgec='black'
        linew=0.5
    plt.scatter(x[ind],y[ind],c=v[ind],s=(dx[ind]*800/rescale)**2,marker="s",vmin=vmin,vmax=vmax,
                cmap=cmap,edgecolor=edgec,linewidth=linew)
    if show_colorbar:
        plt.colorbar(shrink=0.8)
    plt.rcParams['figure.dpi'] = olddpi

def mk_movie(**kwargs):
    '''The function mk_movie() takes 2D data files containing maps and converts them into a sequence of images, 
    before combining them into a movie. It requires a standard set of python packages and the Linux packages
    ffmpeg and convert (ImageMagick).
    
    Args:
    
        start: starting index of the sequence of numpy array you wish to turn into image frames.

        stop: number of arrays you wish to be turned into plots. 
            This will be the variable "snum" for the end product. 
            For now, if you wish to test out the function, 
            you can try out other smaller values to adjust the image for your preferences.

        path: path leading to the directory where your files are stored, Default: "."

        prefix: starting name of a typical file. Ex: if you have 50 files, called "fig01.npy", "fig02.npy" … "fig50.npy", write in "fig".

        fill: This is for the zfill parameter. If your files are standardized into "fig001.npy", "fig002.npy"… "fig100.npy",
            write in 3, for example. If this is not how your files are formatted, write in the number 1.

        suffix: suffix at the end of a file: Ex: ".npy", ".map", etc…

        cmap: write in what color you wish your array to be displayed in (value for cmap). Options include "Reds", "Blues", and more.

        cbar: write "YES" for this parameter if you want your figure to have a colorbar. Write anything else if not.

        cbunit: units of the colormapping to be displayed next to the colorbar: Ex: "Concentration [code units]"
              If you do not plan on using a colorbar, write in any script.

        tunit: units of time displayed by rd_img. Ex: "seconds", "minutes", "hours", "[code units]"

        bunit: units of the box size displayed by rd_img. Ex: "cm", "kpc", "Mpc", "[code units]" …

        fname: starting name of each of your images.

        mvname: what you want your movie to be called.

    Returns:

        info: a string stating that the movie was done.

    Exemple:

        import miniramses as ram
        info = ram.mk_movie(start=100,stop=2000,path="../movie1",prefix="dens_",fill=5,suffix=".map",cmap="Reds", 
                cbar="YES", cbunit="log Density [H/cc]", tunit="Gyr",
                fname="img", mvname="movie", vmin=-1, vmax=6)

    By default, the movie's framerate is 30 frames per second, at a resolution of 420p
    You can edit this function and its parameters according to what fits your model best.

    As it runs, the function will print the files it is currently converting.

    Authors: Thomas Decugis and Romain Teyssier (Princeton University, October 2022)
    '''
    start = kwargs.get("start",1)
    stop = kwargs.get("stop",1)
    prefix = kwargs.get("prefix")
    suffix = kwargs.get("suffix")
    fill = kwargs.get("fill",5)
    path = kwargs.get("path",".")
    vmin = kwargs.get("vmin",None)
    vmax = kwargs.get("vmax",None)
    cmap = kwargs.get("cmap","Reds")
    cbar = kwargs.get("cbar",None)
    cbunit = kwargs.get("cbunit",None)
    tunit = kwargs.get("tunit"," ")
    bsize = kwargs.get("bsize",1)
    bunit = kwargs.get("bunit","[code units]")
    fname = kwargs.get("fname","frame")
    mvname = kwargs.get("mvname","movie")
    
    cmd="curl https://tigress-web.princeton.edu/~rt3504/DAT/logo_essai.jpg --output logo_essai.jpg"
    os.system(cmd)
    concom = "convert logo_essai.jpg -resize 280x200 logo_essai.png"
    os.system(concom)

    for snapshot in range(start, stop + 1): 
        ar = path + "/" + str(prefix) + str(snapshot).zfill(fill) + str(suffix)
        print(ar) #prints file that function is working on.

        map =rd_map(ar)
        time = map.time
        array = map.data
        
        if (not (cbar is None)):
            px = 1/plt.rcParams['figure.dpi']
            fig, ax = plt.subplots(figsize=(1000*px,1000*px))
            plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
            print(np.min(array),np.max(array))
            shw = ax.imshow(array, cmap = cmap, vmin=vmin, vmax=vmax, origin="lower", extent=[0,bsize,0,bsize])
            bar = plt.colorbar(shw,shrink=0.8)
            bar.set_label(cbunit, fontsize=18) 
            bar.ax.tick_params(labelsize=18) 
            plt.ylabel(bunit,fontsize=18)

        else:
            plt.imshow(array, cmap = cmap)#if you wish to graph model in a specific way, modify this program

        ax = plt.gca()
        txt = f't = {time:4.2f}' + tunit
        label = ax.set_xlabel(txt, fontsize = 18, color = "black")
        ax.xaxis.set_label_coords(0.1, 0.95)
        ax.tick_params(axis='both', labelsize=18)
        newname = str(fname)+ str(snapshot).zfill(fill) + ".png"
        print(newname)
        plt.savefig(newname) #saves created images as pngs under the name that was given
        if snapshot == start:
            plt.show()
        plt.close(fig)
        com = "convert logo_essai.png -bordercolor white -border 0.1 " + newname + " +swap -geometry +100+850 -composite " + newname
        os.system(com)
    print("Input files converted into frames: done")
    moviecom = "ffmpeg -y -r 30 -f image2 -s 1000x1000 -start_number " +str(start)+" -i " + str(fname) + "%05d.png" + " -vcodec libx264 -crf 25  -pix_fmt yuv420p " + str(mvname) + ".mp4" 
    os.system(moviecom)
    ok = "Movie: done"
    print(ok)
    return ok 

class ClumpCat:
    """
    This is the class for RAMSES clump catalog.
    """
    def __init__(self):
        """
        This function initialize the clump catalog.
        """
        self.index = np.empty(shape=(0),dtype=int)
        self.parent= np.empty(shape=(0),dtype=int)
        self.halo = np.empty(shape=(0),dtype=int)
        self.ncell = np.empty(shape=(0),dtype=int)
        self.npart = np.empty(shape=(0),dtype=int)
        self.x = np.empty(shape=(0))
        self.y = np.empty(shape=(0))
        self.z = np.empty(shape=(0))
        self.u = np.empty(shape=(0))
        self.v = np.empty(shape=(0))
        self.w = np.empty(shape=(0))
        self.mpatch = np.empty(shape=(0))
        self.mass = np.empty(shape=(0))
        self.dmax = np.empty(shape=(0))
        self.dmin = np.empty(shape=(0))
        self.dsad = np.empty(shape=(0))
        self.dave = np.empty(shape=(0))
        self.r200c = np.empty(shape=(0))
        self.rmax = np.empty(shape=(0))
        self.c200c = np.empty(shape=(0))

def rd_clump(nout,**kwargs):
    """
    This function reads and compiles data for position, mass,
    density, index, etc from the clump catalog.

    Args:
        nout: output file number

    Returns:
        A RAMSES clump catalog with all clump properties

    Authors: Josiah Taylor (Princeton University)
    """
    backup = kwargs.get("backup",False)
    center = kwargs.get("center")
    radius = kwargs.get("radius")
    path = kwargs.get("path","./")
    silent = kwargs.get("silent",False)

    car1 = str(nout).zfill(5)
    info = rd_info(nout,path=path,backup=backup)
    ncpu = info.ncpu
    ndim = info.ndim
    boxlen = info.boxlen

    output = str(nout).zfill(5)
    cat = ClumpCat()
    for i in range(0, ncpu):
        name = str(i+1).zfill(5)
        file_name = path+"/output_%s/clump.%s" % (output,name)
        read_cat = ascii.read(file_name)
        index = read_cat['index']
        parent = read_cat['parent']
        halo = read_cat['halo']
        ncell = read_cat['ncell']
        npart = read_cat['npart']
        x = read_cat['pos_x']
        y = read_cat['pos_y']
        z = read_cat['pos_z']
        u = read_cat['vel_x']
        v = read_cat['vel_y']
        w = read_cat['vel_z']
        dmin = read_cat['rho_min']
        dmax = read_cat['rho_max']
        dsad = read_cat['rho_sad']
        dave = read_cat['rho_ave']
        mpatch = read_cat['mpatch']
        mass = read_cat['mass']
        r200c = read_cat['r200'] #my current version uses 200c as the definition
        rmax = read_cat['rmax']
        c200c = read_cat['c200'] #my current version uses 200c as the definition

        cat.index = np.append(cat.index,index)
        cat.parent = np.append(cat.parent,parent)
        cat.halo = np.append(cat.halo,halo)
        cat.ncell = np.append(cat.ncell,ncell)
        cat.npart = np.append(cat.npart,npart)
        cat.x = np.append(cat.x,x)
        cat.y = np.append(cat.y,y)
        cat.z = np.append(cat.z,z)
        cat.u = np.append(cat.u,u)
        cat.v = np.append(cat.v,v)
        cat.w = np.append(cat.w,w)
        cat.mpatch = np.append(cat.mpatch,mpatch)
        cat.mass = np.append(cat.mass,mass)
        cat.dmax = np.append(cat.dmax,dmax)
        cat.dmin = np.append(cat.dmin,dmin)
        cat.dsad = np.append(cat.dsad,dsad)
        cat.dave = np.append(cat.dave,dave)
        cat.r200c = np.append(cat.r200c,r200c)
        cat.rmax = np.append(cat.rmax,rmax)
        cat.c200c = np.append(cat.c200c,c200c)

    # - - - - - - - -
    # derived fields
    # - - - - - - - -

    # m200c
    omega_mz = info.omega_m * info.aexp**(-3) / (info.omega_m * info.aexp**(-3) + info.omega_l)
    d200c = 200 / omega_mz  # 200 times the critical density

    cat.m200c = (4 * np.pi / 3) * cat.r200c**3 * d200c # (x rho_m = 1.0 in code units)
    cat.v200c = np.sqrt(info.G * cat.m200c / cat.r200c) # in code units

    # rs
    cat.rs = cat.r200c / cat.c200c

    # quantities at rmax
    xmax = cat.rmax / cat.rs

    # mass and circular velocity at rmax
    cat.mmax = (4 * np.pi / 3) * cat.rs**3 * (d200c * cat.c200c**3 / _f(cat.c200c)) * _f(xmax) # x rho_m = 1.0 in code units
    #cat.vmax = np.sqrt(info.G * cat.mmax / cat.rmax) # in code units
    c0 = 2.1626
    cat.vmax = cat.v200c * np.sqrt((cat.c200c/_f(cat.c200c))/(c0/_f(c0)))

    # tidal radius
    cat.rtidal = (cat.mpatch / (4 * np.pi * cat.dsad / 3))**(1/3)

    # most massive central clump in the halo 
    cat.is_central = cat.index == cat.parent 
    cat.is_most_dense = cat.index == cat.halo
    cat.is_satellite = cat.index != cat.parent

    # Filtering clumps
    if ( not (center is None)  and not (radius is None) ):

        # Periodic boundaries
        xx = cat.x-center[0]
        xx[xx>boxlen/2]=xx[xx>boxlen/2]-boxlen
        xx[xx<-boxlen/2]=xx[xx<-boxlen/2]+boxlen
        cat.x = xx+center[0]
        xx = cat.y-center[1]
        xx[xx>boxlen/2]=xx[xx>boxlen/2]-boxlen
        xx[xx<-boxlen/2]=xx[xx<-boxlen/2]+boxlen
        cat.y = xx+center[1]
        xx = cat.z-center[2]
        xx[xx>boxlen/2]=xx[xx>boxlen/2]-boxlen
        xx[xx<-boxlen/2]=xx[xx<-boxlen/2]+boxlen
        cat.z = xx+center[2]

        r = np.sqrt((cat.x-center[0])**2+(cat.y-center[1])**2+(cat.z-center[2])**2)

        for attr, value in cat.__dict__.items():
            if isinstance(value, np.ndarray) and len(value) == len(r):
                setattr(cat, attr, value[r < radius])

    if silent==False:
        txt = "Found "+str(len(cat.index))+" clumps"
        print(txt)
    return cat

def plot_tree(nout,pid,**kwargs):

    # read clump and sink files
    c=rd_clump(nout,**kwargs)
    s=rd_part(nout,prefix='tree',peak=True,**kwargs)

    # collect sinks in chosen clump
    ind=np.where(s.peak_id==pid)
    idp=s.birth_id[ind]
    idm=s.merging_id[ind]
    idt=s.tracking_id[ind]
    tp=s.birth_date[ind]
    tm=s.merging_date[ind]

    # sort sinks according to id
    isort=np.argsort(idp)
    idp=idp[isort]
    idm=idm[isort]
    tp=tp[isort]
    tm=tm[isort]

    plt.plot([idp[0],idp[0]],[tp[0],0],'r')
    print("merger tree particle id=",idp[0])
    for i in range(1,len(idp)):
        plt.plot([idp[i],idp[i]],[tp[i],tm[i]],'b')
        plt.plot([idp[i],idm[i]],[tm[i],tm[i]],'g')

class GraficFile:
    """
    Thid is the empty class for grafic files data
    """

def rd_grafic(filein):
    """This function reads a grafic file (unformatted Fortran binary)
    as produced by the MUSIC code.

    Args:
        filename: the complete path (including the name) of the input grafic file.

    Returns:
        A grafic (class GraficFile) object.

    Example:
        import miniramses as ram
        g = ram.rd_grafic("ic_deltab")
        plt.imshow(g.data[:,:,0],origin="lower")

    Authors: Romain Teyssier (Princeton University, October 2022)
    """
    with FortranFile(filein, 'r') as f:
        recl = ["i4", "i4", "i4", "f4", "f4", "f4", "f4", "f4", "f4", "f4", "f4"] 
        n1, n2, n3, dx, x1, x2, x3, a, omega_m, omega_l, h0 = f.read_record(*recl)
        n1=int(n1[0])
        n2=int(n2[0])
        n3=int(n3[0])
        print("Reading file "+filein)
        print("Found array of size=",n1,n2,n3)
        dat = np.zeros((n3,n2,n1))
        for k in range(n3):
            plane = f.read_reals('f4')
            dat[k,:,:] = plane.reshape((n2,n1))

    out = GraficFile()
    out.n1=n1
    out.n2=n2
    out.n3=n3
    out.dx=dx[0]
    out.x1=x1[0]
    out.x2=x2[0]
    out.x3=x3[0]
    out.omega_m=omega_m[0]
    out.omega_l=omega_l[0]
    out.h0=h0[0]
    out.data=np.array(dat.T)

    return out

def wr_grafic(dat,header1,header2,fileout):
    """This function writes a grafic file (unformatted Fortran binary)
    which is the file format produced e.g. by the MUSIC code.

    Args:
        dat: a 3D numpy array of type "f4"

        header1: a 1D numpy array with 3 elements of type "i4". It should contain the 3 dimensions of the input array.

        header2: a 1D mumpy array with 8 elements of type "f4". It should contain dx, xoff1, xoff2, xoff3 and 4 additional constants,

        filename: the complete path (including the name) of the output grafic file.

    Returns:
        Nothing

    Example:
        import miniramses as ram
        dat = np.zeros((512,512,512),dtype="f4")
        dx = 1./512.
        header1 = np.array([512,512,512],dtype="i4")
        header2 = np.array([dx,0,0,0,0,0,0,0],dtype="f4")
        ram.wr_grafic(dat,header1,header2,"ic_d")

    Authors: Romain Teyssier (Princeton University, October 2022)
    """
    with FortranFile(fileout, 'w') as f:
        f.write_record(header1,header2)
        n3 = int(header1[2])
        for k in range(n3):
            plane = dat[:, :, k]
            f.write_record(plane.T)



# HELPER FUNCTIONS 

def _f(x):
    """
    NFW profile helper function.

    Parameters
    ----------
    x : float or array-like
        Dimensionless radius (r/rs)

    Returns
    -------
    float or array-like
        Value of ln(1+x) - x/(1+x)
    """
    return np.log(1 + x) - x / (1 + x)
