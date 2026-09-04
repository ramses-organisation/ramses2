import numpy as np
import os
import re
import healpy as hp
import astropy.units as u
from astropy.cosmology import LambdaCDM

class LightconeReader:

    @staticmethod
    def rd_metadata(path, verbose=False):
        """
        Read the lightcone shell metadata from the .txt file
        
        Returns:
            Dictionary with keys: 'npart', 'aexp_old', 'aexp'
        """
        if verbose:
            print(f"Reading metadata from {path}")
        with open(path, 'r') as file:
            npart = int(file.readline().strip()) 

            # Read the scale factors
            aexp_old = float(file.readline().strip())
            aexp = float(file.readline().strip()) 
            aexp_center = (aexp + aexp_old)/2 

            # Compute redshifts
            z = 1/aexp - 1 
            z_old = 1/aexp_old - 1 
            z_center = 1/aexp_center - 1

        if verbose:
            print(f"Found {npart} particles")

        return {'npart': npart, 
                'aexp_old': aexp_old, 'aexp': aexp, 'aexp_center': aexp_center, 
                'z': z, 'z_old': z_old, 'z_center': z_center}

    @staticmethod
    def rd_part(path, nproperties=7, verbose=False):
        """
        Read the lightcone shell from the output directory.
        nproperties: number of non-idp properties per particle (default 7 for x,y,z,vx,vy,vz,mass)
        
        Returns:
            idp: numpy array of particle IDs (int32) with shape (npart,)
            properties: numpy array of properties (float32) with shape (nproperties, npart)
                       where rows are x, y, z, vx, vy, vz, mass (depending on nproperties)
        """
        # Construct metadata file path by adding .txt extension
        if verbose:
            print(f"Reading lightcone data from {path}")
        txt_path = path + ".txt"
        metadata = LightconeReader.rd_metadata(txt_path, verbose=verbose)
        
        npart = metadata['npart']
        
        # Read the raw data
        with open(path, 'rb') as f:
            # Read particle IDs first (8 bytes each)
            # idp_data = np.frombuffer(f.read(0 * npart), dtype=np.int32) # use this version when processing old output without idp
            idp_data = np.frombuffer(f.read(4 * npart), dtype=np.int32)
            
            # Read the remaining properties (positions, velocities, masses) (4 bytes each)
            real_data = np.frombuffer(f.read(4 * nproperties * npart), dtype=np.float32)
            real_data = real_data.reshape(nproperties, npart)
        
        return idp_data, real_data

    @staticmethod
    def rd_positions(path, verbose=False):
        """
        Read only the particle positions (x, y, z) from a lightcone shell.

        Seeks past the idp block instead of loading it (saves both the int32
        array and the corresponding file I/O). The on-disk layout is idp
        (int32, npart) followed by the float32 properties stored row-major, so
        the first three property rows are x, y, z.

        Returns:
            positions: float32 array with shape (3, npart), rows x, y, z.
        """
        if verbose:
            print(f"Reading positions from {path}")
        txt_path = path + ".txt"
        metadata = LightconeReader.rd_metadata(txt_path, verbose=verbose)
        npart = metadata['npart']

        with open(path, 'rb') as f:
            f.seek(4 * npart)  # skip idp (int32 per particle)
            real_data = np.frombuffer(f.read(4 * 3 * npart), dtype=np.float32)
            real_data = real_data.reshape(3, npart)

        return real_data

    @staticmethod
    def rd_cell(path, nproperties=8, verbose=False):
        """
        Read the lightcone shell from the output directory.
        nproperties: number of properties per cell (default 9 for x,y,z,rho,phi,accelx,accely,accelz,dphidt)
        
        Returns:
            properties: numpy array of properties (float32) with shape (nproperties, ncell)
                       where rows are x, y, z, rho, phi, accelx, accely, accelz, dphidt (depending on nproperties)
        """
        # Construct metadata file path by adding .txt extension
        if verbose:
            print(f"Reading lightcone data from {path}")
        txt_path = path + ".txt"
        metadata = LightconeReader.rd_metadata(txt_path, verbose=verbose)
        
        ncell = metadata['npart']
        
        # Read the raw data
        with open(path, 'rb') as f:

            # Read the remaining properties (positions, velocities, masses) (4 bytes each)
            real_data = np.frombuffer(f.read(4 * nproperties * ncell), dtype=np.float32)
            real_data = real_data.reshape(nproperties, ncell)
        
        return real_data

    @staticmethod
    def rd_positions_as_healpix(path, nside, verbose=False):
        """
        Read the lightcone shell from the output directory and convert it to a Healpix map.
        nside: Healpix resolution parameter
        """
        # Read only the position data (x, y, z), seeking past idp.
        properties = LightconeReader.rd_positions(path, verbose=verbose)
        x, y, z = properties[0], properties[1], properties[2]

        # Convert Cartesian coordinates to spherical coordinates.
        # x is the depth (cone axis), y and z are the transverse coordinates.
        # Free intermediates as soon as they are dead so peak memory is a
        # single map-sized (npix) array rather than several particle-sized ones.
        r = np.sqrt(x**2 + y**2 + z**2)
        theta = np.arccos(z / r)  # polar angle from x-axis
        del r
        phi = np.arctan2(y, x)
        del properties, x, y, z

        # Convert spherical coordinates to HEALPix pixel indices
        npix = hp.nside2npix(nside)
        pix_indices = hp.ang2pix(nside, theta, phi)
        del theta, phi

        # Count particles per pixel. bincount(minlength=npix) yields exactly the
        # same per-pixel counts as np.unique(return_counts=True) + assignment,
        # but without the O(npart) sort temporaries. Cast to float32 and drop
        # the int64 counts immediately.
        healpix_map = np.bincount(pix_indices, minlength=npix).astype(np.float32)
        del pix_indices

        return healpix_map

    @staticmethod
    def get_shells(path, verbose=False):
        """
        Get all lightcone shell information from the lightcone directory.
        Looks for files named 'part_xxxxx' and 'tree_xxxxx' and returns shell information.
        
        Args:
            path: Path to the lightcone directory
            verbose: Print debug information
            
        Returns:
            List of dictionaries with shell information, sorted by nstep in descending order
            (largest nout corresponds to shell closest to observer)
            
            Each dictionary contains:
            - 'nstep': shell number (int)
            - 'part_file': path to part binary file (str, or None if doesn't exist)
            - 'part_metadata': path to part .txt file (str, or None if doesn't exist)
            - 'part_size': size of part binary file in bytes (int, or None if doesn't exist)
            - 'tree_file': path to tree binary file (str, or None if doesn't exist)
            - 'tree_metadata': path to tree .txt file (str, or None if doesn't exist)
            - 'tree_size': size of tree binary file in bytes (int, or None if doesn't exist)
            - 'grav_file': path to grav binary file (str, or None if doesn't exist)
            - 'grav_metadata': path to grav .txt file (str, or None if doesn't exist)
            - 'grav_size': size of grav binary file in bytes (int, or None if doesn't exist)
        """
        # Check if path exists
        if not os.path.exists(path):
            if verbose:
                print(f"Path {path} does not exist")
            return []
        
        # Define patterns for each file type
        patterns = {
            'part_file': re.compile(r'^part_(\d{5})$'),
            'part_metadata': re.compile(r'^part_(\d{5})\.txt$'),
            'tree_file': re.compile(r'^tree_(\d{5})$'),
            'tree_metadata': re.compile(r'^tree_(\d{5})\.txt$'),
            'grav_file': re.compile(r'^grav_(\d{5})$'),
            'grav_metadata': re.compile(r'^grav_(\d{5})\.txt$')
        }
        
        shells = {}  # Dictionary to collect shell information by nstep
        
        try:
            # Loop over all files in the directory
            for filename in os.listdir(path):
                filepath = os.path.join(path, filename)
                if not os.path.isfile(filepath):
                    continue
                
                # Check each pattern
                for file_type, pattern in patterns.items():
                    match = pattern.match(filename)
                    if match:
                        nstep = int(match.group(1))
                        
                        # Initialize shell entry if needed
                        if nstep not in shells:
                            shells[nstep] = {'nstep': nstep}
                        
                        # Store file path
                        shells[nstep][file_type] = filepath
                        
                        # Store file size for binary files
                        if file_type in ['part_file', 'tree_file']:
                            try:
                                shells[nstep][file_type.replace('_file', '_size')] = os.path.getsize(filepath)
                                if verbose:
                                    print(f"Found {file_type} shell {nstep}: {os.path.getsize(filepath)/1024**2:.2f} MB")
                            except OSError:
                                if verbose:
                                    print(f"Warning: Could not get size for {file_type} shell {nstep}")
                        break
                        
        except OSError as e:
            if verbose:
                print(f"Error reading directory {path}: {e}")
            return []
        
        # Convert to list and sort by nstep in descending order
        shell_list = list(shells.values())
        shell_list.sort(key=lambda x: x['nstep'], reverse=True)
        
        # Fill in None values for missing fields
        for shell in shell_list:
            for field in ['part_file', 'part_metadata', 'part_size', 'tree_file', 'tree_metadata', 'tree_size']:
                if field not in shell:
                    shell[field] = None
        
        # Print statistics only in verbose mode
        if verbose and shell_list:
            part_shells = [s for s in shell_list if s['part_file'] is not None]
            tree_shells = [s for s in shell_list if s['tree_file'] is not None]
            grav_shells = [s for s in shell_list if s['grav_file'] is not None]
            
            print(f"Found {len(shell_list)} total shells ({len(part_shells)} part, {len(tree_shells)} tree)")
            
            if part_shells:
                total_part_size = sum(s['part_size'] for s in part_shells if s['part_size'] is not None)
                print(f"Part files total size: {total_part_size/1024**3:.2f} GB")
            
            if tree_shells:
                total_tree_size = sum(s['tree_size'] for s in tree_shells if s['tree_size'] is not None)
                print(f"Tree files total size: {total_tree_size/1024**3:.2f} GB")

            if grav_shells:
                total_grav_size = sum(s['grav_size'] for s in grav_shells if s['grav_size'] is not None)
                print(f"Grav files total size: {total_grav_size/1024**3:.2f} GB")

        return shell_list

    @staticmethod
    def overdensity_hp_map(shell, nside, info, **kwargs):
        """ 
        This function computes the overdensity map for a given shell using HEALPix.

        Parameters
        ----------
        shell : dict
            A dictionary containing the properties of the shell, including its
            comoving distance, width, and scale factor.
        nside : int
            The HEALPix nside parameter for the output map.
        info : object
            An object containing cosmological parameters, including omega_m and c.

        Returns
        -------
        overdensity : ndarray
            The computed overdensity map as a HEALPix map.
        """

        # Compute the comoving distances in the shell, but only if the caller
        # (e.g. convergence_hp_map) has not already set them on the shell dict.
        if 'comoving_distance' not in shell or 'comoving_width' not in shell:
            cosmo = LambdaCDM(info.H0, info.omega_m, info.omega_l)
            shell['comoving_distance'] = cosmo.comoving_distance(shell['z_center']).value  / info.unit_l.to(u.Mpc).value
            shell['comoving_width'] = (cosmo.comoving_distance(shell['z_old']) - cosmo.comoving_distance(shell['z'])).value / info.unit_l.to(u.Mpc).value

        solid_angle_per_pix = (4 * np.pi / hp.nside2npix(nside)) # in steradians

        # in code units
        rho_bar = 1.0
        m_part = 1 / info.npart

        # float32 per-pixel count map; fold rho = m_part*counts/dOmega/dchi/chi^2
        # and overdensity = rho/rho_bar - 1 as in-place scalar ops so peak memory
        # is a single npix-sized array (no full-size temporaries).
        overdensity = LightconeReader.rd_positions_as_healpix(shell['part_file'], nside)
        scale = np.float32(
            m_part / solid_angle_per_pix
            / shell['comoving_width'] / shell['comoving_distance']**2 / rho_bar
        )
        overdensity *= scale
        overdensity -= np.float32(1.0)

        return overdensity

    @staticmethod
    def _shell_kappa_contribution(shell, nside, info, chi_j):
        """
        Compute a single shell's contribution to the convergence map.

        This is the body of the (embarrassingly parallel) loop in
        convergence_hp_map, factored out into a module-level function so it can be
        dispatched to a multiprocessing.Pool. On Linux the worker processes are
        forked, so the module globals (overdensity_hp_map, info, np, hp) are
        inherited from the parent and don't need to be passed in.
        """
        # overdensity_hp_map reuses the comoving distances that convergence_hp_map
        # already set on the shell dict (no redundant cosmology recompute here).
        delta_k = LightconeReader.overdensity_hp_map(shell, nside, info)
        dchi_k, chi_k = shell['comoving_width'], shell['comoving_distance']
        a_k = shell['aexp_center']
        # Fold the lensing-kernel factor in place on the float32 map.
        kernel = np.float32(dchi_k * chi_k * (chi_j - chi_k) / chi_j / a_k)
        delta_k *= kernel
        return delta_k

    @staticmethod
    def convergence_hp_map(shells, nside, info, **kwargs):
        """
        Compute the convergence map from a list of shells.

        Parameters
        ----------
        shells : list of dict
            Each dict contains the properties of a shell, including its
            comoving distance, width, and scale factor.
        nside : int
            The HEALPix nside parameter for the output map.
        info : object
            An object containing cosmological parameters, including omega_m and c.
        nproc : int, optional
            The number of processes to use for parallel computation. Default is None,
            which uses the number of CPU cores available.
        
        Returns
        -------
        kappa_map : ndarray
            The computed convergence map as a HEALPix map.
        """

        import multiprocessing as mp
        from functools import partial
        from tqdm import tqdm

        # python 3.14+ have default forkserver (aka reimport everything, we just want fork here)
        ctx = mp.get_context('fork')

        # compute the comoving distances in each shell
        cosmo = LambdaCDM(info.H0, info.omega_m, info.omega_l)
        for shell in shells:
            shell['comoving_distance'] = cosmo.comoving_distance(shell['z_center']).value  / info.unit_l.to(u.Mpc).value
            shell['comoving_width'] = (cosmo.comoving_distance(shell['z_old']) - cosmo.comoving_distance(shell['z'])).value / info.unit_l.to(u.Mpc).value

        nproc = kwargs.get('nproc', None)  # None -> os.cpu_count() workers
        n_workers = nproc or os.cpu_count()
        if verbose := kwargs.get('verbose', True):
            print(f"Computing convergence map with nside={nside}, using {n_workers} processes")

        # float64 accumulator for precision; per-shell contributions arrive as
        # float32 maps and are summed in here.
        kappa_map = np.zeros(hp.nside2npix(nside), dtype=np.float64)

        chi_j = shells[-1]['comoving_distance'] # the comoving distance to the last shell in sequence

        worker = partial(LightconeReader._shell_kappa_contribution, nside=nside, info=info, chi_j=chi_j)
        with ctx.Pool(processes=nproc) as pool:
            # chunksize=1 (the imap_unordered default) is what enables dynamic,
            # work-stealing dispatch across the oversubscribed chunk list.
            for contribution in tqdm(pool.imap_unordered(worker, shells),
                                     total=len(shells), desc="Accumulating kappa map"):
                kappa_map += contribution

        kappa_map *= 3 * info.omega_m / 2 / info.c**2  # in code units

        return kappa_map

    @staticmethod 
    def cone_hp_mask(nside, alpha_y, alpha_z):
        """
        Rectangular pencil-beam mask centered on the +x axis (lon=0, lat=0).
        Note: nside = 1024 takes ~ 1 second, O(nside^2) scaling. We compute the 
        corners of the cone by solving the conditions at equality for a unit sphere 
        x^2 + y^2 + z^2 = 1 then 
    
        x / sqrt(x^2 + y^2) = cos(alpha_y) and x / sqrt(x^2 + z^2) = cos(alpha_z). 
    
        This gives the corners needed for the healpix query_polygon function.
    
        Parameters
        -------------
        nside (int): 
            healpix nside parameter
        alpha_y (float): 
            half opening angle in longitude (x-y projection), degrees
        alpha_z (float): 
            half opening angle in latitude  (x-z projection), degrees
            
        Returns 
        ----------
        mask (np.ndarray): 
            healpix map with 1 inside the beam, 0 outside.
        """
        
        ty, tz = np.tan(np.radians(alpha_y)), np.tan(np.radians(alpha_z))
        x = 1.0 / np.sqrt(1 + ty**2 + tz**2)
        corners = np.array([[x,  x*ty,  x*tz],
                            [x, -x*ty,  x*tz],
                            [x, -x*ty, -x*tz],
                            [x,  x*ty, -x*tz]])   # ordered around the boundary
        ipix = hp.query_polygon(nside, corners)
        mask = np.zeros(hp.nside2npix(nside))
        mask[ipix] = 1.0
        return mask


