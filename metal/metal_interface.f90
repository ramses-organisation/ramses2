module metal_interface
  use iso_c_binding
  implicit none

  interface

    ! Initialize Metal device and command queue.
    subroutine mtl_init() bind(C, name="mtl_init")
    end subroutine mtl_init

    function mtl_nsubgrid() bind(C, name="mtl_nsubgrid")
      import c_int
      integer(c_int) :: mtl_nsubgrid
    end function mtl_nsubgrid

    ! Allocate Metal buffers for uold, unew, grid, nbor, hash_key, hash_val.
    ! ncachemax is included so grid/nbor cover the full ngridmax+ncachemax range
    ! needed for AMR ghost-zone caching.
    subroutine mtl_alloc_amr(ngridmax, ncachemax, nvar, twotondim, hash_size) &
        bind(C, name="mtl_alloc_amr")
      import c_int
      integer(c_int), value :: ngridmax, ncachemax, nvar, twotondim, hash_size
    end subroutine mtl_alloc_amr

    ! Allocate AMR refinement device buffers (flag1/2, father, sort arrays,
    ! cache pointers, per-level Hilbert parameters).
    ! Mirrors the flag/sort allocation block in gpu_allocate_amr (gpu_manager.cuf).
    subroutine mtl_alloc_refine(ngridmax, ncachemax, nlevelmax) &
        bind(C, name="mtl_alloc_refine")
      import c_int
      integer(c_int), value :: ngridmax, ncachemax, nlevelmax
    end subroutine mtl_alloc_refine

    ! Copy per-level Hilbert parameters from host to device.
    ! Called once from metal_allocate_amr after init_amr populates the arrays.
    subroutine mtl_upload_level_params(ckey_max, key_off, &
        box_ckey_min, box_ckey_max, periodic, nlevelmax) &
        bind(C, name="mtl_upload_level_params")
      import c_ptr, c_int
      type(c_ptr), value     :: ckey_max       ! integer(c_int)(1:nlevelmax+1)
      type(c_ptr), value     :: key_off        ! integer(c_long)(1:nlevelmax+1)
      type(c_ptr), value     :: box_ckey_min   ! integer(c_int)(1:ndim,1:nlevelmax+1)
      type(c_ptr), value     :: box_ckey_max   ! integer(c_int)(1:ndim,1:nlevelmax+1)
      integer(c_int), intent(in) :: periodic(3)
      integer(c_int), value  :: nlevelmax
    end subroutine mtl_upload_level_params

    ! Copy host arrays into Metal buffers (H->D).
    ! Mirrors r_set_grid_device / cudaMemcpy in the CUDA path.
    ! On Apple Silicon the memcpy stays within DRAM; no PCIe transfer.
    subroutine mtl_set_grid_device(uold, unew, bold, &
        grid, ngridmax, nvar, twotondim) &
        bind(C, name="mtl_set_grid_device")
      import c_ptr, c_int
      type(c_ptr), value  :: uold   ! real(kind=4)(1:twotondim,1:nvar,1:ngridmax)
      type(c_ptr), value  :: unew   ! real(kind=4)(1:twotondim,1:nvar,1:ngridmax)
      type(c_ptr), value  :: bold
      type(c_ptr), value  :: grid   ! type(oct)(1:ngridmax)
      integer(c_int), value :: ngridmax, nvar, twotondim
    end subroutine mtl_set_grid_device

    ! Upload host flag1(8,ngridmax) to device. Mirrors CUDA `flag1 = pst%s%m%flag1`.
    ! Must be called from r_set_grid_device when nlevelmax > levelmin.
    subroutine mtl_upload_flag1(flag1, ngridmax) bind(C, name="mtl_upload_flag1")
      import c_ptr, c_int
      type(c_ptr),    value :: flag1
      integer(c_int), value :: ngridmax
    end subroutine mtl_upload_flag1

    ! Copy Metal uold buffer back to host (D->H), called before I/O.
    ! Mirrors r_transfer_grid_host / cudaMemcpy in the CUDA path.
    subroutine mtl_transfer_grid_host(uold, bold, &
        ngridmax, nvar, twotondim) &
        bind(C, name="mtl_transfer_grid_host")
      import c_ptr, c_int
      type(c_ptr), value  :: uold
      type(c_ptr), value  :: bold
      integer(c_int), value :: ngridmax, nvar, twotondim
    end subroutine mtl_transfer_grid_host

    ! Copy Metal s_grid buffer (oct ckey/refined/lev) back to host (D->H).
    ! Required for AMR runs where metal_refine reorders octs on the device.
    subroutine mtl_transfer_grid_struct_host(grid, ngridmax) &
        bind(C, name="mtl_transfer_grid_struct_host")
      import c_ptr, c_int
      type(c_ptr), value    :: grid
      integer(c_int), value :: ngridmax
    end subroutine mtl_transfer_grid_struct_host

    ! Dispatch set_unew_kernel: unew(:,:,oct) = uold(:,:,oct) for octs at ilevel.
    subroutine mtl_set_unew(head_idx, num_octs) bind(C, name="mtl_set_unew")
      import c_int
      integer(c_int), value :: head_idx   ! 1-based index of first oct
      integer(c_int), value :: num_octs
    end subroutine mtl_set_unew

    ! Dispatch set_uold_kernel: uold(:,:,oct) = unew(:,:,oct) for octs at ilevel.
    subroutine mtl_set_uold(head_idx, num_octs) bind(C, name="mtl_set_uold")
      import c_int
      integer(c_int), value :: head_idx
      integer(c_int), value :: num_octs
    end subroutine mtl_set_uold

    ! Dispatch cmpdt_kernel and return integrated quantities + Courant dt.
    ! Mirrors gpu_cmpdt in gpu_runner.cuf; output args passed by reference.
    subroutine mtl_cmpdt(head_idx, num_octs, dx, gamma, smallr, smallc2, &
        courant_factor, &
#ifdef MHD
        induction, &
#endif
        constant_gravity, mass, ekin, eint, emag, dt) &
        bind(C, name="mtl_cmpdt")
      import c_int, c_float
      integer(c_int), value :: head_idx, num_octs
      real(c_float),  value :: dx, gamma, smallr, smallc2, courant_factor
#ifdef MHD
      integer(c_int), value :: induction
#endif
      real(c_float), intent(in) :: constant_gravity(3)
      real(c_float), intent(out) :: mass, ekin, eint, emag, dt
    end subroutine mtl_cmpdt

    ! Build nbor array from the already-populated hash table.
    ! Mirrors the update_nbor_array loop in r_set_grid_device
    ! (gpu_manager.cuf); a single dispatch replaces those launches.
    subroutine mtl_build_nbor(head_idx, num_subgrids, hash_size, &
        ckey_max_l, key_off_l, &
        box_ckey_min, box_ckey_max, periodic_i) &
        bind(C, name="mtl_build_nbor")
      import c_int, c_long
      integer(c_int), value :: head_idx, num_subgrids
      integer(c_int), value :: hash_size
      integer(c_int), value :: ckey_max_l
      integer(c_long), value :: key_off_l
      integer(c_int), intent(in) :: box_ckey_min(3)
      integer(c_int), intent(in) :: box_ckey_max(3)
      integer(c_int), intent(in) :: periodic_i(3)
    end subroutine mtl_build_nbor

    ! Inclusive int32 prefix scan of s_prefix_sum[offset..offset+n-1].
    ! offset is 0-based (Fortran head_idx - 1).  Three internal phases.
    subroutine mtl_prefix_scan(offset, n) bind(C, name="mtl_prefix_scan")
      import c_int
      integer(c_int), value :: offset, n
    end subroutine mtl_prefix_scan

    ! Read prefix_sum[offset + n - 1] after mtl_prefix_scan.
    ! Mirrors get_total_sum in gpu_scan.cuf.
    function mtl_get_prefix_total(offset, n) bind(C, name="mtl_get_prefix_total")
      import c_int
      integer(c_int)        :: mtl_get_prefix_total
      integer(c_int), value :: offset, n
    end function mtl_get_prefix_total

    ! ---- Flag kernels (AMR refinement flagging) ----------------------------

    ! Batch: reset_flag1 [+ init_flag if noct_fine>0] + count_flag1 in one cmd buf.
    ! Replaces metal_init_flag's 3 separate commit/waits with 1.
    function mtl_init_flag_batch(head_coarse, noct_coarse, head_fine, noct_fine) &
        bind(C, name="mtl_init_flag_batch")
      import c_int
      integer(c_int)        :: mtl_init_flag_batch
      integer(c_int), value :: head_coarse, noct_coarse, head_fine, noct_fine
    end function mtl_init_flag_batch

    ! Batch: poisson_flag + hydro_flag + count_flag1 in one cmd buf.
    ! Replaces metal_user_flag's separate calls with 1 batch.
    function mtl_user_flag_batch(head_idx, num_octs, &
        gamma, smallr, smallc2, err_grad_d, err_grad_p, floor_d, floor_p, &
        mass_sph, m_refine, jeans_refine, factG, dx_loc &
#ifdef MHD
        , err_grad_b2, floor_b2, err_grad_A, floor_A, &
        err_grad_B, floor_B, err_grad_C, floor_C &
#endif
        ) &
        bind(C, name="mtl_user_flag_batch")
      import c_int, c_float
      integer(c_int)        :: mtl_user_flag_batch
      integer(c_int), value :: head_idx, num_octs
      real(c_float),  value :: gamma, smallr, smallc2
      real(c_float),  value :: err_grad_d, err_grad_p, floor_d, floor_p
      real(c_float),  value :: mass_sph, m_refine, jeans_refine, factG, dx_loc
#ifdef MHD
      real(c_float), value :: err_grad_b2, floor_b2, err_grad_A, floor_A
      real(c_float), value :: err_grad_B, floor_B, err_grad_C, floor_C
#endif
    end function mtl_user_flag_batch

    subroutine mtl_enforce_subgrid(head_idx, num_octs) bind(C, name="mtl_enforce_subgrid")
      import c_int
      integer(c_int), value :: head_idx, num_octs
    end subroutine mtl_enforce_subgrid

    ! Batch: 3×(count_neighbors + flag_count) + count_flag1 in one cmd buf.
    ! n_nbor = {1,2,2} hardcoded for NDIM=3.
    ! Replaces metal_smooth_flag's 7 separate commit/waits with 1.
    function mtl_smooth_flag_batch(head_idx, num_octs) &
        bind(C, name="mtl_smooth_flag_batch")
      import c_int
      integer(c_int)        :: mtl_smooth_flag_batch
      integer(c_int), value :: head_idx, num_octs
    end function mtl_smooth_flag_batch

    ! Zero flag1 for octs [head_idx .. head_idx+num_octs-1].
    subroutine mtl_reset_flag1(head_idx, num_octs) bind(C, name="mtl_reset_flag1")
      import c_int
      integer(c_int), value :: head_idx, num_octs
    end subroutine mtl_reset_flag1

    ! Zero flag2 for octs [head_idx .. head_idx+num_octs-1].
    subroutine mtl_reset_flag2(head_idx, num_octs) bind(C, name="mtl_reset_flag2")
      import c_int
      integer(c_int), value :: head_idx, num_octs
    end subroutine mtl_reset_flag2

    ! Flag parent cell when any fine child is refined or flagged (ilevel+1 octs).
    subroutine mtl_init_flag(head_idx, num_octs) bind(C, name="mtl_init_flag")
      import c_int
      integer(c_int), value :: head_idx, num_octs
    end subroutine mtl_init_flag

    ! Atomic-reduce sum of flag1 cells; returns count.
    function mtl_count_flag1(head_idx, num_octs) bind(C, name="mtl_count_flag1")
      import c_int
      integer(c_int)        :: mtl_count_flag1
      integer(c_int), value :: head_idx, num_octs
    end function mtl_count_flag1

    ! Gradient-based hydro and MHD refinement criteria (GRAV=0).
    subroutine mtl_hydro_flag(head_idx, num_octs, &
        gamma, smallr, smallc2, err_grad_d, err_grad_p, floor_d, floor_p &
#ifdef MHD
        , err_grad_b2, floor_b2, err_grad_A, floor_A, &
        err_grad_B, floor_B, err_grad_C, floor_C &
#endif
        ) &
        bind(C, name="mtl_hydro_flag")
      import c_int, c_float
      integer(c_int), value :: head_idx, num_octs
      real(c_float),  value :: gamma, smallr, smallc2
      real(c_float),  value :: err_grad_d, err_grad_p, floor_d, floor_p
#ifdef MHD
      real(c_float), value :: err_grad_b2, floor_b2, err_grad_A, floor_A
      real(c_float), value :: err_grad_B, floor_B, err_grad_C, floor_C
#endif
    end subroutine mtl_hydro_flag

    ! Count flagged face-adjacent neighbours → flag2.
    subroutine mtl_count_neighbors(head_idx, num_octs) &
        bind(C, name="mtl_count_neighbors")
      import c_int
      integer(c_int), value :: head_idx, num_octs
    end subroutine mtl_count_neighbors

    ! Promote flag1 if flag2 >= num_nbors; clear flag2 if flag1 already set.
    subroutine mtl_flag_count(head_idx, num_octs, num_nbors) &
        bind(C, name="mtl_flag_count")
      import c_int
      integer(c_int), value :: head_idx, num_octs, num_nbors
    end subroutine mtl_flag_count

    ! Clear flag1 if any 3×3×3 neighbour slot is absent or in ghost cache.
    subroutine mtl_enforce_rules(head_idx, num_octs, ngridmax) &
        bind(C, name="mtl_enforce_rules")
      import c_int
      integer(c_int), value :: head_idx, num_octs, ngridmax
    end subroutine mtl_enforce_rules

    ! Build father[] array: father[oct] = 1-based parent oct via hash lookup.
    ! ckey_max_l / key_off_l are the Hilbert params for the PARENT level.
    subroutine mtl_build_father(head_idx, num_octs, hash_size, &
        ckey_max_l, key_off_l) bind(C, name="mtl_build_father")
      import c_int, c_long
      integer(c_int), value :: head_idx, num_octs, hash_size, ckey_max_l
      integer(c_long), value :: key_off_l
    end subroutine mtl_build_father

    ! ---- AMR refine / sort / cache bridges --------------------------------

    ! Write/read ifree_dev (device free-oct counter).
    subroutine mtl_set_ifree(val) bind(C, name="mtl_set_ifree")
      import c_int
      integer(c_int), value :: val
    end subroutine mtl_set_ifree

    function mtl_get_ifree() bind(C, name="mtl_get_ifree")
      import c_int
      integer(c_int) :: mtl_get_ifree
    end function mtl_get_ifree

    ! Write/read ifree_cache_dev (device cache pointer).
    subroutine mtl_set_ifree_cache(val) bind(C, name="mtl_set_ifree_cache")
      import c_int
      integer(c_int), value :: val
    end subroutine mtl_set_ifree_cache

    function mtl_get_ifree_cache() bind(C, name="mtl_get_ifree_cache")
      import c_int
      integer(c_int) :: mtl_get_ifree_cache
    end function mtl_get_ifree_cache

    ! Increment ifree_cache_dev in-place (CPU-side, unified memory).
    subroutine mtl_advance_ifree_cache(new_noct) &
        bind(C, name="mtl_advance_ifree_cache")
      import c_int
      integer(c_int), value :: new_noct
    end subroutine mtl_advance_ifree_cache

    ! Create child octs for flagged unrefined cells (refine_kernel).
    subroutine mtl_refine_cells(head_idx, num_octs, hash_size, interpol_var, interpol_type, smallr) &
        bind(C, name="mtl_refine_cells")
      import c_int, c_float
      integer(c_int), value :: head_idx, num_octs, hash_size, interpol_var, interpol_type
      real(c_float), value :: smallr
    end subroutine mtl_refine_cells

    ! Reset hash table every coarse step: zero + re-insert all octs + cache.
    ! Mirrors gpu_reset_hash in gpu_runner.cuf.
    subroutine mtl_reset_hash(ifree, ngridmax, ifree_cache, hash_size) &
        bind(C, name="mtl_reset_hash")
      import c_int
      integer(c_int), value :: ifree, ngridmax, ifree_cache, hash_size
    end subroutine mtl_reset_hash

    ! Free child octs whose parent cell is no longer flagged.
    subroutine mtl_derefine_cells(head_idx, num_octs, hash_size) &
        bind(C, name="mtl_derefine_cells")
      import c_int
      integer(c_int), value :: head_idx, num_octs, hash_size
    end subroutine mtl_derefine_cells

    ! Wipe hash entries for a range of octs (cache cleanup).
    subroutine mtl_free_hash_range(head_idx, num_octs, hash_size) &
        bind(C, name="mtl_free_hash_range")
      import c_int
      integer(c_int), value :: head_idx, num_octs, hash_size
    end subroutine mtl_free_hash_range

    ! Overwrite hash val for each oct after sort rearrangement.
    subroutine mtl_update_hash_range(head_idx, num_octs, hash_size) &
        bind(C, name="mtl_update_hash_range")
      import c_int
      integer(c_int), value :: head_idx, num_octs, hash_size
    end subroutine mtl_update_hash_range

    ! Insert all octs in [head..head+n-1] into hash (reads level from grid).
    subroutine mtl_insert_hash_all(head_idx, num_octs, hash_size) &
        bind(C, name="mtl_insert_hash_all")
      import c_int
      integer(c_int), value :: head_idx, num_octs, hash_size
    end subroutine mtl_insert_hash_all

    ! Initialise swap_global to identity permutation.
    subroutine mtl_init_swap_table(head_idx, num_octs) &
        bind(C, name="mtl_init_swap_table")
      import c_int
      integer(c_int), value :: head_idx, num_octs
    end subroutine mtl_init_swap_table

    ! Set prefix_sum[oct] = (lev != ilevel).
    subroutine mtl_init_prefix_level(head_idx, num_octs, ilevel) &
        bind(C, name="mtl_init_prefix_level")
      import c_int
      integer(c_int), value :: head_idx, num_octs, ilevel
    end subroutine mtl_init_prefix_level

    ! Batch all num_bits Hilbert-sort passes in one command buffer.
    ! Replaces the inner loop over mtl_init_prefix_bit/scan/local_swap/global_swap.
    subroutine mtl_hilbert_sort_level(head_idx, num_octs, num_bits) &
        bind(C, name="mtl_hilbert_sort_level")
      import c_int
      integer(c_int), value :: head_idx, num_octs, num_bits
    end subroutine mtl_hilbert_sort_level

    ! Set prefix_sum[oct] = bit ibit of oct's Hilbert key.
    subroutine mtl_init_prefix_bit(head_idx, num_octs, ibit) &
        bind(C, name="mtl_init_prefix_bit")
      import c_int
      integer(c_int), value :: head_idx, num_octs, ibit
    end subroutine mtl_init_prefix_bit

    ! LSD scatter: write swap_local using prefix_sum + swap_global.
    subroutine mtl_compute_local_swap(head_idx, num_octs) &
        bind(C, name="mtl_compute_local_swap")
      import c_int
      integer(c_int), value :: head_idx, num_octs
    end subroutine mtl_compute_local_swap

    ! Apply local swap to global permutation.
    subroutine mtl_update_global_swap(head_idx, num_octs) &
        bind(C, name="mtl_update_global_swap")
      import c_int
      integer(c_int), value :: head_idx, num_octs
    end subroutine mtl_update_global_swap

    ! Pack lev/ckey/refined-mask from grid (via swap_global) into flag2 scratch.
    subroutine mtl_sort_gather_grid(head_idx, num_octs) &
        bind(C, name="mtl_sort_gather_grid")
      import c_int
      integer(c_int), value :: head_idx, num_octs
    end subroutine mtl_sort_gather_grid

    ! Unpack flag2 scratch back to grid and recompute hkey via Hilbert.
    subroutine mtl_sort_scatter_grid(head_idx, num_octs) &
        bind(C, name="mtl_sort_scatter_grid")
      import c_int
      integer(c_int), value :: head_idx, num_octs
    end subroutine mtl_sort_scatter_grid

    ! flag2[:,oct] = flag1[:,swap_global[oct]-1].
    subroutine mtl_sort_gather_flag(head_idx, num_octs) &
        bind(C, name="mtl_sort_gather_flag")
      import c_int
      integer(c_int), value :: head_idx, num_octs
    end subroutine mtl_sort_gather_flag

    ! flag1[:,oct] = flag2[:,oct].
    subroutine mtl_sort_scatter_flag(head_idx, num_octs) &
        bind(C, name="mtl_sort_scatter_flag")
      import c_int
      integer(c_int), value :: head_idx, num_octs
    end subroutine mtl_sort_scatter_flag

    ! unew[cell,ivar,oct] = uold[cell,ivar,swap_global[oct]-1].
    subroutine mtl_sort_gather_hydro(head_idx, num_octs) &
        bind(C, name="mtl_sort_gather_hydro")
      import c_int
      integer(c_int), value :: head_idx, num_octs
    end subroutine mtl_sort_gather_hydro

    ! Blit unew[head..head+n-1] → uold (replaces cudaMemcpy after sort).
    subroutine mtl_blit_unew_to_uold(head_idx, num_octs) &
        bind(C, name="mtl_blit_unew_to_uold")
      import c_int
      integer(c_int), value :: head_idx, num_octs
    end subroutine mtl_blit_unew_to_uold

    subroutine mtl_sort_gather_force(head_idx, num_octs, idim) &
        bind(C, name="mtl_sort_gather_force")
      import c_int
      integer(c_int), value :: head_idx, num_octs, idim
    end subroutine mtl_sort_gather_force

    subroutine mtl_sort_scatter_force(head_idx, num_octs, idim) &
        bind(C, name="mtl_sort_scatter_force")
      import c_int
      integer(c_int), value :: head_idx, num_octs, idim
    end subroutine mtl_sort_scatter_force

    subroutine mtl_sort_gather_phi(head_idx, num_octs, is_phi_old) &
        bind(C, name="mtl_sort_gather_phi")
      import c_int
      integer(c_int), value :: head_idx, num_octs, is_phi_old
    end subroutine mtl_sort_gather_phi

    subroutine mtl_sort_scatter_phi(head_idx, num_octs, is_phi_old) &
        bind(C, name="mtl_sort_scatter_phi")
      import c_int
      integer(c_int), value :: head_idx, num_octs, is_phi_old
    end subroutine mtl_sort_scatter_phi

    ! Run update_nbor_prefix + all prefix-scan phases in one command buffer.
    ! Returns total number of subgrids needing a cache oct for direction input_ind.
    ! Replaces mtl_update_nbor_prefix + mtl_prefix_scan + mtl_get_prefix_total.
    function mtl_nbor_scan(head_idx, num_subgrids, hash_size, input_ind) &
        bind(C, name="mtl_nbor_scan")
      import c_int
      integer(c_int)        :: mtl_nbor_scan
      integer(c_int), value :: head_idx, num_subgrids, hash_size, input_ind
    end function mtl_nbor_scan

    ! Run compute_cache_swap + make_cache_octs + insert_hash_cache in one command buffer.
    ! Replaces three separate commit/wait calls per cache direction.
    subroutine mtl_cache_fill(head_idx, num_subgrids, hash_size, input_ind, &
        ngridmax, ifree_cache, new_noct, interpol_var, interpol_type, smallr) bind(C, name="mtl_cache_fill")
      import c_int, c_float
      integer(c_int), value :: head_idx, num_subgrids, hash_size, input_ind
      integer(c_int), value :: ngridmax, ifree_cache, new_noct
      integer(c_int), value :: interpol_var, interpol_type
      real(c_float), value :: smallr
    end subroutine mtl_cache_fill

    ! Set nbor[input_ind, sg] and prefix_sum[sg] = (nbor==0 ? 1 : 0).
    subroutine mtl_update_nbor_prefix(head_idx, num_subgrids, &
        hash_size, input_ind) bind(C, name="mtl_update_nbor_prefix")
      import c_int
      integer(c_int), value :: head_idx, num_subgrids, hash_size, input_ind
    end subroutine mtl_update_nbor_prefix

    ! Build swap_local of subgrids that need a new cache oct.
    subroutine mtl_compute_cache_swap(head_idx, num_subgrids) &
        bind(C, name="mtl_compute_cache_swap")
      import c_int
      integer(c_int), value :: head_idx, num_subgrids
    end subroutine mtl_compute_cache_swap

    ! Create ghost octs (cache region), interpolating MHD state when enabled.
    subroutine mtl_make_cache_octs(head_idx, num_subgrids, hash_size, &
        ngridmax, ifree_cache, new_noct, input_ind, interpol_var, interpol_type, smallr) &
        bind(C, name="mtl_make_cache_octs")
      import c_int, c_float
      integer(c_int), value :: head_idx, num_subgrids, hash_size
      integer(c_int), value :: ngridmax, ifree_cache, new_noct, input_ind
      integer(c_int), value :: interpol_var, interpol_type
      real(c_float), value :: smallr
    end subroutine mtl_make_cache_octs

    ! Insert new cache octs into the hash table.
    subroutine mtl_insert_hash_cache_r(hash_size, ngridmax, &
        ifree_cache, new_noct) bind(C, name="mtl_insert_hash_cache_r")
      import c_int
      integer(c_int), value :: hash_size, ngridmax, ifree_cache, new_noct
    end subroutine mtl_insert_hash_cache_r

    ! Restriction (averaging down) for fine→coarse levels (upload_kernel).
    ! Mirrors gpu_upload in gpu_runner.cuf.
    ! internal_energy=1: average internal energy instead of total.
    subroutine mtl_upload(head_idx, num_octs, internal_energy, &
        gamma, smallr, smallc2) bind(C, name="mtl_upload")
      import c_int, c_float
      integer(c_int), value :: head_idx, num_octs, internal_energy
      real(c_float),  value :: gamma, smallr, smallc2
    end subroutine mtl_upload

    ! Flush and wait for all pending Metal work to complete.
    ! Mirrors cudaDeviceSynchronize() in the CUDA path; called by m_timer
    ! before reading the wall clock so GPU time is captured correctly.
    subroutine mtl_device_sync() bind(C, name="mtl_device_sync")
    end subroutine mtl_device_sync

    ! Dispatch hydro_integrator_kernel (MUSCL-Hancock Godunov).
    ! Mirrors gpu_godunov in gpu_runner.cuf.
    ! constant_gravity: 3-element array (all zeros for GRAV=0).
    subroutine mtl_godunov(head_idx, num_subgrids, ngridmax, &
        ilevel, levelmin, levelmax, &
        gamma, smallr, smallc2, dt, dx, slope, &
#ifdef MHD
        slope_mag, riemann, riemann2d, switch_llf_dmin, switch_llf_pmin, induction, etamag, &
#else
        riemann, &
#endif
        constant_gravity) &
        bind(C, name="mtl_godunov")
      import c_int, c_float
      integer(c_int), value :: head_idx, num_subgrids, ngridmax
      integer(c_int), value :: ilevel, levelmin, levelmax
      real(c_float),  value :: gamma, smallr, smallc2, dt, dx
#ifdef MHD
      integer(c_int), value :: slope, slope_mag, riemann, riemann2d, induction
      real(c_float), value :: switch_llf_dmin, switch_llf_pmin, etamag
#else
      integer(c_int), value :: slope, riemann
#endif
      real(c_float), intent(in) :: constant_gravity(3)
    end subroutine mtl_godunov

    ! -------------------------------------------------------------------------
    ! Gravity / Poisson solver bindings
    ! -------------------------------------------------------------------------

    ! Allocate all gravity-related device buffers (AMR + MG).
    subroutine mtl_alloc_grav(ngridmax, ncachemax, ngridmax_mg, ncachemax_mg, &
        hash_size_mg) bind(C, name="mtl_alloc_grav")
      import c_int
      integer(c_int), value :: ngridmax, ncachemax, ngridmax_mg, ncachemax_mg
      integer(c_int), value :: hash_size_mg
    end subroutine mtl_alloc_grav

    ! Upload rho array (cell densities for all octs) to device.
    subroutine mtl_upload_rho(rho_host, head_idx, num_octs) bind(C, name="mtl_upload_rho")
      import c_ptr, c_int
      type(c_ptr), value :: rho_host
      integer(c_int), value :: head_idx, num_octs
    end subroutine mtl_upload_rho

    ! Reset rho to zero (reset_rho_kernel, 2D {8,16}).
    subroutine mtl_reset_rho(head_idx, num_octs) bind(C, name="mtl_reset_rho")
      import c_int
      integer(c_int), value :: head_idx, num_octs
    end subroutine mtl_reset_rho

    ! Compute multipole leaf contributions (multipole_leaf_kernel, 2D {8,16}).
    subroutine mtl_multipole_leaf(head_idx, num_octs, scale) &
        bind(C, name="mtl_multipole_leaf")
      import c_int, c_float
      integer(c_int), value :: head_idx, num_octs
      real(c_float),  value :: scale
    end subroutine mtl_multipole_leaf

    ! Upload multipole to coarser level (multipole_upload_kernel, 1D 256).
    subroutine mtl_multipole_upload(head_idx, num_octs) &
        bind(C, name="mtl_multipole_upload")
      import c_int
      integer(c_int), value :: head_idx, num_octs
    end subroutine mtl_multipole_upload

    ! Reduce total multipole (multipole_tot_kernel, 1D 256); result in *tot4.
    subroutine mtl_multipole_tot(head_idx, num_octs, tot4) &
        bind(C, name="mtl_multipole_tot")
      import c_int, c_float
      integer(c_int), value    :: head_idx, num_octs
      real(c_float), intent(out) :: tot4(4)
    end subroutine mtl_multipole_tot

    ! CIC density deposit (deposit_rho_kernel, 2D {8,16}).
    subroutine mtl_deposit_rho(head_idx, num_octs, dx, vol_loc, &
        m_refine, mass_sph, var_cut_refine, ivar_refine, ngridmax) &
        bind(C, name="mtl_deposit_rho")
      import c_int, c_float
      integer(c_int), value :: head_idx, num_octs
      real(c_float),  value :: dx, vol_loc, m_refine, mass_sph, var_cut_refine
      integer(c_int), value :: ivar_refine, ngridmax
    end subroutine mtl_deposit_rho
    ! Initialise phi from time-extrapolated coarser-level interpolation.
    subroutine mtl_init_phi(head_idx, num_octs, ngridmax, tfrac) &
        bind(C, name="mtl_init_phi")
      import c_int, c_float
      integer(c_int), value :: head_idx, num_octs, ngridmax
      real(c_float),  value :: tfrac
    end subroutine mtl_init_phi

    ! Compute gravitational force gradient (gradient_phi_kernel, 2D {8,16}).
    ! "fine" variant uses AMR phi/f; "mg" uses MG phi/f.
    subroutine mtl_gradient_phi_fine(head_idx, num_octs, dx) &
        bind(C, name="mtl_gradient_phi_fine")
      import c_int, c_float
      integer(c_int), value :: head_idx, num_octs
      real(c_float),  value :: dx
    end subroutine mtl_gradient_phi_fine

    subroutine mtl_gradient_phi_mg(head_idx, num_octs, dx) &
        bind(C, name="mtl_gradient_phi_mg")
      import c_int, c_float
      integer(c_int), value :: head_idx, num_octs
      real(c_float),  value :: dx
    end subroutine mtl_gradient_phi_mg

    ! Compute potential energy (cmp_epot_kernel, 1D 256); result in epot_out.
    subroutine mtl_cmp_epot(head_idx, num_octs, epot_out) &
        bind(C, name="mtl_cmp_epot")
      import c_int, c_float
      integer(c_int), value    :: head_idx, num_octs
      real(c_float), intent(out) :: epot_out
    end subroutine mtl_cmp_epot

    ! Compute maximum density (cmp_rhomax_kernel, 1D 256); result in rhomax_out.
    subroutine mtl_cmp_rhomax(head_idx, num_octs, rhomax_out) &
        bind(C, name="mtl_cmp_rhomax")
      import c_int, c_float
      integer(c_int), value    :: head_idx, num_octs
      real(c_float), intent(out) :: rhomax_out
    end subroutine mtl_cmp_rhomax

    ! Zero MG hash tables (removes tombstones before each multigrid cycle).
    subroutine mtl_clean_mg_hashes() bind(C, name="mtl_clean_mg_hashes")
    end subroutine mtl_clean_mg_hashes

    ! Build the MG level (prefix sum → father octs → hash → nbor).
    ! "fine": source octs from AMR grid; "mg": source octs from MG grid.
    subroutine mtl_build_mg_fine(ifine, ilevel, head_idx, num_octs, &
        head_father, head_mg, new_noct) bind(C, name="mtl_build_mg_fine")
      import c_int
      integer(c_int), value  :: ifine, ilevel, head_idx, num_octs, head_father, head_mg
      integer(c_int), intent(out) :: new_noct
    end subroutine mtl_build_mg_fine

    subroutine mtl_build_mg_mg(ifine, ilevel, head_idx, num_octs, &
        head_father, head_mg, new_noct) bind(C, name="mtl_build_mg_mg")
      import c_int
      integer(c_int), value  :: ifine, ilevel, head_idx, num_octs, head_father, head_mg
      integer(c_int), intent(out) :: new_noct
    end subroutine mtl_build_mg_mg

    ! Reset mask to mask_val in f(cell,3,oct) (reset_mask_kernel, 2D {8,16}).
    ! "fine" binds s_f; "mg" binds s_f_mg.
    subroutine mtl_reset_mask_fine(head_idx, num_octs, mask_val) &
        bind(C, name="mtl_reset_mask_fine")
      import c_int, c_float
      integer(c_int), value :: head_idx, num_octs
      real(c_float),  value :: mask_val
    end subroutine mtl_reset_mask_fine

    subroutine mtl_reset_mask_mg(head_idx, num_octs, mask_val) &
        bind(C, name="mtl_reset_mask_mg")
      import c_int, c_float
      integer(c_int), value :: head_idx, num_octs
      real(c_float),  value :: mask_val
    end subroutine mtl_reset_mask_mg

    ! Restrict mask from fine → coarser MG level (restrict_mask_kernel, 1D 128).
    subroutine mtl_restrict_mask_fine(head_idx, head_father, num_octs) &
        bind(C, name="mtl_restrict_mask_fine")
      import c_int
      integer(c_int), value :: head_idx, head_father, num_octs
    end subroutine mtl_restrict_mask_fine

    subroutine mtl_restrict_mask_mg(head_idx, head_father, num_octs) &
        bind(C, name="mtl_restrict_mask_mg")
      import c_int
      integer(c_int), value :: head_idx, head_father, num_octs
    end subroutine mtl_restrict_mask_mg

    ! Convert volume accumulation to signed mask (volume_to_mask_kernel, 1D 256).
    subroutine mtl_volume_to_mask_fine(head_idx, num_octs, mask_max_out) &
        bind(C, name="mtl_volume_to_mask_fine")
      import c_int, c_float
      integer(c_int), value    :: head_idx, num_octs
      real(c_float), intent(out) :: mask_max_out
    end subroutine mtl_volume_to_mask_fine

    subroutine mtl_volume_to_mask_mg(head_idx, num_octs, mask_max_out) &
        bind(C, name="mtl_volume_to_mask_mg")
      import c_int, c_float
      integer(c_int), value    :: head_idx, num_octs
      real(c_float), intent(out) :: mask_max_out
    end subroutine mtl_volume_to_mask_mg

    ! Compute residual (cmp_residual_kernel, 2D {8,16}).
    subroutine mtl_cmp_residual_fine(head_idx, num_octs, ngridmax, &
        fourpi, offset, oneoverdx2) bind(C, name="mtl_cmp_residual_fine")
      import c_int, c_float
      integer(c_int), value :: head_idx, num_octs, ngridmax
      real(c_float),  value :: fourpi, offset, oneoverdx2
    end subroutine mtl_cmp_residual_fine

    subroutine mtl_cmp_residual_mg(head_idx, num_octs, ngridmax, &
        fourpi, offset, oneoverdx2) bind(C, name="mtl_cmp_residual_mg")
      import c_int, c_float
      integer(c_int), value :: head_idx, num_octs, ngridmax
      real(c_float),  value :: fourpi, offset, oneoverdx2
    end subroutine mtl_cmp_residual_mg

    ! Gauss-Seidel smoother red/black step (gauss_seidel_kernel, 2D {4,32}).
    subroutine mtl_gauss_seidel_fine(head_idx, num_octs, ngridmax, &
        dx2, safe, redstep) bind(C, name="mtl_gauss_seidel_fine")
      import c_int, c_float
      integer(c_int), value :: head_idx, num_octs, ngridmax
      real(c_float),  value :: dx2
      integer(c_int), value :: safe, redstep
    end subroutine mtl_gauss_seidel_fine

    subroutine mtl_gauss_seidel_mg(head_idx, num_octs, ngridmax, &
        dx2, safe, redstep) bind(C, name="mtl_gauss_seidel_mg")
      import c_int, c_float
      integer(c_int), value :: head_idx, num_octs, ngridmax
      real(c_float),  value :: dx2
      integer(c_int), value :: safe, redstep
    end subroutine mtl_gauss_seidel_mg

    ! Restrict residual from fine → coarser MG level (restrict_residual_kernel).
    subroutine mtl_restrict_residual_fine(head_idx, head_father, num_octs) &
        bind(C, name="mtl_restrict_residual_fine")
      import c_int
      integer(c_int), value :: head_idx, head_father, num_octs
    end subroutine mtl_restrict_residual_fine

    subroutine mtl_restrict_residual_mg(head_idx, head_father, num_octs) &
        bind(C, name="mtl_restrict_residual_mg")
      import c_int
      integer(c_int), value :: head_idx, head_father, num_octs
    end subroutine mtl_restrict_residual_mg


    ! Interpolate coarse correction and add to fine grid.
    subroutine mtl_interpolate_correct_fine(head_idx, head_father, num_octs) &
        bind(C, name="mtl_interpolate_correct_fine")
      import c_int
      integer(c_int), value :: head_idx, head_father, num_octs
    end subroutine mtl_interpolate_correct_fine

    subroutine mtl_interpolate_correct_mg(head_idx, head_father, num_octs) &
        bind(C, name="mtl_interpolate_correct_mg")
      import c_int
      integer(c_int), value :: head_idx, head_father, num_octs
    end subroutine mtl_interpolate_correct_mg

    ! L2 residual norm (residual_norm_kernel, 1D 256); result in norm_out.
    subroutine mtl_residual_norm_fine(head_idx, num_octs, norm_out) &
        bind(C, name="mtl_residual_norm_fine")
      import c_int, c_float
      integer(c_int), value    :: head_idx, num_octs
      real(c_float), intent(out) :: norm_out
    end subroutine mtl_residual_norm_fine

    subroutine mtl_rhs_norm_fine(head_idx, num_octs, norm_out) &
        bind(C, name="mtl_rhs_norm_fine")
      import c_int, c_float
      integer(c_int), value    :: head_idx, num_octs
      real(c_float), intent(out) :: norm_out
    end subroutine mtl_rhs_norm_fine

    subroutine mtl_residual_norm_mg(head_idx, num_octs, norm_out) &
        bind(C, name="mtl_residual_norm_mg")
      import c_int, c_float
      integer(c_int), value    :: head_idx, num_octs
      real(c_float), intent(out) :: norm_out
    end subroutine mtl_residual_norm_mg

    ! Save phi → phi_old (save_phi_old_kernel, 2D {8,16}).
    subroutine mtl_save_phi_old_fine(head_idx, num_octs) &
        bind(C, name="mtl_save_phi_old_fine")
      import c_int
      integer(c_int), value :: head_idx, num_octs
    end subroutine mtl_save_phi_old_fine

    ! Reset phi (and f[1..3]) to zero (reset_phi_kernel_mg, 2D {8,16}).
    subroutine mtl_reset_phi_fine(head_idx, num_octs) &
        bind(C, name="mtl_reset_phi_fine")
      import c_int
      integer(c_int), value :: head_idx, num_octs
    end subroutine mtl_reset_phi_fine

    subroutine mtl_reset_phi_mg(head_idx, num_octs) &
        bind(C, name="mtl_reset_phi_mg")
      import c_int
      integer(c_int), value :: head_idx, num_octs
    end subroutine mtl_reset_phi_mg

    ! Set phi to a uniform value (reset_phi_val_kernel, 2D {8,16}).
    subroutine mtl_reset_phi_val_fine(head_idx, num_octs, phi_val) &
        bind(C, name="mtl_reset_phi_val_fine")
      import c_int, c_float
      integer(c_int), value :: head_idx, num_octs
      real(c_float),  value :: phi_val
    end subroutine mtl_reset_phi_val_fine

    subroutine mtl_reset_phi_val_mg(head_idx, num_octs, phi_val) &
        bind(C, name="mtl_reset_phi_val_mg")
      import c_int, c_float
      integer(c_int), value :: head_idx, num_octs
      real(c_float),  value :: phi_val
    end subroutine mtl_reset_phi_val_mg

    ! Set RHS in f(cell,2,oct) (reset_rhs_kernel_mg, 2D {8,16}).
    subroutine mtl_reset_rhs_fine(head_idx, num_octs, ngridmax, &
        fourpi, offset, oneoverdx2) bind(C, name="mtl_reset_rhs_fine")
      import c_int, c_float
      integer(c_int), value :: head_idx, num_octs, ngridmax
      real(c_float),  value :: fourpi, offset, oneoverdx2
    end subroutine mtl_reset_rhs_fine

    subroutine mtl_reset_rhs_mg(head_idx, num_octs, ngridmax, &
        fourpi, offset, oneoverdx2) bind(C, name="mtl_reset_rhs_mg")
      import c_int, c_float
      integer(c_int), value :: head_idx, num_octs, ngridmax
      real(c_float),  value :: fourpi, offset, oneoverdx2
    end subroutine mtl_reset_rhs_mg

    subroutine mtl_upload_cooling_table(n1, n2, nH_tbl, T2_tbl, &
        cool, heat, cool_com, heat_com, metal, cool_prime, &
        heat_prime, cool_com_prime, heat_com_prime, metal_prime) &
        bind(C, name="mtl_upload_cooling_table")
      import c_int, c_ptr
      integer(c_int), value :: n1, n2
      type(c_ptr), value :: nH_tbl, T2_tbl
      type(c_ptr), value :: cool, heat, cool_com, heat_com, metal
      type(c_ptr), value :: cool_prime, heat_prime, cool_com_prime
      type(c_ptr), value :: heat_com_prime, metal_prime
    end subroutine mtl_upload_cooling_table

    subroutine mtl_cooling(head_idx, num_octs, gamma, smallr, smallc2, &
        dtcool, eos_type, eos_T2, eos_nH, eos_index, scale_T2, scale_nH, &
        cooling, metal, imetal, z_ave, self_shielding, X_frac, T2max, isothermal) &
        bind(C, name="mtl_cooling")
      import c_int, c_double, c_float
      integer(c_int), value :: head_idx, num_octs
      real(c_float), value :: gamma, smallr, smallc2
      real(c_double), value :: dtcool
      integer(c_int), value :: eos_type
      real(c_double), value :: eos_T2, eos_nH, eos_index, scale_T2, scale_nH
      integer(c_int), value :: cooling, metal, imetal
      real(c_double), value :: z_ave
      integer(c_int), value :: self_shielding
      real(c_double), value :: X_frac, T2max
      integer(c_int), value :: isothermal
    end subroutine mtl_cooling

    subroutine mtl_sync_hydro(head_idx, num_octs, gamma, smallr, smallc2, dt, constant_gravity) &
        bind(C, name="mtl_sync_hydro")
      import c_int, c_float
      integer(c_int), value :: head_idx, num_octs
      real(c_float),  value :: gamma, smallr, smallc2, dt
      real(c_float),  intent(in) :: constant_gravity(3)
    end subroutine mtl_sync_hydro

    subroutine mtl_grav_hydro(head_idx, num_octs, gamma, smallr, smallc2, dt, constant_gravity) &
        bind(C, name="mtl_grav_hydro")
      import c_int, c_float
      integer(c_int), value :: head_idx, num_octs
      real(c_float),  value :: gamma, smallr, smallc2, dt
      real(c_float),  intent(in) :: constant_gravity(3)
    end subroutine mtl_grav_hydro

    ! Particle interfaces
    subroutine mtl_alloc_part(npartmax) bind(C, name="mtl_alloc_part")
      import c_int
      integer(c_int), value :: npartmax
    end subroutine mtl_alloc_part

    subroutine mtl_upload_part(xp, vp, mp, levelp, sortp, idp, npart) bind(C, name="mtl_upload_part")
      import c_ptr, c_int
      type(c_ptr), value :: xp, vp, mp, levelp, sortp, idp
      integer(c_int), value :: npart
    end subroutine mtl_upload_part

    subroutine mtl_download_part(xp, vp, mp, levelp, sortp, idp, npart) bind(C, name="mtl_download_part")
      import c_ptr, c_int
      type(c_ptr), value :: xp, vp, mp, levelp, sortp, idp
      integer(c_int), value :: npart
    end subroutine mtl_download_part

    subroutine mtl_kick_drift_part(action_part, ilevel, head_idx, num_parts, skip1, skip2, skip3, dx_loc, &
        box_size, periodic, dtnew, dtold) bind(C, name="mtl_kick_drift_part")
      import c_int, c_float, c_ptr
      integer(c_int), value :: action_part, ilevel, head_idx, num_parts
      real(c_float), value :: skip1, skip2, skip3, dx_loc
      type(c_ptr), value :: box_size, periodic, dtnew, dtold
    end subroutine mtl_kick_drift_part

    subroutine mtl_newdt_part(head_idx, num_parts, vmax_out, ekin_out) bind(C, name="mtl_newdt_part")
      import c_int, c_float
      integer(c_int), value :: head_idx, num_parts
      real(c_float), intent(out) :: vmax_out, ekin_out
    end subroutine mtl_newdt_part

    subroutine mtl_split_part(head_idx, num_parts, ilevel, skip1, skip2, skip3, dx_loc, n_fine_out) &
        bind(C, name="mtl_split_part")
      import c_int, c_float
      integer(c_int), value :: head_idx, num_parts, ilevel
      real(c_float), value :: skip1, skip2, skip3, dx_loc
      integer(c_int), intent(out) :: n_fine_out
    end subroutine mtl_split_part

    subroutine mtl_sort_part(head_idx, num_parts, level, shift, dx_inv, skip) bind(C, name="mtl_sort_part")
      import c_int, c_float, c_ptr
      integer(c_int), value :: head_idx, num_parts, level
      real(c_float), value :: shift, dx_inv
      type(c_ptr), value :: skip
    end subroutine mtl_sort_part

    subroutine mtl_cic_part_medium(head_idx, num_parts, skip1, skip2, skip3, dx_loc, vol_loc, mass_sph, &
        star, m_refine_at_level, mass_cut_refine, ilevel) bind(C, name="mtl_cic_part_medium")
      import c_int, c_float
      integer(c_int), value :: head_idx, num_parts, ilevel, star
      real(c_float), value :: skip1, skip2, skip3, dx_loc, vol_loc, mass_sph, m_refine_at_level, mass_cut_refine
    end subroutine mtl_cic_part_medium

    subroutine mtl_multipole_q_part(head_idx, num_parts, leading, q_out) bind(C, name="mtl_multipole_q_part")
      import c_int, c_long, c_float
      integer(c_int),  value :: head_idx, num_parts
      integer(c_long), value :: leading
      real(c_float), intent(out) :: q_out(4)
    end subroutine mtl_multipole_q_part

    subroutine mtl_debug_rho(head_grid, num_octs, twotondim, ilevel) bind(C, name="mtl_debug_rho")
      import c_int
      integer(c_int), value :: head_grid, num_octs, twotondim, ilevel
    end subroutine mtl_debug_rho

    subroutine mtl_debug_xp(head_idx, num_parts, ilevel) bind(C, name="mtl_debug_xp")
      import c_int
      integer(c_int), value :: head_idx, num_parts, ilevel
    end subroutine mtl_debug_xp

  end interface

end module metal_interface
