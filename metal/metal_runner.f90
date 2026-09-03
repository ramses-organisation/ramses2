module metal_runner
  use metal_interface
  use iso_c_binding
  implicit none

contains

integer function metal_nsubgrid()
  metal_nsubgrid = int(mtl_nsubgrid())
end function metal_nsubgrid

!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine metal_allocate_amr(sim)
  use amr_parameters, only: twotondim, ndim
  use hydro_parameters, only: nvar
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t) :: sim
  integer :: ilevel
  integer(c_int) :: periodic_i(3)

  ! Set hash_size
  sim%m%hash_size = 2 * (sim%m%ngridmax + sim%m%ncachemax)

  ! Allocate key offsets
  allocate(sim%m%key_off(1:sim%r%nlevelmax+1))
  sim%m%key_off(1) = 1_8
  do ilevel = 2, sim%r%nlevelmax + 1
     sim%m%key_off(ilevel) = sim%m%key_off(ilevel-1) + sim%m%hkey_max(1, ilevel-1)
  end do

  ! Allocate cache pointers
  allocate(sim%m%head_cache(1:sim%r%nlevelmax))
  allocate(sim%m%tail_cache(1:sim%r%nlevelmax))
  allocate(sim%m%noct_cache(1:sim%r%nlevelmax))
  sim%m%head_cache=1
  sim%m%tail_cache=0
  sim%m%noct_cache=0

  ! Allocate Metal-owned buffers for uold/unew/grid/nbor/hash.
  ! ncachemax is now passed so grid and nbor cover the full
  ! ngridmax+ncachemax range needed for AMR ghost-zone caching.
  call mtl_alloc_amr( &
       int(sim%m%ngridmax,   c_int), &
       int(sim%m%ncachemax,  c_int), &
       int(nvar,             c_int), &
       int(twotondim,        c_int), &
       int(sim%m%hash_size,  c_int))

  ! Allocate AMR refinement device buffers (flag1/2, father, sort, cache,
  ! per-level params).
  call mtl_alloc_refine( &
       int(sim%m%ngridmax,  c_int), &
       int(sim%m%ncachemax, c_int), &
       int(sim%r%nlevelmax, c_int))

  ! Initialise ifree_cache (CUDA path does this in gpu_allocate_amr).
  sim%m%ifree_cache = 1

  ! Upload per-level Hilbert parameters (ckey_max, key_off, box bounds,
  ! periodicity) to device.  These are populated by init_amr before this
  ! subroutine is called, and key_off was just allocated above.
  periodic_i = merge(int(1, c_int), int(0, c_int), sim%r%periodic)
  call mtl_upload_level_params( &
       c_loc(sim%m%ckey_max(1)),          &
       c_loc(sim%m%key_off(1)),           &
       c_loc(sim%m%box_ckey_min(1,1)),    &
       c_loc(sim%m%box_ckey_max(1,1)),    &
       periodic_i,                         &
       int(sim%r%nlevelmax, c_int))

end subroutine metal_allocate_amr

!###########################################################
!###########################################################
!###########################################################
!###########################################################
recursive subroutine r_set_grid_device(pst)
  use mdl_module
  use mdl_parameters, only: MDL_SET_GRID_DEVICE
  use amr_parameters, only: twotondim
  use hydro_parameters, only: nvar
  use ramses_commons, only: pst_t
  implicit none
  type(pst_t) :: pst

  integer :: rID

  if(pst%nLower > 0) then
     rID = mdl_send_request(pst%s%mdl, MDL_SET_GRID_DEVICE, pst%iUpper+1)
     call r_set_grid_device(pst%pLower)
     call mdl_get_reply(pst%s%mdl, rID, 0)
  else
     ! Copy host arrays into Metal buffers (mirrors H->D cudaMemcpy in gpu_manager.cuf).
#ifdef HYDRO
     call mtl_set_grid_device( &
          c_loc(pst%s%m%uold(1,1,1)), &
          c_loc(pst%s%m%unew(1,1,1)), &
#ifdef MHD
          c_loc(pst%s%m%bold(1,1,1)), &
#else
          c_null_ptr, &
#endif
          c_loc(pst%s%m%grid(1)),     &
          int(pst%s%m%ngridmax, c_int), &
          int(nvar,             c_int), &
          int(twotondim,        c_int))
#else
     call mtl_set_grid_device( &
          c_null_ptr, &
          c_null_ptr, &
          c_null_ptr, &
          c_loc(pst%s%m%grid(1)),     &
          int(pst%s%m%ngridmax, c_int), &
          int(nvar,             c_int), &
          int(twotondim,        c_int))
#endif
     if (pst%s%r%nlevelmax > pst%s%r%levelmin) then
        call mtl_upload_flag1( &
             c_loc(pst%s%m%flag1(1,1)), &
             int(pst%s%m%ngridmax, c_int))
     end if

     ! Build nbor on device via insert_hash + build_nbor kernels
     call metal_set_nbor(pst%s, pst%s%r%levelmin)

     ! Copy particles to device
     if (pst%s%r%part) then
        call metal_upload_part(pst%s)
     end if
     pst%s%m%data_on_device = .true.
  endif

end subroutine r_set_grid_device

!###########################################################
!###########################################################
!###########################################################
!###########################################################
recursive subroutine r_transfer_grid_host(pst)
  use mdl_module
  use mdl_parameters, only: MDL_TRANSFER_GRID_HOST
  use amr_parameters, only: twotondim
  use hydro_parameters, only: nvar
  use ramses_commons, only: pst_t
  implicit none
  type(pst_t) :: pst

  integer :: rID

  if(pst%nLower > 0) then
     rID = mdl_send_request(pst%s%mdl, MDL_TRANSFER_GRID_HOST, pst%iUpper+1)
     call r_transfer_grid_host(pst%pLower)
     call mdl_get_reply(pst%s%mdl, rID, 0)
  else
#ifdef HYDRO
     call mtl_transfer_grid_host( &
          c_loc(pst%s%m%uold(1,1,1)), &
#ifdef MHD
          c_loc(pst%s%m%bold(1,1,1)), &
#else
          c_null_ptr, &
#endif
          int(pst%s%m%ngridmax, c_int), &
          int(nvar,             c_int), &
          int(twotondim,        c_int))
#endif
     if (pst%s%r%nlevelmax > pst%s%r%levelmin) then
        call mtl_transfer_grid_struct_host( &
             c_loc(pst%s%m%grid(1)), &
             int(pst%s%m%ngridmax, c_int))
     end if
     if (pst%s%r%part) then
        call metal_download_part(pst%s)
     end if
  endif

end subroutine r_transfer_grid_host

!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine metal_set_unew(sim, ilevel)
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t), intent(inout) :: sim
  integer, intent(in) :: ilevel

  call mtl_set_unew( &
       int(sim%m%head(ilevel), c_int), &
       int(sim%m%noct(ilevel), c_int))

end subroutine metal_set_unew

!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine metal_set_uold(sim, ilevel)
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t), intent(inout) :: sim
  integer, intent(in) :: ilevel

  call mtl_set_uold( &
       int(sim%m%head(ilevel), c_int), &
       int(sim%m%noct(ilevel), c_int))

end subroutine metal_set_uold

!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine metal_cmpdt(sim, ilevel, mass, ekin, eint, emag, dt)
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t), intent(inout) :: sim
  integer,        intent(in)    :: ilevel
  real(kind=8),   intent(out)   :: mass, ekin, eint, emag, dt

  real(c_float) :: dx, gamma, smallr, smallc2, courant_factor
  real(c_float) :: constant_gravity(3)
  real(c_float) :: mass_f, ekin_f, eint_f, emag_f, dt_f
#ifdef MHD
  integer(c_int) :: induction
#endif

  dx               = real(sim%r%boxlen / 2**ilevel, c_float)
  gamma            = real(sim%r%gamma,              c_float)
  smallr           = real(sim%r%smallr,             c_float)
  smallc2          = real(sim%r%smallc**2,          c_float)
  courant_factor   = real(sim%r%courant_factor,     c_float)
  constant_gravity = real(sim%r%constant_gravity,   c_float)
#ifdef MHD
  induction = merge(int(1, c_int), int(0, c_int), sim%r%induction)
#endif

  call mtl_cmpdt(                      &
       int(sim%m%head(ilevel), c_int), &
       int(sim%m%noct(ilevel), c_int), &
       dx, gamma, smallr, smallc2,     &
       courant_factor,                 &
#ifdef MHD
       induction,                      &
#endif
       constant_gravity,               &
       mass_f, ekin_f, eint_f, emag_f, dt_f)

  mass = real(mass_f, 8)
  ekin = real(ekin_f, 8)
  eint = real(eint_f, 8)
  emag = real(emag_f, 8)
  dt   = real(dt_f,   8)

end subroutine metal_cmpdt

!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine metal_godunov(sim, ilevel)
  use ramses_commons, only: ramses_t
#ifdef MHD
  use amr_parameters, only: ndim
  use hydro_parameters, only: solver_llf, solver_hlld, solver_uct_hlld, solver2d_llf, solver2d_hlld
#endif
  implicit none
  type(ramses_t), intent(inout) :: sim
  integer,        intent(in)    :: ilevel

  real(c_float) :: gamma, smallr, smallc2, dt, dx
  real(c_float) :: constant_gravity(3)
#ifdef MHD
  integer(c_int) :: head_idx, num_subgrids, nsubgrid_l, nsubgridtondim_l, induction
#endif

  gamma            = real(sim%r%gamma,              c_float)
  smallr           = real(sim%r%smallr,             c_float)
  smallc2          = real(sim%r%smallc**2,          c_float)
  dt               = real(sim%g%dtnew(ilevel),      c_float)
  dx               = real(sim%r%boxlen / 2**ilevel, c_float)
  constant_gravity = real(sim%r%constant_gravity,   c_float)

#ifdef MHD
  nsubgrid_l = mtl_nsubgrid()
  nsubgridtondim_l = nsubgrid_l**ndim
  if (sim%r%riemann == solver_hlld .or. sim%r%riemann == solver_llf) then
     if (sim%r%riemann2d /= solver2d_hlld .and. sim%r%riemann2d /= solver2d_llf) &
          error stop 'Metal MHD requires riemann2d=hlld or llf'
  else if (sim%r%riemann /= solver_uct_hlld) then
     error stop 'Metal MHD requires riemann=hlld, llf, or uct-hlld'
  end if
  if (mod(sim%m%head(ilevel)-1, nsubgridtondim_l) /= 0) error stop 'Metal MHD oct head is not nsubgrid aligned'
  if (mod(sim%m%noct(ilevel), nsubgridtondim_l) /= 0) error stop 'Metal MHD oct count is not divisible by nsubgrid'
  head_idx = int((sim%m%head(ilevel)-1)/nsubgridtondim_l+1, c_int)
  num_subgrids = int(sim%m%noct(ilevel)/nsubgridtondim_l, c_int)
  induction = merge(int(1, c_int), int(0, c_int), sim%r%induction)
#endif

  if(sim%r%verbose .and. sim%g%myid==1) &
       write(*,'("   Entering metal_godunov for level ",I2)') ilevel

  call mtl_godunov(                    &
#ifdef MHD
       head_idx, num_subgrids,          &
#else
       int(sim%m%head(ilevel), c_int), &
       int(sim%m%noct(ilevel), c_int), &
#endif
       int(sim%m%ngridmax,     c_int), &
       int(ilevel,             c_int), &
       int(sim%r%levelmin,     c_int), &
       int(sim%r%nlevelmax,    c_int), &
       gamma, smallr, smallc2,         &
       dt, dx,                         &
       int(sim%r%slope_type,   c_int), &
#ifdef MHD
       int(sim%r%slope_mag_type, c_int), &
       int(sim%r%riemann,      c_int), &
       int(sim%r%riemann2d,    c_int), &
       real(sim%r%switch_llf_dmin, c_float), &
       real(sim%r%switch_llf_pmin, c_float), &
       induction, real(sim%r%etamag, c_float), &
#else
       int(sim%r%riemann,      c_int), &
#endif
       constant_gravity)

end subroutine metal_godunov

!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine metal_set_nbor(sim, ilevel)
  ! Insert all oct Hilbert keys into the device hash table, then build
  ! the device nbor array for octs at ilevel.
  ! NOTE: father[] is NOT populated here — it is built inside metal_refine
  use amr_parameters, only: ndim
  use ramses_commons,  only: ramses_t
  use iso_c_binding
  implicit none
  type(ramses_t), intent(inout) :: sim
  integer,        intent(in)    :: ilevel

  integer(c_int)  :: hash_size_l, ckey_max_l
  integer(c_long) :: key_off_l
  integer(c_int)  :: bmin(3), bmax(3), periodic_i(3)
#ifdef MHD
  integer(c_int) :: head_idx, num_subgrids, nsubgrid_l, nsubgridtondim_l
#endif

  hash_size_l = int(sim%m%hash_size,                   c_int)
  ckey_max_l  = int(sim%m%ckey_max(ilevel),            c_int)
  key_off_l   = int(sim%m%key_off(ilevel),             c_long)
  bmin        = int(sim%m%box_ckey_min(1:ndim,ilevel), c_int)
  bmax        = int(sim%m%box_ckey_max(1:ndim,ilevel), c_int)
  periodic_i  = merge(int(1, c_int), int(0, c_int), sim%r%periodic)
#ifdef MHD
  nsubgrid_l = mtl_nsubgrid()
  nsubgridtondim_l = nsubgrid_l**ndim
  if (mod(sim%m%head(ilevel)-1, nsubgridtondim_l) /= 0) error stop 'Metal MHD oct head is not nsubgrid aligned'
  if (mod(sim%m%noct(ilevel), nsubgridtondim_l) /= 0) error stop 'Metal MHD oct count is not divisible by nsubgrid'
  head_idx = int((sim%m%head(ilevel)-1)/nsubgridtondim_l+1, c_int)
  num_subgrids = int(sim%m%noct(ilevel)/nsubgridtondim_l, c_int)
#endif

  ! Step 1: insert all allocated octs into the hash table (1 .. ifree-1).
  ! Use insert_hash_all (reads grid[].lev per oct) so fine octs get their own
  ! level's ckey_max/key_off rather than the coarse-level scalars.
  call mtl_insert_hash_all( &
       int(1,             c_int), &
       int(sim%m%ifree-1, c_int), &
       hash_size_l)

  ! Step 2: build nbor for octs at ilevel
  call mtl_build_nbor( &
#ifdef MHD
       head_idx, num_subgrids, &
#else
       int(sim%m%head(ilevel), c_int), &
       int(sim%m%noct(ilevel), c_int), &
#endif
       hash_size_l, ckey_max_l, key_off_l, &
       bmin, bmax, periodic_i)

end subroutine metal_set_nbor

!###########################################################
!###########################################################
!###########################################################
!###########################################################

subroutine metal_init_flag(sim, ilevel, nflag)
  ! 1. Zero flag1 for all cells at ilevel.
  ! 2. For each fine oct (ilevel+1), flag its parent cell if any child
  !    is refined or already flagged (uses father[] built by metal_refine).
  ! 3. Return total flagged-cell count.
  use ramses_commons, only: ramses_t
  use iso_c_binding
  implicit none
  type(ramses_t), intent(inout) :: sim
  integer,        intent(in)    :: ilevel
  integer,        intent(out)   :: nflag

  nflag = int(mtl_init_flag_batch( &
       int(sim%m%head(ilevel), c_int), &
       int(sim%m%noct(ilevel), c_int), &
       int(sim%m%head(ilevel+1), c_int), &
       int(merge(sim%m%noct(ilevel+1), 0, &
       ilevel < sim%r%nlevelmax .and. sim%m%noct(ilevel+1) > 0), c_int)))

end subroutine metal_init_flag

!###########################################################
!###########################################################
!###########################################################
!###########################################################

subroutine metal_user_flag(sim, ilevel, nflag)
  ! Applies gradient-based density/pressure refinement criterion then
  ! returns updated flagged-cell count.
  use ramses_commons, only: ramses_t
  use iso_c_binding
  implicit none
  type(ramses_t), intent(inout) :: sim
  integer,        intent(in)    :: ilevel
  integer,        intent(out)   :: nflag

  real(c_float) :: gamma, smallr, smallc2
  real(c_float) :: err_grad_d, err_grad_p, floor_d, floor_p
#ifdef MHD
  real(c_float) :: err_grad_b2, floor_b2, err_grad_A, floor_A
  real(c_float) :: err_grad_B, floor_B, err_grad_C, floor_C
#endif
  real(kind=8)  :: dx, factG

  dx = sim%r%boxlen/2**ilevel
  factG = 1.0d0
  if (sim%r%cosmo) factG = 3.0d0 / 8.0d0 / acos(-1.0d0) * sim%g%omega_m * sim%g%aexp

  gamma      = real(sim%r%gamma,      c_float)
  smallr     = real(sim%r%smallr,     c_float)
  smallc2    = real(sim%r%smallc**2,  c_float)
  err_grad_d = real(sim%r%err_grad_d, c_float)
  err_grad_p = real(sim%r%err_grad_p, c_float)
  floor_d    = real(sim%r%floor_d,    c_float)
  floor_p    = real(sim%r%floor_p,    c_float)
#ifdef MHD
  err_grad_b2 = real(sim%r%err_grad_b2, c_float)
  floor_b2 = real(sim%r%floor_b2, c_float)
  err_grad_A = real(sim%r%err_grad_A, c_float)
  floor_A = real(sim%r%floor_A, c_float)
  err_grad_B = real(sim%r%err_grad_B, c_float)
  floor_B = real(sim%r%floor_B, c_float)
  err_grad_C = real(sim%r%err_grad_C, c_float)
  floor_C = real(sim%r%floor_C, c_float)
#endif

  nflag = int(mtl_user_flag_batch( &
       int(sim%m%head(ilevel), c_int), &
       int(sim%m%noct(ilevel), c_int), &
       gamma, smallr, smallc2,          &
       err_grad_d, err_grad_p, floor_d, floor_p, &
       real(sim%r%mass_sph, c_float), &
       real(sim%r%m_refine(ilevel), c_float), &
       real(sim%r%jeans_refine(ilevel), c_float), &
       real(factG, c_float), &
       real(dx, c_float) &
#ifdef MHD
       , err_grad_b2, floor_b2, err_grad_A, floor_A, &
       err_grad_B, floor_B, err_grad_C, floor_C &
#endif
       ))

end subroutine metal_user_flag

subroutine metal_enforce_subgrid(sim, ilevel)
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t), intent(inout) :: sim
  integer, intent(in) :: ilevel
  if (sim%m%noct(ilevel) > 0 .and. mtl_nsubgrid() > 1) then
     call mtl_enforce_subgrid(int(sim%m%head(ilevel), c_int), int(sim%m%noct(ilevel), c_int))
  end if
end subroutine metal_enforce_subgrid

!###########################################################
!###########################################################
!###########################################################
!###########################################################

subroutine metal_enforce_rules(sim, ilevel)
  ! Clears flag1 for any oct whose 3x3x3 nbor stencil contains a
  ! missing or ghost-cache entry (enforces 2:1 refinement constraint).
  use ramses_commons, only: ramses_t
  use iso_c_binding
  implicit none
  type(ramses_t), intent(inout) :: sim
  integer,        intent(in)    :: ilevel

  if (sim%m%noct(ilevel) > 0) then
     call mtl_enforce_rules(int(sim%m%head(ilevel), c_int), &
                            int(sim%m%noct(ilevel), c_int), &
                            int(sim%m%ngridmax,     c_int))
  end if

end subroutine metal_enforce_rules

!###########################################################
!###########################################################
!###########################################################
!###########################################################

subroutine metal_smooth_flag(sim, ilevel, nflag)
  ! Performs ndim dilatation steps: each step counts flagged face-adjacent
  ! neighbours (→ flag2) then promotes flag1 when count >= n_nbor(idim).
  ! n_nbor = [1, 2, 2] for NDIM=3.
  ! Returns updated flagged-cell count.
  use amr_parameters, only: ndim
  use ramses_commons,  only: ramses_t
  use iso_c_binding
  implicit none
  type(ramses_t), intent(inout) :: sim
  integer,        intent(in)    :: ilevel
  integer,        intent(out)   :: nflag

  nflag = int(mtl_smooth_flag_batch( &
       int(sim%m%head(ilevel), c_int), &
       int(sim%m%noct(ilevel), c_int)))

end subroutine metal_smooth_flag

!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine metal_refine(sim, ilevel, nmake, nkill)
  ! Steps:
  !  1. Wipe hash entries for existing cache octs.
  !  2. refine_kernel → read back new ifree.
  !  3. insert_hash_all for newly created octs.
  !  4. derefine_kernel per level (nlevelmax → ilevel+1).
  !  5. Level bucket sort: compact [head_child..new_ifree-1] by level;
  !     update sim%m%head/noct/tail and sim%m%ifree.
  !  6. Hilbert sort per level.
  !  7. sort_gather/scatter grid/flag/hydro; blit unew→uold.
  !  8. update_hash for sorted range.
  !  9. update_father per level (mtl_build_father).
  ! 10. Per-level cache rebuild (27 neighbour directions).
  use amr_parameters, only: ndim
  use ramses_commons, only: ramses_t
  use iso_c_binding
  implicit none
  type(ramses_t), intent(inout) :: sim
  integer,        intent(in)    :: ilevel
  integer,        intent(out)   :: nmake, nkill

  integer(c_int) :: hash_size_l, old_ifree, new_ifree
  integer(c_int) :: head_child, n_all, n_rem
  integer(c_int) :: ilev, ind
  integer(c_int) :: cur_head, nones, new_noct
  integer(c_int) :: total_valid, ifree_cache_now
  integer(c_int) :: cache_noct_lev
  integer(c_int) :: new_head, new_tail
  integer(c_int) :: head_cache, head_idx, num_octs
#ifdef MHD
  integer(c_int) :: nsubgrid_l, nsubgridtondim_l, subgridsize_l
  integer(c_int) :: head_subgrid, num_subgrids
#endif

  hash_size_l = int(sim%m%hash_size, c_int)

  ! --- Step 1: wipe old cache hash entries for levels >= ilevel+1 -------------
  head_cache = int(sim%m%head_cache(ilevel+1), c_int)
  head_idx   = head_cache + int(sim%m%ngridmax, c_int)
  num_octs   = int(sim%m%ifree_cache - head_cache, c_int)
  if (num_octs > 0) then
     call mtl_free_hash_range( &
          int(head_idx, c_int), &
          int(num_octs, c_int), &
          hash_size_l)
  end if

  ! Reset cache pointer to head of cache region for level ilevel+1.
  sim%m%ifree_cache = int(head_cache)

  ! --- Step 2: refine_kernel -------------------------------------------------
  old_ifree = int(sim%m%ifree, c_int)
  call mtl_set_ifree(old_ifree)

  call mtl_refine_cells( &
       int(sim%m%head(ilevel),        c_int), &
       old_ifree - int(sim%m%head(ilevel), c_int), &
       hash_size_l, int(sim%r%interpol_var, c_int), &
       int(sim%r%interpol_type, c_int), real(sim%r%smallr, c_float))

  new_ifree = mtl_get_ifree()
  nmake     = int(new_ifree - old_ifree)

  ! --- Step 3: insert new octs into hash table -------------------------------
  call mtl_insert_hash_all( &
       old_ifree, &
       int(nmake, c_int), &
       hash_size_l)

  ! --- Step 4: derefine per level (top-down) ---------------------------------
  do ilev = int(sim%r%nlevelmax, c_int), int(ilevel + 1, c_int), -1
     if (sim%m%noct(ilev) > 0) then
        call mtl_derefine_cells( &
             int(sim%m%head(ilev), c_int), &
             int(sim%m%noct(ilev), c_int), &
             hash_size_l)
     end if
  end do

  ! --- Step 5: level bucket sort on [head_child .. new_ifree-1] --------------
  head_child = int(sim%m%tail(ilevel) + 1, c_int)
  n_all      = new_ifree - head_child   ! total slots in child region

  if (n_all > 0) then
     ! Init identity permutation.
     call mtl_init_swap_table(head_child, n_all)

     cur_head    = head_child
     total_valid = 0
     do ilev = int(ilevel + 1, c_int), int(sim%r%nlevelmax, c_int)
        n_rem = new_ifree - cur_head
        if (n_rem <= 0) exit

        ! Init prefix: bit = (lev != ilev) → 0 for octs at this level.
        call mtl_init_prefix_level(cur_head, n_rem, ilev)
        call mtl_prefix_scan(int(cur_head - 1, c_int), n_rem)
        nones   = mtl_get_prefix_total(int(cur_head - 1, c_int), n_rem)
        new_noct = n_rem - nones   ! octs at level ilev

        ! Scatter octs at ilev to front of [cur_head..cur_head+n_rem-1].
        call mtl_compute_local_swap(cur_head, n_rem)
        call mtl_update_global_swap(cur_head, n_rem)

        sim%m%head(ilev) = int(cur_head)
        sim%m%noct(ilev) = int(new_noct)
        sim%m%tail(ilev) = int(cur_head + new_noct - 1)
        cur_head    = cur_head + new_noct
        total_valid = total_valid + int(new_noct)
     end do

     ! Update nkill: slots consumed but not in any valid level.
     nkill          = int(n_all - total_valid)
     sim%m%noct_used = sim%m%tail(sim%r%nlevelmax)
     sim%m%ifree    = int(cur_head)   ! first free slot after valid octs

     ! --- Step 6: Hilbert sort per level ------------------------------------
     ! All passes for a given level are batched into one command buffer by
     ! mtl_hilbert_sort_level, replacing ndim*ilev × 4 commit/waits with 1.
     do ilev = int(ilevel + 1, c_int), int(sim%r%nlevelmax, c_int)
        if (sim%m%noct(ilev) <= 0) cycle
        call mtl_hilbert_sort_level( &
             int(sim%m%head(ilev), c_int), &
             int(sim%m%noct(ilev), c_int), &
             int(ndim * ilev,      c_int))
     end do

     ! --- Step 7: sort_gather/scatter grid/flag/hydro; blit ----------------
     call mtl_sort_gather_grid(head_child, n_all)
     call mtl_sort_scatter_grid(head_child, n_all)
     call mtl_sort_gather_flag(head_child, n_all)
     call mtl_sort_scatter_flag(head_child, n_all)
     call mtl_sort_gather_hydro(head_child, n_all)
     call mtl_blit_unew_to_uold(head_child, n_all)

     ! Reorder gravity variables if gravity is active
     if (sim%r%poisson) then
        call mtl_sort_gather_force(head_child, n_all, 1_c_int)
        call mtl_sort_scatter_force(head_child, n_all, 1_c_int)
        call mtl_sort_gather_force(head_child, n_all, 2_c_int)
        call mtl_sort_scatter_force(head_child, n_all, 2_c_int)
        call mtl_sort_gather_force(head_child, n_all, 3_c_int)
        call mtl_sort_scatter_force(head_child, n_all, 3_c_int)
        call mtl_sort_gather_phi(head_child, n_all, 0_c_int)
        call mtl_sort_scatter_phi(head_child, n_all, 0_c_int)
        call mtl_sort_gather_phi(head_child, n_all, 1_c_int)
        call mtl_sort_scatter_phi(head_child, n_all, 1_c_int)
     end if

     ! --- Step 8: update hash with new positions ----------------------------
     call mtl_update_hash_range(head_child, n_all, hash_size_l)

     ! --- Step 9: update father array per level (uses existing bridge) ------
     do ilev = int(ilevel + 1, c_int), int(sim%r%nlevelmax, c_int)
        if (sim%m%noct(ilev) <= 0) cycle
        call mtl_build_father( &
             int(sim%m%head(ilev),    c_int), &
             int(sim%m%noct(ilev),    c_int), &
             hash_size_l, &
             int(sim%m%ckey_max(ilev - 1), c_int), &
             int(sim%m%key_off(ilev - 1),  c_long))
     end do

  else
     ! No new or existing child octs — nothing to sort.
     nkill = 0
  end if

  ! --- Step 10: per-level cache rebuild -------------------------------------
  ! For each level above ilevel, walk all 27 neighbour directions.
  ! For each direction, find subgrids missing that nbor, create cache octs.
  ! ifree_cache_now is 1-based: starts at the tail of level ilevel's cache + 1.
  ifree_cache_now = int(sim%m%tail_cache(ilevel) + 1, c_int)
#ifdef MHD
  nsubgrid_l = mtl_nsubgrid()
  nsubgridtondim_l = nsubgrid_l**ndim
  subgridsize_l = (nsubgrid_l + 2)**ndim
#endif

  do ilev = int(ilevel + 1, c_int), int(sim%r%nlevelmax, c_int)
     if (sim%m%noct(ilev) <= 0) cycle

     new_head = ifree_cache_now
#ifdef MHD
     head_subgrid = int((sim%m%head(ilev) - 1) / nsubgridtondim_l + 1, c_int)
     num_subgrids = int(sim%m%noct(ilev) / nsubgridtondim_l, c_int)
     do ind = 1_c_int, subgridsize_l
#else
     do ind = 1_c_int, 27_c_int
#endif
        ! nbor_prefix + scan in one command buffer; returns missing-nbor count.
        cache_noct_lev = mtl_nbor_scan( &
#ifdef MHD
             head_subgrid, num_subgrids, &
#else
             int(sim%m%head(ilev), c_int), &
             int(sim%m%noct(ilev), c_int), &
#endif
             hash_size_l, ind)
        if (cache_noct_lev <= 0) cycle

        ! cache_swap + make_cache_octs + insert_hash in one command buffer.
        call mtl_cache_fill( &
#ifdef MHD
             head_subgrid, num_subgrids, &
#else
             int(sim%m%head(ilev), c_int), &
             int(sim%m%noct(ilev), c_int), &
#endif
             hash_size_l, ind, &
             int(sim%m%ngridmax, c_int), &
             ifree_cache_now, cache_noct_lev, &
             int(sim%r%interpol_var, c_int), int(sim%r%interpol_type, c_int), &
             real(sim%r%smallr, c_float))

        call mtl_advance_ifree_cache(cache_noct_lev)
        ifree_cache_now = ifree_cache_now + cache_noct_lev
     end do
     new_tail = ifree_cache_now - 1
     sim%m%head_cache(ilev) = int(new_head)
     sim%m%tail_cache(ilev) = int(new_tail)
     sim%m%noct_cache(ilev) = int(new_tail - new_head + 1)
  end do

  sim%m%ifree_cache = int(ifree_cache_now)

  ! Reset hash every coarse step so stale cache entries are purged.
  if (ilevel == sim%r%levelmin) then
     call mtl_reset_hash( &
          int(sim%m%ifree,       c_int), &
          int(sim%r%ngridmax,    c_int), &
          int(sim%m%ifree_cache, c_int), &
          hash_size_l)
  end if

end subroutine metal_refine

!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine metal_upload(sim, ilevel)
  ! Restriction: average 8 fine octs (ilevel+1) → coarse parent cells (ilevel).
  use ramses_commons, only: ramses_t
  use iso_c_binding
  implicit none
  type(ramses_t), intent(inout) :: sim
  integer,        intent(in)    :: ilevel

  integer(c_int) :: internal_energy

  if (ilevel >= sim%r%nlevelmax) return
  if (sim%m%noct(ilevel+1) <= 0) return

  internal_energy = int(merge(1, 0, sim%r%interpol_var == 1), c_int)

  call mtl_upload( &
       int(sim%m%head(ilevel+1), c_int), &
       int(sim%m%noct(ilevel+1), c_int), &
       internal_energy,                   &
       real(sim%r%gamma,        c_float), &
       real(sim%r%smallr,       c_float), &
       real(sim%r%smallc**2,    c_float))

end subroutine metal_upload

!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine metal_allocate_grav(sim)
  use amr_parameters, only: twotondim, ndim
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t), intent(inout) :: sim
  integer :: ngridmax_mg, ncachemax_mg, hash_size_mg

  ngridmax_mg  = sim%r%ngridmax / 7
  ncachemax_mg = max(sim%r%ncachemax / 7, 10000)
  sim%m_mg%ngridmax  = ngridmax_mg
  sim%m_mg%ncachemax = ncachemax_mg
  sim%m_mg%hash_size = 2 * (ngridmax_mg + ncachemax_mg)
  hash_size_mg = sim%m_mg%hash_size

  sim%m_mg%head = 1
  sim%m_mg%tail = 0
  sim%m_mg%noct = 0

  call mtl_alloc_grav( &
       int(sim%r%ngridmax,  c_int), &
       int(sim%r%ncachemax, c_int), &
       int(ngridmax_mg,     c_int), &
       int(ncachemax_mg,    c_int), &
       int(hash_size_mg,    c_int))

end subroutine metal_allocate_grav

!###########################################################
!###########################################################
subroutine metal_reset_rho(sim, ilevel)
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t), intent(inout) :: sim
  integer, intent(in) :: ilevel
  if (sim%m%noct(ilevel) <= 0) return
  call mtl_reset_rho( &
       int(sim%m%head(ilevel), c_int), &
       int(sim%m%noct(ilevel), c_int))
end subroutine metal_reset_rho

!###########################################################
!###########################################################
subroutine metal_multipole_leaf(sim, ilevel)
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t), intent(inout) :: sim
  integer, intent(in) :: ilevel
  real(c_float) :: scale
  scale = real(sim%r%boxlen / 2**ilevel, c_float)
  if (sim%m%noct(ilevel) <= 0) return
  call mtl_multipole_leaf( &
       int(sim%m%head(ilevel), c_int), &
       int(sim%m%noct(ilevel), c_int), &
       scale)
end subroutine metal_multipole_leaf

!###########################################################
!###########################################################
subroutine metal_multipole_upload(sim, ilevel)
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t), intent(inout) :: sim
  integer, intent(in) :: ilevel
  if (ilevel >= sim%r%nlevelmax) return
  if (sim%m%noct(ilevel+1) <= 0) return
  call mtl_multipole_upload( &
       int(sim%m%head(ilevel+1), c_int), &
       int(sim%m%noct(ilevel+1), c_int))
end subroutine metal_multipole_upload

!###########################################################
!###########################################################
subroutine metal_multipole_tot(sim, ilevel, tot4)
  use amr_parameters, only: ndim
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t), intent(inout) :: sim
  integer, intent(in) :: ilevel
  real(c_float), intent(out) :: tot4(ndim+1)
  if (sim%m%noct(ilevel) <= 0) then
     tot4 = 0.0_c_float; return
  end if
  call mtl_multipole_tot( &
       int(sim%m%head(ilevel), c_int), &
       int(sim%m%noct(ilevel), c_int), &
       tot4)
end subroutine metal_multipole_tot

!###########################################################
!###########################################################
subroutine metal_deposit_rho(sim, ilevel)
  use amr_parameters, only: ndim
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t), intent(inout) :: sim
  integer, intent(in) :: ilevel
  real(c_float) :: dx_f, vol_f, mref_f, msph_f, vcut_f
  real(kind=8) :: dx
  if (sim%m%noct(ilevel) <= 0) return
  dx     = sim%r%boxlen / 2**ilevel
  dx_f   = real(dx,                      c_float)
  vol_f  = real(dx**ndim,                c_float)
  mref_f = real(sim%r%m_refine(ilevel),  c_float)
  msph_f = real(sim%r%mass_sph,          c_float)
  vcut_f = real(sim%r%var_cut_refine,    c_float)
  call mtl_deposit_rho( &
       int(sim%m%head(ilevel),      c_int), &
       int(sim%m%noct(ilevel),      c_int), &
       dx_f, vol_f, mref_f, msph_f, vcut_f, &
       int(sim%r%ivar_refine,       c_int), &
       int(sim%m%ngridmax,          c_int))
end subroutine metal_deposit_rho

!###########################################################
!###########################################################
subroutine metal_cic_multipole2(sim, ilevel)
  use amr_parameters, only: ndim
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t), intent(inout) :: sim
  integer, intent(in) :: ilevel
  real(c_float) :: dx_f, vol_f, mref_f, msph_f, vcut_f
  real(c_float) :: tot4(ndim+1)
  real(kind=8) :: dx
  if (sim%m%noct(ilevel) <= 0) return
  dx     = sim%r%boxlen / 2**ilevel
  dx_f   = real(dx,                      c_float)
  vol_f  = real(dx**ndim,                c_float)
  mref_f = real(sim%r%m_refine(ilevel),  c_float)
  msph_f = real(sim%r%mass_sph,          c_float)
  vcut_f = real(sim%r%var_cut_refine,    c_float)
  call mtl_deposit_rho( &
       int(sim%m%head(ilevel),      c_int), &
       int(sim%m%noct(ilevel),      c_int), &
       dx_f, vol_f, mref_f, msph_f, vcut_f, &
       int(sim%r%ivar_refine,       c_int), &
       int(sim%m%ngridmax,          c_int))
  if (ilevel == sim%r%levelmin) then
     call mtl_multipole_tot( &
          int(sim%m%head(ilevel), c_int), &
          int(sim%m%noct(ilevel), c_int), &
          tot4)
     sim%g%multipole%q = sim%g%multipole%q + real(tot4, kind=8)
  end if
end subroutine metal_cic_multipole2

!###########################################################
!###########################################################
subroutine metal_upload_rho(sim, ilevel)
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t), intent(inout) :: sim
  integer, intent(in) :: ilevel
#ifdef GRAV
  if (sim%m%noct(ilevel) <= 0) return
  call mtl_upload_rho( &
       c_loc(sim%m%rho(1, sim%m%head(ilevel))), &
       int(sim%m%head(ilevel), c_int),           &
       int(sim%m%noct(ilevel), c_int))
#endif
end subroutine metal_upload_rho

!###########################################################
!###########################################################
subroutine metal_init_phi(sim, ilevel, icount)
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t), intent(inout) :: sim
  integer, intent(in) :: ilevel, icount
  real(kind=8) :: tfrac
  integer :: head_cache, noct_cache
  if (sim%m%noct(ilevel) <= 0) return
  ! Reset phi and f to zero for inner domain octs (initial guess)
  call mtl_reset_phi_fine( &
       int(sim%m%head(ilevel), c_int), &
       int(sim%m%noct(ilevel), c_int))
  if (ilevel > sim%r%levelmin) then
     if (sim%g%dtold(ilevel-1) > 0.0d0) then
        tfrac = sim%g%dtnew(ilevel) / sim%g%dtold(ilevel-1) * (icount - 1)
     else
        tfrac = 0.0d0
     end if
     ! Interpolate phi from coarser level for inner domain octs
     call mtl_init_phi( &
          int(sim%m%head(ilevel), c_int), &
          int(sim%m%noct(ilevel), c_int), &
          int(sim%m%ngridmax,     c_int), &
          real(tfrac, c_float))
     ! Interpolate phi from coarser level for cache (ghost) octs
     ! Note: we preserve f in ghost cells, so no reset here
     noct_cache = sim%m%noct_cache(ilevel)
     if (noct_cache > 0) then
        head_cache = sim%m%ngridmax + sim%m%head_cache(ilevel)
        call mtl_init_phi( &
             int(head_cache, c_int), &
             int(noct_cache, c_int), &
             int(sim%m%ngridmax, c_int), &
             real(tfrac, c_float))
     end if
  end if
end subroutine metal_init_phi

!###########################################################
!###########################################################
subroutine metal_save_phi_old(sim, ilevel)
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t), intent(inout) :: sim
  integer, intent(in) :: ilevel
  if (sim%m%noct(ilevel) <= 0) return
  call mtl_save_phi_old_fine( &
       int(sim%m%head(ilevel), c_int), &
       int(sim%m%noct(ilevel), c_int))
end subroutine metal_save_phi_old

!###########################################################
!###########################################################
subroutine metal_gradient_phi(sim, ilevel)
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t), intent(inout) :: sim
  integer, intent(in) :: ilevel
  real(c_float) :: dx
  dx = real(sim%r%boxlen / 2**ilevel, c_float)
  if (sim%m%noct(ilevel) <= 0) return
  call mtl_gradient_phi_fine( &
       int(sim%m%head(ilevel), c_int), &
       int(sim%m%noct(ilevel), c_int), &
       dx)
end subroutine metal_gradient_phi

!###########################################################
!###########################################################
subroutine metal_cmp_epot(sim, ilevel, epot)
  use amr_parameters, only: ndim
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t), intent(inout) :: sim
  integer, intent(in) :: ilevel
  real(kind=8), intent(out) :: epot
  real(c_float) :: epot_f
  real(kind=8) :: dx, fourpi, fact
  epot = 0.0d0
  if (sim%m%noct(ilevel) <= 0) return
  dx     = sim%r%boxlen / 2**ilevel
  fourpi = 4.0d0 * acos(-1.0d0)
  if (sim%r%cosmo) fourpi = 1.5d0 * sim%g%omega_m * sim%g%aexp
  fact = -dx**ndim / fourpi / 2.0d0
  call mtl_cmp_epot( &
       int(sim%m%head(ilevel), c_int), &
       int(sim%m%noct(ilevel), c_int), &
       epot_f)
  epot = fact * real(epot_f, kind=8)
end subroutine metal_cmp_epot

!###########################################################
!###########################################################
subroutine metal_cmp_rhomax(sim, ilevel, rhomax)
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t), intent(inout) :: sim
  integer, intent(in) :: ilevel
  real(kind=8), intent(out) :: rhomax
  real(c_float) :: rhomax_f
  rhomax = 0.0d0
  if (sim%m%noct(ilevel) <= 0) return
  call mtl_cmp_rhomax( &
       int(sim%m%head(ilevel), c_int), &
       int(sim%m%noct(ilevel), c_int), &
       rhomax_f)
  rhomax = real(rhomax_f, kind=8)
end subroutine metal_cmp_rhomax

!###########################################################
!###########################################################
subroutine metal_make_mask(sim, ilevel)
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t), intent(inout) :: sim
  integer, intent(in) :: ilevel
  if (sim%m%noct(ilevel) <= 0) return
  call mtl_reset_mask_fine( &
       int(sim%m%head(ilevel), c_int), &
       int(sim%m%noct(ilevel), c_int), &
       1.0_c_float)
end subroutine metal_make_mask

!###########################################################
!###########################################################
subroutine metal_make_rhs(sim, ilevel)
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t), intent(inout) :: sim
  integer, intent(in) :: ilevel
  real(c_float) :: fourpi_f, offset_f, oneoverdx2_f
  real(kind=8) :: fourpi, dx, oneoverdx2
  fourpi = 4.0d0 * acos(-1.0d0)
  if (sim%r%cosmo) fourpi = 1.5d0 * sim%g%omega_m * sim%g%aexp
  dx         = sim%r%boxlen / 2**ilevel
  oneoverdx2 = 1.0d0 / (dx * dx)
  fourpi_f     = real(fourpi,         c_float)
  offset_f     = real(sim%g%rho_tot,  c_float)
  if (any(.not. sim%r%periodic(1:3))) offset_f = 0.0_c_float
  oneoverdx2_f = real(oneoverdx2,     c_float)
  if (sim%m%noct(ilevel) <= 0) return
  call mtl_reset_rhs_fine( &
       int(sim%m%head(ilevel), c_int), &
       int(sim%m%noct(ilevel), c_int), &
       int(sim%m%ngridmax,     c_int), &
       fourpi_f, offset_f, oneoverdx2_f)
end subroutine metal_make_rhs

!###########################################################
!###########################################################
subroutine metal_clean_mg(sim)
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t), intent(inout) :: sim
  call mtl_clean_mg_hashes()
  sim%m_mg%head = 1
  sim%m_mg%tail = 0
  sim%m_mg%noct = 0
end subroutine metal_clean_mg

!###########################################################
!###########################################################
subroutine metal_build_mg(sim, ilevel, ifine)
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t), intent(inout) :: sim
  integer, intent(in) :: ilevel, ifine
  integer(c_int) :: head_idx, num_octs, head_father, new_noct

  if (ifine == ilevel) then
     head_idx    = int(sim%m%head(ilevel),    c_int)
     num_octs    = int(sim%m%noct(ilevel),    c_int)
     head_father = 1_c_int
  else
     head_idx    = int(sim%m_mg%head(ifine),  c_int)
     num_octs    = int(sim%m_mg%noct(ifine),  c_int)
     head_father = int(sim%m%noct(ilevel) + sim%m_mg%head(ifine), c_int)
  end if
  if (num_octs == 0) return

  sim%m_mg%head(ifine-1) = sim%m_mg%tail(ifine) + 1

  if (ifine == ilevel) then
     call mtl_build_mg_fine(int(ifine,c_int), int(ilevel,c_int), &
          head_idx, num_octs, head_father, &
          int(sim%m_mg%head(ifine-1), c_int), new_noct)
  else
     call mtl_build_mg_mg(int(ifine,c_int), int(ilevel,c_int), &
          head_idx, num_octs, head_father, &
          int(sim%m_mg%head(ifine-1), c_int), new_noct)
  end if

  sim%m_mg%tail(ifine-1) = sim%m_mg%head(ifine-1) + new_noct - 1
  sim%m_mg%noct(ifine-1) = new_noct


end subroutine metal_build_mg

!###########################################################
!###########################################################
subroutine metal_restrict_mask(sim, ilevel, ifine, allmasked)
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t), intent(inout) :: sim
  integer, intent(in) :: ilevel, ifine
  logical, intent(out) :: allmasked
  integer(c_int) :: head_idx, num_octs, head_father
  real(c_float) :: mask_max_f

  allmasked = .true.

  ! Zero volume fraction at coarse MG level
  call mtl_reset_mask_mg( &
       int(sim%m_mg%head(ifine-1), c_int), &
       int(sim%m_mg%noct(ifine-1), c_int), &
       0.0_c_float)

  ! Restrict from fine to coarse
  if (ifine == ilevel) then
     head_idx    = int(sim%m%head(ilevel),    c_int)
     num_octs    = int(sim%m%noct(ilevel),    c_int)
     head_father = 1_c_int
     if (num_octs > 0) call mtl_restrict_mask_fine(head_idx, head_father, num_octs)
  else
     head_idx    = int(sim%m_mg%head(ifine),  c_int)
     num_octs    = int(sim%m_mg%noct(ifine),  c_int)
     head_father = int(sim%m%noct(ilevel) + sim%m_mg%head(ifine), c_int)
     if (num_octs > 0) call mtl_restrict_mask_mg(head_idx, head_father, num_octs)
  end if

  ! Convert volume fraction → mask, get max
  if (sim%m_mg%noct(ifine-1) > 0) then
     call mtl_volume_to_mask_mg( &
          int(sim%m_mg%head(ifine-1), c_int), &
          int(sim%m_mg%noct(ifine-1), c_int), &
          mask_max_f)
     allmasked = (real(mask_max_f, kind=8) <= 0.0d0)
  end if

end subroutine metal_restrict_mask

!###########################################################
!###########################################################
subroutine metal_cmp_residual(sim, ilevel, ifine)
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t), intent(inout) :: sim
  integer, intent(in) :: ilevel, ifine
  integer(c_int) :: head_idx, num_octs, ngridmax_loc
  real(c_float) :: fourpi_f, offset_f, oneoverdx2_f
  real(kind=8) :: fourpi, dx, oneoverdx2

  fourpi = 4.0d0 * acos(-1.0d0)
  if (sim%r%cosmo) fourpi = 1.5d0 * sim%g%omega_m * sim%g%aexp
  dx         = sim%r%boxlen / 2**ifine
  oneoverdx2 = 1.0d0 / (dx*dx)
  fourpi_f     = real(fourpi,         c_float)
  offset_f     = real(sim%g%rho_tot,  c_float)
  if (any(.not. sim%r%periodic(1:3))) offset_f = 0.0_c_float
  oneoverdx2_f = real(oneoverdx2,     c_float)

  if (ifine == ilevel) then
     head_idx     = int(sim%m%head(ilevel),    c_int)
     num_octs     = int(sim%m%noct(ilevel),    c_int)
     ngridmax_loc = int(sim%m%ngridmax,        c_int)
     if (num_octs <= 0) return
     call mtl_cmp_residual_fine(head_idx, num_octs, ngridmax_loc, &
          fourpi_f, offset_f, oneoverdx2_f)
  else
     head_idx     = int(sim%m_mg%head(ifine),  c_int)
     num_octs     = int(sim%m_mg%noct(ifine),  c_int)
     ngridmax_loc = int(sim%m_mg%ngridmax,     c_int)
     if (num_octs <= 0) return
     call mtl_cmp_residual_mg(head_idx, num_octs, ngridmax_loc, &
          fourpi_f, offset_f, oneoverdx2_f)
  end if

end subroutine metal_cmp_residual

!###########################################################
!###########################################################
subroutine metal_gauss_seidel(sim, ilevel, ifine, safe, redstep)
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t), intent(inout) :: sim
  integer, intent(in) :: ilevel, ifine
  logical, intent(in) :: safe, redstep
  integer(c_int) :: head_idx, num_octs, ngridmax_loc, safe_i, redstep_i
  real(c_float) :: dx2_f

  safe_i    = merge(1_c_int, 0_c_int, safe)
  redstep_i = merge(1_c_int, 0_c_int, redstep)
  dx2_f     = real((sim%r%boxlen / 2**ifine)**2, c_float)

  if (ifine == ilevel) then
     head_idx     = int(sim%m%head(ilevel),    c_int)
     num_octs     = int(sim%m%noct(ilevel),    c_int)
     ngridmax_loc = int(sim%m%ngridmax,        c_int)
     if (num_octs <= 0) return
     call mtl_gauss_seidel_fine(head_idx, num_octs, ngridmax_loc, &
          dx2_f, safe_i, redstep_i)
  else
     head_idx     = int(sim%m_mg%head(ifine),  c_int)
     num_octs     = int(sim%m_mg%noct(ifine),  c_int)
     ngridmax_loc = int(sim%m_mg%ngridmax,     c_int)
     if (num_octs <= 0) return
     call mtl_gauss_seidel_mg(head_idx, num_octs, ngridmax_loc, &
          dx2_f, safe_i, redstep_i)
  end if

end subroutine metal_gauss_seidel

!###########################################################
!###########################################################
subroutine metal_reset_corr(sim, ilevel)
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t), intent(inout) :: sim
  integer, intent(in) :: ilevel
  if (sim%m_mg%noct(ilevel) <= 0) return
  call mtl_reset_phi_val_mg( &
       int(sim%m_mg%head(ilevel), c_int), &
       int(sim%m_mg%noct(ilevel), c_int), &
       0.0_c_float)
end subroutine metal_reset_corr

!###########################################################
!###########################################################
subroutine metal_restrict_residual(sim, ilevel, ifine)
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t), intent(inout) :: sim
  integer, intent(in) :: ilevel, ifine
  integer(c_int) :: head_idx, num_octs, head_father

  if (ifine == ilevel) then
     head_idx    = int(sim%m%head(ilevel),    c_int)
     num_octs    = int(sim%m%noct(ilevel),    c_int)
     head_father = 1_c_int
     if (num_octs <= 0) return
     call mtl_restrict_residual_fine(head_idx, head_father, num_octs)
  else
     head_idx    = int(sim%m_mg%head(ifine),  c_int)
     num_octs    = int(sim%m_mg%noct(ifine),  c_int)
     head_father = int(sim%m%noct(ilevel) + sim%m_mg%head(ifine), c_int)
     if (num_octs <= 0) return
     call mtl_restrict_residual_mg(head_idx, head_father, num_octs)
  end if

end subroutine metal_restrict_residual

!###########################################################
!###########################################################
subroutine metal_interpolate_correct(sim, ilevel, ifine)
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t), intent(inout) :: sim
  integer, intent(in) :: ilevel, ifine
  integer(c_int) :: head_idx, num_octs, head_father

  if (ifine == ilevel) then
     head_idx    = int(sim%m%head(ilevel),    c_int)
     num_octs    = int(sim%m%noct(ilevel),    c_int)
     head_father = 1_c_int
     if (num_octs <= 0) return
     call mtl_interpolate_correct_fine(head_idx, head_father, num_octs)
  else
     head_idx    = int(sim%m_mg%head(ifine),  c_int)
     num_octs    = int(sim%m_mg%noct(ifine),  c_int)
     head_father = int(sim%m%noct(ilevel) + sim%m_mg%head(ifine), c_int)
     if (num_octs <= 0) return
     call mtl_interpolate_correct_mg(head_idx, head_father, num_octs)
  end if

end subroutine metal_interpolate_correct

!###########################################################
!###########################################################
subroutine metal_residual_norm2(sim, ilevel, norm)
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t), intent(inout) :: sim
  integer, intent(in) :: ilevel
  real(kind=8), intent(out) :: norm
  real(c_float) :: norm_f
  real(kind=8) :: dx2
  norm = 0.0d0
  if (sim%m%noct(ilevel) <= 0) return
  dx2 = (sim%r%boxlen / 2**ilevel)**2
  call mtl_residual_norm_fine( &
       int(sim%m%head(ilevel), c_int), &
       int(sim%m%noct(ilevel), c_int), &
       norm_f)
  norm = dx2 * real(norm_f, kind=8)
end subroutine metal_residual_norm2

!###########################################################
!###########################################################
subroutine metal_rhs_norm2(sim, ilevel, norm)
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t), intent(inout) :: sim
  integer, intent(in) :: ilevel
  real(kind=8), intent(out) :: norm
  real(c_float) :: norm_f
  real(kind=8) :: dx2
  norm = 0.0d0
  if (sim%m%noct(ilevel) <= 0) return
  dx2 = (sim%r%boxlen / 2**ilevel)**2
  call mtl_rhs_norm_fine( &
       int(sim%m%head(ilevel), c_int), &
       int(sim%m%noct(ilevel), c_int), &
       norm_f)
  norm = dx2 * real(norm_f, kind=8)
end subroutine metal_rhs_norm2

!###########################################################
!###########################################################
subroutine metal_upload_cooling_table(c)
  use cooling_module, only: cooling_t
  implicit none
  type(cooling_t), intent(in), target :: c
  integer(c_int) :: n1, n2

  n1 = int(c%table%n1, c_int)
  n2 = int(c%table%n2, c_int)

  call mtl_upload_cooling_table(n1, n2, &
       c_loc(c%table%nH(1)), c_loc(c%table%T2(1)), &
       c_loc(c%table%cool(1,1)), c_loc(c%table%heat(1,1)), &
       c_loc(c%table%cool_com(1,1)), c_loc(c%table%heat_com(1,1)), &
       c_loc(c%table%metal(1,1)), c_loc(c%table%cool_prime(1,1)), &
       c_loc(c%table%heat_prime(1,1)), c_loc(c%table%cool_com_prime(1,1)), &
       c_loc(c%table%heat_com_prime(1,1)), c_loc(c%table%metal_prime(1,1)))

end subroutine metal_upload_cooling_table

!###########################################################
!###########################################################
subroutine metal_cooling(sim, ilevel)
  use ramses_commons, only: ramses_t
  use constants, only: rhoc, mH
  implicit none
  type(ramses_t), intent(inout) :: sim
  integer, intent(in) :: ilevel

  integer(c_int) :: head_idx, num_octs
  real(c_float) :: gamma, smallr, smallc2
  real(kind=8) :: dtcool, scale_T2, scale_nH, scale_l, scale_d, scale_t, scale_v
  real(kind=8) :: nH_eos, nCOM
  integer(c_int) :: cooling, metal, imetal, self_shielding, eos_type, isothermal

  gamma = real(sim%r%gamma, c_float)
  smallr = real(sim%r%smallr, c_float)
  smallc2 = real(sim%r%smallc**2, c_float)

  call units(sim%r,sim%g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Reference density for EOS floor
  nH_eos = sim%r%eos_nH
  if (sim%r%cosmo) then
     nCOM = 200.0d0*sim%g%omega_b*rhoc*(sim%g%h0/100)**2/sim%g%aexp**3*sim%cool%X/mH
     nH_eos = max(nCOM, nH_eos)
  end if

  dtcool = sim%g%dtnew(ilevel) * scale_t

  head_idx = int(sim%m%head(ilevel), c_int)
  num_octs = int(sim%m%noct(ilevel), c_int)

  if (num_octs <= 0) return

  cooling = 0_c_int
  if (sim%r%cooling) cooling = 1_c_int

  metal = 0_c_int
  if (sim%r%metal) metal = 1_c_int

  self_shielding = 0_c_int
  if (sim%r%self_shielding) self_shielding = 1_c_int

  isothermal = 0_c_int
  if (sim%r%isothermal) isothermal = 1_c_int

  imetal = int(sim%r%imetal, c_int)
  eos_type = int(sim%r%eos_type, c_int)

  call mtl_cooling(head_idx, num_octs, gamma, smallr, smallc2, &
       dtcool, eos_type, sim%r%eos_T2, nH_eos, sim%r%eos_index, &
       scale_T2, scale_nH, cooling, metal, imetal, sim%r%z_ave, &
       self_shielding, sim%cool%X, sim%r%T2max, isothermal)

end subroutine metal_cooling

!###########################################################
!###########################################################
subroutine metal_sync_hydro(sim, ilevel, dt)
  use ramses_commons, only: ramses_t
  use iso_c_binding
  implicit none
  type(ramses_t), intent(inout) :: sim
  integer,        intent(in)    :: ilevel
  real(kind=8),   intent(in)    :: dt

  integer(c_int) :: head_idx, num_octs
  real(c_float)  :: gamma, smallr, smallc2
  real(c_float)  :: constant_gravity(3)

  if (sim%m%noct(ilevel) <= 0) return

  gamma = real(sim%r%gamma, c_float)
  smallr = real(sim%r%smallr, c_float)
  smallc2 = real(sim%r%smallc**2, c_float)
  constant_gravity = real(sim%r%constant_gravity, c_float)

  head_idx = int(sim%m%head(ilevel), c_int)
  num_octs = int(sim%m%noct(ilevel), c_int)

  call mtl_sync_hydro( &
       head_idx, num_octs, &
       gamma, smallr, smallc2, &
       real(dt, c_float), constant_gravity)
end subroutine metal_sync_hydro

!###########################################################
!###########################################################
subroutine metal_grav_hydro(sim, ilevel)
  use ramses_commons, only: ramses_t
  use iso_c_binding
  implicit none
  type(ramses_t), intent(inout) :: sim
  integer,        intent(in)    :: ilevel

  integer(c_int) :: head_idx, num_octs
  real(c_float)  :: gamma, smallr, smallc2, dt
  real(c_float)  :: constant_gravity(3)

  if (sim%m%noct(ilevel) <= 0) return

  dt = real(0.5d0 * sim%g%dtnew(ilevel), c_float)
  gamma = real(sim%r%gamma, c_float)
  smallr = real(sim%r%smallr, c_float)
  smallc2 = real(sim%r%smallc**2, c_float)
  constant_gravity = real(sim%r%constant_gravity, c_float)

  head_idx = int(sim%m%head(ilevel), c_int)
  num_octs = int(sim%m%noct(ilevel), c_int)

  call mtl_grav_hydro( &
       head_idx, num_octs, &
       gamma, smallr, smallc2, &
       dt, constant_gravity)
end subroutine metal_grav_hydro

!###########################################################
!###########################################################
subroutine metal_allocate_part(sim)
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t), intent(inout) :: sim
  call mtl_alloc_part(int(sim%r%npartmax, c_int))
end subroutine metal_allocate_part

!###########################################################
!###########################################################
subroutine metal_upload_part(sim)
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t), intent(inout), target :: sim
  type(c_ptr) :: xp_ptr, vp_ptr, mp_ptr, levelp_ptr, sortp_ptr, idp_ptr
  xp_ptr = c_null_ptr
  vp_ptr = c_null_ptr
  mp_ptr = c_null_ptr
  levelp_ptr = c_null_ptr
  sortp_ptr = c_null_ptr
  idp_ptr = c_null_ptr
  if (allocated(sim%p%xp)) xp_ptr = c_loc(sim%p%xp(1,1))
  if (allocated(sim%p%vp)) vp_ptr = c_loc(sim%p%vp(1,1))
  if (allocated(sim%p%mp)) mp_ptr = c_loc(sim%p%mp(1))
  if (allocated(sim%p%levelp)) levelp_ptr = c_loc(sim%p%levelp(1))
  if (allocated(sim%p%sortp)) sortp_ptr = c_loc(sim%p%sortp(1))
  if (allocated(sim%p%idp)) idp_ptr = c_loc(sim%p%idp(1))
  call mtl_upload_part(xp_ptr, vp_ptr, mp_ptr, levelp_ptr, sortp_ptr, idp_ptr, int(sim%p%npart, c_int))
end subroutine metal_upload_part

!###########################################################
!###########################################################
subroutine metal_download_part(sim)
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t), intent(inout), target :: sim
  type(c_ptr) :: xp_ptr, vp_ptr, mp_ptr, levelp_ptr, sortp_ptr, idp_ptr
  xp_ptr = c_null_ptr
  vp_ptr = c_null_ptr
  mp_ptr = c_null_ptr
  levelp_ptr = c_null_ptr
  sortp_ptr = c_null_ptr
  idp_ptr = c_null_ptr
  if (allocated(sim%p%xp)) xp_ptr = c_loc(sim%p%xp(1,1))
  if (allocated(sim%p%vp)) vp_ptr = c_loc(sim%p%vp(1,1))
  if (allocated(sim%p%mp)) mp_ptr = c_loc(sim%p%mp(1))
  if (allocated(sim%p%levelp)) levelp_ptr = c_loc(sim%p%levelp(1))
  if (allocated(sim%p%sortp)) sortp_ptr = c_loc(sim%p%sortp(1))
  if (allocated(sim%p%idp)) idp_ptr = c_loc(sim%p%idp(1))
  call mtl_download_part(xp_ptr, vp_ptr, mp_ptr, levelp_ptr, sortp_ptr, idp_ptr, int(sim%p%npart, c_int))
end subroutine metal_download_part

!###########################################################
!###########################################################
subroutine metal_kick_drift_part(sim, ilevel, action_part)
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t), intent(inout), target :: sim
  integer, intent(in) :: ilevel
  integer, intent(in) :: action_part
  integer :: head_idx, tail_idx, num_parts
  real(kind=8) :: dx_loc
  real(c_float), target :: box_size_f(3), dtnew_f(size(sim%g%dtnew)), dtold_f(size(sim%g%dtold))
  integer(c_int), target :: periodic_i(3)
  box_size_f = real(sim%r%boxlen, c_float)
  dtnew_f = real(sim%g%dtnew, c_float)
  dtold_f = real(sim%g%dtold, c_float)
  periodic_i(1:3) = merge(1, 0, sim%r%periodic(1:3))
  dx_loc = sim%r%boxlen/2.0d0**ilevel
  head_idx = sim%p%headp(ilevel)
  tail_idx = sim%p%tailp(ilevel)
  num_parts = tail_idx - head_idx + 1
  if (num_parts > 0) then
     call mtl_kick_drift_part( &
          int(action_part, c_int), &
          int(ilevel, c_int), &
          int(head_idx, c_int), &
          int(num_parts, c_int), &
          real(sim%m%skip(1), c_float), &
          real(sim%m%skip(2), c_float), &
          real(sim%m%skip(3), c_float), &
          real(dx_loc, c_float), &
          c_loc(box_size_f(1)), &
          c_loc(periodic_i(1)), &
          c_loc(dtnew_f(1)), &
          c_loc(dtold_f(1)))
  endif
end subroutine metal_kick_drift_part

!###########################################################
!###########################################################
subroutine metal_newdt_part(sim, ilevel, vmax, ekin)
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t), intent(inout), target :: sim
  integer, intent(in) :: ilevel
  real(kind=8), intent(inout) :: vmax, ekin
  integer :: head_idx, tail_idx, num_parts
  real(c_float) :: vmax_f, ekin_f
  head_idx = sim%p%headp(ilevel)
  tail_idx = sim%p%tailp(ilevel)
  num_parts = tail_idx - head_idx + 1
  if (num_parts > 0) then
     call mtl_newdt_part( &
          int(head_idx, c_int), &
          int(num_parts, c_int), &
          vmax_f, ekin_f)
     vmax = max(vmax, real(vmax_f, kind=8))
     ekin = ekin + real(ekin_f, kind=8)
  endif
end subroutine metal_newdt_part

!###########################################################
!###########################################################
subroutine metal_split_part(sim, ilevel)
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t), intent(inout), target :: sim
  integer, intent(in) :: ilevel
  integer :: head_idx, num_parts, n_fine, n_coarse, ilev
  real(kind=8) :: dx_loc
  integer(c_int) :: n_fine_c
  head_idx = sim%p%headp(ilevel)
  num_parts = sim%p%tailp(sim%r%nlevelmax) - head_idx + 1
  if (num_parts <= 0) then
     sim%p%tailp(ilevel) = sim%p%headp(ilevel) - 1
     do ilev = ilevel + 1, sim%r%nlevelmax
        sim%p%headp(ilev) = sim%p%tailp(ilevel) + 1
        sim%p%tailp(ilev) = sim%p%npart
     end do
     return
  endif
  dx_loc = sim%r%boxlen/2.0d0**ilevel
  call mtl_split_part( &
       int(head_idx, c_int), &
       int(num_parts, c_int), &
       int(ilevel, c_int), &
       real(sim%m%skip(1), c_float), &
       real(sim%m%skip(2), c_float), &
       real(sim%m%skip(3), c_float), &
       real(dx_loc, c_float), &
       n_fine_c)
  n_fine = int(n_fine_c)
  n_coarse = num_parts - n_fine
  sim%p%tailp(ilevel) = sim%p%headp(ilevel) + n_coarse - 1
  do ilev = ilevel + 1, sim%r%nlevelmax
     sim%p%headp(ilev) = sim%p%tailp(ilevel) + 1
     sim%p%tailp(ilev) = sim%p%npart
  end do
end subroutine metal_split_part

!###########################################################
!###########################################################
subroutine metal_sort_part(sim, ilevel)
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t), intent(inout), target :: sim
  integer, intent(in) :: ilevel
  integer :: head_idx, tail_idx, num_parts
  real(kind=8) :: dx, dx_inv, shift
  real(c_float), target :: skip_f(3)
  head_idx = sim%p%headp(ilevel)
  tail_idx = sim%p%tailp(sim%r%nlevelmax)
  num_parts = tail_idx - head_idx + 1
  if (num_parts <= 0) return
  dx = sim%r%boxlen/2.0d0**ilevel
  dx_inv = 1.0d0 / dx
  if (sim%r%part_dep_algo >= 2) then
     shift = 0.5d0
  else
     shift = 0.0d0
  endif
  skip_f = real(sim%m%skip, c_float)
  call mtl_sort_part( &
       int(head_idx, c_int), &
       int(num_parts, c_int), &
       int(ilevel, c_int), &
       real(shift, c_float), &
       real(dx_inv, c_float), &
       c_loc(skip_f(1)))
end subroutine metal_sort_part

!###########################################################
!###########################################################
subroutine metal_cic_part_medium(sim, ilevel, rtype)
  use ramses_commons, only: ramses_t
  use amr_parameters, only: ndim, twotondim
  implicit none
  type(ramses_t), intent(inout) :: sim
  integer, intent(in) :: ilevel
  integer, intent(in) :: rtype
  integer :: head_idx, tail_idx, num_parts, idim
  real(kind=8) :: dx_loc, vol_loc
  real(c_float) :: q_out(4)
  if (rtype /= 0 .and. rtype /= 1) return

  ! At the coarsest level accumulate monopole + dipole on the GPU.
  if (ilevel == sim%r%levelmin .and. sim%p%npart > 0) then
     call mtl_multipole_q_part( &
          int(1, c_int), &
          int(sim%p%npart, c_int), &
          int(sim%r%npartmax, c_long), &
          q_out)
     do idim = 1, ndim + 1
        sim%g%multipole%q(idim) = sim%g%multipole%q(idim) + real(q_out(idim), kind=8)
     end do
  end if

  head_idx = sim%p%headp(ilevel)
  tail_idx = sim%p%tailp(sim%r%nlevelmax)
  num_parts = tail_idx - head_idx + 1
  if (num_parts <= 0) return
  dx_loc = sim%r%boxlen/2.0d0**ilevel
  vol_loc = dx_loc**ndim
  call mtl_cic_part_medium( &
       int(head_idx, c_int), &
       int(num_parts, c_int), &
       real(sim%m%skip(1), c_float), &
       real(sim%m%skip(2), c_float), &
       real(sim%m%skip(3), c_float), &
       real(dx_loc, c_float), &
       real(vol_loc, c_float), &
       real(1.0d0, c_float), &
       int(0, c_int), &
       real(sim%r%m_refine(ilevel), c_float), &
       real(sim%r%mass_cut_refine, c_float), &
       int(ilevel, c_int))
end subroutine metal_cic_part_medium

end module metal_runner
