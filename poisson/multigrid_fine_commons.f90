module multigrid_fine_commons

#ifdef GRAV
  use multigrid_fine_coarse, only: r_cmp_residual_mg, r_cmp_residual_norm2, r_gauss_seidel_mg,&
        r_interpolate_and_correct, r_reset_correction, r_restrict_mask, r_restrict_residual, r_set_scan_flag,&
        double_level_t, level_count_t, gs_step_t
#endif

#ifdef _CUDA
  use gpu_runner, only: gpu_make_mask, gpu_make_rhs, gpu_build_mg, gpu_clean_mg, gpu_rhs_norm2
#endif
#ifdef _METAL
  use metal_runner, only: metal_make_mask, metal_make_rhs, metal_build_mg, metal_clean_mg, metal_rhs_norm2
#endif

contains

! ------------------------------------------------------------------------
! Multigrid Poisson solver for refined AMR levels
! ------------------------------------------------------------------------
! This file contains all generic fine multigrid routines, such as
!   * multigrid iterations @ MG levels
!   * MG workspace building
!
! Used variables:
!     -----------------------------------------------------------------
!     potential            phi
!     physical RHS         rho
!     residual             f(:,1)
!     BC-modified RHS      f(:,2)
!     mask                 f(:,3)
!     scan flag            flag2
!
! ------------------------------------------------------------------------

#ifdef GRAV

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Main multigrid routine, called by amr_step
! ------------------------------------------------------------------------

subroutine multigrid(pst,ilevel,icount)
  use amr_parameters, only: twotondim
  use poisson_parameters, only: ngs_fine, ngs_coarse, ncycles_coarse_safe
  use ramses_commons, only: pst_t
  use phi_fine_cg_module, only: r_make_initial_phi
  implicit none
  type(pst_t)::pst
  integer,intent(in) :: ilevel,icount

  integer,parameter :: MAXITER  = 20
  real(kind=8),parameter :: SAFE_FACTOR = 0.5

  integer :: igrid, ifine, i, iter, allmasked
  integer,dimension(1:4) :: output_array
  real(kind=8) :: res_norm2, i_res_norm2, rhs_norm2
  real(kind=8) :: err, last_err
  real(kind=8) :: i_res_norm2_tot, res_norm2_tot
  type(double_level_t)::double_level
  type(level_count_t)::level_count
  type(gs_step_t)::gs_step

  if(pst%s%r%gravity_type>0)return
  if(pst%s%m%noct_tot(ilevel)==0)return

  if(pst%s%r%verbose) print '(A,I2)','Entering multigrid at level ',ilevel

  ! ---------------------------------------------------------------------
  ! Prepare first guess, mask and BCs at finest level
  ! ---------------------------------------------------------------------
  level_count%ilevel=ilevel
  level_count%icount=icount
  call r_make_initial_phi(pst,level_count,storage_size(level_count)/32) ! Initial guess
  call r_make_mask(pst,ilevel,1) ! Fill the fine level mask
  call r_make_bc_rhs(pst,level_count,storage_size(level_count)/32) ! Fill BC-modified RHS
  call r_cmp_rhs_norm2(pst,ilevel,1,rhs_norm2,2) ! Compute RHS norm

  if(pst%s%r%verbose) print '(A)','Initial guess done '
  if(pst%s%r%verbose) print '(A,1pE10.3)','RHS norm2=',rhs_norm2

  ! ---------------------------------------------------------------------
  ! Initialize Domain Decomposition and Hash Table for Multigrid
  ! ---------------------------------------------------------------------
  call r_init_mg(pst,ilevel,1)

  if(pst%s%r%verbose) print '(A)','Multigrid init done '

  ! ---------------------------------------------------------------------
  ! Build Multigrid hierarchy in memory
  ! ---------------------------------------------------------------------
  double_level%ilevel=ilevel
  do ifine=ilevel,pst%s%r%bound_levelmin+1,-1
     double_level%ifine=ifine
     if(pst%s%r%verbose) print '(A,I2)','Build MG ',ifine
     call r_build_mg(pst,double_level,storage_size(double_level)/32)
  end do

  if(pst%s%r%verbose) print '(A)','Multigrid hierarchy done '

  ! ---------------------------------------------------------------------
  ! Restrict mask up
  ! ---------------------------------------------------------------------
  pst%s%g%levelmin_mg=pst%s%r%bound_levelmin
  double_level%ilevel=ilevel
  do ifine=ilevel,pst%s%r%bound_levelmin+1,-1
     ! Restrict and communicate mask
     double_level%ifine=ifine
     call r_restrict_mask(pst,double_level,storage_size(double_level)/32,allmasked,1)
     if(allmasked==1) then ! Coarser level is fully masked: stop here
        pst%s%g%levelmin_mg=ifine
        exit
     end if
  end do

  if(pst%s%r%verbose) print '(A)','Restrict mask up done '

  ! ---------------------------------------------------------------------
  ! Set scan flag (for optimisation)
  ! ---------------------------------------------------------------------
  double_level%ilevel=ilevel
  do ifine=ilevel,pst%s%g%levelmin_mg,-1
     double_level%ifine=ifine
     call r_set_scan_flag(pst,double_level,storage_size(double_level)/32)
  end do

  if(pst%s%r%verbose) print '(A)','Mask and scan done '

  ! ---------------------------------------------------------------------
  ! Initiate solve at fine level
  ! ---------------------------------------------------------------------

  ! Switch back to normal mode after mg_safe_mode_reset solves in safe mode
  if(pst%s%g%safe_mode(ilevel).and.pst%s%r%mg_safe_mode_reset>0)then
     pst%s%g%safe_mode_countdown(ilevel)=pst%s%g%safe_mode_countdown(ilevel)-1
     if(pst%s%g%safe_mode_countdown(ilevel)<=0)then
        if(pst%s%r%verbose)print *,'CAUTION: Switching back to normal MG mode for level ',ilevel
        pst%s%g%safe_mode(ilevel)=.false.
     end if
  end if

  iter = 0
  err = 1.0d0
  main_iteration_loop: do

     iter=iter+1

     gs_step%ilevel=ilevel
     gs_step%ifine=ilevel
     gs_step%safe=pst%s%g%safe_mode(ilevel)

     ! Pre-smoothing
     do i=1,ngs_fine
        gs_step%redstep=.true.   ! Red step
        call r_gauss_seidel_mg(pst,gs_step,storage_size(gs_step)/32)
        gs_step%redstep=.false.  ! Black step
        call r_gauss_seidel_mg(pst,gs_step,storage_size(gs_step)/32)
     end do

     ! Compute new residual
     double_level%ilevel=ilevel
     double_level%ifine=ilevel
     call r_cmp_residual_mg(pst,double_level,storage_size(double_level)/32)

!!$     ! Compute initial residual norm
!!$     if(iter==1) then
!!$        call r_cmp_residual_norm2(pst,ilevel,1,i_res_norm2,2)
!!$     end if

     ! Recurse only if a coarser multigrid level exists
     if(ilevel>pst%s%g%levelmin_mg) then

        ! Restrict residual to coarser level
        call r_restrict_residual(pst,double_level,storage_size(double_level)/32)

        ! Reset correction from upper level before solve
        call r_reset_correction(pst,ilevel-1,1)

        ! Multigrid-solve the upper level
        call recursive_multigrid(pst,ilevel,ilevel-1,pst%s%g%safe_mode(ilevel))

        ! Interpolate coarse solution and correct fine solution
        call r_interpolate_and_correct(pst,double_level,storage_size(double_level)/32)

     end if

     ! Post-smoothing
     do i=1,ngs_fine
        gs_step%redstep=.true.   ! Red step
        call r_gauss_seidel_mg(pst,gs_step,storage_size(gs_step)/32)
        gs_step%redstep=.false.  ! Black step
        call r_gauss_seidel_mg(pst,gs_step,storage_size(gs_step)/32)
     end do

     ! Update fine residual
     double_level%ilevel=ilevel
     double_level%ifine=ilevel
     call r_cmp_residual_mg(pst,double_level,storage_size(double_level)/32)

     ! Compute residual norm
     call r_cmp_residual_norm2(pst,ilevel,1,res_norm2,2)

     last_err = err
!     err = sqrt(res_norm2/(i_res_norm2+1d-20*pst%s%g%rho_tot**2))
     err = sqrt(res_norm2/(rhs_norm2+1d-20*pst%s%g%rho_tot**2))

     ! Verbosity
     if(pst%s%r%verbose) print '(A,I0,A,1pE10.3)','   ==> Step ',iter,' Error=',err

     ! Converged?
     if(err<pst%s%r%epsilon .or. iter == pst%s%r%nvcycle .or. iter>=MAXITER) exit

     ! Not converged, check error and possibly enable safe mode for the level
     if(err > last_err*SAFE_FACTOR .and. (.not. pst%s%g%safe_mode(ilevel))) then
        if(pst%s%r%verbose)print *,'CAUTION: Switching to safe MG mode for level ',ilevel
        pst%s%g%safe_mode(ilevel) = .true.
        pst%s%g%safe_mode_countdown(ilevel) = pst%s%r%mg_safe_mode_reset
     end if

  end do main_iteration_loop

  print '(A,I0,A,I0,A,1pE10.3)','   ==> Level ',ilevel,' Step ',iter,' Error=',err
  if(iter==MAXITER) print *,'WARN: Fine multigrid Poisson failed to converge...'

  ! ---------------------------------------------------------------------
  ! Cleanup MG levels after solve complete
  ! ---------------------------------------------------------------------
  call r_cleanup_mg(pst)

end subroutine multigrid

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Recursive multigrid routine for coarse MG levels
! ------------------------------------------------------------------------

recursive subroutine recursive_multigrid(pst,ilevel,ifinelevel,safe)
  use amr_parameters, only: twotondim
  use poisson_parameters, only: ngs_fine, ngs_coarse, ncycles_coarse_safe
  use ramses_commons, only: pst_t
  implicit none
  type(pst_t)::pst
  integer,intent(in) :: ilevel
  integer,intent(in) :: ifinelevel
  logical,intent(in) :: safe

  integer :: i, igrid, icycle, ncycle
  type(double_level_t)::double_level
  type(gs_step_t)::gs_step

  ! Set parameter array
  gs_step%ilevel=ilevel
  gs_step%ifine=ifinelevel
  gs_step%safe=safe

  if(ifinelevel<=pst%s%g%levelmin_mg) then

     ! Solve 'directly' :
     do i=1,2*ngs_coarse
        gs_step%redstep=.true.   ! Red step
        call r_gauss_seidel_mg(pst,gs_step,storage_size(gs_step)/32)
        gs_step%redstep=.false.  ! Black step
        call r_gauss_seidel_mg(pst,gs_step,storage_size(gs_step)/32)
     end do

     return
  end if

  if(safe) then
     ncycle=ncycles_coarse_safe
  else
     ncycle=1
  endif

  do icycle=1,ncycle

     ! Pre-smoothing
     do i=1,ngs_coarse
        gs_step%redstep=.true.   ! Red step
        call r_gauss_seidel_mg(pst,gs_step,storage_size(gs_step)/32)
        gs_step%redstep=.false.  ! Black step
        call r_gauss_seidel_mg(pst,gs_step,storage_size(gs_step)/32)
     end do

     ! Compute residual and restrict into upper level RHS
     double_level%ilevel=ilevel
     double_level%ifine=ifinelevel
     call r_cmp_residual_mg(pst,double_level,storage_size(double_level)/32)

     ! Restrict residual to coarser level
     call r_restrict_residual(pst,double_level,storage_size(double_level)/32)

     ! Reset correction from upper level before solve
     call r_reset_correction(pst,ifinelevel-1,1)

     ! Multigrid-solve the upper level
     call recursive_multigrid(pst,ilevel,ifinelevel-1, safe)

     ! Interpolate coarse solution and correct back into fine solution
     call r_interpolate_and_correct(pst,double_level,storage_size(double_level)/32)

     ! Post-smoothing
     do i=1,ngs_coarse
        gs_step%redstep=.true.   ! Red step
        call r_gauss_seidel_mg(pst,gs_step,storage_size(gs_step)/32)
        gs_step%redstep=.false.  ! Black step
        call r_gauss_seidel_mg(pst,gs_step,storage_size(gs_step)/32)
     end do

  end do

end subroutine recursive_multigrid

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Multigrid workspace initialisation
! ------------------------------------------------------------------------

recursive subroutine r_init_mg(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::ilevel

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_INIT_MG,pst%iUpper+1,input_size,0,ilevel)
     call r_init_mg(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call init_mg(pst%s%r,pst%s%m,pst%s%m_mg,ilevel)
  endif

end subroutine r_init_mg

subroutine init_mg(r,m,m_mg,ilevel)
  use amr_parameters, only: nhilbert
  use amr_commons, only: run_t,mesh_t
  use hilbert
  implicit none
  type(run_t)::r
  type(mesh_t)::m,m_mg
  integer::ilevel

  integer::ilev,idom

  ! Compute multigrid Hilbert key tick marks
  call m_mg%domain(ilevel)%copy(m%domain(ilevel))
  do ilev=ilevel-1,1,-1
     call m_mg%domain(ilev)%copy(m_mg%domain(ilev+1))
     do idom=0,m_mg%domain(ilev)%ncpu
        m_mg%domain(ilev)%b(1:nhilbert,idom)=coarsen_key(m_mg%domain(ilev+1)%b(1:nhilbert,idom),ilev)
     end do
  end do

end subroutine init_mg

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ---------------------------------------------------------------------
! Coarse grid MG activation for local grids
! ---------------------------------------------------------------------

recursive subroutine r_build_mg(pst,input,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  type(double_level_t)::input

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_BUILD_MG,pst%iUpper+1,input_size,0,input)
     call r_build_mg(pst%pLower,input,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
#ifdef _CUDA
     call gpu_build_mg(pst%s,input%ilevel,input%ifine)
#elif defined(_METAL)
     call metal_build_mg(pst%s,input%ilevel,input%ifine)
#else
     if(input%ifine==input%ilevel)then
        call build_mg(pst%s,pst%s%m,input%ifine)
     else
        call build_mg(pst%s,pst%s%m_mg,input%ifine)
     end if
#endif
  endif

end subroutine r_build_mg

subroutine build_mg(s,m,ifinelevel)
  use mdl_module
  use amr_parameters, only: nhilbert, ndim, twotondim
  use ramses_commons, only: ramses_t
  use amr_commons, only: mesh_t
  use cache_commons
  use cache
  use nbors_utils
  use hilbert
  use hash
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  type(ramses_t)::s
  type(mesh_t)::m
  integer,intent(in)::ifinelevel

  integer::icoarselevel,igrid,inbor,idim,ichild,grid_cpu,ind,ifather
  integer(kind=8),dimension(0:ndim)::hash_key,hash_father,hash_nbor
  integer(kind=4),dimension(1:ndim)::cart_key
  integer,dimension(1:3,1:8),save::shift_oct=reshape(&
       & (/-1,-1,-1,+1,-1,-1,-1,+1,-1,+1,+1,-1,&
       &   -1,-1,+1,+1,-1,+1,-1,+1,+1,+1,+1,+1/),(/3,8/))
  integer(kind=8),dimension(1:nhilbert)::hk
  integer(kind=8),dimension(1:ndim)::ix
  logical::in_rank,in_domain
  type(msg_small_realdp)::dummy_small_realdp

  associate(r=>s%r,g=>s%g,m_mg=>s%m_mg,mdl=>s%mdl)

  icoarselevel=ifinelevel-1
  m_mg%ifree=m_mg%noct_used+1
  m_mg%head(icoarselevel)=m_mg%ifree
  hash_father(0)=icoarselevel

  call open_cache(mdl, m_mg, pack_size=storage_size(dummy_small_realdp)/32, &
       flush=pack_flush_build_mg, combine=unpack_flush_build_mg)

  ! Loop over fine grids
  do igrid=m%head(ifinelevel),m%tail(ifinelevel)

     hash_key(1:ndim)=m%grid(igrid)%ckey(1:ndim)

     ! Gather twotondim neighboring father grids
     do inbor=1,twotondim

#ifndef WITHOUTMPI
        ! If counter is good, check on incoming messages and perform actions
        if(mdl%mail_counter==32)then
           call check_mail(mdl,MPI_REQUEST_NULL)
           mdl%mail_counter=0
        endif
        mdl%mail_counter=mdl%mail_counter+1
#endif
        ! Get neighboring grid
        hash_nbor(1:ndim)=hash_key(1:ndim)+shift_oct(1:ndim,inbor)

        ! Periodic boundary conditions
        do idim=1,ndim
           if(r%periodic(idim))then
              if(hash_nbor(idim)< m%box_ckey_min(idim,ifinelevel))hash_nbor(idim)=m%box_ckey_max(idim,ifinelevel)-1
              if(hash_nbor(idim)>=m%box_ckey_max(idim,ifinelevel))hash_nbor(idim)=m%box_ckey_min(idim,ifinelevel)
           endif
        enddo
        hash_father(1:ndim)=hash_nbor(1:ndim)/2

        in_domain = .true.
        do idim = 1, ndim
           in_domain = in_domain .and. hash_father(idim) .ge. m_mg%box_ckey_min(idim,icoarselevel) &
                &                .and. hash_father(idim) .lt. m_mg%box_ckey_max(idim,icoarselevel)
        end do

        if(in_domain)then

           ! Access hash table
           ifather=hash_getp(m_mg%grid_dict,hash_father)

           ! If grid does not exist, create it in memory
           if(ifather<=0)then

              ! Compute Cartesian keys of new oct
              cart_key(1:ndim)=int(hash_father(1:ndim),kind=4)

              ! Compute Hilbert keys of new octs
              ix(1:ndim)=cart_key(1:ndim)
              hk(1:nhilbert)=hilbert_key(ix,icoarselevel-1)

              ! Check if grid sits inside processor boundaries
              in_rank = ge_keys(hk,m_mg%domain(icoarselevel)%b(1:nhilbert,mdl_self(mdl)-1)).and. &
                   &    gt_keys(m_mg%domain(icoarselevel)%b(1:nhilbert,mdl_self(mdl)),hk)

              if(in_rank)then

                 ! Set grid index to a virtual grid in local main memory
                 ichild=m_mg%ifree
                 ! Go to next main memory free line
                 m_mg%ifree=m_mg%ifree+1
                 if(m_mg%ifree.GT.m_mg%ngridmax)then
                    write(*,*)'No more free memory'
                    write(*,*)'in multigrid'
                    write(*,*)'Increase ngridmax'
                    call mdl_abort(mdl)
                 end if
                 ! Insert new grid in hash table
                 call hash_setp(m_mg%grid_dict,hash_father,ichild)

              else

                 ! Otherwise, determine parent processor and use the cache
                 grid_cpu = m_mg%domain(icoarselevel)%get_rank(hk)
                 ! If next cache line is occupied, free it.
                 if(m_mg%occupied(m_mg%free_cache))call destage(mdl,m_mg%ngridmax+m_mg%free_cache)
                 ! Set grid index to a virtual grid in local cache memory
                 ichild=m_mg%ngridmax+m_mg%free_cache
                 m_mg%occupied(m_mg%free_cache)=.true.
                 m_mg%parent_cpu(m_mg%free_cache)=grid_cpu
                 m_mg%dirty(m_mg%free_cache)=.true.
                 m_mg%ghost_parent_grid(m_mg%free_cache)=0
                 m_mg%ghost_parent_cell(m_mg%free_cache)=0
                 ! Go to next free cache line
                 m_mg%free_cache=m_mg%free_cache+1
                 m_mg%ncache=m_mg%ncache+1
                 if(m_mg%free_cache.GT.m_mg%ncachemax)m_mg%free_cache=1
                 if(m_mg%ncache.GT.m_mg%ncachemax)m_mg%ncache=m_mg%ncachemax
                 ! Insert new grid in hash table
                 call hash_setp(m_mg%grid_dict,hash_father,ichild)

              endif

              ! Set oct properties
              m_mg%grid(ichild)%lev=icoarselevel
              m_mg%grid(ichild)%ckey(1:ndim)=cart_key(1:ndim)
              m_mg%grid(ichild)%hkey(1:nhilbert)=hk(1:nhilbert)
              m_mg%grid(ichild)%refined(1:twotondim)=.true.

              ! Set flag arrays
              m_mg%flag1(1:twotondim,ichild)=0
              m_mg%flag2(1:twotondim,ichild)=0

              ! Intitialize gravity variables
              do ind=1,twotondim
                 m_mg%f(ind,1:ndim,ichild)=0
                 m_mg%phi(ind,ichild)=0
              enddo

           end if

        end if

     end do
     ! End loop over coarse neighbors

  end do
  ! End loop over grids

  call close_cache(mdl)

  ! Multigrid oct statistics
  m_mg%tail(icoarselevel)=m_mg%ifree-1
  m_mg%noct(icoarselevel)=m_mg%tail(icoarselevel)-m_mg%head(icoarselevel)+1
  m_mg%noct_used=m_mg%tail(icoarselevel)

  end associate

end subroutine build_mg

subroutine pack_flush_build_mg(mesh,igrid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_small_realdp
  type(mesh_t)::mesh
  integer::igrid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind
  type(msg_small_realdp)::msg

#ifdef GRAV
  do ind=1,twotondim
     msg%realdp(ind)=0.0d0
  end do
#endif

  msg_array=transfer(msg,msg_array)

end subroutine pack_flush_build_mg

subroutine unpack_flush_build_mg(mesh,igrid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim, twotondim
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_small_realdp
  type(mesh_t)::mesh
  integer::igrid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::idim,ind
  type(msg_small_realdp)::msg

  mesh%grid(igrid)%lev=hash_key(0)
  mesh%grid(igrid)%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)

  do ind=1,twotondim
     mesh%grid(igrid)%refined(ind)=.true.
  end do

#ifdef GRAV
  do idim=1,ndim
     do ind=1,twotondim
        mesh%f(ind,idim,igrid)=0.0d0
     end do
  end do
  do ind=1,twotondim
     mesh%phi(ind,igrid)=0.0d0
  end do
#endif

end subroutine unpack_flush_build_mg

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Multigrid cleanup
! ------------------------------------------------------------------------

recursive subroutine r_cleanup_mg(pst)
  use mdl_module
  use amr_parameters, only: twotondim
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_CLEANUP_MG,pst%iUpper+1)
     call r_cleanup_mg(pst%pLower)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
#ifdef _CUDA
     call gpu_clean_mg(pst%s)
#elif defined(_METAL)
     call metal_clean_mg(pst%s)
#else
     call cleanup_mg(pst%s%m_mg)
#endif
  endif

end subroutine r_cleanup_mg

subroutine cleanup_mg(m)
  use amr_commons, only: mesh_t
  use hash
  implicit none
  type(mesh_t)::m

  ! Reset the MG hash table
  call reset_entire_hash(m%grid_dict,.false.)

  ! Restore AMR grid array into its original state
  m%head=1
  m%tail=0
  m%noct=0
  m%ifree=1
  m%noct_tot=0
  m%noct_used=0

end subroutine cleanup_mg

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Initialize mask at fine level into f(:,3)
! ------------------------------------------------------------------------

recursive subroutine r_make_mask(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::ilevel

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_MAKE_MASK,pst%iUpper+1,input_size,0,ilevel)
     call r_make_mask(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
#ifdef _CUDA
     call gpu_make_mask(pst%s,ilevel)
#elif defined(_METAL)
     call metal_make_mask(pst%s,ilevel)
#else
     call make_mask(pst%s%m,ilevel)
#endif
  endif

end subroutine r_make_mask

subroutine make_mask(m,ilevel)
  use amr_parameters, only: twotondim
  use amr_commons, only: mesh_t
  implicit none
  type(mesh_t)::m
  integer,intent(in)::ilevel

  integer::igrid,ind

  ! Init mask to 1.0 on all fine level cells :
  do igrid=m%head(ilevel),m%tail(ilevel)
     do ind=1,twotondim
        m%f(ind,3,igrid)=1.0d0
     end do
  end do

end subroutine make_mask

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Preprocess the fine (AMR) level RHS to account for boundary conditions
!
!  _____#_____
! |     #     |      Cell I is INSIDE active domain (mask > 0)
! |  I  #  O  |      Cell O is OUTSIDE (mask <= 0 or nonexistent cell)
! |_____#_____|      # is the boundary
!       #
!
! phi(I) and phi(O) must BOTH be set at call time, if applicable
! phi(#) is computed from phi(I), phi(O) and the mask values
! If AMR cell O does not exist, phi(O) is computed by interpolation
!
! Sets BC-modified RHS    into f(:,2)
!
! ------------------------------------------------------------------------

recursive subroutine r_make_bc_rhs(pst,input,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  type(level_count_t)::input

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_MAKE_BC_RHS,pst%iUpper+1,input_size,0,input)
     call r_make_bc_rhs(pst%pLower,input,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
#ifdef _CUDA
     call gpu_make_rhs(pst%s,input%ilevel)
#elif defined(_METAL)
     call metal_make_rhs(pst%s,input%ilevel)
#else
     call make_bc_rhs(pst%s,input%ilevel,input%icount)
#endif
  endif

end subroutine r_make_bc_rhs

recursive subroutine r_cmp_rhs_norm2(pst,ilevel,input_size,norm2,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer::ilevel
  real(kind=8)::norm2,next_norm2

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_CMP_RHS_NORM2,pst%iUpper+1,input_size,output_size,ilevel)
     call r_cmp_rhs_norm2(pst%pLower,ilevel,input_size,norm2,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_norm2)
     norm2=norm2+next_norm2
  else
#ifdef _CUDA
     call gpu_rhs_norm2(pst%s, ilevel, norm2)
#elif defined(_METAL)
     call metal_rhs_norm2(pst%s, ilevel, norm2)
#else
     call cmp_rhs_norm2(pst%s%r,pst%s%m,ilevel,norm2)
#endif
  endif

end subroutine r_cmp_rhs_norm2

subroutine cmp_rhs_norm2(r,m,ilevel, norm2)
  use amr_parameters, only: ndim, twotondim
  use amr_commons, only: run_t, mesh_t
  implicit none
  type(run_t)::r
  type(mesh_t)::m
  integer,  intent(in)  :: ilevel
  real(kind=8), intent(out) :: norm2

  real(kind=8) :: dx2
  integer  :: ind, igrid

  ! Set constants
  dx2  = (r%boxlen/2**ilevel)**2
  norm2 = 0.0d0

  ! Loop over grids
  do igrid=m%head(ilevel),m%tail(ilevel)

     ! Loop over cells
     do ind=1,twotondim
        if(m%f(ind,3,igrid)<=0.0)cycle      ! Do not count masked cells
        norm2 = norm2 + m%f(ind,2,igrid)**2
     end do
     ! End loop over cells

  end do
  ! End loop over grids

  norm2 = dx2*norm2

end subroutine cmp_rhs_norm2

subroutine make_bc_rhs(s,ilevel,icount)
  use mdl_module
  use amr_parameters, only: ndim, twondim, twotondim, threetondim
  use ramses_commons, only: ramses_t
  use nbors_utils
  use cache_commons
  use cache
  use phi_fine_cg_module, only: pack_fetch_interpol,unpack_fetch_interpol
  use interpol_phi_module, only: interpol_phi
  use boundaries, only: init_bound_phi
  implicit none
  type(ramses_t)::s

  integer, intent(in) :: ilevel,icount

  integer, dimension(1:3,1:2,1:8) :: iii, jjj
  integer::igrid,idim,ind,igridn,inbor,ig,id
  integer,dimension(1:8,1:8)::ccc
  real(kind=8)::aa,bb,cc,dd,tfrac
  real(kind=8),dimension(1:8)::bbb
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer,dimension(1:threetondim)::ind_nbor
  integer,dimension(1:threetondim)::igrid_nbor
  real(kind=8),dimension(1:twotondim,0:twondim)::phi_nbor,dis_nbor
  integer,dimension(1:3,1:6),save::shift=reshape(&
       & (/-1,0,0,1,0,0,0,-1,0,0,1,0,0,0,-1,0,0,1/),(/3,6/))

  real(kind=8) :: dx, oneoverdx2, phi_b, nb_mask, nb_phi, w
  real(kind=8) :: fourpi, offset
  type(msg_three_realdp)::dummy_three_realdp

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

  ! Set constants
  fourpi = 4.D0*ACOS(-1.0D0)
  if(r%cosmo) fourpi = 1.5D0*g%omega_m*g%aexp

  dx  = r%boxlen/2**ilevel
  oneoverdx2 = 1.0d0/(dx*dx)
  offset = g%rho_tot
  if(any(.not. r%periodic(1:ndim)))offset = 0d0

  iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

  ! CIC method constants
  aa = 1.0D0/4.0D0**ndim
  bb = 3.0D0*aa
  cc = 9.0D0*aa
  dd = 27.D0*aa
  bbb(:)  =(/aa ,bb ,bb ,cc ,bb ,cc ,cc ,dd/)

  ! Sampling positions in the 3x3x3 father cell cube
  ccc(:,1)=(/1 ,2 ,4 ,5 ,10,11,13,14/)
  ccc(:,2)=(/3 ,2 ,6 ,5 ,12,11,15,14/)
  ccc(:,3)=(/7 ,8 ,4 ,5 ,16,17,13,14/)
  ccc(:,4)=(/9 ,8 ,6 ,5 ,18,17,15,14/)
  ccc(:,5)=(/19,20,22,23,10,11,13,14/)
  ccc(:,6)=(/21,20,24,23,12,11,15,14/)
  ccc(:,7)=(/25,26,22,23,16,17,13,14/)
  ccc(:,8)=(/27,26,24,23,18,17,15,14/)

  if (icount .ne. 1 .and. icount .ne. 2)then
     write(*,*)'icount has bad value'
     call mdl_abort(mdl)
  endif

  ! Compute fraction of time steps for interpolation
  if (g%dtold(ilevel-1)>0.0)then
     tfrac=g%dtnew(ilevel)/g%dtold(ilevel-1)*(icount-1)
  else
     tfrac=0.0
  end if

  call open_cache(mdl, m, pack_size=storage_size(dummy_three_realdp)/32, &
       pack=pack_fetch_interpol, unpack=unpack_fetch_interpol, bound=init_bound_phi)

  hash_nbor(0)=ilevel

  ! Loop over grids
  do igrid=m%head(ilevel),m%tail(ilevel)

     ! Get central oct potential
     do ind=1,twotondim
        phi_nbor(ind,0)=m%phi(ind,igrid)
        dis_nbor(ind,0)=m%f(ind,3,igrid)
     end do

     ! Get neighboring octs potential
     do inbor=1,twondim

        ! Get neighboring grid
        hash_nbor(1:ndim)=m%grid(igrid)%ckey(1:ndim)+shift(1:ndim,inbor)

        ! Periodic boundary conditions
        do idim=1,ndim
           if(r%periodic(idim))then
              if(hash_nbor(idim)< m%box_ckey_min(idim,ilevel))hash_nbor(idim)=m%box_ckey_max(idim,ilevel)-1
              if(hash_nbor(idim)>=m%box_ckey_max(idim,ilevel))hash_nbor(idim)=m%box_ckey_min(idim,ilevel)
           endif
        enddo

        ! Get neighbouring grid using read-only cache
        call get_grid(s,hash_nbor,igridn,flush_cache=.false.,fetch_cache=.true.)

        ! If grid exists, then copy into array
        if(igridn>0)then
           do ind=1,twotondim
              phi_nbor(ind,inbor)=m%phi(ind,igridn)
              dis_nbor(ind,inbor)=m%f(ind,3,igridn)
           end do

        ! Otherwise interpolate from coarser level
        else

           ! Get 3**ndim neighbouring parent cells using read-only cache
           call get_threetondim_nbor_parent_cell(s,hash_nbor,igrid_nbor,ind_nbor,flush_cache=.false.,fetch_cache=.true.)
           call interpol_phi(m,igrid_nbor,ind_nbor,ccc,bbb,tfrac,phi_nbor(1,inbor))
           do ind=1,threetondim
              call unlock_cache(m,igrid_nbor(ind))
           end do
           do ind=1,twotondim
              dis_nbor(ind,inbor)=-1.0
           end do
        endif

     end do
     ! End loop over neighboring octs

     ! Loop over cells
     do ind=1,twotondim

        ! Init BC-modified RHS to rho - offset :
        m%f(ind,2,igrid) = fourpi*(m%rho(ind,igrid) - offset)

        ! Do not process masked cells
        if(m%f(ind,3,igrid)<=0.0) cycle

        ! Separate directions
        do idim=1,ndim

           ! Loop over the 2 neighbors
           do inbor=1,2
              id=jjj(idim,inbor,ind); ig=iii(idim,inbor,ind)

              nb_mask=dis_nbor(id,ig)
              if(nb_mask>0.0)cycle

              ! phi(#) interpolated with mask:
              nb_phi = phi_nbor(id,ig)
              w = nb_mask/(nb_mask-m%f(ind,3,igrid)) ! Linear parameter
              phi_b = ((1.0d0-w)*nb_phi + w*m%phi(ind,igrid))

              ! Increment correction for current cell
              m%f(ind,2,igrid) = m%f(ind,2,igrid) - 2.0d0*oneoverdx2*phi_b

           end do
        end do

     end do
  end do

  call close_cache(mdl)

  end associate

end subroutine make_bc_rhs

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

#endif

end module multigrid_fine_commons
