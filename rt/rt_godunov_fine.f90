module rt_godunov_fine_module
#if defined(_CUDA) && defined(DO_RT)
  use gpu_runner, only: gpu_rt_godunov, gpu_rt_set_rtunew, gpu_rt_set_rtuold, gpu_rt_set_emissivity
#endif
contains
!###########################################################
!###########################################################
!###########################################################
!###########################################################
recursive subroutine r_rt_godunov_fine(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::ilevel

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_RT_GODUNOV_FINE,pst%iUpper+1,input_size,0,ilevel)
     call r_rt_godunov_fine(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
#if defined(_CUDA) && defined(DO_RT)
     call gpu_rt_godunov(pst%s,ilevel)
#else
     call rt_godunov_fine(pst%s,ilevel)
#endif
  endif

end subroutine r_rt_godunov_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine rt_godunov_fine(s,ilevel)
  use amr_parameters, only: nsubgrid
  use ramses_commons, only: ramses_t
  use cache_commons
  use cache
  use marshal, only: pack_fetch_refine,unpack_fetch_refine
  use boundaries, only: init_bound_refine
  implicit none
  type(ramses_t)::s
  integer::ilevel
  type(msg_large_realdp)::dummy_large_realdp
  !--------------------------------------------------------------------------
  ! This routine is a wrapper to the second order Godunov solver.
  ! Small grids (2x2x2) are gathered from level ilevel and sent to the
  ! hydro solver. On entry, hydro variables are gathered from array rtuold.
  ! On exit, rtunew has been updated.
  !--------------------------------------------------------------------------
  integer::igrid

  if(s%r%verbose.and.s%g%myid==1)write(*,'("   Entering godunov_fine for level ",I2)')ilevel

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

  call open_cache(mdl, m, pack_size=storage_size(dummy_large_realdp)/32, &
       pack=pack_fetch_refine, unpack=unpack_fetch_refine, &
       init=init_flush_rt_godunov, flush=pack_flush_rt_godunov, &
       combine=unpack_flush_rt_godunov, bound=init_bound_refine)

  ! Loop over active grids by blocks of nsubgrid**ndim octs
  do igrid=m%head(ilevel),m%tail(ilevel),nsubgrid**ndim
     call rt_godfine1(s,igrid,ilevel,m%rt_w%kernel)
  end do

  call close_cache(mdl)

  end associate

end subroutine rt_godunov_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
recursive subroutine r_set_emissivity(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::ilevel

  integer::rID, i

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_SET_EMISSIVITY,pst%iUpper+1,input_size,0,ilevel)
     call r_set_emissivity(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
#if defined(_CUDA) && defined(DO_RT)
     call gpu_rt_set_emissivity(pst%s,ilevel)
#elif defined(DO_RT)
     do i=pst%s%m%head(ilevel),pst%s%m%tail(ilevel)
        pst%s%m%emissivity(:,:,i)=0
     end do
#endif
  endif

end subroutine r_set_emissivity
!###########################################################
!###########################################################
!###########################################################
!###########################################################
recursive subroutine r_set_rtunew(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::ilevel

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_SET_RTUNEW,pst%iUpper+1,input_size,0,ilevel)
     call r_set_rtunew(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
#if defined(_CUDA) && defined(DO_RT)
     ! rtunew lives on the device when the photons are advected or a star can turn them on
     if(pst%s%r%rt_advect .or. (pst%s%r%star .and. pst%s%r%rt_star))then
        call gpu_rt_set_rtunew(pst%s,ilevel)
        return
     endif
#endif
     call set_rtunew(pst%s%m,ilevel)
  endif

end subroutine r_set_rtunew
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine set_rtunew(m,ilevel)
  use amr_parameters, only: ndim, twotondim
  use amr_commons, only: run_t, global_t, mesh_t
  implicit none
  type(mesh_t)::m
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine sets array rtunew to its initial value rtuold before calling
  ! the hydro scheme. rtunew is set to zero in virtual boundaries.
  !--------------------------------------------------------------------------
  integer::i

  ! Set rtunew to rtuold for myid cells
#ifdef DO_RT
  do i=m%head(ilevel),m%tail(ilevel)
     m%rtunew(:,:,i)=m%rtuold(:,:,i)
  end do
#endif

end subroutine set_rtunew
!###########################################################
!###########################################################
!###########################################################
!###########################################################
recursive subroutine r_set_rtuold(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::ilevel
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_SET_RTUOLD,pst%iUpper+1,input_size,0,ilevel)
     call r_set_rtuold(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
#if defined(_CUDA) && defined(DO_RT)
     call gpu_rt_set_rtuold(pst%s,ilevel)
#else
     call set_rtuold(pst%s%r, pst%s%g, pst%s%m, ilevel)
#endif
  endif

end subroutine r_set_rtuold
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine set_rtuold(r, g, m, ilevel)
  use rt_parameters, only: nrtgrp, smallNp
  use amr_parameters, only: ndim, twotondim
  use amr_commons, only: run_t, global_t, mesh_t
  implicit none
  type(run_t) :: r
  type(global_t) :: g
  type(mesh_t) :: m
  integer :: ilevel
  real(kind=8)::Npc,fred
  !---------------------------------------------------------
  ! This routine sets array rtuold to its new value rtunew
  ! after the hydro step.
  !---------------------------------------------------------
  integer :: i, j, ig, iN

  ! Add emissivity source term
#ifdef DO_RT
  if(.not.r%neq_chem)then
     do ig=1,nrtgrp
        iN=1+(ig-1)*(ndim+1)
        do i=m%head(ilevel),m%tail(ilevel)
           m%rtunew(1:twotondim,iN,i)=m%rtunew(1:twotondim,iN,i) + &
                & m%emissivity(1:twotondim,ig,i)*g%dtnew(ilevel)
        end do
     end do
  endif
#endif

  if(r%neq_chem.and.r%rt_smooth)return

 ! Set rtuold to rtunew
#ifdef DO_RT
  do i = m%head(ilevel), m%tail(ilevel)
     m%rtuold(:,:,i) = m%rtunew(:,:,i)
  end do

  ! Make a photon conservation fix (prevent light explosions)
  do ig = 1, nrtgrp
      iN = 1 + (ig-1)*(ndim+1)
      do i = m%head(ilevel), m%tail(ilevel)
        do j = 1, twotondim
          ! No negative photon densities:
          m%rtuold(j,iN,i) = max(m%rtuold(j,iN,i),smallNp)
          Npc=m%rtuold(j,iN,i)*g%rt_c(ilevel)
          ! Reduced flux, should always be .le. 1
          fred = sqrt(sum((m%rtuold(j,iN+1:iN+ndim,i))**2))/Npc
          if(fred .gt. 1d0) then ! Too big so normalize flux to one
            m%rtuold(j,iN+1:iN+ndim,i) &
                = m%rtuold(j,iN+1:iN+ndim,i)/fred
          endif
        end do
      end do
  end do
  ! End photon conservation fix

#endif

end subroutine set_rtuold
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine rt_godfine1(s,ind_grid,ilevel,h)
  use mdl_module
  use amr_parameters, only: ndim, twondim, twotondim
  use rt_parameters, only: nrtvar
  use ramses_commons, only: ramses_t
  use nbors_utils
  use rt_commons
  use hash
  use rt_flux_module, only: rt_unsplit
  implicit none
  type(ramses_t)::s
  integer::ind_grid,ilevel
  type(rt_kernel_t)::h
  !-------------------------------------------------------------------
  ! This routine gathers first RT variables from neighboring grids
  ! to set initial conditions in a 6x6x6 grid. It interpolate from
  ! coarser level missing grid variables. It then calls the
  ! Godunov solver that computes fluxes. These fluxes are zeroed at
  ! coarse-fine boundaries, since contribution from finer levels has
  ! already been taken into account. Conservative variables are updated
  ! and stored in array rtunew(:), both at the current level and at the
  ! coarser level if necessary.
  !-------------------------------------------------------------------
  integer::ivar,idim,ind_son,ind_oct
  integer::icell,inbor,ipass
  integer::igrid,ichild
  integer::i0,j0,k0,i1,j1,k1,i2,j2,k2,i3,j3,k3
  integer::ii0,jj0,kk0,ii1,jj1,kk1
  integer::i1min,i1max,j1min,j1max,k1min,k1max
  integer::ii1min,ii1max,jj1min,jj1max,kk1min,kk1max
  integer::i2min,i2max,j2min,j2max,k2min,k2max
  integer::i3min,i3max,j3min,j3max,k3min,k3max
  integer,dimension(1:ndim)::ckey_corner,ckey
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer,dimension(0:twondim)::icell_nbor
  integer,dimension(0:twondim)::igrid_nbor
  real(kind=8)::dx,oneontwotondim,rt_c_diff,face_flux
  real(kind=8),dimension(0:twondim  ,1:nrtvar)::u1
  real(kind=8),dimension(1:twotondim,1:nrtvar)::u2
  logical::okx,oky,okz,oknbor

  i2min=0; i2max=0; j2min=0; j2max=0; k2min=0; k2max=0
  i3min=1; i3max=1; j3min=1; j3max=1; k3min=1; k3max=1
  okx=.true.; oky=.true.; okz=.true.

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

  oneontwotondim = 1.d0/dble(twotondim)

  ! Mesh spacing in that level
  dx=r%boxlen/2**ilevel

  ! Integer constants
  i1min=h%io1; i1max=h%io2; j1min=h%jo1; j1max=h%jo2; k1min=h%ko1; k1max=h%ko2
#if NDIM>0
  i2max=1; i3min=h%iu1+2; i3max=h%iu2-2
#endif
#if NDIM>1
  j2max=1; j3min=h%ju1+2; j3max=h%ju2-2
#endif
#if NDIM>2
  k2max=1; k3min=h%ku1+2; k3max=h%ku2-2
#endif

  !---------------------
  ! Gather hydro stencil
  !---------------------
  hash_nbor(0)=m%grid(ind_grid)%lev
  ckey_corner(1:ndim)=(m%grid(ind_grid)%ckey(1:ndim)/(i1max-1))*(i1max-1)
  ind_oct=ind_grid
  h%inkernel=.false.

  ! Loop over 3x3x3 neighboring father cells using 7 passes
  do ipass = 1, 1+2*ndim

  if(ipass == 1)then
     ii1min = i1min+1; ii1max = i1max-1
     jj1min = j1min; jj1max = j1max
#if NDIM>1
     jj1min = jj1min+1; jj1max = jj1max-1
#endif
     kk1min = k1min; kk1max = k1max;
#if NDIM>2
     kk1min = kk1min+1; kk1max = kk1max-1
#endif
  endif
  if(ipass == 2)then
     ii1min = i1min; ii1max = i1min
  endif
  if(ipass == 3)then
     ii1min = i1max; ii1max = i1max
  endif
  if(ipass == 4)then
     ii1min = i1min+1; ii1max = i1max-1
     jj1min = j1min; jj1max = j1min
  endif
  if(ipass == 5)then
     jj1min = j1max; jj1max = j1max
  endif
  if(ipass == 6)then
     jj1min = j1min+1; jj1max = j1max-1
     kk1min = k1min; kk1max = k1min
  endif
  if(ipass == 7)then
     kk1min = k1max; kk1max = k1max
  endif

  do k1=kk1min,kk1max
#if NDIM>2
     okz=(k1>k1min.and.k1<k1max)
#endif
     do j1=jj1min,jj1max
#if NDIM>1
        oky=(j1>j1min.and.j1<j1max)
#endif
        do i1=ii1min,ii1max
           okx=(i1>i1min.and.i1<i1max)

           !--------------------
           ! For inner octs only
           !--------------------
           if(okx.and.oky.and.okz)then

              ! Compute relative Cartesian key
              ckey(1:ndim)=m%grid(ind_oct)%ckey(1:ndim)-ckey_corner(1:ndim)

              ! Store grid index
              ii1=1; jj1=1; kk1=1
              ii1=ckey(1)+1
#if NDIM>1
              jj1=ckey(2)+1
#endif
#if NDIM>2
              kk1=ckey(3)+1
#endif
              h%childloc(ii1,jj1,kk1)=ind_oct
              h%inkernel(ii1,jj1,kk1)=.true.
              h%gridloc(ii1,jj1,kk1)=0
              h%cellloc(ii1,jj1,kk1)=0

              ! Loop over 2x2x2 cells
              do k2=k2min,k2max
                 do j2=j2min,j2max
                    do i2=i2min,i2max
                       ind_son=1+i2+2*j2+4*k2
                       i3=1; j3=1; k3=1
                       i3=1+2*(ii1-1)+i2
#if NDIM>1
                       j3=1+2*(jj1-1)+j2
#endif
#if NDIM>2
                       k3=1+2*(kk1-1)+k2
#endif
                       ! Gather RT variables
                       do ivar=1,nrtvar
#ifdef DO_RT
                          h%rtuloc(i3,j3,k3,ivar)=m%rtuold(ind_son,ivar,ind_oct)
#endif
                       end do
                       ! Gather refinement flag
                       h%okloc(i3,j3,k3)=m%grid(ind_oct)%refined(ind_son)
                    end do
                 end do
              end do

              ! Go to next inner oct
              ind_oct=ind_oct+1

           !-----------------------
           ! For boundary octs only
           !-----------------------
           else

              ! Compute neighboring grid Cartesian index
              hash_nbor(1)=ckey_corner(1)+i1-1
#if NDIM>1
              hash_nbor(2)=ckey_corner(2)+j1-1
#endif
#if NDIM>2
              hash_nbor(3)=ckey_corner(3)+k1-1
#endif
              ! Periodic boundary conditions
              do idim=1,ndim
                 if(r%periodic(idim))then
                    if(hash_nbor(idim)< m%box_ckey_min(idim,ilevel))hash_nbor(idim)=m%box_ckey_max(idim,ilevel)-1
                    if(hash_nbor(idim)>=m%box_ckey_max(idim,ilevel))hash_nbor(idim)=m%box_ckey_min(idim,ilevel)
                 endif
              enddo

              ! Set grid index to null
              icell=0
              igrid=0
              do inbor=0,twondim
                 igrid_nbor(inbor)=0
              end do
              ! Get neighboring grid index with read-only cache
              call get_grid(s,hash_nbor,ichild,flush_cache=.false.,fetch_cache=.true.,lock=.true.)

              !----------------------------------------------------
              ! If grid does not exist, interpolate RT variables
              !----------------------------------------------------
              if(ichild==0)then

                 ! Get parent father cell with read-write cache
                 call get_parent_cell(s,hash_nbor,igrid,icell,flush_cache=.true.,fetch_cache=.true.,lock=.true.)
                 if(igrid==0)then
                    write(*,*)'RT-GODUNOV: parent_cell should exist'
                    write(*,*)'PE ',g%myid,hash_nbor
                    call mdl_abort(mdl)
                 endif

                 ! Get 2ndim neighboring father cells with read-write cache
                 ! Note that possible cache grids are locked inside the routine
                 call get_twondim_nbor_parent_cell(s,hash_nbor,igrid_nbor,icell_nbor,flush_cache=.true.,fetch_cache=.true.)
                 oknbor=.true.
                 do inbor=0,twondim
                    oknbor=oknbor.and.(igrid_nbor(inbor)>0)
                 end do
                 if(.not. oknbor)then
                    write(*,*)"RT-GODUNOV: parent neighbors should exist"
                    write(*,*)'PE ',g%myid,hash_nbor
                    write(*,*)igrid_nbor(0)
                    do idim=1,ndim
                       write(*,*)igrid_nbor(2*idim-1)
                       write(*,*)igrid_nbor(2*idim)
                    end do
                    call mdl_abort(mdl)
                 endif

                 ! Gather RT variables
                 do inbor=0,twondim
                    do ivar=1,nrtvar
#ifdef DO_RT
                       u1(inbor,ivar)=m%rtuold(icell_nbor(inbor),ivar,igrid_nbor(inbor))
#endif
                    end do
                 end do
                 ! Interpolate using rt variables
                 call interpol_rt(u1,u2,r%interpol_type)
              endif

              ! Store grid index
              h%childloc(i1,j1,k1)=ichild
              h%inkernel(i1,j1,k1)=.true.
              h%gridloc (i1,j1,k1)=igrid
              h%cellloc (i1,j1,k1)=icell
              do inbor=1,twondim
                 h%nborloc(i1,j1,k1,inbor)=igrid_nbor(inbor)
              end do

              ! Loop over 2x2x2 cells
              do k2=k2min,k2max
                 do j2=j2min,j2max
                    do i2=i2min,i2max
                       ind_son=1+i2+2*j2+4*k2
                       i3=1; j3=1; k3=1
                       i3=1+2*(i1-1)+i2
#if NDIM>1
                       j3=1+2*(j1-1)+j2
#endif
#if NDIM>2
                       k3=1+2*(k1-1)+k2
#endif
                       !----------------------------------------------------
                       ! If neighboring grid exists, use refined grid values
                       !----------------------------------------------------
                       if(ichild>0)then

                          ! Gather refinement flag
                          h%okloc(i3,j3,k3)=m%grid(ichild)%refined(ind_son)
                          ! Gather RT variables
                          do ivar=1,nrtvar
#ifdef DO_RT
                             h%rtuloc(i3,j3,k3,ivar)=m%rtuold(ind_son,ivar,ichild)
#endif
                          end do

                       !-----------------------------------------------------------
                       ! If neighboring grid doesn't exist, use interpolated values
                       !-----------------------------------------------------------
                       else

                          ! Gather refinement flag
                          h%okloc(i3,j3,k3)=.false.
                          ! Gather interpolated hydro variables
                          do ivar=1,nrtvar
#ifdef DO_RT
                             h%rtuloc(i3,j3,k3,ivar)=u2(ind_son,ivar)
                             if(mod(ivar,ndim+1)==1) then
                                ! Variable light speed correction for coarser level
                                h%rtuloc(i3,j3,k3,ivar)           &
                                  = h%rtuloc(i3,j3,k3,ivar)       &
                                  * g%rt_c(ilevel-1)/g%rt_c(ilevel)
                             endif
#endif
                          end do
                       end if

                    end do
                 end do
              end do
              ! End loop over 2x2x2 cells
           endif
        end do
     end do
  end do
  ! End over octs
  end do

  !-------------------------------------------------
  ! Compute flux using second-order Godunov method
  !-------------------------------------------------
#ifdef DO_RT
  call rt_unsplit(h%rtuloc,h%rtflux,h%cFlx,        &
       & g%rt_c(ilevel),dx,dx,dx,g%dtnew(ilevel),  &
       & h%iu1,h%iu2,h%ju1,h%ju2,h%ku1,h%ku2,      &
       & h%if1,h%if2,h%jf1,h%jf2,h%kf1,h%kf2)
#endif

  !-------------------------------------------------
  ! Reset flux along direction at refined interfaces
  !-------------------------------------------------
  if(r%rt_nsubcycle.eq.1) then
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     do k3=k3min,k3max+k0
        do j3=j3min,j3max+j0
           do i3=i3min,i3max+i0
              if(h%okloc(i3-i0,j3-j0,k3-k0) .or. h%okloc(i3,j3,k3))then
                 do ivar=1,nrtvar
#ifdef DO_RT
                    h%rtflux(i3,j3,k3,ivar,idim)=0.0d0
#endif
                 end do
              end if
           end do
        end do
     end do
  end do
  endif
  !--------------------------------------
  ! Conservative update at level ilevel
  !--------------------------------------
  ! Loop over dimensions
  do idim=1,ndim
     i0=0; j0=0; k0=0
     ii0=0; jj0=0; kk0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     ii0=1
#if NDIM>1
     jj0=1
#endif
#if NDIM>2
     kk0=1
#endif
     ! Loop over inner octs
     do k1=k1min+kk0,k1max-kk0
        do j1=j1min+jj0,j1max-jj0
           do i1=i1min+ii0,i1max-ii0
              ! Get oct index
              ichild=h%childloc(i1,j1,k1)
              ! Loop over cells
              do k2=k2min,k2max
                 do j2=j2min,j2max
                    do i2=i2min,i2max
                       ind_son=1+i2+2*j2+4*k2
                       i3=1; j3=1; k3=1
                       i3=1+2*(i1-1)+i2
#if NDIM>1
                       j3=1+2*(j1-1)+j2
#endif
#if NDIM>2
                       k3=1+2*(k1-1)+k2
#endif
                       ! Update conservative variables new state vector
                       if(.not. h%okloc(i3,j3,k3)) then
                        do ivar=1,nrtvar
#ifdef DO_RT
                          m%rtunew(ind_son,ivar,ichild)=m%rtunew(ind_son,ivar,ichild)+ &
                               & (h%rtflux(i3   ,j3   ,k3   ,ivar,idim) &
                               & -h%rtflux(i3+i0,j3+j0,k3+k0,ivar,idim))
#endif
                        end do
                       endif
                    end do
                 end do
              end do
           end do
        end do
     end do
  end do

  ! If sitting in coarsest level, skip.
  if(ilevel>r%levelmin)then

  !--------------------------------------
  ! Conservative update at level ilevel-1
  !--------------------------------------
  if(r%rt_nsubcycle.eq.1) then
  ! Loop over dimensions
  do idim=1,ndim
     i0=0; j0=0; k0=0
     ii0=0; jj0=0; kk0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     ii0=1
#if NDIM>1
     jj0=1
#endif
#if NDIM>2
     kk0=1
#endif
     ii1min=i1min+ii0; ii1max=i1max-ii0
     jj1min=j1min+jj0; jj1max=j1max-jj0
     kk1min=k1min+kk0; kk1max=k1max-kk0
     !----------------------
     ! Left flux at boundary
     !----------------------
     if(idim==1)then
        ii1min=i1min; ii1max=i1min
     endif
     if(idim==2)then
        jj1min=j1min; jj1max=j1min
     endif
     if(idim==3)then
        kk1min=k1min; kk1max=k1min
     endif
     ! Loop over outer octs on left face
     do k1=kk1min,kk1max
        do j1=jj1min,jj1max
           do i1=ii1min,ii1max
              ! Get oct index
              ichild=h%childloc(i1,j1,k1)
              ! Check that parent cell is not refined
              if(ichild==0)then
                 ! Get parent cell index
                 igrid=h%gridloc(i1,j1,k1)
                 icell=h%cellloc(i1,j1,k1)
                 ! Sum the fluxes of the inner cell left faces, then update once
                 do ivar=1,nrtvar
                    ! For VSLA, when updating coarser level, rescale radiation flux
                    ! to the expression which would be seen from there.
                    rt_c_diff = g%rt_c(ilevel-1)/g%rt_c(ilevel)
                    if (mod(ivar,ndim+1)==1) then
                       rt_c_diff=1.d0
                    end if
                    face_flux=0.0d0
                    do k2=k2min,k2max-k0
                       do j2=j2min,j2max-j0
                          do i2=i2min,i2max-i0
                             i3=1; j3=1; k3=1
                             i3=1+2*(i1+i0-1)+i2
#if NDIM>1
                             j3=1+2*(j1+j0-1)+j2
#endif
#if NDIM>2
                             k3=1+2*(k1+k0-1)+k2
#endif
#ifdef DO_RT
                             face_flux=face_flux+h%rtflux(i3,j3,k3,ivar,idim)
#endif
                          end do
                       end do
                    end do
                    ! Conservative update of new state variables
#ifdef DO_RT
                    m%rtunew(icell,ivar,igrid)=m%rtunew(icell,ivar,igrid) &
                         & -face_flux * oneontwotondim * rt_c_diff
#endif
                 end do
              endif
           end do
        end do
     end do
     !-----------------------
     ! Right flux at boundary
     !-----------------------
     if(idim==1)then
        ii1min=i1max; ii1max=i1max
     endif
     if(idim==2)then
        jj1min=j1max; jj1max=j1max
     endif
     if(idim==3)then
        kk1min=k1max; kk1max=k1max
     endif
     ! Loop over outer octs on right face
     do k1=kk1min,kk1max
        do j1=jj1min,jj1max
           do i1=ii1min,ii1max
              ! Get oct index
              ichild=h%childloc(i1,j1,k1)
              ! Check that parent cell is not refined
              if(ichild==0)then
                 ! Get parent cell index
                 igrid=h%gridloc(i1,j1,k1)
                 icell=h%cellloc(i1,j1,k1)
                 ! Sum the fluxes of the inner cell right faces, then update once
                 do ivar=1,nrtvar
                    ! Rescale for VSLA, as above
                    rt_c_diff = g%rt_c(ilevel-1)/g%rt_c(ilevel)
                    if (mod(ivar,ndim+1)==1) then
                       rt_c_diff=1.d0
                    end if
                    face_flux=0.0d0
                    do k2=k2min+k0,k2max
                       do j2=j2min+j0,j2max
                          do i2=i2min+i0,i2max
                             i3=1; j3=1; k3=1
                             i3=1+2*(i1-i0-1)+i2
#if NDIM>1
                             j3=1+2*(j1-j0-1)+j2
#endif
#if NDIM>2
                             k3=1+2*(k1-k0-1)+k2
#endif
#ifdef DO_RT
                             face_flux=face_flux+h%rtflux(i3+i0,j3+j0,k3+k0,ivar,idim)
#endif
                          end do
                       end do
                    end do
                    ! Conservative update of new state variables
#ifdef DO_RT
                    m%rtunew(icell,ivar,igrid)=m%rtunew(icell,ivar,igrid) &
                         & +face_flux * oneontwotondim * rt_c_diff
#endif
                 end do
                 ! End loop over faces
              endif
           end do
        end do
     end do
     ! End loop over boundary octs
  end do
  ! End loop over dimensions
  endif ! if(r%rt_nsubcycle.eq.1)

  endif

  !----------------
  ! Unlock all octs
  !----------------
  do k1=k1min,k1max
     do j1=j1min,j1max
        do i1=i1min,i1max
           ! Check if in kernel
           if(.not.h%inkernel(i1,j1,k1))cycle
           ! Get oct index
           ichild=h%childloc(i1,j1,k1)
           ! Check that parent cell is not refined
           if(ichild>0)then
              call unlock_cache(m,ichild)
           else
              ! Get parent cell index
              igrid=h%gridloc(i1,j1,k1)
              call unlock_cache(m,igrid)
              ! Get neighbouring parent oct index
              do inbor=1,twondim
                 igrid=h%nborloc(i1,j1,k1,inbor)
                 call unlock_cache(m,igrid)
              end do
           endif
        end do
     end do
  end do

  end associate

end subroutine rt_godfine1
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine init_flush_rt_godunov(mesh,igrid,hash_key)
  use amr_parameters, only: ndim, twotondim
  use rt_parameters, only: nrtvar
  use amr_commons, only: mesh_t
  type(mesh_t)::mesh
  integer::igrid
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind,ivar

  mesh%grid(igrid)%lev=hash_key(0)
  mesh%grid(igrid)%ckey(1:ndim)=hash_key(1:ndim)

#ifdef DO_RT
  do ivar=1,nrtvar
     do ind=1,twotondim
        mesh%rtunew(ind,ivar,igrid)=0.0d0
     enddo
  enddo
#endif

end subroutine init_flush_rt_godunov
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine pack_flush_rt_godunov(mesh,igrid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use rt_parameters, only: nrtvar
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_large_realdp
  type(mesh_t)::mesh
  integer::igrid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind,ivar
  type(msg_large_realdp)::msg

#ifdef DO_RT
  do ivar=1,nrtvar
     do ind=1,twotondim
        msg%realdp_rt(ind,ivar)=mesh%rtunew(ind,ivar,igrid)
     end do
  end do
#endif

  msg_array=transfer(msg,msg_array)

end subroutine pack_flush_rt_godunov
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine unpack_flush_rt_godunov(mesh,igrid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim, twotondim
  use rt_parameters, only: nrtvar
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_large_realdp
  type(mesh_t)::mesh
  integer::igrid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind,ivar
  type(msg_large_realdp)::msg

  mesh%grid(igrid)%lev=hash_key(0)
  mesh%grid(igrid)%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)

#ifdef DO_RT
  do ivar=1,nrtvar
     do ind=1,twotondim
        mesh%rtunew(ind,ivar,igrid)=mesh%rtunew(ind,ivar,igrid)+msg%realdp_rt(ind,ivar)
     end do
  end do
#endif

end subroutine unpack_flush_rt_godunov
!###########################################################
!###########################################################
!###########################################################
!###########################################################
end module rt_godunov_fine_module
