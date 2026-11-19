module force_fine_module

  use multigrid_fine_coarse, only: level_count_t

#ifdef _CUDA
  use gpu_runner, only: gpu_gradient_phi, gpu_epot, gpu_rhomax, gpu_force_analytic, gpu_check_nan
#endif
#ifdef _METAL
  use metal_runner, only: metal_gradient_phi, metal_cmp_epot, metal_cmp_rhomax
#endif

contains
#ifdef GRAV
!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine m_force_fine(pst,ilevel,icount)
  use amr_parameters, only: ndim, twotondim, nvector
  use ramses_commons, only: pst_t
  implicit none
  type(pst_t)::pst
  integer::ilevel,icount
  !----------------------------------------------------------
  ! This routine computes the gravitational acceleration,
  ! the maximum density rho_max, and the potential energy
  !----------------------------------------------------------
  integer::dummy(2)
  real(kind=8)::rhomax,epot
  type(level_count_t)::level_count

  if(pst%s%m%noct_tot(ilevel)==0)return
  if(pst%s%r%verbose)write(*,'("   Entering force_fine for level ",I2)')ilevel

  if(pst%s%r%gravity_type>0)then
     ! Compute analytical gravity force
     call r_force_analytic(pst,ilevel,1)
  else
     ! Compute gradient of potential
     level_count%ilevel=ilevel
     level_count%icount=icount
     call r_gradient_phi(pst,level_count,2)
     ! Add external acceleration
     if(pst%s%r%gravity_type<0)then
        call r_force_analytic(pst,ilevel,1)
     endif
  endif
  if(pst%s%r%verbose)write(*,'("   Gradient phi done for level ",I2)')ilevel

  ! Compute gravity potential energy
  call r_compute_epot(pst,ilevel,1,epot,2)
  pst%s%g%epot_tot=pst%s%g%epot_tot+epot
  if(pst%s%r%verbose)write(*,'("   Potential energy done for level ",I2)')ilevel

  ! Compute maximum mass density
  call r_compute_rhomax(pst,ilevel,1,rhomax,2)
  pst%s%g%rho_max(ilevel)=rhomax
  if(pst%s%r%verbose)write(*,'("   Maximum density done for level ",I2)')ilevel

end subroutine m_force_fine
!#########################################################
!#########################################################
!#########################################################
!#########################################################
recursive subroutine r_force_analytic(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::ilevel

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_FORCE_ANALYTIC,pst%iUpper+1,input_size,0,ilevel)
     call r_force_analytic(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
#ifdef _CUDA
     call gpu_force_analytic(pst%s,ilevel)
#else
     call force_analytic(pst%s%r,pst%s%g,pst%s%m,ilevel)
#endif
  endif

end subroutine r_force_analytic
!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine force_analytic(r,g,m,ilevel)
  use amr_parameters, only: ndim, twotondim, nvector
  use amr_commons, only: run_t, global_t, mesh_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer::ilevel
  !-------------------------------------
  ! Compute analytical gravity force
  !-------------------------------------
  integer::igrid,ind,i,ngrid,idim,nstride
  real(kind=8)::dx
  real(kind=8),dimension(1:nvector,1:ndim)::xx,ff

  ! Mesh size at level ilevel in code units
  dx=r%boxlen/2**ilevel

  ! Initialize force zero
  if(r%gravity_type>0)then
     m%f(1:twotondim,1:ndim,m%head(ilevel):m%tail(ilevel))=0d0
  endif

  ! Loop over grids by vector sweeps
  do igrid=m%head(ilevel),m%tail(ilevel),nvector
     ngrid=MIN(nvector,m%tail(ilevel)-igrid+1)

     ! Loop over cells
     do ind=1,twotondim

        ! Compute cell centre position in code units
        do idim=1,ndim
           nstride=2**(idim-1)
           do i=1,ngrid
              xx(i,idim)=(2*m%grid(igrid+i-1)%ckey(idim)+MOD((ind-1)/nstride,2)+0.5)*dx-m%skip(idim)
           end do
        end do

        ! Call analytical gravity routine
        ff=0d0
        call gravana(r,g,xx,ff,dx,ngrid)

        ! Scatter variables to main memory
        do idim=1,ndim
           do i=1,ngrid
              m%f(ind,idim,igrid+i-1)=m%f(ind,idim,igrid+i-1)+ff(i,idim)
           end do
        end do

     end do
     ! End loop over cells

  end do
  ! End loop over grid

end subroutine force_analytic
!#########################################################
!#########################################################
!#########################################################
!#########################################################
recursive subroutine r_gradient_phi(pst,input,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  type(level_count_t)::input

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_GRADIENT_PHI,pst%iUpper+1,input_size,0,input)
     call r_gradient_phi(pst%pLower,input,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
#ifdef _CUDA
     call gpu_gradient_phi(pst%s,input%ilevel,input%icount)
     call gpu_check_nan(pst%s, input%ilevel, 'force')
#elif defined(_METAL)
     call metal_gradient_phi(pst%s,input%ilevel)
#else
     call gradient_phi(pst%s,input%ilevel,input%icount)
#endif
  endif

end subroutine r_gradient_phi
!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine gradient_phi(s,ilevel,icount)
  use mdl_module
  use amr_parameters, only: ndim, twondim, twotondim, threetondim, nvector
  use ramses_commons, only: ramses_t
  use nbors_utils
  use cache_commons
  use cache
  use interpol_phi_module, only: interpol_phi
  use phi_fine_cg_module, only: pack_fetch_interpol, unpack_fetch_interpol
  use boundaries, only: init_bound_phi
  implicit none
  type(ramses_t)::s
  integer::ilevel,icount
  !-------------------------------------------------
  ! This routine compute the 3-force for all cells
  ! in the current level grids, using a
  ! 5 nodes kernel (5 points FDA).
  !-------------------------------------------------
  integer::i_nbor,igrid,idim,ind,igridn
  integer::id1,id2,id3,id4
  integer::ig1,ig2,ig3,ig4
  integer,dimension(1:3,1:4,1:8)::ggg,hhh
  integer,dimension(1:8,1:8)::ccc
  integer,dimension(1:threetondim)::igrid_nbor,ind_nbor
  integer,dimension(1:3,1:6)::shift=reshape(&
       & (/-1,0,0,1,0,0,0,-1,0,0,1,0,0,0,-1,0,0,1/),(/3,6/))
  integer(kind=8),dimension(0:ndim)::hash_nbor
  real(kind=8),dimension(1:8)::bbb
  real(kind=8)::dx,a,b,aa,bb,cc,dd,tfrac
  real(kind=8)::phi1,phi2,phi3,phi4
  real(kind=8),dimension(1:twotondim,0:twondim)::phi_nbor
  type(msg_three_realdp)::dummy_three_realdp

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

  ! Mesh size at level ilevel in code units
  dx=r%boxlen/2**ilevel

  ! Rescaling factor
  a=0.50D0*4.0D0/3.0D0/dx
  b=0.25D0*1.0D0/3.0D0/dx
  !   |dim
  !   | |node
  !   | | |cell
  !   v v v
  ggg(1,1,1:8)=(/1,0,1,0,1,0,1,0/); hhh(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,2,1:8)=(/0,2,0,2,0,2,0,2/); hhh(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,3,1:8)=(/1,1,1,1,1,1,1,1/); hhh(1,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(1,4,1:8)=(/2,2,2,2,2,2,2,2/); hhh(1,4,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(2,1,1:8)=(/3,3,0,0,3,3,0,0/); hhh(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,2,1:8)=(/0,0,4,4,0,0,4,4/); hhh(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,3,1:8)=(/3,3,3,3,3,3,3,3/); hhh(2,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(2,4,1:8)=(/4,4,4,4,4,4,4,4/); hhh(2,4,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(3,1,1:8)=(/5,5,5,5,0,0,0,0/); hhh(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,2,1:8)=(/0,0,0,0,6,6,6,6/); hhh(3,2,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,3,1:8)=(/5,5,5,5,5,5,5,5/); hhh(3,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(3,4,1:8)=(/6,6,6,6,6,6,6,6/); hhh(3,4,1:8)=(/1,2,3,4,5,6,7,8/)

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
  if (g%dtold(ilevel-1)>0.0d0)then
     tfrac=g%dtnew(ilevel)/g%dtold(ilevel-1)*(icount-1)
  else
     tfrac=0.0
  end if

  call open_cache(mdl, m, pack_size=storage_size(dummy_three_realdp)/32, &
       pack=pack_fetch_interpol, unpack=unpack_fetch_interpol, &
       bound=init_bound_phi)

  hash_nbor(0)=ilevel

  ! Loop over grids
  do igrid=m%head(ilevel),m%tail(ilevel)

     ! Get central oct potential
     do ind=1,twotondim
        phi_nbor(ind,0)=m%phi(ind,igrid)
     end do

     ! Get neighboring octs potential
     do i_nbor=1,twondim

        ! Get neighboring grid
        hash_nbor(1:ndim)=m%grid(igrid)%ckey(1:ndim)+shift(1:ndim,i_nbor)

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
              phi_nbor(ind,i_nbor)=m%phi(ind,igridn)
           end do

        ! Otherwise interpolate from coarser level
        else
           ! Get 3**ndim parent cell using read-only cache
           call get_threetondim_nbor_parent_cell(s,hash_nbor,igrid_nbor,ind_nbor,flush_cache=.false.,fetch_cache=.true.)
           call interpol_phi(m,igrid_nbor,ind_nbor,ccc,bbb,tfrac,phi_nbor(1,i_nbor))
           do ind=1,threetondim
              call unlock_cache(m,igrid_nbor(ind))
           end do
        endif

     end do
     ! End loop over neighboring octs

     ! Loop over cells
     do ind=1,twotondim

        ! Loop over dimensions
        do idim=1,ndim

           ! Gather nodes indices
           id1=hhh(idim,1,ind); ig1=ggg(idim,1,ind)
           id2=hhh(idim,2,ind); ig2=ggg(idim,2,ind)
           id3=hhh(idim,3,ind); ig3=ggg(idim,3,ind)
           id4=hhh(idim,4,ind); ig4=ggg(idim,4,ind)

           ! Gather potential
           phi1=phi_nbor(id1,ig1)
           phi2=phi_nbor(id2,ig2)
           phi3=phi_nbor(id3,ig3)
           phi4=phi_nbor(id4,ig4)

           ! Compute acceleration
           m%f(ind,idim,igrid)=a*(phi1-phi2)-b*(phi3-phi4)

        end do
        ! End loop over dimensions

     end do
     ! End loop over cells

  end do
  ! End loop over grids

  call close_cache(mdl)

  end associate

end subroutine gradient_phi
!#########################################################
!#########################################################
!#########################################################
!#########################################################
recursive subroutine r_compute_epot(pst,ilevel,input_size,epot,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer::ilevel
  real(kind=8)::epot,next_epot

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_COMPUTE_EPOT,pst%iUpper+1,input_size,output_size,ilevel)
     call r_compute_epot(pst%pLower,ilevel,input_size,epot,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_epot)
     epot=epot+next_epot
  else
#ifdef _CUDA
     call gpu_epot(pst%s,ilevel,epot)
#elif defined(_METAL)
     call metal_cmp_epot(pst%s,ilevel,epot)
#else
     call compute_epot(pst%s%r,pst%s%g,pst%s%m,ilevel,epot)
#endif
  endif

end subroutine r_compute_epot
!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine compute_epot(r,g,m,ilevel,epot)
  use amr_parameters, only: ndim, twotondim
  use amr_commons, only: run_t, global_t, mesh_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer::ilevel
  real(kind=8)::epot
  !----------------------------------------------------------
  ! This routine computes the potential energy
  !----------------------------------------------------------
  integer::igrid,ind,idim
  real(kind=8)::dx,fact,fourpi

  ! Mesh size at level ilevel in code units
  dx=r%boxlen/2**ilevel
  ! Local constants
  fourpi=4.0D0*ACOS(-1.0D0)
  if(r%cosmo)fourpi=1.5D0*g%omega_m*g%aexp
  fact=-dx**ndim/fourpi/2.0D0

  ! Compute gravity potential
  epot=0D0

  ! Loop over myid grids by vector sweeps
  do igrid=m%head(ilevel),m%tail(ilevel)
     ! Loop over cells
     do ind=1,twotondim
        ! Loop over dimensions
        do idim=1,ndim
           if(.not.m%grid(igrid)%refined(ind))then
              epot=epot+fact*m%f(ind,idim,igrid)**2
           endif
        end do
     end do
     ! End loop over cells
  end do
  ! End loop over grids

end subroutine compute_epot
!#########################################################
!#########################################################
!#########################################################
!#########################################################
recursive subroutine r_compute_rhomax(pst,ilevel,input_size,rhomax,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer::ilevel

  real(kind=8)::rhomax,next_rhomax
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_COMPUTE_RHOMAX,pst%iUpper+1,input_size,output_size,ilevel)
     call r_compute_rhomax(pst%pLower,ilevel,input_size,rhomax,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_rhomax)
     rhomax=MAX(rhomax,next_rhomax)
  else
#ifdef _CUDA
     call gpu_rhomax(pst%s,ilevel,rhomax)
#elif defined(_METAL)
     call metal_cmp_rhomax(pst%s,ilevel,rhomax)
#else
     call compute_rhomax(pst%s%r,pst%s%g,pst%s%m,ilevel,rhomax)
#endif
  endif

end subroutine r_compute_rhomax
!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine compute_rhomax(r,g,m,ilevel,rhomax)
  use amr_parameters, only: ndim, twotondim
  use amr_commons, only: run_t, global_t, mesh_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer::ilevel
  real(kind=8)::rhomax
  !----------------------------------------------------------
  ! This routine computes the potential energy
  !----------------------------------------------------------
  integer::igrid,ind,idim

  ! Compute maximum total mass density
  rhomax=0D0

  ! Loop over myid grids by vector sweeps
  do igrid=m%head(ilevel),m%tail(ilevel)
     ! Loop over cells
     do ind=1,twotondim
        rhomax=MAX(rhomax,dble(abs(m%rho(ind,igrid))))
     end do
     ! End loop over cells
  end do
  ! End loop over grids

end subroutine compute_rhomax
!#########################################################
!#########################################################
!#########################################################
!#########################################################
#endif
end module force_fine_module
