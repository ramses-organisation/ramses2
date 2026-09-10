module init_amr_module

contains

!###############################################
!###############################################
!###############################################
!###############################################
recursive subroutine r_set_add(pst,iUpper,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size

  integer::n,iLower,iMiddle,iUpper
  integer::rID

  associate(mdl=>pst%s%mdl)

  iLower = mdl_self(mdl)-1
  n = iUpper - iLower
  iMiddle = (iUpper + iLower) / 2

  if(n>1)then
     pst%iUpper = iMiddle
     pst%nLower = iMiddle - iLower
     pst%nUpper = iUpper - iMiddle
     allocate(pst%pLower)
     pst%pLower%s => pst%s
     rID = mdl_send_request(mdl,MDL_SET_ADD,pst%iUpper+1,input_size,0,iUpper)
     call r_set_add(pst%pLower,iMiddle,input_size)
     call mdl_get_reply(mdl,rID,0)
  end if

  end associate

end subroutine r_set_add
!###############################################
!###############################################
!###############################################
!###############################################
recursive subroutine r_init_amr(pst)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
#ifdef _CUDA
  use gpu_manager
  use clump_device, only: gpu_allocate_clump
#endif
#ifdef _METAL
  use metal_runner, only: metal_allocate_amr, metal_allocate_grav
#endif
  implicit none
  type(pst_t)::pst

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_INIT_AMR,pst%iUpper+1)
     call r_init_amr(pst%pLower)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     allocate(pst%s%m)
     call init_amr(pst%s%r,pst%s%g,pst%s%m,'amr')
     if(pst%s%r%poisson)then
        allocate(pst%s%m_mg)
        call init_amr(pst%s%r,pst%s%g,pst%s%m_mg,'mg')
     endif
     if(pst%s%r%clump_finder .or. pst%s%r%sink_form)then
        allocate(pst%s%c)
     endif
     call init_params(pst%s%mdl,pst%s%r,pst%s%g)
#ifdef _CUDA
     call gpu_allocate_amr(pst%s)
     if(pst%s%r%poisson) call gpu_allocate_mg(pst%s)
     if(pst%s%r%clump_finder .or. pst%s%r%sink_form)then
        call gpu_allocate_clump(pst%s)
     endif
#endif
#ifdef _METAL
     call metal_allocate_amr(pst%s)
     if(pst%s%r%poisson) call metal_allocate_grav(pst%s)
#endif
  endif

end subroutine r_init_amr
!###############################################
!###############################################
!###############################################
!###############################################
subroutine init_amr(r,g,m,type)
  use amr_parameters, ONLY: nhilbert, ndim, twotondim, threetondim, dp
  use hydro_parameters, ONLY: nvar
  use rt_parameters, ONLY: nrtvar, nrtgrp
  use cr_parameters, ONLY: ncruvar
  use amr_commons, ONLY: run_t, global_t, mesh_t
  use hash
  use hilbert
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  character(len=*)::type
  ! Local variables
  integer::idim,ilevel,icpu,igrid,ibound,ilevelmin
  integer(kind=8)::max_key
  real(kind=8)::dx
  integer(kind=8)::ngrid_tot,ikey
  integer(kind=4)::ngrid,nremain
  integer(kind=8),dimension(1:nhilbert)::hk
  integer(kind=8),dimension(1:ndim)::ix

  ! Store size in mesh object
  if(type=='amr')then
     m%ngridmax=r%ngridmax
     m%ncachemax=MAX(r%ncachemax,10000)
  endif
  if(type=='mg')then
     m%ngridmax=r%ngridmax/7
     m%ncachemax=MAX(r%ncachemax/7,10000)
  endif

  ! Allocate main oct array
  allocate(m%grid(1:m%ngridmax+m%ncachemax))
  do igrid=1,m%ngridmax+m%ncachemax
     m%grid(igrid)%lev=0
  end do

  ! Allocate grid arrays
  allocate(m%flag1(1:twotondim,1:m%ngridmax+m%ncachemax))
  allocate(m%flag2(1:twotondim,1:m%ngridmax+m%ncachemax))
  m%flag1=0
  m%flag2=0

  ! Allocate AMR specific arrays
  if(type=='amr')then
#ifdef HYDRO
     allocate(m%uold(1:twotondim,1:nvar,1:m%ngridmax+m%ncachemax))
     allocate(m%unew(1:twotondim,1:nvar,1:m%ngridmax+m%ncachemax))
     m%uold=0d0
     m%unew=0d0
#endif
#ifdef MHD
     allocate(m%bold(1:twotondim,1:6,1:m%ngridmax+m%ncachemax))
     allocate(m%bnew(1:twotondim,1:6,1:m%ngridmax+m%ncachemax))
#endif
#ifdef DO_RT
     allocate(m%rtuold(1:twotondim,1:nrtvar,1:m%ngridmax+m%ncachemax))
     allocate(m%rtunew(1:twotondim,1:nrtvar,1:m%ngridmax+m%ncachemax))
     allocate(m%emissivity(1:twotondim,1:nrtgrp,1:m%ngridmax+m%ncachemax))
#endif
#ifdef DO_CR
     allocate(m%cruold(1:twotondim,1:ncruvar,1:m%ngridmax+m%ncachemax))
     allocate(m%crunew(1:twotondim,1:ncruvar,1:m%ngridmax+m%ncachemax))
#endif
#ifdef TURB
     allocate(m%fturb(1:twotondim,1:3,1:m%ngridmax+m%ncachemax))
#endif
#ifdef GRAV
     allocate(m%rho(1:twotondim,1:m%ngridmax+m%ncachemax))
     allocate(m%phi(1:twotondim,1:m%ngridmax+m%ncachemax))
     allocate(m%nref(1:twotondim,1:m%ngridmax+m%ncachemax))
     allocate(m%f(1:twotondim,1:3,1:m%ngridmax+m%ncachemax))
     allocate(m%phi_old(1:twotondim,1:m%ngridmax+m%ncachemax))
     m%f=0d0
     m%rho=0d0
     m%phi=0d0
     m%nref=0d0
     m%phi_old=0d0
#endif
  endif

  ! Allocate MG solver specific arrays
#ifdef GRAV
  if(type=='mg')then
     allocate(m%phi(1:twotondim,1:m%ngridmax+m%ncachemax))
     allocate(m%f(1:twotondim,1:3,1:m%ngridmax+m%ncachemax))
  endif
#endif

  ! Allocate cache-related arrays
  allocate(m%dirty(1:m%ncachemax))
  allocate(m%locked(1:m%ncachemax))
  allocate(m%occupied(1:m%ncachemax))
  allocate(m%parent_cpu(1:m%ncachemax))
  allocate(m%ghost_parent_grid(1:m%ncachemax))
  allocate(m%ghost_parent_cell(1:m%ncachemax))
  m%dirty=.false.
  m%locked=.false.
  m%occupied=.false.
  m%ghost_parent_grid=0
  m%ghost_parent_cell=0
  m%free_cache=1; m%ncache=0; m%nlocked=0; m%nlocked_max=0

  allocate(m%lev_null(1:m%ncachemax))
  allocate(m%ckey_null(1:ndim,1:m%ncachemax))
  allocate(m%occupied_null(1:m%ncachemax))
  m%occupied_null=.false.
  m%free_null=1; m%nnull=0

  ! Allocate hash table for AMR data
  if(r%verbose.and.g%myid==1)write(*,*)'Initialize empty hash'
  call init_empty_hash(m%grid_dict,2*(m%ngridmax+m%ncachemax),'simple')

  ! Set initial cpu boundaries
  ! Set maximum Cartesian key per level
  if(r%verbose.and.g%myid==1)write(*,*)'Initialize level cpu boundaries'

  allocate(m%ckey_max(1:r%nlevelmax+1))
  allocate(m%hkey_max(1:nhilbert,1:r%nlevelmax+1))
  m%hkey_max=0

  allocate(m%box_ckey_min(1:3,1:r%nlevelmax+1))
  allocate(m%box_ckey_max(1:3,1:r%nlevelmax+1))

  if(r%nbound>0)then
     allocate(m%bound_ckey_min(1:3,1:r%nbound,1:r%nlevelmax+1))
     allocate(m%bound_ckey_max(1:3,1:r%nbound,1:r%nlevelmax+1))
  endif

  allocate(m%domain(1:r%nlevelmax+1))
  do ilevel=1,r%nlevelmax+1
     m%domain(ilevel)%ncpu=0
     call m%domain(ilevel)%create(g%myid,g%ncpu,g%ncpu)
  end do

  ! Make sure that the coarsest level uses only one Hilbert integer
  if(r%levelmin>(levels_per_key(ndim)+1))then
     write(*,*)'levelmin is way too large !'
     write(*,*)'are you crazy ?'
     stop
  endif

  ! Set max. Cartesian and Hilbert keys for coarse levels
  do ilevel=1,r%levelmin
     m%ckey_max(ilevel)=2**(ilevel-1)
     m%hkey_max(1,ilevel)=int(m%ckey_max(ilevel),kind=8)**ndim
  end do

  ! Set max. Cartesian and Hilbert keys for fine levels
  do ilevel=r%levelmin+1,r%nlevelmax+1
     m%ckey_max(ilevel)=2**(ilevel-1)
     m%hkey_max(1:nhilbert,ilevel)=refine_key(m%hkey_max(1:nhilbert,ilevel-1),ilevel)
  end do

  ! Bounding box for computational domain
  ! Examples are: box_xmin=1 box_xmax=-1
  ! Default is:   box_xmin=0 box_xmax=0
  ! This sets the min. and max. Cartesian keys of the box in each direction
  do ilevel=r%bound_levelmin,r%nlevelmax+1
     if(r%box_xmin.GE.0)then
        m%box_ckey_min(1,ilevel)=r%box_xmin*2**(ilevel-r%bound_levelmin)
     else
        m%box_ckey_min(1,ilevel)=(m%ckey_max(r%bound_levelmin)+r%box_xmin)*2**(ilevel-r%bound_levelmin)
     endif
     if(r%box_xmax.LE.0)then
        m%box_ckey_max(1,ilevel)=(m%ckey_max(r%bound_levelmin)+r%box_xmax)*2**(ilevel-r%bound_levelmin)
     else
        m%box_ckey_max(1,ilevel)=r%box_xmax*2**(ilevel-r%bound_levelmin)
     endif
#if NDIM>1
     if(r%box_ymin.GE.0)then
        m%box_ckey_min(2,ilevel)=r%box_ymin*2**(ilevel-r%bound_levelmin)
     else
        m%box_ckey_min(2,ilevel)=(m%ckey_max(r%bound_levelmin)+r%box_ymin)*2**(ilevel-r%bound_levelmin)
     endif
     if(r%box_ymax.LE.0)then
        m%box_ckey_max(2,ilevel)=(m%ckey_max(r%bound_levelmin)+r%box_ymax)*2**(ilevel-r%bound_levelmin)
     else
        m%box_ckey_max(2,ilevel)=r%box_ymax*2**(ilevel-r%bound_levelmin)
     endif
#endif
#if NDIM>2
     if(r%box_zmin.GE.0)then
        m%box_ckey_min(3,ilevel)=r%box_zmin*2**(ilevel-r%bound_levelmin)
     else
        m%box_ckey_min(3,ilevel)=(m%ckey_max(r%bound_levelmin)+r%box_zmin)*2**(ilevel-r%bound_levelmin)
     endif
     if(r%box_zmax.LE.0)then
        m%box_ckey_max(3,ilevel)=(m%ckey_max(r%bound_levelmin)+r%box_zmax)*2**(ilevel-r%bound_levelmin)
     else
        m%box_ckey_max(3,ilevel)=r%box_zmax*2**(ilevel-r%bound_levelmin)
     endif
#endif
  end do

  ! Bounding box for boundary regions
  do ibound=1,r%nbound
     do ilevel=r%bound_levelmin,r%nlevelmax+1
        if(r%bound_xmin(ibound).GE.0)then
           m%bound_ckey_min(1,ibound,ilevel)=r%bound_xmin(ibound)*2**(ilevel-r%bound_levelmin)
        else
           m%bound_ckey_min(1,ibound,ilevel)=(m%ckey_max(r%bound_levelmin)+r%bound_xmin(ibound))*2**(ilevel-r%bound_levelmin)
        endif
        if(r%bound_xmax(ibound).LE.0)then
           m%bound_ckey_max(1,ibound,ilevel)=(m%ckey_max(r%bound_levelmin)+r%bound_xmax(ibound))*2**(ilevel-r%bound_levelmin)
        else
           m%bound_ckey_max(1,ibound,ilevel)=r%bound_xmax(ibound)*2**(ilevel-r%bound_levelmin)
        endif
#if NDIM>1
        if(r%bound_ymin(ibound).GE.0)then
           m%bound_ckey_min(2,ibound,ilevel)=r%bound_ymin(ibound)*2**(ilevel-r%bound_levelmin)
        else
           m%bound_ckey_min(2,ibound,ilevel)=(m%ckey_max(r%bound_levelmin)+r%bound_ymin(ibound))*2**(ilevel-r%bound_levelmin)
        endif
        if(r%bound_ymax(ibound).LE.0)then
           m%bound_ckey_max(2,ibound,ilevel)=(m%ckey_max(r%bound_levelmin)+r%bound_ymax(ibound))*2**(ilevel-r%bound_levelmin)
        else
           m%bound_ckey_max(2,ibound,ilevel)=r%bound_ymax(ibound)*2**(ilevel-r%bound_levelmin)
        endif
#endif
#if NDIM>2
        if(r%bound_zmin(ibound).GE.0)then
           m%bound_ckey_min(3,ibound,ilevel)=r%bound_zmin(ibound)*2**(ilevel-r%bound_levelmin)
        else
           m%bound_ckey_min(3,ibound,ilevel)=(m%ckey_max(r%bound_levelmin)+r%bound_zmin(ibound))*2**(ilevel-r%bound_levelmin)
        endif
        if(r%bound_zmax(ibound).LE.0)then
           m%bound_ckey_max(3,ibound,ilevel)=(m%ckey_max(r%bound_levelmin)+r%bound_zmax(ibound))*2**(ilevel-r%bound_levelmin)
        else
           m%bound_ckey_max(3,ibound,ilevel)=r%bound_zmax(ibound)*2**(ilevel-r%bound_levelmin)
        endif
#endif
     end do
  end do

  ! Number of octs across the grid at levelmin
  m%nx = 1; m%ny = 1; m%nz = 1
  m%nx = m%box_ckey_max(1,r%levelmin)-m%box_ckey_min(1,r%levelmin)
#if NDIM>1
  m%ny = m%box_ckey_max(2,r%levelmin)-m%box_ckey_min(2,r%levelmin)
#endif
#if NDIM>2
  m%nz = m%box_ckey_max(3,r%levelmin)-m%box_ckey_min(3,r%levelmin)
#endif

  ! Size of the box in code units in the x-direction
  if(r%box_size(1)>0)then
     r%boxlen = r%box_size(1)/dble(m%nx)*m%ckey_max(r%levelmin)
  else
     r%box_size(1) = r%boxlen*dble(m%nx)/m%ckey_max(r%levelmin)
  endif
  ! Size of the box in code units in the y- and z-directions
  r%box_size(2) = r%box_size(1)/dble(m%nx)*dble(m%ny)
  r%box_size(3) = r%box_size(1)/dble(m%nx)*dble(m%nz)

  ! Coordinates of lower left corner of the box
  allocate(m%skip(1:ndim))
  dx=r%boxlen/2**(r%levelmin-1)
  do idim=1,ndim
     m%skip(idim)=m%box_ckey_min(idim,r%levelmin)*dx
  end do

  ! Set bounds for Hilbert keys for coarse levels
  do ilevel=1,r%levelmin-1
     max_key = m%hkey_max(1,ilevel)
     do icpu=1,g%ncpu-1
        m%domain(ilevel)%b(1,icpu) = (icpu*max_key)/g%ncpu
     end do
     m%domain(ilevel)%b(1,0) = 0
     m%domain(ilevel)%b(1,g%ncpu) = max_key
  end do

  ! Set bounds for Hilbert keys for levelmin
  ngrid_tot=int(m%nx,kind=8)*int(m%ny,kind=8)*int(m%nz,kind=8)
  if(ngrid_tot.EQ.m%hkey_max(1,r%levelmin))then
     max_key = m%hkey_max(1,r%levelmin)
     do icpu=1,g%ncpu-1
        m%domain(r%levelmin)%b(1,icpu) = (icpu*max_key)/g%ncpu
     end do
     m%domain(r%levelmin)%b(1,0) = 0
     m%domain(r%levelmin)%b(1,g%ncpu) = m%hkey_max(1,r%levelmin)
  else
     ngrid=int(ngrid_tot/int(g%ncpu,kind=8),kind=4)
     nremain=int(ngrid_tot-int(ngrid,kind=8)*int(g%ncpu,kind=8),kind=4)
     igrid=0
     icpu=1
     do ikey=1, m%hkey_max(1,r%levelmin)-1
        hk(1)=ikey
        ix=hilbert_reverse(hk,r%levelmin-1)
        if(ix(1).ge.m%box_ckey_min(1,r%levelmin).and.ix(1).lt.m%box_ckey_max(1,r%levelmin))then
#if NDIM>1
        if(ix(2).ge.m%box_ckey_min(2,r%levelmin).and.ix(2).lt.m%box_ckey_max(2,r%levelmin))then
#endif
#if NDIM>2
        if(ix(3).ge.m%box_ckey_min(3,r%levelmin).and.ix(3).lt.m%box_ckey_max(3,r%levelmin))then
#endif
           igrid=igrid+1
           if(icpu.LE.nremain)then
              if(igrid.EQ.ngrid+1)then
                 m%domain(r%levelmin)%b(1,icpu:g%ncpu) = hk(1)+1
                 igrid=0
                 icpu=icpu+1
              endif
           else
              if(igrid.EQ.ngrid)then
                 m%domain(r%levelmin)%b(1,icpu:g%ncpu) = hk(1)+1
                 igrid=0
                 icpu=icpu+1
              endif
           endif
#if NDIM>2
        endif
#endif
#if NDIM>1
        endif
#endif
        endif
     end do
     m%domain(r%levelmin)%b(1,0) = 0
     m%domain(r%levelmin)%b(1,g%ncpu) = m%hkey_max(1,r%levelmin)
  endif

  ! Set bounds for Hilbert keys for fine levels
  do ilevel=r%levelmin+1,r%nlevelmax+1
     do icpu=0,g%ncpu
        m%domain(ilevel)%b(1:nhilbert,icpu) = refine_key(m%domain(ilevel-1)%b(1:nhilbert,icpu),ilevel)
     end do
  end do

  ! Allocate head, tail and numbers for each level
  if(r%verbose.and.g%myid==1)write(*,*)'Initialize oct decomposition'
  allocate(m%head(1:r%nlevelmax))
  allocate(m%tail(1:r%nlevelmax))
  allocate(m%noct(1:r%nlevelmax))
  allocate(m%noct_min(1:r%nlevelmax))
  allocate(m%noct_max(1:r%nlevelmax))
  allocate(m%noct_tot(1:r%nlevelmax))

  ! Initialize level-based arrays
  m%head=1       ! Head oct in the level
  m%tail=0       ! Tail oct in the level
  m%noct=0       ! Number of oct in the level and in the cpu
  m%noct_tot=0   ! Total number of oct in the level (all cpus)
  m%noct_min=0   ! Minimum number of oct across all cpus
  m%noct_max=0   ! Maximum number of oct across all cpus
  m%noct_used=0  ! Number of oct used across all levels
  m%noct_used_tot=0  ! Total number of oct used (all cpus)

  ! Allocate head, tail, numbers and indice for clean and dirty octs at each level
!!$  allocate(m%head_clean(1:r%nlevelmax))
!!$  allocate(m%tail_clean(1:r%nlevelmax))
!!$  allocate(m%noct_clean(1:r%nlevelmax))
!!$  allocate(m%indx_clean(1:m%ngridmax))
!!$  allocate(m%head_dirty(1:r%nlevelmax))
!!$  allocate(m%tail_dirty(1:r%nlevelmax))
!!$  allocate(m%noct_dirty(1:r%nlevelmax))
!!$  allocate(m%indx_dirty(1:m%ngridmax))

end subroutine init_amr
!###############################################
!###############################################
!###############################################
!###############################################
subroutine init_params(mdl,r,g)
  use mdl_module
  use amr_parameters, ONLY: nhilbert,ndim
  use amr_commons, ONLY: run_t, global_t
  use hash
  use hilbert
  use output_amr_module, only: input_params
  implicit none
  type(mdl_t)::mdl
  type(run_t)::r
  type(global_t)::g

  character(len=5)::nchar
  character(len=80)::file_params
  integer::nfile_file,ncpu_file,levelmin_file,nlevelmax_file
  logical::file_exist

  ! Initial time step for each level
  g%dtold=0.0D0
  g%dtnew=0.0D0

  ! Read parameters from restart file
  if(r%nrestart>0)then
     call title(r%nrestart,nchar)
     file_params='backup_'//TRIM(nchar)//'/params.bin'
     inquire(file=file_params, exist=file_exist)
     if(file_exist)then
        call input_params(mdl,r,g,file_params,nfile_file,ncpu_file,levelmin_file,nlevelmax_file)
        if(g%myid==1)write(*,'(" Restarting from backup number ",I8)')r%nrestart
        if(g%myid==1)write(*,'(" Restart file has ",I8," files")')ncpu_file
     else
        if(g%myid==1)write(*,'(" Could not restart from file ",(A))')'backup_'//TRIM(nchar)
        stop
     endif
  else
  ! Read parameters from ramses output file
     if(r%filetype=='ramses')then
        file_params=TRIM(r%initfile(r%levelmin))//'/params.bin'
        inquire(file=file_params, exist=file_exist)
        if(file_exist)then
           call input_params(mdl,r,g,file_params,nfile_file,ncpu_file,levelmin_file,nlevelmax_file)
           if(g%myid==1)write(*,'(" Starting from ramses output folder ",(A))')r%initfile(r%levelmin)
           if(g%myid==1)write(*,'(" Output folder has ",I8," files")')nfile_file
        else
           if(g%myid==1)write(*,'(" Could not read folder ",(A))')r%initfile(r%levelmin)
           stop
        endif
     endif
  endif

end subroutine init_params
!###############################################
!###############################################
!###############################################
!###############################################
end module init_amr_module
