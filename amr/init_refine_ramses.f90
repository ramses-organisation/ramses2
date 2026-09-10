module init_refine_ramses_module
contains
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_init_refine_ramses(pst)
  use amr_parameters, only: nhilbert
  use ramses_commons, only: pst_t
  use hilbert
  use params_module, only: m_broadcast_params, m_broadcast_global
  use init_refine_basegrid_module, only:r_init_refine_basegrid,r_noct_max,r_noct_min,r_noct_tot,r_noct_used_max
  use load_balance_module, only: r_broadcast_bound_key
#ifdef GRAV
  use rho_fine_module, only: m_rho_fine
#endif
  use output_amr_module, only: input_params
  implicit none
  type(pst_t)::pst
  !--------------------------------------------------------------------
  ! This routine is the master procedure to set the base grid
  ! and initialize all cell-based variables within it.
  !--------------------------------------------------------------------
  character(LEN=80)::file_params,file_amr
  character(LEN=5)::nchar,ncharcpu
  integer(kind=8),dimension(1:nhilbert)::zero_key
  integer(kind=8),dimension(1:nhilbert,0:pst%s%g%ncpu)::bound_key
  integer::i,ilevel,icpu,ilun,ipos,dummy(1)
  integer::levelmin_file,nlevelmax_file
  integer::levelmin_max,nlevelmax_min,noct_tmp
  integer::ncpu_file_dummy,ncpu_file,input_size,output_size
  integer,dimension(:),allocatable::noct_file,noct_skip
  integer,dimension(:),allocatable::input_array,output_array

  associate(r=>pst%s%r,g=>pst%s%g,m=>pst%s%m,p=>pst%s%p,mdl=>pst%s%mdl)

  if(r%verbose)write(*,*)'Building adaptive grid from output file '//TRIM(r%initfile(r%levelmin))

  ! Set local constants
  zero_key=0

  ! Read parameters from output file
  file_params=TRIM(r%initfile(r%levelmin))//'/params.bin'
  call input_params(mdl,r,g,file_params,ncpu_file,ncpu_file_dummy,levelmin_file,nlevelmax_file)
  write(*,'(" Output folder has levelmin=",I4)')levelmin_file
  write(*,'(" Output folder has levelmax=",I4)')nlevelmax_file

  ! Broadcast parameters to all CPUs.
  call m_broadcast_params(pst)

  ! Broadcast global variables to all CPUs.
  call m_broadcast_global(pst)

  if(r%verbose)write(*,*)'Broadcast completed'

  ! Compute the proper level interval
  levelmin_max=MAX(r%levelmin,levelmin_file)
  nlevelmax_min=MIN(r%nlevelmax,nlevelmax_file)

  ! Create base grids if necessary
  if(levelmin_max>r%levelmin)then

     do ilevel=r%levelmin,levelmin_max-1

        ! Set unigrid at coarser levels
        call r_init_refine_basegrid(pst,ilevel,1)

        ! Get total, min and max grid count (only in master)
        call r_noct_tot(pst,ilevel,1,m%noct_tot(ilevel),2)
        call r_noct_min(pst,ilevel,1,m%noct_min(ilevel),1)
        call r_noct_max(pst,ilevel,1,m%noct_max(ilevel),1)

     end do

     ! Get maximum used memory (only in master)
     call r_noct_used_max(pst,r%levelmin,1,m%noct_used_max,1)

  endif

  ! Allocate local variables
  allocate(noct_file(1:ncpu_file))
  allocate(noct_skip(1:ncpu_file))

  ! Loop over relevant levels
  do ilevel=levelmin_max,nlevelmax_min

     ! Count number of octs in each CPU file
     do icpu=1,ncpu_file
        call title(icpu,ncharcpu)
        file_amr=TRIM(r%initfile(r%levelmin))//'/amr.'//TRIM(ncharcpu)
        ilun=10
        noct_skip(icpu)=0
        open(unit=ilun,file=file_amr,access="stream",action="read",form='unformatted')
        do i=levelmin_file,ilevel-1
           ipos=13+4*(i-levelmin_file)
           read(ilun,POS=ipos)noct_tmp
           noct_skip(icpu)=noct_skip(icpu)+noct_tmp
        end do
        ipos=13+4*(ilevel-levelmin_file)
        read(ilun,POS=ipos)noct_file(icpu)
        close(ilun)
     end do

     ! Allocate input array
     input_size=2*ncpu_file+4
     allocate(input_array(1:input_size))
     input_array(1)=ilevel
     input_array(2)=ncpu_file
     input_array(3)=levelmin_file
     input_array(4)=nlevelmax_file
     input_array(5:ncpu_file+4)=noct_file
     input_array(ncpu_file+5:2*ncpu_file+4)=noct_skip

     ! Allocate output array
     output_size=2*nhilbert*(g%ncpu+1)
     allocate(output_array(1:output_size))

     ! Call recursive slave routine
     call r_init_refine_ramses(pst,input_array,input_size,output_array,output_size)

     bound_key=reshape(transfer(output_array,zero_key),[nhilbert,g%ncpu+1])
     deallocate(input_array,output_array)

     ! Get total, min and max grid count (only in master).
     call r_noct_tot(pst,ilevel,1,m%noct_tot(ilevel),2)
     call r_noct_min(pst,ilevel,1,m%noct_min(ilevel),1)
     call r_noct_max(pst,ilevel,1,m%noct_max(ilevel),1)

     ! Finalize new domain decomposition
     bound_key(1:nhilbert,0)=zero_key
     do icpu=1,g%ncpu
        if(gt_keys(bound_key(1:nhilbert,icpu-1),bound_key(1:nhilbert,icpu)))then
           bound_key(1:nhilbert,icpu)=bound_key(1:nhilbert,icpu-1)
        endif
     end do
     bound_key(1:nhilbert,g%ncpu)=m%hkey_max(1:nhilbert,ilevel)

     ! Scatter new domain decomposition to all processors
     input_size=2*nhilbert*(g%ncpu+1)+1
     allocate(input_array(1:input_size))
     input_array(1)=ilevel
     input_array(2:input_size)=transfer(reshape(bound_key,[nhilbert*(g%ncpu+1)]),input_array)
     call r_broadcast_bound_key(pst,input_array,input_size,dummy,0)
     deallocate(input_array)

  end do
  ! End loop over relevant levels

  ! Get maximum used memory (only in master)
  call r_noct_used_max(pst,r%levelmin,1,m%noct_used_max,1)

  ! Deallocate local variables
  deallocate(noct_file,noct_skip)

#ifdef GRAV
  ! Compute total mass density from gas and particles on the grid
  call m_rho_fine(pst,r%levelmin,0)
#endif

  end associate

end subroutine m_init_refine_ramses
!###############################################
!###############################################
!###############################################
!###############################################
recursive subroutine r_init_refine_ramses(pst,input_array,input_size,output_array,output_size)
  use mdl_module
  use amr_parameters, only: nhilbert
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array

  integer,dimension(:),allocatable::next_output_array

  integer::ilevel,ncpu_file
  integer::levelmin_file,nlevelmax_file
  integer,dimension(:),allocatable::noct_file,nskip_file
  integer(kind=8),dimension(1)::dummy
  integer(kind=8),dimension(:,:),allocatable::bound_key
  integer(kind=8),dimension(:,:),allocatable::next_bound_key
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_INIT_REFINE_RAMSES,pst%iUpper+1,input_size,output_size,input_array)
     call r_init_refine_ramses(pst%pLower,input_array,input_size,output_array,output_size)
     allocate(next_output_array(1:output_size))
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_output_array)
     allocate(bound_key(1:nhilbert,0:pst%s%g%ncpu))
     allocate(next_bound_key(1:nhilbert,0:pst%s%g%ncpu))
     bound_key=reshape(transfer(output_array,dummy),[nhilbert,pst%s%g%ncpu+1])
     next_bound_key=reshape(transfer(next_output_array,dummy),[nhilbert,pst%s%g%ncpu+1])
     bound_key=bound_key+next_bound_key
     output_array=transfer(reshape(bound_key,[nhilbert*(pst%s%g%ncpu+1)]),output_array)
     deallocate(bound_key)
     deallocate(next_bound_key)
     deallocate(next_output_array)
  else
     ilevel=input_array(1)
     ncpu_file=input_array(2)
     levelmin_file=input_array(3)
     nlevelmax_file=input_array(4)
     allocate(noct_file(1:ncpu_file))
     noct_file=input_array(5:ncpu_file+4)
     allocate(nskip_file(1:ncpu_file))
     nskip_file=input_array(ncpu_file+5:2*ncpu_file+4)
     allocate(bound_key(1:nhilbert,0:pst%s%g%ncpu))
     bound_key=0
     call init_refine_ramses(pst%s,ilevel,ncpu_file,levelmin_file,nlevelmax_file,noct_file,nskip_file,bound_key)
     output_array=transfer(reshape(bound_key,[nhilbert*(pst%s%g%ncpu+1)]),output_array)
     deallocate(bound_key)
     deallocate(noct_file,nskip_file)
  endif

end subroutine r_init_refine_ramses
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_refine_ramses(s,ilevel,ncpu_file,levelmin_file,nlevelmax_file,noct_file,nskip_file,bound_key_target)
  !--------------------------------------------------------------
  ! This routine builds from a RAMSES output file
  ! the initial AMR grid.
  !--------------------------------------------------------------
  use mdl_module, only: mdl_abort
  use amr_parameters, only: nhilbert, ndim, twotondim, nvector
  use input_hydro_condinit_module, only: input_hydro_vecpot, cons_from_prim
  use hydro_parameters, only: nvar, nprim, nener, ie
  use ramses_commons, only: ramses_t
  use hash
  use hilbert
  implicit none
  type(ramses_t)::s
  integer::ilevel,ncpu_file
  integer::levelmin_file,nlevelmax_file
  integer,dimension(1:ncpu_file)::noct_file,nskip_file
  integer(kind=8),dimension(1:nhilbert,0:s%g%ncpu)::bound_key_target

  ! Local variables
  integer::icpu,iskip_amr=0,iskip_hydro=0,iskip_grav=0,ilun
  integer::i,n,ind,istart,iend,noct_tmp,ilev,ioct,idim,i1,j1,k1
  integer::igrid,igrid_start,nleft,nright,ileft,iright
  character(LEN=80)::file_params,file_amr,file_hydro,file_grav
  character(LEN=5)::nchar,ncharcpu
  logical::clean

  integer,dimension(1:ncpu_file)::noct_cum
  integer,dimension(1:s%g%ncpu)::ntarget_cum

  integer(kind=8),dimension(1:nhilbert)::hk
  integer(kind=8),dimension(1:ndim)::ix
  integer(kind=8),dimension(0:ndim)::hash_key
  integer(kind=8),dimension(1:nhilbert,1:s%r%nlevelmax)::key_ref
  integer(kind=8),dimension(1:nhilbert)::coarse_key,one_key
  integer,dimension(1:s%r%nlevelmax)::n_same,npatch
  integer(kind=8)::ipos

  integer,dimension(1:ndim)::ckey
  integer::refined_int
  logical,dimension(1:twotondim)::refined
  real(kind=4),dimension(1:twotondim,1:nprim)::qout
#if NENER>0
  integer::irad
#endif

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

  ! Set some constants
  one_key=0
  one_key(1)=1

  ! Compute starting grid index at that level
  if(ilevel.EQ.r%levelmin)then
     igrid_start=1
!!$     m%head_clean(ilevel)=1
!!$     m%head_dirty(ilevel)=1
  else
     igrid_start=m%tail(ilevel-1)+1
!!$     m%head_clean(ilevel)=m%head_clean(ilevel-1)+m%noct_clean(ilevel-1)
!!$     m%head_dirty(ilevel)=m%head_dirty(ilevel-1)+m%noct_dirty(ilevel-1)
  endif

  ! Set grid at current level
  igrid=igrid_start-1

  !-------------------------------------------
  ! Count cumulative numbers of octs per file
  !-------------------------------------------
  noct_cum(1)=noct_file(1)
  do icpu=2,ncpu_file
     noct_cum(icpu)=noct_cum(icpu-1)+noct_file(icpu)
  end do

  !------------------------------------------
  ! New cumulative numbers of octs per cpu
  !------------------------------------------
  do icpu=1,g%ncpu
     ntarget_cum(icpu)=int(dble(icpu)*dble(noct_cum(ncpu_file))/dble(g%ncpu))
  end do

  if(g%myid>1)then
     nleft=ntarget_cum(g%myid-1)
  else
     nleft=0
  endif
  nright=ntarget_cum(g%myid)

  !-----------------------------------------------------
  ! Compute interval of file to open for current process
  !-----------------------------------------------------
  ileft=0
  iright=-1
  bound_key_target=0
  if(nright.GT.nleft)then
     do icpu=1,ncpu_file
        if(icpu>1)then
           if(noct_cum(icpu).GT.nleft.AND.noct_cum(icpu-1).LT.nright)then
              if(ileft==0)ileft=icpu
              iright=MAX(icpu,iright)
           endif
        else
           if(noct_cum(icpu).GT.nleft)then
              if(ileft==0)ileft=icpu
              iright=MAX(icpu,iright)
           endif
        endif
     end do
  endif

  !----------------------------
  ! Read octs data in files
  !----------------------------
  ! Loop over relevant files (if any)
  do icpu=ileft,iright
     if(icpu>1)then
        istart=MAX(nleft-noct_cum(icpu-1),0)+1
        iend=MIN(nright-noct_cum(icpu-1),noct_file(icpu))
     else
        istart=nleft+1
        iend=MIN(nright,noct_file(icpu))
     endif
     call title(icpu,ncharcpu)

     ! Prepare reading the AMR file
     file_amr=TRIM(r%initfile(r%levelmin))//'/amr.'//TRIM(ncharcpu)
     open(unit=10,file=file_amr,access="stream",action="read",form='unformatted')
     iskip_amr=13+4*(nlevelmax_file-levelmin_file+1)+(4*ndim+4)*nskip_file(icpu)

     ! Prepare reading the HYDRO file
     if(r%hydro)then
        file_hydro=TRIM(r%initfile(r%levelmin))//'/hydro.'//TRIM(ncharcpu)
        open(unit=11,file=file_hydro,access="stream",action="read",form='unformatted')
        iskip_hydro=17+4*(nlevelmax_file-levelmin_file+1)+(4*twotondim*nprim)*nskip_file(icpu)
     endif

     ! Loop over useful octs in file
     do i=istart,iend

        ! Read values from AMR files
        ipos=iskip_amr+(4*ndim+4)*(i-1)
        read(10,POS=ipos)ckey
        ipos=iskip_amr+(4*ndim+4)*(i-1)+4*ndim
        read(10,POS=ipos)refined_int
        do ind=1,twotondim
           refined(ind)=btest(refined_int,ind-1)
        end do

        ! Read values from HYDRO files
        if(r%hydro)then
           ipos=iskip_hydro+(4*twotondim*nprim)*(i-1)
           read(11,POS=ipos)qout
        endif

        ! Create new oct in memory
        igrid=igrid+1
        if(igrid.GT.m%ngridmax)then
           write(*,*)'No more free memory'
           write(*,*)'Increase ngridmax'
           call mdl_abort(mdl)
        end if
        if(igrid==igrid_start)m%head(ilevel)=igrid
        m%tail(ilevel)=igrid
        m%noct(ilevel)=m%noct(ilevel)+1
        m%noct_used=m%noct_used+1

        ! Fill primitive variables from files
        ! Ignore magnetic field if present
        m%grid(igrid)%lev=ilevel
        m%grid(igrid)%ckey=ckey
        m%grid(igrid)%refined=refined
        if(r%hydro)then
#ifdef HYDRO
           do ind=1,twotondim
              m%uold(ind,1:5,igrid)=qout(ind,1:5) ! Density, velocity, pressure
#if NENER>0
              do irad=1,nener ! Non-thermal pressures
                 m%uold(ind,5+irad,igrid)=qout(ind,ie+irad)
              end do
#endif
#if NVAR>5+NENER
              do n=1,nvar-5-nener ! Passive scalar mass fraction
                 m%uold(ind,5+nener+n,igrid)=qout(ind,ie+nener+n)
              end do
#endif
           end do
#endif
        endif

        ! Set flag1 to preserve refinements
        do ind=1,twotondim
           if(m%grid(igrid)%refined(ind))then
              m%flag1(ind,igrid)=1
           else
              m%flag1(ind,igrid)=0
           endif
        end do

        ! Insert in hash table
        hash_key(0)=ilevel
        hash_key(1:ndim)=ckey
        call hash_setp(m%grid_dict,hash_key,igrid)

        ! Compute Hilbert keys of new octs
        ix(1:ndim)=ckey(1:ndim)
        hk(1:nhilbert)=hilbert_key(ix,ilevel-1)
        m%grid(igrid)%hkey(1:nhilbert)=hk(1:nhilbert)
        bound_key_target(1:nhilbert,g%myid)=hk(1:nhilbert)+one_key
     end do
     close(10)
     if(r%hydro)then
        close(11)
     endif
  end do

  if(r%hydro)then
     ! Compute initial magnetic field
     call input_hydro_vecpot(r,g,m,ilevel)
     ! Convert primitive to conservative
     call cons_from_prim(r,g,m,ilevel)
  endif

  !-----------
  ! Super-octs
  !-----------
  do ilev=1,ilevel
     npatch(ilev)=twotondim**ilev
  end do
  ilev=ilevel
  n_same=0
  key_ref=0
  key_ref(1,1:r%nlevelmax)=-1
  do ioct=m%head(ilev),m%tail(ilev)
     m%grid(ioct)%superoct=1
     coarse_key(1:nhilbert)=m%grid(ioct)%hkey(1:nhilbert)
     do i=1,MIN(ilev-1,r%nsuperoct)
        coarse_key(1:nhilbert)=coarsen_key(coarse_key(1:nhilbert),ilev-1) ! ilev-1 used to speed up only
        if(eq_keys(coarse_key(1:nhilbert),key_ref(1:nhilbert,i)))then
           n_same(i)=n_same(i)+1
        else
           n_same(i)=1
           key_ref(1:nhilbert,i)=coarse_key(1:nhilbert)
        endif
        if(n_same(i).EQ.npatch(i))then
           m%grid(ioct-npatch(i)+1:ioct)%superoct=npatch(i)
        endif
     end do
  end do

!!$  !---------------------
!!$  ! Clean and dirty octs
!!$  !---------------------
!!$  m%noct_clean(ilevel)=0
!!$  m%noct_dirty(ilevel)=0
!!$  hash_key(0)=ilevel
!!$  do ioct=m%head(ilevel),m%tail(ilevel)
!!$     clean=.true.
!!$#if NDIM>2
!!$     do k1=-1,1
!!$     hash_key(3)=m%grid(ioct)%ckey(3)+k1
!!$#endif
!!$#if NDIM>1
!!$     do j1=-1,1
!!$     hash_key(2)=m%grid(ioct)%ckey(2)+j1
!!$#endif
!!$     do i1=-1,1
!!$        hash_key(1)=m%grid(ioct)%ckey(1)+i1
!!$        clean=clean.and.hash_is_clean(m%grid_dict,hash_key)
!!$     end do
!!$#if NDIM>1
!!$     end do
!!$#endif
!!$#if NDIM>2
!!$     end do
!!$#endif
!!$     if(clean)then
!!$        m%indx_clean(m%head_clean(ilevel)+m%noct_clean(ilevel))=ioct
!!$        m%noct_clean(ilevel)=m%noct_clean(ilevel)+1
!!$     else
!!$        m%indx_dirty(m%head_dirty(ilevel)+m%noct_dirty(ilevel))=ioct
!!$        m%noct_dirty(ilevel)=m%noct_dirty(ilevel)+1
!!$     endif
!!$  end do

  end associate

end subroutine init_refine_ramses
!################################################################
!################################################################
!################################################################
!################################################################
end module init_refine_ramses_module
