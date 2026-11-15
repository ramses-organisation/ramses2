module load_balance_module
contains
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_load_balance(pst,ilevel)
  use amr_parameters, only: nhilbert
  use ramses_commons, only: pst_t
  use hilbert
  use init_refine_basegrid_module, only: r_collect_noct
  implicit none
  type(pst_t)::pst
  integer::ilevel
  !--------------------------------------------------------------------
  ! This routine is the master procedure to load balance the AMR grid
  ! for all levels strictly larger than ilevel.
  !--------------------------------------------------------------------
  integer,dimension(1:pst%s%g%ncpu)::noct
  integer,allocatable,dimension(:)::input_array
  integer,allocatable,dimension(:)::output_array
  integer(kind=8),dimension(1:nhilbert)::zero_key=0
  integer(kind=8),dimension(1:nhilbert,0:pst%s%g%ncpu)::bound_key
  integer::ilev,icpu,input_size,output_size,dummy,adummy(1)

  associate(r=>pst%s%r,g=>pst%s%g,m=>pst%s%m,p=>pst%s%p,mdl=>pst%s%mdl)

  if(g%ncpu==1)return
  if(ilevel==r%nlevelmax)return

  if(r%verbose)write(*,111)ilevel
111 format(' Load balancing for all levels greater than ',I2)

  ! Compute the new domain decomposition
  do ilev=ilevel+1,r%nlevelmax

     if(m%noct_tot(ilev)>0)then

        ! Collect number of oct in each cpu for current level
        adummy(1) = ilev
        call r_collect_noct(pst,adummy,1,noct,g%ncpu)

        ! Compute input array
        input_size=g%ncpu+1
        allocate(input_array(1:input_size))
        input_array(1)=ilev
        input_array(2:input_size)=noct(1:g%ncpu)

        ! Allocate output array
        output_size=2*nhilbert*(g%ncpu+1)
        allocate(output_array(1:output_size))

        ! Compute and collect new Hilbert key boundaries for the new domain decomposition
        call r_collect_bound_key(pst,input_array,input_size,output_array,output_size)
        bound_key=reshape(transfer(output_array,zero_key),[nhilbert,g%ncpu+1])
        deallocate(input_array,output_array)

        ! Finalize new domain decomposition
        bound_key(1:nhilbert,0)=zero_key
        do icpu=1,g%ncpu
           if(gt_keys(bound_key(1:nhilbert,icpu-1),bound_key(1:nhilbert,icpu)))then
              bound_key(1:nhilbert,icpu)=bound_key(1:nhilbert,icpu-1)
           endif
        end do
        bound_key(1:nhilbert,g%ncpu)=m%hkey_max(1:nhilbert,ilev)

        ! Scatter new domain decomposition to all processors
        input_size=2*nhilbert*(g%ncpu+1)+1
        allocate(input_array(1:input_size))
        input_array(1)=ilev
        input_array(2:input_size)=transfer(reshape(bound_key,[nhilbert*(g%ncpu+1)]),input_array)
        call r_broadcast_bound_key(pst,input_array,input_size,adummy,0)
        deallocate(input_array)

     else

        ! Level is empty: use the same domain decomposition as the coarser level
        do icpu=0,g%ncpu
           bound_key(1:nhilbert,icpu) = refine_key(m%domain(ilev-1)%b(1:nhilbert,icpu),ilev)
        end do

        ! Scatter new domain decomposition to all processors
        input_size=2*nhilbert*(g%ncpu+1)+1
        allocate(input_array(1:input_size))
        input_array(1)=ilev
        input_array(2:input_size)=transfer(reshape(bound_key,[nhilbert*(g%ncpu+1)]),input_array)
        call r_broadcast_bound_key(pst,input_array,input_size,adummy,0)
        deallocate(input_array)

     endif

  end do
  ! End loop over finer levels

  ! Redistribute the grid across CPU according to the new domains
  call r_load_balance(pst,ilevel,1,dummy,0)

  end associate

end subroutine m_load_balance
!###############################################
!###############################################
!###############################################
!###############################################
recursive subroutine r_broadcast_bound_key(pst,input_array,input_size,output_array,output_size)
  use mdl_module
  use amr_parameters, only: nhilbert
  use ramses_commons, only: pst_t
  use mdl_parameters
  use hilbert
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array

  integer::ilevel
  integer(kind=8),dimension(1)::dummy
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_BROADCAST_BOUND_KEY,pst%iUpper+1,input_size,output_size,input_array)
     call r_broadcast_bound_key(pst%pLower,input_array,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size)
  else
     ilevel=input_array(1)
     pst%s%m%domain(ilevel)%b=reshape(transfer(input_array(2:input_size),dummy),[nhilbert,pst%s%g%ncpu+1])
  endif

end subroutine r_broadcast_bound_key
!###############################################
!###############################################
!###############################################
!###############################################
recursive subroutine r_collect_bound_key(pst,input_array,input_size,output_array,output_size)
  use mdl_module
  use amr_parameters, only: nhilbert
  use ramses_commons, only: pst_t
  use mdl_parameters
  use hilbert
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array

  integer,dimension(:),allocatable::next_output_array

  integer::ilevel
  integer,dimension(:),allocatable::noct
  integer(kind=8),dimension(1)::dummy
  integer(kind=8),dimension(:,:),allocatable::bound_key
  integer(kind=8),dimension(:,:),allocatable::next_bound_key
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_COLLECT_BOUND_KEY,pst%iUpper+1,input_size,output_size,input_array)
     call r_collect_bound_key(pst%pLower,input_array,input_size,output_array,output_size)
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
     allocate(noct(1:pst%s%g%ncpu))
     allocate(bound_key(1:nhilbert,0:pst%s%g%ncpu))
     bound_key=0
     ilevel=input_array(1)
     noct(1:pst%s%g%ncpu)=input_array(2:input_size)
     call compute_new_bound_key(pst%s%r,pst%s%g,pst%s%m,ilevel,noct,bound_key)
     output_array=transfer(reshape(bound_key,[nhilbert*(pst%s%g%ncpu+1)]),output_array)
     deallocate(bound_key)
     deallocate(noct)
  endif

end subroutine r_collect_bound_key
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine compute_new_bound_key(r,g,m,ilevel,noct,bound_key_target)
  use amr_parameters, only: ndim, twotondim, nhilbert, nsubgrid
  use amr_commons, only: run_t, global_t, mesh_t, oct
  use hilbert
  use hash
  use cache_commons
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer::ilevel
  integer,dimension(1:g%ncpu)::noct
  integer(kind=8),dimension(1:nhilbert,0:g%ncpu)::bound_key_target
  !----------------------------------------------------
  ! This routine compute the new Hilbert keys so that
  ! perfect load balancing is enforced.
  !----------------------------------------------------
  integer::idom,ioct,nblock
  integer::ileft,iright,istart
  integer(kind=8)::nleft,nright,nstart
  integer(kind=8),dimension(1:g%ncpu)::noct_cum
  integer(kind=8),dimension(1:g%ncpu)::ntarget_cum
  integer(kind=8),dimension(1:nhilbert)::one_key
  real(kind=8)::xtarget

  ! Hilbert key corresponding to one units
  one_key=0
  one_key(1)=1

  ! Compute cumulative grid counts
  noct_cum(1)=noct(1)
  do idom=2,g%ncpu
     noct_cum(idom)=noct_cum(idom-1)+noct(idom)
  end do

  ! Perfect load balancing, in whole blocks of nsubgrid**ndim octs
  nblock=nsubgrid**ndim
  xtarget=dble(noct_cum(g%ncpu)/nblock)/dble(g%ncpu)

  ! Find left and right domains
  ileft=0
  iright=-1
  bound_key_target=0
  do idom=1,g%ncpu
     ntarget_cum(idom)=int(dble(idom)*xtarget,kind=8)*nblock
     if(g%myid>1)then
        nleft=noct_cum(g%myid-1)
     else
        nleft=0
     endif
     nright=noct_cum(g%myid)
     IF(nright.GT.nleft)then
        if(ntarget_cum(idom).GT.nleft.AND.ntarget_cum(idom).LE.nright)then
           if(ileft==0)ileft=idom
           iright=MAX(idom,iright)
        endif
     endif
  end do

  ! Find corresponding Hilbert keys
  if(iright.GE.ileft)then
     if(g%myid.GT.1)then
        nstart=noct_cum(g%myid-1)
     else
        nstart=0
     endif
     istart=ileft
     do ioct=m%head(ilevel),m%tail(ilevel)
        nstart=nstart+1
        ! Several domains can share the same target
        do while(istart.LE.iright)
           if(nstart.LT.ntarget_cum(istart))exit
           bound_key_target(1:nhilbert,istart)=m%grid(ioct)%hkey(1:nhilbert)+one_key
           istart=istart+1
        end do
        if(istart.GT.iright)exit
     end do
  endif

end subroutine compute_new_bound_key
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_load_balance(pst,ilevel,input_size,output_array,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer::ilevel
  integer::output_array

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_LOAD_BALANCE,pst%iUpper+1,input_size,output_size,ilevel)
     call r_load_balance(pst%pLower,ilevel,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size)
  else
     call load_balance(pst%s,ilevel)
  endif

end subroutine r_load_balance
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine load_balance(s,ilevel)
  use amr_parameters, only: ndim, twotondim, nhilbert
  use init_refine_basegrid_module, only: check_subgrid
  use oct_commons, only: oct
  use ramses_commons, only: ramses_t
  use cr_parameters, only: ncruvar
  use hilbert
  use hash
  use cache_commons
  use cache
  use marshal, only: pack_fetch_refine, unpack_fetch_refine
  use nbors_utils
  implicit none
  type(ramses_t)::s
  integer::ilevel
  !------------------------------------------------
  ! This routine performs parallel load balancing.
  !------------------------------------------------
  integer::i
  integer::grid_cpu,ichild,idom,jdom,lastdom,domains_matched
  integer::nleft,nright,ileft,iright,istart,nstart,noverlaps
  integer::ioct,ilev

  integer::j,ibit,ibucket,inew
  integer::noct_zero,head_zero,indx_zero
  integer::skip_bit,ikey,true_level
  integer(kind=8),dimension(0:ndim)::hash_key
  integer(kind=8),dimension(1:nhilbert)::hk
  integer,dimension(:),allocatable::noct_level,head_level,indx_level
  integer,dimension(:),allocatable::swap_table,swap_tmp
  integer,dimension(0:twotondim-1)::bucket_count,bucket_offset
  integer,dimension(1:twotondim)::flag1_tmp,flag2_tmp
#ifdef HYDRO
  real(kind=8),dimension(1:twotondim,1:nvar)::uold_tmp
#endif
#ifdef MHD
  real(kind=8),dimension(1:twotondim,1:6)::bold_tmp
#endif
#ifdef DO_RT
  real(kind=8),dimension(1:twotondim,1:nrtvar)::rtuold_tmp
#endif
#ifdef DO_CR
  real(kind=8),dimension(1:twotondim,1:ncruvar)::cruold_tmp
#endif
#ifdef GRAV
  real(kind=8),dimension(1:twotondim,1:3)::f_tmp
  real(kind=8),dimension(1:twotondim)::phi_tmp,phi_old_tmp
#endif
  type(oct)::grid_tmp
  type(msg_large_realdp)::dummy_large_realdp

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

  !-----------------------------------------------------
  ! Step 1: dispatch octs and empty slots according to
  ! the new target Hilbert tick marks
  !-----------------------------------------------------
  m%ifree=m%noct_used+1
  do ilev=ilevel+1,r%nlevelmax

     call open_cache(mdl, m, pack_size=storage_size(dummy_large_realdp)/32, &
          pack=pack_fetch_refine, unpack=unpack_fetch_refine, &
          flush=pack_flush_loadbalance, combine=unpack_flush_loadbalance)

     hash_key(0)=ilev
     do ioct=m%head(ilev),m%tail(ilev)

        ! Check if grid sits outside future processor boundaries
        if (.not. m%domain(ilev)%in_rank(m%grid(ioct)%hkey)) then

           ! Determine the future processor
           hk = m%grid(ioct)%hkey(1:nhilbert)
           grid_cpu = m%domain(ilev)%get_rank(hk)

           ! If next cache line is occupied, free it.
           if(m%occupied(m%free_cache))call destage(mdl,m%ngridmax+m%free_cache)
           ! Set grid index to a virtual grid in local cache memory
           ichild=m%ngridmax+m%free_cache
           m%occupied(m%free_cache)=.true.
           m%parent_cpu(m%free_cache)=grid_cpu
           m%dirty(m%free_cache)=.true.
           m%ghost_parent_grid(m%free_cache)=0
           m%ghost_parent_cell(m%free_cache)=0
           ! Go to next free cache line
           m%free_cache=m%free_cache+1
           m%ncache=m%ncache+1
           if(m%free_cache.GT.m%ncachemax)m%free_cache=1
           if(m%ncache.GT.m%ncachemax)m%ncache=m%ncachemax

           ! Copy all data to the cache grid
           m%grid(ichild)=m%grid(ioct)
           m%flag1(:,ichild)=m%flag1(:,ioct)
           m%flag2(:,ichild)=m%flag2(:,ioct)
#ifdef HYDRO
           m%uold(:,:,ichild)=m%uold(:,:,ioct)
#endif
#ifdef MHD
           m%bold(:,:,ichild)=m%bold(:,:,ioct)
#endif
#ifdef DO_RT
           m%rtuold(:,:,ichild)=m%rtuold(:,:,ioct)
#endif
#ifdef DO_CR
           m%cruold(:,:,ichild)=m%cruold(:,:,ioct)
#endif
#ifdef GRAV
           m%f(:,:,ichild)=m%f(:,:,ioct)
           m%phi(:,ichild)=m%phi(:,ioct)
           m%phi_old(:,ichild)=m%phi_old(:,ioct)
#endif
           ! Set grid level to zero
           m%grid(ioct)%lev=0
           ! Free grid from hash table
           hash_key(1:ndim)=m%grid(ioct)%ckey(1:ndim)
           call hash_free(m%grid_dict,hash_key)

           ! Insert new cache grid in hash table
           call hash_setp(m%grid_dict,hash_key,ichild)

        endif

     end do

     call close_cache(mdl)

  end do

  !-----------------------------------------------------
  ! Step 2: sort new octs and empty slots according to
  ! their level (using counting sort algorithm).
  !-----------------------------------------------------
  allocate(noct_level(r%levelmin:r%nlevelmax))
  allocate(head_level(r%levelmin:r%nlevelmax))
  allocate(indx_level(r%levelmin:r%nlevelmax))
  ! Count number of octs per bucket
  noct_level=0
  noct_zero=0

  do ioct=m%tail(ilevel)+1,m%ifree-1
     true_level=m%grid(ioct)%lev
     if(true_level>0)then
        noct_level(true_level)=noct_level(true_level)+1
     else
        noct_zero=noct_zero+1
     end if
  end do
  head_level(ilevel+1)=m%tail(ilevel)+1
  do ilev=ilevel+2,r%nlevelmax
     head_level(ilev)=head_level(ilev-1)+noct_level(ilev-1)
  end do
  head_zero=head_level(r%nlevelmax)+noct_level(r%nlevelmax)

  ! Allocate main swap table
  if(m%ifree.GT.head_level(ilevel+1))then
  allocate(swap_table(head_level(ilevel+1):m%ifree-1))

  ! Build index permutation table
  indx_level=head_level
  indx_zero=head_zero
  do ioct=m%tail(ilevel)+1,m%ifree-1
     true_level=m%grid(ioct)%lev
     if(true_level>0)then
        swap_table(indx_level(true_level))=ioct
        indx_level(true_level)=indx_level(true_level)+1
     else
        swap_table(indx_zero)=ioct
        indx_zero=indx_zero+1
     end if
  end do

  !-----------------------------------------------------
  ! Step 3: sort octs level by level according to their
  ! Hilbert key using LSD Radix Sort algorithm.
  !-----------------------------------------------------
  ! Loop over levels
  do ilev=ilevel+1,r%nlevelmax
     if(noct_level(ilev)>0)then

        ! Allocate temporary swap table just for the level
        allocate(swap_tmp(head_level(ilev):head_level(ilev)+noct_level(ilev)-1))

        ! Loop over useful bits at that level
        do ibit=ilev,1,-1

           ! Get bit and key to read from
           skip_bit=ndim*(ilev-ibit)
           ikey = skip_bit / bits_per_int(ndim) + 1
           skip_bit = mod(skip_bit, bits_per_int(ndim))

           ! Count octs in buckets
           bucket_count=0
           do inew=head_level(ilev),head_level(ilev)+noct_level(ilev)-1
              ioct=swap_table(inew)
              ibucket=ibits(m%grid(ioct)%hkey(ikey),skip_bit,ndim)
              bucket_count(ibucket)=bucket_count(ibucket)+1
           end do

           ! Compute offsets
           bucket_offset(0)=head_level(ilev)
           do ibucket=1,twotondim-1
              bucket_offset(ibucket)=bucket_offset(ibucket-1)+bucket_count(ibucket-1)
           end do

           ! Sort according to Hilbert key
           do inew=head_level(ilev),head_level(ilev)+noct_level(ilev)-1
              ioct=swap_table(inew)
              ibucket=ibits(m%grid(ioct)%hkey(ikey),skip_bit,ndim)
              swap_tmp(bucket_offset(ibucket))=ioct
              bucket_offset(ibucket)=bucket_offset(ibucket)+1
           end do

           ! Store permutations in swap table
           do inew=head_level(ilev),head_level(ilev)+noct_level(ilev)-1
              swap_table(inew)=swap_tmp(inew)
           end do
        end do

        ! Deallocate tmp swap array
        deallocate(swap_tmp)
     endif
  end do

  !-----------------------------------------------------
  ! Step 4: Apply permutations directly in main memory
  ! Remember: swap_table(inew)=iold means:
  ! New data at position inew COMES FROM
  ! Old data at position iold.
  !-----------------------------------------------------
  ! Perform the swap
  do j=head_level(ilevel+1),m%ifree-1
     if(j.NE.swap_table(j))then
        hash_key(0)=m%grid(j)%lev
        hash_key(1:ndim)=m%grid(j)%ckey(1:ndim)
        if(m%grid(j)%lev>0)then
           call hash_free(m%grid_dict,hash_key)
        endif
        grid_tmp=m%grid(j)
        flag1_tmp=m%flag1(:,j)
        flag2_tmp=m%flag2(:,j)
#ifdef HYDRO
        uold_tmp=m%uold(:,:,j)
#endif
#ifdef MHD
        bold_tmp=m%bold(:,:,j)
#endif
#ifdef DO_RT
        rtuold_tmp=m%rtuold(:,:,j)
#endif
#ifdef DO_CR
        cruold_tmp=m%cruold(:,:,j)
#endif
#ifdef GRAV
        f_tmp=m%f(:,:,j)
        phi_tmp=m%phi(:,j)
        phi_old_tmp=m%phi_old(:,j)
#endif
        i=j
        inew=swap_table(j)
        do while(inew.NE.j)
           m%grid(i)=m%grid(inew)
           m%flag1(:,i)=m%flag1(:,inew)
           m%flag2(:,i)=m%flag2(:,inew)
#ifdef HYDRO
           m%uold(:,:,i)=m%uold(:,:,inew)
#endif
#ifdef MHD
           m%bold(:,:,i)=m%bold(:,:,inew)
#endif
#ifdef DO_RT
           m%rtuold(:,:,i)=m%rtuold(:,:,inew)
#endif
#ifdef DO_CR
           m%cruold(:,:,i)=m%cruold(:,:,inew)
#endif
#ifdef GRAV
           m%f(:,:,i)=m%f(:,:,inew)
           m%phi(:,i)=m%phi(:,inew)
           m%phi_old(:,i)=m%phi_old(:,inew)
#endif
           hash_key(0)=m%grid(inew)%lev
           hash_key(1:ndim)=m%grid(inew)%ckey(1:ndim)
           if(m%grid(inew)%lev>0)then
              call hash_free(m%grid_dict,hash_key)
              call hash_setp(m%grid_dict,hash_key,i)
           endif
           swap_table(i)=i
           i=inew
           inew=swap_table(inew)
        end do
        m%grid(i)=grid_tmp
        m%flag1(:,i)=flag1_tmp
        m%flag2(:,i)=flag2_tmp
#ifdef HYDRO
        m%uold(:,:,i)=uold_tmp
#endif
#ifdef MHD
        m%bold(:,:,i)=bold_tmp
#endif
#ifdef DO_RT
        m%rtuold(:,:,i)=rtuold_tmp
#endif
#ifdef DO_CR
        m%cruold(:,:,i)=cruold_tmp
#endif
#ifdef GRAV
        m%f(:,:,i)=f_tmp
        m%phi(:,i)=phi_tmp
        m%phi_old(:,i)=phi_old_tmp
#endif
        hash_key(0)=m%grid(i)%lev
        hash_key(1:ndim)=m%grid(i)%ckey(1:ndim)
        if(m%grid(i)%lev>0)then
           call hash_setp(m%grid_dict,hash_key,i)
        end if
        swap_table(i)=i
     endif
  end do
  deallocate(swap_table)
  endif

  !-----------------------------------------------------
  ! Step 5: Clean up final AMR structure
  !-----------------------------------------------------
  do ilev=ilevel+1,r%nlevelmax
     m%head(ilev)=head_level(ilev)
     m%tail(ilev)=head_level(ilev)+noct_level(ilev)-1
     m%noct(ilev)=noct_level(ilev)
  end do
  m%noct_used=m%tail(r%nlevelmax)
  deallocate(noct_level,head_level,indx_level)

  ! Check that each level is made of complete blocks
  do ilev=ilevel+1,r%nlevelmax
     call check_subgrid(s,ilev)
  end do

  end associate

end subroutine load_balance
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine pack_flush_loadbalance(mesh, igrid, msg_size, msg_array)
  use amr_parameters, only: ndim, twotondim
  use hydro_parameters, only: nvar
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_large_realdp
  use rt_parameters, only: nrtvar
  use cr_parameters, only: ncruvar
  type(mesh_t)::mesh
  integer::igrid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind,ivar,idim
  type(msg_large_realdp)::msg

  do ind=1,twotondim
     if(mesh%grid(igrid)%refined(ind))then
        msg%int4(ind)=1
     else
        msg%int4(ind)=0
     endif
  end do

#ifdef HYDRO
  do ind=1,twotondim
     do ivar=1,nvar
        msg%realdp_hydro(ind,ivar)=mesh%uold(ind,ivar,igrid)
     end do
  end do
#endif

#ifdef MHD
  do ind=1,twotondim
     do ivar=1,6
        msg%realdp_mhd(ind,ivar)=mesh%bold(ind,ivar,igrid)
     end do
  end do
#endif

#ifdef DO_RT
  do ind=1,twotondim
     do ivar=1,nrtvar
        msg%realdp_rt(ind,ivar)=mesh%rtuold(ind,ivar,igrid)
     end do
  end do
#endif

#ifdef DO_CR
  do ind=1,twotondim
     do ivar=1,ncruvar
        msg%realdp_cr(ind,ivar)=mesh%cruold(ind,ivar,igrid)
     end do
  end do
#endif

#ifdef GRAV
  do ind=1,twotondim
     do idim=1,ndim
        msg%realdp_poisson(ind,idim)=mesh%f(ind,idim,igrid)
     end do
     msg%realdp_poisson(ind,ndim+1)=mesh%phi(ind,igrid)
     msg%realdp_poisson(ind,ndim+2)=mesh%phi_old(ind,igrid)
  end do
#endif

  msg_array=transfer(msg,msg_array)

end subroutine pack_flush_loadbalance
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine unpack_flush_loadbalance(mesh,igrid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim,twotondim
  use hydro_parameters, only: nvar
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_large_realdp
  use rt_parameters, only: nrtvar
  use cr_parameters, only: ncruvar
  type(mesh_t)::mesh
  integer::igrid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind,ivar,idim
  type(msg_large_realdp)::msg

  mesh%grid(igrid)%lev=hash_key(0)
  mesh%grid(igrid)%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)

  do ind=1,twotondim
     if(msg%int4(ind)==1)then
        mesh%grid(igrid)%refined(ind)=.true.
     else
        mesh%grid(igrid)%refined(ind)=.false.
     endif
  enddo

#ifdef HYDRO
  do ind=1,twotondim
     do ivar=1,nvar
        mesh%uold(ind,ivar,igrid)=msg%realdp_hydro(ind,ivar)
     end do
  end do
#endif

#ifdef MHD
  do ind=1,twotondim
     do ivar=1,6
        mesh%bold(ind,ivar,igrid)=msg%realdp_mhd(ind,ivar)
     end do
  end do
#endif

#ifdef DO_RT
  do ind=1,twotondim
     do ivar=1,nrtvar
        mesh%rtuold(ind,ivar,igrid)=msg%realdp_rt(ind,ivar)
     end do
  end do
#endif

#ifdef DO_CR
  do ind=1,twotondim
     do ivar=1,ncruvar
        mesh%cruold(ind,ivar,igrid)=msg%realdp_cr(ind,ivar)
     end do
  end do
#endif

#ifdef GRAV
  do ind=1,twotondim
     do idim=1,ndim
        mesh%f(ind,idim,igrid)=msg%realdp_poisson(ind,idim)
     end do
     mesh%phi(ind,igrid)=msg%realdp_poisson(ind,ndim+1)
     mesh%phi_old(ind,igrid)=msg%realdp_poisson(ind,ndim+2)
  end do
#endif

end subroutine unpack_flush_loadbalance
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_balance_part(pst,ilevel,input_size,output_array,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer::ilevel
  integer::output_array

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_BALANCE_PART,pst%iUpper+1,input_size,output_size,ilevel)
     call r_balance_part(pst%pLower,ilevel,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size)
  else
#ifndef WITHOUTMPI
     if(pst%s%r%part)then
        call balance_part(pst%s,pst%s%p   ,ilevel)
     endif
     if(pst%s%r%star)then
        call balance_part(pst%s,pst%s%star,ilevel)
     endif
     if(pst%s%r%sink)then
        call balance_part(pst%s,pst%s%sink,ilevel)
     endif
     if(pst%s%r%core)then
        call balance_part(pst%s,pst%s%core,ilevel)
     endif
     if(pst%s%r%tree)then
        call balance_part(pst%s,pst%s%tree,ilevel)
     endif
     if(pst%s%r%trac)then
        call balance_part(pst%s,pst%s%trac,ilevel)
     endif
     if(pst%s%r%dust)then
        call balance_part(pst%s,pst%s%dust,ilevel)
     endif
#endif
  endif

end subroutine r_balance_part
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
#ifndef WITHOUTMPI
subroutine balance_part(s,p,ilevel)
  use amr_parameters, only: nhilbert, ndim, i8b, dp
  use pm_parameters, only: part_memory
  use pm_commons, only: part_t
  use ramses_commons, only: ramses_t
  use domain_m, only: domain_t
  use rho_fine_module, only: sort_hilbert
  use hilbert
  use mpi
  use mdl_module, only: mdl_abort
  implicit none
  type(ramses_t)::s
  type(part_t)::p
  integer::ilevel
  !---------------------------------------------------------------------
  ! This routine will dispatch particles across processors according to
  ! their Hilbert key and for a given domain decomposition.
  ! It assumes that particles are sorted according to their levels
  ! and within the level, according to their Hilbert key.
  ! It can be used only if routine rho has been called once before.
  !---------------------------------------------------------------------
  integer::info
  integer,dimension(MPI_STATUS_SIZE,s%g%ncpu)::statuses
  integer(kind=8),dimension(1:nhilbert)::hk_ref
  integer(kind=8),dimension(1:ndim)::ix_ref

  integer,dimension(1:ndim)::ix
  integer::i,istart,ipart,jpart,idim,grid_cpu,myid,ncpu
  integer::ilev,idom,icpu,mydom,ndom,count_loc,recv_cnt_tot,send_cnt_tot
  integer::nbuffer,countrecv,countsend,tag=101
  real(kind=8)::dx_loc

  integer,allocatable,dimension(:)::send_cnt,recv_cnt
  integer,allocatable,dimension(:)::send_oft,recv_oft,offset_cpu
  integer,dimension(s%g%ncpu)::reqsend,reqrecv

  real(kind=8),dimension(:),allocatable::x_recv_buf,x_send_buf
  integer(i8b),dimension(:),allocatable::l_recv_buf,l_send_buf
  integer,dimension(:),allocatable::i_recv_buf,i_send_buf

  integer(kind=8)::unbalance
  integer(kind=8),dimension(1:nhilbert)::diff_key
  type(domain_t),allocatable,dimension(:)::domain_part
  integer(kind=8),allocatable,dimension(:,:)::bound_key_target,bound_key_new
  integer(kind=8),allocatable,dimension(:,:)::bound_key_left,bound_key_right
  integer::iter
  integer(kind=8)::npart_lev,npart_lev_tot
  integer(kind=8),dimension(0:s%g%ncpu)::npart_cum,npart_cum_tot
  integer,dimension(1:s%g%ncpu)::npart_cpu,npart_cpu_tot
  real(kind=8)::xpart_target,xcum_target

  real(dp),dimension(1:ndim)::xp_tmp,vp_tmp,fp_tmp,jp_tmp
  real(dp)::mp_tmp,zp_tmp,tp_tmp
  integer::levelp_tmp
  integer(i8b)::idp_tmp

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

  !----------------------------------------------------
  ! Default for particle domains are grid domains
  !----------------------------------------------------
  allocate(domain_part(ilevel:r%nlevelmax+1))
  do ilev=ilevel,r%nlevelmax+1
     call domain_part(ilev)%copy(m%domain(ilev))
  enddo

  !-----------------------------------------------
  ! Determine particle domains if needed
  !-----------------------------------------------
  if(part_memory)then

     !-----------------------------
     ! Allocate temporary work space
     !-----------------------------
     ncpu=g%ncpu
     myid=g%myid
     allocate(bound_key_new(1:nhilbert,0:ncpu))
     allocate(bound_key_left(1:nhilbert,0:ncpu))
     allocate(bound_key_right(1:nhilbert,0:ncpu))
     allocate(bound_key_target(1:nhilbert,0:ncpu))

     ! Loop over levels
     do ilev=ilevel,r%nlevelmax
        dx_loc=r%boxlen/2**ilev

        ! Compute number of particles
        npart_lev=p%tailp(ilev)-p%headp(ilev)+1
        npart_lev_tot=0
        call MPI_ALLREDUCE(npart_lev,npart_lev_tot,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
        if(npart_lev_tot.EQ.0)cycle
        xpart_target=dble(npart_lev_tot)/dble(ncpu)

        npart_cpu=0
        npart_cpu(myid)=npart_lev
        call MPI_ALLREDUCE(npart_cpu,npart_cpu_tot,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
        npart_cpu=npart_cpu_tot

        npart_cum=0
        do icpu=1,ncpu
           npart_cum(icpu)=npart_cum(icpu-1)+npart_cpu(icpu)
        end do

        if(g%myid==1.and.r%verbose)write(*,*)"====================================="
        if(g%myid==1.and.r%verbose)write(*,'("Level=",I4," npart=",I12)')ilev,npart_lev_tot
        if(g%myid==1.and.r%verbose)write(*,'(17(I12,1X))')npart_cum

        !---------------------------------------------------------
        ! Sort particle according to current level Hilbert key
        !---------------------------------------------------------
        do i=p%headp(ilev),p%tailp(ilev)
           p%sortp(i)=i
        end do
        ix=0
        call sort_hilbert(r,g,m,p,p%headp(ilev),p%tailp(ilev),ix,0,1,ilev-1)

        ! Compute first guess domain decomposition
        bound_key_target(1:nhilbert,0:ncpu)=domain_part(ilev)%b(1:nhilbert,0:ncpu)
        bound_key_new=bound_key_target
        bound_key_left=0
        do icpu=0,ncpu
           bound_key_right(1:nhilbert,icpu)=m%hkey_max(1:nhilbert,ilev)
        end do

        unbalance=10
        iter=0

        !---------------------------------------------------------
        ! Find new Hilbert tick marks by dichotomy
        !---------------------------------------------------------
        do while (unbalance.GT.1)
           iter=iter+1

           ! Compute number of particles above tick marks
           npart_cum=0

           ! Loop over particles in Hilbert order
           do i=p%headp(ilev),p%tailp(ilev)
              ipart=p%sortp(i)

              ! Compute Hilbert key of particle parent grid
              ix_ref(1:ndim)=int((p%xp(ipart,1:ndim)+m%skip(1:ndim))/(2*dx_loc))
              hk_ref(1:nhilbert)=hilbert_key(ix_ref,ilev-1)

              do icpu=1,ncpu
                 if(gt_keys(bound_key_target(1:nhilbert,icpu),hk_ref(1:nhilbert)))then
                    npart_cum(icpu)=npart_cum(icpu)+1
                 end if
              enddo
           end do

           ! Compute global histogram
           call MPI_ALLREDUCE(npart_cum,npart_cum_tot,ncpu+1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
           npart_cum=npart_cum_tot

           unbalance=0
           do icpu=1,ncpu-1
              xcum_target=dble(icpu)*xpart_target
              if(npart_cum(icpu)>xcum_target)then
                 bound_key_new(1:nhilbert,icpu)=average_keys(bound_key_left(1:nhilbert,icpu),bound_key_target(1:nhilbert,icpu))
                 bound_key_right(1:nhilbert,icpu)=bound_key_target(1:nhilbert,icpu)
              else
                 bound_key_new(1:nhilbert,icpu)=average_keys(bound_key_right(1:nhilbert,icpu),bound_key_target(1:nhilbert,icpu))
                 bound_key_left(1:nhilbert,icpu)=bound_key_target(1:nhilbert,icpu)
              endif
              ! Note that the right key is always larger than the left key
              ! so the difference is always positive unless it overflows
              diff_key=difference_keys(bound_key_right(1:nhilbert,icpu),bound_key_left(1:nhilbert,icpu))
#if NHILBERT<=1
              unbalance=MAX(unbalance,ABS(diff_key(1)))
#endif
#if NHILBERT==2
              unbalance=MAX(unbalance,ABS(diff_key(1))+ABS(2*diff_key(2)))
#endif
#if NHILBERT==3
              unbalance=MAX(unbalance,ABS(diff_key(1))+ABS(2*diff_key(2))+ABS(2*diff_key(3)))
#endif
           end do

           bound_key_target=bound_key_new

        end do
        if(g%myid==1.and.r%verbose)write(*,'("iter=",I4,1X,17(I12,1X))')iter,npart_cum
        if(g%myid==1.and.r%verbose)write(*,'("iter=",I4,1X,17(I12,1X))')iter,(int(dble(icpu)*xpart_target,kind=8),icpu=0,ncpu)

        !---------------------------------------------------------
        ! Store new Hilbert tick marks after convergence
        !---------------------------------------------------------
        domain_part(ilev)%b(1:nhilbert,0:ncpu)=bound_key_target(1:nhilbert,0:ncpu)

     end do
     if(g%myid==1.and.r%verbose)write(*,*)"====================================="

     !-----------------------------
     ! Deallocate work space
     !-----------------------------
     deallocate(bound_key_target)
     deallocate(bound_key_new)
     deallocate(bound_key_left)
     deallocate(bound_key_right)

  endif

  !---------------------------------
  ! Balance particles across cpus
  !---------------------------------

  !-----------------------------
  ! Allocate work space
  !-----------------------------
  allocate(send_cnt(1:g%ncpu))
  allocate(recv_cnt(1:g%ncpu))
  allocate(recv_oft(1:g%ncpu))
  allocate(send_oft(1:g%ncpu))
  allocate(offset_cpu(1:g%ncpu))

  !-------------------------------------
  ! Compute number of particles to send
  !-------------------------------------
  send_cnt=0
  count_loc=0

  ! Loop over levels
  do ilev=ilevel,r%nlevelmax
     ix_ref=-1
     dx_loc=r%boxlen/2**ilev

     ! Loop over particles
     do ipart=p%headp(ilev),p%tailp(ilev)

        ! Determine in which cpu particle should sit.
        ix = int((p%xp(ipart,1:ndim)+m%skip(1:ndim))/(2*dx_loc))
        if(.NOT. ALL(ix.EQ.ix_ref(1:ndim)))then
           ix_ref(1:ndim)=ix(1:ndim)
           grid_cpu=g%myid

           ! Compute Hilbert key of particle parent grid
           hk_ref(1:nhilbert)=hilbert_key(ix_ref,ilev-1)

           ! Check if grid sits outside future processor boundaries
           if (.not. domain_part(ilev)%in_rank(hk_ref)) then
              ! Determine the future processor
              grid_cpu=domain_part(ilev)%get_rank(hk_ref)
           endif
        endif

        if(grid_cpu.EQ.g%myid)then
           ! Update count for local cpu
           count_loc=count_loc+1
        else
           ! Update count for remote cpu
           send_cnt(grid_cpu)=send_cnt(grid_cpu)+1
        end if

     end do
     ! End loop over particles
  end do
  ! End loop over levels

  !-------------------------------------
  ! Compute number of particles to receive
  !-------------------------------------
  call MPI_ALLTOALL(send_cnt(1),1,MPI_INTEGER,recv_cnt(1),1,MPI_INTEGER,MPI_COMM_WORLD,info)
  send_cnt_tot=SUM(send_cnt)
  recv_cnt_tot=SUM(recv_cnt)

  !-------------------------------------
  ! Compute offsets
  !-------------------------------------
  recv_oft=0; send_oft=0
  offset_cpu(1)=p%headp(ilevel)-1+count_loc
  do icpu=2,g%ncpu
     offset_cpu(icpu)=offset_cpu(icpu-1)+send_cnt(icpu-1)
     recv_oft(icpu)=recv_oft(icpu-1)+recv_cnt(icpu-1)
     send_oft(icpu)=send_oft(icpu-1)+send_cnt(icpu-1)
  end do

  !-------------------------------------
  ! Shift to the right particles to send
  !-------------------------------------
  send_cnt=0
  count_loc=0

  ! Loop over levels
  do ilev=ilevel,r%nlevelmax
     ix_ref=-1
     dx_loc=r%boxlen/2**ilev

     ! Loop over particles
     do ipart=p%headp(ilev),p%tailp(ilev)

        ! Determine in which cpu particle should sit.
        ix = int((p%xp(ipart,1:ndim)+m%skip(1:ndim))/(2*dx_loc))
        if(.NOT. ALL(ix.EQ.ix_ref(1:ndim)))then
           ix_ref(1:ndim)=ix(1:ndim)
           grid_cpu=g%myid

           ! Compute Hilbert key of particle parent grid
           hk_ref(1:nhilbert)=hilbert_key(ix_ref,ilev-1)

           ! Check if grid sits outside future processor boundaries
           if (.not. domain_part(ilev)%in_rank(hk_ref)) then
              ! Determine the future processor
              grid_cpu=domain_part(ilev)%get_rank(hk_ref)
           endif
        endif

        if(grid_cpu.EQ.g%myid)then
           ! Update count for local cpu
           count_loc=count_loc+1
           p%workp(ipart)=p%headp(ilevel)-1+count_loc
        else
           ! Update count for remote cpu
           send_cnt(grid_cpu)=send_cnt(grid_cpu)+1
           p%workp(ipart)=offset_cpu(grid_cpu)+send_cnt(grid_cpu)
        end if

     end do
     ! End loop over particles
  end do
  ! End loop over levels

  !-------------------------------------
  ! Swap particles using new index table
  !-------------------------------------
  do ipart=p%headp(ilevel),p%tailp(r%nlevelmax)
     do while(p%workp(ipart).NE.ipart)
        ! Swap new index
        jpart=p%workp(ipart)
        p%workp(ipart)=p%workp(jpart)
        p%workp(jpart)=jpart
        ! Swap positions
        xp_tmp(1:ndim)=p%xp(ipart,1:ndim)
        p%xp(ipart,1:ndim)=p%xp(jpart,1:ndim)
        p%xp(jpart,1:ndim)=xp_tmp(1:ndim)
        ! Swap velocities
        vp_tmp(1:ndim)=p%vp(ipart,1:ndim)
        p%vp(ipart,1:ndim)=p%vp(jpart,1:ndim)
        p%vp(jpart,1:ndim)=vp_tmp(1:ndim)
        ! Swap masses
        mp_tmp=p%mp(ipart)
        p%mp(ipart)=p%mp(jpart)
        p%mp(jpart)=mp_tmp
        ! Swap potential
        if(allocated(p%phip))then
           mp_tmp=p%phip(ipart)
           p%phip(ipart)=p%phip(jpart)
           p%phip(jpart)=mp_tmp
        endif
        ! Swap metallicity
        if(allocated(p%zp))then
           zp_tmp=p%zp(ipart)
           p%zp(ipart)=p%zp(jpart)
           p%zp(jpart)=zp_tmp
        endif
        ! Swap acceleration
        if(allocated(p%fp))then
           fp_tmp(1:ndim)=p%fp(ipart,1:ndim)
           p%fp(ipart,1:ndim)=p%fp(jpart,1:ndim)
           p%fp(jpart,1:ndim)=fp_tmp(1:ndim)
        endif
        ! Swap angular momentum (spin)
        if(allocated(p%jp))then
           jp_tmp(1:ndim)=p%jp(ipart,1:ndim)
           p%jp(ipart,1:ndim)=p%jp(jpart,1:ndim)
           p%jp(jpart,1:ndim)=jp_tmp(1:ndim)
        end if
        ! Swap birth time
        if(allocated(p%tp))then
           tp_tmp=p%tp(ipart)
           p%tp(ipart)=p%tp(jpart)
           p%tp(jpart)=tp_tmp
        endif
        ! Swap merging time
        if(allocated(p%tm))then
           tp_tmp=p%tm(ipart)
           p%tm(ipart)=p%tm(jpart)
           p%tm(jpart)=tp_tmp
        endif
        ! Swap specific energy
        if(allocated(p%up))then
           tp_tmp=p%up(ipart)
           p%up(ipart)=p%up(jpart)
           p%up(jpart)=tp_tmp
        endif
        ! Swap levels
        levelp_tmp=p%levelp(ipart)
        p%levelp(ipart)=p%levelp(jpart)
        p%levelp(jpart)=levelp_tmp
        ! Swap last accretion steps
        if(allocated(p%nacc))then
           levelp_tmp=p%nacc(ipart)
           p%nacc(ipart)=p%nacc(jpart)
           p%nacc(jpart)=levelp_tmp
        endif
        ! Swap ids
        idp_tmp=p%idp(ipart)
        p%idp(ipart)=p%idp(jpart)
        p%idp(jpart)=idp_tmp
        ! Swap merging ids
        if(allocated(p%idm))then
           idp_tmp=p%idm(ipart)
           p%idm(ipart)=p%idm(jpart)
           p%idm(jpart)=idp_tmp
        endif
        ! Swap tracking ids
        if(allocated(p%idt))then
           idp_tmp=p%idt(ipart)
           p%idt(ipart)=p%idt(jpart)
           p%idt(jpart)=idp_tmp
        endif
        ! Swap size
        if(allocated(p%size))then
           mp_tmp=p%size(ipart)
           p%size(ipart)=p%size(jpart)
           p%size(jpart)=mp_tmp
        endif
        ! Swap charge
        if(allocated(p%charge))then
           mp_tmp=p%charge(ipart)
           p%charge(ipart)=p%charge(jpart)
           p%charge(jpart)=mp_tmp
        endif
     end do
  end do

  !-------------------------------------------------------------------
  ! Set new number of particles in local processor
  !-------------------------------------------------------------------
  p%npart=p%headp(ilevel)-1+count_loc+recv_cnt_tot
  p%tailp(r%nlevelmax)=p%npart

  ! Stop if this cpu would hold more particles than its arrays can store
  if(p%npart>size(p%mp))then
     write(*,'(A,I0,A,I0,A,I0,A,I0)') &
          & 'Error: cpu ',g%myid,' particle type ',p%type, &
          & ' would hold ',p%npart,' particles in arrays sized ',size(p%mp)
     write(*,*)'Increase npartmax, nstarmax or nsinkmax (a maximum per cpu)'
     call mdl_abort(mdl)
  endif

  !-------------------------------------------------------------------
  ! Swap particles positions, velocities and masses between processors
  !-------------------------------------------------------------------
  allocate(x_recv_buf(1:recv_cnt_tot))
  allocate(x_send_buf(1:send_cnt_tot))

  !-------------------------
  ! Swap positions
  !-------------------------
  do idim=1,ndim

  countrecv=0
  do icpu=1,g%ncpu
     nbuffer=recv_cnt(icpu)
     if(nbuffer>0)then
        countrecv=countrecv+1
        istart=recv_oft(icpu)+1
        call MPI_IRECV(x_recv_buf(istart),nbuffer,MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
     endif
  end do

  do i=1,send_cnt_tot
     ipart=p%headp(ilevel)-1+count_loc+i
     x_send_buf(i)=p%xp(ipart,idim)
  end do

  countsend=0
  do icpu=1,g%ncpu
     nbuffer=send_cnt(icpu)
     if(nbuffer>0) then
        countsend=countsend+1
        istart=send_oft(icpu)+1
        call MPI_ISEND(x_send_buf(istart),nbuffer,MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
     end if
  end do

  ! Wait for full completion of receives
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)

  do i=1,recv_cnt_tot
     ipart=p%headp(ilevel)-1+count_loc+i
     p%xp(ipart,idim)=x_recv_buf(i)
  end do

  ! Wait for full completion of sends
  call MPI_WAITALL(countsend,reqsend,statuses,info)

  end do

  !-------------------------
  ! Swap velocities
  !-------------------------
  do idim=1,ndim

  countrecv=0
  do icpu=1,g%ncpu
     nbuffer=recv_cnt(icpu)
     if(nbuffer>0)then
        countrecv=countrecv+1
        istart=recv_oft(icpu)+1
        call MPI_IRECV(x_recv_buf(istart),nbuffer,MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
     endif
  end do

  do i=1,send_cnt_tot
     ipart=p%headp(ilevel)-1+count_loc+i
     x_send_buf(i)=p%vp(ipart,idim)
  end do

  countsend=0
  do icpu=1,g%ncpu
     nbuffer=send_cnt(icpu)
     if(nbuffer>0) then
        countsend=countsend+1
        istart=send_oft(icpu)+1
        call MPI_ISEND(x_send_buf(istart),nbuffer,MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
     end if
  end do

  ! Wait for full completion of receives
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)

  do i=1,recv_cnt_tot
     ipart=p%headp(ilevel)-1+count_loc+i
     p%vp(ipart,idim)=x_recv_buf(i)
  end do

  ! Wait for full completion of sends
  call MPI_WAITALL(countsend,reqsend,statuses,info)

  end do

  !-------------------------
  ! Swap masses
  !-------------------------
  countrecv=0
  do icpu=1,g%ncpu
     nbuffer=recv_cnt(icpu)
     if(nbuffer>0)then
        countrecv=countrecv+1
        istart=recv_oft(icpu)+1
        call MPI_IRECV(x_recv_buf(istart),nbuffer,MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
     endif
  end do

  do i=1,send_cnt_tot
     ipart=p%headp(ilevel)-1+count_loc+i
     x_send_buf(i)=p%mp(ipart)
  end do

  countsend=0
  do icpu=1,g%ncpu
     nbuffer=send_cnt(icpu)
     if(nbuffer>0) then
        countsend=countsend+1
        istart=send_oft(icpu)+1
        call MPI_ISEND(x_send_buf(istart),nbuffer,MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
     end if
  end do

  ! Wait for full completion of receives
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)

  do i=1,recv_cnt_tot
     ipart=p%headp(ilevel)-1+count_loc+i
     p%mp(ipart)=x_recv_buf(i)
  end do

  ! Wait for full completion of sends
  call MPI_WAITALL(countsend,reqsend,statuses,info)

  !-------------------------
  ! Swap potentials
  !-------------------------
  if(allocated(p%phip))then

     countrecv=0
     do icpu=1,g%ncpu
        nbuffer=recv_cnt(icpu)
        if(nbuffer>0)then
           countrecv=countrecv+1
           istart=recv_oft(icpu)+1
           call MPI_IRECV(x_recv_buf(istart),nbuffer,MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
        endif
     end do

     do i=1,send_cnt_tot
        ipart=p%headp(ilevel)-1+count_loc+i
        x_send_buf(i)=p%phip(ipart)
     end do

     countsend=0
     do icpu=1,g%ncpu
        nbuffer=send_cnt(icpu)
        if(nbuffer>0) then
           countsend=countsend+1
           istart=send_oft(icpu)+1
           call MPI_ISEND(x_send_buf(istart),nbuffer,MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
        end if
     end do

     ! Wait for full completion of receives
     call MPI_WAITALL(countrecv,reqrecv,statuses,info)

     do i=1,recv_cnt_tot
        ipart=p%headp(ilevel)-1+count_loc+i
        p%phip(ipart)=x_recv_buf(i)
     end do

     ! Wait for full completion of sends
     call MPI_WAITALL(countsend,reqsend,statuses,info)

  endif

  !-------------------------
  ! Swap metallicities
  !-------------------------
  if(allocated(p%zp))then

     countrecv=0
     do icpu=1,g%ncpu
        nbuffer=recv_cnt(icpu)
        if(nbuffer>0)then
           countrecv=countrecv+1
           istart=recv_oft(icpu)+1
           call MPI_IRECV(x_recv_buf(istart),nbuffer,MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
        endif
     end do

     do i=1,send_cnt_tot
        ipart=p%headp(ilevel)-1+count_loc+i
        x_send_buf(i)=p%zp(ipart)
     end do

     countsend=0
     do icpu=1,g%ncpu
        nbuffer=send_cnt(icpu)
        if(nbuffer>0) then
           countsend=countsend+1
           istart=send_oft(icpu)+1
           call MPI_ISEND(x_send_buf(istart),nbuffer,MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
        end if
     end do

     ! Wait for full completion of receives
     call MPI_WAITALL(countrecv,reqrecv,statuses,info)

     do i=1,recv_cnt_tot
        ipart=p%headp(ilevel)-1+count_loc+i
        p%zp(ipart)=x_recv_buf(i)
     end do

     ! Wait for full completion of sends
     call MPI_WAITALL(countsend,reqsend,statuses,info)

  endif

  !-------------------------
  ! Swap accelerations
  !-------------------------
  if(allocated(p%fp))then

     do idim=1,ndim

     countrecv=0
     do icpu=1,g%ncpu
        nbuffer=recv_cnt(icpu)
        if(nbuffer>0)then
           countrecv=countrecv+1
           istart=recv_oft(icpu)+1
           call MPI_IRECV(x_recv_buf(istart),nbuffer,MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
        endif
     end do

     do i=1,send_cnt_tot
        ipart=p%headp(ilevel)-1+count_loc+i
        x_send_buf(i)=p%fp(ipart,idim)
     end do

     countsend=0
     do icpu=1,g%ncpu
        nbuffer=send_cnt(icpu)
        if(nbuffer>0) then
           countsend=countsend+1
           istart=send_oft(icpu)+1
           call MPI_ISEND(x_send_buf(istart),nbuffer,MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
        end if
     end do

     ! Wait for full completion of receives
     call MPI_WAITALL(countrecv,reqrecv,statuses,info)

     do i=1,recv_cnt_tot
        ipart=p%headp(ilevel)-1+count_loc+i
        p%fp(ipart,idim)=x_recv_buf(i)
     end do

     ! Wait for full completion of sends
     call MPI_WAITALL(countsend,reqsend,statuses,info)

     end do

  endif

  !-------------------------
  ! Swap Angular momenta (spins)
  !-------------------------
  if(allocated(p%jp))then
     do idim=1,ndim

        countrecv=0
        do icpu=1,g%ncpu
           nbuffer=recv_cnt(icpu)
           if(nbuffer>0)then
              countrecv=countrecv+1
              istart=recv_oft(icpu)+1
              call MPI_IRECV(x_recv_buf(istart),nbuffer,MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
           endif
        end do

        do i=1,send_cnt_tot
           ipart=p%headp(ilevel)-1+count_loc+i
           x_send_buf(i)=p%jp(ipart,idim)
        end do

        countsend=0
        do icpu=1,g%ncpu
           nbuffer=send_cnt(icpu)
           if(nbuffer>0) then
              countsend=countsend+1
              istart=send_oft(icpu)+1
              call MPI_ISEND(x_send_buf(istart),nbuffer,MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
           end if
        end do

        ! Wait for full completion of receives
        call MPI_WAITALL(countrecv,reqrecv,statuses,info)

        do i=1,recv_cnt_tot
           ipart=p%headp(ilevel)-1+count_loc+i
           p%jp(ipart,idim)=x_recv_buf(i)
        end do

        ! Wait for full completion of sends
        call MPI_WAITALL(countsend,reqsend,statuses,info)

     end do
  end if

  !-------------------------
  ! Swap birth ages
  !-------------------------
  if(allocated(p%tp))then

     countrecv=0
     do icpu=1,g%ncpu
        nbuffer=recv_cnt(icpu)
        if(nbuffer>0)then
           countrecv=countrecv+1
           istart=recv_oft(icpu)+1
           call MPI_IRECV(x_recv_buf(istart),nbuffer,MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
        endif
     end do

     do i=1,send_cnt_tot
        ipart=p%headp(ilevel)-1+count_loc+i
        x_send_buf(i)=p%tp(ipart)
     end do

     countsend=0
     do icpu=1,g%ncpu
        nbuffer=send_cnt(icpu)
        if(nbuffer>0) then
           countsend=countsend+1
           istart=send_oft(icpu)+1
           call MPI_ISEND(x_send_buf(istart),nbuffer,MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
        end if
     end do

     ! Wait for full completion of receives
     call MPI_WAITALL(countrecv,reqrecv,statuses,info)

     do i=1,recv_cnt_tot
        ipart=p%headp(ilevel)-1+count_loc+i
        p%tp(ipart)=x_recv_buf(i)
     end do

     ! Wait for full completion of sends
     call MPI_WAITALL(countsend,reqsend,statuses,info)

  endif

  !-------------------------
  ! Swap merging ages
  !-------------------------
  if(allocated(p%tm))then

     countrecv=0
     do icpu=1,g%ncpu
        nbuffer=recv_cnt(icpu)
        if(nbuffer>0)then
           countrecv=countrecv+1
           istart=recv_oft(icpu)+1
           call MPI_IRECV(x_recv_buf(istart),nbuffer,MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
        endif
     end do

     do i=1,send_cnt_tot
        ipart=p%headp(ilevel)-1+count_loc+i
        x_send_buf(i)=p%tm(ipart)
     end do

     countsend=0
     do icpu=1,g%ncpu
        nbuffer=send_cnt(icpu)
        if(nbuffer>0) then
           countsend=countsend+1
           istart=send_oft(icpu)+1
           call MPI_ISEND(x_send_buf(istart),nbuffer,MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
        end if
     end do

     ! Wait for full completion of receives
     call MPI_WAITALL(countrecv,reqrecv,statuses,info)

     do i=1,recv_cnt_tot
        ipart=p%headp(ilevel)-1+count_loc+i
        p%tm(ipart)=x_recv_buf(i)
     end do

     ! Wait for full completion of sends
     call MPI_WAITALL(countsend,reqsend,statuses,info)

  endif

  !-------------------------
  ! Swap specific energies
  !-------------------------
  if(allocated(p%up))then

     countrecv=0
     do icpu=1,g%ncpu
        nbuffer=recv_cnt(icpu)
        if(nbuffer>0)then
           countrecv=countrecv+1
           istart=recv_oft(icpu)+1
           call MPI_IRECV(x_recv_buf(istart),nbuffer,MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
        endif
     end do

     do i=1,send_cnt_tot
        ipart=p%headp(ilevel)-1+count_loc+i
        x_send_buf(i)=p%up(ipart)
     end do

     countsend=0
     do icpu=1,g%ncpu
        nbuffer=send_cnt(icpu)
        if(nbuffer>0) then
           countsend=countsend+1
           istart=send_oft(icpu)+1
           call MPI_ISEND(x_send_buf(istart),nbuffer,MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
        end if
     end do

     ! Wait for full completion of receives
     call MPI_WAITALL(countrecv,reqrecv,statuses,info)

     do i=1,recv_cnt_tot
        ipart=p%headp(ilevel)-1+count_loc+i
        p%up(ipart)=x_recv_buf(i)
     end do

     ! Wait for full completion of sends
     call MPI_WAITALL(countsend,reqsend,statuses,info)

  endif

  !-------------------------
  ! Swap sizes
  !-------------------------
  if(allocated(p%size))then

     countrecv=0
     do icpu=1,g%ncpu
        nbuffer=recv_cnt(icpu)
        if(nbuffer>0)then
           countrecv=countrecv+1
           istart=recv_oft(icpu)+1
           call MPI_IRECV(x_recv_buf(istart),nbuffer,MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
        endif
     end do

     do i=1,send_cnt_tot
        ipart=p%headp(ilevel)-1+count_loc+i
        x_send_buf(i)=p%size(ipart)
     end do

     countsend=0
     do icpu=1,g%ncpu
        nbuffer=send_cnt(icpu)
        if(nbuffer>0) then
           countsend=countsend+1
           istart=send_oft(icpu)+1
           call MPI_ISEND(x_send_buf(istart),nbuffer,MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
        end if
     end do

     ! Wait for full completion of receives
     call MPI_WAITALL(countrecv,reqrecv,statuses,info)

     do i=1,recv_cnt_tot
        ipart=p%headp(ilevel)-1+count_loc+i
        p%size(ipart)=x_recv_buf(i)
     end do

     ! Wait for full completion of sends
     call MPI_WAITALL(countsend,reqsend,statuses,info)

  endif

  !-------------------------
  ! Swap charges
  !-------------------------
  if(allocated(p%charge))then

     countrecv=0
     do icpu=1,g%ncpu
        nbuffer=recv_cnt(icpu)
        if(nbuffer>0)then
           countrecv=countrecv+1
           istart=recv_oft(icpu)+1
           call MPI_IRECV(x_recv_buf(istart),nbuffer,MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
        endif
     end do

     do i=1,send_cnt_tot
        ipart=p%headp(ilevel)-1+count_loc+i
        x_send_buf(i)=p%charge(ipart)
     end do

     countsend=0
     do icpu=1,g%ncpu
        nbuffer=send_cnt(icpu)
        if(nbuffer>0) then
           countsend=countsend+1
           istart=send_oft(icpu)+1
           call MPI_ISEND(x_send_buf(istart),nbuffer,MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
        end if
     end do

     ! Wait for full completion of receives
     call MPI_WAITALL(countrecv,reqrecv,statuses,info)

     do i=1,recv_cnt_tot
        ipart=p%headp(ilevel)-1+count_loc+i
        p%charge(ipart)=x_recv_buf(i)
     end do

     ! Wait for full completion of sends
     call MPI_WAITALL(countsend,reqsend,statuses,info)

  endif

  deallocate(x_recv_buf,x_send_buf)

  allocate(i_recv_buf(1:recv_cnt_tot))
  allocate(i_send_buf(1:send_cnt_tot))

  !-------------------------
  ! Swap levels
  !-------------------------
  countrecv=0
  do icpu=1,g%ncpu
     nbuffer=recv_cnt(icpu)
     if(nbuffer>0)then
        countrecv=countrecv+1
        istart=recv_oft(icpu)+1
        call MPI_IRECV(i_recv_buf(istart),nbuffer,MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
     endif
  end do

  do i=1,send_cnt_tot
     ipart=p%headp(ilevel)-1+count_loc+i
     i_send_buf(i)=p%levelp(ipart)
  end do

  countsend=0
  do icpu=1,g%ncpu
     nbuffer=send_cnt(icpu)
     if(nbuffer>0) then
        countsend=countsend+1
        istart=send_oft(icpu)+1
        call MPI_ISEND(i_send_buf(istart),nbuffer,MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
     end if
  end do

  ! Wait for full completion of receives
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)

  do i=1,recv_cnt_tot
     ipart=p%headp(ilevel)-1+count_loc+i
     p%levelp(ipart)=i_recv_buf(i)
  end do

  ! Wait for full completion of sends
  call MPI_WAITALL(countsend,reqsend,statuses,info)

  !-------------------------
  ! Swap last accretion steps
  !-------------------------
  if(allocated(p%nacc))then

     countrecv=0
     do icpu=1,g%ncpu
        nbuffer=recv_cnt(icpu)
        if(nbuffer>0)then
           countrecv=countrecv+1
           istart=recv_oft(icpu)+1
           call MPI_IRECV(i_recv_buf(istart),nbuffer,MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
        endif
     end do

     do i=1,send_cnt_tot
        ipart=p%headp(ilevel)-1+count_loc+i
        i_send_buf(i)=p%nacc(ipart)
     end do

     countsend=0
     do icpu=1,g%ncpu
        nbuffer=send_cnt(icpu)
        if(nbuffer>0) then
           countsend=countsend+1
           istart=send_oft(icpu)+1
           call MPI_ISEND(i_send_buf(istart),nbuffer,MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
        end if
     end do

     ! Wait for full completion of receives
     call MPI_WAITALL(countrecv,reqrecv,statuses,info)

     do i=1,recv_cnt_tot
        ipart=p%headp(ilevel)-1+count_loc+i
        p%nacc(ipart)=i_recv_buf(i)
     end do

     ! Wait for full completion of sends
     call MPI_WAITALL(countsend,reqsend,statuses,info)

  endif

  deallocate(i_send_buf,i_recv_buf)

  allocate(l_recv_buf(1:recv_cnt_tot))
  allocate(l_send_buf(1:send_cnt_tot))

  !-------------------------
  ! Swap ids
  !-------------------------
  countrecv=0
  do icpu=1,g%ncpu
     nbuffer=recv_cnt(icpu)
     if(nbuffer>0)then
        countrecv=countrecv+1
        istart=recv_oft(icpu)+1
#ifndef LONGINT
        call MPI_IRECV(l_recv_buf(istart),nbuffer,MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
#else
        call MPI_IRECV(l_recv_buf(istart),nbuffer,MPI_INTEGER8,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
#endif
     endif
  end do

  do i=1,send_cnt_tot
     ipart=p%headp(ilevel)-1+count_loc+i
     l_send_buf(i)=p%idp(ipart)
  end do

  countsend=0
  do icpu=1,g%ncpu
     nbuffer=send_cnt(icpu)
     if(nbuffer>0) then
        countsend=countsend+1
        istart=send_oft(icpu)+1
#ifndef LONGINT
        call MPI_ISEND(l_send_buf(istart),nbuffer,MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
#else
        call MPI_ISEND(l_send_buf(istart),nbuffer,MPI_INTEGER8,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
#endif
     end if
  end do

  ! Wait for full completion of receives
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)

  do i=1,recv_cnt_tot
     ipart=p%headp(ilevel)-1+count_loc+i
     p%idp(ipart)=l_recv_buf(i)
  end do

  ! Wait for full completion of sends
  call MPI_WAITALL(countsend,reqsend,statuses,info)

  !-------------------------
  ! Swap merging ids
  !-------------------------
  if(allocated(p%idm))then

     countrecv=0
     do icpu=1,g%ncpu
        nbuffer=recv_cnt(icpu)
        if(nbuffer>0)then
           countrecv=countrecv+1
           istart=recv_oft(icpu)+1
#ifndef LONGINT
           call MPI_IRECV(l_recv_buf(istart),nbuffer,MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
#else
           call MPI_IRECV(l_recv_buf(istart),nbuffer,MPI_INTEGER8,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
#endif
        endif
     end do

     do i=1,send_cnt_tot
        ipart=p%headp(ilevel)-1+count_loc+i
        l_send_buf(i)=p%idm(ipart)
     end do

     countsend=0
     do icpu=1,g%ncpu
        nbuffer=send_cnt(icpu)
        if(nbuffer>0) then
           countsend=countsend+1
           istart=send_oft(icpu)+1
#ifndef LONGINT
           call MPI_ISEND(l_send_buf(istart),nbuffer,MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
#else
           call MPI_ISEND(l_send_buf(istart),nbuffer,MPI_INTEGER8,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
#endif
        end if
     end do

     ! Wait for full completion of receives
     call MPI_WAITALL(countrecv,reqrecv,statuses,info)

     do i=1,recv_cnt_tot
        ipart=p%headp(ilevel)-1+count_loc+i
        p%idm(ipart)=l_recv_buf(i)
     end do

     ! Wait for full completion of sends
     call MPI_WAITALL(countsend,reqsend,statuses,info)

  endif

  !-------------------------
  ! Swap tracking ids
  !-------------------------
  if(allocated(p%idt))then

     countrecv=0
     do icpu=1,g%ncpu
        nbuffer=recv_cnt(icpu)
        if(nbuffer>0)then
           countrecv=countrecv+1
           istart=recv_oft(icpu)+1
#ifndef LONGINT
           call MPI_IRECV(l_recv_buf(istart),nbuffer,MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
#else
           call MPI_IRECV(l_recv_buf(istart),nbuffer,MPI_INTEGER8,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
#endif
        endif
     end do

     do i=1,send_cnt_tot
        ipart=p%headp(ilevel)-1+count_loc+i
        l_send_buf(i)=p%idt(ipart)
     end do

     countsend=0
     do icpu=1,g%ncpu
        nbuffer=send_cnt(icpu)
        if(nbuffer>0) then
           countsend=countsend+1
           istart=send_oft(icpu)+1
#ifndef LONGINT
           call MPI_ISEND(l_send_buf(istart),nbuffer,MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
#else
           call MPI_ISEND(l_send_buf(istart),nbuffer,MPI_INTEGER8,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
#endif
        end if
     end do

     ! Wait for full completion of receives
     call MPI_WAITALL(countrecv,reqrecv,statuses,info)

     do i=1,recv_cnt_tot
        ipart=p%headp(ilevel)-1+count_loc+i
        p%idt(ipart)=l_recv_buf(i)
     end do

     ! Wait for full completion of sends
     call MPI_WAITALL(countsend,reqsend,statuses,info)

  endif

  deallocate(l_recv_buf,l_send_buf)

  !----------------------------
  ! Deallocate work space
  !----------------------------
  deallocate(send_cnt)
  deallocate(recv_cnt)
  deallocate(recv_oft)
  deallocate(send_oft)
  deallocate(offset_cpu)
  do ilev=ilevel,r%nlevelmax+1
    call domain_part(ilev)%destroy
  end do
  deallocate(domain_part)

  !----------------------------------
  ! Put all particles in level ilevel
  !----------------------------------
  p%tailp(ilevel)=p%npart
  do ilev=ilevel+1,r%nlevelmax
     p%headp(ilev)=p%npart+1
     p%tailp(ilev)=p%npart
  end do

  call MPI_ALLREDUCE(p%npart,p%npart_max,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)

  end associate

end subroutine balance_part
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
#endif
end module load_balance_module
