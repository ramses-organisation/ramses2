module cache
  use amr_parameters, only: ndim, twotondim, nhilbert
  use amr_commons, only: mesh_t
  use clfind_commons, only: clump_t
  use cache_commons
  use hash
  use hilbert
  use mdl_module
#ifndef WITHOUTMPI
  use mpi
#endif

  type cache_oct
     integer(kind=8),dimension(0:ndim)::key
     integer::igrid
  end type cache_oct

contains
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine get_tile(m,igrid,tile,ntile)
  implicit none
  type(mesh_t)::m
  integer,intent(in)::igrid
  type(cache_oct),dimension(1:ntilemax),intent(out)::tile
  integer,intent(out)::ntile

  integer::ilevel,itile,ipos,i

  ilevel=m%grid(igrid)%lev
  if(igrid.GE.m%head(ilevel).and.igrid.LE.m%tail(ilevel))then
     itile=(igrid-m%head(ilevel))/ntilemax
     ntile=MIN(m%tail(ilevel)-itile*ntilemax-m%head(ilevel)+1,ntilemax)
     do i=1,ntile
        ipos=m%head(ilevel)+itile*ntilemax+i-1
        tile(i)%key(0)=ilevel
        tile(i)%key(1:ndim)=m%grid(ipos)%ckey(1:ndim)
        tile(i)%igrid=ipos
     end do
  else
     ntile=1
     tile(1)%key(0)=ilevel
     tile(1)%key(1:ndim)=m%grid(igrid)%ckey(1:ndim)
     tile(1)%igrid=igrid
  endif

end subroutine get_tile
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine close_cache(mdl)
  implicit none
  type(mdl_t)::mdl
  !
  ! This routine closes all cache operations.
  ! It purges all remaining flush messages.
  !
  integer::info,icache,igrid,icpu,ibuf,iskip,ipeak
  integer::send_flush_id,send_flush_id_clump,nflush
  integer::dummy_int,close_tag=7,close_id
  integer::i,ilevel,child,flush_guard,flush_clump_guard
  integer::local_peak_id
  integer(kind=8)::global_peak_id
  integer(kind=8),dimension(1:ndim)::ix
  integer(kind=8),dimension(1:nhilbert)::hk
  integer(kind=8),dimension(0:ndim)::hash_child
#ifndef WITHOUTMPI
  integer,dimension(MPI_STATUS_SIZE)::reply_status,request_status,flush_status
  integer,dimension(MPI_STATUS_SIZE)::reply_status_clump,request_status_clump,flush_status_clump
  logical::flush_cancelled,flush_clump_cancelled
#endif

  associate(m=>mdl%m,c=>mdl%c)

  ! EMPTY AND CLEAN THE GRID CACHE
  if(mdl%cache_opened)then
     do icache=1,m%ncache
        igrid=m%ngridmax+icache
        m%locked(icache)=.false.
        if(m%occupied(icache))call destage(mdl,igrid)
        m%occupied(icache)=.false.
        m%dirty(icache)=.false.
        m%ghost_parent_grid(icache)=0
        m%ghost_parent_cell(icache)=0
     end do
     m%free_cache=1
     m%ncache=0
     m%nlocked=0
     do icache=1,m%nnull
        if(m%occupied_null(icache))then
           hash_child(0)=m%lev_null(icache)
           hash_child(1:ndim)=m%ckey_null(1:ndim,icache)
           call hash_free(m%grid_dict,hash_child)
        endif
        m%occupied_null(icache)=.false.
     end do
     m%free_null=1
     m%nnull=0
  endif

  ! EMPTY AND CLEAN THE CLUMP CACHE
  if(mdl%cache_opened_clump)then
     do icache=1,c%ncache
        ipeak=c%npeak+icache
        c%locked(icache)=.false.
        if(c%occupied(icache))call destage_clump(mdl,ipeak)
        c%occupied(icache)=.false.
        c%dirty(icache)=.false.
     end do
     c%free_cache=1
     c%ncache=0
  endif

#ifndef WITHOUTMPI

  if(mdl%cache_opened)then
     ! COMPLETE THE LAST GRID FLUSH
     do icpu=1,mdl_threads(mdl)
        ibuf=mdl%cpu2buf_flush(icpu)
        if(ibuf>0)then
           iskip=1
           nflush=mdl%send_flush(ibuf)%array(iskip)
           if(nflush>0)then
              ! Post send
              call MPI_ISSEND(mdl%send_flush(ibuf)%array(iskip),mdl%size_flush_array,MPI_INTEGER,icpu-1,flush_tag,MPI_COMM_WORLD,send_flush_id,info)
              ! While waiting for completion, check on incoming messages and perform actions
              call check_mail(mdl,send_flush_id)
              mdl%send_flush(ibuf)%array(iskip)=0
           endif
        endif
     end do
  endif

  if(mdl%cache_opened_clump)then
     ! COMPLETE THE LAST CLUMP FLUSH
     do icpu=1,mdl_threads(mdl)
        ibuf=mdl%cpu2buf_flush_clump(icpu)
        if(ibuf>0)then
           iskip=1
           nflush=mdl%send_flush_clump(ibuf)%array(iskip)
           if(nflush>0)then
              ! Post send
              call MPI_ISSEND(mdl%send_flush_clump(ibuf)%array(iskip),mdl%size_flush_array_clump,MPI_INTEGER,icpu-1,flush_tag_clump,MPI_COMM_WORLD,send_flush_id_clump,info)
              ! While waiting for completion, check on incoming messages and perform actions
              call check_mail(mdl,send_flush_id_clump)
              mdl%send_flush_clump(ibuf)%array(iskip)=0
           endif
        endif
     end do
  endif

  ! CHECK-IN CHECK-OUT
  if(mdl_self(mdl).NE.1)then
     call MPI_ISEND(dummy_int,1,MPI_INTEGER,0,close_tag,MPI_COMM_WORLD,close_id,info)
     call check_mail(mdl,close_id)
     call MPI_IRECV(dummy_int,1,MPI_INTEGER,0,close_tag,MPI_COMM_WORLD,close_id,info)
     call check_mail(mdl,close_id)
  else
     do icpu=2,mdl_threads(mdl)
        call MPI_IRECV(dummy_int,1,MPI_INTEGER,MPI_ANY_SOURCE,close_tag,MPI_COMM_WORLD,close_id,info)
        call check_mail(mdl,close_id)
     end do
     do icpu=2,mdl_threads(mdl)
        call MPI_ISEND(dummy_int,1,MPI_INTEGER,icpu-1,close_tag,MPI_COMM_WORLD,close_id,info)
        call check_mail(mdl,close_id)
     end do
  endif

  ! Barrier to get the last flush message
  call MPI_BARRIER(MPI_COMM_WORLD,info)
  call check_mail(mdl,MPI_REQUEST_NULL)

  ! Finally CANCEL THE 2 GRID RECV
  if(mdl%cache_opened)then

     call MPI_CANCEL(mdl%request_id,info)
     call MPI_WAIT(mdl%request_id,request_status,info)

     flush_guard=0
     do
        call MPI_CANCEL(mdl%flush_id,info)
        call MPI_WAIT(mdl%flush_id,flush_status,info)
        call MPI_TEST_CANCELLED(flush_status,flush_cancelled,info)
        if(flush_cancelled)exit

        ! Cancellation failed: a flush message was matched and received.
        ! Process it before reposting a receive and trying to close again.
        if(mdl%combiner_rule.eq.COMBINER_EXIST)then
           iskip=1
           nflush=mdl%recv_flush_array(iskip)
           iskip=iskip+1

           do i=1,nflush
              ilevel=mdl%recv_flush_array(iskip)
              hash_child(0)=ilevel
              hash_child(1:ndim)=mdl%recv_flush_array(iskip+1:iskip+ndim)
              iskip=iskip+ndim+1
              child=hash_getp(m%grid_dict,hash_child)
              if(child>0)then
                 call unpack_flush%proc(m,child,mdl%size_msg_array, &
                      mdl%recv_flush_array(iskip:iskip+mdl%size_msg_array-1),hash_child)
              endif
              iskip=iskip+mdl%size_msg_array
           end do
        endif

        if(mdl%combiner_rule.eq.COMBINER_CREATE)then
           iskip=1
           nflush=mdl%recv_flush_array(iskip)
           iskip=iskip+1

           do i=1,nflush
              ilevel=mdl%recv_flush_array(iskip)
              hash_child(0)=ilevel
              hash_child(1:ndim)=mdl%recv_flush_array(iskip+1:iskip+ndim)
              iskip=iskip+ndim+1

              if(hash_getp(m%grid_dict,hash_child)<=0)then
                 ix(1:ndim)=hash_child(1:ndim)
                 hk(1:nhilbert)=hilbert_key(ix,ilevel-1)

                 child=m%ifree
                 m%ifree=m%ifree+1
                 if(m%ifree.GT.m%ngridmax)then
                    write(*,*)'No more free memory'
                    write(*,*)'Increase ngridmax'
                    call mdl_abort(mdl)
                 endif
                 m%grid(child)%hkey(1:nhilbert)=hk(1:nhilbert)
                 m%flag1(1:twotondim,child)=0
                 m%flag2(1:twotondim,child)=0
                 call hash_setp(m%grid_dict,hash_child,child)

                 call unpack_flush%proc(m,child,mdl%size_msg_array, &
                      mdl%recv_flush_array(iskip:iskip+mdl%size_msg_array-1),hash_child)
              endif

              iskip=iskip+mdl%size_msg_array
           end do
        endif

        call MPI_IRECV(mdl%recv_flush_array,mdl%size_flush_array,MPI_INTEGER, &
             MPI_ANY_SOURCE,flush_tag,MPI_COMM_WORLD,mdl%flush_id,info)

        flush_guard=flush_guard+1
        if(flush_guard>100000)then
           write(*,'(A,I0)')'CACHE_CLOSE_FLUSH_DRAIN_GUARD rank=',mdl_self(mdl)
           call mdl_abort(mdl)
        endif
     end do

     ! Test to free memory in corresponding MPI buffer
     do icpu=1,mdl_threads(mdl)
        call MPI_WAIT(mdl%reply_id(icpu),reply_status,info)
     end do

     ! Reset cpu mapping to flush and fetch buffers
     do icpu=1,mdl_threads(mdl)
        mdl%cpu2buf_fetch(icpu)=0
        mdl%cpu2buf_flush(icpu)=0
     end do
     mdl%ibuffer_fetch=0
     mdl%ibuffer_flush=0
  endif

  ! Finally CANCEL THE 2 CLUMP RECV
  if(mdl%cache_opened_clump)then

     call MPI_CANCEL(mdl%request_id_clump,info)
     call MPI_WAIT(mdl%request_id_clump,request_status_clump,info)

     flush_clump_guard=0
     do
        call MPI_CANCEL(mdl%flush_id_clump,info)
        call MPI_WAIT(mdl%flush_id_clump,flush_status_clump,info)
        call MPI_TEST_CANCELLED(flush_status_clump,flush_clump_cancelled,info)
        if(flush_clump_cancelled)exit

        ! Cancellation failed: a clump flush was received and must be applied.
        iskip=1
        nflush=mdl%recv_flush_array_clump(iskip)
        iskip=iskip+1
        do i=1,nflush
           global_peak_id=transfer(mdl%recv_flush_array_clump(iskip:iskip+1),global_peak_id)
           iskip=iskip+2
           local_peak_id=int(global_peak_id-c%npeak_cum(mdl_self(mdl)-1),kind=4)
           call unpack_flush_clump%proc(c,local_peak_id,mdl%size_msg_array_clump, &
                mdl%recv_flush_array_clump(iskip:iskip+mdl%size_msg_array_clump-1))
           iskip=iskip+mdl%size_msg_array_clump
        end do

        call MPI_IRECV(mdl%recv_flush_array_clump,mdl%size_flush_array_clump,MPI_INTEGER, &
             MPI_ANY_SOURCE,flush_tag_clump,MPI_COMM_WORLD,mdl%flush_id_clump,info)

        flush_clump_guard=flush_clump_guard+1
        if(flush_clump_guard>100000)then
           write(*,'(A,I0)')'CACHE_CLOSE_CLUMP_FLUSH_DRAIN_GUARD rank=',mdl_self(mdl)
           call mdl_abort(mdl)
        endif
     end do

     ! Test to free memory in corresponding MPI buffer
     do icpu=1,mdl_threads(mdl)
        call MPI_WAIT(mdl%reply_id_clump(icpu),reply_status_clump,info)
     end do

     ! Reset cpu mapping to flush and fetch buffers
     do icpu=1,mdl_threads(mdl)
        mdl%cpu2buf_fetch_clump(icpu)=0
        mdl%cpu2buf_flush_clump(icpu)=0
     end do
     mdl%ibuffer_fetch_clump=0
     mdl%ibuffer_flush_clump=0
  endif

  ! Barrier to prevent interference with the next cache
  call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif

  mdl%cache_opened=.false.
  nullify(mdl%m)
  mdl%cache_opened_clump=.false.
  nullify(mdl%c)

  end associate

end subroutine close_cache
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine open_cache(mdl, m, pack_size, pack, unpack, init, flush, combine, bound)
  implicit none
  type(mdl_t)::mdl
  type(mesh_t),target::m
  integer::pack_size
  procedure(cache_function_unpack),optional::unpack
  procedure(cache_function),optional::pack
  procedure(cache_function),optional::flush
  procedure(cache_function_init),optional::init
  procedure(cache_function_unpack),optional::combine
  procedure(cache_function_bound),optional::bound
  integer::info,icpu,iskip

  mdl%size_msg_array = pack_size

  mdl%m => m

  pack_fetch%proc => null()
  init_flush%proc => null()
  pack_flush%proc => null()
  init_bound%proc => null()
  unpack_fetch%proc => null()
  unpack_flush%proc => null()
  if (present(pack))      pack_fetch%proc => pack
  if (present(init))      init_flush%proc => init
  if (present(flush))     pack_flush%proc => flush
  if (present(bound))     init_bound%proc => bound
  if (present(unpack))  unpack_fetch%proc => unpack
  if (present(combine)) unpack_flush%proc => combine

  mdl%cache_opened=.true.

#ifndef WITHOUTMPI
  do icpu=1,mdl_threads(mdl)
     mdl%reply_id(icpu)=MPI_REQUEST_NULL
  end do

  mdl%mail_counter=0

  if (.not.present(init) .and. present(flush)) then
     mdl%combiner_rule = COMBINER_CREATE
  else
     mdl%combiner_rule = COMBINER_EXIST
  end if
  mdl%size_request_array=1+ndim
  mdl%size_flush_array=1+(1+ndim+mdl%size_msg_array)*nflushmax
  mdl%size_fetch_array=2+(1+ndim+mdl%size_msg_array)*ntilemax

  ! Set communication counters to zero
  do icpu=1,mdl%nbuffer_flush
     mdl%send_flush(icpu)%array(1)=0
  end do

  ! Post the first RECV for request
  call MPI_IRECV(mdl%recv_request_array,mdl%size_request_array,MPI_INTEGER,MPI_ANY_SOURCE,request_tag,MPI_COMM_WORLD,mdl%request_id,info)

  ! Post the first RECV for flush
  call MPI_IRECV(mdl%recv_flush_array,mdl%size_flush_array,MPI_INTEGER,MPI_ANY_SOURCE,flush_tag,MPI_COMM_WORLD,mdl%flush_id,info)
#endif

end subroutine open_cache
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine open_cache_clump(mdl, c, pack_size, pack, unpack, init, flush, combine)
  implicit none
  type(mdl_t)::mdl
  type(clump_t),target::c
  integer::pack_size
  procedure(cache_function_clump),optional::pack
  procedure(cache_function_unpack_clump),optional::unpack
  procedure(cache_function_clump),optional::flush
  procedure(cache_function_init_clump),optional::init
  procedure(cache_function_unpack_clump),optional::combine
  integer::info,icpu,iskip

  mdl%size_msg_array_clump = pack_size

  mdl%c => c

  pack_fetch_clump%proc => null()
  unpack_fetch_clump%proc => null()
  init_flush_clump%proc => null()
  pack_flush_clump%proc => null()
  unpack_flush_clump%proc => null()
  if (present(pack))      pack_fetch_clump%proc => pack
  if (present(unpack))  unpack_fetch_clump%proc => unpack
  if (present(init))      init_flush_clump%proc => init
  if (present(flush))     pack_flush_clump%proc => flush
  if (present(combine)) unpack_flush_clump%proc => combine

  mdl%cache_opened_clump = .true.

#ifndef WITHOUTMPI
  do icpu=1,mdl_threads(mdl)
     mdl%reply_id_clump(icpu)=MPI_REQUEST_NULL
  end do

  mdl%mail_counter=0

  mdl%size_request_array_clump=2
  mdl%size_flush_array_clump=1+(2+mdl%size_msg_array_clump)*nflushmax
  mdl%size_fetch_array_clump=2+(2+mdl%size_msg_array_clump)*ntilemax

  ! Set communication counters to zero
  do icpu=1,mdl%nbuffer_flush_clump
     mdl%send_flush_clump(icpu)%array(1)=0
  end do

  ! Post the first RECV for request
  call MPI_IRECV(mdl%recv_request_array_clump,mdl%size_request_array_clump,MPI_INTEGER,MPI_ANY_SOURCE,request_tag_clump,MPI_COMM_WORLD,mdl%request_id_clump,info)

  ! Post the first RECV for flush
  call MPI_IRECV(mdl%recv_flush_array_clump,mdl%size_flush_array_clump,MPI_INTEGER,MPI_ANY_SOURCE,flush_tag_clump,MPI_COMM_WORLD,mdl%flush_id_clump,info)
#endif

end subroutine open_cache_clump
!##############################################################
!##############################################################
!##############################################################
!##############################################################
end module cache
