module clump_merger_module
contains
!################################################################
!################################################################
!################################################################
!################################################################
subroutine get_peak(s,global_peak_id,local_peak_id,flush_cache,fetch_cache,lock)
  use amr_commons
  use ramses_commons, only: ramses_t
  use clfind_commons
  use cache_commons
  use hash
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  type(ramses_t)::s
  integer(kind=8)::global_peak_id
  integer::local_peak_id
  logical::flush_cache
  logical::fetch_cache
  logical,optional::lock
  !-----------------------------------------------------------
  ! This routine implement the cache request for a peak patch
  ! that is in distant MPI memory.
  !-----------------------------------------------------------
  logical::failed_request
#ifndef WITHOUTMPI
  integer,dimension(MPI_STATUS_SIZE)::send_request_status_clump
  integer::info
#endif
  integer::i,iskip,ntile_response,icounter,peak_cpu
  integer::send_request_id_clump,response_id_clump
  integer(kind=8)::gpid
  integer::lpid

  associate(g=>s%g,c=>s%c,m=>s%m,mdl=>s%mdl)

#ifndef WITHOUTMPI
  ! If counter is good, check on incoming messages and perform actions
  if(mdl%mail_counter==32)then
     call check_mail(mdl,MPI_REQUEST_NULL)
     mdl%mail_counter=0
  endif
  mdl%mail_counter=mdl%mail_counter+1
#endif

  ! Peak is in the local processor memory
  if(    global_peak_id > c%npeak_cum(g%myid-1) .and. &
       & global_peak_id <= c%npeak_cum(g%myid) )then

     local_peak_id=global_peak_id-c%npeak_cum(g%myid-1)
     return

  else

#ifndef WITHOUTMPI
     ! Get position in cache memory from hash table
     local_peak_id = hash_getp_simple(c%peak_dict,global_peak_id)

     if(present(lock))then
        if(lock)call lock_cache_clump(c,local_peak_id)
     endif

     ! If peak exists, all good
     if(local_peak_id > 0)return

     ! Get remote cpu to which peak belongs
     call get_global_peak_cpu(s,global_peak_id,peak_cpu)

     if(fetch_cache)then

        ! Send a request to the relevant cpu
        mdl%send_request_array_clump(1:2)=transfer(global_peak_id,local_peak_id,2)

        ! Post RECV for the expected response
        call MPI_IRECV(mdl%recv_fetch_array_clump,mdl%size_fetch_array_clump,MPI_INTEGER,peak_cpu-1,msg_tag_clump,MPI_COMM_WORLD,response_id_clump,info)

        ! Post SEND for the request
        call MPI_ISEND(mdl%send_request_array_clump,mdl%size_request_array_clump,MPI_INTEGER,peak_cpu-1,request_tag_clump,MPI_COMM_WORLD,send_request_id_clump,info)

        ! While waiting for reply, check on incoming messages and perform actions
        call check_mail(mdl,response_id_clump)

        ! Wait for ISEND completion to free memory in corresponding MPI buffer
        call MPI_WAIT(send_request_id_clump,send_request_status_clump,info)

        ! Check header for type of response
        iskip=1
        failed_request=(mdl%recv_fetch_array_clump(iskip).NE.1)
        iskip=iskip+1

        ! Number of tiles in the response buffer
        ntile_response=mdl%recv_fetch_array_clump(iskip)
        iskip=iskip+1

        ! Loop over tiles
        do i=1,ntile_response

           ! Find next available cache line
           if(c%locked(c%free_cache))then
              icounter=0
              do while(c%locked(c%free_cache))
                 c%free_cache=c%free_cache+1
                 icounter=icounter+1
                 if(c%free_cache>c%ncachemax)c%free_cache=1
                 if(icounter>c%ncachemax)then
                    write(*,*)'PE ',g%myid,'clump cache entirely locked'
                    stop
                 endif
              end do
           end if

           ! Next available peak in memory
           lpid=c%npeak+c%free_cache

           ! If cache line is occupied, free it.
           if(c%occupied(c%free_cache))call destage_clump(mdl,lpid)

           ! Get global peak id from message header
           gpid=transfer(mdl%recv_fetch_array_clump(iskip:iskip+1),gpid)
           iskip=iskip+2

           ! Insert local peak id in hash table
           call hash_setp_simple(c%peak_dict,gpid,lpid)

           c%gid(c%free_cache)=gpid
           c%parent_cpu(c%free_cache)=peak_cpu
           c%occupied(c%free_cache)=.true.
           c%dirty(c%free_cache)=.false.

           ! Set the request peak local peak id
           if(global_peak_id==gpid)local_peak_id=lpid

           ! Unpack response to fetch request
           call unpack_fetch_clump%proc(c,lpid,mdl%size_msg_array_clump,mdl%recv_fetch_array_clump(iskip:iskip+mdl%size_msg_array_clump-1))

           ! If we also have also a flush cache...
           ! This is for combined read-write cache operations
           if(flush_cache)then

              c%dirty(c%free_cache)=.true.

              ! Set initialisation rule for combiner operations
              call init_flush_clump%proc(c,lpid)

           endif

           ! Go to next free cache line
           c%free_cache=c%free_cache+1
           c%ncache=c%ncache+1
           if(c%free_cache.GT.c%ncachemax)then
              c%free_cache=1
           endif
           if(c%ncache.GT.c%ncachemax)c%ncache=c%ncachemax

           ! Go to next tile
           iskip=iskip+mdl%size_msg_array_clump

        end do

     else if(flush_cache)then

	! Find next available cache line
        if(c%locked(c%free_cache))then
           icounter=0
           do while(c%locked(c%free_cache))
              c%free_cache=c%free_cache+1
              icounter=icounter+1
              if(c%free_cache>c%ncachemax)c%free_cache=1
              if(icounter>c%ncachemax)then
                 write(*,*)'PE ',g%myid,'clump cache entirely locked'
                 stop
              endif
           end do
        end if

        ! Next available peak in memory
        local_peak_id=c%npeak+c%free_cache

        ! If cache line is occupied, free it.
        if(c%occupied(c%free_cache))call destage_clump(mdl,local_peak_id)

        ! Insert local peak id in hash table
        call hash_setp_simple(c%peak_dict,global_peak_id,local_peak_id)

        c%gid(c%free_cache)=global_peak_id
        c%parent_cpu(c%free_cache)=peak_cpu
        c%occupied(c%free_cache)=.true.
        c%dirty(c%free_cache)=.true.

        ! Set initialisation rule for combiner operations
        call init_flush_clump%proc(c,local_peak_id)

        ! Go to next free cache line
        c%free_cache=c%free_cache+1
        c%ncache=c%ncache+1
        if(c%free_cache.GT.c%ncachemax)then
           c%free_cache=1
        endif
        if(c%ncache.GT.c%ncachemax)c%ncache=c%ncachemax

     endif

     if(present(lock))then
        if(lock)call lock_cache_clump(c,local_peak_id)
     endif

#endif
  endif

  end associate

end subroutine get_peak
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine lock_cache_clump(c,local_peak_id)
  use clfind_commons, only: clump_t
  implicit none
  type(clump_t)::c
  integer::local_peak_id
  !----------------------------------------------------
  ! This routine locks a cache line because
  ! it will be updated later.
  !----------------------------------------------------
  integer::icache
  if(local_peak_id<=c%npeak)return
  icache=local_peak_id-c%npeak
  c%locked(icache)=.true.
end subroutine lock_cache_clump
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine unlock_cache_clump(c,local_peak_id)
  use clfind_commons, only: clump_t
  implicit none
  type(clump_t)::c
  integer::local_peak_id
  !----------------------------------------------------
  ! This routine unlocks a cache line because
  ! it has been updated and can be flushed.
  !----------------------------------------------------
  integer::icache
  if(local_peak_id<=c%npeak)return
  icache=local_peak_id-c%npeak
  c%locked(icache)=.false.
end subroutine unlock_cache_clump
!################################################################
!################################################################
!################################################################
!################################################################
subroutine get_global_peak_cpu(s,global_peak_id,peak_cpu)
  use amr_commons
  use clfind_commons
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t)::s
  integer(kind=8)::global_peak_id
  integer::peak_cpu
  !----------------------------------------------------
  ! This routine returns the mpi domain of a peak with
  ! input argument the peak global id
  !----------------------------------------------------
  integer::icpu
  associate(g=>s%g,c=>s%c)
  peak_cpu = g%ncpu
  do icpu = 1,g%ncpu
     if(    global_peak_id > c%npeak_cum(icpu-1) .and. &
          & global_peak_id <= c%npeak_cum(icpu))then
        peak_cpu = icpu
     endif
  end do
  end associate
end subroutine get_global_peak_cpu
!################################################################
!################################################################
!################################################################
!################################################################
subroutine allocate_peak_patch_arrays(s)
  use amr_parameters, ONLY: ndim, nbin
  use clfind_commons
  use ramses_commons, ONLY: ramses_t
  implicit none
  type(ramses_t)::s
  !----------------------------------------------------
  ! This routine allocate all arrays that are needed at
  ! different steps of the clump finder.
  !----------------------------------------------------

  associate(r=>s%r,g=>s%g,m=>s%m,c=>s%c)

  !------------------------------------------
  ! Compute the max size of peak-based arrays
  !------------------------------------------
  c%npeak_max=c%npeak+c%ncachemax

  !-------------------------------
  ! Allocate peak-patch properties
  !-------------------------------
  allocate(c%saddle_dens(1:c%npeak_max))
  allocate(c%saddle_nbor(1:c%npeak_max))

  allocate(c%lev_peak(1:c%npeak_max))
  allocate(c%new_peak(c%npeak_max))
  allocate(c%relevance(1:c%npeak_max))
  allocate(c%tidal_dens(1:c%npeak_max))

  allocate(c%peak_pos(1:c%npeak_max,1:ndim))
  allocate(c%peak_vel(1:c%npeak_max,1:ndim))
  allocate(c%peak_acc(1:c%npeak_max,1:ndim))
  allocate(c%peak_com(1:c%npeak_max,1:ndim))

  allocate(c%min_dens(1:c%npeak_max))
  allocate(c%n_cells(1:c%npeak_max))
  allocate(c%clump_mass(1:c%npeak_max))
  allocate(c%clump_vol(1:c%npeak_max))
  allocate(c%clump_rad(1:c%npeak_max))

  allocate(c%npart(1:c%npeak_max))
  allocate(c%particle_mass(1:c%npeak_max))

  allocate(c%ind_final(1:c%npeak_max))

  !-------------------------------
  ! Allocate halo-patch properties
  !-------------------------------
  if(c%saddle_threshold>0)then
     allocate(c%ind_halo(1:c%npeak_max))
     allocate(c%halo_mass(1:c%npeak_max))
     allocate(c%mass_bin(1:c%npeak_max,1:nbin))
     allocate(c%phi(1:c%npeak_max,1:nbin))
  endif

  !------------------------------------------
  ! Allocate merger tree particles properties
  !------------------------------------------
  if(r%tree)then
     allocate(c%ntree(1:c%npeak_max))
     allocate(c%form_tree(1:c%npeak_max))
     allocate(c%min_tree_id(1:c%npeak_max))
  endif

  !-----------------------------------
  ! Allocate sink particles properties
  !-----------------------------------
  if(r%sink)then
     allocate(c%nsink(1:c%npeak_max))
     allocate(c%form_sink(1:c%npeak_max))
  endif

  !--------------------
  ! Allocate hash table
  !--------------------
  call init_empty_hash_simple(c%peak_dict,c%ncachemax)

  !----------------------
  ! Allocate cache memory
  !----------------------
  allocate(c%dirty(1:c%ncachemax))
  allocate(c%locked(1:c%ncachemax))
  allocate(c%occupied(1:c%ncachemax))
  allocate(c%parent_cpu(1:c%ncachemax))
  allocate(c%gid(1:c%ncachemax))
  c%dirty=.false.
  c%locked=.false.
  c%occupied=.false.
  c%free_cache=1
  c%ncache=0

  !---------------------
  ! Allocate cache comms
  !---------------------
  call init_cache_clump(s%mdl)

  !-------------------------------------------------
  ! Allocate various particle types peak and halo id
  !-------------------------------------------------
  if(r%part)then
     allocate(s%p%pid(1:s%p%npart))
     allocate(s%p%hid(1:s%p%npart))
  endif
  if(r%star)then
     allocate(s%star%pid(1:s%star%npart))
     allocate(s%star%hid(1:s%star%npart))
  endif
  if(r%sink)then
     allocate(s%sink%pid(1:s%sink%npart))
     allocate(s%sink%hid(1:s%sink%npart))
  endif
  if(r%tree)then
     allocate(s%tree%pid(1:s%tree%npart))
     allocate(s%tree%hid(1:s%tree%npart))
  endif

  end associate

end subroutine allocate_peak_patch_arrays
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_deallocate_clump(pst)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_CLUMP_DEALLOC,pst%iUpper+1)
     call r_deallocate_clump(pst%pLower)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call deallocate_peak_patch_arrays(pst%s)
  endif

end subroutine r_deallocate_clump
!################################################################
!################################################################
!################################################################
!################################################################
subroutine deallocate_peak_patch_arrays(s)
  use clfind_commons
  use ramses_commons, only: ramses_t
  use hash, only: reset_entire_hash_simple
  implicit none
  type(ramses_t)::s
  !----------------------------------------------------
  ! This routine deallocates all arrays that habe been
  ! previously allocated by the clump finder.
  !----------------------------------------------------
  associate(r=>s%r,g=>s%g,m=>s%m,c=>s%c)

  ! Deallocate test particle arrays
  if(c%ntest>0)then
     deallocate(c%cell)
     deallocate(c%grid)
     deallocate(c%level)
     deallocate(c%hash)
     c%ntest=0
  endif
  if(c%ntest_tot==0)return

  ! Deallocate cumulative peak count CPU array
  deallocate(c%npeak_cum)
  if(c%npeak_tot==0)return

  ! Deallocate peak-patch arrays
  c%npeak_tot=0
  deallocate(c%peak_cell)
  deallocate(c%peak_grid)
  deallocate(c%peak_level)
  deallocate(c%max_dens)

  deallocate(c%saddle_dens)
  deallocate(c%saddle_nbor)

  deallocate(c%lev_peak)
  deallocate(c%new_peak)
  deallocate(c%relevance)
  deallocate(c%tidal_dens)

  deallocate(c%peak_pos)
  deallocate(c%peak_vel)
  deallocate(c%peak_acc)
  deallocate(c%peak_com)

  deallocate(c%min_dens)
  deallocate(c%n_cells)
  deallocate(c%clump_mass)
  deallocate(c%clump_vol)
  deallocate(c%clump_rad)

  deallocate(c%npart)
  deallocate(c%particle_mass)

  deallocate(c%ind_final)

  ! Deallocate halo-patch arrays
  if(c%saddle_threshold>0)then
     deallocate(c%ind_halo)
     deallocate(c%halo_mass)
     deallocate(c%mass_bin)
     deallocate(c%phi)
  endif

  ! Deallocate tree arrays
  if(r%tree)then
     deallocate(c%ntree)
     deallocate(c%form_tree)
     deallocate(c%min_tree_id)
  endif

  ! Deallocate sink arrays
  if(r%sink)then
     deallocate(c%nsink)
     deallocate(c%form_sink)
  endif

  ! Deallocate hash table
  call reset_entire_hash_simple(c%peak_dict)

  ! Dellocate cache memory
  deallocate(c%dirty)
  deallocate(c%locked)
  deallocate(c%occupied)
  deallocate(c%parent_cpu)
  deallocate(c%gid)

  ! Deallocate cache comms
  call kill_cache_clump(s%mdl)

  ! Deallocate particle pid and hid
  if(r%part)then
     deallocate(s%p%pid)
     deallocate(s%p%hid)
  endif
  if(r%star)then
     deallocate(s%star%pid)
     deallocate(s%star%hid)
  endif
  if(r%sink)then
     deallocate(s%sink%pid)
     deallocate(s%sink%hid)
  endif
  if(r%tree)then
     deallocate(s%tree%pid)
     deallocate(s%tree%hid)
  endif

  end associate

end subroutine deallocate_peak_patch_arrays
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine collect_saddle(s)
  use amr_parameters, only: twotondim, ndim
  use ramses_commons, only: ramses_t
  use cache_commons
  use cache
  use nbors_utils
  implicit none
  type(ramses_t)::s
  !-------------------------------------------------------------------
  ! This is the clump finder routine for collecting the densest saddle
  ! points between each patch and its corresponding neighbor.
  ! Written by Ziyong Wu (mini-ramses version December 2023).
  !-------------------------------------------------------------------
  type(msg_int4_small_realdp)::dummy_int4_small_realdp
  type(msg_saddle_clump)::dummy_saddle_clump
  integer::ilevel,icpu,next_level,now_level,icelln,igridn,idim,j,jpeak,k
  integer::ipart,jpart,ip,i,ipeak,itest,igrid,ind,peak_cen,peak_nbor
  integer(kind=8),dimension(1:s%g%ncpu)::npeak_cpu,npeak_cpu_all
  integer,dimension(1:ndim)::ckey,ckey_nbor
  integer(kind=8),dimension(0:ndim)::hash_cell,hash_nbor
  real(kind=8)::dens_cen,dens_ave,dens_nbor,x,y,z
  real(kind=8),dimension(1:ndim)::xcen,xnei
  integer,parameter::nSnei=56
  real(kind=8),dimension(1:ndim,1:nSnei)::xSnei
  integer,dimension(1:nSnei)::igrid_nbor
  integer,dimension(1:nSnei)::icell_nbor
  integer(kind=8)::global_peak_id
  logical::ok

  associate(r=>s%r,g=>s%g,m=>s%m,c=>s%c,mdl=>s%mdl)

  !--------------------------------------------------------
  ! Arrays to define neighbors (center=[0,0,0])
  ! normalized to dx = 1 = size of the central leaf cell
  ! from -0.75 to 0.75
  !--------------------------------------------------------
  ind=0
  do k=1,4
     do j=1,4
        do i=1,4
           ok=.true.
           if((i==2.or.i==3).and.(j==2.or.j==3).and.(k==2.or.k==3)) ok=.false. ! centre
           if(ok)then
              ind = ind+1
              x = (i-1)+0.5d0 - 2
              y = (j-1)+0.5d0 - 2
              z = (k-1)+0.5d0 - 2
              xSnei(1,ind) = x/2d0
#if NDIM>1
              xSnei(2,ind) = y/2d0
#endif
#if NDIM>2
              xSnei(3,ind) = z/2d0
#endif
           endif
        enddo
     enddo
  enddo

  !----------------------------------------------------
  ! Compute densest saddle point and associated peak id
  !----------------------------------------------------
  call open_cache(mdl, m, pack_size=storage_size(dummy_int4_small_realdp)/32, &
       pack=pack_fetch_saddle, unpack=unpack_fetch_saddle)

  call open_cache_clump(mdl, c, pack_size=storage_size(dummy_saddle_clump)/32, &
       init=init_flush_saddle,flush=pack_flush_saddle,combine=unpack_flush_saddle)

  c%saddle_dens = c%density_threshold
  c%saddle_nbor = 0

  do itest=1,c%ntest
     ilevel=c%level(itest)
     igrid=c%grid(itest)
     ind=c%cell(itest)

     peak_cen = m%flag1(ind,igrid)
#ifdef GRAV
     dens_cen = m%rho(ind,igrid)
#endif
     ! Set grid index to null
     icelln=0
     igridn=0
     do j=1,nSnei
        igrid_nbor(j)=0
     end do

     xcen(1)=2*m%grid(igrid)%ckey(1)+MOD((ind-1)  ,2)+0.5
#if NDIM>1
     xcen(2)=2*m%grid(igrid)%ckey(2)+MOD((ind-1)/2,2)+0.5
#endif
#if NDIM>2
     xcen(3)=2*m%grid(igrid)%ckey(3)+MOD((ind-1)/4,2)+0.5
#endif
     ! Collect all neighboring cell from hash table
     do j=1,nSnei

        ! Compute neighboring cell coordinates
        xnei(1:ndim)=xcen(1:ndim)+xSnei(1:ndim,j)
        ! Periodic boundary conditions
        do idim=1,ndim
           if(r%periodic(idim))then
              if(xnei(idim)< m%box_ckey_min(idim,ilevel+1))xnei(idim)=xnei(idim)-m%box_ckey_min(idim,ilevel+1)+m%box_ckey_max(idim,ilevel+1)
              if(xnei(idim)>=m%box_ckey_max(idim,ilevel+1))xnei(idim)=xnei(idim)+m%box_ckey_min(idim,ilevel+1)-m%box_ckey_max(idim,ilevel+1)
           endif
        end do

        ! Get neighboring cell at ilevel
        ckey_nbor(1:ndim)=int(xnei(1:ndim))
        hash_nbor(0)=ilevel+1
        hash_nbor(1:ndim)=ckey_nbor(1:ndim)
        call get_parent_cell(s,hash_nbor,igridn,icelln,flush_cache=.false.,fetch_cache=.true.)

        ! If missing, get neighboring cell at ilevel-1
        if(igridn==0)then
           ckey_nbor(1:ndim)=int(xnei(1:ndim)/2.0)
           hash_nbor(0)=ilevel
           hash_nbor(1:ndim)=ckey_nbor(1:ndim)
           call get_parent_cell(s,hash_nbor,igridn,icelln,flush_cache=.false.,fetch_cache=.true.)

           ! If refined, get neighboring cell at ilevel+1
        else if (m%grid(igridn)%refined(icelln))then
           ckey_nbor(1:ndim)=int(xnei(1:ndim)*2.0)
           hash_nbor(0)=ilevel+2
           hash_nbor(1:ndim)=ckey_nbor(1:ndim)
           call get_parent_cell(s,hash_nbor,igridn,icelln,flush_cache=.false.,fetch_cache=.true.)
        endif

        ! Lock grid in cache
        call lock_cache(m,igridn)

        igrid_nbor(j) = igridn
        icell_nbor(j) = icelln

     end do

     do j=1,nSnei
        igridn = igrid_nbor(j) ! Gather neighboring grid
        icelln = icell_nbor(j)

        if(igridn>0)then
           peak_nbor = m%flag1(icelln,igridn)
#ifdef GRAV
           dens_nbor = m%rho(icelln,igridn)
#endif
           ok = peak_cen/=0
           ok = ok .and. peak_nbor/=0
           ok = ok .and. peak_cen/=peak_nbor

           ! If saddle density is larger, set new densest saddle point
           if(ok)then
              dens_ave = 0.5*(dens_cen+dens_nbor)
              global_peak_id = peak_cen
              call get_peak(s,global_peak_id,ipeak,flush_cache=.true.,fetch_cache=.false.)
              if(dens_ave>c%saddle_dens(ipeak))then
                 c%saddle_dens(ipeak)=dens_ave
                 c%saddle_nbor(ipeak)=peak_nbor
              end if
           endif
        endif
     end do

     ! Unlock neighboring grids
     do j=1,nSnei
        igridn = igrid_nbor(j)
        call unlock_cache(m,igridn)
     end do

  end do

  call close_cache(mdl)

  end associate

end subroutine collect_saddle
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_fetch_saddle(mesh,igrid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_int4_small_realdp
  type(mesh_t)::mesh
  integer::igrid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind
  type(msg_int4_small_realdp)::msg

  do ind=1,twotondim
     msg%flg(ind)=mesh%flag1(ind,igrid)
  end do

  do ind=1,twotondim
     if(mesh%grid(igrid)%refined(ind))then
        msg%ref(ind)=1
     else
        msg%ref(ind)=0
     endif
  end do

#ifdef GRAV
  do ind=1,twotondim
     msg%realdp(ind)=mesh%rho(ind,igrid)
  end do
#endif

  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_saddle
!################################################################
!################################################################
!################################################################
!################################################################
subroutine unpack_fetch_saddle(mesh,igrid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim, twotondim
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_int4_small_realdp
  type(mesh_t)::mesh
  integer::igrid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind
  type(msg_int4_small_realdp)::msg

  mesh%grid(igrid)%lev=hash_key(0)
  mesh%grid(igrid)%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)

  do ind=1,twotondim
     mesh%flag1(ind,igrid)=msg%flg(ind)
  end do

  do ind=1,twotondim
     if(msg%ref(ind)==1)then
        mesh%grid(igrid)%refined(ind)=.true.
     else
        mesh%grid(igrid)%refined(ind)=.false.
     endif
  end do

#ifdef GRAV
  do ind=1,twotondim
     mesh%rho(ind,igrid)=msg%realdp(ind)
  end do
#endif

end subroutine unpack_fetch_saddle
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_flush_saddle(c,local_peak_id)
  use clfind_commons, only: clump_t
  type(clump_t)::c
  integer::local_peak_id

  c%saddle_dens(local_peak_id)=c%density_threshold
  c%saddle_nbor(local_peak_id)=0

end subroutine init_flush_saddle
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_flush_saddle(c,local_peak_id,msg_size,msg_array)
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_saddle_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_saddle_clump)::msg

  msg%nbor=c%saddle_nbor(local_peak_id)
  msg%dens=c%saddle_dens(local_peak_id)

  msg_array=transfer(msg,msg_array)

end subroutine pack_flush_saddle
!################################################################
!################################################################
!################################################################
!################################################################
subroutine unpack_flush_saddle(c,local_peak_id,msg_size,msg_array)
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_saddle_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_saddle_clump)::msg

  msg=transfer(msg_array,msg)

  if(msg%dens>c%saddle_dens(local_peak_id))then
     c%saddle_dens(local_peak_id)=msg%dens
     c%saddle_nbor(local_peak_id)=msg%nbor
  endif

end subroutine unpack_flush_saddle
!################################################################
!################################################################
!################################################################
!################################################################
subroutine merge_clumps(s,action)
  use amr_commons, only: ndim
  use ramses_commons, only: ramses_t
  use cache_commons, only: msg_merge_clump, msg_prop_clump, msg_halo_clump
  use cache
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  type(ramses_t)::s
  character(len=9)::action
  !---------------------------------------------------------
  ! This routine merges the irrelevant clumps
  ! - clumps are sorted by ascending max density
  ! - irrelevent clumps are merged to most relevant neighbor
  !---------------------------------------------------------
#ifndef WITHOUTMPI
  integer::info
#endif
  integer::j,i,ipart,igrid,ind,itest
  integer::current,nmove,ipeak,jpeak,iter
  integer::nsurvive,nzero,idepth
  integer::ilev,ilevel
  integer(kind=8)::global_peak_id,merge_to
  real(kind=8)::value_iij,zero=0,relevance_peak
  real(kind=8)::d,dx_loc,vol
  integer,dimension(1:s%c%npeak_max)::alive,ind_sort
  real(kind=8),dimension(1:s%c%npeak_max)::peakd
  logical::do_merge=.false.
  type(msg_merge_clump)::dummy_merge_clump
  type(msg_prop_clump)::dummy_prop_clump
  type(msg_halo_clump)::dummy_halo_clump
#ifndef WITHOUTMPI
  integer::nmove_all,nsurvive_all,nzero_all
#endif

  associate(g=>s%g,r=>s%r,m=>s%m,c=>s%c,mdl=>s%mdl)

  if (r%verbose.and.g%myid==1)then
     if(action.EQ.'relevance')then
        write(*,*)'Now merging irrelevant clumps'
     endif
     if(action.EQ.'saddleden')then
        write(*,*)'Now merging clumps into halos'
     endif
  endif

  ! Initialize new_peak array to global peak id
  ! All peaks are alive at the start
  do i=1,c%npeak
     c%new_peak(i)=i+c%npeak_cum(g%myid-1)
     c%tidal_dens(i)=c%saddle_threshold
     if(action.EQ.'relevance')then
        alive(i)=1
        c%lev_peak(i)=-1
     endif
     if(action.EQ.'saddleden')then
        if(c%relevance(i)>c%relevance_threshold)then
           alive(i)=1
        else
           alive(i)=0
        endif
        c%lev_peak(i)=-1
     endif
  end do

  ! Sort peaks by maximum peak density in ascending order
  do i=1,c%npeak
     peakd(i)=c%max_dens(i)
     ind_sort(i)=i
  end do
  call quick_sort_dp(peakd,ind_sort,c%npeak)

  ! Loop over peak levels
  nzero=c%npeak_tot
  idepth=0
  do while(nzero>0)

     ! Merge peaks
     nmove=c%npeak_tot
     iter=0
     do while(nmove>0)

        call open_cache_clump(mdl, c, pack_size=storage_size(dummy_merge_clump)/32, &
             pack=pack_fetch_merge, unpack=unpack_fetch_merge)

        nmove=0
        do i=c%npeak,1,-1
           ipeak=ind_sort(i)
           merge_to=c%new_peak(ipeak)
           if(alive(ipeak)>0)then
              if(action.EQ.'relevance')then
                 relevance_peak=c%max_dens(ipeak)/c%saddle_dens(ipeak)
                 do_merge=relevance_peak<c%relevance_threshold
              endif
              if(action.EQ.'saddleden')then
                 do_merge=(c%saddle_dens(ipeak)>c%saddle_threshold)
              endif
              if(do_merge)then
                 if(c%saddle_nbor(ipeak)>0)then
                    global_peak_id=c%saddle_nbor(ipeak)
                    call get_peak(s,global_peak_id,jpeak,flush_cache=.false.,fetch_cache=.true.)
                    if(c%max_dens(jpeak)>c%max_dens(ipeak))then
                       merge_to=c%new_peak(jpeak)
                    else if(c%max_dens(jpeak)==c%max_dens(ipeak))then
                       merge_to=MIN(c%new_peak(ipeak),c%new_peak(jpeak))
                    endif
                 endif
              endif
           endif
           if(c%new_peak(ipeak).NE.merge_to)then
              nmove=nmove+1
              c%new_peak(ipeak)=merge_to
              c%tidal_dens(ipeak)=c%saddle_dens(ipeak)
           endif
        end do

        call close_cache(mdl)

        iter=iter+1
#ifndef WITHOUTMPI
        call MPI_ALLREDUCE(nmove,nmove_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
        nmove=nmove_all
#endif
        if(g%myid==1.and.r%verbose)write(*,*)'niter=',iter,'nmove=',nmove
     end do

     ! Update flag1 field
     call open_cache_clump(mdl, c, pack_size=storage_size(dummy_merge_clump)/32, &
          pack=pack_fetch_merge, unpack=unpack_fetch_merge)
     do itest=1,c%ntest
        igrid=c%grid(itest)
        ind=c%cell(itest)
        global_peak_id=m%flag1(ind,igrid)
        if (global_peak_id>0)then
           call get_peak(s,global_peak_id,ipeak,flush_cache=.false.,fetch_cache=.true.)
           m%flag1(ind,igrid)=c%new_peak(ipeak)
        end if
     end do
     call close_cache(mdl)

     ! Compute new saddle points and corresponding nboring peaks.
     call collect_saddle(s)

     ! Set alive to zero for newly merged peaks
     nzero=0
     nsurvive=0
     do ipeak=1,c%npeak
        if(alive(ipeak)>0)then
           merge_to=c%new_peak(ipeak)
           if(merge_to.NE.ipeak+c%npeak_cum(g%myid-1))then
              alive(ipeak)=0
              c%lev_peak(ipeak)=idepth
              nzero=nzero+1
           else
              nsurvive=nsurvive+1
           end if
        endif
     end do

#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(nzero,nzero_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
     nzero=nzero_all
     call MPI_ALLREDUCE(nsurvive,nsurvive_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
     nsurvive=nsurvive_all
#endif
     if(r%verbose.and.g%myid==1)write(*,*)'level=',idepth,'nmove=',nzero,'survived=',nsurvive
     idepth=idepth+1

  end do
  ! End loop over peak levels

  ! Last level has no more clumps, also idepth=idepth+1 still happens on last level.
  c%merge_levelmax=idepth-2

  ! Count surviving peaks
  nsurvive=0
  do ipeak=1,c%npeak
     if(alive(ipeak)>0)then
        nsurvive=nsurvive+1
     endif
  end do
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(nsurvive,nsurvive_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  nsurvive=nsurvive_all
#endif
  if(g%myid==1)then
     if(action.EQ.'relevance')then
        write(*,*)'Found',nsurvive,' relevant peaks'
     endif
     if(action.EQ.'saddleden')then
        write(*,*)'Found',nsurvive,' halos'
     endif
  endif

  if(action.EQ.'relevance')then

     ! Compute relevance
     do ipeak=1,c%npeak
        if(alive(ipeak)>0)then
           c%relevance(ipeak)=c%max_dens(ipeak)/c%saddle_dens(ipeak)
        else
           c%relevance(ipeak)=0
        endif
     end do

     ! Merge all peaks to deepest level
     do ilev=idepth-2,0,-1
        call open_cache_clump(mdl, c, pack_size=storage_size(dummy_merge_clump)/32, &
             pack=pack_fetch_merge, unpack=unpack_fetch_merge)
        do ipeak=1,c%npeak
           if(c%lev_peak(ipeak)==ilev)then
              merge_to=c%new_peak(ipeak)
              call get_peak(s,merge_to,jpeak,flush_cache=.false.,fetch_cache=.true.)
              c%new_peak(ipeak)=c%new_peak(jpeak)
           endif
        end do
        call close_cache(mdl)
     end do

  endif

  if(action.EQ.'saddleden')then

     ! Compute halo index
     do ipeak=1,c%npeak
        c%ind_halo(ipeak)=c%new_peak(ipeak)
     end do
     do ilev=idepth-2,0,-1
        call open_cache_clump(mdl, c, pack_size=storage_size(dummy_halo_clump)/32, &
             pack=pack_fetch_halo, unpack=unpack_fetch_halo)
        do ipeak=1,c%npeak
           if(c%lev_peak(ipeak)==ilev)then
              merge_to=c%ind_halo(ipeak)
              call get_peak(s,merge_to,jpeak,flush_cache=.false.,fetch_cache=.true.)
              c%ind_halo(ipeak)=c%ind_halo(jpeak)
           endif
        end do
        call close_cache(mdl)
     end do

     ! Set merging level of halos to idepth-1
     do ipeak=1,c%npeak
        if(c%ind_halo(ipeak)==ipeak+c%npeak_cum(g%myid-1))then
           c%lev_peak(ipeak)=idepth-1
        endif
     end do

     ! Compute halo mass
     c%halo_mass=0
     call open_cache_clump(mdl, c, pack_size=storage_size(dummy_halo_clump)/32, &
          init=init_flush_halo, flush=pack_flush_halo, combine=unpack_flush_halo)
     do ipeak=1,c%npeak
        merge_to=c%ind_halo(ipeak)
        call get_peak(s,merge_to,jpeak,flush_cache=.true.,fetch_cache=.false.)
        c%halo_mass(jpeak)=c%halo_mass(jpeak)+c%clump_mass(ipeak)
     end do
     call close_cache(mdl)

     ! Assign back halo mass to peak
     call open_cache_clump(mdl, c, pack_size=storage_size(dummy_halo_clump)/32, &
          pack=pack_fetch_halo, unpack=unpack_fetch_halo)
     do ipeak=1,c%npeak
        merge_to=c%ind_halo(ipeak)
        call get_peak(s,merge_to,jpeak,flush_cache=.false.,fetch_cache=.true.)
        c%halo_mass(ipeak)=c%halo_mass(jpeak)
     end do
     call close_cache(mdl)

     ! Compute hierarchical clump properties
     do ilev=0,idepth-2
        call open_cache_clump(mdl, c, pack_size=storage_size(dummy_prop_clump)/32, &
             init=init_flush_prop, flush=pack_flush_prop, combine=unpack_flush_prop)
        do ipeak=1,c%npeak
           if(c%lev_peak(ipeak)==ilev)then
              ! If clump is too massive, then it does not merge
              if(c%clump_mass(ipeak).GT.c%fraction_threshold*c%halo_mass(ipeak))then
                 c%new_peak(ipeak)=ipeak+c%npeak_cum(g%myid-1)
                 c%lev_peak(ipeak)=idepth-1
              endif
              merge_to=c%new_peak(ipeak)
              if(merge_to.NE.ipeak+c%npeak_cum(g%myid-1))then
                 call get_peak(s,merge_to,jpeak,flush_cache=.true.,fetch_cache=.false.)
                 c%clump_mass(jpeak)=c%clump_mass(jpeak)+c%clump_mass(ipeak)
                 c%clump_vol(jpeak)=c%clump_vol(jpeak)+c%clump_vol(ipeak)
                 c%n_cells(jpeak)=c%n_cells(jpeak)+c%n_cells(ipeak)
                 c%min_dens(jpeak)=min(c%min_dens(jpeak),c%min_dens(ipeak))
              endif
           endif
        end do
        call close_cache(mdl)
     end do

  endif

  end associate

end subroutine merge_clumps
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_fetch_merge(c,local_peak_id,msg_size,msg_array)
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_merge_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_merge_clump)::msg

  msg%npeak=c%new_peak(local_peak_id)
  msg%mdens=c%max_dens(local_peak_id)

  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_merge
!################################################################
!################################################################
!################################################################
!################################################################
subroutine unpack_fetch_merge(c,local_peak_id,msg_size,msg_array)
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_merge_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_merge_clump)::msg

  msg=transfer(msg_array,msg)

  c%max_dens(local_peak_id)=msg%mdens
  c%new_peak(local_peak_id)=msg%npeak

end subroutine unpack_fetch_merge
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_fetch_halo(c,local_peak_id,msg_size,msg_array)
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_halo_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_halo_clump)::msg

  msg%ihalo=c%ind_halo(local_peak_id)
  msg%mhalo=c%halo_mass(local_peak_id)

  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_halo
!################################################################
!################################################################
!################################################################
!################################################################
subroutine unpack_fetch_halo(c,local_peak_id,msg_size,msg_array)
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_halo_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_halo_clump)::msg

  msg=transfer(msg_array,msg)

  c%ind_halo(local_peak_id)=msg%ihalo
  c%halo_mass(local_peak_id)=msg%mhalo

end subroutine unpack_fetch_halo
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_flush_halo(c,local_peak_id)
  use clfind_commons, only: clump_t
  type(clump_t)::c
  integer::local_peak_id

  c%halo_mass(local_peak_id)=0d0

end subroutine init_flush_halo
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_flush_halo(c,local_peak_id,msg_size,msg_array)
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_halo_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_halo_clump)::msg

  msg%mhalo=c%halo_mass(local_peak_id)

  msg_array=transfer(msg,msg_array)

end subroutine pack_flush_halo
!################################################################
!################################################################
!################################################################
!################################################################
subroutine unpack_flush_halo(c,local_peak_id,msg_size,msg_array)
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_halo_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_halo_clump)::msg

  msg=transfer(msg_array,msg)

  c%halo_mass(local_peak_id)=c%halo_mass(local_peak_id)+msg%mhalo

end subroutine unpack_flush_halo
!################################################################
!################################################################
!################################################################
!################################################################
subroutine compute_clump_properties(s,rtype)
  use amr_commons, only: ndim
  use clfind_commons
  use ramses_commons, only: ramses_t
  use cache_commons, only: msg_prop_clump
  use cache
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
  type(ramses_t)::s
  integer::rtype
  !----------------------------------------------------------------------------
  ! This subroutine performs a loop over all cells above the threshold and
  ! collects the  relevant information. After some MPI communications,
  ! all necessary peak-patch properties are computed
  !----------------------------------------------------------------------------
  type(msg_prop_clump)::dummy_prop_clump
  integer(kind=8)::global_peak_id
  integer::ipart,grid,peak_nr,ilevel,ipeak,plevel,igrid,itest,idim,ind
  real(kind=8),dimension(1:ndim)::xcell,accel
  real(kind=8)::dx_loc,tot_mass
  real(kind=8)::zero=0
  real(kind=8)::d=0, vol=0, nref=0, pi
#ifndef WITHOUTMPI
  integer::i
  real(kind=8)::tot_mass_tot
#endif

  associate(g=>s%g,r=>s%r,m=>s%m,c=>s%c,mdl=>s%mdl)

  if(g%myid==1.and.r%verbose)write(*,*)'Entering compute clump properties'

  ! Constants
  pi=ACOS(-1.0D0)

  !-----------------------------------------------------------------------
  ! Loop over local peaks and compute peak cell coordinates, velocities...
  !-----------------------------------------------------------------------
  c%peak_pos=0d0; c%peak_com=0d0; c%peak_vel=0d0; c%peak_acc=0d0

  do ipeak=1,c%npeak
     ilevel=c%peak_level(ipeak)
     igrid=c%peak_grid(ipeak)
     ind=c%peak_cell(ipeak)
     dx_loc=r%boxlen/2**ilevel
     ! Peak cell coordinates and acceleration
     xcell(1)=(2*m%grid(igrid)%ckey(1)+MOD((ind-1)  ,2)+0.5)*dx_loc-m%skip(1)
#if NDIM>1
     xcell(2)=(2*m%grid(igrid)%ckey(2)+MOD((ind-1)/2,2)+0.5)*dx_loc-m%skip(2)
#endif
#if NDIM>2
     xcell(3)=(2*m%grid(igrid)%ckey(3)+MOD((ind-1)/4,2)+0.5)*dx_loc-m%skip(3)
#endif
     c%peak_pos(ipeak,1:ndim)=xcell(1:ndim)
#ifdef HYDRO
     if (r%hydro.AND.rtype.eq.4)then
        c%peak_vel(ipeak,1:ndim)=m%uold(ind,2:ndim+1,igrid)/m%uold(ind,1,igrid)
        c%peak_com(ipeak,1:ndim)=xcell(1:ndim)
     endif
#endif
#ifdef GRAV
     c%peak_acc(ipeak,1:ndim)=m%f(ind,1:ndim,igrid)
#endif
  end do

  !--------------------------------------------------------------------------
  ! Loop over all cells above the threshold and compute mass, volume...
  !--------------------------------------------------------------------------
  c%min_dens=huge(zero); c%n_cells=0
  c%clump_mass=0d0; c%clump_vol=0d0

  call open_cache_clump(mdl, c, storage_size(dummy_prop_clump)/32, &
       init=init_flush_prop, flush=pack_flush_prop, combine=unpack_flush_prop)
  do itest=1,c%ntest
     ilevel=c%level(itest)
     igrid=c%grid(itest)
     ind=c%cell(itest)
     global_peak_id=m%flag1(ind,igrid)

     ! Save peak patch id into flag2 because flag1 will become halo patch id
     m%flag2(ind,igrid)=m%flag1(ind,igrid)

     if (global_peak_id /=0 ) then
        call get_peak(s,global_peak_id,peak_nr,flush_cache=.true.,fetch_cache=.false.)

        ! Cell density
#ifdef GRAV
        d=m%rho(ind,igrid)
        nref=m%nref(ind,igrid)
#endif
        ! Cell volume
        dx_loc=r%boxlen/2**ilevel
        vol=dx_loc**ndim

        ! Number of leaf cells per clump
        c%n_cells(peak_nr)=c%n_cells(peak_nr)+1

        ! Clump min density
        c%min_dens(peak_nr)=min(c%min_dens(peak_nr),d)

        ! Clump mass
        c%clump_mass(peak_nr)=c%clump_mass(peak_nr)+vol*nref

        ! Clump volume
        c%clump_vol(peak_nr)=c%clump_vol(peak_nr)+vol

     end if
  end do
  call close_cache(mdl)

  ! Calculate clump tidal radius
  do ipeak=1,c%npeak
     c%clump_rad(ipeak)=(c%clump_mass(ipeak)*3d0/4d0/pi/c%tidal_dens(ipeak))**(1d0/3d0)
  end do

  ! Calculate total mass above threshold
  tot_mass=sum(c%clump_mass(1:c%npeak))

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(tot_mass,tot_mass_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  tot_mass=tot_mass_tot
#endif

  end associate

end subroutine compute_clump_properties
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_flush_prop(c,local_peak_id)
  use amr_commons, only: ndim
  use clfind_commons, only: clump_t
  type(clump_t)::c
  integer::local_peak_id
  real(kind=8)::zero

  c%n_cells(local_peak_id)=0
  c%min_dens(local_peak_id)=huge(zero)
  c%clump_mass(local_peak_id)=0d0
  c%clump_vol(local_peak_id)=0d0

end subroutine init_flush_prop
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_flush_prop(c,local_peak_id,msg_size,msg_array)
  use amr_commons, only: ndim
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_prop_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_prop_clump)::msg

  msg%ncell=c%n_cells(local_peak_id)
  msg%dens=c%min_dens(local_peak_id)
  msg%mass=c%clump_mass(local_peak_id)
  msg%vol=c%clump_vol(local_peak_id)

  msg_array=transfer(msg,msg_array)

end subroutine pack_flush_prop
!################################################################
!################################################################
!################################################################
!################################################################
subroutine unpack_flush_prop(c,local_peak_id,msg_size,msg_array)
  use amr_commons, only: ndim
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_prop_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_prop_clump)::msg

  msg=transfer(msg_array,msg)

  c%n_cells(local_peak_id)=c%n_cells(local_peak_id)+msg%ncell
  c%min_dens(local_peak_id)=min(c%min_dens(local_peak_id),msg%dens)
  c%clump_mass(local_peak_id)=c%clump_mass(local_peak_id)+msg%mass
  c%clump_vol(local_peak_id)=c%clump_vol(local_peak_id)+msg%vol

end subroutine unpack_flush_prop
!################################################################
!################################################################
!################################################################
!################################################################
subroutine trim_clumps(s)
#ifndef WITHOUTMPI
  use mpi
#endif
  use amr_commons, only: ndim
  use clfind_commons
  use ramses_commons, only: ramses_t
  use cache_commons, only: msg_prop_clump
  use cache
  implicit none
  type(ramses_t)::s
  !----------------------------------------------------------------------------
  ! This subroutine remove all clumps and halos that are considered irrelevant.
  ! They are removed because their relevance (or peakiness) is below the
  ! relevance threshold or because their mass is too small.
  ! The flag1 and flag2 arrays are modified accordingly.
  !----------------------------------------------------------------------------
#ifndef WITHOUTMPI
  integer::info
  integer,dimension(1:s%g%ncpu)::nfinal_cpu_all
#endif
  type(msg_prop_clump)::dummy_prop_clump
  integer(kind=8)::global_peak_id,global_halo_id
  integer::nfinal
  integer(kind=8)::nfinal_tot
  integer,dimension(1:s%g%ncpu)::nfinal_cpu
  integer(kind=8),dimension(0:s%g%ncpu)::nfinal_cum
  integer::icpu,ipeak,jpeak,igrid,ilevel,ind,itest

  associate(g=>s%g,r=>s%r,m=>s%m,c=>s%c,mdl=>s%mdl)

  if(g%myid==1.and.r%verbose)write(*,*)'Entering trim_clumps'

  ! Change peak indexing to keep only relevant ones
  nfinal=0
  c%ind_final=0
  do ipeak=1,c%npeak
     if( c%clump_mass(ipeak) > c%mass_threshold.AND. &
          & c%relevance(ipeak) > c%relevance_threshold.AND. &
          & c%npart(ipeak) > 0 ) then
        nfinal=nfinal+1
        c%ind_final(ipeak)=nfinal
     endif
  end do
  nfinal_cpu=0
  nfinal_cpu(g%myid)=nfinal
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(nfinal_cpu,nfinal_cpu_all,g%ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  nfinal_cpu=nfinal_cpu_all
#endif
  nfinal_cum=0
  do icpu=1,g%ncpu
     nfinal_cum(icpu)=nfinal_cum(icpu-1)+int(nfinal_cpu(icpu),kind=8)
  end do
#ifdef WITHOUTMPI
  nfinal_tot=int(nfinal,kind=8)
#else
  nfinal_tot=nfinal_cum(g%ncpu)
#endif
!!$  if (g%myid==1)then
!!$    write(*,*)'Found',int(nfinal_tot,kind=4),' clumps'
!!$  endif
  do ipeak=1,c%npeak
     if( c%clump_mass(ipeak) > c%mass_threshold.AND. &
          & c%relevance(ipeak) > c%relevance_threshold.AND. &
          & c%npart(ipeak) > 0 ) then
        c%ind_final(ipeak)=c%ind_final(ipeak)+nfinal_cum(g%myid-1)
     endif
  end do

  ! Change other important peak indices
  call open_cache_clump(mdl, c, storage_size(dummy_prop_clump)/32, &
       pack=pack_fetch_prop, unpack=unpack_fetch_prop)
  do ipeak=1,c%npeak
     if( c%clump_mass(ipeak) > c%mass_threshold.AND. &
          & c%relevance(ipeak) > c%relevance_threshold.AND. &
          & c%npart(ipeak) > 0 ) then
        global_peak_id=c%new_peak(ipeak)
        call get_peak(s,global_peak_id,jpeak,flush_cache=.false.,fetch_cache=.true.)
        c%new_peak(ipeak)=c%ind_final(jpeak)
        global_peak_id=c%ind_halo(ipeak)
        call get_peak(s,global_peak_id,jpeak,flush_cache=.false.,fetch_cache=.true.)
        c%ind_halo(ipeak)=c%ind_final(jpeak)
     endif
  end do
  call close_cache(mdl)

  ! Modify particle peak ids accordingly
  if(s%r%part)then
     call trim_particles(s,s%p)
  endif
  if(s%r%star)then
     call trim_particles(s,s%star)
  endif
  if(s%r%sink)then
     call trim_particles(s,s%sink)
  endif
  if(s%r%tree)then
     call trim_particles(s,s%tree)
  endif

  ! Modify flag1 and flag2 accordingly
  call open_cache_clump(mdl, c, storage_size(dummy_prop_clump)/32, &
       pack=pack_fetch_prop, unpack=unpack_fetch_prop)
  do itest=1,c%ntest
     ilevel=c%level(itest)
     igrid=c%grid(itest)
     ind=c%cell(itest)
     global_peak_id=m%flag2(ind,igrid)
     if (global_peak_id /=0 ) then
        call get_peak(s,global_peak_id,ipeak,flush_cache=.false.,fetch_cache=.true.)
        if( c%clump_mass(ipeak) > c%mass_threshold.AND. &
             & c%relevance(ipeak) > c%relevance_threshold.AND. &
             & c%npart(ipeak) > 0 ) then
           m%flag2(ind,igrid)=c%ind_final(ipeak)
        else
           m%flag2(ind,igrid)=0
        endif
     endif
     if(c%saddle_threshold>0)then
        global_halo_id=m%flag1(ind,igrid)
        if (global_halo_id /=0 ) then
           call get_peak(s,global_halo_id,ipeak,flush_cache=.false.,fetch_cache=.true.)
           if( c%clump_mass(ipeak) > c%mass_threshold.AND. &
                & c%relevance(ipeak) > c%relevance_threshold.AND. &
                & c%npart(ipeak) > 0 ) then
              m%flag1(ind,igrid)=c%ind_final(ipeak)
           else
              m%flag1(ind,igrid)=0
           endif
        endif
     endif
  end do
  call close_cache(mdl)

  end associate

end subroutine trim_clumps
!################################################################
!################################################################
!################################################################
!################################################################
subroutine trim_particles(s,p)
  use amr_parameters, only: ndim, nbin, twotondim
  use ramses_commons, only: ramses_t
  use pm_commons, only: part_t
  use cache_commons
  use cache
  implicit none
  type(ramses_t)::s
  type(part_t)::p
  !------------------------------------------------------------------
  ! This routine modifies the peak id of the input particles
  ! to accound for the final clump index.
  ! Written by Romain Teyssier (mini-ramses version in June 2024).
  !------------------------------------------------------------------
  type(msg_prop_clump)::dummy_prop_clump
  integer::ipeak,ipart
  integer(kind=8)::global_peak_id

  associate(r=>s%r,g=>s%g,m=>s%m,c=>s%c,mdl=>s%mdl)

  call open_cache_clump(mdl, c, storage_size(dummy_prop_clump)/32, &
       pack=pack_fetch_prop, unpack=unpack_fetch_prop)
  do ipart=1,p%npart
     global_peak_id=p%pid(ipart)
     if (global_peak_id /=0 ) then
        call get_peak(s,global_peak_id,ipeak,flush_cache=.false.,fetch_cache=.true.)
        if( c%clump_mass(ipeak) > c%mass_threshold.AND. &
             & c%relevance(ipeak) > c%relevance_threshold.AND. &
             & c%npart(ipeak) > 0 ) then
           p%pid(ipart)=c%ind_final(ipeak)
        else
           p%pid(ipart)=0
        endif
     endif
  end do
  call close_cache(mdl)

  end associate

end subroutine trim_particles
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_fetch_prop(c,local_peak_id,msg_size,msg_array)
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_prop_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_prop_clump)::msg

  msg%ind=c%ind_final(local_peak_id)
  msg%ncell=c%npart(local_peak_id)
  msg%dens=c%relevance(local_peak_id)
  msg%mass=c%clump_mass(local_peak_id)

  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_prop
!################################################################
!################################################################
!################################################################
!################################################################
subroutine unpack_fetch_prop(c,local_peak_id,msg_size,msg_array)
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_prop_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_prop_clump)::msg

  msg=transfer(msg_array,msg)

  c%ind_final(local_peak_id)=msg%ind
  c%npart(local_peak_id)=msg%ncell
  c%relevance(local_peak_id)=msg%dens
  c%clump_mass(local_peak_id)=msg%mass

end subroutine unpack_fetch_prop
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine particle_clump_properties(s,p)
  use amr_parameters, only: ndim, nbin, twotondim
  use ramses_commons, only: ramses_t
  use pm_commons, only: part_t
  use cache_commons
  use cache
  implicit none
  type(ramses_t)::s
  type(part_t)::p
  !------------------------------------------------------------------
  ! This routine computes various particle-based clump properties.
  ! The particle peak id must be available in p%pid.
  ! Note that this routine must be called prior to unbinding.
  ! Written by Romain Teyssier (mini-ramses version in June 2024).
  !------------------------------------------------------------------
  type(msg_prop_clump)::dummy_prop_clump
  integer::i,ipeak,ipart,ind,idim,ibin,ilevel
  integer(kind=8)::global_peak_id
  integer::halo_nr,peak_nr
  real(kind=8)::xx,rad,dx_loc,r2,core_radius
  real(kind=8),dimension(1:ndim)::xpart

  associate(r=>s%r,g=>s%g,m=>s%m,c=>s%c,mdl=>s%mdl)

  !--------------------------------------------
  ! Sort particles according to global clump id
  !--------------------------------------------
  do ipart=1,p%npart
     p%sortp(ipart)=ipart
     p%workp(ipart)=p%pid(ipart)
  end do
  call quick_sort_int_int(p%workp(1),p%sortp(1),p%npart)

  !-----------------------------------------
  ! Compute peak velocity based on particles
  !-----------------------------------------
  c%particle_mass=0
  c%peak_vel=0
  c%peak_com=0
  c%npart=0
  call open_cache_clump(mdl, c, storage_size(dummy_prop_clump)/32, &
       pack=pack_fetch_part, unpack=unpack_fetch_part, &
       init=init_flush_part, flush=pack_flush_part, combine=unpack_flush_part)
  do i=1+p%norphan_peak,p%npart
     ! Get peak id
     ipart=p%sortp(i)
     global_peak_id=p%workp(i)
     if (global_peak_id /=0 ) then
        call get_peak(s,global_peak_id,peak_nr,flush_cache=.true.,fetch_cache=.true.)
        ! Compute distance to cell center
        xpart(1:ndim)=p%xp(ipart,1:ndim)-c%peak_pos(peak_nr,1:ndim)
        r2=0d0
        do idim=1,ndim
           ! In case of periodic boundaries
           if(r%periodic(idim))then
              if(xpart(idim)> r%box_size(idim)*0.5)xpart(idim)=xpart(idim)-r%box_size(idim)
              if(xpart(idim)<-r%box_size(idim)*0.5)xpart(idim)=xpart(idim)+r%box_size(idim)
           endif
           r2=r2+xpart(idim)**2
        end do
        ! Keep only particles within 4-cell or 10% of tidal radius
        ilevel=c%peak_level(peak_nr)
        dx_loc=r%boxlen/2**ilevel
        core_radius=2d0*dx_loc
        if(r2<=core_radius**2)then
           c%peak_com(peak_nr,1:ndim)=c%peak_com(peak_nr,1:ndim)+p%mp(ipart)*xpart(1:ndim)
           c%peak_vel(peak_nr,1:ndim)=c%peak_vel(peak_nr,1:ndim)+p%mp(ipart)*p%vp(ipart,1:ndim)
           c%particle_mass(peak_nr)=c%particle_mass(peak_nr)+p%mp(ipart)
           c%npart(peak_nr)=c%npart(peak_nr)+1
        endif
     endif
  end do
  call close_cache(mdl)

  ! Compute specific quantities
  do ipeak=1,c%npeak
     if (c%particle_mass(ipeak)>0)then
        c%peak_vel(ipeak,1:ndim)=c%peak_vel(ipeak,1:ndim)/c%particle_mass(ipeak)
        c%peak_com(ipeak,1:ndim)=c%peak_com(ipeak,1:ndim)/c%particle_mass(ipeak)
     end if
     c%peak_com(ipeak,1:ndim)=c%peak_com(ipeak,1:ndim)+c%peak_pos(ipeak,1:ndim)
     ! Periodic boundary conditions
     do idim=1,ndim
        if(r%periodic(idim))then
           if(c%peak_com(ipeak,idim)< 0.0d0           )c%peak_com(ipeak,idim)=c%peak_com(ipeak,idim)+r%box_size(idim)
           if(c%peak_com(ipeak,idim)>=r%box_size(idim))c%peak_com(ipeak,idim)=c%peak_com(ipeak,idim)-r%box_size(idim)
        end if
     end do
  end do

  end associate

end subroutine particle_clump_properties
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_fetch_part(c,local_peak_id,msg_size,msg_array)
  use amr_commons, only: ndim
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_prop_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_prop_clump)::msg

  msg%pos(1:ndim)=c%peak_pos(local_peak_id,1:ndim)
  msg%ind=c%peak_level(local_peak_id)
  msg%rad=c%clump_rad(local_peak_id)

  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_part
!################################################################
!################################################################
!################################################################
!################################################################
subroutine unpack_fetch_part(c,local_peak_id,msg_size,msg_array)
  use amr_commons, only: ndim
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_prop_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_prop_clump)::msg

  msg=transfer(msg_array,msg)

  c%peak_pos(local_peak_id,1:ndim)=msg%pos(1:ndim)
  c%peak_level(local_peak_id)=msg%ind
  c%clump_rad(local_peak_id)=msg%rad

end subroutine unpack_fetch_part
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_flush_part(c,local_peak_id)
  use amr_commons, only: ndim
  use clfind_commons, only: clump_t
  type(clump_t)::c
  integer::local_peak_id

  c%npart(local_peak_id)=0
  c%particle_mass(local_peak_id)=0
  c%peak_vel(local_peak_id,1:ndim)=0d0
  c%peak_com(local_peak_id,1:ndim)=0d0

end subroutine init_flush_part
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_flush_part(c,local_peak_id,msg_size,msg_array)
  use amr_commons, only: ndim
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_prop_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_prop_clump)::msg

  msg%ncell=c%npart(local_peak_id)
  msg%mass=c%particle_mass(local_peak_id)
  msg%vel(1:ndim)=c%peak_vel(local_peak_id,1:ndim)
  msg%pos(1:ndim)=c%peak_com(local_peak_id,1:ndim)

  msg_array=transfer(msg,msg_array)

end subroutine pack_flush_part
!################################################################
!################################################################
!################################################################
!################################################################
subroutine unpack_flush_part(c,local_peak_id,msg_size,msg_array)
  use amr_commons, only: ndim
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_prop_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_prop_clump)::msg

  msg=transfer(msg_array,msg)

  c%npart(local_peak_id)=c%npart(local_peak_id)+msg%ncell
  c%particle_mass(local_peak_id)=c%particle_mass(local_peak_id)+msg%mass
  c%peak_vel(local_peak_id,1:ndim)=c%peak_vel(local_peak_id,1:ndim)+msg%vel(1:ndim)
  c%peak_com(local_peak_id,1:ndim)=c%peak_com(local_peak_id,1:ndim)+msg%pos(1:ndim)

end subroutine unpack_flush_part
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine particle_potential(s,p)
  use amr_parameters, only: ndim, nbin, twotondim
  use ramses_commons, only: ramses_t
  use pm_commons, only: part_t
  use cache_commons
  use cache
  implicit none
  type(ramses_t)::s
  type(part_t)::p
  !------------------------------------------------------------------
  ! This routine computes the particle self-potential of each clump.
  ! Particle are detached hierarchically from children clumps
  ! to their parent clumps in the saddle point merging hierarchy.
  ! The halo-patch is used as a garbage colector.
  ! Written by Romain Teyssier (mini-ramses version in June 2024).
  !------------------------------------------------------------------
  type(msg_mbin_clump)::dummy_mbin_clump
  integer::i,ipart,ind,idim,ibin,ilevel
  integer(kind=8)::global_peak_id
  integer::ipeak
  real(kind=8)::pi,grav,rad,dist,dr,bound
  real(kind=8),dimension(1:ndim)::xpart,vpart

  associate(r=>s%r,g=>s%g,m=>s%m,c=>s%c,mdl=>s%mdl)

  ! Constants
  pi=ACOS(-1.0D0)
  grav=1d0
  if(s%r%cosmo)grav=3d0/8d0/pi*s%g%omega_m*s%g%aexp

  ! Sort particles according to global peak id
  do ipart=1,p%npart
     p%sortp(ipart)=ipart
     p%workp(ipart)=p%pid(ipart)
  end do
  call quick_sort_int_int(p%workp(1),p%sortp(1),p%npart)

  ! Loop over merging hierarchy levels
  c%mass_bin=0d0
  c%npart=0
  do ilevel=0,c%merge_levelmax+1
     ! Open cache
     call open_cache_clump(mdl, c, storage_size(dummy_mbin_clump)/32, &
          pack=pack_fetch_mbin, unpack=unpack_fetch_mbin, &
          init=init_flush_mbin, flush=pack_flush_mbin, combine=unpack_flush_mbin)
     do i=1+p%norphan_peak,p%npart
        ! Get peak id
        ipart=p%sortp(i)
        global_peak_id=p%workp(i)
        call get_peak(s,global_peak_id,ipeak,flush_cache=.true.,fetch_cache=.true.)
        if(c%lev_peak(ipeak)==ilevel)then
           ! Get clump tidal radius
           rad=c%clump_rad(ipeak)
           ! Compute particle radius
           xpart(1:ndim)=p%xp(ipart,1:ndim)-c%peak_pos(ipeak,1:ndim)
           dist=0d0
           do idim=1,ndim
              ! In case of periodic boundaries
              if(r%periodic(idim))then
                 if(xpart(idim)> r%box_size(idim)*0.5)xpart(idim)=xpart(idim)-r%box_size(idim)
                 if(xpart(idim)<-r%box_size(idim)*0.5)xpart(idim)=xpart(idim)+r%box_size(idim)
              endif
              dist=dist+xpart(idim)**2
           end do
           dist=sqrt(dist)
           dr=2d0*rad/dble(nbin)
           do ibin=1,nbin
              ! We use a simple linear binning as the mass is usually propto r
              if(dist<=dble(ibin)*dr)then
                 c%mass_bin(ipeak,ibin)=c%mass_bin(ipeak,ibin)+p%mp(ipart)
                 c%npart(ipeak)=c%npart(ipeak)+1
                 exit
              endif
           end do
           ! Assign particle to next peak in hierarchy
           p%workp(i)=c%new_peak(ipeak)
           p%pid(ipart)=c%new_peak(ipeak)
        endif
     end do
     call close_cache(mdl)
  end do
  ! End loop over levels

  !------------------------
  ! Compute cumulative mass
  !------------------------
  do ipeak=1,c%npeak
     if(   c%clump_mass(ipeak) > c%mass_threshold.AND. &
          & c%relevance(ipeak) > c%relevance_threshold)then
        do ibin=1,nbin-1
           c%mass_bin(ipeak,ibin+1)=c%mass_bin(ipeak,ibin+1)+c%mass_bin(ipeak,ibin)
        end do
     endif
  end do

  !-----------------------
  ! Compute self-potential
  !-----------------------
  c%phi=0d0
  do ipeak=1,c%npeak
     if(   c%clump_mass(ipeak) > c%mass_threshold.AND. &
          & c%relevance(ipeak) > c%relevance_threshold)then
        ! Get clump tidal radius
        rad=c%clump_rad(ipeak)
        ! Convert mass profile into potential energy
        dr=2d0*rad/dble(nbin)
        dist=(dble(nbin-1)+0.5d0)*dr
        c%phi(ipeak,nbin)=grav*c%mass_bin(ipeak,nbin)/dist**2*dr
        do ibin=nbin-1,1,-1
           dist=(dble(ibin-1)+0.5d0)*dr
           c%phi(ipeak,ibin)=c%phi(ipeak,ibin+1)+grav*c%mass_bin(ipeak,ibin)/dist**2*dr
        end do
     endif
  end do

  end associate

end subroutine particle_potential
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine particle_unbind(s,p)
  use amr_parameters, only: ndim, nbin, twotondim
  use ramses_commons, only: ramses_t
  use pm_commons, only: part_t
  use cache_commons
  use cache
  implicit none
  type(ramses_t)::s
  type(part_t)::p
  !------------------------------------------------------------------
  ! This routine unbinds particle hierarchically between children
  ! clumps and parent clumps in the saddle point merging hierarchy.
  ! The halo-patch is used as a garbage colector.
  ! Written by Romain Teyssier (mini-ramses version in June 2024).
  !------------------------------------------------------------------
  type(msg_unbind_clump)::dummy_unbind_clump
  integer::i,ipart,ind,idim,ibin,ilevel
  integer(kind=8)::global_peak_id
  integer::ipeak
  real(kind=8)::pi,rad,vel,bound
  real(kind=8),dimension(1:ndim)::xpart,vpart

  associate(r=>s%r,g=>s%g,m=>s%m,c=>s%c,mdl=>s%mdl)

  ! Constants
  pi=ACOS(-1.0D0)

  ! Sort particles according to global peak id
  do ipart=1,p%npart
     p%sortp(ipart)=ipart
     p%workp(ipart)=p%pid(ipart)
  end do
  call quick_sort_int_int(p%workp(1),p%sortp(1),p%npart)

  ! Loop over merging hierarchy levels
  do ilevel=0,c%merge_levelmax
     ! Demote particles to parent clump if unbound
     call open_cache_clump(mdl, c, storage_size(dummy_unbind_clump)/32, &
          pack=pack_fetch_unbind, unpack=unpack_fetch_unbind)
     do i=1+p%norphan_peak,p%npart
        ! Get peak id
        ipart=p%sortp(i)
        global_peak_id=p%workp(i)
        if(global_peak_id > 0)then !ignore orphans
            call get_peak(s,global_peak_id,ipeak,flush_cache=.false.,fetch_cache=.true.)
            if(c%lev_peak(ipeak)==ilevel)then
                ! Get clump tidal radius
                rad=c%clump_rad(ipeak)
                ! Compute total energy
                bound=total_energy(dble(p%xp(ipart,1:ndim)),c%peak_pos(ipeak,1:ndim), &
                    &             dble(p%vp(ipart,1:ndim)),c%peak_vel(ipeak,1:ndim), &
                    &             c%phi(ipeak,1:nbin),rad,r%box_size,r%periodic)
                ! If unbound, assign to next peak in hierarchy
                if(bound.GE.0d0.or.c%clump_mass(ipeak).LE.c%mass_threshold)then
                    ! if central, unbind particles to the void
                    if(global_peak_id==c%new_peak(ipeak))then
                        p%workp(i)=0
                        p%pid(ipart)=0
                    else ! if satellite, then offer unbound particles to the next peak in the hierarchy
                        p%workp(i)=c%new_peak(ipeak)
                        p%pid(ipart)=c%new_peak(ipeak)
                    endif
                endif
            endif
        endif
     end do
     call close_cache(mdl)
  end do
  ! End loop over levels

  end associate

end subroutine particle_unbind
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_fetch_unbind(c,local_peak_id,msg_size,msg_array)
  use amr_commons, only: ndim,nbin
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_unbind_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_unbind_clump)::msg

  msg%lev=c%lev_peak(local_peak_id)
  msg%ind=c%new_peak(local_peak_id)
  msg%rad=c%clump_rad(local_peak_id)
  msg%mass=c%clump_mass(local_peak_id)
  msg%vel=c%peak_vel(local_peak_id,1:ndim)
  msg%pos=c%peak_pos(local_peak_id,1:ndim)
  msg%mbin=c%phi(local_peak_id,1:nbin)

  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_unbind
!################################################################
!################################################################
!################################################################
!################################################################
subroutine unpack_fetch_unbind(c,local_peak_id,msg_size,msg_array)
  use amr_commons, only: ndim, nbin
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_unbind_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_unbind_clump)::msg

  msg=transfer(msg_array,msg)

  c%lev_peak(local_peak_id)=msg%lev
  c%new_peak(local_peak_id)=msg%ind
  c%clump_rad(local_peak_id)=msg%rad
  c%clump_mass(local_peak_id)=msg%mass
  c%peak_vel(local_peak_id,1:ndim)=msg%vel
  c%peak_pos(local_peak_id,1:ndim)=msg%pos
  c%phi(local_peak_id,1:nbin)=msg%mbin

end subroutine unpack_fetch_unbind
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine mass_profile(s,p)
  use amr_parameters, only: ndim, nbin, twotondim
  use ramses_commons, only: ramses_t
  use pm_commons, only: part_t
  use cache_commons
  use cache
  implicit none
  type(ramses_t)::s
  type(part_t)::p
  !------------------------------------------------------------------
  ! This routine compute the cumulative mass profiles of input
  ! particles p around each subhalo. The particle peak id must
  ! have been computed before by the unbinding routine.
  ! Written by Romain Teyssier (mini-ramses version in June 2024).
  !------------------------------------------------------------------
  type(msg_mbin_clump)::dummy_mbin_clump
  integer::i,ipart,ind,idim,ibin,ilevel
  integer(kind=8)::global_peak_id
  integer::ipeak,jpeak
  real(kind=8)::pi,rad,dist,dr
  real(kind=8),dimension(1:ndim)::xpart,vpart

  associate(r=>s%r,g=>s%g,m=>s%m,c=>s%c,mdl=>s%mdl)

  ! Constants
  pi=ACOS(-1.0D0)

  !-------------------------------------------
  ! Sort particles according to global peak id
  !-------------------------------------------
  do ipart=1,p%npart
     p%sortp(ipart)=ipart
     p%workp(ipart)=p%pid(ipart)
  end do
  call quick_sort_int_int(p%workp(1),p%sortp(1),p%npart)

  !--------------------------------------------------------
  ! Compute particle mass profile around their central peak
  !--------------------------------------------------------
  call open_cache_clump(mdl, c, storage_size(dummy_mbin_clump)/32, &
       pack=pack_fetch_mbin, unpack=unpack_fetch_mbin, &
       init=init_flush_mbin, flush=pack_flush_mbin, combine=unpack_flush_mbin)
  c%mass_bin=0d0
  c%npart=0
  do i=1+p%norphan_peak,p%npart
     ! Get global peak id
     ipart=p%sortp(i)
     global_peak_id=p%workp(i)
     if (global_peak_id /=0 ) then
        call get_peak(s,global_peak_id,ipeak,flush_cache=.true.,fetch_cache=.true.)
        ! Get clump tidal radius
        rad=c%clump_rad(ipeak)
        ! Compute particle radius
        xpart(1:ndim)=p%xp(ipart,1:ndim)-c%peak_pos(ipeak,1:ndim)
        dist=0d0
        do idim=1,ndim
           ! In case of periodic boundaries
           if(r%periodic(idim))then
              if(xpart(idim)> r%box_size(idim)*0.5)xpart(idim)=xpart(idim)-r%box_size(idim)
              if(xpart(idim)<-r%box_size(idim)*0.5)xpart(idim)=xpart(idim)+r%box_size(idim)
           endif
           dist=dist+xpart(idim)**2
        end do
        dist=sqrt(dist)
        dr=2d0*rad/dble(nbin)
        do ibin=1,nbin
           ! We use a simple linear binning as the mass is usually propto r
           if(dist<=dble(ibin)*dr)then
              c%mass_bin(ipeak,ibin)=c%mass_bin(ipeak,ibin)+p%mp(ipart)
              c%npart(ipeak)=c%npart(ipeak)+1
              exit
           endif
        end do
     endif
  end do
  call close_cache(mdl)

  !------------------------
  ! Compute cumulative mass
  !------------------------
  do ipeak=1,c%npeak
     if(   c%clump_mass(ipeak) > c%mass_threshold.AND. &
          & c%relevance(ipeak) > c%relevance_threshold)then
        do ibin=1,nbin-1
           c%mass_bin(ipeak,ibin+1)=c%mass_bin(ipeak,ibin+1)+c%mass_bin(ipeak,ibin)
        end do
     endif
  end do

  end associate

end subroutine mass_profile
!################################################################
!################################################################
!################################################################
!################################################################
function cmp_distance(x1,x2,v1,v2,radius,velocity,box_size,periodic)
  use amr_parameters, only: ndim
  real(kind=8),dimension(1:ndim)::x1,x2,v1,v2
  real(kind=8)::radius,velocity
  real(kind=8),dimension(1:3)::box_size
  real(kind=8)::cmp_distance
  logical,dimension(1:3)::periodic
  !-----------------------------------------------------------
  ! This function computes the phase-space Euclidian distance
  !-----------------------------------------------------------
  integer::idim
  real(kind=8)::xdist,vdist
  real(kind=8),dimension(1:ndim)::xpart,vpart
  xdist=0d0
  xpart(1:ndim)=x1(1:ndim)-x2(1:ndim)
  ! In case of periodic boundaries
  do idim=1,ndim
     if(periodic(idim))then
        if(xpart(idim)> box_size(idim)*0.5)xpart(idim)=xpart(idim)-box_size(idim)
        if(xpart(idim)<-box_size(idim)*0.5)xpart(idim)=xpart(idim)+box_size(idim)
     endif
     xdist=xdist+xpart(idim)**2
  end do
  ! Rescale distance in configuration space
  xdist=xdist/radius**2
  ! Compute Euclidian distance in velocity space
  vdist=0d0
  vpart(1:ndim)=v1(1:ndim)-v2(1:ndim)
  do idim=1,ndim
     vdist=vdist+vpart(idim)**2
  end do
  ! Rescale distance in velocity space
  vdist=vdist/velocity**2
  cmp_distance=xdist+vdist
end function cmp_distance
!################################################################
!################################################################
!################################################################
!################################################################
function total_energy(x1,x2,v1,v2,phi,radius,box_size,periodic)
  use amr_parameters, only: ndim, nbin
  real(kind=8),dimension(1:ndim)::x1,x2,v1,v2
  real(kind=8),dimension(1:nbin)::phi
  real(kind=8)::radius
  real(kind=8),dimension(1:3)::box_size
  real(kind=8)::total_energy
  logical,dimension(1:3)::periodic
  !-----------------------------------------------------------
  ! This function computes the phase-space Euclidian distance
  !-----------------------------------------------------------
  integer::idim,ibin,ileft,iright
  real(kind=8)::epot,ekin,xbin,dr,r,r2
  real(kind=8),dimension(1:ndim)::xpart,vpart
  r2=0d0
  xpart(1:ndim)=x1(1:ndim)-x2(1:ndim)
  do idim=1,ndim  ! In case of periodic boundaries
     if(periodic(idim))then
        if(xpart(idim)> box_size(idim)*0.5)xpart(idim)=xpart(idim)-box_size(idim)
        if(xpart(idim)<-box_size(idim)*0.5)xpart(idim)=xpart(idim)+box_size(idim)
     endif
     r2=r2+xpart(idim)**2
  end do
  r=sqrt(r2)
  dr=2d0*radius/dble(nbin)
  xbin=min(dble(r/dr),dble(nbin))
  ibin=int(xbin)
  ileft=ibin+1
  iright=ileft+1
  epot=0d0
  if(ileft<nbin)then
     epot=phi(ileft)*(ileft-xbin)+phi(iright)*(xbin-ibin)
  elseif(ileft==nbin)then
     epot=phi(ileft)*(ileft-xbin)
  endif
  ekin=0d0
  vpart(1:ndim)=v1(1:ndim)-v2(1:ndim)
  do idim=1,ndim
     ekin=ekin+0.5d0*vpart(idim)**2
  end do
  total_energy=ekin-epot
end function total_energy
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_fetch_mbin(c,local_peak_id,msg_size,msg_array)
  use amr_commons, only: ndim
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_mbin_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_mbin_clump)::msg

  msg%pos(1:ndim)=c%peak_pos(local_peak_id,1:ndim)
  msg%rad=c%clump_rad(local_peak_id)
  msg%lev=c%lev_peak(local_peak_id)
  msg%ind=c%new_peak(local_peak_id)

  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_mbin
!################################################################
!################################################################
!################################################################
!################################################################
subroutine unpack_fetch_mbin(c,local_peak_id,msg_size,msg_array)
  use amr_commons, only: ndim
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_mbin_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_mbin_clump)::msg

  msg=transfer(msg_array,msg)

  c%peak_pos(local_peak_id,1:ndim)=msg%pos(1:ndim)
  c%clump_rad(local_peak_id)=msg%rad
  c%lev_peak(local_peak_id)=msg%lev
  c%new_peak(local_peak_id)=msg%ind

end subroutine unpack_fetch_mbin
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_flush_mbin(c,local_peak_id)
  use amr_commons, only: nbin
  use clfind_commons, only: clump_t
  type(clump_t)::c
  integer::local_peak_id

  c%mass_bin(local_peak_id,1:nbin)=0d0
  c%npart(local_peak_id)=0

end subroutine init_flush_mbin
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_flush_mbin(c,local_peak_id,msg_size,msg_array)
  use amr_commons, only: nbin
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_mbin_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_mbin_clump)::msg

  msg%mbin(1:nbin)=c%mass_bin(local_peak_id,1:nbin)
  msg%npart=c%npart(local_peak_id)

  msg_array=transfer(msg,msg_array)

end subroutine pack_flush_mbin
!################################################################
!################################################################
!################################################################
!################################################################
subroutine unpack_flush_mbin(c,local_peak_id,msg_size,msg_array)
  use amr_commons, only: nbin
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_mbin_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_mbin_clump)::msg

  msg=transfer(msg_array,msg)

  c%mass_bin(local_peak_id,1:nbin)=c%mass_bin(local_peak_id,1:nbin)+msg%mbin(1:nbin)
  c%npart(local_peak_id)=c%npart(local_peak_id)+msg%npart

end subroutine unpack_flush_mbin
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine particle_peak_id(s,p)
  use amr_parameters, only: ndim, nbin, twotondim
  use ramses_commons, only: ramses_t
  use pm_commons, only: part_t
  use nbors_utils
  use cache_commons
  use cache
  use marshal, only: pack_fetch_flag2, unpack_fetch_flag2
  use hilbert
  implicit none
  type(ramses_t)::s
  type(part_t)::p
  !-------------------------------------------------------------------
  ! This routine reads from the grid peak map (flag2) the peak id
  ! of the input particle object.
  ! Written by Romain Teyssier (mini-ramses version in June 2024).
  !-------------------------------------------------------------------
  integer,dimension(1:ndim)::ckey
  integer(kind=8),dimension(0:ndim)::hash_cell
  integer::i,ipart,icell,igrid,ind,idim,ibin,ilevel
  integer(kind=8)::global_peak_id
  integer::local_peak_id,ipeak,jpeak,merge_to
  integer::halo_nr,peak_nr
  real(kind=8)::dx_loc,rmin,rmax,xx
  type(msg_int4)::dummy_int4
  logical::ok_level,ok_leaf
  integer::no_peak

  associate(r=>s%r,g=>s%g,m=>s%m,c=>s%c,mdl=>s%mdl)

  ! Open cache for array flag1 (fetch)
  call open_cache(mdl, m, pack_size=storage_size(dummy_int4)/32, &
       pack=pack_fetch_flag2, unpack=unpack_fetch_flag2)

  ! Loop over particles
  no_peak=0
  p%pid=0
  do ilevel=r%levelmin,r%nlevelmax

     ! Mesh spacing in that level
     dx_loc=r%boxlen/2**ilevel

     do ipart=p%headp(ilevel),p%tailp(ilevel)

        ok_level=.true.

        ! Find parent cell at level ilevel
        do idim=1,ndim
           ckey(idim)=int(p%xp(ipart,idim)/dx_loc)
        end do

        ! Get parent cell at level ilevel using cache
        hash_cell(0)=ilevel+1
        hash_cell(1:ndim)=ckey(1:ndim)
        call get_parent_cell(s,hash_cell,igrid,icell,flush_cache=.false.,fetch_cache=.true.)

        ! Read flag2 value
        global_peak_id=0
        if(igrid>0)global_peak_id=m%flag2(icell,igrid)
        if(global_peak_id==0)no_peak=no_peak+1

        ! Store global peak id in pid array
        p%pid(ipart)=global_peak_id

     end do
     ! End loop over particles
  end do
  ! End loop over levels

  call close_cache(mdl)

  p%norphan_peak=no_peak

  end associate

end subroutine particle_peak_id
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine particle_halo_id(s,p)
  use amr_parameters, only: ndim, nbin, twotondim
  use ramses_commons, only: ramses_t
  use pm_commons, only: part_t
  use nbors_utils
  use cache_commons
  use cache
  use boundaries, only: init_bound_flag
  use marshal, only: pack_fetch_flag, unpack_fetch_flag
  use hilbert
  implicit none
  type(ramses_t)::s
  type(part_t)::p
  !------------------------------------------------------------------
  ! This routine reads from the grid peak map (flag1) the peak id
  ! of the input particle object.
  ! Written by Romain Teyssier (mini-ramses version in June 2024).
  !------------------------------------------------------------------
  integer,dimension(1:ndim)::ckey
  integer(kind=8),dimension(0:ndim)::hash_cell
  integer::i,ipart,icell,igrid,ind,idim,ibin,ilevel
  integer(kind=8)::global_halo_id
  integer::local_peak_id,ipeak,jpeak,merge_to
  integer::halo_nr,peak_nr
  real(kind=8)::dx_loc,rmin,rmax,xx
  type(msg_int4)::dummy_int4
  logical::ok_level,ok_leaf
  integer::no_halo

  associate(r=>s%r,g=>s%g,m=>s%m,c=>s%c,mdl=>s%mdl)

  ! Open cache for array flag1 (fetch)
  call open_cache(mdl, m, pack_size=storage_size(dummy_int4)/32, &
       pack=pack_fetch_flag, unpack=unpack_fetch_flag, bound=init_bound_flag)

  ! Loop over particles
  no_halo=0
  p%hid=0
  do ilevel=r%levelmin,r%nlevelmax

     ! Mesh spacing in that level
     dx_loc=r%boxlen/2**ilevel

     do ipart=p%headp(ilevel),p%tailp(ilevel)

        ok_level=.true.

        ! Find parent cell at level ilevel
        do idim=1,ndim
           ckey(idim)=int(p%xp(ipart,idim)/dx_loc)
        end do

        ! Get parent cell at level ilevel using cache
        hash_cell(0)=ilevel+1
        hash_cell(1:ndim)=ckey(1:ndim)
        call get_parent_cell(s,hash_cell,igrid,icell,flush_cache=.false.,fetch_cache=.true.)

        ! Read flag1 value
        global_halo_id=0
        if(igrid>0)global_halo_id=m%flag1(icell,igrid)
        if(global_halo_id==0)no_halo=no_halo+1

        ! Store global halo id in hid array
        p%hid(ipart)=global_halo_id

     end do
     ! End loop over particles
  end do
  ! End loop over levels

  call close_cache(mdl)

  p%norphan_halo=no_halo

  end associate

end subroutine particle_halo_id
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
end module clump_merger_module
