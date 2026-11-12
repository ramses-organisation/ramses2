module refine_utils
#ifdef _CUDA
  use gpu_runner, only: gpu_refine
#elif defined(_METAL)
  use metal_runner, only: metal_refine
#endif
  type out_refine_fine_t
    integer::make,kill
  end type out_refine_fine_t
contains
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_refine_fine(pst,ilevel)
  use ramses_commons, only: pst_t
  use init_refine_basegrid_module, only:r_noct_max,r_noct_min,r_noct_tot,r_noct_used_max
  use load_balance_module, only: m_load_balance, r_balance_part
  use mdl_module
  implicit none
  type(pst_t)::pst
  integer::ilevel
  !--------------------------------------------------------------------
  ! This routine is the master procedure to refine the AMR grid
  ! from level ilevel to nlevelmax.
  !--------------------------------------------------------------------
  integer::ilev,dummy
  double precision::ttend,ttstart
  type(out_refine_fine_t)::out_refine_fine

  associate(s=>pst%s)

  if(ilevel==s%r%nlevelmax)return
  if(s%m%noct_tot(ilevel)==0)return

  if(s%r%verbose)write(*,111)ilevel
111 format('   Entering refine_fine for level ',I2)

  ! Create new octs and destroy unecessary octs
  call r_refine_fine(pst,ilevel,1,out_refine_fine,2)

  if(s%r%verbose)write(*,112)out_refine_fine%make
112 format('   ==> Make ',i7,' sub-grids')

  if(s%r%verbose)write(*,113)out_refine_fine%kill
113 format('   ==> Kill ',i7,' sub-grids')

  ! Get total, min and max grid count (only in master)
  do ilev=ilevel+1,s%r%nlevelmax
     call r_noct_tot(pst,ilev,1,s%m%noct_tot(ilev),2)
     call r_noct_min(pst,ilev,1,s%m%noct_min(ilev),1)
     call r_noct_max(pst,ilev,1,s%m%noct_max(ilev),1)
  end do

  ! Get maximum used memory (only in master)
  call r_noct_used_max(pst,ilevel,1,s%m%noct_used_max,1)

  ! Load balance all levels across cpus
  call m_load_balance(pst,ilevel)

  ! Find clean and dirty octs
!  call r_clean_dirty(pst,ilevel,1)

  ! Get total, min and max grid count (only in master).
  do ilev=ilevel+1,s%r%nlevelmax
     call r_noct_tot(pst,ilev,1,s%m%noct_tot(ilev),2)
     call r_noct_min(pst,ilev,1,s%m%noct_min(ilev),1)
     call r_noct_max(pst,ilev,1,s%m%noct_max(ilev),1)
  end do

  ! Get maximum used memory (only in master)
  call r_noct_used_max(pst,ilevel,1,s%m%noct_used_max,1)

  ! Balance particles across cpus every nremap coarse steps (never if nremap<=0)
  if(mdl_threads(s%mdl)>1.AND.ilevel==s%r%levelmin.AND.s%r%nremap>0)then
     if(s%r%pic.AND.mod(s%g%nstep_coarse,s%r%nremap)==0.AND.s%g%nstep_coarse>0)then
        write(*,*)'Load balancing particle distribution'
        ttstart = mdl_wtime(s%mdl)
        call r_balance_part(pst,ilevel,1,dummy,0)
        ttend = mdl_wtime(s%mdl)
        print '(A,F0.7)',' Time elapsed load balancing: ',ttend-ttstart
     endif
  endif

  end associate

end subroutine m_refine_fine
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_refine_fine(pst,ilevel,input_size,output,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer::ilevel
  type(out_refine_fine_t)::output

  type(out_refine_fine_t)::next_output
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_REFINE_FINE,pst%iUpper+1,input_size,output_size,ilevel)
     call r_refine_fine(pst%pLower,ilevel,input_size,output,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_output)
     output%make = output%make + next_output%make
     output%kill = output%kill + next_output%kill
  else
#ifdef _CUDA
     if(pst%s%m%data_on_device)then
        call gpu_refine(pst%s,ilevel,output%make,output%kill)
     else
        call refine_fine(pst%s,ilevel,output%make,output%kill)
     endif
#elif defined(_METAL)
     if(pst%s%m%data_on_device)then
        call metal_refine(pst%s,ilevel,output%make,output%kill)
     else
        call refine_fine(pst%s,ilevel,output%make,output%kill)
     endif
#else
     call refine_fine(pst%s,ilevel,output%make,output%kill)
#endif
  endif

end subroutine r_refine_fine
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine refine_fine(s,ilevel,ncreate,nkill)
  use amr_parameters, only: ndim, nhilbert, twotondim
  use init_refine_basegrid_module, only: check_subgrid
  use flag_utils, only: ensure_ref_rules
  use ramses_commons, only: ramses_t
  use oct_commons, only: oct
  use marshal, only: pack_fetch_refine, unpack_fetch_refine, pack_fetch_flag, unpack_fetch_flag
  use boundaries, only: init_bound_refine
  use cache_commons
  use cache
  use hash
  use hilbert
  use call_back, only: cache_f
  use nbors_utils
  use hydro_parameters, only: nvar, nion
  use rt_parameters, only: nrtvar
  use cr_parameters, only: ncruvar
#ifndef DO_RTZ
  use init_xion_module, only: calc_equilibrium_xion
#endif
  implicit none
  type(ramses_t)::s
  integer::ilevel
  integer::ncreate,nkill
  !---------------------------------------------------------
  ! This routine refines cells at level ilevel if cells
  ! are flagged for refinement and are not already refined.
  ! This routine destroys refinements for cells that are
  ! not flagged for refinement and are already refined.
  ! For single time-stepping, numerical rules are
  ! automatically satisfied. For adaptive time-stepping,
  ! numerical rules are checked before refining any cell.
  !---------------------------------------------------------
  integer::icell,i,j,ibit,ibucket,ilev,ind,inew,ioct
  integer::noct_zero,head_zero,indx_zero
  integer::skip_bit,ikey,true_level
  integer::ind_cell,igrid,ind_parent
  integer(kind=8),dimension(0:ndim)::hash_key
  integer,dimension(:),allocatable::noct_level,head_level,indx_level
  integer,dimension(:),allocatable::swap_table,swap_tmp
  integer,dimension(0:twotondim-1)::bucket_count,bucket_offset
  logical::ok
  type(msg_large_realdp)::dummy_large_realdp
  type(msg_int4)::dummy_int4
  integer,dimension(1:twotondim)::flag1_tmp,flag2_tmp
#ifdef HYDRO
  real(kind=8),dimension(1:nion)::xion
  real(kind=8),dimension(1:nvar)::uold
  real(kind=8),dimension(1:twotondim,1:nvar)::uold_tmp
#endif
#ifdef MHD
  real(kind=8),dimension(1:6)::bold
  real(kind=8),dimension(1:twotondim,1:6)::bold_tmp
#endif
#ifdef DO_RT
  real(kind=8),dimension(1:nrtvar)::rtuold
  real(kind=8),dimension(1:twotondim,1:nrtvar)::rtuold_tmp
#endif
#ifdef DO_CR
  real(kind=8),dimension(1:nrtvar)::cruold
  real(kind=8),dimension(1:twotondim,1:ncruvar)::cruold_tmp
#endif
#ifdef GRAV
  real(kind=8),dimension(1:twotondim,1:3)::f_tmp
  real(kind=8),dimension(1:twotondim)::phi_tmp,phi_old_tmp
#endif
  type(oct)::oct_tmp

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

  !---------------------------------------------------
  ! Step 1: if a cell is flagged for refinement and
  ! if it is not already refined, create its son grid.
  !---------------------------------------------------
  g%ncreate=0
  m%ifree=m%noct_used+1
  do ilev=ilevel,r%nlevelmax-1

     ! Finer level flags may rely on coarser buffer cells cancelled by
     ! ensure_ref_rules (adaptive time steps): check them now that the
     ! octs of level ilev created in this pass exist
     if(ilev>ilevel)call ensure_ref_rules(s,ilev)

     call open_cache(mdl, m, pack_size=storage_size(dummy_large_realdp)/32, &
          pack=pack_fetch_refine, unpack=unpack_fetch_refine, &
          flush=pack_flush_refine, combine=unpack_flush_refine, &
          bound=init_bound_refine)

     do ioct=m%head(ilev),m%tail(ilev)
        do ind=1,twotondim
           ok   = m%flag1(ind,ioct)==1 .and. &
                & .not.m%grid(ioct)%refined(ind)
           if(ok)then
              ind_parent=ioct
              ind_cell=ind
              call make_new_oct(s,ind_parent,ind_cell,ilev+1)
              g%ncreate=g%ncreate+1
           endif
        end do
     end do

     call close_cache(mdl)

     ! Set status of parent cell to "refined"
     do ioct=m%head(ilev),m%tail(ilev)
        do ind=1,twotondim
           ok   = m%flag1(ind,ioct)==1 .and. &
                & .not.m%grid(ioct)%refined(ind)
           if(ok)then
              m%grid(ioct)%refined(ind)=.true.
           endif
        end do
     end do

  end do
  ncreate=g%ncreate

#ifndef DO_RTZ
#ifdef HYDRO
  if(r%neq_chem .and. r%upload_equilibrium_x) then
     ! Enforce equilibrium on ionization states when derefining, to
     ! prevent unnatural values (e.g when merging hot and cold cells).
     ! Skip this during grid initialization (i.e. nstep_coarse=0)
     do ilev=ilevel,r%nlevelmax-1
        do ioct=m%head(ilev),m%tail(ilev)
           do ind=1,twotondim
              ok   = m%flag1(ind,ioct)==0 .and. &
                   & m%grid(ioct)%refined(ind)
              if(ok)then
                 uold=m%uold(ind,1:nvar,ioct)
#ifdef MHD
                 bold=m%bold(ind,1:6,ioct)
#endif
#ifdef DO_RT
                 rtuold=m%rtuold(ind,1:nrtvar,ioct)
#endif
                 call calc_equilibrium_xion(s, uold, &
#ifdef MHD
                      & bold, &
#endif
#ifdef DO_RT
                      & rtuold, &
#endif
                      & ilev, xion)
                 m%uold(ind,r%iIons:r%iIons+nion-1,ioct)=xion*m%uold(ind,1,ioct)
              endif
           end do
        end do
     end do
  endif
#endif
#endif

  !----------------------------------------------------------
  ! Step 2: if the parent cell is not flagged for refinement,
  ! but it is refined, then destroy the child grid.
  !----------------------------------------------------------
  g%nkill=0
  do ilev=r%nlevelmax,ilevel+1,-1

     call open_cache(mdl, m, pack_size=storage_size(dummy_int4)/32, &
          pack=pack_fetch_flag, unpack=unpack_fetch_flag, &
          init=init_flush_derefine, flush=pack_flush_derefine, combine=unpack_flush_derefine)

     hash_key(0)=ilev
     do ioct=m%head(ilev),m%tail(ilev)
        hash_key(1:ndim)=m%grid(ioct)%ckey(1:ndim)
        ! Get parent cell using a read-write cache
        call get_parent_cell(s,hash_key,igrid,icell,flush_cache=.true.,fetch_cache=.true.)
        if (igrid==0) then
          write(*,*) 'FATAL: no parent',hash_key
          stop
        endif
        ok   = m%flag1(icell,igrid)==0 .and. &
             & m%grid(igrid)%refined(icell)
        if(ok)then
           ! Set grid level to zero
           m%grid(ioct)%lev=0
           ! Set parent cell to "unrefined" status
           m%grid(igrid)%refined(icell)=.false.
           ! Free grid from hash table
           call hash_free(m%grid_dict,hash_key)
           g%nkill=g%nkill+1
        end if
     end do

     call close_cache(mdl)

  end do
  nkill=g%nkill

  !-----------------------------------------------------
  ! Step 3: sort new octs and empty slots according to
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
     if(true_level.NE.0)then
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
     if(true_level.NE.0)then
        swap_table(indx_level(true_level))=ioct
        indx_level(true_level)=indx_level(true_level)+1
     else
        swap_table(indx_zero)=ioct
        indx_zero=indx_zero+1
     end if
  end do

  !-----------------------------------------------------
  ! Step 4: sort octs level by level according to their
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
              ibucket=int(ibits(m%grid(ioct)%hkey(ikey),skip_bit,ndim),kind=4)
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
              ibucket=int(ibits(m%grid(ioct)%hkey(ikey),skip_bit,ndim),kind=4)
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
  ! Step 5: Apply permutations directly in main memory
  ! Remember: swap_table(inew)=iold means:
  ! New data at position inew COMES FROM
  ! Old data at position iold.
  !-----------------------------------------------------
  ! Perform the swap
  do j=head_level(ilevel+1),m%ifree-1
     if(j.NE.swap_table(j))then
        hash_key(0)=m%grid(j)%lev
        hash_key(1:ndim)=m%grid(j)%ckey(1:ndim)
        if(m%grid(j)%lev>0)call hash_free(m%grid_dict,hash_key)
        oct_tmp=m%grid(j)
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
        m%grid(i)=oct_tmp
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
  ! Step 6: Clean up final AMR structure
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

end subroutine refine_fine
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine pack_flush_refine(mesh,igrid,msg_size,msg_array)
  use amr_parameters, only: ndim, twotondim
  use hydro_parameters, only: nvar
  use rt_parameters, only: nrtvar
  use cr_parameters, only: ncruvar
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_large_realdp
  type(mesh_t)::mesh
  integer::igrid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind,ivar,idim
  type(msg_large_realdp)::msg

#ifdef HYDRO
  do ivar=1,nvar
     do ind=1,twotondim
        msg%realdp_hydro(ind,ivar)=mesh%uold(ind,ivar,igrid)
     end do
  end do
#endif

#ifdef MHD
  do ivar=1,6
     do ind=1,twotondim
        msg%realdp_mhd(ind,ivar)=mesh%bold(ind,ivar,igrid)
     end do
  end do
#endif

#ifdef GRAV
  do idim=1,ndim
     do ind=1,twotondim
        msg%realdp_poisson(ind,idim)=mesh%f(ind,idim,igrid)
     end do
  end do
  do ind=1,twotondim
     msg%realdp_poisson(ind,ndim+1)=mesh%phi(ind,igrid)
     msg%realdp_poisson(ind,ndim+2)=mesh%phi_old(ind,igrid)
  end do
#endif

#ifdef DO_RT
  do ivar=1,nrtvar
     do ind=1,twotondim
        msg%realdp_rt(ind,ivar)=mesh%rtuold(ind,ivar,igrid)
     end do
  end do
#endif

#ifdef DO_CR
  do ivar=1,ncruvar
     do ind=1,twotondim
        msg%realdp_cr(ind,ivar)=mesh%cruold(ind,ivar,igrid)
     end do
  end do
#endif

  msg_array=transfer(msg,msg_array)

end subroutine pack_flush_refine
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine unpack_flush_refine(mesh,igrid,msg_size,msg_array,hash_key)
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
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind,ivar,idim
  type(msg_large_realdp)::msg

  mesh%grid(igrid)%lev=hash_key(0)
  mesh%grid(igrid)%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)

  do ind=1,twotondim
     mesh%grid(igrid)%refined(ind)=.false.
  end do

#ifdef HYDRO
  do ivar=1,nvar
     do ind=1,twotondim
        mesh%uold(ind,ivar,igrid)=msg%realdp_hydro(ind,ivar)
     end do
  end do
#endif

#ifdef MHD
  do ivar=1,6
     do ind=1,twotondim
        mesh%bold(ind,ivar,igrid)=msg%realdp_mhd(ind,ivar)
     end do
  end do
#endif

#ifdef GRAV
  do idim=1,ndim
     do ind=1,twotondim
        mesh%f(ind,idim,igrid)=msg%realdp_poisson(ind,idim)
     end do
  end do
  do ind=1,twotondim
     mesh%phi(ind,igrid)=msg%realdp_poisson(ind,ndim+1)
     mesh%phi_old(ind,igrid)=msg%realdp_poisson(ind,ndim+2)
  end do
#endif

#ifdef DO_RT
  do ivar=1,nrtvar
     do ind=1,twotondim
        mesh%rtuold(ind,ivar,igrid)=msg%realdp_rt(ind,ivar)
     end do
  end do
#endif

#ifdef DO_CR
  do ivar=1,ncruvar
     do ind=1,twotondim
        mesh%cruold(ind,ivar,igrid)=msg%realdp_cr(ind,ivar)
     end do
  end do
#endif

end subroutine unpack_flush_refine
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine init_flush_derefine(mesh,igrid,hash_key)
  use amr_parameters, only: ndim, twotondim
  use amr_commons, only: mesh_t
  type(mesh_t)::mesh
  integer::igrid
  integer(kind=8),dimension(0:ndim)::hash_key

  mesh%grid(igrid)%lev=hash_key(0)
  mesh%grid(igrid)%ckey(1:ndim)=hash_key(1:ndim)
  mesh%grid(igrid)%refined(1:twotondim)=.true.

end subroutine init_flush_derefine
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine pack_flush_derefine(mesh,igrid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_int4
  type(mesh_t)::mesh
  integer::igrid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind
  type(msg_int4)::msg

  do ind=1,twotondim
     if(mesh%grid(igrid)%refined(ind))then
        msg%int4(ind)=1
     else
        msg%int4(ind)=0
     endif
  end do
  msg_array=transfer(msg,msg_array)

end subroutine pack_flush_derefine
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine unpack_flush_derefine(mesh,igrid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim, twotondim
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_int4
  type(mesh_t)::mesh
  integer::igrid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind
  type(msg_int4)::msg

  mesh%grid(igrid)%lev=hash_key(0)
  mesh%grid(igrid)%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)
  do ind=1,twotondim
     if(mesh%grid(igrid)%refined(ind))then
        if(msg%int4(ind)==0)then
           mesh%grid(igrid)%refined(ind)=.false.
        endif
     endif
  end do

end subroutine unpack_flush_derefine
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine make_new_oct(s,iparent,icell,ilevel)
  use mdl_module
  use amr_parameters, only: ndim, nhilbert, twotondim, twondim, nvector
  use hydro_parameters, only: nvar
  use ramses_commons, only: ramses_t
#ifdef DO_CR
  use cr_parameters, only: ncruvar
#endif
  use nbors_utils
  use cache_commons
  use hilbert
  use hash
#ifndef WITHOUTMPI
  use mpi
#endif
  use rt_parameters, only: nrtvar, smallnp
  implicit none
  type(ramses_t)::s
  integer::iparent,icell,ilevel
  !--------------------------------------------------------------
  ! This routine creates a children oct at level ilevel.
  ! ilevel is thus the level of the new children oct.
  ! The new oct is labelled using the new index ichild.
  ! The parent cell is labeled with the parent oct index iparent
  ! and the cell index icell (from 1 to 8).
  !--------------------------------------------------------------
  integer::idim,ivar,ind,ichild,inbor,nstride,grid_cpu
  integer(kind=8),dimension(1:nhilbert)::hk
  integer(kind=8),dimension(1:ndim)::ix
  integer(kind=8),dimension(1:ndim)::cart_key
  integer(kind=8),dimension(0:ndim)::hash_key
  integer,dimension(0:twondim)::ind_nbor
  integer,dimension(0:twondim)::igrid_nbor
  logical::ok
#ifdef HYDRO
  real(kind=8),dimension(0:twondim,1:nvar)::u1
  real(kind=8),dimension(1:twotondim,1:nvar)::u2
#endif
#ifdef MHD
  real(kind=8),dimension(0:twondim,1:6)::b1
  real(kind=8),dimension(1:twotondim,1:6)::b2
  real(kind=8),dimension(1:twondim,1:twotondim,1:6)::b3
  logical,dimension(1:twondim)::refined
  integer,dimension(1:3,1:6),save::shift=reshape(&
       & (/-1,0,0,1,0,0,0,-1,0,0,1,0,0,0,-1,0,0,1/),(/3,6/))
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer::igridn
#endif
#ifdef DO_RT
  real(kind=8),dimension(0:twondim  ,1:nrtvar)::rtu1
  real(kind=8),dimension(1:twotondim,1:nrtvar)::rtu2
#endif
#ifdef DO_CR
  real(kind=8),dimension(0:twondim  ,1:ncruvar)::cru1
  real(kind=8),dimension(1:twotondim,1:ncruvar)::cru2
#endif

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

#ifndef WITHOUTMPI
  ! If counter is good, check on incoming messages and perform actions
  if(mdl%mail_counter==32)then
     call check_mail(mdl,MPI_REQUEST_NULL)
     mdl%mail_counter=0
  endif
  mdl%mail_counter=mdl%mail_counter+1
#endif

  !=================================
  ! Create new octs into main memory
  !=================================
  ! Compute Cartesian keys of new octs
  do idim=1,ndim
     nstride=2**(idim-1)
     cart_key(idim)=2*m%grid(iparent)%ckey(idim)+MOD((icell-1)/nstride,2)
  end do
  hash_key(0)=ilevel
  hash_key(1:ndim)=cart_key(1:ndim)

  ! Compute Hilbert keys of new octs
  ix(1:ndim)=cart_key(1:ndim)
  hk(1:nhilbert)=hilbert_key(ix,ilevel-1)

  ! Check if grid sits inside processor boundaries
  if (m%domain(ilevel)%in_rank(hk)) then

     ! Set grid index to a virtual grid in local main memory
     ichild=m%ifree

     ! Insert new grid in hash table
     call hash_setp(m%grid_dict,hash_key,ichild)

     ! Go to next main memory free line
     m%ifree=m%ifree+1
     if(m%ifree.GT.m%ngridmax)then
        write(*,*)'No more free memory'
        write(*,*)'Increase ngridmax'
        call mdl_abort(mdl)
     end if

  ! Otherwise, determine parent processor and use the cache
  else

     grid_cpu=m%domain(ilevel)%get_rank(hk)

     ! If next cache line is occupied, free it.
     if(m%occupied(m%free_cache))call destage(mdl,m%ngridmax+m%free_cache)

     ! Set grid index to a virtual grid in local cache memory
     ichild=m%ngridmax+m%free_cache
     m%occupied(m%free_cache)=.true.
     m%parent_cpu(m%free_cache)=grid_cpu
     m%dirty(m%free_cache)=.true.
     m%ghost_parent_grid(m%free_cache)=0
     m%ghost_parent_cell(m%free_cache)=0

     ! Insert new grid in hash table
     call hash_setp(m%grid_dict,hash_key,ichild)

     ! Go to next free cache line
     m%free_cache=m%free_cache+1
     m%ncache=m%ncache+1
     if(m%free_cache.GT.m%ncachemax)m%free_cache=1
     if(m%ncache.GT.m%ncachemax)m%ncache=m%ncachemax

  endif

  ! Set oct properties
  m%grid(ichild)%lev=ilevel
  m%grid(ichild)%ckey(1:ndim)=int(cart_key(1:ndim),kind=4)
  m%grid(ichild)%hkey(1:nhilbert)=hk(1:nhilbert)
  m%grid(ichild)%refined(1:twotondim)=.false.

  ! Set flag arrays to zero
  m%flag1(1:twotondim,ichild)=0
  m%flag2(1:twotondim,ichild)=0

  !===================================
  ! Interpolate parent hydro variables
  !===================================
#ifdef HYDRO
  ! Get 2ndim neighboring father cells with read-only cache
  call get_twondim_nbor_parent_cell(s,hash_key,igrid_nbor,ind_nbor,flush_cache=.false.,fetch_cache=.true.)
  ok=.true.
  do inbor=0,twondim
     ok=ok.and.(igrid_nbor(inbor)>0)
  end do
  if(.not. ok)then
     write(*,*)"OUPS parent neighbors should exist"
     write(*,*)hash_key
     write(*,*)igrid_nbor(0)
     do idim=1,ndim
        write(*,*)igrid_nbor(2*idim-1)
        write(*,*)igrid_nbor(2*idim)
     end do
     call mdl_abort(mdl)
  endif
  ! Store parent cell hydro variables
  do inbor=0,twondim
     do ivar=1,nvar
        u1(inbor,ivar)=m%uold(ind_nbor(inbor),ivar,igrid_nbor(inbor))
     end do
  end do
#ifdef DO_RT
  do inbor=0,twondim
     do ivar=1,nrtvar
        rtu1(inbor,ivar)=m%rtuold(ind_nbor(inbor),ivar,igrid_nbor(inbor))
     end do
  end do
#endif
#ifdef DO_CR
  do inbor=0,twondim
     do ivar=1,ncruvar
        cru1(inbor,ivar)=m%cruold(ind_nbor(inbor),ivar,igrid_nbor(inbor))
     end do
  end do
#endif
#ifdef MHD
  ! Store parent cell MHD variables
  do inbor=0,twondim
     do ivar=1,6
        b1(inbor,ivar)=m%bold(ind_nbor(inbor),ivar,igrid_nbor(inbor))
     end do
  end do
  ! Get neighboring children grids
  hash_nbor(0)=ilevel
  do inbor=1,twondim
     hash_nbor(1:ndim)=hash_key(1:ndim)+shift(1:ndim,inbor)
     ! Periodic boundary conditions
     do idim=1,ndim
        if(r%periodic(idim))then
           if(hash_nbor(idim)< m%box_ckey_min(idim,ilevel))hash_nbor(idim)=m%box_ckey_max(idim,ilevel)-1
           if(hash_nbor(idim)>=m%box_ckey_max(idim,ilevel))hash_nbor(idim)=m%box_ckey_min(idim,ilevel)
        endif
     enddo
     igridn=0
     if(m%grid(igrid_nbor(inbor))%refined(ind_nbor(inbor)))then
        call get_grid(s,hash_nbor,igridn,flush_cache=.false.,fetch_cache=.true.)
     endif
     refined(inbor)=(igridn>0)
     if(refined(inbor))then
        do ind=1,twotondim
           do ivar=1,6
              b3(inbor,ind,ivar)=m%bold(ind,ivar,igridn)
           end do
        end do
     endif
  end do
  ! Interpolate using MHD variables
  call interpol_mhd(u1,u2,b1,b2,b3,refined,r%interpol_var,r%interpol_type,r%smallr)
#else
  ! Interpolate using hydro variables
  call interpol_hydro(u1,u2,r%interpol_var,r%interpol_type,r%smallr)
#endif
  ! Store children cell hydro variables
  do ivar=1,nvar
     do ind=1,twotondim
        m%uold(ind,ivar,ichild)=u2(ind,ivar)
     end do
  end do
#ifdef MHD
  ! Store children cell MHD variables
  do ivar=1,6
     do ind=1,twotondim
        m%bold(ind,ivar,ichild)=b2(ind,ivar)
     end do
  end do
#endif
#ifdef DO_RT
  ! Interpolate using rt variables
  call interpol_rt(rtu1,rtu2,r%interpol_type)
  ! Store children cell rt variables
  do ivar=1,nrtvar
     do ind=1,twotondim
        m%rtuold(ind,ivar,ichild)=rtu2(ind,ivar)
        ! Rescale according to speed of light difference
        if(mod(ivar,ndim+1).eq.1)then
           m%rtuold(ind,ivar,ichild) = m%rtuold(ind,ivar,ichild) * g%rt_c(ilevel-1)/g%rt_c(ilevel)
        end if
     end do
  end do
#endif
#ifdef DO_CR
  ! Interpolate using CR variables
  call interpol_cr(cru1,cru2,r%interpol_type)
  ! Store children cell CR variables
  do ivar=1,ncruvar
     do ind=1,twotondim
        m%cruold(ind,ivar,ichild)=cru2(ind,ivar)
     end do
  end do
#endif
  do inbor=1,twondim
     call unlock_cache(m,igrid_nbor(inbor))
  end do
#endif

  !================================
  ! Inject parent gravity variables
  !================================
#ifdef GRAV
  do ind=1,twotondim
     m%f(ind,1:ndim,ichild)=m%f(icell,1:ndim,iparent)
     m%phi(ind,ichild)=m%phi(icell,iparent)
     m%phi_old(ind,ichild)=m%phi_old(icell,iparent)
  end do
#endif

  end associate

end subroutine make_new_oct
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_clean_dirty(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::ilevel

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_CLEAN_DIRTY,pst%iUpper+1,input_size,0,ilevel)
     call r_clean_dirty(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call clean_dirty(pst%s,ilevel)
  endif

end subroutine r_clean_dirty
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine clean_dirty(s,ilevel)
  use amr_parameters, only: ndim, twotondim, nhilbert
  use ramses_commons, only: ramses_t
  use hilbert
  use hash
  implicit none
  type(ramses_t)::s
  integer::ilevel
  !-------------------------------------------------
  ! This routine identifies clean and dirty octs for
  ! optimization purposes. Dirty octs require MPI or
  ! inter-level communications. Clean octs have all
  ! their 26 neighbors around them.
  !-------------------------------------------------
  integer::ioct,ilev,i1
#if NDIM>1
  integer::j1
#endif
#if NDIM>2
  integer::k1
#endif
  logical::clean

  integer(kind=8),dimension(0:ndim)::hash_key

  associate(r=>s%r,g=>s%g,m=>s%m)

  !---------------------
  ! Clean and dirty octs
  !---------------------
  do ilev=ilevel+1,r%nlevelmax
     m%head_clean(ilev)=m%head_clean(ilev-1)+m%noct_clean(ilev-1)
     m%head_dirty(ilev)=m%head_dirty(ilev-1)+m%noct_dirty(ilev-1)
     m%noct_clean(ilev)=0
     m%noct_dirty(ilev)=0
     hash_key(0)=ilev
     do ioct=m%head(ilev),m%tail(ilev)
        clean=.true.
#if NDIM>2
        do k1=-1,1
        hash_key(3)=m%grid(ioct)%ckey(3)+k1
#endif
#if NDIM>1
        do j1=-1,1
        hash_key(2)=m%grid(ioct)%ckey(2)+j1
#endif
        do i1=-1,1
           hash_key(1)=m%grid(ioct)%ckey(1)+i1
           clean=clean.and.hash_is_clean(m%grid_dict,hash_key)
        end do
#if NDIM>1
        end do
#endif
#if NDIM>2
        end do
#endif
        if(clean)then
           m%indx_clean(m%head_clean(ilev)+m%noct_clean(ilev))=ioct
           m%noct_clean(ilev)=m%noct_clean(ilev)+1
        else
           m%indx_dirty(m%head_dirty(ilev)+m%noct_dirty(ilev))=ioct
           m%noct_dirty(ilev)=m%noct_dirty(ilev)+1
        endif
     end do
  end do

  end associate

end subroutine clean_dirty
!###############################################################
!###############################################################
!###############################################################
!###############################################################
end module refine_utils
