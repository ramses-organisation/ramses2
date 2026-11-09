module rt_upload_module
#if defined(_CUDA) && defined(DO_RT)
  use gpu_runner, only: gpu_rt_upload_fine
#endif
contains
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_rt_upload_fine(pst,ilevel)
  use ramses_commons, only: pst_t
  implicit none
  type(pst_t)::pst
  integer::ilevel
  !--------------------------------------------------------------------
  ! This routine is the master procedure to upload HYDRO variables
  ! from level ilevel+1 to ilevel (averaging down or restriction).
  !--------------------------------------------------------------------
  if(ilevel==pst%s%r%nlevelmax)return
  if(pst%s%m%noct_tot(ilevel)==0)return
  if(pst%s%m%noct_tot(ilevel+1)==0)return
  if(pst%s%r%verbose)write(*,111)ilevel
111 format(' Entering upload_fine for level',i2)

  call r_rt_upload_fine(pst,ilevel,1)

end subroutine m_rt_upload_fine
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_rt_upload_fine(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size

  integer::ilevel
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_RT_UPLOAD_FINE,pst%iUpper+1,input_size,0,ilevel)
     call r_rt_upload_fine(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
#if defined(_CUDA) && defined(DO_RT)
     ! With the mesh on the device, the restriction runs there
     if(pst%s%m%data_on_device)then
        call gpu_rt_upload_fine(pst%s,ilevel)
        return
     endif
#endif
     call rt_upload_fine(pst%s,ilevel)
  endif

end subroutine r_rt_upload_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine rt_upload_fine(s,ilevel)
  use mdl_module
  use rt_parameters, only: nrtvar
  use amr_parameters, only: ndim, twotondim
  use ramses_commons, only: ramses_t
  use nbors_utils
  use cache_commons
  use cache
  use rt_flag_module, only: pack_fetch_rt, unpack_fetch_rt
  implicit none
  type(ramses_t)::s
  integer::ilevel
  !----------------------------------------------------------------------
  ! This routine performs a restriction operation (averaging down)
  ! for the RT variables.
  !----------------------------------------------------------------------
  integer::ioct,ind,ivar,icell,igrid
  integer(kind=8),dimension(0:ndim)::hash_key
  integer,dimension(1:6,1:4)::hh
  real(kind=8)::average
  type(msg_realdp)::dummy_realdp

  hh(1,1:4)=(/1,3,5,7/)
  hh(2,1:4)=(/2,4,6,8/)
  hh(3,1:4)=(/1,2,5,6/)
  hh(4,1:4)=(/3,4,7,8/)
  hh(5,1:4)=(/1,2,3,4/)
  hh(6,1:4)=(/5,6,7,8/)

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

  ! Set conservative variable to zero in refined cells
  do ioct=m%head(ilevel),m%tail(ilevel)
     do ivar=1,nrtvar
        do ind=1,twotondim
           if(m%grid(ioct)%refined(ind))then
#ifdef DO_RT
              m%rtuold(ind,ivar,ioct)=0.0
#endif
           endif
        end do
     end do
  end do

  call open_cache(mdl, m, pack_size=storage_size(dummy_realdp)/32, &
       pack=pack_fetch_rt, unpack=unpack_fetch_rt, &
       init=init_flush_upload_rt, flush=pack_flush_upload_rt, combine=unpack_flush_upload_rt)

  ! Loop over finer level grids
  hash_key(0)=ilevel+1
  do ioct=m%head(ilevel+1),m%tail(ilevel+1)

     ! Get parent cell and grid index
     hash_key(1:ndim)=m%grid(ioct)%ckey(1:ndim)
     call get_parent_cell(s,hash_key,igrid,icell,flush_cache=.true.,fetch_cache=.false.)

     ! Average conservative variables
     do ivar=1,nrtvar
        average=0.0d0
        do ind=1,twotondim
#ifdef DO_RT
           average=average+m%rtuold(ind,ivar,ioct)
#endif
        end do
        ! Scatter result to parent cell
#ifdef DO_RT
        m%rtuold(icell,ivar,igrid)=average/dble(twotondim)
        ! Rescale according to light speed difference between levels
        if (mod(ivar,ndim+1).eq.1) &
            m%rtuold(icell,ivar,igrid) = m%rtuold(icell,ivar,igrid) * g%rt_c(ilevel+1)/g%rt_c(ilevel)
#endif
     end do

  end do

  call close_cache(mdl)

  end associate

end subroutine rt_upload_fine
!##########################################################################
!##########################################################################
!##########################################################################
!##########################################################################
subroutine init_flush_upload_rt(mesh,igrid,hash_key)
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
        mesh%rtuold(ind,ivar,igrid)=0.0d0
     end do
  end do
#endif

end subroutine init_flush_upload_rt
!##########################################################################
!##########################################################################
!##########################################################################
!##########################################################################
subroutine pack_flush_upload_rt(mesh,igrid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use rt_parameters, only: nrtvar
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_realdp
  type(mesh_t)::mesh
  integer::igrid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind,ivar
  type(msg_realdp)::msg

#ifdef DO_RT
  do ivar=1,nrtvar
     do ind=1,twotondim
        msg%realdp_rt(ind,ivar)=mesh%rtuold(ind,ivar,igrid)
     end do
  end do
#endif

  msg_array=transfer(msg,msg_array)

end subroutine pack_flush_upload_rt
!##########################################################################
!##########################################################################
!##########################################################################
!##########################################################################
subroutine unpack_flush_upload_rt(mesh,igrid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim, twotondim
  use rt_parameters, only: nrtvar
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_realdp
  type(mesh_t)::mesh
  integer::igrid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind,ivar
  type(msg_realdp)::msg

  mesh%grid(igrid)%lev=hash_key(0)
  mesh%grid(igrid)%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)

#ifdef DO_RT
  do ivar=1,nrtvar
     do ind=1,twotondim
        if(mesh%grid(igrid)%refined(ind))then
           mesh%rtuold(ind,ivar,igrid)=mesh%rtuold(ind,ivar,igrid)+msg%realdp_rt(ind,ivar)
        endif
     end do
  end do
#endif

end subroutine unpack_flush_upload_rt
!##########################################################################
!##########################################################################
!##########################################################################
!##########################################################################
end module rt_upload_module
