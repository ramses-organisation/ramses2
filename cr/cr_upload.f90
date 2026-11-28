module cr_upload_module
#if defined(_CUDA) && defined(DO_CR)
  use gpu_runner, only: gpu_cr_upload_fine
#endif
contains
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_cr_upload_fine(pst,ilevel)
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

  call r_cr_upload_fine(pst,ilevel,1)

end subroutine m_cr_upload_fine
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_cr_upload_fine(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size

  integer::ilevel
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_CR_UPLOAD_FINE,pst%iUpper+1,input_size,0,ilevel)
     call r_cr_upload_fine(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
#if defined(_CUDA) && defined(DO_CR)
     ! With the mesh on the device, the restriction runs there
     if(pst%s%m%data_on_device)then
        call gpu_cr_upload_fine(pst%s,ilevel)
        return
     endif
#endif
     call cr_upload_fine(pst%s,ilevel)
  endif

end subroutine r_cr_upload_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cr_upload_fine(s,ilevel)
  use mdl_module
  use cr_parameters, only: ncruvar
  use amr_parameters, only: ndim, twotondim
  use ramses_commons, only: ramses_t
  use nbors_utils
  use cache_commons
  use cache
  implicit none
  type(ramses_t)::s
  integer::ilevel
  !----------------------------------------------------------------------
  ! This routine performs a restriction operation (averaging down)
  ! for the CR variables.
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
     do ivar=1,ncruvar
        do ind=1,twotondim
           if(m%grid(ioct)%refined(ind))then
#ifdef DO_CR
              m%cruold(ind,ivar,ioct)=0.0
#endif
           endif
        end do
     end do
  end do

  call open_cache(mdl, m, pack_size=storage_size(dummy_realdp)/32, &
       pack=pack_fetch_cr, unpack=unpack_fetch_cr, &
       init=init_flush_upload_cr, flush=pack_flush_upload_cr, combine=unpack_flush_upload_cr)

  ! Loop over finer level grids
  hash_key(0)=ilevel+1
  do ioct=m%head(ilevel+1),m%tail(ilevel+1)

     ! Get parent cell and grid index
     hash_key(1:ndim)=m%grid(ioct)%ckey(1:ndim)
     call get_parent_cell(s,hash_key,igrid,icell,flush_cache=.true.,fetch_cache=.false.)

     ! Average conservative variables
     do ivar=1,ncruvar
        average=0.0d0
        do ind=1,twotondim
#ifdef DO_CR
           average=average+m%cruold(ind,ivar,ioct)
#endif
        end do
        ! Scatter result to parent cell
#ifdef DO_CR
        m%cruold(icell,ivar,igrid)=average/dble(twotondim)
#endif
     end do

  end do

  call close_cache(mdl)

  end associate

end subroutine cr_upload_fine
!##########################################################################
!##########################################################################
!##########################################################################
!##########################################################################
subroutine init_flush_upload_cr(mesh,igrid,hash_key)
  use amr_parameters, only: ndim, twotondim
  use cr_parameters, only: ncruvar
  use amr_commons, only: mesh_t
  type(mesh_t)::mesh
  integer::igrid
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind,ivar

  mesh%grid(igrid)%lev=hash_key(0)
  mesh%grid(igrid)%ckey(1:ndim)=hash_key(1:ndim)

  do ivar=1,ncruvar
     do ind=1,twotondim
#ifdef DO_CR
        mesh%cruold(ind,ivar,igrid)=0.0d0
#endif
     end do
  end do

end subroutine init_flush_upload_cr
!##########################################################################
!##########################################################################
!##########################################################################
!##########################################################################
subroutine pack_flush_upload_cr(mesh,igrid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use cr_parameters, only: ncruvar
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_realdp
  type(mesh_t)::mesh
  integer::igrid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind,ivar
  type(msg_realdp)::msg

  do ivar=1,ncruvar
     do ind=1,twotondim
#ifdef DO_CR
        msg%realdp_cr(ind,ivar)=mesh%cruold(ind,ivar,igrid)
#endif
     end do
  end do

  msg_array=transfer(msg,msg_array)

end subroutine pack_flush_upload_cr
!##########################################################################
!##########################################################################
!##########################################################################
!##########################################################################
subroutine unpack_flush_upload_cr(mesh,igrid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim, twotondim
  use cr_parameters, only: ncruvar
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

  do ivar=1,ncruvar
     do ind=1,twotondim
#ifdef DO_CR
        if(mesh%grid(igrid)%refined(ind))then
           mesh%cruold(ind,ivar,igrid)=mesh%cruold(ind,ivar,igrid)+msg%realdp_cr(ind,ivar)
        endif
#endif
     end do
  end do

end subroutine unpack_flush_upload_cr
!##########################################################################
!##########################################################################
!##########################################################################
!##########################################################################
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine pack_fetch_cr(mesh,igrid,msg_size,msg_array)
  use amr_parameters, only: ndim, twotondim
  use hydro_parameters, only: nvar
  use cr_parameters, only: ncruvar
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_realdp
  type(mesh_t)::mesh
  integer::igrid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind,ivar
  type(msg_realdp)::msg

  do ind=1,twotondim
     if(mesh%grid(igrid)%refined(ind))then
        msg%int4(ind)=1
     else
        msg%int4(ind)=0
     endif
  end do

#ifdef HYDRO
  do ivar=1,nvar
     do ind=1,twotondim
        msg%realdp(ind,ivar)=mesh%uold(ind,ivar,igrid)
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

end subroutine pack_fetch_cr
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine unpack_fetch_cr(mesh,igrid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim, twotondim
  use hydro_parameters, only: nvar
  use cr_parameters, only: ncruvar
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

  do ind=1,twotondim
     if(msg%int4(ind)==1)then
        mesh%grid(igrid)%refined(ind)=.true.
     else
        mesh%grid(igrid)%refined(ind)=.false.
     endif
  end do

#ifdef HYDRO
  do ivar=1,nvar
     do ind=1,twotondim
        mesh%uold(ind,ivar,igrid)=msg%realdp(ind,ivar)
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

end subroutine unpack_fetch_cr
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
end module cr_upload_module
