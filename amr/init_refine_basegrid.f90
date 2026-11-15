module init_refine_basegrid_module
contains
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_init_refine_basegrid(pst)
  use ramses_commons, only: pst_t
  use flag_utils, only: m_flag_fine
#ifdef GRAV
  use rho_fine_module, only: m_rho_fine
#endif
  use input_hydro_grafic_module, only: r_input_refmap_grafic
  implicit none
  type(pst_t)::pst
  !--------------------------------------------------------------------
  ! This routine is the master procedure to set the base grid
  ! and initialize all cell-based variables within it.
  !--------------------------------------------------------------------
  associate(r=>pst%s%r,g=>pst%s%g,m=>pst%s%m,p=>pst%s%p,mdl=>pst%s%mdl)

  if(r%verbose)write(*,*)'Entering init_refine_basegrid'

  write(*,'(" Building initial base grid at level ",I0)')r%levelmin

  ! Call recursive slave routine
  call r_init_refine_basegrid(pst,r%levelmin,1)

  ! Get total, min and max grid count (only in master).
  call r_noct_tot(pst,r%levelmin,1,m%noct_tot(r%levelmin),2)
  call r_noct_min(pst,r%levelmin,1,m%noct_min(r%levelmin),1)
  call r_noct_max(pst,r%levelmin,1,m%noct_max(r%levelmin),1)
  call r_noct_used_max(pst,r%levelmin,1,m%noct_used_max,1)

  ! Initialize hydro variables on the base grid
  if(r%hydro)call m_init_flow_fine(pst,r%levelmin)

  ! Initialize rt variables on the base grid
#ifdef DO_RT
  if(r%rt)call m_rt_init_flow_fine(pst,r%levelmin)
#endif

  ! Initialize cr variables on the base grid
#ifdef DO_CR
  if(r%cr)call m_cr_init_flow_fine(pst,r%levelmin)
#endif

  ! Initialize refinement map on the base grid
  if(pst%s%r%filetype=='grafic_zoom'.and.pst%s%r%ivar_refine==0)then
     call r_input_refmap_grafic(pst,r%levelmin,1)
  endif

  ! Compute total mass density from gas and particles on the base grid
#ifdef GRAV
  if(pst%s%r%filetype.NE.'grafic_zoom')then
     call m_rho_fine(pst,r%levelmin,0)
  endif
#endif

  ! Flag coarse level cells for refinement
  call m_flag_fine(pst,r%levelmin,2)

  end associate

end subroutine m_init_refine_basegrid
!###############################################
!###############################################
!###############################################
!###############################################
recursive subroutine r_collect_noct(pst,ilevel,input_size,noct,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer,dimension(1:input_size)::ilevel
  integer,dimension(1:output_size)::noct

  integer,dimension(1:output_size)::next_noct
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_COLLECT_NOCT,pst%iUpper+1,input_size,output_size,ilevel)
     call r_collect_noct(pst%pLower,ilevel,input_size,noct,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_noct)
     noct=noct+next_noct
  else
     noct=0
     noct(pst%s%g%myid)=pst%s%m%noct(ilevel(1))
  endif

end subroutine r_collect_noct
!###############################################
!###############################################
!###############################################
!###############################################
recursive subroutine r_noct_tot(pst,ilevel,input_size,noct_tot,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer::ilevel
  integer(kind=8)::noct_tot

  integer(kind=8)::next_noct_tot
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_NOCT_TOT,pst%iUpper+1,input_size,output_size,ilevel)
     call r_noct_tot(pst%pLower,ilevel,input_size,noct_tot,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_noct_tot)
     noct_tot=noct_tot+next_noct_tot
  else
     noct_tot=int(pst%s%m%noct(ilevel),kind=8)
  endif

end subroutine r_noct_tot
!###############################################
!###############################################
!###############################################
!###############################################
recursive subroutine r_noct_max(pst,ilevel,input_size,noct_max,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer::ilevel
  integer::noct_max

  integer::next_noct_max
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_NOCT_MAX,pst%iUpper+1,input_size,output_size,ilevel)
     call r_noct_max(pst%pLower,ilevel,input_size,noct_max,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_noct_max)
     noct_max=MAX(noct_max,next_noct_max)
  else
     noct_max=pst%s%m%noct(ilevel)
  endif

end subroutine r_noct_max
!###############################################
!###############################################
!###############################################
!###############################################
recursive subroutine r_noct_used_max(pst,ilevel,input_size,noct_used_max,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer::ilevel
  integer::noct_used_max

  integer::next_noct_used_max
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_NOCT_USED_MAX,pst%iUpper+1,input_size,output_size,ilevel)
     call r_noct_used_max(pst%pLower,ilevel,input_size,noct_used_max,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_noct_used_max)
     noct_used_max=MAX(noct_used_max,next_noct_used_max)
  else
     noct_used_max=pst%s%m%noct_used
  endif

end subroutine r_noct_used_max
!###############################################
!###############################################
!###############################################
!###############################################
recursive subroutine r_noct_min(pst,ilevel,input_size,noct_min,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer::ilevel
  integer::noct_min

  integer::next_noct_min
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_NOCT_MIN,pst%iUpper+1,input_size,output_size,ilevel)
     call r_noct_min(pst%pLower,ilevel,input_size,noct_min,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_noct_min)
     noct_min=MIN(noct_min,next_noct_min)
  else
     noct_min=pst%s%m%noct(ilevel)
  endif

end subroutine r_noct_min
!###############################################
!###############################################
!###############################################
!###############################################
recursive subroutine r_init_refine_basegrid(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::ilevel

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_INIT_REFINE_BASEGRID,pst%iUpper+1,input_size,0,ilevel)
     call r_init_refine_basegrid(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call init_refine_basegrid(pst%s,ilevel)
  endif

end subroutine r_init_refine_basegrid
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_refine_basegrid(s,ilevel)
  use mdl_module, only: mdl_abort
  use amr_parameters, only: nhilbert,ndim,twotondim
  use ramses_commons, only: ramses_t
  use hilbert
  use hash, only: hash_setp, hash_is_clean, hash_stats
  implicit none
  type(ramses_t)::s
  integer::ilevel
  !-------------------------------------------------------
  ! This routine builds a fully refined Cartesian grid
  ! at level ilevel. Always starts at levelmin.
  !-------------------------------------------------------
  logical::clean
  integer::igrid,istart,i1,j1,k1
  integer(kind=8)::ikey
  integer(kind=8),dimension(1:nhilbert)::hk
  integer(kind=8),dimension(1:ndim)::ix
  integer(kind=8),dimension(0:ndim)::hash_key

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

  ! Compute starting grid index at that level
  if(ilevel.EQ.r%levelmin)then
     istart=1
  else
     istart=m%tail(ilevel-1)+1
  endif
  hk=0
  ! New grid in current level
  igrid=istart-1

  ! Loop over the Cartesian grid in Hilbert order
  do ikey=m%domain(ilevel)%b(1,g%myid-1), m%domain(ilevel)%b(1,g%myid)-1
     ! Compute Cartesian index from Hilbert index
     hk(1)=ikey
     ix=hilbert_reverse(hk,ilevel-1)
     if(ix(1).ge.m%box_ckey_min(1,ilevel).and.ix(1).lt.m%box_ckey_max(1,ilevel))then
#if NDIM>1
     if(ix(2).ge.m%box_ckey_min(2,ilevel).and.ix(2).lt.m%box_ckey_max(2,ilevel))then
#endif
#if NDIM>2
     if(ix(3).ge.m%box_ckey_min(3,ilevel).and.ix(3).lt.m%box_ckey_max(3,ilevel))then
#endif
        ! Insert new grid in main array
        igrid=igrid+1
        if(igrid.GT.m%ngridmax)then
           write(*,*)'No more free memory'
           write(*,*)'Increase ngridmax'
           call mdl_abort(mdl)
        end if
        if(igrid==istart)m%head(ilevel)=istart
        m%tail(ilevel)=igrid
        m%noct(ilevel)=m%noct(ilevel)+1
        m%noct_used=m%noct_used+1
        m%grid(igrid)%lev=ilevel
        m%grid(igrid)%ckey(1:ndim)=int(ix(1:ndim),kind=4)
        m%grid(igrid)%hkey(1:nhilbert)=hk(1:nhilbert)
        m%grid(igrid)%refined(1:twotondim)=.false.
        ! Insert new grid in hash table
        hash_key(0)=ilevel
        hash_key(1:ndim)=ix(1:ndim)
        call hash_setp(m%grid_dict,hash_key,igrid)
     endif
#if NDIM>1
     endif
#endif
#if NDIM>2
     endif
#endif
  end do
  m%ifree=m%noct_used+1

  ! Check that the level is made of complete blocks
  call check_subgrid(s,ilevel)

!!$  !---------------------
!!$  ! Clean and dirty octs
!!$  !---------------------
!!$  m%head_clean(ilev)=1
!!$  m%head_dirty(ilev)=1
!!$  m%noct_clean(ilev)=0
!!$  m%noct_dirty(ilev)=0
!!$  hash_key(0)=ilev
!!$  do ioct=m%head(ilev),m%tail(ilev)
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
!!$        m%indx_clean(m%head_clean(ilev)+m%noct_clean(ilev))=ioct
!!$        m%noct_clean(ilev)=m%noct_clean(ilev)+1
!!$     else
!!$        m%indx_dirty(m%head_dirty(ilev)+m%noct_dirty(ilev))=ioct
!!$        m%noct_dirty(ilev)=m%noct_dirty(ilev)+1
!!$     endif
!!$  end do

  end associate

end subroutine init_refine_basegrid
!################################################################
!################################################################
!################################################################
!################################################################
subroutine check_subgrid(s,ilevel)
  use amr_parameters, only: ndim, nsubgrid
  use ramses_commons, only: ramses_t
  use mdl_module, only: mdl_abort
  implicit none
  type(ramses_t)::s
  integer::ilevel
  !--------------------------------------------------------------
  ! Check that the local octs of this level form complete blocks
  ! of nsubgrid**ndim consecutive octs. Abort otherwise.
  !--------------------------------------------------------------
  integer::igrid,iskip,nblock
  logical::ok

  associate(m=>s%m,mdl=>s%mdl)

  if(nsubgrid==1)return
  nblock=nsubgrid**ndim
  ok=mod(m%noct(ilevel),nblock)==0
  if(ok.and.m%noct(ilevel)>0)then
     do igrid=m%head(ilevel),m%tail(ilevel),nblock
        do iskip=1,nblock-1
           ok=ok.and.all(m%grid(igrid+iskip)%ckey(1:ndim)/nsubgrid==m%grid(igrid)%ckey(1:ndim)/nsubgrid)
        end do
     end do
  endif
  if(.not.ok)then
     write(*,'(" Level ",I0," is not made of complete blocks of ",I0," octs (NSUBGRID=",I0,")")')ilevel,nblock,nsubgrid
     write(*,*)'Was the restart or initial file written with a different NSUBGRID?'
     call mdl_abort(mdl)
  endif

  end associate

end subroutine check_subgrid
!################################################################
!################################################################
!################################################################
!################################################################
end module init_refine_basegrid_module
