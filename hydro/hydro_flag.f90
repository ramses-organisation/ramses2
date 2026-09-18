module hydro_flag_module
contains
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine hydro_flag(s,ilevel)
  use amr_parameters, only: ndim, twotondim, twondim
  use ramses_commons, only: ramses_t
  use hydro_parameters, only: nvar, nener
  use cache_commons
  use cache
  use nbors_utils
  use boundaries, only: init_bound_refine
  implicit none
  type(ramses_t)::s
  integer::ilevel
  ! -------------------------------------------------------------------
  ! This routine flag for refinement cells that satisfies
  ! some user-defined physical criteria at the level ilevel.
  ! -------------------------------------------------------------------
  integer,dimension(1:3,1:8),save::iii=reshape(&
       & (/0,0,0,1,0,0,0,1,0,1,1,0,0,0,1,1,0,1,0,1,1,1,1,1/),(/3,8/))
  integer,dimension(1:3,1:6),save::shift=reshape(&
       & (/-1,0,0,1,0,0,0,-1,0,0,1,0,0,0,-1,0,0,1/),(/3,6/))
  integer::igrid,ind,idim,ivar,i_nbor
  integer::igridd,igridg,indd,indg,igridp,icellp
  integer,dimension(1:twondim)::igridn,icelln
  integer(kind=8),dimension(0:ndim)::hash_key,hash_nbor
  real(kind=8),dimension(1:nvar)::uug,uum,uud
#ifdef MHD
  real(kind=8),dimension(1:6)::bbg,bbm,bbd
#endif
  logical::ok
  type(msg_realdp)::dummy_realdp
#ifdef HYDRO


  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

  if(    r%err_grad_d==-1.0.and.&
#ifdef MHD
       & r%err_grad_A==-1.0.and.&
       & r%err_grad_B==-1.0.and.&
       & r%err_grad_C==-1.0.and.&
       & r%err_grad_B2==-1.0.and.&
#endif
       & r%err_grad_p==-1.0.and.&
       & r%err_grad_u==-1.0.and.&
#if NENER>0
       & sum(r%err_grad_prad(1:NENER))==-nener .and.&
#endif
       & r%err_grad_xHI==-1.0.and.&
       & r%err_grad_xHII==-1.0)return

  hash_key(0)=ilevel+1

  call open_cache(mdl, m, pack_size=storage_size(dummy_realdp)/32, &
       pack=pack_fetch_hydro, unpack=unpack_fetch_hydro, bound=init_bound_refine)

  ! Loop over active grids
  do igrid=m%head(ilevel),m%tail(ilevel)

     ! Loop over cells
     do ind=1,twotondim

        ! Compute cell hash key
        hash_key(1:ndim)=2*m%grid(igrid)%ckey(1:ndim)+iii(1:ndim,ind)

        ! Initialize refinement to false
        ok=.false.

        ! If a neighbor cell does not exist,
        ! replace it by its father cell
        do i_nbor=1,twondim
           hash_nbor(0)=hash_key(0)
           ! Periodic boundary conditions
           do idim=1,ndim
              hash_nbor(idim)=hash_key(idim)+shift(idim,i_nbor)
              if(r%periodic(idim))then
                 if(hash_nbor(idim)< m%box_ckey_min(idim,ilevel+1))hash_nbor(idim)=m%box_ckey_max(idim,ilevel+1)-1
                 if(hash_nbor(idim)>=m%box_ckey_max(idim,ilevel+1))hash_nbor(idim)=m%box_ckey_min(idim,ilevel+1)
              endif
           enddo
           call get_parent_cell(s,hash_nbor,igridp,icellp,flush_cache=.false.,fetch_cache=.true.,lock=.true.)
           if(igridp>0)then
              igridn(i_nbor)=igridp
              icelln(i_nbor)=icellp
           else
              hash_nbor(0)=hash_nbor(0)-1
              hash_nbor(1:ndim)=hash_nbor(1:ndim)/2
              call get_parent_cell(s,hash_nbor,igridp,icellp,flush_cache=.false.,fetch_cache=.true.,lock=.true.)
              igridn(i_nbor)=igridp
              icelln(i_nbor)=icellp
           endif
        end do

        ! Loop over dimensions
        do idim=1,ndim
           ! Gather hydro variables
           do ivar=1,nvar
              indg=icelln(2*idim-1)
              indd=icelln(2*idim  )
              igridg=igridn(2*idim-1)
              igridd=igridn(2*idim  )
              uug(ivar)=m%uold(indg,ivar,igridg)
              uum(ivar)=m%uold(ind,ivar,igrid)
              uud(ivar)=m%uold(indd,ivar,igridd)
           end do
#ifdef MHD
           ! Gather MHD variables
           do ivar=1,6
              indg=icelln(2*idim-1)
              indd=icelln(2*idim  )
              igridg=igridn(2*idim-1)
              igridd=igridn(2*idim  )
              bbg(ivar)=m%bold(indg,ivar,igridg)
              bbm(ivar)=m%bold(ind,ivar,igrid)
              bbd(ivar)=m%bold(indd,ivar,igridd)
           end do
           call hydro_refine(r,uug,uum,uud,bbg,bbm,bbd,ok)
#else
           call hydro_refine(r,uug,uum,uud,ok)
#endif
        end do

        do i_nbor=1,twondim
           call unlock_cache(m,igridn(i_nbor))
        end do

        ! Count only newly flagged cells
        if(m%flag1(ind,igrid)==0.and.ok)g%nflag=g%nflag+1
        if(ok)m%flag1(ind,igrid)=1

     end do
     ! End loop over cells
  end do
  ! End loop over grids

  call close_cache(mdl)

  end associate

#endif

end subroutine hydro_flag
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine pack_fetch_hydro(mesh,igrid,msg_size,msg_array)
  use amr_parameters, only: ndim,twotondim
  use hydro_parameters, only: nvar
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

#ifdef MHD
  do ivar=1,6
     do ind=1,twotondim
        msg%realdp_mhd(ind,ivar)=mesh%bold(ind,ivar,igrid)
     end do
  end do
#endif

  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_hydro
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine unpack_fetch_hydro(mesh,igrid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim, twotondim
  use hydro_parameters, only: nvar
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

#ifdef MHD
  do ivar=1,6
     do ind=1,twotondim
        mesh%bold(ind,ivar,igrid)=msg%realdp_mhd(ind,ivar)
     end do
  end do
#endif

end subroutine unpack_fetch_hydro
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
end module hydro_flag_module
