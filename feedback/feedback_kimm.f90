module feedback_kimm_module
  !==================================================================
  ! Mechanical supernova feedback with one deposit per host cell.
  ! Written by Taysun Kimm and Daniel Han (mini-ramses version in November 2026).
  !==================================================================
  use amr_parameters, only: ndim, twotondim
  implicit none

  ! Slots of a host cell record
  integer,parameter::Q_M=1              ! ejecta mass
  integer,parameter::Q_VX=2,Q_VY=3,Q_VZ=4 ! ejecta velocity
  integer,parameter::Q_Z=5              ! ejecta metal fraction
  integer,parameter::NQ=5

  ! Host cell records of the current call
  integer,save::sn_ncell=0,sn_cap=0
  integer,allocatable,dimension(:),save::sn_igrid,sn_icell
  real(kind=8),allocatable,dimension(:,:),save::sn_q

  ! Neighbor cells with a deposit in the coarse step: uold plus deposits (nb_u),
  ! last energy and entropy (nb_e, nb_s), uold at entry (nb_o), unsent deposits (nb_p)
  integer,parameter::NB_CAP_MAX=2**23
  integer,save::nb_cap=0,nb_n=0,nb_step=-1
  integer(kind=8),allocatable,dimension(:,:),save::nb_key
  real(kind=8),allocatable,dimension(:,:),save::nb_u,nb_o,nb_p
  real(kind=8),allocatable,dimension(:),save::nb_e,nb_s

  ! Run parameters used by the owner-side combine, and a full table seen by it
  real(kind=8),save::kb_smallr=0d0,kb_gamma=0d0
  logical,save::kb_entropy=.false.,kb_full=.false.
  integer,save::kb_ientropy=0

  ! Records of a cache grid, flushed to its owner
  type msg_kimm_record
     real(kind=8),dimension(1:twotondim,1:NQ)::q
  end type msg_kimm_record

  ! Density and momentum deposited into a cache grid, flushed to its owner
  type msg_kimm_mech
     real(kind=8),dimension(1:twotondim,1:4)::du
  end type msg_kimm_mech

  ! Stencil and blast wave parameters
  integer,parameter::nSNnei=48
  integer,parameter::nSNcen=4
  real(kind=8),parameter::A_SN=2.5d5
  real(kind=8),parameter::expN_SN=-2d0/17d0
  real(kind=8),parameter::expE_SN=+16d0/17d0
  real(kind=8),parameter::expZ_SN=-0.14d0
  real(kind=8),parameter::f_LOAD=real(nSNnei,kind=8)/real(nSNcen+nSNnei,kind=8) ! mass loading of the ejecta
  real(kind=8),parameter::f_ESN=0.676d0 ! Blondin et al. (98) at t=trad

contains
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine mechanical_feedback_kimm(s,p,ilevel,texp_old,msn_loc)
  use ramses_commons, only: ramses_t
  use pm_commons, only: part_t
  use nbors_utils
  use cache_commons, only: msg_large_realdp
  use cache, only: open_cache, close_cache
  use marshal, only: pack_fetch_refine, unpack_fetch_refine
  use mdl_module, only: mdl_abort
  implicit none
  type(ramses_t)::s
  type(part_t)::p
  integer::ilevel
  real(kind=8)::texp_old,msn_loc
  !==================================================================
  ! Collects the supernovae of level ilevel per host cell and deposits them.
  ! The routine modifies hydro variables unew.
  !==================================================================
  integer,dimension(1:ndim)::ckey
  integer(kind=8),dimension(0:ndim)::hash_cell
  integer::ipart,idim,igrid,icell
  real(kind=8)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(kind=8)::dx_loc,birth_time,t_SN,mejecta
  real(kind=8),dimension(1:NQ)::q
  logical::ok_level,ok_face
  type(msg_large_realdp)::dummy_large_realdp

#ifdef HYDRO
#if NDIM==3
  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Mesh spacing in that level
  dx_loc=r%boxlen/2**ilevel

  ! Supernovae progenitors life time from Myr to proper time in code units
  t_SN=r%t_SNII*1d6*(365.*24.*3600.)/(scale_t/g%aexp**2)

  ! Total mass of ejecta
  msn_loc=0d0

  ! Open cache for host cells; a cache grid flushes its records to its owner
  sn_ncell=0
  if(.not.allocated(sn_igrid))then
     sn_cap=128
     allocate(sn_igrid(1:sn_cap),sn_icell(1:sn_cap),sn_q(1:NQ,1:sn_cap))
  endif
  call open_cache(mdl, m, pack_size=storage_size(dummy_large_realdp)/32, &
       pack=pack_fetch_refine, unpack=unpack_fetch_refine, &
       init=init_flush_kimm, flush=pack_flush_kimm, combine=unpack_flush_kimm)

  ! Loop over particles in Hilbert order
  do ipart=p%headp(ilevel),p%tailp(ilevel)

     ! Select only stars whose supernova time falls within this level step
     birth_time=p%tp(ipart) ! Proper time
     if(birth_time.lt.(texp_old-t_SN))cycle ! Already exploded
     if(birth_time.ge.(g%texp-t_SN))cycle ! Not old enough

     ok_level=.true.

     ! Find parent cell at level ilevel
     do idim=1,ndim
        ckey(idim)=int((p%xp(ipart,idim)+m%skip(idim))/dx_loc)
     end do
     hash_cell(0)=ilevel+1
     hash_cell(1:ndim)=ckey(1:ndim)
     call get_parent_cell(s,hash_cell,igrid,icell,flush_cache=.true.,fetch_cache=.true.)

     ! If cell does not exist at current level, then find cell at coarser level
     if(igrid==0)then
        do idim=1,ndim
           ckey(idim)=int((p%xp(ipart,idim)+m%skip(idim))/dx_loc/2)
        end do
        hash_cell(0)=ilevel
        hash_cell(1:ndim)=ckey(1:ndim)
        call get_parent_cell(s,hash_cell,igrid,icell,flush_cache=.true.,fetch_cache=.true.)
        if(igrid==0)ok_level=.false.
     end if

     ! No feedback from a host cell touching a non-periodic face
     ok_face=.true.
     do idim=1,ndim
        if(.not.r%periodic(idim))then
           if(ckey(idim)==m%box_ckey_min(idim,hash_cell(0)))ok_face=.false.
           if(ckey(idim)==m%box_ckey_max(idim,hash_cell(0))-1)ok_face=.false.
        endif
     end do
     if(.not.ok_face)cycle

     if(.not.ok_level)then
        write(*,*)'Error: mechanical_feedback_kimm found no host cell at this level or the coarser one'
        call mdl_abort(mdl)
     endif
     if(m%grid(igrid)%refined(icell))then
        write(*,*)'Error: mechanical_feedback_kimm found a refined host cell'
        call mdl_abort(mdl)
     endif

     ! Supernova ejecta into the host cell record
     mejecta=r%eta_SNII*p%mp(ipart)
     q(Q_M)=mejecta
     q(Q_VX)=p%vp(ipart,1)
     q(Q_VY)=p%vp(ipart,2)
     q(Q_VZ)=p%vp(ipart,3)
     q(Q_Z)=r%yield_SNII+(1d0-r%yield_SNII)*p%zp(ipart)
     call add_record(igrid,icell,q)

     ! Update particle mass
     p%mp(ipart)=p%mp(ipart)-mejecta

     ! Update total mass of ejecta
     msn_loc=msn_loc+mejecta

  end do
  ! End loop over particles

  call close_cache(mdl)

  ! Deposit the records of the cells this rank owns
  call deposit_kimm(s)

  end associate
#endif
#endif
end subroutine mechanical_feedback_kimm
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine add_record(igrid,icell,q)
  implicit none
  integer,intent(in)::igrid,icell
  real(kind=8),dimension(1:NQ),intent(in)::q
  !------------------------------------------------------------------
  ! Adds a supernova to the record of its host cell (mass-weighted velocity
  ! and metal fraction).
  !------------------------------------------------------------------
  integer::i,old_cap
  integer,allocatable,dimension(:)::igrid_old,icell_old
  real(kind=8),allocatable,dimension(:,:)::q_old
  real(kind=8)::m1,m2

  if(q(Q_M)<=0d0)return

  do i=1,sn_ncell
     if(sn_igrid(i)==igrid.and.sn_icell(i)==icell)then
        m1=sn_q(Q_M,i); m2=q(Q_M)
        sn_q(Q_VX:Q_VZ,i)=(sn_q(Q_VX:Q_VZ,i)*m1+q(Q_VX:Q_VZ)*m2)/(m1+m2)
        sn_q(Q_Z,i)=(sn_q(Q_Z,i)*m1+q(Q_Z)*m2)/(m1+m2)
        sn_q(Q_M,i)=m1+m2
        return
     endif
  end do

  if(sn_ncell==sn_cap)then
     old_cap=sn_cap
     call move_alloc(sn_igrid,igrid_old)
     call move_alloc(sn_icell,icell_old)
     call move_alloc(sn_q,q_old)
     sn_cap=2*old_cap
     allocate(sn_igrid(1:sn_cap),sn_icell(1:sn_cap),sn_q(1:NQ,1:sn_cap))
     sn_igrid(1:old_cap)=igrid_old
     sn_icell(1:old_cap)=icell_old
     sn_q(:,1:old_cap)=q_old
     deallocate(igrid_old,icell_old,q_old)
  endif
  sn_ncell=sn_ncell+1
  sn_igrid(sn_ncell)=igrid
  sn_icell(sn_ncell)=icell
  sn_q(:,sn_ncell)=q
end subroutine add_record
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine init_flush_kimm(mesh,igrid,hash_key)
  use amr_commons, only: mesh_t
  type(mesh_t)::mesh
  integer::igrid
  integer(kind=8),dimension(0:ndim)::hash_key
  mesh%grid(igrid)%lev=int(hash_key(0),kind=4)
  mesh%grid(igrid)%ckey(1:ndim)=int(hash_key(1:ndim),kind=4)
end subroutine init_flush_kimm
!##############################################################################
subroutine pack_flush_kimm(mesh,igrid,msg_size,msg_array)
  use amr_commons, only: mesh_t
  type(mesh_t)::mesh
  integer::igrid,msg_size
  integer,dimension(1:msg_size),optional::msg_array
  ! Records of a cache grid, removed from the list once sent
  integer::i,nw
  type(msg_kimm_record)::msg
  nw=storage_size(msg)/32
  msg%q=0d0
  do i=1,sn_ncell
     if(sn_igrid(i)==igrid)then
        msg%q(sn_icell(i),1:NQ)=sn_q(1:NQ,i)
        sn_igrid(i)=0
     endif
  end do
  msg_array(1:nw)=transfer(msg,msg_array(1:nw))
  if(msg_size>nw)msg_array(nw+1:msg_size)=0
end subroutine pack_flush_kimm
!##############################################################################
subroutine unpack_flush_kimm(mesh,igrid,msg_size,msg_array,hash_key)
  use amr_commons, only: mesh_t
  type(mesh_t)::mesh
  integer::igrid,msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key
  integer::ind,nw
  type(msg_kimm_record)::msg
  nw=storage_size(msg)/32
  mesh%grid(igrid)%lev=int(hash_key(0),kind=4)
  mesh%grid(igrid)%ckey(1:ndim)=int(hash_key(1:ndim),kind=4)
  msg=transfer(msg_array(1:nw),msg)
  do ind=1,twotondim
     if(msg%q(ind,Q_M)>0d0)call add_record(igrid,ind,msg%q(ind,1:NQ))
  end do
end subroutine unpack_flush_kimm
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine nb_begin(nhost,nstep)
  implicit none
  integer,intent(in)::nhost,nstep
  !------------------------------------------------------------------
  ! Empties the neighbor table at the first call of a coarse step.
  !------------------------------------------------------------------
  integer::cap
  if(allocated(nb_key).and.nstep==nb_step)return
  cap=256
  do while(cap<64*nhost.and.cap<2**20)
     cap=2*cap
  end do
  if(cap>nb_cap.or.4*cap<nb_cap)then
     if(allocated(nb_key))deallocate(nb_key,nb_u,nb_o,nb_p,nb_e,nb_s)
     allocate(nb_key(0:ndim,0:cap-1),nb_u(1:4,0:cap-1),nb_o(1:5,0:cap-1),nb_p(1:4,0:cap-1))
     allocate(nb_e(0:cap-1),nb_s(0:cap-1))
     nb_cap=cap
  endif
  nb_key(0,:)=-1
  nb_n=0
  nb_step=nstep
end subroutine nb_begin
!##############################################################################
integer function nb_probe(key)
  implicit none
  integer(kind=8),dimension(0:ndim),intent(in)::key
  ! Slot of key, or the empty slot for it
  integer(kind=8),parameter,dimension(0:3)::mult=(/2654435761_8,73856093_8,19349663_8,83492791_8/)
  integer(kind=8)::h
  integer::idim,i
  h=0
  do idim=0,ndim
     h=ieor(h,key(idim)*mult(idim))
  end do
  i=int(iand(h,int(nb_cap-1,kind=8)))
  do
     if(nb_key(0,i)<0)exit
     if(all(nb_key(:,i)==key))exit
     i=iand(i+1,nb_cap-1)
  end do
  nb_probe=i
end function nb_probe
!##############################################################################
subroutine nb_lookup(mesh,key,it,found)
  use amr_commons, only: mesh_t
  implicit none
  type(mesh_t)::mesh
  integer(kind=8),dimension(0:ndim),intent(in)::key
  integer,intent(out)::it
  logical,intent(out)::found
  !------------------------------------------------------------------
  ! Entry of a cell, made if missing; it=-1 when the table is full.
  ! A half-full table first drops the entries no deposit can use
  ! again, then doubles.
  !------------------------------------------------------------------
  if(2*(nb_n+1)>nb_cap)then
     call nb_rebuild(mesh,nb_cap,.true.)
     if(2*(nb_n+1)>nb_cap)then
        if(2*nb_cap>NB_CAP_MAX)then
           it=-1
           found=.false.
           return
        endif
        call nb_rebuild(mesh,2*nb_cap,.false.)
     endif
  endif
  it=nb_probe(key)
  found=nb_key(0,it)>=0
  if(.not.found)then
     nb_key(:,it)=key
     nb_n=nb_n+1
  endif
end subroutine nb_lookup
!##############################################################################
subroutine nb_rebuild(mesh,cap,purge)
  use amr_commons, only: mesh_t
  use hash, only: hash_getp
  implicit none
  type(mesh_t)::mesh
  integer,intent(in)::cap
  logical,intent(in)::purge
  !------------------------------------------------------------------
  ! Rehash into cap slots. With purge, drop the entries of this rank's
  ! cells whose uold has changed (their level has finished its step)
  ! and the entries of other ranks' cells with nothing left to send.
  !------------------------------------------------------------------
  integer(kind=8),allocatable,dimension(:,:)::key_old
  real(kind=8),allocatable,dimension(:,:)::u_old,o_old,p_old
  real(kind=8),allocatable,dimension(:)::e_old,s_old
  integer(kind=8),dimension(0:ndim)::key_grid
  integer::i,it,cap_old,igrid,ind,idim,ivar
  logical::keep
  cap_old=nb_cap
  call move_alloc(nb_key,key_old)
  call move_alloc(nb_u,u_old)
  call move_alloc(nb_o,o_old)
  call move_alloc(nb_p,p_old)
  call move_alloc(nb_e,e_old)
  call move_alloc(nb_s,s_old)
  nb_cap=cap
  allocate(nb_key(0:ndim,0:nb_cap-1),nb_u(1:4,0:nb_cap-1),nb_o(1:5,0:nb_cap-1),nb_p(1:4,0:nb_cap-1))
  allocate(nb_e(0:nb_cap-1),nb_s(0:nb_cap-1))
  nb_key(0,:)=-1
  nb_n=0
  do i=0,cap_old-1
     if(key_old(0,i)<0)cycle
     if(purge)then
        key_grid(0)=key_old(0,i)-1
        key_grid(1:ndim)=key_old(1:ndim,i)/2
        ind=1
        do idim=1,ndim
           ind=ind+2**(idim-1)*int(key_old(idim,i)-2*key_grid(idim))
        end do
        igrid=hash_getp(mesh%grid_dict,key_grid)
        if(igrid>0.and.igrid<=mesh%ngridmax)then
           keep=.true.
           do ivar=1,5
              if(o_old(ivar,i)/=mesh%uold(ind,ivar,igrid))keep=.false.
           end do
        else
           keep=any(p_old(:,i)/=0d0)
        endif
        if(.not.keep)cycle
     endif
     it=nb_probe(key_old(:,i))
     nb_key(:,it)=key_old(:,i)
     nb_u(:,it)=u_old(:,i)
     nb_o(:,it)=o_old(:,i)
     nb_p(:,it)=p_old(:,i)
     nb_e(it)=e_old(i)
     nb_s(it)=s_old(i)
     nb_n=nb_n+1
  end do
  deallocate(key_old,u_old,o_old,p_old,e_old,s_old)
end subroutine nb_rebuild
!##############################################################################
integer function nb_find(key)
  implicit none
  integer(kind=8),dimension(0:ndim),intent(in)::key
  ! Slot of key, or -1
  integer::i
  nb_find=-1
  if(.not.allocated(nb_key))return
  i=nb_probe(key)
  if(nb_key(0,i)>=0)nb_find=i
end function nb_find
!##############################################################################
subroutine cell_key(mesh,igrid,ind,key)
  use amr_commons, only: mesh_t
  implicit none
  type(mesh_t)::mesh
  integer,intent(in)::igrid,ind
  integer(kind=8),dimension(0:ndim),intent(out)::key
  ! Hash key of cell ind of grid igrid
  integer::idim
  key(0)=mesh%grid(igrid)%lev+1
  do idim=1,ndim
     key(idim)=2*mesh%grid(igrid)%ckey(idim)+mod((ind-1)/2**(idim-1),2)
  end do
end subroutine cell_key
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine init_flush_kimm_mech(mesh,igrid,hash_key)
  use amr_commons, only: mesh_t
  type(mesh_t)::mesh
  integer::igrid
  integer(kind=8),dimension(0:ndim)::hash_key
  mesh%grid(igrid)%lev=int(hash_key(0),kind=4)
  mesh%grid(igrid)%ckey(1:ndim)=int(hash_key(1:ndim),kind=4)
end subroutine init_flush_kimm_mech
!##############################################################################
subroutine pack_flush_kimm_mech(mesh,igrid,msg_size,msg_array)
  use amr_commons, only: mesh_t
  type(mesh_t)::mesh
  integer::igrid,msg_size
  integer,dimension(1:msg_size),optional::msg_array
  ! Density and momentum deposited into the cells of a cache grid
  integer(kind=8),dimension(0:ndim)::key
  integer::ind,it
  type(msg_kimm_mech)::msg
  msg%du=0d0
  do ind=1,twotondim
     call cell_key(mesh,igrid,ind,key)
     it=nb_find(key)
     if(it>=0)then
        msg%du(ind,1:4)=nb_p(1:4,it)
        nb_p(1:4,it)=0d0
     endif
  end do
  msg_array=transfer(msg,msg_array)
end subroutine pack_flush_kimm_mech
!##############################################################################
subroutine unpack_flush_kimm_mech(mesh,igrid,msg_size,msg_array,hash_key)
  use amr_commons, only: mesh_t
  type(mesh_t)::mesh
  integer::igrid,msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key
  !------------------------------------------------------------------
  ! The owner adds the deposits made on other ranks to the entry of each
  ! cell and sets the cell energy from the entry.
  !------------------------------------------------------------------
  type(msg_kimm_mech)::msg
  integer(kind=8),dimension(0:ndim)::key
  integer::ind,it,ivar
#ifdef MHD
  integer::idim
#endif
  logical::found
  real(kind=8)::d_nei,ekk_nei,emag_nei,eth_nei,d_new,u_new,v_new,w_new,ekk_new,eth_new,e_new,s_new
  mesh%grid(igrid)%lev=int(hash_key(0),kind=4)
  mesh%grid(igrid)%ckey(1:ndim)=int(hash_key(1:ndim),kind=4)
  msg=transfer(msg_array,msg)
  do ind=1,twotondim
     if(all(msg%du(ind,1:4)==0d0))cycle
     call cell_key(mesh,igrid,ind,key)
     call nb_lookup(mesh,key,it,found)
     if(it<0)then
        kb_full=.true.
        cycle
     endif
     if(found)then
        do ivar=1,5
           if(nb_o(ivar,it)/=mesh%uold(ind,ivar,igrid))found=.false.
        end do
     endif
     if(.not.found)then
        do ivar=1,4
           nb_u(ivar,it)=mesh%uold(ind,ivar,igrid)
        end do
        do ivar=1,5
           nb_o(ivar,it)=mesh%uold(ind,ivar,igrid)
        end do
        nb_e(it)=mesh%uold(ind,5,igrid)
        nb_s(it)=0d0
        if(kb_entropy)nb_s(it)=mesh%uold(ind,kb_ientropy,igrid)
        nb_p(:,it)=0d0
     endif
     d_nei=max(dble(mesh%uold(ind,1,igrid)),kb_smallr)
     ekk_nei=0.5d0*(mesh%uold(ind,2,igrid)**2+mesh%uold(ind,3,igrid)**2+mesh%uold(ind,4,igrid)**2)/d_nei
     emag_nei=0d0
#ifdef MHD
     do idim=1,3
        emag_nei=emag_nei+0.125d0*(mesh%bold(ind,idim,igrid)+mesh%bold(ind,idim+3,igrid))**2
     end do
#endif
     eth_nei=mesh%uold(ind,5,igrid)-ekk_nei-emag_nei
     nb_u(1,it)=nb_u(1,it)+msg%du(ind,1)
     nb_u(2:4,it)=nb_u(2:4,it)+msg%du(ind,2:4)
     d_new=max(nb_u(1,it),kb_smallr)
     u_new=nb_u(2,it)/d_new
     v_new=nb_u(3,it)/d_new
     w_new=nb_u(4,it)/d_new
     ekk_new=0.5d0*d_new*(u_new**2+v_new**2+w_new**2)
     eth_new=eth_nei*(d_new/d_nei)
     e_new=ekk_new+eth_new+emag_nei
     mesh%unew(ind,5,igrid)=(mesh%unew(ind,5,igrid)-nb_e(it))+e_new
     nb_e(it)=e_new
     if(kb_entropy)then
        s_new=eth_new/d_new**(kb_gamma-1d0)*(kb_gamma-1d0)
        mesh%unew(ind,kb_ientropy,igrid)=(mesh%unew(ind,kb_ientropy,igrid)-nb_s(it))+s_new
        nb_s(it)=s_new
     endif
  end do
end subroutine unpack_flush_kimm_mech
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine deposit_kimm(s)
  use hydro_parameters, only: nvar
  use cr_parameters, only: ncrgrp
  use ramses_commons, only: ramses_t
  use nbors_utils
  use cache_commons, only: msg_large_realdp
  use cache, only: open_cache, close_cache
  use godunov_fine_module, only: init_flush_godunov, pack_flush_godunov, unpack_flush_godunov
  use marshal, only: pack_fetch_refine, unpack_fetch_refine
  use boundaries, only: init_bound_refine
  use mdl_module, only: mdl_abort
  implicit none
  type(ramses_t)::s
  !==================================================================
  ! Deposits the record of every host cell of this rank into its stencil.
  !==================================================================
  integer::isn,igrid,icell,igridn,icelln,i,j,k,jj,ind,idim,ivar,igrp,it,nsend
  integer(kind=8),allocatable,dimension(:,:)::key_send
  integer,dimension(1:ndim)::ckey,ckey_nbor
  integer(kind=8),dimension(0:ndim)::hash_cell,hash_nbor
  integer(kind=8),dimension(0:ndim,1:nSNnei)::key_nbor
  integer,dimension(1:nSNnei)::igrid_nbor,icell_nbor,level_nbor
  ! Points per host-level cell, first point of the host-level cell (vrep) and
  ! of the actual cell (srep), points per actual cell (nsub)
  integer,dimension(1:nSNnei)::wSNnei_amr,vrep,srep,nsub
  integer,dimension(1:ndim,1:nSNnei)::vkey_nbor
  real(kind=8),dimension(1:3,1:nSNnei)::xSNnei,vSNnei,vSNnei_amr
  real(kind=8),dimension(1:3)::xcen,xnei
  real(kind=8),dimension(1:nvar)::q
  real(kind=8),dimension(1:ncrgrp)::e_SN_crs
  real(kind=8)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(kind=8)::dx_loc,vol_cell,w_nei,w_sub,vol_act,d_act,zsum,fsolid,chi
  real(kind=8)::mejecta,e_SN,m_SN,num_sn,fecrtot,up,vp,wp,zloss
  real(kind=8)::dloss,ekloss,dzloss,dm_ejecta,dm_load
  real(kind=8)::d_nei,Z_nei,nH_nei,Zdepen,f_w_cell,f_w_crit,vload,vload_rad,f_esn2,p_solid
  real(kind=8),dimension(1:nSNnei)::d_nei_pt,f_w_cell_pt,f_w_crit_pt
  real(kind=8)::ekk_nei,eth_nei,emag_nei,d_new,u_new,v_new,w_new,ekk_new,eth_new,emag_new,e_new,s_new
  real(kind=8)::drho_load,chi_sum,w_sum,eth_add,dold,x,y,zz,rr
  real(kind=8),dimension(1:3)::dmom_load,dmom_shell
  logical::ok,is_matched,dump_thermal,chi_undefined,found
  type(msg_large_realdp)::dummy_large_realdp
  type(msg_kimm_mech)::dummy_mech

#ifdef HYDRO
#if NDIM==3
  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

  ! Arrays to define neighbors (center=[0,0,0])
  ind=0
  do k=1,4
     do j=1,4
        do i=1,4
           ok=.true.
           if((i==1.or.i==4).and.(j==1.or.j==4).and.(k==1.or.k==4)) ok=.false. ! edge
           if((i==2.or.i==3).and.(j==2.or.j==3).and.(k==2.or.k==3)) ok=.false. ! centre
           if(ok)then
              ind=ind+1
              x=(i-1)+0.5d0-2d0
              y=(j-1)+0.5d0-2d0
              zz=(k-1)+0.5d0-2d0
              rr=sqrt(x**2+y**2+zz**2)
              xSNnei(1,ind)=x/2d0
              xSNnei(2,ind)=y/2d0
              xSNnei(3,ind)=zz/2d0
              vSNnei(1,ind)=x/rr
              vSNnei(2,ind)=y/rr
              vSNnei(3,ind)=zz/rr
           endif
        end do
     end do
  end do

  ! Conversion factor from user units to cgs units
  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Supernovae specific energy from cgs to code units
  e_SN=r%E_SNII/(r%M_SNII*2d33)/scale_v**2

  ! Supernovae progenitor mass in code units
  m_SN=r%M_SNII*2d33/(scale_d*scale_l**3)

  if(r%cr)then
     ! Reserve part of the SN energy for cosmic rays
     fecrtot=0.0
     e_SN_crs(:)=0.0
     do igrp=1,ncrgrp
        e_SN_crs(igrp)=e_SN*r%fecr(igrp)
        fecrtot=fecrtot+r%fecr(igrp)
     enddo
     e_SN=e_SN*(1.0-fecrtot)
  endif

  ! Neighbor table of the coarse step, and run parameters for the combine
  call nb_begin(sn_ncell,g%nstep_coarse)
  kb_smallr=r%smallr
  kb_gamma=r%gamma
  kb_entropy=r%entropy.and.r%dual_energy>=0
  kb_ientropy=0
  if(kb_entropy)kb_ientropy=r%ientropy
  kb_full=.false.

  ! Open cache for array uold (fetch) and unew (flush)
  call open_cache(mdl, m, pack_size=storage_size(dummy_large_realdp)/32, &
       pack=pack_fetch_refine, unpack=unpack_fetch_refine, &
       init=init_flush_godunov, flush=pack_flush_godunov, &
       combine=unpack_flush_godunov, bound=init_bound_refine)

  ! Loop over host cell records
  do isn=1,sn_ncell
     igrid=sn_igrid(isn)
     icell=sn_icell(isn)

     ! Records of cache grids are deposited by their owner
     if(igrid==0.or.igrid>m%ngridmax)cycle
     if(m%grid(igrid)%refined(icell))then
        write(*,*)'Error: mechanical_feedback_kimm found a refined host cell'
        call mdl_abort(mdl)
     endif

     ! Host cell key, position and volume
     hash_cell(0)=m%grid(igrid)%lev+1
     do idim=1,ndim
        ckey(idim)=2*m%grid(igrid)%ckey(idim)+mod((icell-1)/2**(idim-1),2)
        hash_cell(idim)=ckey(idim)
        xcen(idim)=ckey(idim)+0.5d0
     end do
     dx_loc=r%boxlen/2**m%grid(igrid)%lev
     vol_cell=dx_loc**ndim

     ! Supernova properties of the cell
     mejecta=sn_q(Q_M,isn)
     up=sn_q(Q_VX,isn)
     vp=sn_q(Q_VY,isn)
     wp=sn_q(Q_VZ,isn)
     zloss=sn_q(Q_Z,isn)
     num_sn=mejecta/m_SN
     dloss=mejecta/vol_cell
     ekloss=dloss*0.5d0*(up**2+vp**2+wp**2)
     dzloss=dloss*zloss

     call lock_cache(m,igrid)

     ! Collect all neighboring cells from hash table
     do j=1,nSNnei
        xnei(1:ndim)=xcen(1:ndim)+xSNnei(1:ndim,j)
        ! Periodic boundary conditions
        do idim=1,ndim
           if(r%periodic(idim))then
              if(xnei(idim)< m%box_ckey_min(idim,hash_cell(0)))xnei(idim)=xnei(idim)-m%box_ckey_min(idim,hash_cell(0))+m%box_ckey_max(idim,hash_cell(0))
              if(xnei(idim)>=m%box_ckey_max(idim,hash_cell(0)))xnei(idim)=xnei(idim)+m%box_ckey_min(idim,hash_cell(0))-m%box_ckey_max(idim,hash_cell(0))
           endif
        end do
        ! Get neighboring cell at current level
        ckey_nbor(1:ndim)=int(xnei(1:ndim))
        vkey_nbor(1:ndim,j)=ckey_nbor(1:ndim)
        hash_nbor(0)=hash_cell(0)
        hash_nbor(1:ndim)=ckey_nbor(1:ndim)
        call get_parent_cell(s,hash_nbor,igridn,icelln,flush_cache=.true.,fetch_cache=.true.)
        ! If missing, get neighboring cell at current level - 1
        if(igridn==0)then
           ckey_nbor(1:ndim)=int(xnei(1:ndim)/2.0)
           hash_nbor(0)=hash_cell(0)-1
           hash_nbor(1:ndim)=ckey_nbor(1:ndim)
           call get_parent_cell(s,hash_nbor,igridn,icelln,flush_cache=.true.,fetch_cache=.true.)
        ! If refined, get neighboring cell at current level + 1
        else if(m%grid(igridn)%refined(icelln))then
           ckey_nbor(1:ndim)=int(xnei(1:ndim)*2.0)
           hash_nbor(0)=hash_cell(0)+1
           hash_nbor(1:ndim)=ckey_nbor(1:ndim)
           call get_parent_cell(s,hash_nbor,igridn,icelln,flush_cache=.true.,fetch_cache=.true.)
        endif
        if(igridn==0)then
           write(*,*)'Error: mechanical_feedback_kimm found no stencil cell at any level'
           call mdl_abort(mdl)
        endif
        call lock_cache(m,igridn)
        igrid_nbor(j)=igridn
        icell_nbor(j)=icelln
        level_nbor(j)=hash_nbor(0)-1
        key_nbor(0:ndim,j)=hash_nbor(0:ndim)

        ! Merge the point into the first one in the same cell of the host level
        if(j==1)then
           wSNnei_amr=0
           vSNnei_amr=0d0
        endif
        is_matched=.false.
        do k=1,j-1
           if(all(vkey_nbor(1:ndim,k)==vkey_nbor(1:ndim,j)))then
              wSNnei_amr(k)=wSNnei_amr(k)+1
              vSNnei_amr(1:3,k)=vSNnei_amr(1:3,k)+vSNnei(1:3,j)
              vrep(j)=k
              is_matched=.true.
              exit
           endif
        end do
        if(.not.is_matched)then
           wSNnei_amr(j)=1
           vSNnei_amr(1:3,j)=vSNnei(1:3,j)
           vrep(j)=j
        endif
        ! Group the points of the host-level cell per actual cell
        srep(j)=j
        nsub(j)=1
        do k=vrep(j),j-1
           if(vrep(k)==vrep(j).and.srep(k)==k.and.igrid_nbor(k)==igridn.and. &
                & icell_nbor(k)==icelln.and.level_nbor(k)==level_nbor(j))then
              nsub(k)=nsub(k)+1
              srep(j)=k
              nsub(j)=0
              exit
           endif
        end do
     end do
     do j=1,nSNnei
        if(wSNnei_amr(j)>0)then
           vSNnei_amr(1:3,j)=vSNnei_amr(1:3,j)/sqrt(dot_product(vSNnei_amr(1:3,j),vSNnei_amr(1:3,j)))
        endif
     end do

     ! Update unew in central cell
     m%unew(icell,1,igrid)=m%unew(icell,1,igrid)+dloss-dloss*f_LOAD
     m%unew(icell,2,igrid)=m%unew(icell,2,igrid)+dloss*up-dloss*up*f_LOAD
     m%unew(icell,3,igrid)=m%unew(icell,3,igrid)+dloss*vp-dloss*vp*f_LOAD
     m%unew(icell,4,igrid)=m%unew(icell,4,igrid)+dloss*wp-dloss*wp*f_LOAD
     m%unew(icell,5,igrid)=m%unew(icell,5,igrid)+ekloss-ekloss*f_LOAD
     if(r%metal)m%unew(icell,r%imetal,igrid)=m%unew(icell,r%imetal,igrid)+dzloss-dzloss*f_LOAD

     ! Update passive scalars so that they don't change
     do ivar=6,nvar
        if(r%metal.and.ivar==r%imetal)cycle
        q(ivar)=m%uold(icell,ivar,igrid)/max(dble(m%uold(icell,1,igrid)),r%smallr)
        m%unew(icell,ivar,igrid)=m%unew(icell,ivar,igrid)+(dloss-dloss*f_LOAD)*q(ivar)
     end do

     if(r%cr)then
        ! Inject cosmic rays
#if NENER>0
        do igrp=1,ncrgrp
           m%unew(icell,r%iEcr+igrp-1,igrid)=m%unew(icell,r%iEcr+igrp-1,igrid)+e_SN_crs(igrp)*dloss
           m%unew(icell,5,igrid)=m%unew(icell,5,igrid)+e_SN_crs(igrp)*dloss
        enddo
#endif
     endif

     ! Mass given to each stencil point
     dm_ejecta=dloss*f_LOAD/dble(nSNnei)
     dm_load=dm_ejecta

     ! Thermal dump if the mean f_w_cell/f_w_crit is below mechanical_r_thermal_dump
     chi_sum=0d0
     w_sum=0d0
     chi_undefined=.false.
     do j=1,nSNnei
        if(wSNnei_amr(j)==0)cycle
        w_nei=dble(wSNnei_amr(j))
        igridn=igrid_nbor(j)
        icelln=icell_nbor(j)
        Z_nei=r%z_ave*0.02d0
        if(nsub(j)==wSNnei_amr(j))then
           d_nei=max(dble(m%uold(icelln,1,igridn)),r%smallr)
           if(r%metal)Z_nei=m%uold(icelln,r%imetal,igridn)/d_nei
        else
           d_nei=0d0
           zsum=0d0
           do jj=j,nSNnei
              if(vrep(jj)/=j)cycle
              d_act=max(dble(m%uold(icell_nbor(jj),1,igrid_nbor(jj))),r%smallr)
              d_nei=d_nei+d_act
              if(r%metal)zsum=zsum+m%uold(icell_nbor(jj),r%imetal,igrid_nbor(jj))/d_act
           end do
           d_nei=d_nei/w_nei
           if(r%metal)Z_nei=zsum/w_nei
        endif
        f_w_cell=(dm_load+d_nei/8d0)/dm_ejecta-1d0
        nH_nei=d_nei*scale_nH
        Zdepen=(max(0.01d0,Z_nei/0.02d0))**(expZ_SN*2d0)
        f_w_crit=(A_SN/1d4)**2d0/(f_ESN*r%M_SNII)*num_sn**((expE_SN-1d0)*2d0)*(r%E_SNII/1d51)**(2d0*expE_SN-1d0)*nH_nei**(expN_SN*2d0)*Zdepen-1d0
        f_w_crit=max(0d0,f_w_crit)
        d_nei_pt(j)=d_nei
        f_w_cell_pt(j)=f_w_cell
        f_w_crit_pt(j)=f_w_crit
        if(f_w_crit>0d0)then
           chi_sum=chi_sum+f_w_cell/f_w_crit*w_nei
        else
           chi_undefined=.true.
        endif
        w_sum=w_sum+w_nei
     end do
     dump_thermal=.false.
     if(w_sum>0d0.and..not.chi_undefined)dump_thermal=chi_sum/w_sum<r%mechanical_r_thermal_dump
     if(dump_thermal)then
        ! The central cell share of the blast as heat
        m%unew(icell,5,igrid)=m%unew(icell,5,igrid)+(1d0-f_LOAD)*dloss*e_SN
     endif

     ! Entropy of the central cell from its new energy
     if(r%entropy.and.r%dual_energy>=0)then
        d_new=max(dble(m%unew(icell,1,igrid)),r%smallr)
        u_new=m%unew(icell,2,igrid)/d_new
        v_new=m%unew(icell,3,igrid)/d_new
        w_new=m%unew(icell,4,igrid)/d_new
        ekk_new=0.5d0*d_new*(u_new**2+v_new**2+w_new**2)
        emag_new=0d0
#ifdef MHD
        do idim=1,3
           emag_new=emag_new+0.125d0*(m%bold(icell,idim,igrid)+m%bold(icell,idim+3,igrid))**2
        end do
#endif
        eth_new=m%unew(icell,5,igrid)-ekk_new-emag_new
        m%unew(icell,r%ientropy,igrid)=eth_new/d_new**(r%gamma-1d0)*(r%gamma-1d0)
     endif

     ! Loop over the cells of the host level, each with w_nei of the 48 points
     do j=1,nSNnei

        if(wSNnei_amr(j)==0)cycle
        w_nei=dble(wSNnei_amr(j))

        ! Density, mass ratio and critical mass ratio of the host-level cell
        d_nei=d_nei_pt(j)
        f_w_cell=f_w_cell_pt(j)
        f_w_crit=f_w_crit_pt(j)

        p_solid=0d0
        if(.not.dump_thermal)then
           ! Compute SN terminal momentum
           vload_rad=sqrt(2d0*f_ESN*e_SN*(1d0+f_w_crit))/(1d0+f_w_cell)/f_LOAD
           if(f_w_cell.ge.f_w_crit)then ! radiative phase
              vload=vload_rad
              p_solid=(1d0+f_w_cell)*dm_ejecta*vload*w_nei
           else ! adiabatic phase: the whole cell mass moves
              chi=f_w_cell/f_w_crit
              f_esn2=min(1d0,-2.6495d0*chi**3+6.2807d0*chi**2-5.2005d0*chi+1.7882d0)
              fsolid=w_nei/dble(nSNnei)
              vload=sqrt(2d0*e_SN*f_esn2)/sqrt(d_nei*vol_cell/(mejecta*fsolid)+f_LOAD)
              if(vload>vload_rad)vload=vload_rad
              p_solid=(f_LOAD*mejecta*fsolid+d_nei*vol_cell)*vload/vol_cell
           endif
        endif

        ! Each actual cell takes the share of its points over its volume
        do jj=j,nSNnei
           if(vrep(jj)/=j.or.srep(jj)/=jj)cycle
           igridn=igrid_nbor(jj)
           icelln=icell_nbor(jj)
           w_sub=dble(nsub(jj))
           vol_act=dble(twotondim)**(hash_cell(0)-1-level_nbor(jj))
           d_act=max(dble(m%uold(icelln,1,igridn)),r%smallr)
           ekk_nei=0.5d0*(m%uold(icelln,2,igridn)**2+m%uold(icelln,3,igridn)**2+m%uold(icelln,4,igridn)**2)/d_act
           emag_nei=0d0
#ifdef MHD
           do idim=1,3
              emag_nei=emag_nei+0.125d0*(m%bold(icelln,idim,igridn)+m%bold(icelln,idim+3,igridn))**2
           end do
#endif
           eth_nei=m%uold(icelln,5,igridn)-ekk_nei-emag_nei

           ! Add mass and momentum coming from central cell loading
           drho_load=dloss*f_LOAD*w_sub/dble(nSNnei)/vol_act
           dmom_load(1)=dloss*up*f_LOAD*w_sub/dble(nSNnei)/vol_act
           dmom_load(2)=dloss*vp*f_LOAD*w_sub/dble(nSNnei)/vol_act
           dmom_load(3)=dloss*wp*f_LOAD*w_sub/dble(nSNnei)/vol_act
           m%unew(icelln,1,igridn)=m%unew(icelln,1,igridn)+drho_load
           m%unew(icelln,2,igridn)=m%unew(icelln,2,igridn)+dmom_load(1)
           m%unew(icelln,3,igridn)=m%unew(icelln,3,igridn)+dmom_load(2)
           m%unew(icelln,4,igridn)=m%unew(icelln,4,igridn)+dmom_load(3)
           if(r%metal)then
              m%unew(icelln,r%imetal,igridn)=m%unew(icelln,r%imetal,igridn)+dzloss*f_LOAD*w_sub/dble(nSNnei)/vol_act
           endif
           do ivar=6,nvar
              if(r%metal.and.ivar==r%imetal)cycle
              if(r%entropy.and.r%dual_energy>=0.and..not.dump_thermal.and.ivar==r%ientropy)cycle
              m%unew(icelln,ivar,igridn)=m%unew(icelln,ivar,igridn)+dloss*f_LOAD*q(ivar)*w_sub/dble(nSNnei)/vol_act
           end do

           if(dump_thermal)then
              ! The cell share of the blast as heat
              eth_add=dloss*f_LOAD*w_sub/dble(nSNnei)/vol_act*e_SN
              m%unew(icelln,5,igridn)=m%unew(icelln,5,igridn)+ekloss*f_LOAD*w_sub/dble(nSNnei)/vol_act+eth_add
              if(r%entropy.and.r%dual_energy>=0)then
                 dold=max(dble(m%uold(icelln,1,igridn)),r%smallr)
                 m%unew(icelln,r%ientropy,igridn)=m%unew(icelln,r%ientropy,igridn)+eth_add/dold**(r%gamma-1d0)*(r%gamma-1d0)
              endif
           else
              ! Add momentum coming from the cold shell along the mean direction
              if(nsub(jj)==wSNnei_amr(j))then
                 dmom_shell(1:3)=p_solid*vSNnei_amr(1:3,j)/vol_act
              else
                 dmom_shell(1:3)=(p_solid*(w_sub/w_nei))*vSNnei_amr(1:3,j)/vol_act
              endif
              m%unew(icelln,2,igridn)=m%unew(icelln,2,igridn)+dmom_shell(1)
              m%unew(icelln,3,igridn)=m%unew(icelln,3,igridn)+dmom_shell(2)
              m%unew(icelln,4,igridn)=m%unew(icelln,4,igridn)+dmom_shell(3)

              ! Energy from uold plus all deposits, at the temperature of uold;
              ! a cache copy keeps its deposit for the owner
              call nb_lookup(m,key_nbor(0:ndim,jj),it,found)
              if(it<0)then
                 write(*,*)'Error: mechanical_feedback_kimm neighbor table is full'
                 call mdl_abort(mdl)
              endif
              if(igridn>m%ngridmax)then
                 if(.not.found)nb_p(:,it)=0d0
                 nb_p(1,it)=nb_p(1,it)+drho_load
                 nb_p(2:4,it)=(nb_p(2:4,it)+dmom_load(1:3))+dmom_shell(1:3)
                 cycle
              endif
              if(found)then
                 do ivar=1,5
                    if(nb_o(ivar,it)/=m%uold(icelln,ivar,igridn))found=.false.
                 end do
              endif
              if(.not.found)then
                 nb_p(:,it)=0d0
                 do ivar=1,4
                    nb_u(ivar,it)=m%uold(icelln,ivar,igridn)
                 end do
                 do ivar=1,5
                    nb_o(ivar,it)=m%uold(icelln,ivar,igridn)
                 end do
                 nb_e(it)=m%uold(icelln,5,igridn)
                 nb_s(it)=0d0
                 if(r%entropy.and.r%dual_energy>=0)nb_s(it)=m%uold(icelln,r%ientropy,igridn)
              endif
              nb_u(1,it)=nb_u(1,it)+drho_load
              nb_u(2:4,it)=(nb_u(2:4,it)+dmom_load(1:3))+dmom_shell(1:3)
              d_new=max(nb_u(1,it),r%smallr)
              u_new=nb_u(2,it)/d_new
              v_new=nb_u(3,it)/d_new
              w_new=nb_u(4,it)/d_new
              ekk_new=0.5d0*d_new*(u_new**2+v_new**2+w_new**2)
              eth_new=eth_nei*(d_new/d_act)
              e_new=ekk_new+eth_new+emag_nei
              m%unew(icelln,5,igridn)=(m%unew(icelln,5,igridn)-nb_e(it))+e_new
              nb_e(it)=e_new
              if(r%entropy.and.r%dual_energy>=0)then
                 s_new=eth_new/d_new**(r%gamma-1d0)*(r%gamma-1d0)
                 m%unew(icelln,r%ientropy,igridn)=(m%unew(icelln,r%ientropy,igridn)-nb_s(it))+s_new
                 nb_s(it)=s_new
              endif
           endif
        end do
        ! End loop over the actual cells

     end do
     ! End loop over the cells of the host level

     ! Unlock all octs
     call unlock_cache(m,igrid)
     do j=1,nSNnei
        call unlock_cache(m,igrid_nbor(j))
     end do

  end do
  ! End loop over host cell records

  call close_cache(mdl)
  sn_ncell=0

  ! Send the deposits made into cache copies to their owners
  nsend=0
  do it=0,nb_cap-1
     if(nb_key(0,it)>=0.and.nb_p(1,it)/=0d0)nsend=nsend+1
  end do
  allocate(key_send(0:ndim,1:max(nsend,1)))
  nsend=0
  do it=0,nb_cap-1
     if(nb_key(0,it)>=0.and.nb_p(1,it)/=0d0)then
        nsend=nsend+1
        key_send(0:ndim,nsend)=nb_key(0:ndim,it)
     endif
  end do
  call open_cache(mdl, m, pack_size=storage_size(dummy_mech)/32, &
       init=init_flush_kimm_mech, flush=pack_flush_kimm_mech, combine=unpack_flush_kimm_mech)
  do i=1,nsend
     hash_nbor(0:ndim)=key_send(0:ndim,i)
     call get_parent_cell(s,hash_nbor,igridn,icelln,flush_cache=.true.,fetch_cache=.false.)
  end do
  call close_cache(mdl)
  deallocate(key_send)
  if(kb_full)then
     write(*,*)'Error: mechanical_feedback_kimm neighbor table is full'
     call mdl_abort(mdl)
  endif

  end associate
#endif
#endif
end subroutine deposit_kimm
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
end module feedback_kimm_module
