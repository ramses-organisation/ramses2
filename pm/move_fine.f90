module move_fine_module
  use rho_fine_module, only: cic_weight, cic_index, tsc_weight, tsc_index, pcs_weight, pcs_index
#ifdef _CUDA
  use part_device, only: gpu_kick_drift_part, gpu_kick_drift_star, gpu_kick_drift_sink, gpu_kick_drift_tree
  use part_device, only: gpu_kick_drift_core
#elif defined(_METAL)
  use metal_runner, only: metal_kick_drift_part
#endif
  use rng
  implicit none
  ! Module-level tracer RNG state
  type(RngStream)::tracer_rng
  logical::tracer_rng_ready=.false.
contains
!################################################################
!################################################################
!################################################################
!################################################################
subroutine m_kick_drift_part(pst,ilevel,action_part)
  use amr_parameters, only: ndim, twotondim
  use ramses_commons, only: pst_t
  implicit none
  type(pst_t)::pst
  integer::ilevel
  integer::action_part
  !--------------------------------------------------------------
  ! Move particles according to kick-drift leap frog scheme.
  !--------------------------------------------------------------
  integer,dimension(1:2)::input_array
  integer::dummy(2)

  if(pst%s%m%noct_tot(ilevel)==0)return
  if(pst%s%r%verbose)write(*,'("   Entering kick_drift_part for level",i2," and action=",i2)')ilevel,action_part

  input_array(1)=ilevel
  input_array(2)=action_part
  call r_kick_drift_part(pst,input_array,2,dummy,0)

end subroutine m_kick_drift_part
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_kick_drift_part(pst,input_array,input_size,output_array,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array

  integer::ilevel
  integer::action_part
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_KICK_DRIFT_PART,pst%iUpper+1,input_size,output_size,input_array)
     call r_kick_drift_part(pst%pLower,input_array,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size)
  else
     ilevel=input_array(1)
     action_part=input_array(2)
     ! Force interpolation for various components (part, star, sink, tree...)
     ! based on their respective deposition schemes (CIC 1, TSC 2 or PCS 3)
#ifdef _CUDA
     if (pst%s%r%part) call gpu_kick_drift_part(pst%s, ilevel, action_part)
     if (pst%s%r%star) call gpu_kick_drift_star(pst%s, ilevel, action_part)
     if (pst%s%r%sink) call gpu_kick_drift_sink(pst%s, ilevel, action_part)
     if (pst%s%r%core) call gpu_kick_drift_core(pst%s, ilevel, action_part)
     if (pst%s%r%tree) call gpu_kick_drift_tree(pst%s, ilevel, action_part)
#elif defined(_METAL)
     if (pst%s%r%part) call metal_kick_drift_part(pst%s, ilevel, action_part)
#else
     if(pst%s%r%part)then
        if(pst%s%r%part_force_interpolation_scheme==1)then
           call cic_kick_drift_part(pst%s,pst%s%p   ,ilevel,action_part)
        elseif(pst%s%r%part_force_interpolation_scheme==2)then
           call tsc_kick_drift_part(pst%s,pst%s%p   ,ilevel,action_part)
        elseif(pst%s%r%part_force_interpolation_scheme==3)then
           call pcs_kick_drift_part(pst%s,pst%s%p   ,ilevel,action_part)
        endif
     endif
     if(pst%s%r%star)then
        if(pst%s%r%star_force_interpolation_scheme==1)then
           call cic_kick_drift_part(pst%s,pst%s%star,ilevel,action_part)
        elseif(pst%s%r%star_force_interpolation_scheme==2)then
           call tsc_kick_drift_part(pst%s,pst%s%star,ilevel,action_part)
        elseif(pst%s%r%star_force_interpolation_scheme==3)then
           call pcs_kick_drift_part(pst%s,pst%s%star,ilevel,action_part)
        endif
     endif
     if(pst%s%r%sink)then
        if(pst%s%r%sink_force_interpolation_scheme==1)then
           call cic_kick_drift_part(pst%s,pst%s%sink,ilevel,action_part)
        elseif(pst%s%r%sink_force_interpolation_scheme==2)then
           call tsc_kick_drift_part(pst%s,pst%s%sink,ilevel,action_part)
        elseif(pst%s%r%sink_force_interpolation_scheme==3)then
           call pcs_kick_drift_part(pst%s,pst%s%sink,ilevel,action_part)
        endif
     endif
     if(pst%s%r%core)then
        if(pst%s%r%sink_force_interpolation_scheme==1)then
           call cic_kick_drift_part(pst%s,pst%s%core,ilevel,action_part)
        elseif(pst%s%r%sink_force_interpolation_scheme==2)then
           call tsc_kick_drift_part(pst%s,pst%s%core,ilevel,action_part)
        elseif(pst%s%r%sink_force_interpolation_scheme==3)then
           call pcs_kick_drift_part(pst%s,pst%s%core,ilevel,action_part)
        endif
     endif
     if(pst%s%r%tree)then
        if(pst%s%r%tree_force_interpolation_scheme==1)then
           call cic_kick_drift_part(pst%s,pst%s%tree,ilevel,action_part)
        elseif(pst%s%r%tree_force_interpolation_scheme==2)then
           call tsc_kick_drift_part(pst%s,pst%s%tree,ilevel,action_part)
        elseif(pst%s%r%tree_force_interpolation_scheme==3)then
           call pcs_kick_drift_part(pst%s,pst%s%tree,ilevel,action_part)
        endif
     endif
     if(pst%s%r%trac)then
        if(pst%s%r%trac_interpolation_scheme==0)then
           call mc_trace_gas_part(pst%s,pst%s%trac,ilevel,action_part) ! Classical Monte Carlo (may not work with AMR)
        elseif(pst%s%r%trac_interpolation_scheme==1)then
           call cic_trace_gas_part(pst%s,pst%s%trac,ilevel,action_part)
        elseif(pst%s%r%trac_interpolation_scheme==2)then
           call tsc_trace_gas_part(pst%s,pst%s%trac,ilevel,action_part)
        elseif(pst%s%r%trac_interpolation_scheme==3)then
           call pcs_trace_gas_part(pst%s,pst%s%trac,ilevel,action_part)
        elseif(pst%s%r%trac_interpolation_scheme==4)then
           call cic_trace_gas_part_ito_mc(pst%s,pst%s%trac,ilevel,action_part) ! Ito formulation of the flux-based Monte Carlo tracer
        elseif(pst%s%r%trac_interpolation_scheme==5)then
           call tsc_trace_gas_part_ito_mc(pst%s,pst%s%trac,ilevel,action_part) ! Ito MC tracer with TSC
        elseif(pst%s%r%trac_interpolation_scheme==6)then
           call cic_trace_gas_part_sgs_turb(pst%s,pst%s%trac,ilevel,action_part) ! SGS turbulent diffusion tracer with CIC
        elseif(pst%s%r%trac_interpolation_scheme==7)then
           call tsc_trace_gas_part_sgs_turb(pst%s,pst%s%trac,ilevel,action_part) ! SGS turbulent diffusion tracer with TSC
        endif
     endif
     if(pst%s%r%dust)then
        if(pst%s%r%dust_force_interpolation_scheme==1)then
           call cic_kick_drift_dust(pst%s,pst%s%dust,ilevel,action_part)
        elseif(pst%s%r%dust_force_interpolation_scheme==2)then
           call tsc_kick_drift_dust(pst%s,pst%s%dust,ilevel,action_part)
        elseif(pst%s%r%dust_force_interpolation_scheme==3)then
           call pcs_kick_drift_dust(pst%s,pst%s%dust,ilevel,action_part)
        elseif(pst%s%r%dust_force_interpolation_scheme==4)then
           call cic_kick_drift_dust_ito_mc(pst%s,pst%s%dust,ilevel,action_part) ! Asymptotically approaches Ito MC tracer limit
        elseif(pst%s%r%dust_force_interpolation_scheme==5)then
           call tsc_kick_drift_dust_ito_mc(pst%s,pst%s%dust,ilevel,action_part) ! Asymptotically approaches Ito MC tracer limit
        elseif(pst%s%r%dust_force_interpolation_scheme==6)then
           call cic_kick_drift_dust_guiding_center(pst%s,pst%s%dust,ilevel,action_part) ! CIC guiding center
        elseif(pst%s%r%dust_force_interpolation_scheme==7)then
           call tsc_kick_drift_dust_guiding_center(pst%s,pst%s%dust,ilevel,action_part) ! TSC guiding center
        endif
     endif
#endif
  endif

end subroutine r_kick_drift_part
!################################################################
!################################################################
!################################################################
!################################################################
subroutine cic_kick_drift_part(s,p,ilevel,action_part)
  use amr_parameters, only: ndim, twotondim, nvector
  use pm_parameters
  use pm_commons, only: part_t
  use ramses_commons, only: ramses_t
  use nbors_utils
  use cache_commons
  use cache
  use boundaries, only: init_bound_grav
  implicit none
  type(ramses_t)::s
  type(part_t)::p
  integer::ilevel
  integer::action_part
  !
  !
  real(kind=8),dimension(1:ndim)::x,dr,dl
  integer,dimension(1:ndim)::ir,il
  real(kind=8),dimension(1:twotondim)::vol
  integer,dimension(1:ndim,1:twotondim)::ckey
  integer,dimension(1:twotondim)::icell,igrid
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer::ipart,ind,idim
  real(kind=8)::dx_loc,vol_loc,dteff
  real(kind=8)::gamma,norm2,fnorm,delta
  real(kind=8),dimension(1:ndim)::ff
  logical::ok_level
  real(kind=8),dimension(1:nvector,1:ndim)::xana
  real(kind=8),dimension(1:nvector,1:ndim)::fana
  type(msg_four_realdp)::dummy_four_realdp

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

  if(p%static)then
     ! We still need to set the particle levels correctly, even if they are not moved
     do ipart=p%headp(ilevel),p%tailp(ilevel)
        p%levelp(ipart)=ilevel
     end do
     return
  end if

  ! Mesh spacing in that level
  dx_loc=r%boxlen/2**ilevel
  vol_loc=dx_loc**ndim

  ! Deal with particles that left the computational domain
  if(ilevel==r%levelmin.and.ANY(.not.r%periodic(1:ndim)))then

     ! Loop over particles outside the box
     do ipart=p%headp(ilevel-1),p%tailp(ilevel-1)

        ! Get particle position
        do idim=1,ndim
           x(idim)=p%xp(ipart,idim)
        end do

        ! Without the Poisson solver, particles outside the box have no force
        if(.not.r%poisson)then
           ff(1:ndim)=0d0
        else
           ! Call analytical acceleration routine
           xana(1,1:ndim)=x(1:ndim)
           call gravana(r,g,xana,fana,dx_loc,1)
           ff(1:ndim)=fana(1,1:ndim)
        endif

        ! Perform kick, or drift, or both
        if(action_part==action_kick_drift)then

           ! Update velocity (use levelmin time step)
           p%vp(ipart,1:ndim)=p%vp(ipart,1:ndim)+ff(1:ndim)*0.5d0*g%dtnew(ilevel)

           ! Update position (use levelmin time step)
           p%xp(ipart,1:ndim)=p%xp(ipart,1:ndim)+p%vp(ipart,1:ndim)*g%dtnew(ilevel)

        else if(action_part.EQ.action_kick_only)then

           ! Compute proper time step for second kick
           dteff=g%dtnew(p%levelp(ipart))

           ! Update level to levelmin
           p%levelp(ipart)=ilevel

           ! Update velocity
           p%vp(ipart,1:ndim)=p%vp(ipart,1:ndim)+ff(1:ndim)*0.5d0*dteff

        endif

     end do
     ! End loop over particles

  endif

  ! Open read-only cache
  call open_cache(mdl, m, pack_size=storage_size(dummy_four_realdp)/32, &
       pack=pack_fetch_kick, unpack=unpack_fetch_kick, bound=init_bound_grav)

  ! Loop over particles
  do ipart=p%headp(ilevel),p%tailp(ilevel)

     ! Rescale particle position at level ilevel
     do idim=1,ndim
        x(idim)=(p%xp(ipart,idim)+m%skip(idim))/dx_loc
     end do

     ! CIC at level ilevel (dr: right cloud boundary; dl: left cloud boundary)
     do idim=1,ndim
        dr(idim)=x(idim)+0.5D0
        ir(idim)=int(dr(idim))
        dr(idim)=dr(idim)-ir(idim)
        dl(idim)=1.0D0-dr(idim)
        il(idim)=ir(idim)-1
     end do

     ! Periodic boundary conditions
     do idim=1,ndim
        if(r%periodic(idim))then
           if(il(idim)< m%box_ckey_min(idim,ilevel+1))il(idim)=m%box_ckey_max(idim,ilevel+1)-1
           if(ir(idim)>=m%box_ckey_max(idim,ilevel+1))ir(idim)=m%box_ckey_min(idim,ilevel+1)
        endif
     enddo

     ! Compute cells Cartesian key
     ckey = cic_index(il,ir)

     ! Get parent cell at level ilevel using read-only cache
     ok_level=.true.
     hash_nbor(0)=ilevel+1
     icell=0
     do ind=1,twotondim
        hash_nbor(1:ndim)=ckey(1:ndim,ind)
        call get_parent_cell(s,hash_nbor,igrid(ind),icell(ind),flush_cache=.false.,fetch_cache=.true.,lock=.true.)
        if(igrid(ind)==0)then
           ok_level=.false.
        end if
     end do
     do ind=1,twotondim
        call unlock_cache(m,igrid(ind))
     end do

     ! If cloud is not fully inside level ilevel, re-do CIC at coarser level
     if(.not. ok_level)then

        ! Rescale particle position at level ilevel
        do idim=1,ndim
           x(idim)=x(idim)/2.0d0
        end do

        ! CIC at level ilevel-1 (dr: right cloud boundary; dl: left cloud boundary)
        do idim=1,ndim
           dr(idim)=x(idim)+0.5D0
           ir(idim)=int(dr(idim))
           dr(idim)=dr(idim)-ir(idim)
           dl(idim)=1.0D0-dr(idim)
           il(idim)=ir(idim)-1
        end do

        ! Periodic boundary conditions
        do idim=1,ndim
           if(r%periodic(idim))then
              if(il(idim)< m%box_ckey_min(idim,ilevel))il(idim)=m%box_ckey_max(idim,ilevel)-1
              if(ir(idim)>=m%box_ckey_max(idim,ilevel))ir(idim)=m%box_ckey_min(idim,ilevel)
           endif
        enddo

        ! Compute cells Cartesian key
        ckey = cic_index(il,ir)

        ! Get parent cell at level ilevel-1 using read-only cache
        ok_level=.true.
        hash_nbor(0)=ilevel
        icell=0
        do ind=1,twotondim
           hash_nbor(1:ndim)=ckey(1:ndim,ind)
           call get_parent_cell(s,hash_nbor,igrid(ind),icell(ind),flush_cache=.false.,fetch_cache=.true.,lock=.true.)
           if(igrid(ind)==0)then
              ok_level=.false.
           end if
        end do
        do ind=1,twotondim
           call unlock_cache(m,igrid(ind))
        end do
     end if

     ! Compute cloud volumes
     vol = cic_weight(dl,dr)

     ! Gather 3-force
     ff(1:ndim)=0.0
     if(allocated(p%phip))p%phip(ipart)=0.0
     if(ok_level)then
#ifdef GRAV
        do ind=1,twotondim
           ff(1:ndim)=ff(1:ndim)+m%f(icell(ind),1:ndim,igrid(ind))*vol(ind)
        end do
        if(allocated(p%phip))then ! Store potential
           do ind=1,twotondim
              p%phip(ipart)=p%phip(ipart)+m%phi(icell(ind),igrid(ind))*vol(ind)
           end do
        endif
#endif
     endif

     ! Perform kick, or drift, or both
     if(action_part==action_kick_drift)then

        ! Update velocity
        p%vp(ipart,1:ndim)=p%vp(ipart,1:ndim)+ff(1:ndim)*0.5d0*g%dtnew(ilevel)

        ! Update position
        p%xp(ipart,1:ndim)=p%xp(ipart,1:ndim)+p%vp(ipart,1:ndim)*g%dtnew(ilevel)

        ! For sink particle only
        if(p%type==SINK_TYPE.and.r%sink_descent)then

           ! Compute gradient descent coefficients
           gamma = 0d0
           norm2 = 0d0
           fnorm = 0d0
           do idim=1,ndim
              gamma = gamma + p%vp(ipart,idim)*(ff(idim)-p%fp(ipart,idim))
              norm2 = norm2 + (ff(idim)-p%fp(ipart,idim))**2
              fnorm = fnorm + ff(idim)**2
           enddo
           gamma = gamma*g%dtnew(ilevel) ! cm2 s-2
           fnorm = sqrt(fnorm)
           delta = 0d0
           if(norm2>0)then
              delta = MIN(r%fudge_descent*g%dtnew(ilevel)*sqrt(abs(gamma)/norm2)*fnorm,0.5d0*dx_loc)
           endif
!!$           delta = MIN(r%fudge_descent*g%dtnew(ilevel)*sqrt(abs(gamma)),0.5d0*dx_loc)

           ! Update particle positions
           if(fnorm>0)then
              p%xp(ipart,1:ndim) = p%xp(ipart,1:ndim) + ff(1:ndim)/fnorm*delta
           endif
!!$           if(norm2>0)then
!!$              p%xp(ipart,1:ndim) = p%xp(ipart,1:ndim) + (ff(1:ndim)-p%fp(ipart,1:ndim))/sqrt(norm2)*delta
!!$           endif
        endif

        ! Store force for sink particles
        if(p%type==SINK_TYPE)p%fp(ipart,1:ndim)=ff(1:ndim)
        if(p%type==CORE_TYPE)p%fp(ipart,1:ndim)=ff(1:ndim)

     else if(action_part.EQ.action_kick_only)then

        ! Compute proper time step for second kick
        if (p%levelp(ipart)>=ilevel)then
           dteff=g%dtnew(p%levelp(ipart))
        else
           dteff=g%dtold(p%levelp(ipart))
        endif

        ! Update level
        p%levelp(ipart)=ilevel

        ! Update velocity
        p%vp(ipart,1:ndim)=p%vp(ipart,1:ndim)+ff(1:ndim)*0.5d0*dteff

     endif

  end do
  ! End loop over particles

  call close_cache(mdl)

  ! Periodic boundary conditions
  if(action_part==action_kick_drift)then
     do ipart=p%headp(ilevel),p%tailp(ilevel)
        do idim=1,ndim
           if(r%periodic(idim))then
              if(p%xp(ipart,idim)< 0.0d0           )p%xp(ipart,idim)=p%xp(ipart,idim)+r%box_size(idim)
              if(p%xp(ipart,idim)>=r%box_size(idim))p%xp(ipart,idim)=p%xp(ipart,idim)-r%box_size(idim)
           end if
        end do
     end do
  end if

  end associate

end subroutine cic_kick_drift_part
!################################################################
!################################################################
!################################################################
!################################################################
subroutine tsc_kick_drift_part(s,p,ilevel,action_part)
  use amr_parameters, only: ndim, threetondim
  use pm_parameters
  use pm_commons, only: part_t
  use ramses_commons, only: ramses_t
  use nbors_utils
  use cache_commons
  use cache
  use boundaries, only: init_bound_grav
  implicit none
  type(ramses_t)::s
  type(part_t)::p
  integer::ilevel
  integer::action_part
  !
  !
  real(kind=8),dimension(1:ndim)::x,wl,wc,wr
  integer,dimension(1:ndim)::cl,cc,cr
  real(kind=8),dimension(1:threetondim)::vol
  integer,dimension(1:ndim,1:threetondim)::ckey
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer::ipart,icell,igrid,ind,idim
  real(kind=8)::xl,xc,xr
  real(kind=8)::dx_loc,vol_loc,dteff
  real(kind=8)::gamma,norm2,fnorm,delta
  real(kind=8),dimension(1:ndim)::ff
  logical::ok_level
  type(msg_four_realdp)::dummy_four_realdp

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

  if(p%static)then
     ! We still need to set the particle levels correctly, even if they are not moved
     do ipart=p%headp(ilevel),p%tailp(ilevel)
        p%levelp(ipart)=ilevel
     end do
     return
  end if

  ! Mesh spacing in that level
  dx_loc=r%boxlen/2**ilevel
  vol_loc=dx_loc**ndim

  ! Open read-only cache
  call open_cache(mdl, m, pack_size=storage_size(dummy_four_realdp)/32, &
       pack=pack_fetch_kick, unpack=unpack_fetch_kick, bound=init_bound_grav)

  ! Loop over particles
  do ipart=p%headp(ilevel),p%tailp(ilevel)

     ! Rescale particle position at level ilevel
     do idim=1,ndim
        x(idim)=(p%xp(ipart,idim)+m%skip(idim))/dx_loc
     end do

     ! TSC at level ilevel; a particle contributes to 3 cells in each direction
     do idim=1,ndim
        cl(idim)=int(x(idim))-1 ! cell index
        cc(idim)=int(x(idim))
        cr(idim)=int(x(idim))+1
        xl=dble(cl(idim))+0.5D0 ! cell coordinate
        xc=dble(cc(idim))+0.5D0
        xr=dble(cr(idim))+0.5D0
        wl(idim)=0.5D0*(1.5D0-abs(x(idim)-xl))**2 ! weight
        wc(idim)=0.75D0-         (x(idim)-xc) **2
        wr(idim)=0.5D0*(1.5D0-abs(x(idim)-xr))**2
     end do

     ! Periodic boundary conditions
     do idim=1,ndim
        if(r%periodic(idim))then
           if(cl(idim)< m%box_ckey_min(idim,ilevel+1))cl(idim)=m%box_ckey_max(idim,ilevel+1)-1
           if(cr(idim)>=m%box_ckey_max(idim,ilevel+1))cr(idim)=m%box_ckey_min(idim,ilevel+1)
        endif
     enddo

     ! Compute cells Cartesian key
     ckey = tsc_index(cl,cc,cr)

     ! Compute cloud volumes
     vol = tsc_weight(wl,wc,wr)

     ! Gather 3-force
     hash_nbor(0)=ilevel+1
     ff(1:ndim)=0.0
     if(allocated(p%phip))p%phip(ipart)=0.0
     do ind=1,threetondim
        hash_nbor(1:ndim)=ckey(1:ndim,ind)
        ! Get parent cell at level ilevel using read-only cache
        call get_parent_cell(s,hash_nbor,igrid,icell,flush_cache=.false.,fetch_cache=.true.)
#ifdef GRAV
        if(igrid>0)then
           ff(1:ndim)=ff(1:ndim)+m%f(icell,1:ndim,igrid)*vol(ind)
           if(allocated(p%phip))then ! Store potential
              p%phip(ipart)=p%phip(ipart)+m%phi(icell,igrid)*vol(ind)
           endif
        end if
#endif
     end do

     ! Perform kick, or drift, or both
     if(action_part==action_kick_drift)then

        ! Update velocity
        p%vp(ipart,1:ndim)=p%vp(ipart,1:ndim)+ff(1:ndim)*0.5d0*g%dtnew(ilevel)

        ! Update position
        p%xp(ipart,1:ndim)=p%xp(ipart,1:ndim)+p%vp(ipart,1:ndim)*g%dtnew(ilevel)

        ! For sink particle only
        if(p%type==SINK_TYPE.and.r%sink_descent)then

           ! Compute gradient descent coefficients
           gamma = 0d0
           norm2 = 0d0
           fnorm = 0d0
           do idim=1,ndim
              gamma = gamma + p%vp(ipart,idim)*(ff(idim)-p%fp(ipart,idim))
              norm2 = norm2 + (ff(idim)-p%fp(ipart,idim))**2
              fnorm = fnorm + ff(idim)**2
           enddo
           gamma = gamma*g%dtnew(ilevel) ! cm2 s-2
           fnorm = sqrt(fnorm)
           delta = 0d0
           if(norm2>0)then
              delta = MIN(r%fudge_descent*g%dtnew(ilevel)*sqrt(abs(gamma)/norm2)*fnorm,0.5d0*dx_loc)
           endif
!!$           delta = MIN(r%fudge_descent*g%dtnew(ilevel)*sqrt(abs(gamma)),0.5d0*dx_loc)

           ! Update particle positions
           if(fnorm>0)then
              p%xp(ipart,1:ndim) = p%xp(ipart,1:ndim) + ff(1:ndim)/fnorm*delta
           endif
!!$           if(norm2>0)then
!!$              p%xp(ipart,1:ndim) = p%xp(ipart,1:ndim) + (ff(1:ndim)-p%fp(ipart,1:ndim))/sqrt(norm2)*delta
!!$           endif
        endif

        ! Store force for sink particles
        if(p%type==SINK_TYPE)p%fp(ipart,1:ndim)=ff(1:ndim)
        if(p%type==CORE_TYPE)p%fp(ipart,1:ndim)=ff(1:ndim)

     else if(action_part.EQ.action_kick_only)then

        ! Compute proper time step for second kick
        if (p%levelp(ipart)>=ilevel)then
           dteff=g%dtnew(p%levelp(ipart))
        else
           dteff=g%dtold(p%levelp(ipart))
        endif

        ! Update level
        p%levelp(ipart)=ilevel

        ! Update velocity
        p%vp(ipart,1:ndim)=p%vp(ipart,1:ndim)+ff(1:ndim)*0.5d0*dteff

     endif

  end do
  ! End loop over particles

  call close_cache(mdl)

  ! Periodic boundary conditions
  if(action_part==action_kick_drift)then
     do ipart=p%headp(ilevel),p%tailp(ilevel)
        do idim=1,ndim
           if(r%periodic(idim))then
              if(p%xp(ipart,idim)< 0.0d0           )p%xp(ipart,idim)=p%xp(ipart,idim)+r%box_size(idim)
              if(p%xp(ipart,idim)>=r%box_size(idim))p%xp(ipart,idim)=p%xp(ipart,idim)-r%box_size(idim)
           end if
        end do
     end do
  end if

  end associate

end subroutine tsc_kick_drift_part
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pcs_kick_drift_part(s,p,ilevel,action_part)
  use amr_parameters, only: ndim, fourtondim
  use pm_parameters
  use pm_commons, only: part_t
  use ramses_commons, only: ramses_t
  use nbors_utils
  use cache_commons
  use cache
  use boundaries, only: init_bound_grav
  implicit none
  type(ramses_t)::s
  type(part_t)::p
  integer::ilevel
  integer::action_part
  !
  !
  real(kind=8),dimension(1:ndim)::x,wll,wl,wr,wrr
  integer,dimension(1:ndim)::cll,cl,cr,crr
  real(kind=8),dimension(1:fourtondim)::vol
  integer,dimension(1:ndim,1:fourtondim)::ckey
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer::ipart,icell,igrid,ind,idim
  real(kind=8)::xll,xl,xr,xrr
  real(kind=8)::dx_loc,vol_loc,dteff
  real(kind=8)::gamma,norm2,fnorm,delta
  real(kind=8),dimension(1:ndim)::ff
  logical::ok_level
  type(msg_four_realdp)::dummy_four_realdp

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

  if(p%static)then
     ! We still need to set the particle levels correctly, even if they are not moved
     do ipart=p%headp(ilevel),p%tailp(ilevel)
        p%levelp(ipart)=ilevel
     end do
     return
  end if

  ! Mesh spacing in that level
  dx_loc=r%boxlen/2**ilevel
  vol_loc=dx_loc**ndim

  ! Open read-only cache
  call open_cache(mdl, m, pack_size=storage_size(dummy_four_realdp)/32, &
       pack=pack_fetch_kick, unpack=unpack_fetch_kick, bound=init_bound_grav)

  ! Loop over particles
  do ipart=p%headp(ilevel),p%tailp(ilevel)

     ! Rescale particle position at level ilevel
     do idim=1,ndim
        x(idim)=(p%xp(ipart,idim)+m%skip(idim))/dx_loc
     end do

     ! PCS at level ilevel; a particle contributes to 4 cells in each direction
     do idim=1,ndim
        crr(idim)=int(x(idim)+1.5D0) ! rightermost cell index
        cr (idim)=crr(idim)-1
        cl (idim)=crr(idim)-2
        cll(idim)=crr(idim)-3
        xll=dble(cll(idim))+0.5D0 ! cell coordinate
        xl =dble(cl (idim))+0.5D0
        xr =dble(cr (idim))+0.5D0
        xrr=dble(crr(idim))+0.5D0
        wll(idim)=(2D0                        -abs(x(idim)-xll))**3/6D0 ! weight
        wl (idim)=(4D0-6D0*(x(idim)-xl)**2+3d0*abs(x(idim)-xl )**3)/6D0
        wr (idim)=(4D0-6D0*(x(idim)-xr)**2+3d0*abs(x(idim)-xr )**3)/6D0
        wrr(idim)=(2D0                        -abs(x(idim)-xrr))**3/6D0
     end do

     ! Periodic boundary conditions
     do idim=1,ndim
        if(r%periodic(idim))then
           if(cll(idim)< m%box_ckey_min(idim,ilevel+1))cll(idim)=cll(idim)-m%box_ckey_min(idim,ilevel+1)+m%box_ckey_max(idim,ilevel+1)
           if(cl (idim)< m%box_ckey_min(idim,ilevel+1))cl (idim)=cl (idim)-m%box_ckey_min(idim,ilevel+1)+m%box_ckey_max(idim,ilevel+1)
           if(cr (idim)>=m%box_ckey_max(idim,ilevel+1))cr (idim)=cr (idim)+m%box_ckey_min(idim,ilevel+1)-m%box_ckey_max(idim,ilevel+1)
           if(crr(idim)>=m%box_ckey_max(idim,ilevel+1))crr(idim)=crr(idim)+m%box_ckey_min(idim,ilevel+1)-m%box_ckey_max(idim,ilevel+1)
        endif
     enddo

     ! Compute cells Cartesian key
     ckey = pcs_index(cll,cl,cr,crr)

     ! Compute cloud volumes
     vol = pcs_weight(wll,wl,wr,wrr)

     ! Gather 3-force
     hash_nbor(0)=ilevel+1
     ff(1:ndim)=0.0
     if(allocated(p%phip))p%phip(ipart)=0.0
     do ind=1,fourtondim
        hash_nbor(1:ndim)=ckey(1:ndim,ind)
        ! Get parent cell at level ilevel using read-only cache
        call get_parent_cell(s,hash_nbor,igrid,icell,flush_cache=.false.,fetch_cache=.true.)
#ifdef GRAV
        if(igrid>0)then
           ff(1:ndim)=ff(1:ndim)+m%f(icell,1:ndim,igrid)*vol(ind)
           if(allocated(p%phip))then ! Store potential
              p%phip(ipart)=p%phip(ipart)+m%phi(icell,igrid)*vol(ind)
           endif
        end if
#endif
     end do

     ! Perform kick, or drift, or both
     if(action_part==action_kick_drift)then

        ! Update velocity
        p%vp(ipart,1:ndim)=p%vp(ipart,1:ndim)+ff(1:ndim)*0.5d0*g%dtnew(ilevel)

        ! Update position
        p%xp(ipart,1:ndim)=p%xp(ipart,1:ndim)+p%vp(ipart,1:ndim)*g%dtnew(ilevel)

        ! For sink particle only
        if(p%type==SINK_TYPE.and.r%sink_descent)then

           ! Compute gradient descent coefficients
           gamma = 0d0
           norm2 = 0d0
           fnorm = 0d0
           do idim=1,ndim
              gamma = gamma + p%vp(ipart,idim)*(ff(idim)-p%fp(ipart,idim))
              norm2 = norm2 + (ff(idim)-p%fp(ipart,idim))**2
              fnorm = fnorm + ff(idim)**2
           enddo
           gamma = gamma*g%dtnew(ilevel) ! cm2 s-2
           fnorm = sqrt(fnorm)
           delta = 0d0
           if(norm2>0)then
              delta = MIN(r%fudge_descent*g%dtnew(ilevel)*sqrt(abs(gamma)/norm2)*fnorm,0.5d0*dx_loc)
           endif
!!$           delta = MIN(r%fudge_descent*g%dtnew(ilevel)*sqrt(abs(gamma)),0.5d0*dx_loc)

           ! Update particle positions
           if(fnorm>0)then
              p%xp(ipart,1:ndim) = p%xp(ipart,1:ndim) + ff(1:ndim)/fnorm*delta
           endif
!!$           if(norm2>0)then
!!$              p%xp(ipart,1:ndim) = p%xp(ipart,1:ndim) + (ff(1:ndim)-p%fp(ipart,1:ndim))/sqrt(norm2)*delta
!!$           endif
        endif

        ! Store force for sink particles
        if(p%type==SINK_TYPE)p%fp(ipart,1:ndim)=ff(1:ndim)
        if(p%type==CORE_TYPE)p%fp(ipart,1:ndim)=ff(1:ndim)

     else if(action_part.EQ.action_kick_only)then

        ! Compute proper time step for second kick
        if (p%levelp(ipart)>=ilevel)then
           dteff=g%dtnew(p%levelp(ipart))
        else
           dteff=g%dtold(p%levelp(ipart))
        endif

        ! Update level
        p%levelp(ipart)=ilevel

        ! Update velocity
        p%vp(ipart,1:ndim)=p%vp(ipart,1:ndim)+ff(1:ndim)*0.5d0*dteff

     endif

  end do
  ! End loop over particles

  call close_cache(mdl)

  ! Periodic boundary conditions
  if(action_part==action_kick_drift)then
     do ipart=p%headp(ilevel),p%tailp(ilevel)
        do idim=1,ndim
           if(r%periodic(idim))then
              if(p%xp(ipart,idim)< 0.0d0           )p%xp(ipart,idim)=p%xp(ipart,idim)+r%box_size(idim)
              if(p%xp(ipart,idim)>=r%box_size(idim))p%xp(ipart,idim)=p%xp(ipart,idim)-r%box_size(idim)
           end if
        end do
     end do
  end if

  end associate

end subroutine pcs_kick_drift_part
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine pack_fetch_kick(mesh,igrid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_four_realdp
  type(mesh_t)::mesh
  integer::igrid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind
  type(msg_four_realdp)::msg

#ifdef GRAV
  do ind=1,twotondim
     msg%force(ind,1:3)=mesh%f(ind,1:3,igrid)
#ifdef OUTPUT_PARTICLE_POTENTIAL
     msg%phi(ind)=mesh%phi(ind,igrid)
#endif
  end do
#endif

  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_kick
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine unpack_fetch_kick(mesh,igrid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim, twotondim
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_four_realdp
  type(mesh_t)::mesh
  integer::igrid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind
  type(msg_four_realdp)::msg

  mesh%grid(igrid)%lev=hash_key(0)
  mesh%grid(igrid)%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)

#ifdef GRAV
  do ind=1,twotondim
     mesh%f(ind,1:3,igrid)=msg%force(ind,1:3)
#ifdef OUTPUT_PARTICLE_POTENTIAL
     mesh%phi(ind,igrid)=msg%phi(ind)
#endif
  end do
#endif

end subroutine unpack_fetch_kick
!#########################################################################
!#########################################################################
! Tracer-only cache pack/unpack (hydro only) to keep default paths intact
! In general, different particle types may need different cache packs/unpacks
!#########################################################################
!#########################################################################
subroutine pack_fetch_kick_trac(mesh,igrid,msg_size,msg_array)
  use amr_parameters, only: twotondim, ndim
  use hydro_parameters, only: nvar
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_hydro_mflux
  type(mesh_t)::mesh
  integer::igrid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind, ivar
  type(msg_hydro_mflux)::msg

#ifdef HYDRO
  do ivar=1,nvar
     do ind=1,twotondim
        msg%realdp_hydro(ind,ivar)=mesh%uold(ind,ivar,igrid)
     end do
  end do
#ifdef TRCFLX
  do ind=1,twotondim
     msg%realdp_mflux(ind,1:2*ndim+1)=mesh%mflux(ind,1:2*ndim+1,igrid)
  end do
#endif
#endif

  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_kick_trac
!#########################################################################
!#########################################################################
subroutine unpack_fetch_kick_trac(mesh,igrid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim, twotondim
  use hydro_parameters, only: nvar
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_hydro_mflux
  type(mesh_t)::mesh
  integer::igrid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key
  integer::ind,ivar
  type(msg_hydro_mflux)::msg

  mesh%grid(igrid)%lev=hash_key(0)
  mesh%grid(igrid)%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)

#ifdef HYDRO
  do ivar=1,nvar
     do ind=1,twotondim
        mesh%uold(ind,ivar,igrid)=msg%realdp_hydro(ind,ivar)
     end do
  end do
#ifdef TRCFLX
  do ind=1,twotondim
     mesh%mflux(ind,1:2*ndim+1,igrid)=msg%realdp_mflux(ind,1:2*ndim+1)
  end do
#endif
#endif
end subroutine unpack_fetch_kick_trac
!#########################################################################
!#########################################################################
subroutine pack_fetch_kick_dust(mesh,igrid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use hydro_parameters, only: nvar
  use rt_parameters, only: nrtvar
  use cr_parameters, only: ncruvar
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_large_realdp
  type(mesh_t)::mesh
  integer::igrid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind, ivar
  type(msg_large_realdp)::msg

#ifdef HYDRO
  do ind=1,twotondim
     do ivar=1,nvar
        msg%realdp_hydro(ind,ivar)=mesh%uold(ind,ivar,igrid)
     end do
  end do
#endif
#ifdef GRAV
  do ind=1,twotondim
     msg%realdp_poisson(ind,1)=mesh%f(ind,1,igrid)
     msg%realdp_poisson(ind,2)=mesh%f(ind,2,igrid)
     msg%realdp_poisson(ind,3)=mesh%f(ind,3,igrid)
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
#ifdef MHD
  do ivar=1,6
     do ind=1,twotondim
        msg%realdp_mhd(ind,ivar)=mesh%bold(ind,ivar,igrid)
     end do
  end do
#endif
  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_kick_dust
!#########################################################################
!#########################################################################
subroutine unpack_fetch_kick_dust(mesh,igrid,msg_size,msg_array,hash_key)
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
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind,ivar
  type(msg_large_realdp)::msg

  mesh%grid(igrid)%lev=hash_key(0)
  mesh%grid(igrid)%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)

#ifdef HYDRO
  do ind=1,twotondim
     do ivar=1,nvar
        mesh%uold(ind,ivar,igrid)=msg%realdp_hydro(ind,ivar)
     end do
  end do
#endif
#ifdef GRAV
  do ind=1,twotondim
     mesh%f(ind,1,igrid)=msg%realdp_poisson(ind,1)
     mesh%f(ind,2,igrid)=msg%realdp_poisson(ind,2)
     mesh%f(ind,3,igrid)=msg%realdp_poisson(ind,3)
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
#ifdef MHD
  do ivar=1,6
     do ind=1,twotondim
        mesh%bold(ind,ivar,igrid)=msg%realdp_mhd(ind,ivar)
     end do
  end do
#endif

end subroutine unpack_fetch_kick_dust
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cic_trace_gas_part(s,p,ilevel,action_part)
  use amr_parameters, only: ndim, twotondim, dp
  use pm_parameters
  use pm_commons, only: part_t
  use ramses_commons, only: ramses_t
  use nbors_utils
  use cache_commons
  use cache
  implicit none
  type(ramses_t)::s
  type(part_t)::p
  integer::ilevel
  integer::action_part
  real(kind=8),dimension(1:ndim)::x,x_mid,dr,dl,dr2,dl2
  integer,dimension(1:ndim)::ir,il
  integer,dimension(1:ndim)::ir2,il2
  real(kind=8),dimension(1:twotondim)::vol,vol2
  integer,dimension(1:ndim,1:twotondim)::ckey,ckey2
  integer::icell,icell2,igrid
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer::ipart,ind,idim
  real(kind=8)::dx_loc,vol_loc
  real(kind=8),dimension(1:ndim)::ff,v_pred
  logical :: ok_level
  type(msg_three_realdp)::dummy_three_realdp
  type(msg_hydro_mflux)::dummy_hydro_mflux

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

  if(p%static)return
  if (p%type/=TRAC_TYPE) return

  dx_loc=r%boxlen/2**ilevel
  vol_loc=dx_loc**ndim

  ! Tracer hydro cache (uold, mflux)
  call open_cache(mdl, m, pack_size=storage_size(dummy_hydro_mflux)/32, &
       pack=pack_fetch_kick_trac, unpack=unpack_fetch_kick_trac)
  do ipart=p%headp(ilevel),p%tailp(ilevel)
     ! Position in cell units at current level
     do idim=1,ndim
        x(idim)=(p%xp(ipart,idim)+m%skip(idim))/dx_loc
     end do
     do idim=1,ndim
        if(r%periodic(idim))then
           if(x(idim)< dble(m%box_ckey_min(idim,ilevel+1)))x(idim)=x(idim)+dble(m%box_ckey_max(idim,ilevel+1)-m%box_ckey_min(idim,ilevel+1))
           if(x(idim)>=dble(m%box_ckey_max(idim,ilevel+1)))x(idim)=x(idim)-dble(m%box_ckey_max(idim,ilevel+1)-m%box_ckey_min(idim,ilevel+1))
        endif
     end do
     ! Gather velocity v = mom/rho at x^n using CIC
     do idim=1,ndim
        dr(idim)=x(idim)+0.5D0
        ir(idim)=int(dr(idim))
        dr(idim)=dr(idim)-ir(idim)
        dl(idim)=1.0D0-dr(idim)
        il(idim)=ir(idim)-1
     end do
     do idim=1,ndim
        if(r%periodic(idim))then
           if(il(idim)< m%box_ckey_min(idim,ilevel+1))il(idim)=m%box_ckey_max(idim,ilevel+1)-1
           if(ir(idim)>=m%box_ckey_max(idim,ilevel+1))ir(idim)=m%box_ckey_min(idim,ilevel+1)
        endif
     enddo
     ckey = cic_index(il,ir)
     vol = cic_weight(dl,dr)
     ff(1:ndim)=0.0
     ok_level=.true.
     hash_nbor(0)=ilevel+1
     do ind=1,twotondim
        hash_nbor(1:ndim)=ckey(1:ndim,ind)
        call get_parent_cell(s,hash_nbor,igrid,icell,flush_cache=.false.,fetch_cache=.true.)
#ifdef HYDRO
        if(igrid>0)then
           ff(1:ndim)=ff(1:ndim)+real(m%uold(icell,2:ndim+1,igrid),kind=8)/max(real(m%uold(icell,1,igrid),kind=8),real(r%smallr,kind=8))*vol(ind)
        else
           ok_level=.false.
        end if
#endif
     end do
     if(.not.ok_level)then
        do idim=1,ndim
           x(idim)=x(idim)/2.0d0
        end do
        do idim=1,ndim
           dr(idim)=x(idim)+0.5D0
           ir(idim)=int(dr(idim))
           dr(idim)=dr(idim)-ir(idim)
           dl(idim)=1.0D0-dr(idim)
           il(idim)=ir(idim)-1
        end do
        do idim=1,ndim
           if(r%periodic(idim))then
              if(il(idim)< m%box_ckey_min(idim,ilevel))il(idim)=m%box_ckey_max(idim,ilevel)-1
              if(ir(idim)>=m%box_ckey_max(idim,ilevel))ir(idim)=m%box_ckey_min(idim,ilevel)
           endif
        enddo
        ckey = cic_index(il,ir)
        vol = cic_weight(dl,dr)
        ff(1:ndim)=0.0
        hash_nbor(0)=ilevel
        do ind=1,twotondim
           hash_nbor(1:ndim)=ckey(1:ndim,ind)
           call get_parent_cell(s,hash_nbor,igrid,icell,flush_cache=.false.,fetch_cache=.true.)
#ifdef HYDRO
           if(igrid>0)then
              ff(1:ndim)=ff(1:ndim)+real(m%uold(icell,2:ndim+1,igrid),kind=8)/max(real(m%uold(icell,1,igrid),kind=8),real(r%smallr,kind=8))*vol(ind)
           end if
#endif
        end do
     end if

     if(action_part==action_kick_only)then
        ! RK2 step 1 (early call): stash v^n at x^n, no move
        p%vp(ipart,1:ndim)=ff(1:ndim)
        p%levelp(ipart)=ilevel

     else if(action_part==action_kick_drift)then
        ! RK2 step 2 (late call): predict midpoint and correct
        ! Use stored v^n if available; fallback to current ff at step 0
        if (g%nstep>0) then
           v_pred(1:ndim)=p%vp(ipart,1:ndim)
        else
           v_pred(1:ndim)=ff(1:ndim)
        endif
        ! Predict x_mid in cell units
        do idim=1,ndim
           x_mid(idim)=x(idim)+0.5d0*g%dtnew(ilevel)*v_pred(idim)/dx_loc
        end do
        do idim=1,ndim
           if(r%periodic(idim))then
              if(x_mid(idim)< dble(m%box_ckey_min(idim,ilevel+1)))x_mid(idim)=x_mid(idim)+dble(m%box_ckey_max(idim,ilevel+1)-m%box_ckey_min(idim,ilevel+1))
              if(x_mid(idim)>=dble(m%box_ckey_max(idim,ilevel+1)))x_mid(idim)=x_mid(idim)-dble(m%box_ckey_max(idim,ilevel+1)-m%box_ckey_min(idim,ilevel+1))
           endif
        end do
        ! Gather v^{n+1} at x_mid using CIC
        do idim=1,ndim
           dr2(idim)=x_mid(idim)+0.5D0
           ir2(idim)=int(dr2(idim))
           dr2(idim)=dr2(idim)-ir2(idim)
           dl2(idim)=1.0D0-dr2(idim)
           il2(idim)=ir2(idim)-1
        end do
        do idim=1,ndim
           if(r%periodic(idim))then
              if(il2(idim)< m%box_ckey_min(idim,ilevel+1))il2(idim)=m%box_ckey_max(idim,ilevel+1)-1
              if(ir2(idim)>=m%box_ckey_max(idim,ilevel+1))ir2(idim)=m%box_ckey_min(idim,ilevel+1)
           endif
        enddo
        ckey2 = cic_index(il2,ir2)
        vol2 = cic_weight(dl2,dr2)
        ff(1:ndim)=0.0
        ok_level=.true.
        hash_nbor(0)=ilevel+1
        do ind=1,twotondim
           hash_nbor(1:ndim)=ckey2(1:ndim,ind)
           call get_parent_cell(s,hash_nbor,igrid,icell2,flush_cache=.false.,fetch_cache=.true.)
#ifdef HYDRO
           if(igrid>0)then
              ff(1:ndim)=ff(1:ndim)+real(m%uold(icell2,2:ndim+1,igrid),kind=8)/max(real(m%uold(icell2,1,igrid),kind=8),real(r%smallr,kind=8))*vol2(ind)
           else
              ok_level=.false.
           end if
#endif
        end do
        if(.not.ok_level)then
           do idim=1,ndim
              x_mid(idim)=x_mid(idim)/2.0d0
           end do
           do idim=1,ndim
              dr2(idim)=x_mid(idim)+0.5D0
              ir2(idim)=int(dr2(idim))
              dr2(idim)=dr2(idim)-ir2(idim)
              dl2(idim)=1.0D0-dr2(idim)
              il2(idim)=ir2(idim)-1
           end do
           do idim=1,ndim
              if(r%periodic(idim))then
                 if(il2(idim)< m%box_ckey_min(idim,ilevel))il2(idim)=m%box_ckey_max(idim,ilevel)-1
                 if(ir2(idim)>=m%box_ckey_max(idim,ilevel))ir2(idim)=m%box_ckey_min(idim,ilevel)
              endif
           enddo
           ckey2 = cic_index(il2,ir2)
           vol2 = cic_weight(dl2,dr2)
           ff(1:ndim)=0.0
           hash_nbor(0)=ilevel
           do ind=1,twotondim
              hash_nbor(1:ndim)=ckey2(1:ndim,ind)
              call get_parent_cell(s,hash_nbor,igrid,icell2,flush_cache=.false.,fetch_cache=.true.)
#ifdef HYDRO
              if(igrid>0)then
                 ff(1:ndim)=ff(1:ndim)+real(m%uold(icell2,2:ndim+1,igrid),kind=8)/max(real(m%uold(icell2,1,igrid),kind=8),real(r%smallr,kind=8))*vol2(ind)
              end if
#endif
           end do
        end if
        ! Set time-centered velocity and drift
        p%vp(ipart,1:ndim)=ff(1:ndim)
        p%xp(ipart,1:ndim)=p%xp(ipart,1:ndim)+p%vp(ipart,1:ndim)*g%dtnew(ilevel)
     endif
  end do
  call close_cache(mdl)
  if(action_part==action_kick_drift)then
     do ipart=p%headp(ilevel),p%tailp(ilevel)
        do idim=1,ndim
           if(r%periodic(idim))then
              if(p%xp(ipart,idim)< 0.0d0           )p%xp(ipart,idim)=p%xp(ipart,idim)+r%box_size(idim)
              if(p%xp(ipart,idim)>=r%box_size(idim))p%xp(ipart,idim)=p%xp(ipart,idim)-r%box_size(idim)
           endif
        end do
     end do
  end if
  end associate
end subroutine cic_trace_gas_part

subroutine tsc_trace_gas_part(s,p,ilevel,action_part)
  use amr_parameters, only: ndim, threetondim
  use pm_parameters
  use pm_commons, only: part_t
  use ramses_commons, only: ramses_t
  use nbors_utils
  use cache_commons
  use cache
  implicit none
  type(ramses_t)::s
  type(part_t)::p
  integer::ilevel
  integer::action_part
  real(kind=8),dimension(1:ndim)::x,wl,wc,wr
  real(kind=8),dimension(1:ndim)::x_mid,wl2,wc2,wr2
  integer,dimension(1:ndim)::cl,cc,cr
  integer,dimension(1:ndim)::cl2,cc2,cr2
  real(kind=8),dimension(1:threetondim)::vol,vol2
  integer,dimension(1:ndim,1:threetondim)::ckey,ckey2
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer::ipart,icell,icell2,igrid,ind,idim
  real(kind=8)::xl,xc,xr
  real(kind=8)::dx_loc
  real(kind=8),dimension(1:ndim)::ff,v_pred
  logical::ok_level
  type(msg_three_realdp)::dummy_three_realdp
  type(msg_hydro_mflux)::dummy_hydro_mflux

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

  if(p%static)return
  if (p%type/=TRAC_TYPE) return

  dx_loc=r%boxlen/2**ilevel

  call open_cache(mdl, m, pack_size=storage_size(dummy_hydro_mflux)/32, &
       pack=pack_fetch_kick_trac, unpack=unpack_fetch_kick_trac)
  do ipart=p%headp(ilevel),p%tailp(ilevel)
     do idim=1,ndim
        x(idim)=(p%xp(ipart,idim)+m%skip(idim))/dx_loc
     end do
     do idim=1,ndim
        if(r%periodic(idim))then
           if(x(idim)< dble(m%box_ckey_min(idim,ilevel+1)))x(idim)=x(idim)+dble(m%box_ckey_max(idim,ilevel+1)-m%box_ckey_min(idim,ilevel+1))
           if(x(idim)>=dble(m%box_ckey_max(idim,ilevel+1)))x(idim)=x(idim)-dble(m%box_ckey_max(idim,ilevel+1)-m%box_ckey_min(idim,ilevel+1))
        endif
     end do
     do idim=1,ndim
        cl(idim)=int(x(idim))-1
        cc(idim)=int(x(idim))
        cr(idim)=int(x(idim))+1
        xl=dble(cl(idim))+0.5D0
        xc=dble(cc(idim))+0.5D0
        xr=dble(cr(idim))+0.5D0
        wl(idim)=0.5D0*(1.5D0-abs(x(idim)-xl))**2
        wc(idim)=0.75D0-         (x(idim)-xc) **2
        wr(idim)=0.5D0*(1.5D0-abs(x(idim)-xr))**2
     end do
     do idim=1,ndim
        if(r%periodic(idim))then
           if(cl(idim)< m%box_ckey_min(idim,ilevel+1))cl(idim)=m%box_ckey_max(idim,ilevel+1)-1
           if(cr(idim)>=m%box_ckey_max(idim,ilevel+1))cr(idim)=m%box_ckey_min(idim,ilevel+1)
        endif
     enddo
     ckey = tsc_index(cl,cc,cr)
     vol = tsc_weight(wl,wc,wr)
     hash_nbor(0)=ilevel+1
     ff(1:ndim)=0.0
     do ind=1,threetondim
        hash_nbor(1:ndim)=ckey(1:ndim,ind)
        call get_parent_cell(s,hash_nbor,igrid,icell,flush_cache=.false.,fetch_cache=.true.)
#ifdef HYDRO
        if(igrid>0)then
           ff(1:ndim)=ff(1:ndim)+real(m%uold(icell,2:ndim+1,igrid),kind=8)/max(real(m%uold(icell,1,igrid),kind=8),real(r%smallr,kind=8))*vol(ind)
        end if
#endif
     end do

     if(action_part==action_kick_only)then
        ! RK2 step 1: stash v^n at x^n
        p%vp(ipart,1:ndim)=ff(1:ndim)
        p%levelp(ipart)=ilevel

     else if(action_part==action_kick_drift)then
        ! RK2 step 2: predict with v^n, sample v^{n+1} at midpoint
        if (g%nstep>0) then
           v_pred(1:ndim)=p%vp(ipart,1:ndim)
        else
           v_pred(1:ndim)=ff(1:ndim)
        endif
        do idim=1,ndim
           x_mid(idim)=x(idim)+0.5d0*g%dtnew(ilevel)*v_pred(idim)/dx_loc
        end do
        do idim=1,ndim
           if(r%periodic(idim))then
              if(x_mid(idim)< dble(m%box_ckey_min(idim,ilevel+1)))x_mid(idim)=x_mid(idim)+dble(m%box_ckey_max(idim,ilevel+1)-m%box_ckey_min(idim,ilevel+1))
              if(x_mid(idim)>=dble(m%box_ckey_max(idim,ilevel+1)))x_mid(idim)=x_mid(idim)-dble(m%box_ckey_max(idim,ilevel+1)-m%box_ckey_min(idim,ilevel+1))
           endif
        end do
        do idim=1,ndim
           cl2(idim)=int(x_mid(idim))-1
           cc2(idim)=int(x_mid(idim))
           cr2(idim)=int(x_mid(idim))+1
           xl=dble(cl2(idim))+0.5D0
           xc=dble(cc2(idim))+0.5D0
           xr=dble(cr2(idim))+0.5D0
           wl2(idim)=0.5D0*(1.5D0-abs(x_mid(idim)-xl))**2
           wc2(idim)=0.75D0-         (x_mid(idim)-xc) **2
           wr2(idim)=0.5D0*(1.5D0-abs(x_mid(idim)-xr))**2
        end do
        do idim=1,ndim
           if(r%periodic(idim))then
              if(cl2(idim)< m%box_ckey_min(idim,ilevel+1))cl2(idim)=m%box_ckey_max(idim,ilevel+1)-1
              if(cr2(idim)>=m%box_ckey_max(idim,ilevel+1))cr2(idim)=m%box_ckey_min(idim,ilevel+1)
           endif
        enddo
        ckey2 = tsc_index(cl2,cc2,cr2)
        vol2 = tsc_weight(wl2,wc2,wr2)
        hash_nbor(0)=ilevel+1
        ff(1:ndim)=0.0
        do ind=1,threetondim
           hash_nbor(1:ndim)=ckey2(1:ndim,ind)
           call get_parent_cell(s,hash_nbor,igrid,icell2,flush_cache=.false.,fetch_cache=.true.)
#ifdef HYDRO
           if(igrid>0)then
              ff(1:ndim)=ff(1:ndim)+real(m%uold(icell2,2:ndim+1,igrid),kind=8)/max(real(m%uold(icell2,1,igrid),kind=8),real(r%smallr,kind=8))*vol2(ind)
           end if
#endif
        end do
        p%vp(ipart,1:ndim)=ff(1:ndim)
        p%xp(ipart,1:ndim)=p%xp(ipart,1:ndim)+p%vp(ipart,1:ndim)*g%dtnew(ilevel)
     endif
  end do
  call close_cache(mdl)
  if(action_part==action_kick_drift)then
     do ipart=p%headp(ilevel),p%tailp(ilevel)
        do idim=1,ndim
           if(r%periodic(idim))then
              if(p%xp(ipart,idim) <0.0d0           )p%xp(ipart,idim)=p%xp(ipart,idim)+r%box_size(idim)
              if(p%xp(ipart,idim)>=r%box_size(idim))p%xp(ipart,idim)=p%xp(ipart,idim)-r%box_size(idim)
           endif
        end do
     end do
  end if
  end associate
end subroutine tsc_trace_gas_part

subroutine pcs_trace_gas_part(s,p,ilevel,action_part)
  use amr_parameters, only: ndim, fourtondim
  use pm_parameters
  use pm_commons, only: part_t
  use ramses_commons, only: ramses_t
  use nbors_utils
  use cache_commons
  use cache
  implicit none
  type(ramses_t)::s
  type(part_t)::p
  integer::ilevel
  integer::action_part
  real(kind=8),dimension(1:ndim)::x,wll,wl,wr,wrr
  real(kind=8),dimension(1:ndim)::x_mid,wll2,wl2,wr2,wrr2
  integer,dimension(1:ndim)::cll,cl,cr,crr
  integer,dimension(1:ndim)::cll2,cl2,cr2,crr2
  real(kind=8),dimension(1:fourtondim)::vol,vol2
  integer,dimension(1:ndim,1:fourtondim)::ckey,ckey2
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer::ipart,icell,icell2,igrid,ind,idim
  real(kind=8)::xll,xl,xr,xrr
  real(kind=8)::dx_loc
  real(kind=8),dimension(1:ndim)::ff,v_pred
  type(msg_large_realdp)::dummy_large_realdp
  type(msg_hydro_mflux)::dummy_hydro_mflux

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

  if(p%static)return
  if (p%type/=TRAC_TYPE) return

  dx_loc=r%boxlen/2**ilevel

  call open_cache(mdl, m, pack_size=storage_size(dummy_hydro_mflux)/32, &
       pack=pack_fetch_kick_trac, unpack=unpack_fetch_kick_trac)
  do ipart=p%headp(ilevel),p%tailp(ilevel)
     do idim=1,ndim
        x(idim)=(p%xp(ipart,idim)+m%skip(idim))/dx_loc
     end do
     do idim=1,ndim
        if(r%periodic(idim))then
           if(x(idim)< dble(m%box_ckey_min(idim,ilevel+1)))x(idim)=x(idim)+dble(m%box_ckey_max(idim,ilevel+1)-m%box_ckey_min(idim,ilevel+1))
           if(x(idim)>=dble(m%box_ckey_max(idim,ilevel+1)))x(idim)=x(idim)-dble(m%box_ckey_max(idim,ilevel+1)-m%box_ckey_min(idim,ilevel+1))
        endif
     end do
     do idim=1,ndim
        crr(idim)=int(x(idim)+1.5D0)
        cr (idim)=crr(idim)-1
        cl (idim)=crr(idim)-2
        cll(idim)=crr(idim)-3
        xll=dble(cll(idim))+0.5D0
        xl =dble(cl (idim))+0.5D0
        xr =dble(cr (idim))+0.5D0
        xrr=dble(crr(idim))+0.5D0
        wll(idim)=(2D0-abs(x(idim)-xll))**3/6D0
        wl (idim)=(4D0-6D0*(x(idim)-xl)**2+3d0*abs(x(idim)-xl )**3)/6D0
        wr (idim)=(4D0-6D0*(x(idim)-xr)**2+3d0*abs(x(idim)-xr )**3)/6D0
        wrr(idim)=(2D0-abs(x(idim)-xrr))**3/6D0
     end do
     do idim=1,ndim
        if(r%periodic(idim))then
           if(cll(idim)< m%box_ckey_min(idim,ilevel+1))cll(idim)=cll(idim)-m%box_ckey_min(idim,ilevel+1)+m%box_ckey_max(idim,ilevel+1)
           if(cl (idim)< m%box_ckey_min(idim,ilevel+1))cl (idim)=cl (idim)-m%box_ckey_min(idim,ilevel+1)+m%box_ckey_max(idim,ilevel+1)
           if(cr (idim)>=m%box_ckey_max(idim,ilevel+1))cr (idim)=cr (idim)+m%box_ckey_min(idim,ilevel+1)-m%box_ckey_max(idim,ilevel+1)
           if(crr(idim)>=m%box_ckey_max(idim,ilevel+1))crr(idim)=crr(idim)+m%box_ckey_min(idim,ilevel+1)-m%box_ckey_max(idim,ilevel+1)
        endif
     enddo
     ckey = pcs_index(cll,cl,cr,crr)
     vol = pcs_weight(wll,wl,wr,wrr)
     hash_nbor(0)=ilevel+1
     ff(1:ndim)=0.0
     do ind=1,fourtondim
        hash_nbor(1:ndim)=ckey(1:ndim,ind)
        call get_parent_cell(s,hash_nbor,igrid,icell,flush_cache=.false.,fetch_cache=.true.)
#ifdef HYDRO
        if(igrid>0)then
           ff(1:ndim)=ff(1:ndim)+real(m%uold(icell,2:ndim+1,igrid),kind=8)/max(real(m%uold(icell,1,igrid),kind=8),real(r%smallr,kind=8))*vol(ind)
        end if
#endif
     end do

     if(action_part==action_kick_only)then
        ! RK2 step 1: stash v^n at x^n
        p%vp(ipart,1:ndim)=ff(1:ndim)
        p%levelp(ipart)=ilevel

     else if(action_part==action_kick_drift)then
        ! RK2 step 2: predict with v^n, sample v^{n+1} at midpoint
        if (g%nstep>0) then
           v_pred(1:ndim)=p%vp(ipart,1:ndim)
        else
           v_pred(1:ndim)=ff(1:ndim)
        endif
        do idim=1,ndim
           x_mid(idim)=x(idim)+0.5d0*g%dtnew(ilevel)*v_pred(idim)/dx_loc
        end do
        do idim=1,ndim
           if(r%periodic(idim))then
              if(x_mid(idim)< dble(m%box_ckey_min(idim,ilevel+1)))x_mid(idim)=x_mid(idim)+dble(m%box_ckey_max(idim,ilevel+1)-m%box_ckey_min(idim,ilevel+1))
              if(x_mid(idim)>=dble(m%box_ckey_max(idim,ilevel+1)))x_mid(idim)=x_mid(idim)-dble(m%box_ckey_max(idim,ilevel+1)-m%box_ckey_min(idim,ilevel+1))
           endif
        end do
        do idim=1,ndim
           crr2(idim)=int(x_mid(idim)+1.5D0)
           cr2 (idim)=crr2(idim)-1
           cl2 (idim)=crr2(idim)-2
           cll2(idim)=crr2(idim)-3
           xll=dble(cll2(idim))+0.5D0
           xl =dble(cl2 (idim))+0.5D0
           xr =dble(cr2 (idim))+0.5D0
           xrr=dble(crr2(idim))+0.5D0
           wll2(idim)=(2D0-abs(x_mid(idim)-xll))**3/6D0
           wl2 (idim)=(4D0-6D0*(x_mid(idim)-xl)**2+3d0*abs(x_mid(idim)-xl )**3)/6D0
           wr2 (idim)=(4D0-6D0*(x_mid(idim)-xr)**2+3d0*abs(x_mid(idim)-xr )**3)/6D0
           wrr2(idim)=(2D0-abs(x_mid(idim)-xrr))**3/6D0
        end do
        do idim=1,ndim
           if(r%periodic(idim))then
              if(cll2(idim)< m%box_ckey_min(idim,ilevel+1))cll2(idim)=cll2(idim)-m%box_ckey_min(idim,ilevel+1)+m%box_ckey_max(idim,ilevel+1)
              if(cl2 (idim)< m%box_ckey_min(idim,ilevel+1))cl2 (idim)=cl2 (idim)-m%box_ckey_min(idim,ilevel+1)+m%box_ckey_max(idim,ilevel+1)
              if(cr2 (idim)>=m%box_ckey_max(idim,ilevel+1))cr2 (idim)=cr2 (idim)+m%box_ckey_min(idim,ilevel+1)-m%box_ckey_max(idim,ilevel+1)
              if(crr2(idim)>=m%box_ckey_max(idim,ilevel+1))crr2(idim)=crr2(idim)+m%box_ckey_min(idim,ilevel+1)-m%box_ckey_max(idim,ilevel+1)
           endif
        enddo
        ckey2 = pcs_index(cll2,cl2,cr2,crr2)
        vol2 = pcs_weight(wll2,wl2,wr2,wrr2)
        hash_nbor(0)=ilevel+1
        ff(1:ndim)=0.0
        do ind=1,fourtondim
           hash_nbor(1:ndim)=ckey2(1:ndim,ind)
           call get_parent_cell(s,hash_nbor,igrid,icell2,flush_cache=.false.,fetch_cache=.true.)
#ifdef HYDRO
           if(igrid>0)then
              ff(1:ndim)=ff(1:ndim)+real(m%uold(icell2,2:ndim+1,igrid),kind=8)/max(real(m%uold(icell2,1,igrid),kind=8),real(r%smallr,kind=8))*vol2(ind)
           end if
#endif
        end do
        p%vp(ipart,1:ndim)=ff(1:ndim)
        p%xp(ipart,1:ndim)=p%xp(ipart,1:ndim)+p%vp(ipart,1:ndim)*g%dtnew(ilevel)
     endif
  end do
  call close_cache(mdl)
  if(action_part==action_kick_drift)then
     do ipart=p%headp(ilevel),p%tailp(ilevel)
        do idim=1,ndim
           if(r%periodic(idim))then
              if(p%xp(ipart,idim)< 0.0d0           )p%xp(ipart,idim)=p%xp(ipart,idim)+r%box_size(idim)
              if(p%xp(ipart,idim)>=r%box_size(idim))p%xp(ipart,idim)=p%xp(ipart,idim)-r%box_size(idim)
           endif
        end do
     end do
  end if
  end associate
end subroutine pcs_trace_gas_part
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cic_trace_gas_part_sgs_turb(s,p,ilevel,action_part)
   use amr_parameters, only: ndim, twotondim, dp
   use pm_parameters
   use pm_commons, only: part_t
   use ramses_commons, only: ramses_t
   use rng
   use nbors_utils
   use cache_commons
   use cache
   implicit none
   type(ramses_t)::s
   type(part_t)::p
   integer::ilevel
   integer::action_part
   real(kind=8),dimension(1:ndim)::x,x_mid,v_pred,u_eff,u_mid,disp,xi,momentum
   real(kind=8),dimension(1:ndim)::grad_at_part
   real(kind=8),dimension(1:twotondim)::vol,kappa_cells
   integer,dimension(1:ndim)::il,ir
   integer,dimension(1:ndim,1:twotondim)::ckey
   integer(kind=8),dimension(0:ndim)::hash_nbor
   real(kind=8),dimension(1:ndim,1:2)::w1d,dw1d
   real(kind=8)::dx_loc,dt_level,kappa_mid,noise_amp,rho
   logical :: use_sgs
   type(msg_hydro_mflux)::dummy_hydro_mflux
   type(RngStream),external::RngStream_CreateStream
   real(kind=8),external::RngStream_RandUni
   integer(kind=8)::stream_skip
   external :: RngStream_SetPackageSeed, RngStream_AdvanceState, gaussdev
   integer :: ipart,idim,ind,icell,igrid,ii
   character(LEN=80) :: filename,fileloc
   character(LEN=5) :: nchar

   associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)
   if(p%static)return
   if (p%type/=TRAC_TYPE) return

   dx_loc=r%boxlen/2**ilevel

   dt_level=g%dtnew(ilevel)
   use_sgs = r%sgs_turb .and. (r%iturb>0)

   if(action_part==action_kick_only)then
      do ipart=p%headp(ilevel),p%tailp(ilevel)
         p%levelp(ipart)=ilevel
      end do
      return
   endif

   if(use_sgs .and. .not. tracer_rng_ready)then
      call RngStream_SetPackageSeed(r%seed)
      tracer_rng = RngStream_CreateStream('tracer_sgs')
      stream_skip = int(2*g%myid,kind=8)
      call RngStream_AdvanceState(tracer_rng,0_8,stream_skip)
      tracer_rng_ready = .true.
   end if

   call open_cache(mdl, m, pack_size=storage_size(dummy_hydro_mflux)/32, &
        pack=pack_fetch_kick_trac, unpack=unpack_fetch_kick_trac)

   do ipart=p%headp(ilevel),p%tailp(ilevel)

      do idim=1,ndim
         x(idim)=(p%xp(ipart,idim)+m%skip(idim))/dx_loc
      end do
      call wrap_cell_coords(s,x,ilevel+1)

      if (g%nstep>0) then
         v_pred(1:ndim)=p%vp(ipart,1:ndim)
      else
         call gather_cic_state(s,x,ilevel,dx_loc,.false.,u_eff,kappa_mid)
         v_pred(1:ndim)=u_eff(1:ndim)
      endif

      do idim=1,ndim
         x_mid(idim)=x(idim)+0.5d0*dt_level*v_pred(idim)/dx_loc
      end do
      call wrap_cell_coords(s,x_mid,ilevel+1)

      call cic_weights_and_derivs(x_mid, w1d, dw1d, il, ir, vol)
      do idim=1,ndim
         if(r%periodic(idim))then
            if(il(idim)< m%box_ckey_min(idim,ilevel+1))il(idim)=m%box_ckey_max(idim,ilevel+1)-1
            if(ir(idim)>=m%box_ckey_max(idim,ilevel+1))ir(idim)=m%box_ckey_min(idim,ilevel+1)
         endif
      end do
      ckey = cic_index(il,ir)

      momentum(1:ndim)=0.d0
      rho=0.d0
      kappa_mid=0.d0
      kappa_cells(1:twotondim)=0.d0
      hash_nbor(0)=ilevel+1
      do ind=1,twotondim
         hash_nbor(1:ndim)=ckey(1:ndim,ind)
         call get_parent_cell(s,hash_nbor,igrid,icell,flush_cache=.false.,fetch_cache=.true.)
#ifdef HYDRO
         if(igrid>0)then
            momentum(1:ndim)=momentum(1:ndim)+m%uold(icell,2:ndim+1,igrid)*vol(ind)
            rho=rho+m%uold(icell,1,igrid)*vol(ind)
            if(use_sgs)then
               kappa_cells(ind)=tracer_cell_kappa(m%uold(icell,1,igrid),m%uold(icell,r%iturb,igrid),real(dx_loc,dp),real(r%smallr,dp))
               kappa_mid=kappa_mid+kappa_cells(ind)*vol(ind)
            end if
         end if
#endif
      end do
      if(rho>r%smallr)then
         u_mid(1:ndim)=momentum(1:ndim)/rho
      else
         u_mid(1:ndim)=0.d0
      end if

      grad_at_part(1:ndim)=0.d0
      if(use_sgs)then
         call compute_gradient_cic_scalar(w1d, dw1d, kappa_cells, grad_at_part)
      end if

      disp(1:ndim)=(u_mid(1:ndim) + grad_at_part(1:ndim)/dx_loc)*dt_level

      if(use_sgs .and. kappa_mid>0.0d0)then
         if(r%tracer_kick_pdf=='gaussian')then
            call sample_tracer_gaussian(xi)
         elseif(r%tracer_kick_pdf=='uniform')then
            call sample_tracer_uniform(xi)
         else
            call sample_tracer_uniform(xi)
         endif
         noise_amp = sqrt(2.0d0*kappa_mid*dt_level)
         disp(1:ndim)=disp(1:ndim)+noise_amp*xi(1:ndim)
      end if

      p%levelp(ipart)=ilevel
      p%vp(ipart,1:ndim)=u_mid(1:ndim)
      p%xp(ipart,1:ndim)=p%xp(ipart,1:ndim)+disp(1:ndim)
      if(s%r%ntrajectories>0)then
         do ii=1,s%r%ntrajectories
            if(s%r%trajectories(ii)==p%idp(ipart))then
               call title(g%myid,nchar)
               filename='trajectory.dat'
               fileloc=TRIM(filename)//TRIM(nchar)
               open(25+g%myid,file=fileloc,status='unknown',access='append')
               write(25+g%myid,*) g%t, p%idp(ipart), &
                    (p%xp(ipart,idim), idim=1,ndim), &
                    (p%vp(ipart,idim), idim=1,ndim)
               close(25+g%myid)
               exit
            endif
         end do
      endif

   end do

   call close_cache(mdl)

   if(action_part==action_kick_drift)then
      do ipart=p%headp(ilevel),p%tailp(ilevel)
         do idim=1,ndim
            if(r%periodic(idim))then
               if(p%xp(ipart,idim)< 0.0d0           )p%xp(ipart,idim)=p%xp(ipart,idim)+r%box_size(idim)
               if(p%xp(ipart,idim)>=r%box_size(idim))p%xp(ipart,idim)=p%xp(ipart,idim)-r%box_size(idim)
            endif
         end do
      end do
   end if

   end associate
 end subroutine cic_trace_gas_part_sgs_turb

 subroutine tsc_trace_gas_part_sgs_turb(s,p,ilevel,action_part)
   use amr_parameters, only: ndim, threetondim, dp
   use pm_parameters
   use pm_commons, only: part_t
   use ramses_commons, only: ramses_t
   use rng
   use nbors_utils
   use cache_commons
   use cache
   use rho_fine_module, only: tsc_index
   implicit none
   type(ramses_t)::s
   type(part_t)::p
   integer::ilevel
   integer::action_part
   real(kind=8),dimension(1:ndim)::x,x_mid,v_pred,u_eff,u_mid,disp,xi,momentum
   real(kind=8),dimension(1:ndim)::grad_at_part
   real(kind=8),dimension(1:threetondim)::vol,kappa_cells
   integer,dimension(1:ndim)::il,ic,ir
   integer,dimension(1:ndim,1:threetondim)::ckey
   integer(kind=8),dimension(0:ndim)::hash_nbor
   real(kind=8),dimension(1:ndim,1:3)::w1d,dw1d
   real(kind=8)::dx_loc,dt_level,kappa_mid,noise_amp,rho
   logical :: use_sgs
   type(msg_hydro_mflux)::dummy_hydro_mflux
   type(RngStream),external::RngStream_CreateStream
   real(kind=8),external::RngStream_RandUni
   integer(kind=8)::stream_skip
   external :: RngStream_SetPackageSeed, RngStream_AdvanceState, gaussdev
   integer :: ipart,idim,ind,icell,igrid,ii
   character(LEN=80) :: filename,fileloc
   character(LEN=5) :: nchar

   associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)
   if(p%static)return
   if (p%type/=TRAC_TYPE) return

   dx_loc=r%boxlen/2**ilevel

   dt_level=g%dtnew(ilevel)
   use_sgs = r%sgs_turb .and. (r%iturb>0)

   if(action_part==action_kick_only)then
      do ipart=p%headp(ilevel),p%tailp(ilevel)
         p%levelp(ipart)=ilevel
      end do
      return
   endif

   if(use_sgs .and. .not. tracer_rng_ready)then
      call RngStream_SetPackageSeed(r%seed)
      tracer_rng = RngStream_CreateStream('tracer_sgs')
      stream_skip = int(2*g%myid,kind=8)
      call RngStream_AdvanceState(tracer_rng,0_8,stream_skip)
      tracer_rng_ready = .true.
   end if

   call open_cache(mdl, m, pack_size=storage_size(dummy_hydro_mflux)/32, &
        pack=pack_fetch_kick_trac, unpack=unpack_fetch_kick_trac)

   do ipart=p%headp(ilevel),p%tailp(ilevel)

      do idim=1,ndim
         x(idim)=(p%xp(ipart,idim)+m%skip(idim))/dx_loc
      end do
      call wrap_cell_coords(s,x,ilevel+1)

      if (g%nstep>0) then
         v_pred(1:ndim)=p%vp(ipart,1:ndim)
      else
         call gather_tsc_state(s,x,ilevel,dx_loc,.false.,u_eff,kappa_mid)
         v_pred(1:ndim)=u_eff(1:ndim)
      endif

      do idim=1,ndim
         x_mid(idim)=x(idim)+0.5d0*dt_level*v_pred(idim)/dx_loc
      end do
      call wrap_cell_coords(s,x_mid,ilevel+1)

      call tsc_weights_and_derivs(x_mid, w1d, dw1d, il, ic, ir, vol)
      do idim=1,ndim
         if(r%periodic(idim))then
            if(il(idim)< m%box_ckey_min(idim,ilevel+1))il(idim)=il(idim)-m%box_ckey_min(idim,ilevel+1)+m%box_ckey_max(idim,ilevel+1)
            if(ir(idim)>=m%box_ckey_max(idim,ilevel+1))ir(idim)=ir(idim)+m%box_ckey_min(idim,ilevel+1)-m%box_ckey_max(idim,ilevel+1)
         endif
      end do
      ckey = tsc_index(il,ic,ir)

      momentum(1:ndim)=0.d0
      rho=0.d0
      kappa_mid=0.d0
      kappa_cells(1:threetondim)=0.d0
      hash_nbor(0)=ilevel+1
      do ind=1,threetondim
         hash_nbor(1:ndim)=ckey(1:ndim,ind)
         call get_parent_cell(s,hash_nbor,igrid,icell,flush_cache=.false.,fetch_cache=.true.)
#ifdef HYDRO
         if(igrid>0)then
            momentum(1:ndim)=momentum(1:ndim)+m%uold(icell,2:ndim+1,igrid)*vol(ind)
            rho=rho+m%uold(icell,1,igrid)*vol(ind)
            if(use_sgs)then
               kappa_cells(ind)=tracer_cell_kappa(m%uold(icell,1,igrid),m%uold(icell,r%iturb,igrid),real(dx_loc,dp),real(r%smallr,dp))
               kappa_mid=kappa_mid+kappa_cells(ind)*vol(ind)
            end if
         end if
#endif
      end do
      if(rho>r%smallr)then
         u_mid(1:ndim)=momentum(1:ndim)/rho
      else
         u_mid(1:ndim)=0.d0
      end if

      grad_at_part(1:ndim)=0.d0
      if(use_sgs)then
         call compute_gradient_tsc_scalar(w1d, dw1d, kappa_cells, grad_at_part)
      end if

      disp(1:ndim)=(u_mid(1:ndim) + grad_at_part(1:ndim)/dx_loc)*dt_level

      if(use_sgs .and. kappa_mid>0.0d0)then
         if(r%tracer_kick_pdf=='gaussian')then
            call sample_tracer_gaussian(xi)
         elseif(r%tracer_kick_pdf=='uniform')then
            call sample_tracer_uniform(xi)
         else
            call sample_tracer_uniform(xi)
         endif
         noise_amp = sqrt(2.0d0*kappa_mid*dt_level)
         disp(1:ndim)=disp(1:ndim)+noise_amp*xi(1:ndim)
      end if

      p%levelp(ipart)=ilevel
      p%vp(ipart,1:ndim)=u_mid(1:ndim)
      p%xp(ipart,1:ndim)=p%xp(ipart,1:ndim)+disp(1:ndim)
      if(s%r%ntrajectories>0)then
         do ii=1,s%r%ntrajectories
            if(s%r%trajectories(ii)==p%idp(ipart))then
               call title(g%myid,nchar)
               filename='trajectory.dat'
               fileloc=TRIM(filename)//TRIM(nchar)
               open(25+g%myid,file=fileloc,status='unknown',access='append')
               write(25+g%myid,*) g%t, p%idp(ipart), &
                    (p%xp(ipart,idim), idim=1,ndim), &
                    (p%vp(ipart,idim), idim=1,ndim)
               close(25+g%myid)
               exit
            endif
         end do
      endif

   end do

   call close_cache(mdl)

   if(action_part==action_kick_drift)then
      do ipart=p%headp(ilevel),p%tailp(ilevel)
         do idim=1,ndim
            if(r%periodic(idim))then
               if(p%xp(ipart,idim)< 0.0d0           )p%xp(ipart,idim)=p%xp(ipart,idim)+r%box_size(idim)
               if(p%xp(ipart,idim)>=r%box_size(idim))p%xp(ipart,idim)=p%xp(ipart,idim)-r%box_size(idim)
            endif
         end do
      end do
   end if

   end associate
 end subroutine tsc_trace_gas_part_sgs_turb

 subroutine mc_trace_gas_part(s,p,ilevel,action_part)
   ! Classic Monte Carlo tracer (scheme 0)
   use amr_parameters, only: ndim, twotondim
   use pm_parameters
   use pm_commons, only: part_t
   use ramses_commons, only: ramses_t
   use rng
   use nbors_utils
   use cache_commons
   use cache
   implicit none
   type(ramses_t)::s
   type(part_t)::p
   integer::ilevel
   integer::action_part
   real(kind=8),dimension(1:ndim)::x,vel
   integer,dimension(1:ndim)::icell_idx
   integer(kind=8),dimension(0:ndim)::hash_nbor
   real(kind=8),dimension(1:2*ndim)::prob_face
   real(kind=8)::out_sum,scale,stay_prob,u,cum
   integer::ipart,idim,iface,selected,icell
   real(kind=8)::rho_cell,denom,dx_loc,dist_to_face,tol_corner,dist_to_corner
   logical::near_corner
   integer,dimension(1:ndim)::corner_idx
   integer,dimension(1:ndim,1:twotondim)::corner_nbor_idx
   real(kind=8),dimension(1:twotondim)::corner_weight
   real(kind=8)::weight_sum
   integer::ind,bit,selected_corner,corner_associated
   type(msg_hydro_mflux)::dummy_hydro_mflux
   type(RngStream),external::RngStream_CreateStream
   real(kind=8),external::RngStream_RandUni
   integer(kind=8)::stream_skip
   external :: RngStream_SetPackageSeed, RngStream_AdvanceState
   integer :: igrid,ii
   character(LEN=80) :: filename,fileloc
   character(LEN=5) :: nchar

   associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)
   if(p%static)return
   if (p%type/=TRAC_TYPE) return

   dx_loc=r%boxlen/2**ilevel

   tol_corner=0.05d0

   if(action_part==action_kick_only)then
      do ipart=p%headp(ilevel),p%tailp(ilevel)
         p%levelp(ipart)=ilevel
      end do
      return
   endif

   if(.not. tracer_rng_ready)then
      call RngStream_SetPackageSeed(r%seed)
      tracer_rng = RngStream_CreateStream('tracer_mc')
      stream_skip = int(2*g%myid,kind=8)
      call RngStream_AdvanceState(tracer_rng,0_8,stream_skip)
      tracer_rng_ready = .true.
   end if

   call open_cache(mdl, m, pack_size=storage_size(dummy_hydro_mflux)/32, &
        pack=pack_fetch_kick_trac, unpack=unpack_fetch_kick_trac)

   do ipart=p%headp(ilevel),p%tailp(ilevel)
      do idim=1,ndim
         x(idim)=(p%xp(ipart,idim)+m%skip(idim))/dx_loc
      end do
      call wrap_cell_coords(s,x,ilevel+1)

      near_corner=.true.
      do idim=1,ndim
         dist_to_corner=abs(x(idim)-dble(nint(x(idim))))
         if(dist_to_corner>=tol_corner)near_corner=.false.
      end do

      if(near_corner)then
         do idim=1,ndim
            corner_idx(idim)=nint(x(idim))
         end do
         corner_weight=0.d0
         corner_associated=0
         do ind=1,twotondim
            do idim=1,ndim
               bit=merge(1,0,btest(ind-1,idim-1))
               corner_nbor_idx(idim,ind)=corner_idx(idim)+bit-1
               if(r%periodic(idim))then
                  if(corner_nbor_idx(idim,ind)< m%box_ckey_min(idim,ilevel+1))corner_nbor_idx(idim,ind)=m%box_ckey_max(idim,ilevel+1)-1
                  if(corner_nbor_idx(idim,ind)>=m%box_ckey_max(idim,ilevel+1))corner_nbor_idx(idim,ind)=m%box_ckey_min(idim,ilevel+1)
               end if
            end do
            hash_nbor(0)=ilevel+1
            hash_nbor(1:ndim)=corner_nbor_idx(1:ndim,ind)
            call get_parent_cell(s,hash_nbor,igrid,icell,flush_cache=.false.,fetch_cache=.true.)
#ifdef HYDRO
            if(igrid>0)then
               corner_associated=corner_associated+1
               corner_weight(ind)=max(real(m%uold(icell,1,igrid),kind=8),real(r%smallr,kind=8))
            else
               corner_weight(ind)=0.d0
            end if
#else
            if(igrid>0)then
               corner_associated=corner_associated+1
               corner_weight(ind)=1.d0
            else
               corner_weight(ind)=0.d0
            end if
#endif
         end do
         if(corner_associated<twotondim)near_corner=.false.
         weight_sum=0.d0
         do ind=1,twotondim
            weight_sum=weight_sum+corner_weight(ind)
         end do
         if(near_corner.and.weight_sum>0.d0)then
            u=RngStream_RandUni(tracer_rng)
            cum=0.d0
            selected_corner=1
            do ind=1,twotondim
               cum=cum+corner_weight(ind)/weight_sum
               if(u<=cum)then
                  selected_corner=ind
                  exit
               end if
            end do
            do idim=1,ndim
               x(idim)=dble(corner_nbor_idx(idim,selected_corner))+0.5d0
            end do
            call wrap_cell_coords(s,x,ilevel+1)
            do idim=1,ndim
               p%xp(ipart,idim)=x(idim)*dx_loc-m%skip(idim)
            end do
         else
            near_corner=.false.
         end if
      end if

      if(.not.near_corner)then
         do idim=1,ndim
            dist_to_face=abs(x(idim)-dble(nint(x(idim))))
            if(dist_to_face<tol_corner)then
               if(RngStream_RandUni(tracer_rng)<0.5d0)then
                  x(idim)=dble(nint(x(idim)))-0.5d0
               else
                  x(idim)=dble(nint(x(idim)))+0.5d0
               endif
            else
               x(idim)=dble(int(x(idim)))+0.5d0
            endif
         end do
         call wrap_cell_coords(s,x,ilevel+1)
         do idim=1,ndim
            p%xp(ipart,idim)=x(idim)*dx_loc-m%skip(idim)
         end do
      end if

      do idim=1,ndim
         icell_idx(idim)=int(x(idim))
         if(r%periodic(idim))then
            if(icell_idx(idim)< m%box_ckey_min(idim,ilevel+1))icell_idx(idim)=m%box_ckey_max(idim,ilevel+1)-1
            if(icell_idx(idim)>=m%box_ckey_max(idim,ilevel+1))icell_idx(idim)=m%box_ckey_min(idim,ilevel+1)
         endif
      enddo

      hash_nbor(0)=ilevel+1
      hash_nbor(1:ndim)=icell_idx(1:ndim)
      call get_parent_cell(s,hash_nbor,igrid,icell,flush_cache=.false.,fetch_cache=.true.)
#ifdef HYDRO
      if(igrid<=0)cycle
#ifdef TRCFLX
      rho_cell=max(m%mflux(icell,1,igrid),r%smallr)
      denom=max(rho_cell,r%smallr)
      vel=0.d0
      vel(1:ndim)=m%uold(icell,2:ndim+1,igrid)/max(m%uold(icell,1,igrid),r%smallr)

      prob_face=0.d0
      do idim=1,ndim
         prob_face(2*idim-1)=max(-m%mflux(icell,1+idim,igrid),0.d0)/denom
         prob_face(2*idim  )=max( m%mflux(icell,1+idim+ndim,igrid),0.d0)/denom
      end do
#endif
#else
      if(igrid<=0)cycle
      prob_face=0.d0
      vel=0.d0
      denom=1.d0
#endif

      out_sum=0.d0
      do iface=1,2*ndim
         if(prob_face(iface)>0.d0)out_sum=out_sum+prob_face(iface)
      end do

      if(out_sum>1.d0)then
         scale=1.d0/out_sum
         prob_face=prob_face*scale
         out_sum=1.d0
      end if

      stay_prob=max(0.d0,1.d0-out_sum)

      selected=0
      u=RngStream_RandUni(tracer_rng)
      cum=0.d0
      do iface=1,2*ndim
         cum=cum+prob_face(iface)
         if(u<cum)then
            selected=iface
            exit
         endif
      end do

      p%vp(ipart,1:ndim)=vel(1:ndim)
      p%levelp(ipart)=ilevel

      if(selected>0)then
         idim=(selected+1)/2
         p%xp(ipart,idim)=p%xp(ipart,idim)+(-1)**selected*dx_loc
      endif
      if(s%r%ntrajectories>0)then
         do ii=1,s%r%ntrajectories
            if(s%r%trajectories(ii)==p%idp(ipart))then
               call title(g%myid,nchar)
               filename='trajectory.dat'
               fileloc=TRIM(filename)//TRIM(nchar)
               open(25+g%myid,file=fileloc,status='unknown',access='append')
               write(25+g%myid,*) g%t, p%idp(ipart), &
                    (p%xp(ipart,idim), idim=1,ndim), &
                    (p%vp(ipart,idim), idim=1,ndim)
               close(25+g%myid)
               exit
            endif
         end do
      endif
   end do

   call close_cache(mdl)

   if(action_part==action_kick_drift)then
      do ipart=p%headp(ilevel),p%tailp(ilevel)
         do idim=1,ndim
            if(r%periodic(idim))then
               if(p%xp(ipart,idim)< 0.0d0           )p%xp(ipart,idim)=p%xp(ipart,idim)+r%box_size(idim)
               if(p%xp(ipart,idim)>=r%box_size(idim))p%xp(ipart,idim)=p%xp(ipart,idim)-r%box_size(idim)
            endif
         end do
      end do
   end if

   end associate
 end subroutine mc_trace_gas_part

 subroutine cic_trace_gas_part_ito_mc(s,p,ilevel,action_part)
   ! Ito MC flux-based tracer with CIC (scheme 4)
   use amr_parameters, only: ndim, twotondim
   use pm_parameters
   use pm_commons, only: part_t
   use ramses_commons, only: ramses_t
   use rng
   use nbors_utils
   use cache_commons
   use cache
   implicit none
   type(ramses_t)::s
   type(part_t)::p
   integer::ilevel
   integer::action_part
   real(kind=8),dimension(1:ndim)::x,disp,xi,u_eff,kappa_num
   integer,dimension(1:ndim)::il,ir
   integer,dimension(1:ndim,1:twotondim)::ckey
   integer(kind=8),dimension(0:ndim)::hash_nbor
   real(kind=8),dimension(1:twotondim)::vol
   real(kind=8),dimension(1:ndim,1:twotondim)::u_cells,kappa_num_cells,skew_cells
   real(kind=8),dimension(1:ndim)::skewness_eff
   real(kind=8),dimension(1:ndim,1:2)::w1d,dw1d
   real(kind=8)::dx_loc,dt_level,rho_cell,denom,fluxL,fluxR,noise_amp,cfl_plus,cfl_minus,pr,pl
   type(msg_hydro_mflux)::dummy_hydro_mflux
   type(RngStream),external::RngStream_CreateStream
   real(kind=8),external::RngStream_RandUni
   integer(kind=8)::stream_skip
   external :: RngStream_SetPackageSeed, RngStream_AdvanceState
   integer :: ipart,idim,ind,icell,igrid,ii
   character(LEN=80) :: filename,fileloc
   character(LEN=5) :: nchar

   associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

   if(p%static)return
   if (p%type/=TRAC_TYPE) return

   dx_loc=r%boxlen/2**ilevel

   if(action_part==action_kick_only)then
      do ipart=p%headp(ilevel),p%tailp(ilevel)
         p%levelp(ipart)=ilevel
      end do
      return
   endif

   dt_level=g%dtnew(ilevel)

   if(.not. tracer_rng_ready)then
      call RngStream_SetPackageSeed(r%seed)
      tracer_rng = RngStream_CreateStream('tracer_ito_mc')
      stream_skip = int(2*g%myid,kind=8)
      call RngStream_AdvanceState(tracer_rng,0_8,stream_skip)
      tracer_rng_ready = .true.
   end if

   call open_cache(mdl, m, pack_size=storage_size(dummy_hydro_mflux)/32, &
        pack=pack_fetch_kick_trac, unpack=unpack_fetch_kick_trac)

   do ipart=p%headp(ilevel),p%tailp(ilevel)
      do idim=1,ndim
         x(idim)=(p%xp(ipart,idim)+m%skip(idim))/dx_loc
      end do
      call wrap_cell_coords(s,x,ilevel+1)

      call cic_weights_and_derivs(x, w1d, dw1d, il, ir, vol)
      do idim=1,ndim
         if(r%periodic(idim))then
            if(il(idim)< m%box_ckey_min(idim,ilevel+1))il(idim)=m%box_ckey_max(idim,ilevel+1)-1
            if(ir(idim)>=m%box_ckey_max(idim,ilevel+1))ir(idim)=m%box_ckey_min(idim,ilevel+1)
         endif
      end do
      ckey = cic_index(il,ir)

      u_cells=0.d0
      kappa_num_cells=0.d0
      skew_cells=0.d0

      hash_nbor(0)=ilevel+1
      do ind=1,twotondim
         hash_nbor(1:ndim)=ckey(1:ndim,ind)
         call get_parent_cell(s,hash_nbor,igrid,icell,flush_cache=.false.,fetch_cache=.true.)
#ifdef HYDRO
#ifdef TRCFLX
         if(igrid>0)then
            rho_cell = max(m%mflux(icell,1,igrid), r%smallr)
            denom = rho_cell
            do idim=1,ndim
               fluxL = m%mflux(icell,1+idim,igrid)
               fluxR = m%mflux(icell,1+idim+ndim,igrid)
               pr = max(fluxR,0.d0)/denom
               pl = max(-fluxL,0.d0)/denom
               cfl_plus = pr+pl
               cfl_minus = pr-pl
               u_cells(idim,ind) = cfl_minus*dx_loc/dt_level
               kappa_num_cells(idim,ind) = 0.5d0*(cfl_plus - cfl_minus**2.d0)*dx_loc**2.d0/dt_level
               skew_cells(idim,ind) = mc_kernel_skewness(pr,pl)
            end do
         end if
#endif
#endif
      end do

      u_eff=0.d0
      kappa_num=0.d0
      skewness_eff=0.d0
      do ind=1,twotondim
         do idim=1,ndim
            u_eff(idim)=u_eff(idim)+u_cells(idim,ind)*vol(ind)
            kappa_num(idim)=kappa_num(idim)+kappa_num_cells(idim,ind)*vol(ind)
            skewness_eff(idim)=skewness_eff(idim)+skew_cells(idim,ind)*vol(ind)
         end do
      end do

      if(trim(r%tracer_kick_pdf)=='gaussian')then
         call sample_tracer_gaussian(xi)
      elseif(trim(r%tracer_kick_pdf)=='piecewise_skew_uniform')then
         call sample_tracer_piecewise_skew_uniform(xi,skewness_eff)
      else
         call sample_tracer_uniform(xi)
      endif

      do idim=1,ndim
         disp(idim)=u_eff(idim)*dt_level
         noise_amp = sqrt(max(0.d0,2.d0*kappa_num(idim)*dt_level))
         disp(idim)=disp(idim)+noise_amp*xi(idim)
      end do

      p%levelp(ipart) = ilevel
      p%vp(ipart,1:ndim)=u_eff(1:ndim)
      p%xp(ipart,1:ndim)=p%xp(ipart,1:ndim)+disp(1:ndim)
      if(s%r%ntrajectories>0)then
         do ii=1,s%r%ntrajectories
            if(s%r%trajectories(ii)==p%idp(ipart))then
               call title(g%myid,nchar)
               filename='trajectory.dat'
               fileloc=TRIM(filename)//TRIM(nchar)
               open(25+g%myid,file=fileloc,status='unknown',access='append')
               write(25+g%myid,*) g%t, p%idp(ipart), &
                    (p%xp(ipart,idim), idim=1,ndim), &
                    (p%vp(ipart,idim), idim=1,ndim)
               close(25+g%myid)
               exit
            endif
         end do
      endif
   end do

   call close_cache(mdl)

   if(action_part==action_kick_drift)then
      do ipart=p%headp(ilevel),p%tailp(ilevel)
         do idim=1,ndim
            if(r%periodic(idim))then
               if(p%xp(ipart,idim)< 0.0d0           )p%xp(ipart,idim)=p%xp(ipart,idim)+r%box_size(idim)
               if(p%xp(ipart,idim)>=r%box_size(idim))p%xp(ipart,idim)=p%xp(ipart,idim)-r%box_size(idim)
            endif
         end do
      end do
   end if

   end associate
 end subroutine cic_trace_gas_part_ito_mc

 subroutine tsc_trace_gas_part_ito_mc(s,p,ilevel,action_part)
   ! Ito MC flux-based tracer with TSC (scheme 5) - same logic as cic_trace_gas_part_ito_mc
   use amr_parameters, only: ndim, threetondim
   use pm_parameters
   use pm_commons, only: part_t
   use ramses_commons, only: ramses_t
   use rng
   use nbors_utils
   use cache_commons
   use cache
   use rho_fine_module, only: tsc_index
   implicit none
   type(ramses_t)::s
   type(part_t)::p
   integer::ilevel
   integer::action_part
   real(kind=8),dimension(1:ndim)::x,disp,xi,u_eff,kappa_num
   integer,dimension(1:ndim)::il,ic,ir
   integer,dimension(1:ndim,1:threetondim)::ckey
   integer(kind=8),dimension(0:ndim)::hash_nbor
   real(kind=8),dimension(1:threetondim)::vol
   real(kind=8),dimension(1:ndim,1:threetondim)::u_cells,kappa_num_cells,skew_cells
   real(kind=8),dimension(1:ndim)::skewness_eff
   real(kind=8),dimension(1:ndim,1:3)::w1d,dw1d
   real(kind=8)::dx_loc,dt_level,rho_cell,denom,fluxL,fluxR,noise_amp,cfl_plus,cfl_minus,pr,pl
   type(msg_hydro_mflux)::dummy_hydro_mflux
   type(RngStream),external::RngStream_CreateStream
   real(kind=8),external::RngStream_RandUni
   integer(kind=8)::stream_skip
   external :: RngStream_SetPackageSeed, RngStream_AdvanceState
   integer :: ipart,idim,ind,icell,igrid,ii
   character(LEN=80) :: filename,fileloc
   character(LEN=5) :: nchar

   associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)
   if(p%static)return
   if (p%type/=TRAC_TYPE) return

   dx_loc=r%boxlen/2**ilevel

   if(action_part==action_kick_only)then
      do ipart=p%headp(ilevel),p%tailp(ilevel)
         p%levelp(ipart)=ilevel
      end do
      return
   endif

   dt_level=g%dtnew(ilevel)

   if(.not. tracer_rng_ready)then
      call RngStream_SetPackageSeed(r%seed)
      tracer_rng = RngStream_CreateStream('tracer_ito_mc')
      stream_skip = int(2*g%myid,kind=8)
      call RngStream_AdvanceState(tracer_rng,0_8,stream_skip)
      tracer_rng_ready = .true.
   end if

   call open_cache(mdl, m, pack_size=storage_size(dummy_hydro_mflux)/32, &
        pack=pack_fetch_kick_trac, unpack=unpack_fetch_kick_trac)

   do ipart=p%headp(ilevel),p%tailp(ilevel)
      do idim=1,ndim
         x(idim)=(p%xp(ipart,idim)+m%skip(idim))/dx_loc
      end do
      call wrap_cell_coords(s,x,ilevel+1)

      call tsc_weights_and_derivs(x, w1d, dw1d, il, ic, ir, vol)
      do idim=1,ndim
         if(r%periodic(idim))then
            if(il(idim)< m%box_ckey_min(idim,ilevel+1))il(idim)=m%box_ckey_max(idim,ilevel+1)-1
            if(ir(idim)>=m%box_ckey_max(idim,ilevel+1))ir(idim)=m%box_ckey_min(idim,ilevel+1)
         endif
      end do
      ckey = tsc_index(il,ic,ir)

      u_cells=0.d0
      kappa_num_cells=0.d0
      skew_cells=0.d0

      hash_nbor(0)=ilevel+1
      do ind=1,threetondim
         hash_nbor(1:ndim)=ckey(1:ndim,ind)
         call get_parent_cell(s,hash_nbor,igrid,icell,flush_cache=.false.,fetch_cache=.true.)
#ifdef HYDRO
#ifdef TRCFLX
         if(igrid>0)then
            rho_cell = max(m%mflux(icell,1,igrid), r%smallr)
            denom = rho_cell
            do idim=1,ndim
               fluxL = m%mflux(icell,1+idim,igrid)
               fluxR = m%mflux(icell,1+idim+ndim,igrid)
               pr = max(fluxR,0.d0)/denom
               pl = max(-fluxL,0.d0)/denom
               cfl_plus = pr+pl
               cfl_minus = pr-pl
               u_cells(idim,ind) = cfl_minus*dx_loc/dt_level
               kappa_num_cells(idim,ind) = 0.5d0*(cfl_plus - cfl_minus**2.d0)*dx_loc**2.d0/dt_level
               skew_cells(idim,ind) = mc_kernel_skewness(pr,pl)
            end do
         end if
#endif
#endif
      end do

      u_eff=0.d0
      kappa_num=0.d0
      skewness_eff=0.d0
      do ind=1,threetondim
         do idim=1,ndim
            u_eff(idim)=u_eff(idim)+u_cells(idim,ind)*vol(ind)
            kappa_num(idim)=kappa_num(idim)+kappa_num_cells(idim,ind)*vol(ind)
            skewness_eff(idim)=skewness_eff(idim)+skew_cells(idim,ind)*vol(ind)
         end do
      end do

      if(trim(r%tracer_kick_pdf)=='gaussian')then
         call sample_tracer_gaussian(xi)
      elseif(trim(r%tracer_kick_pdf)=='piecewise_skew_uniform')then
         call sample_tracer_piecewise_skew_uniform(xi,skewness_eff)
      else
         call sample_tracer_uniform(xi)
      endif

      do idim=1,ndim
         disp(idim)=u_eff(idim)*dt_level
         noise_amp = sqrt(max(0.d0,2.d0*kappa_num(idim)*dt_level))
         disp(idim)=disp(idim)+noise_amp*xi(idim)
      end do

      p%levelp(ipart) = ilevel
      p%vp(ipart,1:ndim)=u_eff(1:ndim)
      p%xp(ipart,1:ndim)=p%xp(ipart,1:ndim)+disp(1:ndim)
      if(s%r%ntrajectories>0)then
         do ii=1,s%r%ntrajectories
            if(s%r%trajectories(ii)==p%idp(ipart))then
               call title(g%myid,nchar)
               filename='trajectory.dat'
               fileloc=TRIM(filename)//TRIM(nchar)
               open(25+g%myid,file=fileloc,status='unknown',access='append')
               write(25+g%myid,*) g%t, p%idp(ipart), &
                    (p%xp(ipart,idim), idim=1,ndim), &
                    (p%vp(ipart,idim), idim=1,ndim)
               close(25+g%myid)
               exit
            endif
         end do
      endif
   end do

   call close_cache(mdl)

   if(action_part==action_kick_drift)then
      do ipart=p%headp(ilevel),p%tailp(ilevel)
         do idim=1,ndim
            if(r%periodic(idim))then
               if(p%xp(ipart,idim)< 0.0d0           )p%xp(ipart,idim)=p%xp(ipart,idim)+r%box_size(idim)
               if(p%xp(ipart,idim)>=r%box_size(idim))p%xp(ipart,idim)=p%xp(ipart,idim)-r%box_size(idim)
            endif
         end do
      end do
   end if

   end associate
 end subroutine tsc_trace_gas_part_ito_mc

 !#########################################################################
 ! Stub routines for Ito MC dust (simplified versions)
 !#########################################################################
 subroutine cic_kick_drift_dust_ito_mc(s,p,ilevel,action_part)
   use ramses_commons, only: ramses_t
   use pm_commons, only: part_t
   implicit none
   type(ramses_t)::s
   type(part_t)::p
   integer::ilevel
   integer::action_part
   ! For now, fall back to regular CIC dust
   call cic_kick_drift_dust(s,p,ilevel,action_part)
 end subroutine cic_kick_drift_dust_ito_mc

 subroutine tsc_kick_drift_dust_ito_mc(s,p,ilevel,action_part)
   use ramses_commons, only: ramses_t
   use pm_commons, only: part_t
   implicit none
   type(ramses_t)::s
   type(part_t)::p
   integer::ilevel
   integer::action_part
   ! For now, fall back to regular TSC dust
   call tsc_kick_drift_dust(s,p,ilevel,action_part)
 end subroutine tsc_kick_drift_dust_ito_mc

 !#########################################################################
 ! Stub routines for guiding center dust (not implemented)
 !#########################################################################
 subroutine cic_kick_drift_dust_guiding_center(s,p,ilevel,action_part)
   use ramses_commons, only: ramses_t
   use pm_commons, only: part_t
   implicit none
   type(ramses_t)::s
   type(part_t)::p
   integer::ilevel
   integer::action_part
   write(*,*)'ERROR: cic_kick_drift_dust_guiding_center not implemented'
   stop
 end subroutine cic_kick_drift_dust_guiding_center

 subroutine tsc_kick_drift_dust_guiding_center(s,p,ilevel,action_part)
   use ramses_commons, only: ramses_t
   use pm_commons, only: part_t
   implicit none
   type(ramses_t)::s
   type(part_t)::p
   integer::ilevel
   integer::action_part
   write(*,*)'ERROR: tsc_kick_drift_dust_guiding_center not implemented'
   stop
 end subroutine tsc_kick_drift_dust_guiding_center

!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cic_kick_drift_dust(s,p,ilevel,action_part)
  use amr_parameters, only: ndim, twotondim, dp
  use hydro_parameters, only: nener
  use pm_parameters
  use pm_commons, only: part_t
  use ramses_commons, only: ramses_t
  use nbors_utils
  use cache_commons
  use cache
  implicit none
  ! Intrinsic functions for drag calculations
  intrinsic :: sinh, cosh
  type(ramses_t)::s
  type(part_t)::p
  integer::ilevel
  integer::action_part
  real(kind=8),dimension(1:ndim)::x,x_mid,dr,dl,dr2,dl2
  integer,dimension(1:ndim)::ir,il
  integer,dimension(1:ndim)::ir2,il2
  real(kind=8),dimension(1:twotondim)::vol,vol2
  integer,dimension(1:ndim,1:twotondim)::ckey,ckey2
  integer::icell,icell2,igrid
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer::ipart,ind,idim,irad
  real(kind=8)::dx_loc,vol_loc
  real(kind=8),dimension(1:ndim)::ff,v_pred
  real(kind=8),dimension(1:ndim)::uu
#ifdef MHD
  real(kind=8),dimension(1:3)::bb
  real(kind=8)::emag
#endif
  real(kind=8)::rho_gas,c_sound,eint,coeff,wdrift2
  real(kind=8)::nu_stop,dens,etot,ekin,erad,cs2,pi=4.0d0*atan(1.0d0)
  real(kind=8),dimension(1:ndim)::what,wdrift ! drift velocity unit vector
  logical :: ok_level
  type(msg_large_realdp)::dummy_large_realdp

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

  if(p%static)return
  if (p%type/=DUST_TYPE) return

  dx_loc=r%boxlen/2**ilevel
  vol_loc=dx_loc**ndim

  coeff=9.0d0*pi*r%gamma/128.0d0
  ! Dust hydro+gravity cache
  call open_cache(mdl, m, pack_size=storage_size(dummy_large_realdp)/32, &
       pack=pack_fetch_kick_dust, unpack=unpack_fetch_kick_dust)
#if NDIM==3
  do ipart=p%headp(ilevel),p%tailp(ilevel)
     ! Position in cell units and wrap
     do idim=1,ndim
        x(idim)=(p%xp(ipart,idim)+m%skip(idim))/dx_loc
     end do
     do idim=1,ndim
        if(r%periodic(idim))then
           if(x(idim)< dble(m%box_ckey_min(idim,ilevel+1)))x(idim)=x(idim)+dble(m%box_ckey_max(idim,ilevel+1)-m%box_ckey_min(idim,ilevel+1))
           if(x(idim)>=dble(m%box_ckey_max(idim,ilevel+1)))x(idim)=x(idim)-dble(m%box_ckey_max(idim,ilevel+1)-m%box_ckey_min(idim,ilevel+1))
        endif
     end do

     ! CIC weights/indices at x
     do idim=1,ndim
        dr(idim)=x(idim)+0.5D0
        ir(idim)=int(dr(idim))
        dr(idim)=dr(idim)-ir(idim)
        dl(idim)=1.0D0-dr(idim)
        il(idim)=ir(idim)-1
     end do
     do idim=1,ndim
        if(r%periodic(idim))then
           if(il(idim)< m%box_ckey_min(idim,ilevel+1))il(idim)=m%box_ckey_max(idim,ilevel+1)-1
           if(ir(idim)>=m%box_ckey_max(idim,ilevel+1))ir(idim)=m%box_ckey_min(idim,ilevel+1)
        endif
     enddo
     ckey = cic_index(il,ir)
     vol = cic_weight(dl,dr)

     if(action_part==action_kick_only)then
        p%levelp(ipart)=ilevel
     else if(action_part==action_kick_drift)then
        ! Leapfrog-style sampling: half-step shift only on the first coarse step
        v_pred(1:ndim)=p%vp(ipart,1:ndim)
        if (g%nstep==0) then
           do idim=1,ndim
              x_mid(idim)=x(idim)+0.5d0*g%dtnew(ilevel)*v_pred(idim)/dx_loc
           end do
        else
           do idim=1,ndim
              x_mid(idim)=x(idim)
           end do
        endif
        do idim=1,ndim
           if(r%periodic(idim))then
              if(x_mid(idim)< dble(m%box_ckey_min(idim,ilevel+1)))x_mid(idim)=x_mid(idim)+dble(m%box_ckey_max(idim,ilevel+1)-m%box_ckey_min(idim,ilevel+1))
              if(x_mid(idim)>=dble(m%box_ckey_max(idim,ilevel+1)))x_mid(idim)=x_mid(idim)-dble(m%box_ckey_max(idim,ilevel+1)-m%box_ckey_min(idim,ilevel+1))
           endif
        end do

        ! CIC weights/indices at x_mid
        do idim=1,ndim
           dr2(idim)=x_mid(idim)+0.5D0
           ir2(idim)=int(dr2(idim))
           dr2(idim)=dr2(idim)-ir2(idim)
           dl2(idim)=1.0D0-dr2(idim)
           il2(idim)=ir2(idim)-1
        end do
        do idim=1,ndim
           if(r%periodic(idim))then
              if(il2(idim)< m%box_ckey_min(idim,ilevel+1))il2(idim)=m%box_ckey_max(idim,ilevel+1)-1
              if(ir2(idim)>=m%box_ckey_max(idim,ilevel+1))ir2(idim)=m%box_ckey_min(idim,ilevel+1)
           endif
        enddo
        ckey2 = cic_index(il2,ir2)
        vol2 = cic_weight(dl2,dr2)

        ! Gather hydro and forces at x_mid
        ff(1:ndim)=0.0
        uu(1:ndim)=0.0
        rho_gas = 0.0
        eint = 0.0
#ifdef MHD
        bb(1:3)=0.0
        emag=0.0d0
#endif
        hash_nbor(0)=ilevel+1
        do ind=1,twotondim
           hash_nbor(1:ndim)=ckey2(1:ndim,ind)
           call get_parent_cell(s,hash_nbor,igrid,icell2,flush_cache=.false.,fetch_cache=.true.)
#ifdef HYDRO
           if(igrid>0)then
#ifdef GRAV
              ff(1:ndim)=ff(1:ndim)+m%f(icell2,1:ndim,igrid)*vol2(ind)
#endif
              rho_gas = rho_gas + m%uold(icell2,1,igrid)*vol2(ind)
              uu(1:ndim)=uu(1:ndim)+real(m%uold(icell2,2:ndim+1,igrid),kind=8)/max(real(m%uold(icell2,1,igrid),kind=8),real(r%smallr,kind=8))*vol2(ind)
#ifdef MHD
              bb(1:3)=bb(1:3)+0.5d0*(m%bold(icell2,1:3,igrid)+m%bold(icell2,4:6,igrid))*vol2(ind)
#endif
              dens = max(dble(m%uold(icell2,1,igrid)), r%smallr)
              etot = m%uold(icell2,5,igrid)
              ekin = 0.0d0
              do idim=1,ndim
                 ekin = ekin + 0.5d0 * m%uold(icell2,1+idim,igrid)**2 / dens
              end do
              erad = 0.0d0
#if NENER>0
              do irad=1,nener
                 erad = erad + m%uold(icell2,5+irad,igrid)
              end do
#endif
#ifdef MHD
              emag=0.0d0
              do idim=1,3
                 emag = emag + 0.125d0*(m%bold(icell2,idim,igrid)+m%bold(icell2,3+idim,igrid))**2
              end do
              eint = eint - emag*vol2(ind)
#endif
              eint = eint + (etot - ekin - erad)*vol2(ind)
           end if
#endif
        end do

        cs2 = r%gamma * (r%gamma-1.0d0) * max(eint, r%smallc**2) / max(rho_gas, r%smallr)
        c_sound = max(sqrt(cs2), r%smallc)
        nu_stop = coeff*c_sound*rho_gas/p%size(ipart)
        wdrift(1:ndim)= v_pred(1:ndim) - uu(1:ndim)

        ! Operator-split: drag(half) -> forces -> drag(half)
        call compute_drag_step(wdrift, c_sound, 0.5d0*g%dtnew(ilevel), nu_stop, coeff, r%analytic_dust_force)
#ifdef GRAV
        wdrift(1:ndim)=wdrift(1:ndim)+ff(1:ndim)*0.5d0*g%dtnew(ilevel)
#endif
#ifdef MHD
        call compute_lorentz_step(wdrift, bb(1:3), g%dtnew(ilevel), p%charge(ipart), r%analytic_dust_force)
#endif
#ifdef GRAV
        wdrift(1:ndim)=wdrift(1:ndim)+ff(1:ndim)*0.5d0*g%dtnew(ilevel)
#endif
        call compute_drag_step(wdrift, c_sound, 0.5d0*g%dtnew(ilevel), nu_stop, coeff, r%analytic_dust_force)

        p%vp(ipart,1:ndim)=uu(1:ndim)+wdrift(1:ndim)
        p%xp(ipart,1:ndim)=p%xp(ipart,1:ndim)+p%vp(ipart,1:ndim)*g%dtnew(ilevel)
     endif
  end do
#endif
  call close_cache(mdl)

  ! Periodic boundary conditions
  if(action_part==action_kick_drift)then
     do ipart=p%headp(ilevel),p%tailp(ilevel)
        do idim=1,ndim
           if(r%periodic(idim))then
              if(p%xp(ipart,idim)< 0.0d0           )p%xp(ipart,idim)=p%xp(ipart,idim)+r%box_size(idim)
              if(p%xp(ipart,idim)>=r%box_size(idim))p%xp(ipart,idim)=p%xp(ipart,idim)-r%box_size(idim)
           end if
        end do
     end do
  end if
  end associate
end subroutine cic_kick_drift_dust

subroutine tsc_kick_drift_dust(s,p,ilevel,action_part)
  use amr_parameters, only: ndim, threetondim, dp
  use hydro_parameters, only: nener
  use pm_parameters
  use pm_commons, only: part_t
  use ramses_commons, only: ramses_t
  use nbors_utils
  use cache_commons
  use cache
  implicit none
  type(ramses_t)::s
  type(part_t)::p
  integer::ilevel
  integer::action_part
  real(kind=8),dimension(1:ndim)::x,wl,wc,wr
  real(kind=8),dimension(1:ndim)::x_mid,wl2,wc2,wr2
  integer,dimension(1:ndim)::cl,cc,cr
  integer,dimension(1:ndim)::cl2,cc2,cr2
  real(kind=8),dimension(1:threetondim)::vol,vol2
  integer,dimension(1:ndim,1:threetondim)::ckey,ckey2
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer::ipart,icell,icell2,igrid,ind,idim,irad
  real(kind=8)::xl,xc,xr
  real(kind=8)::dx_loc
  real(kind=8),dimension(1:ndim)::ff,uu,v_pred,wdrift
  real(kind=8),dimension(1:3)::bb
  real(kind=8)::rho_gas,c_sound,eint,coeff
  real(kind=8)::nu_stop,dens,etot,ekin,erad,emag,cs2,pi
  integer :: ii
  character(LEN=80)::filename,fileloc
  character(LEN=5)::nchar
  type(msg_large_realdp)::dummy_large_realdp

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

  if(p%static)return
  if (p%type/=DUST_TYPE) return

  dx_loc=r%boxlen/2**ilevel

  pi=4.0d0*atan(1.0d0)
  coeff=9.0d0*pi*r%gamma/128.0d0

  call open_cache(mdl, m, pack_size=storage_size(dummy_large_realdp)/32, &
       pack=pack_fetch_kick_dust, unpack=unpack_fetch_kick_dust)
#if NDIM==3
  do ipart=p%headp(ilevel),p%tailp(ilevel)
     ! particle position in cell units and periodic wrap
     do idim=1,ndim
        x(idim)=(p%xp(ipart,idim)+m%skip(idim))/dx_loc
     end do
     do idim=1,ndim
        if(r%periodic(idim))then
           if(x(idim)< dble(m%box_ckey_min(idim,ilevel+1)))x(idim)=x(idim)+dble(m%box_ckey_max(idim,ilevel+1)-m%box_ckey_min(idim,ilevel+1))
           if(x(idim)>=dble(m%box_ckey_max(idim,ilevel+1)))x(idim)=x(idim)-dble(m%box_ckey_max(idim,ilevel+1)-m%box_ckey_min(idim,ilevel+1))
        endif
     end do

     if(action_part==action_kick_only)then
        p%levelp(ipart)=ilevel

     else if(action_part==action_kick_drift)then
        ! Leapfrog-style sampling: do a half-step advance only on the first step,
        ! then sample at the current position thereafter (we are already time-centered).
        v_pred(1:ndim)=p%vp(ipart,1:ndim)

        if (g%nstep==0) then
           ! First step: predict position at midpoint t + dt/2 for sampling
           do idim=1,ndim
              x_mid(idim)=x(idim)+0.5d0*g%dtnew(ilevel)*v_pred(idim)/dx_loc
           end do
        else
           ! Subsequent steps: positions are already at mid-time, no shift needed
           do idim=1,ndim
              x_mid(idim)=x(idim)
           end do
        endif
        do idim=1,ndim
           if(r%periodic(idim))then
              if(x_mid(idim)< dble(m%box_ckey_min(idim,ilevel+1)))x_mid(idim)=x_mid(idim)+dble(m%box_ckey_max(idim,ilevel+1)-m%box_ckey_min(idim,ilevel+1))
              if(x_mid(idim)>=dble(m%box_ckey_max(idim,ilevel+1)))x_mid(idim)=x_mid(idim)-dble(m%box_ckey_max(idim,ilevel+1)-m%box_ckey_min(idim,ilevel+1))
           endif
        end do
        do idim=1,ndim
           cl2(idim)=int(x_mid(idim))-1
           cc2(idim)=int(x_mid(idim))
           cr2(idim)=int(x_mid(idim))+1
           xl=dble(cl2(idim))+0.5D0
           xc=dble(cc2(idim))+0.5D0
           xr=dble(cr2(idim))+0.5D0
           wl2(idim)=0.5D0*(1.5D0-abs(x_mid(idim)-xl))**2
           wc2(idim)=0.75D0-         (x_mid(idim)-xc) **2
           wr2(idim)=0.5D0*(1.5D0-abs(x_mid(idim)-xr))**2
        end do
        do idim=1,ndim
           if(r%periodic(idim))then
              if(cl2(idim)< m%box_ckey_min(idim,ilevel+1))cl2(idim)=m%box_ckey_max(idim,ilevel+1)-1
              if(cr2(idim)>=m%box_ckey_max(idim,ilevel+1))cr2(idim)=m%box_ckey_min(idim,ilevel+1)
           endif
        enddo
        ckey2 = tsc_index(cl2,cc2,cr2)
        vol2 = tsc_weight(wl2,wc2,wr2)

        ff(1:ndim)=0.0d0
        uu(1:3)=0.0d0
        bb(1:3)=0.0d0
        rho_gas=0.0d0
        eint=0.0d0
        emag=0.0d0
        hash_nbor(0)=ilevel+1
        do ind=1,threetondim
           hash_nbor(1:ndim)=ckey2(1:ndim,ind)
           call get_parent_cell(s,hash_nbor,igrid,icell2,flush_cache=.false.,fetch_cache=.true.)
#ifdef HYDRO
           if(igrid>0)then
#ifdef GRAV
              ff(1:ndim)=ff(1:ndim)+m%f(icell2,1:ndim,igrid)*vol2(ind)
#endif
              rho_gas = rho_gas + m%uold(icell2,1,igrid)*vol2(ind)
              uu(1:ndim)=uu(1:ndim)+real(m%uold(icell2,2:ndim+1,igrid),kind=8)/max(real(m%uold(icell2,1,igrid),kind=8),real(r%smallr,kind=8))*vol2(ind)
#ifdef MHD
              bb(1:3)=bb(1:3)+0.5d0*(m%bold(icell2,1:3,igrid)+m%bold(icell2,4:6,igrid))*vol2(ind)
#endif
              dens = max(dble(m%uold(icell2,1,igrid)), r%smallr)
              etot = m%uold(icell2,5,igrid)
              ekin = 0.0d0
              do idim=1,ndim
                 ekin = ekin + 0.5d0 * m%uold(icell2,1+idim,igrid)**2 / dens
              end do
              erad = 0.0d0
#if NENER>0
              do irad=1,nener
                 erad = erad + m%uold(icell2,5+irad,igrid)
              end do
#endif
#ifdef MHD
              emag=0.0d0
              do idim=1,3
                 emag = emag + 0.125d0*(m%bold(icell2,idim,igrid)+m%bold(icell2,3+idim,igrid))**2
              end do
#endif
              eint = eint + (etot - ekin - erad - emag) * vol2(ind)
           end if

#endif
        end do
        cs2 = r%gamma * (r%gamma-1.0d0) * max(eint, r%smallc**2) / max(rho_gas, r%smallr)
        c_sound = max(sqrt(cs2), r%smallc)
        nu_stop = coeff*c_sound*rho_gas/p%size(ipart)
        wdrift(1:ndim)= v_pred(1:ndim) - uu(1:ndim)
        ! This is where we want to split off the different physics
        ! There will also be an if(gyro_pic) gate, since everything must be computed very differently
        ! in the gyro case.

        call compute_drag_step(wdrift, c_sound, 0.5*g%dtnew(ilevel), nu_stop, coeff, r%analytic_dust_force)
        ! Gather ff, and apply to either side of the Lorentz force as a half-step
        wdrift(1:ndim)=wdrift(1:ndim)+ff(1:ndim)*0.5d0*g%dtnew(ilevel) ! External force half-step (includes gravity)
#ifdef MHD
        call compute_lorentz_step(wdrift, bb(1:3), g%dtnew(ilevel), p%charge(ipart), r%analytic_dust_force) !Somehow introducing nonphysical oscillations.
#endif
        wdrift(1:ndim)=wdrift(1:ndim)+ff(1:ndim)*0.5d0*g%dtnew(ilevel) ! Second external force half-step
        call compute_drag_step(wdrift, c_sound, 0.5*g%dtnew(ilevel), nu_stop, coeff, r%analytic_dust_force)
        ! Routine will return an intermediate drift velocity.
        p%vp(ipart,1:ndim)=uu(1:ndim)+wdrift(1:ndim)
        ! Leapfrog drift: advance positions using time-centered velocity
        p%xp(ipart,1:ndim)=p%xp(ipart,1:ndim)+p%vp(ipart,1:ndim)*g%dtnew(ilevel)
        ! Trajectory output for selected particles
        if(s%r%ntrajectories>0)then
           do ii=1,s%r%ntrajectories
              if(s%r%trajectories(ii)==p%idp(ipart))then
                 call title(g%myid,nchar)
                 filename='trajectory.dat'
                 fileloc=TRIM(filename)//TRIM(nchar)
                 open(25+g%myid,file=fileloc,status='unknown',access='append')
                 write(25+g%myid,'(1PE15.7,1X,I12,1X,6(1PE15.7,1X),1PE15.7,1X,1PE15.7)') &
                      g%t, p%idp(ipart), &
                      p%xp(ipart,1),p%xp(ipart,2),p%xp(ipart,3), &
                      p%vp(ipart,1),p%vp(ipart,2),p%vp(ipart,3), &
                      p%charge(ipart), p%size(ipart)
                 close(25+g%myid)
                 exit
              endif
           end do
        endif
     endif
  end do
#endif
  call close_cache(mdl)
  if(action_part==action_kick_drift)then
     do ipart=p%headp(ilevel),p%tailp(ilevel)
        do idim=1,ndim
           if(r%periodic(idim))then
              if(p%xp(ipart,idim)< 0.0d0           )p%xp(ipart,idim)=p%xp(ipart,idim)+r%box_size(idim)
              if(p%xp(ipart,idim)>=r%box_size(idim))p%xp(ipart,idim)=p%xp(ipart,idim)-r%box_size(idim)
           end if
        end do
     end do
  end if
  end associate
end subroutine tsc_kick_drift_dust

subroutine pcs_kick_drift_dust(s,p,ilevel,action_part)
  use amr_parameters, only: ndim, fourtondim, dp
  use hydro_parameters, only: nener
  use pm_parameters
  use pm_commons, only: part_t
  use ramses_commons, only: ramses_t
  use nbors_utils
  use cache_commons
  use cache
  implicit none
  type(ramses_t)::s
  type(part_t)::p
  integer::ilevel
  integer::action_part
  real(kind=8),dimension(1:ndim)::x,wll,wl,wr,wrr
  real(kind=8),dimension(1:ndim)::x_mid,wll2,wl2,wr2,wrr2
  integer,dimension(1:ndim)::cll,cl,cr,crr
  integer,dimension(1:ndim)::cll2,cl2,cr2,crr2
  real(kind=8),dimension(1:fourtondim)::vol,vol2
  integer,dimension(1:ndim,1:fourtondim)::ckey,ckey2
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer::ipart,icell,icell2,igrid,ind,idim
  real(kind=8)::xll,xl,xr,xrr
  real(kind=8)::dx_loc
  real(kind=8),dimension(1:ndim)::ff,v_pred,uu,what,wdrift
#ifdef MHD
  real(kind=8),dimension(1:3)::bb
  real(kind=8)::emag
#endif
  real(kind=8)::rho_gas,c_sound,eint,coeff,wdrift2
  real(kind=8)::nu_stop,dens,etot,ekin,erad,cs2,pi
  integer::irad
  type(msg_large_realdp)::dummy_large_realdp

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

  if(p%static)return
  if (p%type/=DUST_TYPE) return

  dx_loc=r%boxlen/2**ilevel

  pi=4.0d0*atan(1.0d0)
  coeff=9.0d0*pi*r%gamma/128.0d0

  call open_cache(mdl, m, pack_size=storage_size(dummy_large_realdp)/32, &
       pack=pack_fetch_kick_dust, unpack=unpack_fetch_kick_dust)
#if NDIM==3
  do ipart=p%headp(ilevel),p%tailp(ilevel)
     do idim=1,ndim
        x(idim)=(p%xp(ipart,idim)+m%skip(idim))/dx_loc
     end do
     do idim=1,ndim
        if(r%periodic(idim))then
           if(x(idim)< dble(m%box_ckey_min(idim,ilevel+1)))x(idim)=x(idim)+dble(m%box_ckey_max(idim,ilevel+1)-m%box_ckey_min(idim,ilevel+1))
           if(x(idim)>=dble(m%box_ckey_max(idim,ilevel+1)))x(idim)=x(idim)-dble(m%box_ckey_max(idim,ilevel+1)-m%box_ckey_min(idim,ilevel+1))
        endif
     end do
     do idim=1,ndim
        crr(idim)=int(x(idim)+1.5D0)
        cr (idim)=crr(idim)-1
        cl (idim)=crr(idim)-2
        cll(idim)=crr(idim)-3
        xll=dble(cll(idim))+0.5D0
        xl =dble(cl (idim))+0.5D0
        xr =dble(cr (idim))+0.5D0
        xrr=dble(crr(idim))+0.5D0
        wll(idim)=(2D0-abs(x(idim)-xll))**3/6D0
        wl (idim)=(4D0-6D0*(x(idim)-xl)**2+3d0*abs(x(idim)-xl )**3)/6D0
        wr (idim)=(4D0-6D0*(x(idim)-xr)**2+3d0*abs(x(idim)-xr )**3)/6D0
        wrr(idim)=(2D0-abs(x(idim)-xrr))**3/6D0
     end do
     do idim=1,ndim
        if(r%periodic(idim))then
           if(cll(idim)< m%box_ckey_min(idim,ilevel+1))cll(idim)=cll(idim)-m%box_ckey_min(idim,ilevel+1)+m%box_ckey_max(idim,ilevel+1)
           if(cl (idim)< m%box_ckey_min(idim,ilevel+1))cl (idim)=cl (idim)-m%box_ckey_min(idim,ilevel+1)+m%box_ckey_max(idim,ilevel+1)
           if(cr (idim)>=m%box_ckey_max(idim,ilevel+1))cr (idim)=cr (idim)+m%box_ckey_min(idim,ilevel+1)-m%box_ckey_max(idim,ilevel+1)
           if(crr(idim)>=m%box_ckey_max(idim,ilevel+1))crr(idim)=crr(idim)+m%box_ckey_min(idim,ilevel+1)-m%box_ckey_max(idim,ilevel+1)
        endif
     enddo
     ckey = pcs_index(cll,cl,cr,crr)
     vol = pcs_weight(wll,wl,wr,wrr)
     hash_nbor(0)=ilevel+1
     ff(1:ndim)=0.0
     do ind=1,fourtondim
        hash_nbor(1:ndim)=ckey(1:ndim,ind)
        call get_parent_cell(s,hash_nbor,igrid,icell,flush_cache=.false.,fetch_cache=.true.)
#ifdef HYDRO
        if(igrid>0)then
           ff(1:ndim)=ff(1:ndim)+real(m%uold(icell,2:ndim+1,igrid),kind=8)/max(real(m%uold(icell,1,igrid),kind=8),real(r%smallr,kind=8))*vol(ind)
        end if
#endif
     end do

    if(action_part==action_kick_only)then
       ! RK2 step 1: stash state; keep velocity as-is for drag scheme
        p%levelp(ipart)=ilevel

    else if(action_part==action_kick_drift)then
       ! RK2 step 2: predict with v^n, sample gas state at midpoint
       if (g%nstep>0) then
          v_pred(1:ndim)=p%vp(ipart,1:ndim)
       else
          v_pred(1:ndim)=ff(1:ndim)
       endif
       do idim=1,ndim
          x_mid(idim)=x(idim)+0.5d0*g%dtnew(ilevel)*v_pred(idim)/dx_loc
       end do
       do idim=1,ndim
          if(r%periodic(idim))then
             if(x_mid(idim)< dble(m%box_ckey_min(idim,ilevel+1)))x_mid(idim)=x_mid(idim)+dble(m%box_ckey_max(idim,ilevel+1)-m%box_ckey_min(idim,ilevel+1))
             if(x_mid(idim)>=dble(m%box_ckey_max(idim,ilevel+1)))x_mid(idim)=x_mid(idim)-dble(m%box_ckey_max(idim,ilevel+1)-m%box_ckey_min(idim,ilevel+1))
          endif
       end do
       do idim=1,ndim
          crr2(idim)=int(x_mid(idim)+1.5D0)
          cr2 (idim)=crr2(idim)-1
          cl2 (idim)=crr2(idim)-2
          cll2(idim)=crr2(idim)-3
          xll=dble(cll2(idim))+0.5D0
          xl =dble(cl2 (idim))+0.5D0
          xr =dble(cr2 (idim))+0.5D0
          xrr=dble(crr2(idim))+0.5D0
          wll2(idim)=(2D0-abs(x_mid(idim)-xll))**3/6D0
          wl2 (idim)=(4D0-6D0*(x_mid(idim)-xl)**2+3d0*abs(x_mid(idim)-xl )**3)/6D0
          wr2 (idim)=(4D0-6D0*(x_mid(idim)-xr)**2+3d0*abs(x_mid(idim)-xr )**3)/6D0
          wrr2(idim)=(2D0-abs(x_mid(idim)-xrr))**3/6D0
       end do
       do idim=1,ndim
          if(r%periodic(idim))then
             if(cll2(idim)< m%box_ckey_min(idim,ilevel+1))cll2(idim)=cll2(idim)-m%box_ckey_min(idim,ilevel+1)+m%box_ckey_max(idim,ilevel+1)
             if(cl2 (idim)< m%box_ckey_min(idim,ilevel+1))cl2 (idim)=cl2 (idim)-m%box_ckey_min(idim,ilevel+1)+m%box_ckey_max(idim,ilevel+1)
             if(cr2 (idim)>=m%box_ckey_max(idim,ilevel+1))cr2 (idim)=cr2 (idim)+m%box_ckey_min(idim,ilevel+1)-m%box_ckey_max(idim,ilevel+1)
             if(crr2(idim)>=m%box_ckey_max(idim,ilevel+1))crr2(idim)=crr2(idim)+m%box_ckey_min(idim,ilevel+1)-m%box_ckey_max(idim,ilevel+1)
          endif
       enddo
       ckey2 = pcs_index(cll2,cl2,cr2,crr2)
       vol2 = pcs_weight(wll2,wl2,wr2,wrr2)
       hash_nbor(0)=ilevel+1
       ff(1:ndim)=0.0
       uu(1:ndim)=0.0
       rho_gas=0.0d0
       eint=0.0d0
#ifdef MHD
       bb(1:3)=0.0d0
       emag=0.0d0
#endif
       do ind=1,fourtondim
          hash_nbor(1:ndim)=ckey2(1:ndim,ind)
          call get_parent_cell(s,hash_nbor,igrid,icell2,flush_cache=.false.,fetch_cache=.true.)
#ifdef HYDRO
          if(igrid>0)then
#ifdef GRAV
             ff(1:ndim)=ff(1:ndim)+m%f(icell2,1:ndim,igrid)*vol2(ind)
#endif
             rho_gas = rho_gas + m%uold(icell2,1,igrid)*vol2(ind)
             uu(1:ndim)=uu(1:ndim)+real(m%uold(icell2,2:ndim+1,igrid),kind=8)/max(real(m%uold(icell2,1,igrid),kind=8),real(r%smallr,kind=8))*vol2(ind)
#ifdef MHD
             bb(1:3)=bb(1:3)+0.5d0*(m%bold(icell2,1:3,igrid)+m%bold(icell2,4:6,igrid))*vol2(ind)
#endif
             dens = max(dble(m%uold(icell2,1,igrid)),r%smallr)
             etot = m%uold(icell2,5,igrid)
             ekin = 0.0d0
             do idim=1,ndim
                ekin = ekin + 0.5d0 * m%uold(icell2,1+idim,igrid)**2 / dens
             end do
             erad = 0.0d0
#if NENER>0
             do irad=1,nener
                erad = erad + m%uold(icell2,5+irad,igrid)
             end do
#endif
#ifdef MHD
             emag=0.0d0
             do idim=1,3
                emag = emag + 0.125d0*(m%bold(icell2,idim,igrid)+m%bold(icell2,3+idim,igrid))**2
             end do
             eint = eint - emag*vol2(ind)
#endif
             eint = eint + (etot - ekin - erad) * vol2(ind)
          end if
#endif
       end do
       cs2 = r%gamma * (r%gamma-1.0d0) * max(eint, r%smallc**2) / max(rho_gas, r%smallr)
       c_sound = max(sqrt(cs2), r%smallc)
       nu_stop = coeff*c_sound*rho_gas/p%size(ipart)
       wdrift(1:ndim)= v_pred(1:ndim) - uu(1:ndim)
       ! Operator-split: drag(half) -> forces -> drag(half)
       call compute_drag_step(wdrift, c_sound, 0.5d0*g%dtnew(ilevel), nu_stop, coeff, r%analytic_dust_force)
#ifdef GRAV
       wdrift(1:ndim)=wdrift(1:ndim)+ff(1:ndim)*0.5d0*g%dtnew(ilevel)
#endif
#ifdef MHD
       call compute_lorentz_step(wdrift, bb(1:3), g%dtnew(ilevel), p%charge(ipart), r%analytic_dust_force)
#endif
#ifdef GRAV
       wdrift(1:ndim)=wdrift(1:ndim)+ff(1:ndim)*0.5d0*g%dtnew(ilevel)
#endif
       call compute_drag_step(wdrift, c_sound, 0.5d0*g%dtnew(ilevel), nu_stop, coeff, r%analytic_dust_force)
       p%vp(ipart,1:ndim)=uu(1:ndim)+wdrift(1:ndim)
       p%xp(ipart,1:ndim)=p%xp(ipart,1:ndim)+p%vp(ipart,1:ndim)*g%dtnew(ilevel)
    endif
  end do
#endif
  call close_cache(mdl)
  if(action_part==action_kick_drift)then
     do ipart=p%headp(ilevel),p%tailp(ilevel)
        do idim=1,ndim
           if(r%periodic(idim))then
              if(p%xp(ipart,idim)< 0.0d0           )p%xp(ipart,idim)=p%xp(ipart,idim)+r%box_size(idim)
              if(p%xp(ipart,idim)>=r%box_size(idim))p%xp(ipart,idim)=p%xp(ipart,idim)-r%box_size(idim)
           end if
        end do
     end do
  end if
  end associate
end subroutine pcs_kick_drift_dust
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine compute_drag_step(wdrift, c_sound, dt, nu_stop, coeff, analytic_dust_force)
  use amr_parameters, only: ndim
  implicit none
  real(kind=8), intent(inout) :: wdrift(1:ndim)
  real(kind=8), intent(in)    :: c_sound, dt, nu_stop, coeff
  logical, intent(in)         :: analytic_dust_force
  real(kind=8) :: wdrift2
  real(kind=8) :: what(1:ndim)

  wdrift2 = dot_product(wdrift(1:ndim), wdrift(1:ndim))
  if (wdrift2 > 0.0d0) then
     what(1:ndim) = wdrift(1:ndim) / sqrt(wdrift2)
     if (analytic_dust_force) then
        wdrift(1:ndim) = c_sound * what(1:ndim) * &
        & exact_drag(dt, nu_stop, coeff, sqrt(wdrift2)/c_sound)
     else
        wdrift(1:ndim) = c_sound * what(1:ndim) * &
        & fully_implicit_drag(dt, nu_stop, coeff, sqrt(wdrift2)/c_sound)
     end if
  else
     wdrift(1:ndim) = 0.0d0
  end if

end subroutine compute_drag_step

function exact_drag(dt, nu, eta, w0) result(drag_value)
  ! Computes drag coefficient based on the analytic solution for Epstein-Baines drag:
  ! drag(dt) = sqrt(((sinh(nu*dt) + sqrt(1 + eta*w0^2)*cosh(nu*dt)) /
  !                  (cosh(nu*dt) + sqrt(1 + nu*w0^2)*sinh(nu*dt)))^2 - 1)
  implicit none
  real(kind=8), intent(in) :: dt, nu, eta, w0
  real(kind=8) :: drag_value

  real(kind=8) :: sinh_term, cosh_term, sqrt_term, numerator, denominator, ratio

  ! Compute hyperbolic functions
  sinh_term = sinh(nu * dt)
  cosh_term = cosh(nu * dt)

  ! Compute sqrt(1 + eta*w0^2) - note: using eta as in the original formula
  sqrt_term = sqrt(1.0d0 + eta * w0**2)

  ! Compute numerator: sinh(nu*dt) + sqrt(1 + eta*w0^2)*cosh(nu*dt)
  numerator = sinh_term + sqrt_term * cosh_term

  ! Compute denominator: cosh(nu*t) + sqrt(1 + eta*w0^2)*sinh(nu*t)
  denominator = cosh_term + sqrt(1.0d0 + eta * w0**2) * sinh_term

  ! Avoid division by zero
  if (abs(denominator) < 1.0d-15) then
     drag_value = 0.0d0
     return
  end if

  ! Compute ratio and final result
  ratio = numerator / denominator
  drag_value = sqrt(max(0.0d0, ratio**2 - 1.0d0))/sqrt(eta)

end function exact_drag

!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
function fully_implicit_drag(dt, nu, eta, w0) result(drag_value)
  ! A second order fully implicit drag update
  implicit none
  real(kind=8), intent(in) :: dt, nu, eta, w0
  real(kind=8) :: drag_value, nu_stop

  real(kind=8) :: sqrt_term, ratio


  ! Compute sqrt(1 + eta*w0^2) w0^2 is actually w0^2/c_sound^2
  sqrt_term = sqrt(1.0d0 + eta * w0**2)
  nu_stop = nu * sqrt_term

  drag_value = w0/(1+nu_stop*dt+0.5d0*nu*dt**2)
  ! Not a typo. Comes from a series expansion of the analytic solution in dt


end function fully_implicit_drag

!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine compute_lorentz(driftvel, bfield, dt, charge_parameter)
  use amr_parameters, only: ndim, dp
  implicit none
  real(kind=8), dimension(1:ndim), intent(inout) :: driftvel
  real(kind=8), dimension(1:ndim), intent(in)    :: bfield
  real(kind=8), intent(in)                       :: dt
  real(dp), intent(in)                           :: charge_parameter
  real(kind=8) :: det, bsquared, dteff
  real(kind=8) :: v1_new, v2_new, v3_new
  real(kind=8), dimension(1:3,1:3) :: matrix

#if NDIM==3
  dteff = -1.0d0 * dt * charge_parameter ! Accidentally flipped sign in the original code.
  bsquared = dot_product(bfield(1:ndim), bfield(1:ndim))

  det = 1 + 0.25d0 * bsquared * dteff**2

  ! Matrix components from Mathematica
  ! Inverse[(e - 0.5*dteff*m)].(e+0.5*dteff*m).v^n = v^(n+1)
  ! Row 1
  matrix(1,1) = 1.0d0 + 0.25d0 * (bfield(1)**2 - bfield(2)**2 - bfield(3)**2) * dteff**2
  matrix(1,2) = 0.5d0 * dteff * (2.0d0 * bfield(3) + bfield(1) * bfield(2) * dteff)
  matrix(1,3) = 0.5d0 * dteff * (-2.0d0 * bfield(2) + bfield(1) * bfield(3) * dteff)

  ! Row 2
  matrix(2,1) = 0.5d0 * dteff * (-2.0d0 * bfield(3) + bfield(1) * bfield(2) * dteff)
  matrix(2,2) = 1.0d0 - 0.25d0 * (bfield(1)**2 - bfield(2)**2 + bfield(3)**2) * dteff**2
  matrix(2,3) = 0.5d0 * dteff * (2.0d0 * bfield(1) + bfield(2) * bfield(3) * dteff)

  ! Row 3
  matrix(3,1) = 0.5d0 * dteff * (2.0d0 * bfield(2) + bfield(1) * bfield(3) * dteff)
  matrix(3,2) = 0.5d0 * dteff * (-2.0d0 * bfield(1) + bfield(2) * bfield(3) * dteff)
  matrix(3,3) = 1.0d0 - 0.25d0 * (bfield(1)**2 + bfield(2)**2 - bfield(3)**2) * dteff**2

  ! Apply matrix/det to driftvel (hard-coded for efficiency)
  if (ndim == 3) then
     v1_new = (matrix(1,1)*driftvel(1) + matrix(1,2)*driftvel(2) + matrix(1,3)*driftvel(3)) / det
     v2_new = (matrix(2,1)*driftvel(1) + matrix(2,2)*driftvel(2) + matrix(2,3)*driftvel(3)) / det
     v3_new = (matrix(3,1)*driftvel(1) + matrix(3,2)*driftvel(2) + matrix(3,3)*driftvel(3)) / det
     driftvel(1) = v1_new
     driftvel(2) = v2_new
     driftvel(3) = v3_new
  end if
#endif

end subroutine compute_lorentz

!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine compute_lorentz_step(driftvel, bfield, dt, charge_parameter, analytic_dust_force)
  use amr_parameters, only: ndim, dp
  implicit none
  real(kind=8), dimension(1:ndim), intent(inout) :: driftvel
  real(kind=8), dimension(1:ndim), intent(in)    :: bfield
  real(kind=8), intent(in)                       :: dt
  real(dp), intent(in)                           :: charge_parameter
  logical, intent(in)                            :: analytic_dust_force
#ifdef MHD
  if (.not. analytic_dust_force) then
     call compute_lorentz(driftvel, bfield, dt, charge_parameter)
  else
     call compute_lorentz_analytic(driftvel, bfield, dt, charge_parameter)
  end if
#endif
end subroutine compute_lorentz_step
!#########################################################################
!#########################################################################
subroutine compute_lorentz_analytic(driftvel, bfield, dt, charge_parameter)
  use amr_parameters, only: ndim, dp
  implicit none
  real(kind=8), dimension(1:ndim), intent(inout) :: driftvel
  real(kind=8), dimension(1:ndim), intent(in)    :: bfield
  real(kind=8), intent(in)                       :: dt
  real(dp), intent(in)                           :: charge_parameter
  real(kind=8) :: dteff, bsquared, bnorm, theta, costh, sinth, vdotb, t2
  real(kind=8), dimension(1:ndim) :: bhat, v, kxv
#if NDIM==3
  dteff = dt * charge_parameter
  if (dteff == 0.0d0) return

  bsquared = dot_product(bfield(1:ndim), bfield(1:ndim))
  if (bsquared <= 0.0d0) return

  bnorm = sqrt(bsquared)
  bhat(1:ndim) = bfield(1:ndim) / bnorm
  theta = dteff * bnorm
  ! It's worth noting that our gyro_factor limits theta to be less than 0.63
  t2 = theta*theta
  if (abs(theta) < 1.0d-1) then
     ! Horner-form small-angle approximations for better performance/accuracy
     costh = 1.0d0 + t2*(-0.5d0 + t2*( 1.0d0/24.0d0 + t2*( -1.0d0/720.0d0 )))
     sinth = theta * (1.0d0 + t2*( -1.0d0/6.0d0 + t2*( 1.0d0/120.0d0 + t2*( -1.0d0/5040.0d0 ))))
  else
     costh = cos(theta)
     sinth = sin(theta)
  end if

  v(1:ndim) = driftvel(1:ndim)

  if (ndim == 3) then
     ! Rodrigues' rotation formula: rotate v around unit axis bhat by angle theta
     kxv(1) = bhat(2)*v(3) - bhat(3)*v(2)
     kxv(2) = bhat(3)*v(1) - bhat(1)*v(3)
     kxv(3) = bhat(1)*v(2) - bhat(2)*v(1)
     vdotb  = bhat(1)*v(1) + bhat(2)*v(2) + bhat(3)*v(3)

     driftvel(1) = v(1)*costh + kxv(1)*sinth + bhat(1)*vdotb*(1.0d0 - costh)
     driftvel(2) = v(2)*costh + kxv(2)*sinth + bhat(2)*vdotb*(1.0d0 - costh)
     driftvel(3) = v(3)*costh + kxv(3)*sinth + bhat(3)*vdotb*(1.0d0 - costh)
  end if
#endif
end subroutine compute_lorentz_analytic

!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################

! Weight computation subroutines for CIC and TSC interpolation
!#########################################################################
subroutine cic_weights_and_derivs(x, w1d, dw1d, il_out, ir_out, vol_out)
  !--------------------------------------------------------------
  ! Compute 1D CIC weights, derivatives, and cell indices from position.
  ! x(1:ndim) is the cell-centered coordinate.
  ! w1d(dim,1) = left weight (dl), w1d(dim,2) = right weight (dr)
  ! dw1d(dim,1) = -1, dw1d(dim,2) = +1 (CIC derivatives are constant)
  ! il_out, ir_out = cell indices (optional)
  ! vol_out(1:twotondim) = 3D volume weights (optional)
  !--------------------------------------------------------------
  use amr_parameters, only: ndim, twotondim
  implicit none
  real(kind=8),intent(in)::x(1:ndim)
  real(kind=8),intent(out)::w1d(1:ndim,1:2),dw1d(1:ndim,1:2)
  integer,intent(out),optional::il_out(1:ndim),ir_out(1:ndim)
  real(kind=8),intent(out),optional::vol_out(1:twotondim)
  integer::idim,ind,ix,iy,iz
  real(kind=8)::xd
#if NDIM==3
  do idim=1,ndim
     xd = x(idim) + 0.5d0
     if(present(ir_out)) ir_out(idim) = int(xd)
     xd = xd - int(xd)  ! fractional offset in [0,1)
     w1d(idim,1) = 1.d0 - xd   ! left weight (dl)
     w1d(idim,2) = xd          ! right weight (dr)
     dw1d(idim,1) = -1.d0      ! derivative of left weight
     dw1d(idim,2) = +1.d0      ! derivative of right weight
     if(present(il_out).and.present(ir_out)) il_out(idim) = ir_out(idim) - 1
  end do
  ! 3D volume weights (optional output)
  if(present(vol_out))then
     ind = 0
     do iz=1,2
        do iy=1,2
           do ix=1,2
              ind = ind+1
              vol_out(ind) = w1d(1,ix)*w1d(2,iy)*w1d(3,iz)
           end do
        end do
     end do
  end if
#else
  w1d=0
  dw1d=0
#endif
end subroutine cic_weights_and_derivs

subroutine tsc_weights_and_derivs(x, w1d, dw1d, il_out, ic_out, ir_out, vol_out)
  !--------------------------------------------------------------
  ! Compute 1D TSC weights, derivatives, and cell indices from position.
  ! x(1:ndim) is the cell-centered coordinate.
  ! w1d(dim,1:3) = left, center, right weights
  ! dw1d(dim,1:3) = corresponding derivatives
  ! il_out, ic_out, ir_out = cell indices (optional)
  ! vol_out(1:threetondim) = 3D volume weights (optional)
  !--------------------------------------------------------------
  use amr_parameters, only: ndim, threetondim
  implicit none
  real(kind=8),intent(in)::x(1:ndim)
  real(kind=8),intent(out)::w1d(1:ndim,1:3),dw1d(1:ndim,1:3)
  integer,intent(out),optional::il_out(1:ndim),ic_out(1:ndim),ir_out(1:ndim)
  real(kind=8),intent(out),optional::vol_out(1:threetondim)
  integer::idim,il,ic,ir,ind,ix,iy,iz
  real(kind=8)::xd,xl,xc,xr
#if NDIM==3
  do idim=1,ndim
     xd = x(idim)
     il = int(xd)-1
     ic = int(xd)
     ir = int(xd)+1
     xl = dble(il)+0.5d0
     xc = dble(ic)+0.5d0
     xr = dble(ir)+0.5d0
     ! TSC weights
     w1d(idim,1) = 0.5d0*(1.5d0-abs(xd-xl))**2
     w1d(idim,2) = 0.75d0-(xd-xc)**2
     w1d(idim,3) = 0.5d0*(1.5d0-abs(xd-xr))**2
     ! TSC derivatives
     dw1d(idim,1) = -(1.5d0-abs(xd-xl))*sign(1.d0,xd-xl)
     dw1d(idim,2) = -2.d0*(xd-xc)
     dw1d(idim,3) = -(1.5d0-abs(xd-xr))*sign(1.d0,xd-xr)
     ! Cell indices (optional output)
     if(present(il_out)) il_out(idim) = il
     if(present(ic_out)) ic_out(idim) = ic
     if(present(ir_out)) ir_out(idim) = ir
  end do
  ! 3D volume weights (optional output)
  if(present(vol_out))then
     ind = 0
     do iz=1,3
        do iy=1,3
           do ix=1,3
              ind = ind+1
              vol_out(ind) = w1d(1,ix)*w1d(2,iy)*w1d(3,iz)
           end do
        end do
     end do
  end if
#else
  w1d=0
  dw1d=0
#endif
end subroutine tsc_weights_and_derivs

!#########################################################################
! Gradient computation subroutines for CIC and TSC interpolation
!#########################################################################
subroutine compute_gradient_cic_scalar(w1d, dw1d, field, grad_out)
  !--------------------------------------------------------------
  ! Compute gradient of a scalar field using CIC (2x2x2) weights.
  !--------------------------------------------------------------
  use amr_parameters, only: ndim, twotondim
  implicit none
  real(kind=8),intent(in)::w1d(1:ndim,1:2),dw1d(1:ndim,1:2)
  real(kind=8),intent(in)::field(1:twotondim)
  real(kind=8),intent(out)::grad_out(1:ndim)
  integer::ix,iy,iz,ind
  grad_out = 0.d0
#if NDIM==3
  ind = 0
  do iz=1,2
     do iy=1,2
        do ix=1,2
           ind = ind+1
           grad_out(1) = grad_out(1) + dw1d(1,ix)*w1d(2,iy)*w1d(3,iz)*field(ind)
           grad_out(2) = grad_out(2) + w1d(1,ix)*dw1d(2,iy)*w1d(3,iz)*field(ind)
           grad_out(3) = grad_out(3) + w1d(1,ix)*w1d(2,iy)*dw1d(3,iz)*field(ind)
        end do
     end do
  end do
#endif
end subroutine compute_gradient_cic_scalar

subroutine compute_gradient_cic(w1d, dw1d, field, grad_out)
  !--------------------------------------------------------------
  ! Compute gradient of a vector field using CIC (2x2x2) weights.
  ! grad_i = d(field_i)/dx_i.
  !--------------------------------------------------------------
  use amr_parameters, only: ndim, twotondim
  implicit none
  real(kind=8),intent(in)::w1d(1:ndim,1:2),dw1d(1:ndim,1:2)
  real(kind=8),intent(in)::field(1:ndim,1:twotondim)
  real(kind=8),intent(out)::grad_out(1:ndim)
  integer::ix,iy,iz,ind
  grad_out = 0.d0
#if NDIM==3
  ind = 0
  do iz=1,2
     do iy=1,2
        do ix=1,2
           ind = ind+1
           grad_out(1) = grad_out(1) + dw1d(1,ix)*w1d(2,iy)*w1d(3,iz)*field(1,ind)
           grad_out(2) = grad_out(2) + w1d(1,ix)*dw1d(2,iy)*w1d(3,iz)*field(2,ind)
           grad_out(3) = grad_out(3) + w1d(1,ix)*w1d(2,iy)*dw1d(3,iz)*field(3,ind)
        end do
     end do
  end do
#endif
end subroutine compute_gradient_cic

subroutine compute_gradient_tsc_scalar(w1d, dw1d, field, grad_out)
  !--------------------------------------------------------------
  ! Compute gradient of a scalar field using TSC (3x3x3) weights.
  !--------------------------------------------------------------
  use amr_parameters, only: ndim, threetondim
  implicit none
  real(kind=8),intent(in)::w1d(1:ndim,1:3),dw1d(1:ndim,1:3)
  real(kind=8),intent(in)::field(1:threetondim)
  real(kind=8),intent(out)::grad_out(1:ndim)
  integer::ix,iy,iz,ind
  grad_out = 0.d0
#if NDIM==3
  ind = 0
  do iz=1,3
     do iy=1,3
        do ix=1,3
           ind = ind+1
           grad_out(1) = grad_out(1) + dw1d(1,ix)*w1d(2,iy)*w1d(3,iz)*field(ind)
           grad_out(2) = grad_out(2) + w1d(1,ix)*dw1d(2,iy)*w1d(3,iz)*field(ind)
           grad_out(3) = grad_out(3) + w1d(1,ix)*w1d(2,iy)*dw1d(3,iz)*field(ind)
        end do
     end do
  end do
#endif
end subroutine compute_gradient_tsc_scalar

subroutine compute_gradient_tsc(w1d, dw1d, field, grad_out)
  !--------------------------------------------------------------
  ! Compute gradient of a vector field using TSC (3x3x3) weights.
  !--------------------------------------------------------------
  use amr_parameters, only: ndim, threetondim
  implicit none
  real(kind=8),intent(in)::w1d(1:ndim,1:3),dw1d(1:ndim,1:3)
  real(kind=8),intent(in)::field(1:ndim,1:threetondim)
  real(kind=8),intent(out)::grad_out(1:ndim)
  integer::ix,iy,iz,ind
  grad_out = 0.d0
#if NDIM==3
  ind = 0
  do iz=1,3
     do iy=1,3
        do ix=1,3
           ind = ind+1
           grad_out(1) = grad_out(1) + dw1d(1,ix)*w1d(2,iy)*w1d(3,iz)*field(1,ind)
           grad_out(2) = grad_out(2) + w1d(1,ix)*dw1d(2,iy)*w1d(3,iz)*field(2,ind)
           grad_out(3) = grad_out(3) + w1d(1,ix)*w1d(2,iy)*dw1d(3,iz)*field(3,ind)
        end do
     end do
  end do
#endif
end subroutine compute_gradient_tsc

!#########################################################################
! Helper routines for tracer particles
!#########################################################################
subroutine wrap_cell_coords(st,x_cell,levelp1)
  use amr_parameters, only: ndim
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t),intent(in)::st
  real(kind=8),intent(inout)::x_cell(1:ndim)
  integer,intent(in)::levelp1
  integer::jd
  real(kind=8)::range
#if NDIM==3
  do jd=1,ndim
     if(st%r%periodic(jd))then
        range=dble(st%m%box_ckey_max(jd,levelp1)-st%m%box_ckey_min(jd,levelp1))
        if(range<=0.d0)cycle
        if(x_cell(jd)< dble(st%m%box_ckey_min(jd,levelp1)))x_cell(jd)=x_cell(jd)+range
        if(x_cell(jd)>=dble(st%m%box_ckey_max(jd,levelp1)))x_cell(jd)=x_cell(jd)-range
     endif
  end do
#endif
end subroutine wrap_cell_coords

real(dp) function tracer_cell_kappa(dens_in,eturb_in,dx_in,smallr_in) result(kappa_val)
  use amr_parameters, only: dp
  implicit none
  real(dp),intent(in)::dens_in,eturb_in,dx_in,smallr_in
  real(dp)::rho_eff,sigma_sq
#if NDIM==3
  rho_eff = max(dens_in,smallr_in)
  sigma_sq = max(2.0_dp*max(eturb_in,0.0_dp)/rho_eff,0.0_dp)
  if(sigma_sq>0.0_dp)then
     kappa_val = dx_in*sqrt(sigma_sq)
  else
     kappa_val = 0.0_dp
  end if
#else
  kappa_val = 0.0_dp
#endif
end function tracer_cell_kappa

subroutine sample_tracer_gaussian(vec)
  use amr_parameters, only: ndim
  implicit none
  real(kind=8),intent(out)::vec(1:ndim)
  integer :: jd
  real(kind=8)::u_rand,tmp
  real(kind=8), external :: RngStream_RandUni
  external :: gaussdev
  vec=0.0d0
#if NDIM==3
  do jd=1,ndim
     u_rand = RngStream_RandUni(tracer_rng)
     call gaussdev(u_rand,tmp)
     vec(jd)=tmp
  end do
#endif
end subroutine sample_tracer_gaussian

subroutine sample_tracer_uniform(vec)
  use amr_parameters, only: ndim
  implicit none
  real(kind=8),intent(out)::vec(1:ndim)
  integer :: jd
  real(kind=8)::u_rand
  real(kind=8), external :: RngStream_RandUni
  vec=0.0d0
#if NDIM==3
  do jd=1,ndim
     u_rand = RngStream_RandUni(tracer_rng)
     ! Uniform distribution [-sqrt(3), sqrt(3)] with variance 1
     vec(jd) = (2.0d0*u_rand - 1.0d0) * sqrt(3.0d0)
  end do
#endif
end subroutine sample_tracer_uniform

subroutine sample_tracer_piecewise_skew_uniform(vec, gamma1_vec)
  !------------------------------------------------------------------------
  ! Sample from a two-piece uniform distribution with mean 0, variance 1,
  ! and skewness gamma1 (one value per dimension).
  !------------------------------------------------------------------------
  use amr_parameters, only: ndim
  implicit none
  real(kind=8),intent(out)::vec(1:ndim)
  real(kind=8),intent(in)::gamma1_vec(1:ndim)
  integer :: jd
  real(kind=8)::u_rand,gamma1,k,a,b,p_left,apb
  real(kind=8), external :: RngStream_RandUni
  vec=0.0d0
#if NDIM==3
  do jd=1,ndim
     gamma1 = gamma1_vec(jd)
     ! Compute two-piece uniform parameters from skewness
     k = 4.0d0 * gamma1 / 3.0d0
     a = 0.5d0 * (-k + sqrt(k*k + 12.0d0))
     b = a + k
     apb = a + b
     p_left = b / apb
     ! Sample using inverse CDF
     u_rand = RngStream_RandUni(tracer_rng)
     if (u_rand < p_left) then
        vec(jd) = -a + a * u_rand * apb / b
     else
        vec(jd) = b * (u_rand * apb - b) / a
     endif
  end do
#endif
end subroutine sample_tracer_piecewise_skew_uniform

real(kind=8) function mc_kernel_skewness(pr, pl) result(gamma1)
  !------------------------------------------------------------------------
  ! Compute the standardized skewness of the MC discrete kernel.
  !------------------------------------------------------------------------
  implicit none
  real(kind=8),intent(in)::pr, pl
  real(kind=8)::mu3, variance
#if NDIM==3
  mu3 = (pr - pl) * (1.0d0 - 3.0d0*pl + 4.0d0*pl**2 - 2.0d0*pl**3 &
       & - 3.0d0*pr - 8.0d0*pl*pr + 2.0d0*pl**2*pr &
       & + 4.0d0*pr**2 + 2.0d0*pl*pr**2 - 2.0d0*pr**3)
  variance = pl - pl**2 + pr + 2.0d0*pl*pr - pr**2
  if (variance > 0.0d0) then
     gamma1 = mu3 / (variance**1.5d0)
  else
     gamma1 = 0.0d0
  endif
#else
  gamma1 = 0.0d0
#endif
end function mc_kernel_skewness

subroutine gather_cic_state(st,x_cell,level_in,dx_cell,use_sgs_in,vel_out,kappa_out)
  use amr_parameters, only: ndim, twotondim, dp
  use ramses_commons, only: ramses_t
  use nbors_utils
  use cache
  implicit none
  type(ramses_t),intent(in)::st
  real(kind=8),intent(in)::x_cell(1:ndim)
  integer,intent(in)::level_in
  real(kind=8),intent(in)::dx_cell
  logical,intent(in)::use_sgs_in
  real(kind=8),intent(out)::vel_out(1:ndim)
  real(kind=8),intent(out)::kappa_out
  integer,dimension(1:ndim)::il,ir
  real(kind=8),dimension(1:twotondim)::vol
  integer,dimension(1:ndim,1:twotondim)::ckey
  real(kind=8),dimension(1:ndim,1:2)::w1d,dw1d
  real(kind=8),dimension(1:ndim)::momentum
  real(kind=8)::rho,kappa_sum
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer::ind,icell,jd,igrid
  vel_out=0.d0
  momentum=0.d0
  rho=0.d0
  kappa_sum=0.d0
  kappa_out=0.d0
#if NDIM==3
  ! Build CIC weights, indices, and volume weights
  call cic_weights_and_derivs(x_cell, w1d, dw1d, il, ir, vol)
  do jd=1,ndim
     if(st%r%periodic(jd))then
        if(il(jd)< st%m%box_ckey_min(jd,level_in+1))il(jd)=st%m%box_ckey_max(jd,level_in+1)-1
        if(ir(jd)>=st%m%box_ckey_max(jd,level_in+1))ir(jd)=st%m%box_ckey_min(jd,level_in+1)
     endif
  end do
  ckey = cic_index(il,ir)

  hash_nbor(0)=level_in+1
  do ind=1,twotondim
     hash_nbor(1:ndim)=ckey(1:ndim,ind)
     call get_parent_cell(st,hash_nbor,igrid,icell,flush_cache=.false.,fetch_cache=.true.)
#ifdef HYDRO
     if(igrid>0)then
        momentum(1:ndim)=momentum(1:ndim)+st%m%uold(icell,2:ndim+1,igrid)*vol(ind)
        rho=rho+st%m%uold(icell,1,igrid)*vol(ind)
        if(use_sgs_in)then
           kappa_sum=kappa_sum+tracer_cell_kappa(st%m%uold(icell,1,igrid),st%m%uold(icell,st%r%iturb,igrid),real(dx_cell,dp),real(st%r%smallr,dp))*vol(ind)
        end if
     end if
#endif
  end do

  if(rho>0.d0)then
     vel_out(1:ndim)=momentum(1:ndim)/max(rho,st%r%smallr)
  else
     vel_out(1:ndim)=0.d0
  end if

  if(use_sgs_in)then
     kappa_out=max(kappa_sum,0.d0)
  else
     kappa_out=0.d0
  end if
#endif
end subroutine gather_cic_state

subroutine gather_tsc_state(st,x_cell,level_in,dx_cell,use_sgs_in,vel_out,kappa_out)
  use amr_parameters, only: ndim, threetondim, dp
  use ramses_commons, only: ramses_t
  use nbors_utils
  use cache
  use rho_fine_module, only: tsc_index
  implicit none
  type(ramses_t),intent(in)::st
  real(kind=8),intent(in)::x_cell(1:ndim)
  integer,intent(in)::level_in
  real(kind=8),intent(in)::dx_cell
  logical,intent(in)::use_sgs_in
  real(kind=8),intent(out)::vel_out(1:ndim)
  real(kind=8),intent(out)::kappa_out
  integer,dimension(1:ndim)::il,ic,ir
  real(kind=8),dimension(1:threetondim)::vol
  integer,dimension(1:ndim,1:threetondim)::ckey
  real(kind=8),dimension(1:ndim,1:3)::w1d,dw1d
  real(kind=8),dimension(1:ndim)::momentum
  real(kind=8)::rho,kappa_sum
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer::ind,icell,jd,igrid
  vel_out=0.d0
  momentum=0.d0
  rho=0.d0
  kappa_sum=0.d0
  kappa_out=0.d0
#if NDIM==3
  ! Build TSC weights, indices, and volume weights
  call tsc_weights_and_derivs(x_cell, w1d, dw1d, il, ic, ir, vol)
  do jd=1,ndim
     if(st%r%periodic(jd))then
        if(il(jd)< st%m%box_ckey_min(jd,level_in+1))il(jd)=il(jd)-st%m%box_ckey_min(jd,level_in+1)+st%m%box_ckey_max(jd,level_in+1)
        if(ir(jd)>=st%m%box_ckey_max(jd,level_in+1))ir(jd)=ir(jd)+st%m%box_ckey_min(jd,level_in+1)-st%m%box_ckey_max(jd,level_in+1)
     endif
  end do
  ckey = tsc_index(il,ic,ir)

  hash_nbor(0)=level_in+1
  do ind=1,threetondim
     hash_nbor(1:ndim)=ckey(1:ndim,ind)
     call get_parent_cell(st,hash_nbor,igrid,icell,flush_cache=.false.,fetch_cache=.true.)
#ifdef HYDRO
     if(igrid>0)then
        momentum(1:ndim)=momentum(1:ndim)+st%m%uold(icell,2:ndim+1,igrid)*vol(ind)
        rho=rho+st%m%uold(icell,1,igrid)*vol(ind)
        if(use_sgs_in)then
           kappa_sum=kappa_sum+tracer_cell_kappa(st%m%uold(icell,1,igrid),st%m%uold(icell,st%r%iturb,igrid),real(dx_cell,dp),real(st%r%smallr,dp))*vol(ind)
        end if
     end if
#endif
  end do

  if(rho>0.d0)then
     vel_out(1:ndim)=momentum(1:ndim)/max(rho,st%r%smallr)
  else
     vel_out(1:ndim)=0.d0
  end if

  if(use_sgs_in)then
     kappa_out=max(kappa_sum,0.d0)
  else
     kappa_out=0.d0
  end if
#endif
end subroutine gather_tsc_state

!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
end module move_fine_module
