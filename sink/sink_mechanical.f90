module sink_mechanical_module
  use sink_evolution_module, only: out_accretion_t, bondi_alpha, dynamical_friction, &
       & dump_sink_data_fine, dump_sink_data_fine_AGN
#ifdef _CUDA
  use sink_device, only: gpu_sink_mechanical
#endif
  type :: in_sink_mechanical_t
     integer::ilevel
     real(kind=8)::dt
     real(kind=8)::texp_old
  end type in_sink_mechanical_t
contains
  !##############################################################################
  !##############################################################################
  !##############################################################################
  !##############################################################################
  subroutine m_sink_mechanical(pst,ilevel,macc)
    use mdl_module
    use ramses_commons, only: pst_t
    use mdl_parameters
    implicit none
    type(pst_t)::pst
    integer::ilevel
    real(kind=8)::macc
    ! This is the master routine calling the 56-point sink evolution.
    ! Note that it loops over all levels that are time-synchronized,
    ! as m_mechanical_feedback does.
    integer::i
    type(in_sink_mechanical_t)::input_acc
    type(out_accretion_t)::output_acc

    macc=0d0
    input_acc%dt=pst%s%g%dtnew(ilevel)
    do i=ilevel,pst%s%r%levelmin,-1
       ! The time step of a coarser level is the sum of the sub-steps of the
       ! finer level, as set at the end of amr_step (dtnew(i) is not yet updated)
       if(i<ilevel)then
          if(pst%s%r%nsubcycle(i)==2)input_acc%dt=pst%s%g%dtold(i+1)+input_acc%dt
       endif
       input_acc%ilevel=i
       input_acc%texp_old=pst%s%g%texp_old(i)
       call r_sink_mechanical(pst,input_acc,storage_size(input_acc)/32,output_acc,2)
       macc=macc+output_acc%mass
       if(i>pst%s%r%levelmin)then
          if(pst%s%r%nsubcycle(i-1)==2.and.pst%s%g%isubcycle(i)==1)exit
       endif
    end do

  end subroutine m_sink_mechanical
  !##############################################################################
  !##############################################################################
  !##############################################################################
  !##############################################################################
  recursive subroutine r_sink_mechanical(pst,input,input_size,output,output_size)
    use mdl_module
    use ramses_commons, only: pst_t
    use mdl_parameters
    implicit none
    type(pst_t)::pst
    integer,VALUE::input_size
    integer::output_size
    type(in_sink_mechanical_t)::input
    type(out_accretion_t)::output,next_output

    integer::rID

    if(pst%nLower>0)then
       rID = mdl_send_request(pst%s%mdl,MDL_SINK_MECHANICAL,pst%iUpper+1,input_size,output_size,input)
       call r_sink_mechanical(pst%pLower,input,input_size,output,output_size)
       call mdl_get_reply(pst%s%mdl,rID,output_size,next_output)
       output%mass=output%mass+next_output%mass
    else
#ifdef _CUDA
       call gpu_sink_mechanical(pst%s,input%ilevel,input%dt,output%mass)
#else
       call sink_mechanical(pst%s,pst%s%sink,input%ilevel,input%dt,input%texp_old,output%mass)
#endif
    endif

  end subroutine r_sink_mechanical
  !##############################################################################
  !##############################################################################
  !##############################################################################
  !##############################################################################
  subroutine sink_mechanical(s,p,ilevel,dt_loc,texp_old,macc_loc)
    use constants
    use amr_parameters, only: ndim, twotondim
    use hydro_parameters, only: nvar, nener
    use ramses_commons, only: ramses_t
    use pm_commons, only: part_t, cross
    use params_module
    use nbors_utils
    use cache_commons
    use cache
    use marshal, only: pack_fetch_refine,unpack_fetch_refine
    use boundaries, only: init_bound_refine
    use godunov_fine_module, only: init_flush_godunov,pack_flush_godunov,unpack_flush_godunov
    use hilbert
    implicit none
    type(ramses_t)::s
    type(part_t)::p
    integer::ilevel
    real(kind=8)::dt_loc,texp_old,macc_loc
    !==================================================================
    ! This is the RAMSES routine for sink evolution using the 56-point
    ! stencil of the mechanical supernova feedback. Gas properties are
    ! sampled in uold and accretion and AGN feedback are applied to unew
    ! across levels ilevel-1, ilevel and ilevel+1, so that sinks do not
    ! need sink_refine=.true.
    !==================================================================
    ! Parameters
    ! Number of sample points (4x4x4 half cells minus the 8 corners)
    integer,parameter::nSKnei=56
    ! Local variables
    real(kind=8),dimension(1:3,1:nSKnei)::xSKnei
    integer,dimension(1:nSKnei)::igrid_nbor,icell_nbor,level_nbor
    real(kind=8),dimension(1:nSKnei)::d_nbor,e_nbor,vol_nbor,w_nbor,mu_nbor
    real(kind=8),dimension(1:3,1:nSKnei)::v_nbor
    real(kind=8)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,factG,unit_yr
    real(kind=8)::dx_loc,vol_loc,dx_cell,vol_cell,vol_pt,current_time
    real(kind=8)::lambda_sonic,r_acc
    real(kind=8)::x,y,z
    integer,dimension(1:ndim)::ckey,ckey_nbor
    integer(kind=8),dimension(0:ndim)::hash_cell,hash_nbor
    integer::i,j,k,ind,ipart,idim,ivar,istep,istep_old,igrid,icell,igridn,icelln,ilev_cell
    real(kind=8),dimension(1:3)::xcen,xnei,xsink,x_rel,vv,v_rel,jet_direction,p_jet
    real(kind=8)::d,e,ethermal,cs,dm,dvol
    real(kind=8)::rho_gas,cs_gas,m_gas,rho_av_all,rho_inf
    real(kind=8),dimension(1:ndim)::vel_gas
    real(kind=8)::velocity,v_bondi,bondi_mass,r2_sink
    real(kind=8)::dMBH_overdt,dMED_overdt
    real(kind=8)::m_acc
    real(kind=8),dimension(1:3)::x_acc,p_acc,l_acc
    real(kind=8)::acc_ratio,fbk_ener_agn,fbk_mom_agn
    real(kind=8)::w_plus,w_minus,de
    real(kind=8)::v_max,e_max,u_rel,lambda_load,eta,m_amb,m_wind,m_bh,m_pt
    real(kind=8),dimension(1:3)::vp_old,jp_new,v_new
    real(kind=8),dimension(1:nvar)::q_acc,q_wind
    logical::quasar_mode,ok_level,ok_leaf,ok_face,output_file
    type(msg_large_realdp)::dummy_large_realdp
#if NENER>0
    integer::irad
#endif
#ifdef MHD
    real(kind=8)::bx,by,bz,emag
#endif

#ifdef HYDRO
#if NDIM==3
    associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

    if(g%myid==1.and.r%verbose)write(*,*)'Entering sink_mechanical...'

    ! Arrays to define the sample points (center=[0,0,0])
    ! normalized to dx = 1 = size of the host cell in which the sink sits
    ! from -0.75 to 0.75
    ind=0
    do k=1,4
       do j=1,4
          do i=1,4
             if((i==1.or.i==4).and.(j==1.or.j==4).and.(k==1.or.k==4))cycle ! edge
             ind=ind+1
             x=(i-1)+0.5d0-2
             y=(j-1)+0.5d0-2
             z=(k-1)+0.5d0-2
             xSKnei(1,ind)=x/2d0
             xSKnei(2,ind)=y/2d0
             xSKnei(3,ind)=z/2d0
          end do
       end do
    end do

    ! Conversion factor from user units to cgs units
    call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

    ! Gravitational constant
    factG=1
    if(r%cosmo)factG=3d0/4d0/twopi*g%omega_m*g%aexp

    ! Bondi sonic constant
    if(r%use_bondi_lambda)then
       if(abs(r%gamma - 1).le.0.01)then
          lambda_sonic = 0.25d0*exp(1.5d0)
       else if(abs(r%gamma - 5.0d0/3.0d0).le.0.01)then
          lambda_sonic = 0.25d0
       else
          lambda_sonic = 0.5d0**((r%gamma + 1)/(2d0*(r%gamma - 1))) &
               & * (0.25d0*(5d0-3d0*r%gamma))**(-(5d0-3d0*r%gamma)/(2d0*(r%gamma - 1)))
       end if
    else
       lambda_sonic = 1.0d0
    end if

    ! Maximum velocity kick (radio mode) and sound speed (quasar mode) of the AGN feedback in code units
    v_max=r%agn_jet_vmax*1d5/scale_v

    ! Mesh spacing in that level
    dx_loc=r%boxlen/2**ilevel
    vol_loc=dx_loc**ndim

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Check if high frequency dump
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Dump the sinks of this level if a multiple of sink_delta_tout falls
    ! within the time step of this level, as for supernova events, so that
    ! each sink is dumped once whatever its level
    output_file = .false.
    if(r%sink_delta_tout>0)then
       unit_yr = 3600*24*365.25
       ! proper time in Myr (careful in cosmology it is negative)
       current_time = g%texp * scale_t / g%aexp**2 / unit_yr / 1e6
       istep = floor(current_time / r%sink_delta_tout)
       current_time = texp_old * scale_t / g%aexp**2 / unit_yr / 1e6
       istep_old = floor(current_time / r%sink_delta_tout)
       if(istep > istep_old)then
          p%step_counter = istep
          output_file = .true.
       endif
    endif

    ! Open cache for array uold (fetch) and unew (flush)
    call open_cache(mdl, m, pack_size=storage_size(dummy_large_realdp)/32, &
         pack=pack_fetch_refine, unpack=unpack_fetch_refine, &
         init=init_flush_godunov, flush=pack_flush_godunov, &
         combine=unpack_flush_godunov, bound=init_bound_refine)

    ! Loop over sink particles
    macc_loc=0d0
    do ipart=p%headp(ilevel),p%tailp(ilevel)

       ! Skip zero-mass (merged/defunct) sinks
       if(p%mp(ipart)<=0.0d0)cycle

       ok_level=.true.

       ! Find parent cell at level ilevel
       do idim=1,ndim
          ckey(idim)=int((p%xp(ipart,idim)+m%skip(idim))/dx_loc)
          xcen(idim)=ckey(idim)+0.5
       end do

       ! Cell size at level ilevel
       dx_cell=dx_loc

       ! Get parent cell at level ilevel using cache
       hash_cell(0)=ilevel+1
       hash_cell(1:ndim)=ckey(1:ndim)
       call get_parent_cell(s,hash_cell,igrid,icell,flush_cache=.true.,fetch_cache=.true.)

       ! If cell does not exist at current level, then find cell at coarser level
       if(igrid==0)then

          ! NGP at level ilevel-1
          do idim=1,ndim
             ckey(idim)=int((p%xp(ipart,idim)+m%skip(idim))/dx_loc/2)
             xcen(idim)=ckey(idim)+0.5
          end do

          ! Cell size at level ilevel-1
          dx_cell=dx_loc*2

          ! Get parent cell at level ilevel-1 using cache
          hash_cell(0)=ilevel
          hash_cell(1:ndim)=ckey(1:ndim)
          call get_parent_cell(s,hash_cell,igrid,icell,flush_cache=.true.,fetch_cache=.true.)
          if(igrid==0)ok_level=.false.

       end if

       ! No accretion or feedback from a host cell touching a non-periodic face
       ok_face=.true.
       do idim=1,ndim
          if(.not.r%periodic(idim))then
             if(ckey(idim)==m%box_ckey_min(idim,hash_cell(0)))ok_face=.false.
             if(ckey(idim)==m%box_ckey_max(idim,hash_cell(0))-1)ok_face=.false.
          endif
       end do
       if(.not.ok_face)cycle

       if(.not. ok_level)then
          write(*,*)"Something went wrong in sink_mechanical"
          write(*,*)"Current level grid and coarser grid both dont exist..."
          stop
       endif

       ! Lock grid in cache
       call lock_cache(m,igrid)

       ok_leaf = .not. m%grid(igrid)%refined(icell)
       if(.not. ok_leaf)then
          write(*,*)"Something went wrong in sink_mechanical"
          write(*,*)"Cell should be a leaf cell..."
          stop
       endif

       ! Host cell properties
       ilev_cell=hash_cell(0)-1
       vol_cell=dx_cell**ndim
       vol_pt=vol_cell/dble(twotondim)
       xsink(1:ndim)=(p%xp(ipart,1:ndim)+m%skip(1:ndim))/dx_cell

       ! Initialize outputs so drag, feedback and dumps
       ! read defined values even when accretion is disabled
       dMBH_overdt=0.0d0; dMED_overdt=0.0d0; rho_inf=0.0d0
       fbk_mom_agn=0.0d0; fbk_ener_agn=0.0d0

       !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       ! Collect all sample points from hash table
       !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
       rho_gas=0d0; vel_gas=0d0; cs_gas=0d0; m_gas=0d0; rho_av_all=0d0
       do j=1,nSKnei

          ! Compute sample point coordinates
          xnei(1:ndim)=xcen(1:ndim)+xSKnei(1:ndim,j)
          ! Periodic boundary conditions
          do idim=1,ndim
             if(r%periodic(idim))then
                if(xnei(idim)< m%box_ckey_min(idim,hash_cell(0)))xnei(idim)=xnei(idim)-m%box_ckey_min(idim,hash_cell(0))+m%box_ckey_max(idim,hash_cell(0))
                if(xnei(idim)>=m%box_ckey_max(idim,hash_cell(0)))xnei(idim)=xnei(idim)+m%box_ckey_min(idim,hash_cell(0))-m%box_ckey_max(idim,hash_cell(0))
             endif
          end do

          ! Get neighboring cell at current level
          ckey_nbor(1:ndim)=int(xnei(1:ndim))
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
          else if (m%grid(igridn)%refined(icelln))then
             ckey_nbor(1:ndim)=int(xnei(1:ndim)*2.0)
             hash_nbor(0)=hash_cell(0)+1
             hash_nbor(1:ndim)=ckey_nbor(1:ndim)
             call get_parent_cell(s,hash_nbor,igridn,icelln,flush_cache=.true.,fetch_cache=.true.)
          endif

          ! Lock grid in cache
          call lock_cache(m,igridn)

          igrid_nbor(j) = igridn
          icell_nbor(j) = icelln
          level_nbor(j) = hash_nbor(0)-1
          vol_nbor(j) = vol_cell*dble(twotondim)**(ilev_cell-level_nbor(j))

          ! Get physical information
          d     = max(dble(m%uold(icelln,1,igridn)),r%smallr)
          vv(1) = m%uold(icelln,2,igridn)/d
          vv(2) = m%uold(icelln,3,igridn)/d
          vv(3) = m%uold(icelln,4,igridn)/d
          e     = m%uold(icelln,5,igridn)

          ! We need to remove all non-thermal energies
#ifdef MHD
          bx=0.5d0*(m%bold(icelln,1,igridn) + m%bold(icelln,4,igridn))
          by=0.5d0*(m%bold(icelln,2,igridn) + m%bold(icelln,5,igridn))
          bz=0.5d0*(m%bold(icelln,3,igridn) + m%bold(icelln,6,igridn))
          emag=0.5d0*(bx**2+by**2+bz**2)
          e = e - emag
#endif
#if NENER>0
          do irad=1,nener
             e = e - m%uold(icelln,5+irad,igridn)
          end do
#endif
          ethermal = (e - 0.5d0*d*sum(vv(:)**2)) / d
          cs       = sqrt(max(dble(r%gamma)*(dble(r%gamma)-1.0d0)*ethermal,dble(r%smallc)**2) &
               &   * dble(r%acc_sink_boost)**(-2d0/3d0))

          d_nbor(j)     = d
          e_nbor(j)     = ethermal
          v_nbor(1:3,j) = vv(1:3)

          ! Add to average information (each point has the same volume)
          rho_gas         = rho_gas         + d          / dble(nSKnei)
          vel_gas(1:ndim) = vel_gas(1:ndim) + vv(1:ndim) / dble(nSKnei)
          cs_gas          = cs_gas          + cs         / dble(nSKnei)
          m_gas           = m_gas           + d          * vol_pt
#ifdef GRAV
          rho_av_all      = rho_av_all      + m%rho(icelln,igridn) / dble(nSKnei)
#endif
       end do

       !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       ! Compute overall accretion rate
       !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
       if(r%accretion_type==1)then

          ! Bondi velocity
          v_rel(1:ndim) = vel_gas(1:ndim) - p%vp(ipart,1:ndim)
          if(r%bondi_use_vrel)then
             velocity   = norm2(v_rel(1:ndim))
             v_bondi    = sqrt(velocity**2 + cs_gas**2)
          else
             v_bondi    = cs_gas
          end if

          ! Bondi mass
          if(r%bondi_use_gas_mass)then
             bondi_mass = p%mp(ipart) + m_gas
          else
             bondi_mass = p%mp(ipart)
          end if

          ! Bondi radius
          r2_sink     = (factG * bondi_mass / v_bondi**2)**2

          ! Density at infinity (using extrapolation from the radius
          ! of the sphere with the stencil volume 7 dx^3)
          if(r%use_rho_inf)then
             r_acc   = (21d0/(4d0*pi))**(1d0/3d0)*dx_cell
             rho_inf = rho_gas / (bondi_alpha(r_acc/(r2_sink+tiny(0.0d0))**0.5d0))
          else
             rho_inf = rho_gas
          end if

          ! Bondi-Hoyle-Lyttleton accretion rate
          dMBH_overdt = 4.0d0 * pi * rho_inf * r2_sink * v_bondi * lambda_sonic
          if(r%verbose_sink)write(*,*)'Bondi: ',dMBH_overdt

       ! Compute threshold accretion rate
       else if(r%accretion_type==3)then
          ! Accrete in one time step all the gas above the effective peak density for sink formation
          dMBH_overdt = (rho_gas - r%sink_relevance_threshold*r%sink_density_threshold)*vol_pt*dble(nSKnei) &
               & / dt_loc
          dMBH_overdt = max(dMBH_overdt, 0.0d0)
          if(r%verbose_sink)write(*,*)'Density: ',dMBH_overdt
       end if

       !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       ! Limit overall accretion rate
       !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
       if(r%accretion_type>0)then

          ! Eddington accretion rate, which introduces an optional cap
          dMED_overdt = 4.0d0 * pi * factG_in_cgs * p%mp(ipart) * mH / (0.1d0 * sigma_T * c_cgs) * scale_t
          if(r%eddington_cap>0)dMBH_overdt = min(dMBH_overdt, dMED_overdt*r%eddington_cap)

          ! If the accretion rate is too low, do nothing
          if((r%eddington_floor>0).and.(dMBH_overdt/dMED_overdt<r%eddington_floor))dMBH_overdt=0d0

          ! Limiting total accreted mass to 75% of the mass in the stencil
          dMBH_overdt = min(dMBH_overdt, 0.75d0*m_gas / dt_loc)

          ! tp is a proper time: compare it with texp, not with the code time t
          if(((g%texp - p%tp(ipart)).lt.r%t_start_black_hole).and.r%t_start_black_hole.gt.0.0)then
             dMBH_overdt = exp(g%texp - p%tp(ipart) - r%t_start_black_hole) * dMBH_overdt
          end if

          if(r%manual_accretion_rate.gt.0.0d0)then
             dMBH_overdt = r%manual_accretion_rate * dMED_overdt
          end if

       end if

       !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       ! Accrete from the sample points
       !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
       m_acc=0.0d0; x_acc=0.0d0; p_acc=0.0d0; l_acc=0.0d0; q_acc=0.0d0
       if(dMBH_overdt>0.0d0)then
          do j=1,nSKnei

             igridn = igrid_nbor(j)
             icelln = icell_nbor(j)

             ! Get physical information
             d          = d_nbor(j)
             vv(1:ndim) = v_nbor(1:ndim,j)
             e          = m%uold(icelln,5,igridn)/d

             ! We need to remove all non-thermal energies as they are not accreted
#ifdef MHD
             bx = 0.5d0*(m%bold(icelln,1,igridn) + m%bold(icelln,4,igridn))
             by = 0.5d0*(m%bold(icelln,2,igridn) + m%bold(icelln,5,igridn))
             bz = 0.5d0*(m%bold(icelln,3,igridn) + m%bold(icelln,6,igridn))
             emag = 0.5d0*(bx**2 + by**2 + bz**2)
             e = e - emag/d
#endif
#if NENER>0
             do irad = 1, nener
                e = e - m%uold(icelln,5+irad,igridn)/d
             end do
#endif
             ! Get accreted mass for this sample point
             dm = dMBH_overdt * dt_loc / dble(nSKnei)
             if (r%mass_weighting) dm = dm * d / rho_gas

             ! Ensure that the accreted mass is positive and never exceeds
             ! 75% of the gas content of this sample point (avoids negative densities)
             dm = max(dm, 0.0d0)
             dm = min(dm, 0.75d0*d*vol_pt)

             ! Accrete from the cell
             dvol = dm / vol_nbor(j)
             m%unew(icelln,1,igridn)          = m%unew(icelln,1,igridn)          - dvol
             m%unew(icelln,2:(ndim+1),igridn) = m%unew(icelln,2:(ndim+1),igridn) - dvol * vv(1:ndim)
             m%unew(icelln,5,igridn)          = m%unew(icelln,5,igridn)          - dvol * e

             ! Accrete passive scalars
#if NVAR>5+NENER
             do ivar = 6 + nener, nvar
                m%unew(icelln,ivar,igridn) = m%unew(icelln,ivar,igridn) - dvol*m%uold(icelln,ivar,igridn)/d
                q_acc(ivar) = q_acc(ivar) + dm*m%uold(icelln,ivar,igridn)/d
             end do
#endif
             ! Accreted mass
             m_acc = m_acc + dm

             if(r%momentum_conserving)then
                ! Sample point position and velocity relative to the sink
                x_rel(1:ndim) = (xcen(1:ndim) + xSKnei(1:ndim,j) - xsink(1:ndim)) * dx_cell
                v_rel(1:ndim) = vv(1:ndim) - p%vp(ipart,1:ndim)

                ! Accreted relative center of mass
                x_acc(1:ndim) = x_acc(1:ndim) + dm * x_rel(1:ndim)

                ! Accreted relative momentum
                p_acc(1:ndim) = p_acc(1:ndim) + dm * v_rel(1:ndim)

                ! Accreted relative angular momentum
                l_acc(1:ndim) = l_acc(1:ndim) + dm * cross(x_rel(1:ndim), v_rel(1:ndim))
             endif

          end do
       end if

       !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       ! AGN mode and mass loading
       !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
       m_wind=0.0d0; q_wind=0.0d0
       if(r%agn.and.dMBH_overdt>0.0d0)then

          ! Compute chi (fraction of Eddington)
          acc_ratio = 0.0d0
          if(dMED_overdt>0.0d0)acc_ratio = dMBH_overdt/dMED_overdt

          ! Direction of the sink angular momentum after accretion
          jp_new(1:ndim) = p%mp(ipart) * p%jp(ipart,1:ndim) + l_acc(1:ndim)

          ! Check if quasar mode
          quasar_mode = (acc_ratio.gt.r%agn_fbk_mode_switch_threshold) .or. (norm2(jp_new(1:ndim)) .le. 0.0d0)

          ! Compute the jet direction
          jet_direction(1:ndim) = jp_new(1:ndim) / (norm2(jp_new(1:ndim)) + tiny(0.0d0))

          ! Fraction of the feedback given to each sample point: all points in
          ! quasar mode, half of the feedback on either side of the jet equatorial
          ! plane through the host cell centre in radio mode (none in the plane)
          w_nbor=0d0; mu_nbor=0d0
          do j=1,nSKnei
             if(.not.quasar_mode)then
                mu_nbor(j) = dot_product(xSKnei(1:ndim,j),jet_direction(1:ndim))
                ! Points in the plane (tolerance: FMA round-off when |x|=|z|)
                if(abs(mu_nbor(j))<1.0d-6)then
                   mu_nbor(j)=0.0d0; cycle
                endif
             endif
             w_nbor(j) = 1d0
             if(r%agn_use_mass_weighting)w_nbor(j) = d_nbor(j)
          end do
          if(quasar_mode)then
             w_nbor = w_nbor / sum(w_nbor)
          else
             w_plus  = sum(w_nbor, mask=mu_nbor>0d0)
             w_minus = sum(w_nbor, mask=mu_nbor<0d0)
             where(mu_nbor>0d0)w_nbor = 0.5d0 * w_nbor / w_plus
             where(mu_nbor<0d0)w_nbor = 0.5d0 * w_nbor / w_minus
          endif

          ! Ambient mass seen by the feedback. It is reduced by the thermal energy
          ! (quasar mode) or the velocity along the jet relative to the sink (radio
          ! mode) already in the gas, so that the sound speed after heating or the
          ! velocity after the kick does not exceed agn_jet_vmax (gas already
          ! hotter or faster does not get hotter or faster)
          e_max = v_max**2 / ( dble(r%gamma) * (dble(r%gamma) - 1.0d0) )
          m_amb = huge(1.0d0)
          do j=1,nSKnei
             if(w_nbor(j)==0.0d0)cycle
             if(quasar_mode)then
                m_amb = min(m_amb, d_nbor(j) * vol_pt / w_nbor(j) * ( 1.0d0 - e_nbor(j) / e_max ))
             else
                u_rel = sign(1.0d0,mu_nbor(j)) * dot_product(v_nbor(1:ndim,j) - p%vp(ipart,1:ndim), jet_direction(1:ndim))
                m_amb = min(m_amb, d_nbor(j) * vol_pt / w_nbor(j) * ( 1.0d0 - u_rel / v_max ))
             endif
          end do
          m_amb = max(m_amb, 0.0d0)

          ! Mass loading: put back into the gas the fraction eta/(1+eta) of the
          ! accreted mass, so that the radio mode velocity kick or the quasar mode
          ! sound speed does not exceed agn_jet_vmax
          if(quasar_mode)then
             lambda_load = dble(r%gamma) * (dble(r%gamma) - 1.0d0) * r%epsilon_quasar * r%epsilon_rad &
                  &      * ( (c_cgs/scale_v) / v_max )**2
          else
             lambda_load = r%epsilon_radio * r%epsilon_rad * (c_cgs/scale_v) / v_max
          endif
          eta = max(0.0d0, (lambda_load * m_acc - m_amb) / (m_amb + m_acc))
          m_wind = eta / (1.0d0 + eta) * m_acc
          q_wind = eta / (1.0d0 + eta) * q_acc

       end if

       ! Mass reaching the black hole
       m_bh = m_acc - m_wind

       ! The wind mass is ejected with the sink velocity before accretion
       vp_old(1:ndim) = p%vp(ipart,1:ndim)

       ! Add accreted properties to sink variables
       if(.not.p%static)p%xp(ipart,1:ndim) = p%xp(ipart,1:ndim) + x_acc(1:ndim) / ( p%mp(ipart) + m_bh )
       p%vp(ipart,1:ndim)                  = p%vp(ipart,1:ndim) + p_acc(1:ndim) / ( p%mp(ipart) + m_bh )
       p%jp(ipart,1:ndim)                  = ( p%mp(ipart) * p%jp(ipart,1:ndim) + l_acc(1:ndim) ) / ( p%mp(ipart) + m_bh )
       if(.not.r%sink_fix_mass)p%mp(ipart) = p%mp(ipart) + m_bh

       ! Save accreted mass to total
       macc_loc = macc_loc + m_bh

       !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       ! Dynamics
       !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
       if(r%drag_sink)then
          call dynamical_friction(s,p,dt_loc,ipart,rho_av_all,cs_gas,vel_gas,factG)
       endif

       !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       ! Sink/AGN Feedback
       !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
       if(r%agn.and.dMBH_overdt>0.0d0)then

          if(quasar_mode)then
             ! Quasar mode thermal energy, powered by the mass reaching the black hole
             fbk_ener_agn = r%epsilon_quasar * r%epsilon_rad * m_bh * (c_cgs/scale_v)**2
          else
             ! Radio mode momentum in units of L/c, powered by the mass reaching the black hole
             fbk_mom_agn  = r%epsilon_radio * r%epsilon_rad * m_bh * (c_cgs/scale_v)
          endif

          ! Wind mass, momentum (radio mode, along the jet direction, opposite
          ! on either side of the equatorial plane) and energy of each sample point
          do j=1,nSKnei
             if(w_nbor(j)==0.0d0)cycle
             igridn = igrid_nbor(j)
             icelln = icell_nbor(j)
             dm = w_nbor(j) * m_wind
             p_jet(1:ndim) = sign(w_nbor(j),mu_nbor(j)) * fbk_mom_agn * jet_direction(1:ndim)

             ! Energy: kinetic energy change of the gas of the sample point, so that
             ! its thermal energy is unchanged, plus the quasar mode heating
             m_pt = d_nbor(j) * vol_pt
             v_new(1:ndim) = ( m_pt * v_nbor(1:ndim,j) + dm * vp_old(1:ndim) + p_jet(1:ndim) ) / ( m_pt + dm )
             de = 0.5d0 * ( m_pt + dm ) * sum(v_new(1:ndim)**2) - 0.5d0 * m_pt * sum(v_nbor(1:ndim,j)**2) &
                  & + w_nbor(j) * fbk_ener_agn

             m%unew(icelln,1,igridn)   = m%unew(icelln,1,igridn)   + dm / vol_nbor(j)
             m%unew(icelln,2:4,igridn) = m%unew(icelln,2:4,igridn) + ( dm * vp_old(1:ndim) + p_jet(1:ndim) ) / vol_nbor(j)
             m%unew(icelln,5,igridn)   = m%unew(icelln,5,igridn)   + de / vol_nbor(j)
#if NVAR>5+NENER
             do ivar = 6 + nener, nvar
                m%unew(icelln,ivar,igridn) = m%unew(icelln,ivar,igridn) + w_nbor(j) * q_wind(ivar) / vol_nbor(j)
             end do
#endif
          end do

       end if

       ! Unlock all octs
       call unlock_cache(m,igrid)
       do j=1,nSKnei
          call unlock_cache(m,igrid_nbor(j))
       end do

       !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       ! Periodic box
       !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       do idim=1,ndim
          if(r%periodic(idim))then
             if(p%xp(ipart,idim)< 0.0d0           )p%xp(ipart,idim)=p%xp(ipart,idim)+r%box_size(idim)
             if(p%xp(ipart,idim)>=r%box_size(idim))p%xp(ipart,idim)=p%xp(ipart,idim)-r%box_size(idim)
          endif
       end do

       !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       ! Save sink data at a high cadence if needed
       !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       if(output_file .and. p%idp(ipart) < 100000)then
          if(r%agn)then
             call dump_sink_data_fine_AGN(s,p,ipart,dt_loc,scale_l,scale_t,scale_d, &
                  &                       dMBH_overdt,dMED_overdt,rho_inf,cs_gas, &
                  &                       fbk_mom_agn,fbk_ener_agn)
          else
             call dump_sink_data_fine(s,p,ipart,dt_loc,scale_l,scale_t,scale_d, &
                  &                   dMBH_overdt,dMED_overdt,rho_inf,cs_gas)
          end if
       end if

    end do ! End loop over ipart

    call close_cache(mdl)

    end associate
#endif
#endif
  end subroutine sink_mechanical
  !##############################################################################
  !##############################################################################
  !##############################################################################
  !##############################################################################
end module sink_mechanical_module
