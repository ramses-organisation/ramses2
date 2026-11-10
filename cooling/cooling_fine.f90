module cooling_fine_module
  use cooling_module, only: set_table, cooling_needs_update
#ifdef _CUDA
  use gpu_runner, only: gpu_cooling, gpu_check_nan
  use cooling_device, only: gpu_upload_cooling_table
#ifndef DO_RTZ
  use neq_cooling_device, only: gpu_neq_cooling
#endif
#elif defined(_METAL)
  use metal_runner, only: metal_cooling, metal_upload_cooling_table
#endif
contains
!###########################################################
!###########################################################
!###########################################################
!###########################################################
recursive subroutine r_cooling_fine(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::ilevel

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_COOLING_FINE,pst%iUpper+1,input_size,0,ilevel)
     call r_cooling_fine(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
#ifdef _CUDA
#ifndef DO_RTZ
     if(pst%s%r%neq_chem)then
        call gpu_neq_cooling(pst%s, ilevel)
     else
        call gpu_cooling(pst%s, ilevel)
     endif
#else
     call gpu_cooling(pst%s, ilevel)
#endif
     call gpu_check_nan(pst%s, ilevel, 'cooling')
#elif defined(_METAL)
     call metal_cooling(pst%s, ilevel)
#else
     call cooling_fine(pst%s%r,pst%s%g,pst%s%m,pst%s%cool,pst%s%tables,ilevel)
#endif

     ! Compute new cooling table for cosmo runs
#ifndef DO_RTZ
     if(pst%s%r%cooling.and.ilevel==pst%s%r%levelmin.and.pst%s%r%cosmo)then
        if (cooling_needs_update(pst%s%cool, dble(pst%s%g%aexp))) then
           if(pst%s%g%myid==1)write(*,*)'Computing new cooling table'
           call set_table(pst%s%cool,dble(pst%s%g%aexp))
#ifdef _CUDA
           call gpu_upload_cooling_table(pst%s%cool)
#elif defined(_METAL)
           call metal_upload_cooling_table(pst%s%cool)
#endif
        endif
     endif
#endif

  endif

end subroutine r_cooling_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cooling_fine(r,g,m,c,tables,ilevel)
  use constants
  use amr_parameters, only: ndim, nvector, twotondim
  use hydro_parameters, only: nener, nion
  use rt_parameters, only: nrtgrp, smallnp
  use amr_commons, only: run_t, global_t, mesh_t
  use cooling_module, only: cooling_t, solve_cooling, T2_min_fix
  use coolrates_module, only: neq_cooling_t
#ifdef DO_RTZ
  use rtz_cooling_module, only: rtz_solve_cooling
  use rtz_module, only: n_elements, elements
#else
  use neq_cooling_module, only: neq_solve_cooling
#endif
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(cooling_t)::c
  type(neq_cooling_t)::tables
  integer::ilevel
  !-------------------------------------------------------------------
  ! Compute cooling for leaf cells at level ilevel
  !-------------------------------------------------------------------
  integer::i,ii,jj,ind,igrid,idim,ngrid,nleaf,counter,e_counter
  real(kind=8)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(kind=8)::dtcool,nH_eos,nCOM,damp_factor,cooling_switch,t_blast
  integer,dimension(1:nvector)::ind_leaf
  real(kind=8),dimension(1:nvector)::nH,T2,delta_T2,ekk,err,emag
  real(kind=8),dimension(1:nvector)::T2min,Zsolar,boost
  real(kind=8)::factor1, factor2, factor3
  !  logical,dimension(1:nvector)::cooling_on=.true.
#ifdef DO_RTZ
  real(kind=8),dimension(1:n_elements, 1:n_elements, 1:nvector):: xion
#else
  real(kind=8),dimension(nion, 1:nvector):: xion
#endif
#ifdef DO_RT
  integer::ig,iNp
  real(kind=8),dimension(1:ndim)::Fpnew
  real(kind=8),dimension(nrtgrp, 1:nvector):: Np, dNpdt=0d0
  real(kind=8),dimension(ndim, nrtgrp, 1:nvector):: Fp, dFpdt=0
  real(kind=8),dimension(ndim, 1:nvector):: p_gas
  real(kind=8)::scale_Np,scale_Fp,Npnew
#endif
#ifdef DO_RTZ
  real(kind=8), dimension(n_elements, 1:nvector):: nElement
  real(kind=8):: dx_SS_H2
#endif
#if NENER>0
  integer::irad
#endif

#ifdef HYDRO
  if(r%verbose.and.g%myid==1)write(*,'("   Entering cooling_fine for level ",I2)')ilevel

  ! Conversion factor from user units to cgs units
  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
#ifdef DO_RT
  call rt_units(r,g,scale_Np,scale_Fp)
#endif
  ! Density for isothermal/polytropic EOS in H/cc
  nH_eos = r%eos_nH

  ! This is used only for cosmological runs
  if(r%cosmo)then
     nCOM = 200.0*g%omega_b*rhoc*(g%h0/100)**2/g%aexp**3*c%X/mH
     nH_eos = MAX(nCOM,nH_eos)
  endif

  ! Loop over cells
  do ind=1,twotondim
     ! Loop over octs with vector sweeps
     do igrid=m%head(ilevel),m%tail(ilevel),nvector

        ! Collect a vector of leaf cells
        ngrid=MIN(nvector,m%tail(ilevel)-igrid+1)
        nleaf=0
        do i=1,ngrid
           if(.NOT. m%grid(igrid+i-1)%refined(ind))then
              nleaf=nleaf+1
              ind_leaf(nleaf)=igrid+i-1
           endif
        end do
        if(nleaf.eq.0)cycle

        ! Compute rho
        do i=1,nleaf
           nH(i)=MAX(dble(m%uold(ind,1,ind_leaf(i))),r%smallr)
        end do

        ! Compute metallicity in solar units
        if(r%metal)then
           do i=1,nleaf
              Zsolar(i)=m%uold(ind,r%imetal,ind_leaf(i))/nH(i)/0.02d0
           end do
        else
           do i=1,nleaf
              Zsolar(i)=r%z_ave
           end do
        endif

        ! Total energy
        do i=1,nleaf
           T2(i)=m%uold(ind,5,ind_leaf(i))
        end do

        ! Kinetic energy
        do i=1,nleaf
           ekk(i)=0.0d0
        end do
        do idim=1,3
           do i=1,nleaf
              ekk(i)=ekk(i)+0.5d0*m%uold(ind,idim+1,ind_leaf(i))**2/nH(i)
           end do
        end do

        ! Non-thermal energies
        do i=1,nleaf
           err(i)=0.0d0
        end do
#if NENER>0
        do irad=1,nener
           do i=1,nleaf
              err(i)=err(i)+m%uold(ind,5+irad,ind_leaf(i))
           end do
        end do
#endif
        ! Magnetic energy
        do i=1,nleaf
           emag(i)=0.0d0
        end do
#ifdef MHD
        do idim=1,3
           do i=1,nleaf
              emag(i)=emag(i)+0.125d0*(m%bold(ind,idim,ind_leaf(i)) &
                   &                  +m%bold(ind,idim+3,ind_leaf(i)))**2
           end do
        end do
#endif
        ! Gas thermal pressure
        do i=1,nleaf
           T2(i)=(r%gamma-1.0d0)*(T2(i)-ekk(i)-err(i)-emag(i))
        end do

        ! Compute T2=T/mu in Kelvin
        do i=1,nleaf
           T2(i)=T2(i)/nH(i)*scale_T2
        end do

        ! Compute nH in H/cc
        do i=1,nleaf
           nH(i)=nH(i)*scale_nH
        end do

        ! Compute radiation self-shielding factor
        if(r%self_shielding)then
           do i=1,nleaf
              boost(i)=MAX(exp(-nH(i)/0.01d0),1.0D-20)
           end do
        else
           do i=1,nleaf
              boost(i)=1
           end do
        endif

        ! Compute ionization fraction
#ifdef DO_RTZ
        counter = 0
        e_counter = 0
        do ii=1,n_elements ! loop over elements
           if (elements(ii)%atomic_number.gt.0) then
              do jj=1,elements(ii)%n_ions ! loop over ions
                 do i=1,nleaf !loop over leaf cells
                    xion(ii,jj,i) = m%uold(ind,r%iIons+counter,ind_leaf(i)) &
                                   /m%uold(ind,1,ind_leaf(i))
                    if (jj.eq.1) then
                       ! This gives us a number density [Atoms/cm^3]
                       nElement(ii,i) = m%uold(ind,r%ichem+e_counter,ind_leaf(i)) &
                                       *scale_nH / elements(ii)%atomic_mass
                    end if
                 end do ! end loop over leaf cells
                 counter = counter + 1 ! increment ionization counter
              end do ! end loop over ions
              e_counter = e_counter + 1 ! increment element counter
           end if
        end do ! end loop over elements

        ! deal with molecules separately
        if (elements(1)%atomic_number.gt.0 .and. r%isH2_rtz) then
           do i=1,nleaf !loop over leaf cells
              xion(1,3,i) = m%uold(ind,r%iIons+counter,ind_leaf(i)) &
                           /m%uold(ind,1,ind_leaf(i))
           end do ! end loop over leaf cells
           counter = counter + 1
        endif
#else
        if(r%neq_chem) then
           do ii=0,nIon-1
              do i=1,nleaf
                 xion(1+ii,i) = m%uold(ind,r%iIons+ii,ind_leaf(i)) &
                      &        /m%uold(ind,1,ind_leaf(i))
              end do
           end do
        endif
#endif

        ! Get photon densities and flux magnitudes
#ifdef DO_RT
        do ig=1,nrtgrp
           iNp=1+(ig-1)*(ndim+1)
           do i=1,nleaf
              Np(ig,i)          = scale_Np * m%rtuold(ind,iNp,ind_leaf(i))
              Fp(1:ndim, ig, i) = scale_Fp * m%rtuold(ind,iNp+1:iNp+ndim,ind_leaf(i))
           enddo
        end do
#endif

        ! Compute gas momentum for radiation force
#ifdef DO_RT
        do i=1,nleaf
           p_gas(1:ndim,i) = m%uold(ind,2:1+ndim,ind_leaf(i)) * scale_d * scale_v
        end do
#endif
        !==========================================
        ! Compute minimum temperature  for cooling
        ! floor or for isothermal/polytropic EOS.
        ! You can put your own polytrope EOS here.
        !==========================================
        if(r%eos_type==1)then ! Strictly isothermal
           do i=1,nleaf
              T2min(i) = r%eos_T2
           end do
        else if(r%eos_type==2)then ! Power law
           do i=1,nleaf
              T2min(i) = r%eos_T2*(nH(i)/nH_eos)**(r%eos_index-1.0d0)
           end do
        else if(r%eos_type==3)then ! Smooth transition from isothermal to power law
           do i=1,nleaf
              T2min(i) = r%eos_T2*(1+(nH(i)/nH_eos)**(r%eos_index-1.0d0))
           end do
        else if(r%eos_type==4)then ! Sharp transition from isothermal to power law
           do i=1,nleaf
              if(nH(i)<nH_eos)then
                 T2min(i) = r%eos_T2
              else
                 T2min(i) = r%eos_T2*(nH(i)/nH_eos)**(r%eos_index-1.0d0)
              endif
           end do
        else if(r%eos_type==5)then ! Second collapse
           do i=1,nleaf
            factor1 = sqrt(1 + (nH(i)/(3.866301516d-15*c%X/mH))**(2*0.4))
            factor2 = (1 + (nH(i)/(3.866301516d-10*c%X/mH)))**(-0.3)
            factor3 = (1 + (nH(i)/(3.866301516d-5*c%X/mH)))**0.56667
            T2min(i) = r%eos_T2 * factor1 * factor2 * factor3
           end do
        endif

        ! Compute thermal temperature by subtracting the polytrope
        if(r%cooling.or.r%cooling_ism.or.r%neq_chem)then
           do i=1,nleaf
              T2(i) = min(max(T2(i)-T2min(i),T2_min_fix),dble(r%T2max))
           end do
        endif

        ! Compute cooling time step in second
        dtcool = g%dtnew(ilevel)*scale_t

#ifdef DO_RT
        ! Isotropic emission (star and sink particles) in cgs units
        do ig=1,nrtgrp
           do i=1,nleaf
              dNpdt(ig,i) = m%emissivity(ind,ig,ind_leaf(i)) * scale_Np / scale_t
              dFpdt(:,ig,i) = 0
           end do
        end do

        ! Smooth RT update
        if(r%rt_smooth) then
           do ig=1,nrtgrp
              iNp=1+(ig-1)*(ndim+1)
              do i=1,nleaf ! Calc addition per sec to Np, Fp for current dt
                 Npnew = scale_Np * m%rtunew(ind,iNp,ind_leaf(i))
                 Fpnew = scale_Fp * m%rtunew(ind,iNp+1:iNp+ndim,ind_leaf(i))
                 dNpdt(ig,i) = dNpdt(ig,i) + (Npnew - Np(ig,i)) / dtcool
                 dFpdt(:,ig,i) = dFpdt(:,ig,i) + (Fpnew - Fp(:,ig,i)) / dtcool
              end do
           end do
        end if
#endif
        !==========================================
        ! Compute net cooling at constant nH
        !==========================================
        if(r%cooling)then
           ! Use classical ramses cooling
           call solve_cooling(c,nH,T2,Zsolar,boost,dtcool,delta_T2,nleaf)
        else if(r%cooling_ism)then
           ! Use cooling from cooling_module_frig described in Audit & Hennebelle 2005
           call solve_cooling_ism(nH,T2,dtcool,delta_T2,r%gamma,r%mu_mol,nleaf)
#ifdef DO_RTZ
        else if(r%neq_chem.and.r%rtz_cooling) then
           ! If both non-equilibrium chemistry and rtz_cooling are turned on
           ! we use a detailed model for the chemistry

           if (r%rtz_equilibrium_test.eq.2) then
               ! for now, just assume some density solar metallicity
               nElement(1:n_elements,1:nleaf)  = 0.d0  ! Initialize to zero
               nElement(1,1:nleaf)  = 1.d-1                          ! Hydrogen
               nElement(2,1:nleaf)  = nElement(1,1:nleaf) * 8.51d-02 ! Helium
               nElement(6,1:nleaf)  = nElement(1,1:nleaf) * 2.69d-04 ! Carbon
               nElement(7,1:nleaf)  = nElement(1,1:nleaf) * 6.76d-05 ! Nitrogen
               nElement(8,1:nleaf)  = nElement(1,1:nleaf) * 4.90d-04 ! Oxygen
               nElement(10,1:nleaf) = nElement(1,1:nleaf) * 8.51d-05 ! Neon
               nElement(12,1:nleaf) = nElement(1,1:nleaf) * 3.98d-05 ! Magnesium
               nElement(14,1:nleaf) = nElement(1,1:nleaf) * 3.24d-05 ! Silicon
               nElement(16,1:nleaf) = nElement(1,1:nleaf) * 1.32d-05 ! Sulfur
               nElement(26,1:nleaf) = nElement(1,1:nleaf) * 3.16d-05 ! Iron
           end if

           ! Compute the cell length in cm if needed
           dx_SS_H2 = 0.d0
           if (r%isH2_rtz) then
              dx_SS_H2 = r%boxlen/2**ilevel * scale_l
           endif

           call rtz_solve_cooling(r, tables, T2, g%aexp, xion, nElement, &
#ifdef DO_RT
                & Np, Fp, p_gas, dNpdt, dFpdt, ilevel, &
#endif
                & dtcool, nleaf, dx_SS_H2)
#else
        else if(r%neq_chem)then
           call neq_solve_cooling(r, tables, T2, xion, nH, Zsolar, &
#ifdef DO_RT
                & Np, Fp, p_gas, dNpdt, dFpdt, ilevel, &
#endif
                & dtcool, nleaf)
#endif
        endif

        ! Compute rho in code units
        do i=1,nleaf
           nH(i) = nH(i)/scale_nH
        end do

        ! Update fluid momentum and kinetic energy due to radiation force
#ifdef DO_RT
        do i=1,nleaf
           m%uold(ind,2:1+ndim,ind_leaf(i)) = p_gas(1:ndim,i) / scale_d / scale_v
        end do
        do i=1,nleaf
           ekk(i)=0.0d0
        end do
        do idim=1,3
           do i=1,nleaf
              ekk(i)=ekk(i)+0.5d0*m%uold(ind,idim+1,ind_leaf(i))**2/nH(i)
           end do
        end do
#endif
        ! Compute radiated energy in code units
        if(r%cooling.or.r%cooling_ism)then
           ! Compute net energy sink
           do i=1,nleaf
              delta_T2(i) = delta_T2(i)*nH(i)/scale_T2/(r%gamma-1.0d0)
           end do
        endif

        ! Compute fluid internal energy in code units
        if(r%cooling.or.r%cooling_ism.or.r%neq_chem)then
           do i=1,nleaf
              T2(i) = T2(i)*nH(i)/scale_T2/(r%gamma-1.0d0)
           end do
        endif

        ! Compute polytrope internal energy in code units
        do i=1,nleaf
           T2min(i) = T2min(i)*nH(i)/scale_T2/(r%gamma-1.0d0)
        end do

        ! Update fluid internal energy
        if(r%cooling.or.r%cooling_ism)then
           do i=1,nleaf
              T2(i) = T2(i) + delta_T2(i)
           end do
        endif

        ! Update ionization fraction
#ifdef DO_RTZ
        counter = 0
        do ii=1,n_elements ! loop over elements
           if (elements(ii)%atomic_number.gt.0) then
              do jj=1,elements(ii)%n_ions ! loop over ions
                 do i=1,nleaf !loop over leaf cells
                    m%uold(ind,r%iIons+counter,ind_leaf(i)) = xion(ii,jj,i)*nH(i)
                 end do ! end loop over leaf cells
                 counter = counter + 1
              end do ! end loop over ions
           end if
        end do ! end loop over elements

        ! deal with molecules separately
        if (elements(1)%atomic_number.gt.0 .and. r%isH2_rtz) then
           do i=1,nleaf !loop over leaf cells
              m%uold(ind,r%iIons+counter,ind_leaf(i)) = xion(1,3,i)*nH(i)
           end do ! end loop over leaf cells
           counter = counter + 1
        endif
#else
        if(r%neq_chem) then
           do ii=0,nion-1
              do i=1,nleaf
                 m%uold(ind,r%iIons+ii,ind_leaf(i)) = xion(1+ii,i)*nH(i)
              end do
           end do
        endif
#endif

        ! Update entropy if dual energy scheme is activated
        if(r%entropy.and.r%dual_energy.GE.0)then
           if(r%isothermal)then ! use only polytrope energy
              do i=1,nleaf
                 m%uold(ind,r%ientropy,ind_leaf(i)) = (T2min(i))/nH(i)**(r%gamma-1)*(r%gamma-1)
              end do
           else if(r%cooling.or.r%cooling_ism.or.r%neq_chem)then
              do i=1,nleaf
                 m%uold(ind,r%ientropy,ind_leaf(i)) = (T2(i) + T2min(i))/nH(i)**(r%gamma-1)*(r%gamma-1)
              end do
           endif
        endif

        ! Update total fluid energy
        if(r%isothermal)then ! use only polytrope energy
           do i=1,nleaf
              m%uold(ind,5,ind_leaf(i)) = T2min(i) + ekk(i) + err(i) + emag(i)
           end do
        else if(r%cooling.or.r%cooling_ism.or.r%neq_chem)then ! add polytrope to thermal energy
           do i=1,nleaf
              m%uold(ind,5,ind_leaf(i)) = T2(i) + T2min(i) + ekk(i) + err(i) + emag(i)
           end do
        endif

        ! Update photon densities and flux magnitudes
#ifdef DO_RT
        do ig=1,nrtgrp
           iNp=1+(ig-1)*(ndim+1)
           do i=1,nleaf
              m%rtuold(ind,iNp,ind_leaf(i)) = max(Np(ig,i)/scale_Np,smallNp)
              m%rtuold(ind,iNp+1:iNp+ndim,ind_leaf(i)) = Fp(1:ndim,ig,i)/scale_Fp
           enddo
        end do
#endif
     end do
     ! End loop over grid
  end do
  ! End loop over cells
#endif

end subroutine cooling_fine

end module cooling_fine_module
