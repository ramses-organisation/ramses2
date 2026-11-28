module cr_source_terms_module
#if defined(_CUDA) && defined(DO_CR)
  use gpu_runner, only: gpu_cr_source_terms
#endif
contains
!###########################################################
!###########################################################
!###########################################################
!###########################################################
recursive subroutine r_cr_source_terms(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::ilevel

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_CR_SOURCE_TERMS,pst%iUpper+1,input_size,0,ilevel)
     call r_cr_source_terms(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
#if defined(_CUDA) && defined(DO_CR)
     call gpu_cr_source_terms(pst%s,ilevel)
#else
     call cr_source_terms(pst%s,ilevel)
#endif
  endif

end subroutine r_cr_source_terms
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cr_source_terms(s,ilevel)
  use amr_parameters, only: ndim, nvector, twondim, twotondim
  use boundaries, only: init_bound_refine
  use cache
  use cache_commons, only: msg_large_realdp
#ifdef DO_CR
  use cr_flux_module, only: invrotatevec, rotatevec
#endif
  use cr_parameters, only: ncrgrp, smallecr
  use marshal, only: pack_fetch_refine,unpack_fetch_refine
  use nbors_utils
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t)::s
  integer::ilevel
#ifdef DO_CR
  integer::i,ind,igrid,idim,ngrid,nleaf,iGrp,i_nbor,icelld,icellg,icellp
  integer::igridd,igridg,igridp
  integer,dimension(1:nvector)::ind_leaf
  integer,dimension(1:3,1:8),save::iii=reshape(&
       & (/0,0,0,1,0,0,0,1,0,1,1,0,0,0,1,1,0,1,0,1,1,1,1,1/),(/3,8/))
  integer,dimension(1:3,1:6),save::shift=reshape(&
       & (/-1,0,0,1,0,0,0,-1,0,0,1,0,0,0,-1,0,0,1/),(/3,6/))
  integer(kind=8),dimension(0:ndim)::hash_key,hash_nbor
  integer,dimension(1:twondim)::igridn,icelln,refined
  type(msg_large_realdp)::dummy_large_realdp
  real(kind=8),dimension(1:nvector,1:ndim,1:ncrgrp),save::gradecr_loc,gradpcr_loc
  real(kind=8),dimension(1:3)::vs_loc,bloc

  real(kind=8)::dx, dx_g,dx_d
  real(kind=8)::pcrg,pcrd,sqrt3,fred,frotx,froty,frotz,sint,cost,sinp,cosp
  real(kind=8)::cr_c_two,dt, va_loc, bdotgradE_loc, norm,bxby, f1, f2, f3
  real(kind=8)::coef_11, coef_12, coef_13, coef_14, coef_21, coef_22
  real(kind=8)::coef_31, coef_33, coef_41, coef_44
  real(kind=8)::e_coef, new_ec, old_ec, sigma_x, sigma_y, sigma_z, sigma_stream
  real(kind=8)::Ecr, rhs2, rhs3, rhs4, v1, v2, v3
  real(kind=8)::g_grp, smallp, three_gmone
  real(kind=8)::scale_nH, scale_T2, scale_l, scale_d, scale_t, scale_v
  real(kind=8)::dtcool, lambda_cr, lambda_rad, nH
  real(kind=8),parameter:: zeta_cr=7.51d-16, fneut = 0.875d0, ne = 1.0d-3
  ! -------------------------------------------------------------------
  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)
  if(r%verbose.and.g%myid==1)write(*,'("   Entering cr_source_terms for level ",I2)')ilevel
  dx = r%boxlen/2**ilevel
  cr_c_two = g%cr_c(ilevel)**2
  smallp = r%smallc**2/r%gamma
  dt = g%dtnew(ilevel)
  if(r%cr_isotropic_pressure)then
     sqrt3=sqrt(3d0)
  else
     sqrt3=1.0d0
  endif

  hash_key(0)=ilevel+1

  call open_cache(mdl, m, pack_size=storage_size(dummy_large_realdp)/32, &
       pack=pack_fetch_refine, unpack=unpack_fetch_refine, &
       init=init_flush_cr_source_terms, flush=pack_flush_cr_source_terms, &
       combine=unpack_flush_cr_source_terms, bound=init_bound_refine)

  if(r%cr_cooling) then
      call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
      dtcool = g%dtnew(ilevel) * scale_t
      ! Fitz Axen + 24. This can easily be improved by using actual ne and neutral fraction
      lambda_cr=zeta_cr*(1d0 + 0.22d0*ne + 0.125d0*fneut)   ! [cm^3 s^-1]
      lambda_rad = zeta_cr*(1d0-0.667d0) ! Radiative cooling [cm^3 s^-1]
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

        ! Loop over leaf cells to compute their CR energy gradients for each group and in each dimension
        do i=1,nleaf
          ! Compute cell hash key
          hash_key(1:ndim)=2*m%grid(ind_leaf(i))%ckey(1:ndim)+iii(1:ndim,ind)

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
              refined(i_nbor) = 0
            else
              hash_nbor(0)=hash_nbor(0)-1
              hash_nbor(1:ndim)=hash_nbor(1:ndim)/2
              call get_parent_cell(s,hash_nbor,igridp,icellp,flush_cache=.false.,fetch_cache=.true.,lock=.true.)
              igridn(i_nbor)=igridp
              icelln(i_nbor)=icellp
              refined(i_nbor) = 1
            endif
          end do

         ! Loop over dimensions
          do idim=1,ndim
            ! Gather neighbouring e_cr values
            do igrp=1,ncrgrp
              icellg=icelln(2*idim-1)
              icelld=icelln(2*idim  )
              igridg=igridn(2*idim-1)
              igridd=igridn(2*idim  )
              pcrg   = max(m%uold(icellg,r%iEcr+igrp-1,igridg),smallecr)
              pcrd   = max(m%uold(icelld,r%iEcr+igrp-1,igridd),smallecr)
              dx_g   = dx+refined(2*idim-1)*dx*0.5
              dx_d   = dx+refined(2*idim)*dx*0.5
              gradecr_loc(i,idim,iGrp) = (pcrd-pcrg)/(dx_g+dx_d)
              gradpcr_loc(i,idim,iGrp) = gradecr_loc(i,idim,iGrp) * (r%gamma_rad(r%iEcr-r%inener+igrp)-1d0)
            end do
          end do

          do i_nbor=1,twondim
            call unlock_cache(m,igridn(i_nbor))
          end do

        end do ! loop over nleaf

        ! Loop again over leaf cells, updating their CR energies and fluxes
        do i=1,nleaf
          vs_loc(:)=0d0
          va_loc=0d0
          norm=0d0
          do idim=1,3
            bloc(idim) = 0.5 * ( m%bold(ind,idim,ind_leaf(i)) + m%bold(ind,idim+3,ind_leaf(i)) )
            vs_loc(idim) = bloc(idim)/sqrt(m%uold(ind,1,ind_leaf(i)))
            va_loc = va_loc + vs_loc(idim)**2
            norm = norm+bloc(idim)**2
          enddo
          va_loc = sqrt(va_loc)
          norm = max(sqrt(norm),1d-30)

          if(r%cr_v_alfven.gt.0d0) then
            vs_loc(:) = 0d0
            vs_loc(1) = r%cr_v_alfven
            va_loc = r%cr_v_alfven
          endif

         bxby = sqrt(bloc(1)**2+bloc(2)**2)
          if(norm.gt.1e-10) then
            sint = bxby/norm
            cost = bloc(3)/norm
          else
            sint = 1d0
            cost = 0d0
          endif
          if(bxby.gt.1e-10) then
            sinp = bloc(2)/bxby
            cosp = bloc(1)/bxby
          else
            sinp = 0d0
            cosp = 1d0
          endif
          do igrp=1,ncrgrp
            g_grp = r%gamma_rad(r%iEcr-r%inener+igrp)
            three_gmone = 3d0*(g_grp - 1d0)
            bdotgradE_loc=0.
            ! bdotgrade is needed for eq 3 in Jiang & Oh
            do idim=1,ndim
              bdotgradE_loc = bdotgradE_loc + bloc(idim)*gradecr_loc(i,idim,iGrp)
            enddo
            ! Local streaming velocity, needed for eq 18 in Jiang & Oh
            vs_loc(:) = -vs_loc(:) * bdotgradE_loc / max(1d-50,abs(bdotgradE_loc)) !! eq 3 in PO17

            Ecr = max(m%unew(ind,r%iEcr+igrp-1,ind_leaf(i)),smallecr) ! Ecr

            ! Flux update ... first rotate flux vector onto B-field,
            ! i.e. describing the flux in the B coordinate system
            f2=0. ; f3=0.
            f1 = m%crunew(ind,1+(igrp-1)*ndim,ind_leaf(i))
#if NDIM>1
            f2 = m%crunew(ind,2+(igrp-1)*ndim,ind_leaf(i))
#endif
#if NDIM>2
            f3 = m%crunew(ind,3+(igrp-1)*ndim,ind_leaf(i))
#endif

            ! Make sure that |F|<=cE/sqrt3 (sqrt3=1 or sqrt(3) for M1 or P1 resp.)
            fred = sqrt(f1**2+f2**2+f3**2)/(g%cr_c(ilevel)*Ecr)*sqrt3
            if(fred .gt. 1d0) then
               f1 = f1/fred ; f2 = f2/fred ; f3 = f3/fred
            endif

            frotx=f1; froty=f2; frotz=f3
            call rotatevec(sint, cost, sinp, cosp, frotx, froty, frotz)
            v1=0. ; v2=0. ; v3=0.
            if(r%cr_streaming_heating) then
              v1 = vs_loc(1) ; v2 = vs_loc(2) ; v3 = vs_loc(3)
            endif

            ! Rotate streaming velocity
            call rotatevec(sint, cost, sinp, cosp, v1, v2, v3)

            rhs2 = frotx ; rhs3 = froty ; rhs4 = frotz

            if(r%cr_streaming_diffusion) then
               sigma_stream = max(1./r%cr_dmax_code &
                  & ,abs(bdotgradE_loc)/3d0 / norm / va_loc / g_grp / Ecr)
               sigma_x = 1./(r%cr_d_code(igrp) + 1./sigma_stream)
            else
               sigma_x = 1./r%cr_d_code(igrp)
            endif
            sigma_y = 1./(r%cr_d_code(iGrp)*r%cr_d_perp_factors(iGrp))
            sigma_z = 1./(r%cr_d_code(iGrp)*r%cr_d_perp_factors(iGrp))

            coef_11 = 1.0
            coef_12 = dt * sigma_x * v1 *three_gmone
            coef_13 = dt * sigma_y * v2 *three_gmone
            coef_14 = dt * sigma_z * v3 *three_gmone

            coef_21 = 0.0
            coef_22 = 1.0 + dt * sigma_x * cr_c_two

            coef_31 = 0.0
            coef_33 = 1.0 + dt * sigma_y * cr_c_two

            coef_41 = 0.0
            coef_44 = 1.0 + dt * sigma_z * cr_c_two

            ! newfr1 = (rhs2 - coef21 * newEc)/coef22
            ! newfr2= (rhs3 - coef31 * newEc)/coef33
            ! newfr3 = (rhs4 - coef41 * newEc)/coef44
            ! coef11 - coef21 * coef12 /coef22 - coef13 * coef31 /coef33 - coef41 * coef14 /coef44)* newec
            !    =Ecr - coef12 *rhs2/coef22 - coef13 * rhs3/coef33 - coef14 * rhs4/coef44

            e_coef = coef_11 - coef_12 * coef_21/coef_22 - coef_13 * coef_31/coef_33 &
                            - coef_14 * coef_41/coef_44
            new_ec = Ecr - coef_12 * rhs2/coef_22 - coef_13 * rhs3/coef_33 &
                         - coef_14 * rhs4/coef_44
            new_ec = new_ec / e_coef

            old_ec = m%unew(ind,r%iEcr+igrp-1,ind_leaf(i))
            m%unew(ind,r%iEcr+igrp-1,ind_leaf(i)) = &
              & m%unew(ind,r%iEcr+igrp-1,ind_leaf(i)) + (new_ec-old_ec)

            ! Floor the CR energy and update total energy if necessary
            m%unew(ind,r%iEcr+igrp-1,ind_leaf(i)) = max(m%unew(ind,r%iEcr+igrp-1,ind_leaf(i)), smallecr)

            frotx = (rhs2 - coef_21 * new_ec)/coef_22
            froty = (rhs3 - coef_31 * new_ec)/coef_33
            frotz = (rhs4 - coef_41 * new_ec)/coef_44

            ! Make sure that |F|<=cE/sqrt3 (sqrt3=1 or sqrt(3) for M1 or P1 resp.)
            fred = sqrt(frotx**2+froty**2+frotz**2)/(g%cr_c(ilevel)*new_ec)*sqrt3
            if(fred .gt. 1d0) then
               frotx = frotx/fred ; froty = froty/fred ; frotz = frotz/fred
            endif

            ! Rotate the flux back to the simulation coordinate system
            call invrotatevec(sint, cost, sinp, cosp, frotx, froty, frotz)
            m%crunew(ind,1+(igrp-1)*ndim,ind_leaf(i)) = frotx
#if NDIM>1
            m%crunew(ind,2+(igrp-1)*ndim,ind_leaf(i)) = froty
#endif
#if NDIM>2
            m%crunew(ind,3+(igrp-1)*ndim,ind_leaf(i)) = frotz
#endif

            ! Add Coulomb+hadronic cooling losses from Fitz Axen + 24,
            ! using exponential decay of the CR energy and flux at the rate lambda_cr * n_H.
            if(r%cr_cooling) then
              nH = m%unew(ind,1,ind_leaf(i)) * scale_nH
              old_ec = m%unew(ind,r%iEcr+igrp-1,ind_leaf(i))
              ! Radiative losses need to be subtracted from total energy
              m%unew(ind,5,ind_leaf(i)) = &
                & m%unew(ind,5,ind_leaf(i)) - old_ec * (1d0 - EXP(-lambda_rad * nH * dtcool))
              ! Update E_cr and make sure it doesn't go below smallecr
              m%unew(ind,r%iEcr+igrp-1,ind_leaf(i)) = old_ec * EXP(-lambda_cr * nH * dtcool)
              m%unew(ind,r%iEcr+igrp-1,ind_leaf(i)) = &
                & max(m%unew(ind,r%iEcr+igrp-1,ind_leaf(i)), smallecr)
              do idim=1,ndim
                ! CR 'flux' cooling:
                m%crunew(ind,idim+(igrp-1)*ndim,ind_leaf(i)) = &
                  & m%crunew(ind,idim+(igrp-1)*ndim,ind_leaf(i)) &
                  & *EXP(-lambda_cr * nH * dtcool * r%cr_c_fraction**2)
              end do
            endif

          end do ! End loop over groups
        end do ! End loop over leaves
     end do ! End loop over octs
  end do ! End loop over cells

  call close_cache(mdl)

  end associate
#endif
end subroutine cr_source_terms
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine init_flush_cr_source_terms(mesh,igrid,hash_key)
  use amr_commons, only: mesh_t
  use amr_parameters, only: ndim, twotondim
  use cache_commons, only: msg_large_realdp
  use cr_parameters, only: ncruvar
  use hydro_parameters, only: nvar
  type(mesh_t)::mesh
  integer::igrid
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind,ivar

  mesh%grid(igrid)%lev=hash_key(0)
  mesh%grid(igrid)%ckey(1:ndim)=hash_key(1:ndim)
#ifdef HYDRO
  do ivar=1,nvar ! Joki: not sure this is needed. hydro vars are updated though...
     do ind=1,twotondim
        mesh%unew(ind,ivar,igrid)=0.0d0
     enddo
  enddo
#endif
#ifdef DO_CR
  do ivar=1,ncruvar
     do ind=1,twotondim
        mesh%crunew(ind,ivar,igrid)=0.0d0
     enddo
  enddo
#endif

end subroutine init_flush_cr_source_terms
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine pack_flush_cr_source_terms(mesh,igrid,msg_size,msg_array)
  use amr_commons, only: mesh_t
  use amr_parameters, only: twotondim
  use cache_commons, only: msg_large_realdp
  use cr_parameters, only: ncruvar
  use hydro_parameters, only: nvar
  type(mesh_t)::mesh
  integer::igrid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind,ivar
  type(msg_large_realdp)::msg

#ifdef HYDRO
  do ivar=1,nvar
     do ind=1,twotondim
        msg%realdp_hydro(ind,ivar)=mesh%unew(ind,ivar,igrid)
     end do
  end do
#endif

  do ivar=1,ncruvar
     do ind=1,twotondim
#ifdef DO_CR
        msg%realdp_cr(ind,ivar)=mesh%crunew(ind,ivar,igrid)
#endif
     end do
  end do

  msg_array=transfer(msg,msg_array)

end subroutine pack_flush_cr_source_terms
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine unpack_flush_cr_source_terms(mesh,igrid,msg_size,msg_array,hash_key)
  use amr_commons, only: mesh_t
  use amr_parameters, only: ndim, twotondim
  use cache_commons, only: msg_large_realdp
  use cr_parameters, only: ncruvar
  use hydro_parameters, only: nvar
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
  do ivar=1,nvar
     do ind=1,twotondim
        mesh%unew(ind,ivar,igrid)=mesh%unew(ind,ivar,igrid)+msg%realdp_hydro(ind,ivar)
     end do
  end do
#endif

  do ivar=1,ncruvar
     do ind=1,twotondim
#ifdef DO_CR
        mesh%crunew(ind,ivar,igrid)=mesh%crunew(ind,ivar,igrid)+msg%realdp_cr(ind,ivar)
#endif
     end do
  end do

end subroutine unpack_flush_cr_source_terms
!###########################################################
!###########################################################
!###########################################################
!###########################################################
end module cr_source_terms_module
