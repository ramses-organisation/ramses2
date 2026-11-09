!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_rt_init_flow_fine(pst,ilevel)
  use amr_commons, only: ndim
  use amr_parameters, only: nvector, twotondim, flen
  use ramses_commons, only: pst_t
  use rt_input_condinit_module, only: r_rt_input_condinit
  use rt_parameters, only: nrtgrp
  implicit none
  type(pst_t)::pst
  integer::ilevel
  !--------------------------------------------------------------------
  ! This routine is the master procedure to input a given initial
  ! condition for RT variables into the exisiting AMR structure from an
  ! analytical model given as an external routine.
  !--------------------------------------------------------------------
  character(len=flen)::filename
  logical::ok_file1,ok_file2,ok_file

  associate(s=>pst%s,r=>pst%s%r,m=>pst%s%m)

  if(s%m%noct_tot(ilevel)==0)return
  if(s%r%verbose)write(*,111)ilevel
111 format(' Entering init_rt_fine for level ',I2)

  ! Use internal-defined or user-defined functions.
  ! We always call the condinit routine, even in cosmological simulations,
  ! just to initialisethe RT variables (to small values)
  if(s%r%verbose)write(*,*)'Initialising RT variables'
  call r_rt_input_condinit(pst,ilevel,1)

  filename=TRIM(s%r%initfile(ilevel))//'/ic_d'
  INQUIRE(file=filename,exist=ok_file1)
  filename=TRIM(s%r%initfile(ilevel))//'/ic_deltab'
  INQUIRE(file=filename,exist=ok_file2)
  ok_file=ok_file1.or.ok_file2
  if(ok_file)then
     ! No initialization necessary for photons
  else if (s%r%filetype=='gadget')then
     ! No initialization necessary for photons
  endif

  end associate

end subroutine m_rt_init_flow_fine
!###############################################
!###############################################
!###############################################
!###############################################
module update_rt_c_module
#if defined(_CUDA) && defined(DO_RT)
  use gpu_manager, only: gpu_star_to_host
#endif
  type :: in_broadcast_rt_c
     real(kind=8)::rt_c
  end type in_broadcast_rt_c
contains
!##############################################################
!##############################################################
!##############################################################
!##############################################################
recursive subroutine r_rt_neq_updates(pst, nstep_coarse, input_size)
  use mdl_module
  use coolrates_module, only: update_rt_c, update_coolrates_tables
#ifdef DO_RTZ
  use rtz_cooling_module, only: rtz_updateRTGroups_CoolConstants
#else
  use neq_cooling_module, only: updateRTGroups_CoolConstants, update_metal_cooling
#ifdef _CUDA
  use neq_cooling_device, only: gpu_upload_neq_tables
#endif
#endif
  use ramses_commons, only: pst_t
  use output_rt_module, only: output_photon_emission_stats
  use SED_module, only: update_SED_group_props
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::nstep_coarse, rID

  associate(s=>pst%s,r=>pst%s%r,m=>pst%s%m,g=>pst%s%g)

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_RT_NEQ_UPDATES,pst%iUpper+1,input_size,0,nstep_coarse)
     call r_rt_neq_updates(pst%pLower, nstep_coarse, input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     ! Update reduced speed of light
     if(r%rt .and. r%cosmo) call update_rt_c(r, g, pst%s%tables)

     ! Update Compton heating
     if(r%neq_chem .and. r%cosmo) &
        call update_coolrates_tables(r, pst%s%tables, dble(g%aexp))

     if(r%rt .and. r%star .and. r%rt_advect .and. r%rt_star &
        .and. r%sedprops_update>0)then
        ! The interval is tested first: mod(n,0) is undefined
        if(mod(nstep_coarse,r%sedprops_update)==0)then
#if defined(_CUDA) && defined(DO_RT)
           call gpu_star_to_host(s)
#endif
           ! Update cross sections based on evolving star properties
           call update_SED_group_props(r, g, s%SED, s%star)
        endif
     endif

     ! Update radiation heating and cooling constants
     if(r%cosmo.or.(r%rt))then
#ifdef DO_RTZ
        call rtz_updateRTGroups_CoolConstants(r, s%tables)
#else
        call updateRTGroups_CoolConstants(r, s%tables)
#endif
     endif

     ! Update UV background constants for metal cooling
#ifndef DO_RTZ
     if(r%neq_chem .and. r%cosmo) &
        call update_metal_cooling(s%tables, dble(g%aexp))
#endif

     if(r%rt_advect .and. r%rt_emission_stats .and. r%rt_star) &
        call output_photon_emission_stats(r,g)

#ifdef _CUDA
#ifndef DO_RTZ
     ! Copy the updated tables to the GPU chemistry
     if(r%neq_chem) call gpu_upload_neq_tables(s)
#endif
#endif

  endif

  end associate

end subroutine r_rt_neq_updates
!##############################################################
!##############################################################
!##############################################################
!##############################################################
end module update_rt_c_module
