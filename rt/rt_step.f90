module rt_step_module
contains
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine m_rt_step(pst,ilevel)
  use ramses_commons, only: pst_t
  use cooling_fine_module, only: r_cooling_fine
  use rt_godunov_fine_module, only: r_rt_godunov_fine,r_set_rtunew,r_set_rtuold
  use rt_upload_module, only: m_rt_upload_fine
  use rt_input_condinit_module, only: r_rt_input_source_regions
  use cooling_fine_module, only: r_cooling_fine
  use rt_star_feedback, only: r_star_RT_feedback
  type(pst_t)::pst
  integer::ilevel

  real(kind=8) :: dt_save, t_save
  real(kind=8) :: dt_rad, t_rad
  integer  :: i, i_substep

  associate(r=>pst%s%r,g=>pst%s%g,m=>pst%s%m,mdl=>pst%s%mdl)

  ! Store hydro timestep length and hydro time
  dt_save = g%dtnew(ilevel)
  t_save = g%t

  ! We shift the time backwards one hydro-dt, to get evolution of stellar
  ! ages within the hydro timestep, in the case of rt subcycling:
  t_rad = t_save - dt_save

  ! Get RT courant time step at coarse level
  call get_rt_courant_dt(r,g,dt_rad,ilevel)
  ! Compute RT courant time step at fine level
  dt_rad = MIN(dt_rad,dt_save)

  ! Compute number of subcycles
  i_substep = ceiling(dt_save/dt_rad-1e-12)
  dt_rad = dt_save/dble(i_substep)

  ! Compute isotropic emissivity from stars
  if(r%star .and. r%rt_star) call r_star_RT_feedback(pst,ilevel,1)

  ! Compute isotropic emissivity from sinks
  ! if(r%rt_sink) call r_sink_RT_feedback(pst,ilevel,1)

  ! RT sub-cycle loop
  do i = 1, i_substep

     ! Shift the RT time forwards one dt_rad
     t_rad = t_rad + dt_rad

     ! Update time variables
     call m_rt_update_time(pst,ilevel,t_rad,dt_rad)

     ! Set rtunew equal to rtuold
     if(i>1)call r_set_rtunew(pst,ilevel,1)

     ! Hyperbolic RT solver
     if(r%rt_advect)call r_rt_godunov_fine(pst,ilevel,1)

     ! Add anisotropic radiation from other sources
     if(r%rt_nsource>0)call r_rt_input_source_regions(pst,ilevel,1)

     ! Set rtuold equal to rtunew if not rt_smooth
     ! Add isotropic sources if not neq_chem
     call r_set_rtuold(pst,ilevel,1)

     ! Source terms for photo-chemistry
     if(r%neq_chem)call r_cooling_fine(pst,ilevel,1)

  end do
  ! End RT subcycle loop

  ! Restore original hydro timestep and time
  call m_rt_update_time(pst,ilevel,t_save,dt_save)

  ! Restriction operator to update coarser level split cells
  call m_rt_upload_fine(pst,ilevel)

  if (g%myid==1 .and. r%rt_nsubcycle .gt. 1) write(*,901) ilevel, i_substep
901 format (' Performed level', I3, ' RT-step with ', I5, ' subcycles')

  end associate

end subroutine m_rt_step
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine m_rt_update_time(pst,ilevel,t,dt)
  use amr_parameters, only: n_frw
  use ramses_commons, only: pst_t
  use update_time_module, only: in_broadcast_aexp_t, r_broadcast_aexp
  use newdt_fine_module, only: in_broadcast_dt_t, r_broadcast_dt
  use mdl_module
  implicit none
  type(pst_t)::pst
  integer::ilevel
  real(kind=8)::dt,t

  type(in_broadcast_dt_t)::in_broadcast_dt
  type(in_broadcast_aexp_t)::in_broadcast_aexp
  integer::i

  associate(r=>pst%s%r,g=>pst%s%g,m=>pst%s%m,p=>pst%s%p,mdl=>pst%s%mdl)

  ! Store time step to global variable
  g%dtnew(ilevel)=dt

  ! Broadcast dtnew to all CPUs
  in_broadcast_dt%ilevel=ilevel
  in_broadcast_dt%dtnew=g%dtnew(ilevel)
  in_broadcast_dt%dtold=g%dtold(ilevel)
  call r_broadcast_dt(pst,in_broadcast_dt,storage_size(in_broadcast_dt)/32)

  ! Store time to global variable
  g%t=t

  ! Compute expansion factor, hubble rate and proper time
  if(r%cosmo)then
     ! Find neighboring times
     i=1
     do while(g%tau_frw(i)>g%t.and.i<n_frw)
        i=i+1
     end do
     ! Interpolate expansion factor
     g%aexp = g%aexp_frw(i  )*(g%t-g%tau_frw(i-1))/(g%tau_frw(i  )-g%tau_frw(i-1))+ &
            & g%aexp_frw(i-1)*(g%t-g%tau_frw(i  ))/(g%tau_frw(i-1)-g%tau_frw(i  ))
     g%hexp = g%hexp_frw(i  )*(g%t-g%tau_frw(i-1))/(g%tau_frw(i  )-g%tau_frw(i-1))+ &
            & g%hexp_frw(i-1)*(g%t-g%tau_frw(i  ))/(g%tau_frw(i-1)-g%tau_frw(i  ))
     g%texp =    g%t_frw(i  )*(g%t-g%tau_frw(i-1))/(g%tau_frw(i  )-g%tau_frw(i-1))+ &
            &    g%t_frw(i-1)*(g%t-g%tau_frw(i  ))/(g%tau_frw(i-1)-g%tau_frw(i  ))
  else
     g%aexp = 1.0
     g%hexp = 0.0
     g%texp = g%t
  end if

  ! Broadcast t, aexp, texp and hexp to all CPUs
  in_broadcast_aexp%t=g%t
  in_broadcast_aexp%texp=g%texp
  in_broadcast_aexp%aexp=g%aexp
  in_broadcast_aexp%aexp_old=g%aexp_old
  in_broadcast_aexp%hexp=g%hexp
  in_broadcast_aexp%nstep=g%nstep
  in_broadcast_aexp%nstep_coarse=g%nstep_coarse
  call r_broadcast_aexp(pst,in_broadcast_aexp,storage_size(in_broadcast_aexp)/32)

  end associate

end subroutine m_rt_update_time
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
end module rt_step_module
