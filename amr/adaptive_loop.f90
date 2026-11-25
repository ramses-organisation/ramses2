subroutine adaptive_loop(pst)
  use mdl_module
  use ramses_commons, only: pst_t
  use init_amr_module, only: r_init_amr
  use params_module, only: m_read_params
  use init_time_module, only: r_init_time
  use init_hydro_module, only: r_init_hydro
  use init_rt_module, only: r_init_rt
  use init_cr_module, only: r_init_cr
  use init_part_module, only: r_init_part
  use init_xion_module, only: m_init_xion
  use input_part_module, only: m_input_part
  use init_refine_basegrid_module, only: m_init_refine_basegrid
  use init_refine_adaptive_module, only: m_init_refine_adaptive
  use init_refine_restart_module, only: m_init_refine_restart
  use init_refine_ramses_module, only: m_init_refine_ramses
  use turb_init_module, only: r_init_turb
  use amr_step, only: m_amr_step
  use update_time_module, only: getmem, writemem, r_hash_stats
  use load_balance_module, only: r_balance_part
  use clump_finder_module, only: m_clump_finder
  use refine_utils, only: m_refine_fine
  use color_module, only: c_bold, c_dim, c_green, c_mark, c_reset
#ifdef _CUDA
  use gpu_manager, only: r_set_grid_device
  use nvtx
#endif
#ifdef _METAL
  use metal_runner, only: r_set_grid_device
#endif

  implicit none
  type(pst_t)::pst
  logical::done

  ! Local variables
  integer::ilevel, dummy
  double precision::tt1,tt2
  real(kind=4)::core_mem
#ifdef _METAL
  real(kind=8), external :: wallclock
#endif
#ifdef _CUDA
  character(len=32) :: str_step
  integer step
  step = 0
#endif

  associate(r=>pst%s%r,g=>pst%s%g,mdl=>pst%s%mdl)

  tt1 = mdl_wtime(mdl)

  ! Read run parameters
  call m_read_params(pst)

  ! Initialize grid variables
  call r_init_amr(pst)

  ! Initialize time variables and cooling tables
  call r_init_time(pst)

  ! Initialize hydro kernel workspace
  if(r%hydro)call r_init_hydro(pst)

  ! Initialize rt kernel workspace
  if(r%rt)call r_init_rt(pst)

  ! Initialize rt kernel workspace
  if(r%cr)call r_init_cr(pst)

  ! Initialize particle variables
  if(r%pic)call r_init_part(pst)

  ! Initialize turbulent driveing
  if(r%turb)call r_init_turb(pst)

  ! Read initial particle properties from files
  if(r%pic)call m_input_part(pst)

  ! Build initial AMR grid
  if(r%nrestart==0)then
     if(r%filetype=='ramses'.and.r%hydro)then
        call m_init_refine_ramses(pst) ! Build AMR grid from output file
     else
        call m_init_refine_basegrid(pst) ! Build coarse grid
        call m_init_refine_adaptive(pst) ! Build adaptive grid
     endif
  else
     call m_init_refine_restart(pst) ! Build AMR grid from restart file
     g%first_coarse_restart = .true. ! Indicate the initial course step post restart
  endif

  ! Initialization of ionization fractions
  if(r%neq_chem .and. r%is_init_xion) call m_init_xion(pst)

  ! Timing since startup
  tt2 = mdl_wtime(mdl)
  print '(A,F0.7,A)',c_dim()//' Time elapsed since startup: ',tt2-tt1,c_reset()

  ! Output mesh structure
  do ilevel=r%levelmin,r%nlevelmax
     if(pst%s%m%noct_tot(ilevel)>0)write(*,999)c_bold(),&
          & ilevel,pst%s%m%noct_tot(ilevel),c_reset()//c_dim(),pst%s%m%noct_min(ilevel),pst%s%m%noct_max(ilevel),pst%s%m%noct_tot(ilevel)/mdl_threads(mdl),c_reset()
  end do
999 format(' ',A,'Level ',I0,' has ',I0,' grids ',A,'(',I0,',',I0,',',I0,')',A)

  g%nstep_coarse_old=g%nstep_coarse
  if(r%nrestart>0)g%nstep_coarse_old=g%nstep_coarse-1

#ifdef _CUDA
  ! Copy entire grid from host to device
  call r_set_grid_device(pst)
#endif
#ifdef _METAL
  call r_set_grid_device(pst)
#endif

  ! Just in case we only do clump finding
  if(r%clump_only)then
     if(mdl_threads(mdl)>1)then
        write(*,*)'Load balancing particle distribution'
        tt1 = mdl_wtime(mdl)
        call r_balance_part(pst,r%levelmin,1,dummy,0)
        tt2 = mdl_wtime(mdl)
        print '(A,F0.7)',' Time elapsed load balancing: ',tt2-tt1
     endif
     call m_refine_fine(pst,r%levelmin)
     call m_clump_finder(pst,.true.,.false.)
     return
  endif

  ! Static mesh runs never refine
  if(r%static_mesh)call m_refine_fine(pst,r%levelmin)

  write(*,'(1x,a)')c_bold()//c_green()//c_mark()//'Starting time integration'//c_reset()

  done = .false.
  do while(.not.done) ! Main time loop

#ifdef _CUDA
     write(str_step,'(A,I0)'),"step_",step
     call nvtxStartRange(trim(str_step), color=5)
     step = step + 1
#endif

#ifdef _METAL
     tt1 = wallclock()
#else
     tt1 = mdl_wtime(mdl)
#endif

     if(r%verbose)write(*,*)'Entering amr_step_coarse'

     g%epot_tot=0.0D0  ! Reset total potential energy
     g%ekin_tot=0.0D0  ! Reset total kinetic energy
     g%mass_tot=0.0D0  ! Reset total mass
     g%eint_tot=0.0D0  ! Reset total internal energy
     g%emag_tot=0.0D0  ! Reset total magnetic energy

     ! Call base level
     call m_amr_step(pst,r%levelmin,1,done)

     ! New coarse time-step
     g%nstep_coarse=g%nstep_coarse+1

#ifdef _METAL
     tt2 = wallclock()
#else
     tt2 = mdl_wtime(mdl)
#endif
     if(g%control_due)then
        if(.not. done)print '(A,F0.7,A)',c_dim()//' Time elapsed since last coarse step: ',tt2-tt1,c_reset()
     endif

     call getmem(core_mem)
     call writemem(core_mem)

#ifdef _CUDA
     call nvtxEndRange()
#endif
  end do

  call m_output_timer(.false.,'dummy')

!  call r_hash_stats(pst)

  return

  end associate

end subroutine adaptive_loop
