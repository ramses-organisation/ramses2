module init_neq_chem_module
contains
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine init_neq_chem(r,g,tables)
  use amr_commons, only: run_t, global_t
  use hydro_parameters, only: nion
#ifdef DO_RTZ
  use rtz_cooling_module, only: rtz_set_model
#else
  use neq_cooling_module, only: neq_set_model, update_metal_cooling, HsurH0, update_UV_rates
#ifndef DO_RT
  use photoionization_UVB_module, only: load_UVB_data
#endif
#endif
  use coolrates_module, only: neq_cooling_t, update_rt_c, update_coolrates_tables
#ifdef DO_RTZ
  use cosmic_ray_ionization_module, only: initialize_cr_rates
  use photoionization_UVB_module, only: load_UVB_data, update_UVB
  use charge_exchange_module, only: load_ct_rates
  use rtz_coolrates_module, only: initialize_high_temperature_metal_cooling, initialize_fine_structure_tables
#endif
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(neq_cooling_t)::tables

  real(kind=8)::T2_sim
  real(kind=8),dimension(1:nion)::xion_sim
  real(kind=8)::xHII_sim

  if(g%myid==1)write(*,*)'Computing non-equilibrium chemistry model'

  call update_rt_c(r, g, tables)

  if(r%cosmo)then
#ifdef DO_RTZ
     call rtz_set_model(r,tables,dble(g%h0/100.),dble(g%omega_b),dble(g%omega_m),dble(g%omega_l),dble(g%aexp_ini),T2_sim,xHII_sim)
     g%xHII_start=xHII_sim
#else
     call neq_set_model(r,tables,dble(g%h0/100.),dble(g%omega_b),dble(g%omega_m),dble(g%omega_l),dble(g%aexp_ini),T2_sim,xion_sim)
     g%xion_start=xion_sim
#endif
     g%T2_start=T2_sim
     if(r%nrestart==0)then
        if(g%myid==1)write(*,*)'Starting with T/mu (K) = ',g%T2_start
     end if
  else
#ifdef DO_RTZ
     call rtz_set_model(r,tables,dble(70./100.),dble(0.049),dble(0.3),dble(0.7),dble(r%aexp_ini),T2_sim,xHII_sim)
#else
     call neq_set_model(r,tables,dble(70./100.),dble(0.049),dble(0.3),dble(0.7),dble(r%aexp_ini),T2_sim,xion_sim)
#endif
  endif

#ifndef DO_RTZ
  ! Cosmological rates at the epoch the run starts from (restart epoch if any)
  if(r%cosmo)call update_coolrates_tables(r, tables, dble(g%aexp),          &
       dble(g%h0)*3.24077927d-20*HsurH0(1d0/dble(g%aexp)-1d0,              &
       dble(g%omega_m), dble(g%omega_l), 1.0d0-dble(g%omega_m)-dble(g%omega_l)))
  ! UV suppression of metal cooling at the same epoch (zero without cosmology)
  tables%phi=0d0
  if(r%cosmo)call update_metal_cooling(tables, dble(g%aexp))
  ! UV background rates at the same epoch (zero with DO_RT),
  ! at z=0 without cosmology (g%aexp=1), as for equilibrium cooling
#ifndef DO_RT
  if(r%haardt_madau)call load_UVB_data()
#endif
  call update_UV_rates(r, tables, dble(g%aexp))
#endif

#ifdef DO_RTZ
  ! If we are running with RTZ than we have numerous other initializations

  ! Initialize cosmic ray data
  call initialize_cr_rates()

  ! Initialize the UV background data
  call load_UVB_data()

  if(r%cosmo)then
     call update_UVB((1.d0/dble(g%aexp_ini))-1.d0)
  else
     call update_UVB((1.d0/dble(r%aexp_ini))-1.d0)
  end if

  !Initialize the charge transfer rates
  call load_ct_rates()

  ! Initialize high temperature cooling tables
  call initialize_high_temperature_metal_cooling()

  ! Initialize the low temperature cooling tables
  call initialize_fine_structure_tables()

#endif

end subroutine init_neq_chem
!###########################################################
!###########################################################
!###########################################################
!###########################################################
end module init_neq_chem_module
