module star_formation_common_module
#if defined(_CUDA) && defined(DO_RT)
  use gpu_manager, only: gpu_star_to_host
#endif
  implicit none
contains

subroutine finish_star_formation(s)
  use ramses_commons, only: ramses_t
  use SED_module, only: update_SED_group_props
  implicit none
  type(ramses_t)::s

  ! Check if radiation advection should be turned on
  if(s%r%rt .and. s%r%rt_star .and. .not. s%r%rt_advect .and. s%star%npart_tot .gt. 0) then
     if(s%g%myid==1) then
       write(*,*) '*****************************************'
       if(s%r%cosmo) then
         write(*,*) 'Stellar RT turned on at a=',s%g%aexp
       else
         write(*,*) 'Stellar RT turned on at t=',s%g%t
       endif
       write(*,*) '*****************************************'
     endif
     s%r%rt_advect=.true.
     ! Update cross sections based on star properties
     if(s%r%sedprops_update .gt. 0) then
#if defined(_CUDA) && defined(DO_RT)
       call gpu_star_to_host(s)
#endif
       call update_SED_group_props(s%r, s%g, s%sed, s%star)
     endif
  endif

  ! Check if CR advection should be turned on
  if(s%r%cr .and. .not. s%r%cr_advect .and. s%star%npart_tot .gt. 0) then
     if(s%g%myid==1) then
       write(*,*) '*****************************************'
       if(s%r%cosmo) then
         write(*,*) 'CR advection turned on at a=',s%g%aexp
       else
         write(*,*) 'CR advection turned on at t=',s%g%t
       endif
       write(*,*) '*****************************************'
     endif
     s%r%cr_advect=.true.
  endif

end subroutine finish_star_formation

end module star_formation_common_module
