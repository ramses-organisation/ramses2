module newdt_fine_module
#ifdef _CUDA
  use part_device, only: gpu_newdt_part, gpu_newdt_star, gpu_newdt_sink, gpu_newdt_tree
  use part_device, only: gpu_newdt_core
#elif defined(_METAL)
  use metal_runner, only: metal_newdt_part
#endif

type :: out_newdt_part_t
  real(kind=8)::ekin,vmax
end type out_newdt_part_t

type :: out_max_bq_t
  real(kind=8)::max_b, max_q
end type out_max_bq_t

type :: in_broadcast_dt_t
  integer::ilevel
  real(kind=8)::dtnew,dtold
end type in_broadcast_dt_t

contains
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine m_newdt_fine(pst,ilevel)
  use amr_parameters, only: nvector, n_frw
  use ramses_commons, only: pst_t
  use courant_fine_module, only: r_courant_fine, out_courant_fine_t
  use constants, only: twopi
  use update_cr_c_module
  implicit none
  type(pst_t)::pst
  integer::ilevel
  !-----------------------------------------------------------
  ! This routine compute the time step using 3 constraints:
  ! 1- a Courant-type condition using particle velocity
  ! 2- the gravity free-fall time
  ! 3- 10% maximum variation for aexp
  ! This routine also compute the particle kinetic energy.
  !-----------------------------------------------------------
  real(kind=8)::dx,tff,fourpi,threepi2
  real(kind=8)::ekin,vmax
  real(kind=8)::dt_gyro, dt_courant, max_b, max_q, cr_c
  real(kind=8)::tnext,trem,dtest
  integer::nrem,i
  type(out_courant_fine_t)::out_courant_fine
  type(out_newdt_part_t)::out_newdt_part
  type(in_broadcast_dt_t)::in_broadcast_dt
  type(in_broadcast_cr_c_t)::in_broadcast_cr_c

  associate(r=>pst%s%r,g=>pst%s%g,m=>pst%s%m,p=>pst%s%p,mdl=>pst%s%mdl)

  if(m%noct_tot(ilevel)==0)return
  if(r%verbose)write(*,'("   Entering newdt_fine for level ",I2)')ilevel

  ! Save old time step
  g%dtold(ilevel)=g%dtnew(ilevel)

  ! Compute local cell spacing
  dx=r%boxlen/2**ilevel

  ! Maximum time step
  g%dtnew(ilevel)=dx/r%smallc

  ! Gravity-based Courant condition
  if(r%poisson.and.r%gravity_type<=0)then
     fourpi=4.0d0*ACOS(-1.0d0)
     if(r%cosmo)fourpi=1.5d0*g%omega_m*g%aexp
     threepi2=3.0d0*ACOS(-1.0d0)**2
     if(g%rho_max(ilevel)>0)then
        tff=sqrt(threepi2/8./fourpi/g%rho_max(ilevel))
        g%dtnew(ilevel)=MIN(g%dtnew(ilevel),r%courant_factor*tff)
     endif
  end if

  ! Cosmic-expansion-based Courant condition
  if(r%cosmo)then
     g%dtnew(ilevel)=MIN(g%dtnew(ilevel),0.1/g%hexp)
  end if

  ! Turbulence driving condition
  if(r%turb)then
     g%dtnew(ilevel)=MIN(g%dtnew(ilevel),pst%s%turb%turb_dt)
  end if

  ! Particle-based Courant condition
  if(r%pic)then
     if(r%verbose)write(*,'("   Entering newdt_part for level ",I2)')ilevel
     call r_newdt_part(pst,ilevel,1,out_newdt_part,storage_size(out_newdt_part)/32)
     ekin=out_newdt_part%ekin
     g%ekin_tot=g%ekin_tot+ekin
     vmax=out_newdt_part%vmax
     if(vmax>0.0d0)then
        g%dtnew(ilevel)=MIN(real(g%dtnew(ilevel),kind=8),r%courant_factor*dx/vmax)
     endif
  endif

  ! Hydro-based Courant condition
  if(r%hydro)then
     if(r%verbose)write(*,'("   Entering newdt_hydro for level ",I2)')ilevel
     call r_courant_fine(pst,ilevel,1,out_courant_fine,storage_size(out_courant_fine)/32)
     g%mass_tot=g%mass_tot+out_courant_fine%mass
     g%ekin_tot=g%ekin_tot+out_courant_fine%ekin
     g%eint_tot=g%eint_tot+out_courant_fine%eint
     g%emag_tot=g%emag_tot+out_courant_fine%emag
     g%dtnew(ilevel)=MIN(real(g%dtnew(ilevel),kind=8),out_courant_fine%dt)
  endif

  ! Dust gyro-frequency timestep constraint (MHD + dust only)
#ifdef MHD
  if(r%dust)then
     ! Compute maximum |B| at this level and maximum particle charge |q|
     block
       type(out_max_bq_t)::bq
       call r_max_B_and_Q(pst, ilevel, 1, bq, storage_size(bq)/32)
       max_b = bq%max_b
       max_q = bq%max_q
     end block
     if(max_b>0.0d0 .and. max_q>0.0d0)then
        dt_gyro = r%dust_gyro_factor * twopi / (max_q * max_b)
        g%dtnew(ilevel)=MIN(real(g%dtnew(ilevel),kind=8), dt_gyro)
     endif
  endif
#endif

  if(r%rt.and.r%rt_advect)then
     call get_rt_courant_dt(r,g,dt_courant,ilevel)
     g%dtnew(ilevel)=MIN(real(g%dtnew(ilevel),kind=8),r%rt_nsubcycle*dt_courant)
  endif

  if(r%cr.and.r%cr_advect)then
     call get_cr_c(r, g, cr_c)
     g%cr_c(ilevel) = cr_c
     ! Broadcast the new value of cr_c
     in_broadcast_cr_c%ilevel=ilevel
     in_broadcast_cr_c%cr_c=g%cr_c(ilevel)
     call r_broadcast_cr_c(pst,in_broadcast_cr_c,storage_size(in_broadcast_cr_c)/32)

     ! Limit timestep to Nsubcycle times dt_cr
     call get_cr_courant_dt(r,g,dt_courant,ilevel)
     g%dtnew(ilevel)=MIN(real(g%dtnew(ilevel),kind=8),r%cr_nsubcycle*dt_courant)
  endif

  ! Adaptive time step condition
  if(ilevel>r%levelmin)then
     g%dtnew(ilevel)=MIN(g%dtnew(ilevel-1)/real(r%nsubcycle(ilevel-1)),g%dtnew(ilevel))
  end if

  ! Land the coarse step on the next output time: within nlanding coarse
  ! steps of it, use dt=trem/nrem, nrem=ceiling(trem/dt), so that the time
  ! step reduction is spread over the last nrem steps (at most 1/nlanding each)
  ! With cosmo, the output time is the conformal time of aout, obtained by
  ! inverting the Friedmann table interpolation used in update_time
  if(ilevel==r%levelmin.and.r%nlanding>0.and.g%iout<=r%noutput)then
     if(r%cosmo)then
        i=1
        do while(g%aexp_frw(i)>r%aout(g%iout).and.i<n_frw)
           i=i+1
        end do
        tnext=g%tau_frw(i  )*(r%aout(g%iout)-g%aexp_frw(i-1))/(g%aexp_frw(i  )-g%aexp_frw(i-1))+ &
            & g%tau_frw(i-1)*(r%aout(g%iout)-g%aexp_frw(i  ))/(g%aexp_frw(i-1)-g%aexp_frw(i  ))
     else
        tnext=r%tout(g%iout)
     endif
     trem=tnext-g%t
     if(trem>0)then
        ! The previous coarse step includes the finer level constraints
        dtest=g%dtnew(ilevel)
        if(g%dtold(ilevel)>0)dtest=MIN(dtest,g%dtold(ilevel))
        if(trem<r%nlanding*dtest)then
           ! Tolerance: trem/dtest=3+round-off must give 3 steps, not 4
           nrem=max(1,ceiling(trem/dtest-1d-6))
           g%dtnew(ilevel)=trem/dble(nrem)
           ! Last step: overshoot by round-off only, so that t>=tout or aexp>=aout
           ! (trem included, since the conformal time is 0 at aout=1)
           if(nrem==1)g%dtnew(ilevel)=trem+1d-10*max(trem,abs(tnext))
        endif
     endif
  endif
  if(r%verbose)write(*,'("   New time step for level ",I2," is ",1PE12.5)')ilevel,g%dtnew(ilevel)

  ! Broadcast new and old time steps to all CPUs
  in_broadcast_dt%ilevel=ilevel
  in_broadcast_dt%dtnew=g%dtnew(ilevel)
  in_broadcast_dt%dtold=g%dtold(ilevel)
  call r_broadcast_dt(pst,in_broadcast_dt,storage_size(in_broadcast_dt)/32)

  end associate

end subroutine m_newdt_fine
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
recursive subroutine r_max_B_and_Q(pst, ilevel, input_size, output, output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  use amr_parameters, only: ndim, twotondim
  implicit none
  type(pst_t)::pst
  integer::ilevel
  integer,VALUE::input_size
  integer::output_size
  type(out_max_bq_t)::output, next_output

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_MAX_B_AND_Q,pst%iUpper+1,input_size,output_size,ilevel)
     call r_max_B_and_Q(pst%pLower,ilevel,input_size,output,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_output)
     output%max_b = MAX(output%max_b, next_output%max_b)
     output%max_q = MAX(output%max_q, next_output%max_q)
  else
     call max_B_and_Q(pst%s%r,pst%s%g,pst%s%m,pst%s%dust,ilevel,output%max_b,output%max_q)
  endif

end subroutine r_max_B_and_Q
!#####################################################################
!#####################################################################
subroutine max_B_and_Q(r,g,m,p,ilevel,max_b,max_q)
  use amr_parameters, only: ndim, twotondim
  use amr_commons, only: run_t, global_t, mesh_t
  use pm_commons, only: part_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
  integer::ilevel
  real(kind=8)::max_b, max_q

  integer::igrid, ind, idim, ipart
  real(kind=8)::bx, by, bz, bmag

  max_b = 0.0d0
  max_q = 0.0d0

#ifdef MHD
  do igrid=m%head(ilevel), m%tail(ilevel)
     do ind=1,twotondim
        bx=0.5d0*(m%bold(ind,1,igrid)+m%bold(ind,4,igrid))
        by=0.5d0*(m%bold(ind,2,igrid)+m%bold(ind,5,igrid))
        bz=0.5d0*(m%bold(ind,3,igrid)+m%bold(ind,6,igrid))
        bmag = sqrt(bx*bx + by*by + bz*bz)
        if(bmag>max_b) max_b=bmag
     end do
  end do
#endif

  if(r%dust)then
     do ipart=p%headp(ilevel), p%tailp(ilevel)
        if (abs(p%charge(ipart))>max_q) max_q=abs(p%charge(ipart))
     end do
  end if

end subroutine max_B_and_Q
recursive subroutine r_newdt_part(pst,ilevel,input_size,output,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  type(out_newdt_part_t)::output,next_output
  integer::ilevel

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_NEWDT_PART,pst%iUpper+1,input_size,output_size,ilevel)
     call r_newdt_part(pst%pLower,ilevel,input_size,output,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_output)
     output%ekin=output%ekin+next_output%ekin
     output%vmax=MAX(output%vmax,next_output%vmax)
  else
     output%vmax=0.0d0
     output%ekin=0.0d0
#ifdef _CUDA
     if (pst%s%r%part) call gpu_newdt_part(pst%s, ilevel, output%vmax, output%ekin)
     if (pst%s%r%star) call gpu_newdt_star(pst%s, ilevel, output%vmax, output%ekin)
     if (pst%s%r%sink) call gpu_newdt_sink(pst%s, ilevel, output%vmax, output%ekin)
     if (pst%s%r%core) call gpu_newdt_core(pst%s, ilevel, output%vmax, output%ekin)
     if (pst%s%r%tree) call gpu_newdt_tree(pst%s, ilevel, output%vmax, output%ekin)
     return
#elif defined(_METAL)
     if (pst%s%r%part) call metal_newdt_part(pst%s, ilevel, output%vmax, output%ekin)
     return
#endif
     if (pst%s%r%part) call newdt_part(pst%s%r,pst%s%g,pst%s%p   ,ilevel,output%ekin,output%vmax)
     if (pst%s%r%star) call newdt_part(pst%s%r,pst%s%g,pst%s%star,ilevel,output%ekin,output%vmax)
     if (pst%s%r%sink) call newdt_part(pst%s%r,pst%s%g,pst%s%sink,ilevel,output%ekin,output%vmax)
     if (pst%s%r%tree) call newdt_part(pst%s%r,pst%s%g,pst%s%tree,ilevel,output%ekin,output%vmax)
     if (pst%s%r%trac) call newdt_part(pst%s%r,pst%s%g,pst%s%trac,ilevel,output%ekin,output%vmax)
     if (pst%s%r%dust) call newdt_part(pst%s%r,pst%s%g,pst%s%dust,ilevel,output%ekin,output%vmax)
     if (pst%s%r%core) call newdt_part(pst%s%r,pst%s%g,pst%s%core,ilevel,output%ekin,output%vmax)
  endif

end subroutine r_newdt_part
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine newdt_part(r,g,p,ilevel,ekin,vmax)
  use amr_parameters, only: ndim
  use amr_commons, only: run_t, global_t
  use pm_commons, only: part_t
  use pm_parameters, only: TREE_TYPE, TRAC_TYPE
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(part_t)::p
  integer::ilevel
  real(kind=8)::ekin,vmax

  integer::ipart,idim

  ! Compute maximum particle velocity
  do idim = 1, ndim
     do ipart = p%headp(ilevel), p%tailp(ilevel)
        vmax = MAX(vmax, ABS(dble(p%vp(ipart, idim))))
     end do
  end do

  if(p%type==TREE_TYPE .or. p%type==TRAC_TYPE)return

  ! Compute kinetic energy
  do idim = 1, ndim
     do ipart = p%headp(ilevel), p%tailp(ilevel)
        ekin = ekin + 0.5D0 * p%mp(ipart) * p%vp(ipart, idim)**2
     end do
  end do

end subroutine newdt_part
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
recursive subroutine r_broadcast_dt(pst,input,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  type(in_broadcast_dt_t)::input

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_BROADCAST_DT,pst%iUpper+1,input_size,0,input)
     call r_broadcast_dt(pst%pLower,input,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     pst%s%g%dtnew(input%ilevel)=input%dtnew
     pst%s%g%dtold(input%ilevel)=input%dtold
  endif

end subroutine r_broadcast_dt
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
end module newdt_fine_module
