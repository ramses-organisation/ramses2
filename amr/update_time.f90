module update_time_module
  type :: in_broadcast_aexp_t
    real(kind=8)::t,texp,aexp,hexp
    real(kind=8)::aexp_old
    integer::nstep,nstep_coarse
  end type in_broadcast_aexp_t
contains
!################################################################
!################################################################
!################################################################
!################################################################
subroutine m_update_time(pst,ilevel,done)
  use amr_parameters, only: n_frw
  use ramses_commons, only: pst_t
  use mdl_module
  use update_rt_c_module, only: r_rt_neq_updates
  use update_cr_c_module, only: r_cr_updates
  use turb_update_module, only: r_update_turb
  implicit none
  type(pst_t)::pst
  integer::ilevel
  logical::done

  ! Local variables
  double precision::ttend
  double precision,save::ttstart=0.0
  real(kind=8)::dt,econs,mcons
  integer::i,itest
  type(in_broadcast_aexp_t)::in_broadcast_aexp
#ifdef _METAL
  real(kind=8), external :: wallclock
#endif

  associate(r=>pst%s%r,g=>pst%s%g,m=>pst%s%m,p=>pst%s%p,mdl=>pst%s%mdl)

  ! Local constants
  dt=g%dtnew(ilevel)
  itest=0

#ifdef _METAL
  if(ttstart.eq.0.0) ttstart = wallclock()
#else
  if(ttstart.eq.0.0) ttstart = mdl_wtime(mdl)
#endif

  if(g%first_coarse_restart)g%first_coarse_restart = .false.

  !-------------------------------------------------------------
  ! At this point, IF nstep_coarse has JUST changed, all levels
  ! are synchronised, and all new refinements have been done.
  !-------------------------------------------------------------
  if(g%nstep_coarse .ne. g%nstep_coarse_old)then

     !----------------------------------------------
     ! RT, non-equilibrium chemistry and CR updates
     !----------------------------------------------
     call r_rt_neq_updates(pst, pst%s%g%nstep_coarse, 1)
     call r_cr_updates(pst, pst%s%g%nstep_coarse, 1)

     !--------------------------
     ! Check mass conservation
     !--------------------------
     if(g%mass_tot_0==0.0D0)then
        g%mass_tot_0=g%mass_tot+g%mass_star_tot+g%mass_sink_tot
        mcons=0.0D0
     else
        mcons=(g%mass_tot+g%mass_star_tot+g%mass_sink_tot-g%mass_tot_0)/g%mass_tot_0
     end if

     !----------------------------
     ! Check energy conservation
     !----------------------------
     if(g%epot_tot_old.ne.0)then
        g%epot_tot_int=g%epot_tot_int+0.5D0*(g%epot_tot_old+g%epot_tot)*log(g%aexp/g%aexp_old)
     end if
     g%epot_tot_old=g%epot_tot
     g%aexp_old=g%aexp
     if(g%const==0.0D0)then
        g%const=g%epot_tot+g%ekin_tot  ! initial total energy
        econs=0.0D0
     else
        econs=(g%ekin_tot+g%epot_tot-g%epot_tot_int-g%const) / &
             &(-(g%epot_tot-g%epot_tot_int-g%const)+g%ekin_tot)
     end if

     if(g%control_due.or.g%output_done)then
        g%control_due=.false.

        !-------------------------------
        ! Output AMR structure to screen
        !-------------------------------
        write(*,*)'Mesh structure'
        do i=r%levelmin,r%nlevelmax
           if(m%noct_tot(i)>0)write(*,999)i,m%noct_tot(i),m%noct_min(i),m%noct_max(i),m%noct_tot(i)/g%ncpu
        end do
999     format(' Level ',I0,' has ',I0,' grids (',I0,',',I0,',',I0,')')

        !------------------------
        ! Output timing data
        !------------------------
        if(r%verbose)call m_output_timer(.false.,'dummy')

        !----------------------------------------------
        ! Output mass and energy conservation to screen
        !----------------------------------------------
        if(r%hydro)then
#ifdef MHD
           write(*,779)g%nstep_coarse,mcons,econs,g%epot_tot,g%ekin_tot,g%eint_tot,g%emag_tot
#else
           write(*,778)g%nstep_coarse,mcons,econs,g%epot_tot,g%ekin_tot,g%eint_tot
#endif
        else
           write(*,777)g%nstep_coarse,mcons,econs,g%epot_tot,g%ekin_tot
        end if
        if(r%star)write(*,'(" Total mass in stars=",1PE14.7)')g%mass_star_tot
        if(r%sink.or.r%core)write(*,'(" Total mass in sinks=",1PE14.7)')g%mass_sink_tot
777     format(' Main step= ',i0,' mcons=',1pe9.2,' econs=',1pe9.2,' epot=',1pe9.2,' ekin=',1pe9.2)
778     format(' Main step= ',i0,' mcons=',1pe9.2,' econs=',1pe9.2,' epot=',1pe9.2,' ekin=',1pe9.2,' eint=',1pe9.2)
779     format(' Main step= ',i0,' mcons=',1pe9.2,' econs=',1pe9.2,' epot=',1pe9.2,' ekin=',1pe9.2,' eint=',1pe9.2,' emag=',1pe9.2)

        !----------------------------------------------
        ! Output fine step information and used memory
        !----------------------------------------------
        if(r%part)then
#if defined(_CUDA) || defined(_METAL)
           write(*,888)g%nstep,g%t,dt,g%aexp,real(100.0D0*dble(m%noct_used_max)/dble(m%ngridmax)),&
                & real(100.0D0*dble(m%ifree_cache)/dble(m%ncachemax)),&
                & real(100.0D0*dble(p%npart_max)/dble(r%npartmax+1))
#else
           write(*,888)g%nstep,g%t,dt,g%aexp,real(100.0D0*dble(m%noct_used_max)/dble(m%ngridmax)),&
                & real(100.0D0*dble(p%npart_max)/dble(r%npartmax+1))
#endif
        else
#if defined(_CUDA) || defined(_METAL)
           write(*,888)g%nstep,g%t,dt,g%aexp,real(100.0D0*dble(m%noct_used_max)/dble(m%ngridmax)),&
                & real(100.0D0*dble(m%ifree_cache)/dble(m%ncachemax))
#else
                write(*,888)g%nstep,g%t,dt,g%aexp,real(100.0D0*dble(m%noct_used_max)/dble(m%ngridmax))
#endif
        endif
        itest=1
     end if
     g%output_done=.false.

     !---------------
     ! Exit program
     !---------------
     if(g%t>=r%tout(r%noutput).or.g%aexp>=r%aout(r%noutput).or.g%nstep_coarse>=r%nstepmax)then
        write(*,*)'Run completed'
#ifdef _METAL
        ttend = wallclock()
#else
        ttend = mdl_wtime(mdl)
#endif
        print '(A,F0.7)',' Total elapsed time: ',ttend-ttstart
        done=.true.
        return
        !call mdl_abort(mdl)
     end if

  end if
  g%nstep_coarse_old=g%nstep_coarse

  !----------------------------
  ! Output controls to screen
  !----------------------------
  if(mod(g%nstep,r%ncontrol)==0)then
     if(itest==0)then
        g%control_due=.true.
        if(r%part)then
#if defined(_CUDA) || defined(_METAL)
           write(*,888)g%nstep,g%t,dt,g%aexp,real(100.0D0*dble(m%noct_used_max)/dble(m%ngridmax)),&
                & real(100.0D0*dble(m%ifree_cache)/dble(m%ncachemax)),&
                & real(100.0D0*dble(p%npart_max)/dble(r%npartmax+1))
#else
           write(*,888)g%nstep,g%t,dt,g%aexp,real(100.0D0*dble(m%noct_used_max)/dble(m%ngridmax)),&
                & real(100.0D0*dble(p%npart_max)/dble(r%npartmax+1))
#endif
        else
#if defined(_CUDA) || defined(_METAL)
           write(*,888)g%nstep,g%t,dt,g%aexp,real(100.0D0*dble(m%noct_used_max)/dble(m%ngridmax)),&
                & real(100.0D0*dble(m%ifree_cache)/dble(m%ncachemax))
#else
                write(*,888)g%nstep,g%t,dt,g%aexp,real(100.0D0*dble(m%noct_used_max)/dble(m%ngridmax))
#endif
        endif
     end if
  end if
888 format(' Fine step= ',i0,' t=',1pe12.5,' dt=',1pe10.3,' a=',1pe10.3,' mem=',0pF0.1,'% ',0pF0.1,'% ',0pF0.1,'%')

  !------------------------
  ! Update time variables
  !------------------------
  g%t=g%t+dt
  g%nstep=g%nstep+1
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

  ! Broadcast t, aexp, aexp_old, texp and hexp to all CPUs
  in_broadcast_aexp%t=g%t
  in_broadcast_aexp%texp=g%texp
  in_broadcast_aexp%aexp=g%aexp
  in_broadcast_aexp%aexp_old = g%aexp_old
  in_broadcast_aexp%hexp=g%hexp
  in_broadcast_aexp%nstep=g%nstep
  in_broadcast_aexp%nstep_coarse=g%nstep_coarse
  call r_broadcast_aexp(pst,in_broadcast_aexp,storage_size(in_broadcast_aexp)/32)

  ! Update turbulent driving field
  call m_timer('hydro - turbulence','start')
  if(r%turb)call r_update_turb(pst)

end associate

end subroutine m_update_time
!##############################################################
!##############################################################
!##############################################################
!##############################################################
recursive subroutine r_broadcast_aexp(pst,input,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  type(in_broadcast_aexp_t)::input

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_BROADCAST_AEXP,pst%iUpper+1,input_size,0,input)
     call r_broadcast_aexp(pst%pLower,input,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     pst%s%g%t   =input%t
     pst%s%g%texp=input%texp
     pst%s%g%aexp=input%aexp
     pst%s%g%aexp_old = input%aexp_old
     pst%s%g%hexp=input%hexp
     pst%s%g%nstep=input%nstep
     pst%s%g%nstep_coarse=input%nstep_coarse
  endif

end subroutine r_broadcast_aexp
!##############################################################
!##############################################################
!##############################################################
!##############################################################
recursive subroutine r_hash_stats(pst)
  use mdl_module
  use ramses_commons, only: pst_t
  use hash, only: hash_stats
  use mdl_parameters
  implicit none
  type(pst_t)::pst

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_HASH_STATS,pst%iUpper+1)
     call r_hash_stats(pst%pLower)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call hash_stats(pst%s%m%grid_dict)
  endif

end subroutine r_hash_stats
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine writemem(usedmem)
  real::usedmem

  usedmem=usedmem*4096

  if(usedmem>1024.**3.)then
     write(*,999)usedmem/1024.**3.
  else if (usedmem>1024.**2.) then
     write(*,998)usedmem/1024.**2
  else if (usedmem>1024.) then
     write(*,997)usedmem/1024.
  endif

997 format(' Used memory: ',F0.1,' kb')
998 format(' Used memory: ',F0.1,' Mb')
999 format(' Used memory: ',F0.1,' Gb')

end subroutine writemem

subroutine getmem(outmem)
  use iso_c_binding, only: c_long
  real::outmem
  character(len=300) :: dir, dir2, file
  integer::ind,j,nmem,read_status
  logical::file_exists
#ifdef _METAL
  interface
    function getmem_mac() bind(c, name='getmem_mac')
      import c_long
      integer(c_long) :: getmem_mac
    end function getmem_mac
  end interface
  outmem = real(getmem_mac(), kind=4)
  return
#endif

  file='/proc/self/stat'
  inquire(file=file, exist=file_exists)
  if (file_exists) then
     open(unit=1,file=file,form='formatted')
     read(1,'(A300)',IOSTAT=read_status)dir
     close(1)
  else
     read_status=-1000
  endif

  if (read_status < 0)then
     outmem=0.
     if (read_status .ne. -1000)write(*,*)'Problem in checking free memory'
  else
     ind=300
     j=0
     do while (j<23)
        ind=index(dir,' ')
        dir2=dir(ind+1:300)
        j=j+1
        dir=dir2
     end do
     ind=index(dir,' ')
     dir2=dir(1:ind)
     read(dir2,'(I12)')nmem
     outmem=real(nmem,kind=4)
  end if

end subroutine getmem
end module update_time_module
