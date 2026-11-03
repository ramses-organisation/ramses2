module input_part_module
contains
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_input_part(pst)
  use mdl_module
  use ramses_commons, only: pst_t
  use input_part_grafic_module, only: m_input_part_grafic
  use input_part_ascii_module, only: m_input_part_ascii
  use input_part_restart_module, only: m_input_part_restart
  use input_part_ramses_module, only: m_input_part_ramses
  use input_part_gadget_module, only: m_input_part_gadget
  implicit none
  type(pst_t)::pst
  !--------------------------------------------------------------------
  ! This routine is the master procedure to read and dispatch particles
  ! from many different initial conditions file formats.
  !--------------------------------------------------------------------
  real(kind=8)::mp_min

  ! Input particle properties from files
  select case(pst%s%r%filetype)
  case ('grafic')
     call m_input_part_grafic(pst)
  case ('grafic_zoom')
     return
  case('ascii')
     call m_input_part_ascii(pst)
  case('gadget')
     call m_input_part_gadget(pst)
  case('ramses')
     call m_input_part_ramses(pst)
  case('restart')
     call m_input_part_restart(pst)
  case DEFAULT
     write(*,*) 'Unsupported format file ' // pst%s%r%filetype
     call mdl_abort(pst%s%mdl)
  end select

  ! Computing maximum particle count (only in master)
  call r_npart_max(pst,pst%s%r%levelmin,1,pst%s%p%npart_max,1)

  ! Compute minimum particle mass and initial coarse particle level
  call r_mass_min_part(pst,pst%s%r%levelmin,1,mp_min,2)
  call r_broadcast_mp_min(pst,mp_min,2)

  ! Check if we should turn on radiation or CR advection
  ! and update rt groups
  call r_check_part_emission(pst)

end subroutine m_input_part
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
recursive subroutine r_mass_min_part(pst,ilevel,input_size,mp_min,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer::ilevel
  real(kind=8)::mp_min,next_mp_min
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_MASS_MIN_PART,pst%iUpper+1,input_size,output_size,ilevel)
     call r_mass_min_part(pst%pLower,ilevel,input_size,mp_min,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_mp_min)
     mp_min=MIN(mp_min,next_mp_min)
  else
     mp_min=MINVAL(pst%s%p%mp(1:pst%s%p%npart))
  endif

end subroutine r_mass_min_part
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
recursive subroutine r_broadcast_mp_min(pst,mp_min,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size

  real(kind=8)::mp_min,mm1
  integer::rID,ilevel

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_BROADCAST_MP_MIN,pst%iUpper+1,input_size,0,mp_min)
     call r_broadcast_mp_min(pst%pLower,mp_min,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     pst%s%g%mp_min=mp_min
     if(pst%s%r%aexp_lock_refine>0d0)then
        ilevel = 1
        do while(.true.)
           mm1 = 0.5d0**(3*ilevel)*(1.0d0-pst%s%g%omega_b/pst%s%g%omega_m)
           if((mm1 > 0.90d0*mp_min).AND.(mm1 < 1.10d0*mp_min))then
              pst%s%g%nlevelmax_part = ilevel
              exit
           endif
           ilevel = ilevel+1
        enddo
     endif
  endif

end subroutine r_broadcast_mp_min
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
recursive subroutine r_npart_max(pst,ilevel,input_size,npart_max,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer::ilevel
  integer::npart_max

  integer::next_npart_max
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_NPART_MAX,pst%iUpper+1,input_size,output_size,ilevel)
     call r_npart_max(pst%pLower,ilevel,input_size,npart_max,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_npart_max)
     npart_max=MAX(npart_max,next_npart_max)
  else
     npart_max=pst%s%p%npart
  endif

end subroutine r_npart_max
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
recursive subroutine r_check_part_emission(pst)

  ! Check if to turn on radiation or cr advection from stellar particles
  use mdl_module
  use ramses_commons, only: pst_t
  use SED_module, only: update_SED_group_props
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::rID

  associate(s=>pst%s,r=>pst%s%r,m=>pst%s%m,g=>pst%s%g)

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_CHECK_PART_EMISSION,pst%iUpper+1)
     call r_check_part_emission(pst%pLower)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else

     ! Check if radiadion advection should be turned on
     if(r%rt .and. r%star .and. .not. r%rt_advect                         &
        .and. pst%s%star%npart_tot .gt. 0) then
        if(g%myid==1) then
          write(*,*) 'Stellar radiation is emitted and advected'
        endif
        r%rt_advect=.true.
        ! Update cross sections based on star properties
        if(r%sedprops_update .gt. 0 .and. r%nrestart .eq. 0) then
           call update_SED_group_props(r, g, pst%s%sed, pst%s%star)
        endif
     endif

     ! Check if CR advection should be turned on
     if(r%cr .and. sum(r%fecr).gt.0.0 .and. .not. r%cr_advect            &
        .and. pst%s%star%npart_tot .gt. 0) then
        if(g%myid==1) then
          write(*,*) 'Cosmic rays are emitted and advected'
        endif
        r%cr_advect=.true.
     endif

  endif

  end associate

end subroutine r_check_part_emission

end module input_part_module
