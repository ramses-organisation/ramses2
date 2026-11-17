module input_part_ascii_module

contains
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_input_part_ascii(pst)
  use ramses_commons, only: pst_t
  use amr_parameters, only: flen
  implicit none
  type(pst_t)::pst
  !--------------------------------------------------------------------
  ! This routine is the master procedure to read and dispatch particles
  ! from an ascii file.
  !--------------------------------------------------------------------
  real(kind=8)::xx1,xx2,xx3,vv1,vv2,vv3,mm1,zz1,tt1,jj1,jj2,jj3
  integer(kind=8)::npart_tot,nstar_tot
  character(LEN=flen)::filename
  integer(kind=8)::ncore_tot
  integer,allocatable,dimension(:)::input_array
  logical::file_exist

  associate(s=>pst%s)

  if(s%r%nrestart>0)return
  if(s%r%verbose)write(*,*)'Entering init_part_ascii'

  if(s%r%part)then

     ! Compute total number of particles in file
     if(TRIM(s%r%initfile(s%r%levelmin)).NE.' ')then
        filename=TRIM(s%r%initfile(s%r%levelmin))//'/ic_part'
        inquire(file=filename, exist=file_exist)
        if(file_exist)then
           write(*,*)'Opening file '//TRIM(filename)
           open(10,file=filename,form='formatted')
           npart_tot=0
           do
#if NDIM==1
              read(10,*,end=101)xx1,vv1,mm1
              if(ABS(xx1)<s%r%box_size(1)/2.0d0)then
#endif
#if NDIM==2
              read(10,*,end=101)xx1,xx2,vv1,vv2,mm1
              if(ABS(xx1)<s%r%box_size(1)/2.0d0.AND.ABS(xx2)<s%r%box_size(2)/2.0d0)then
#endif
#if NDIM==3
              read(10,*,end=101)xx1,xx2,xx3,vv1,vv2,vv3,mm1
              if(ABS(xx1)<s%r%box_size(1)/2.0d0.AND.ABS(xx2)<s%r%box_size(2)/2.0d0.AND.ABS(xx3)<s%r%box_size(3)/2.0d0)then
#endif
                 npart_tot=npart_tot+1
              endif
           end do
101        continue
           s%p%npart_tot=npart_tot
           write(*,*)'Found npart_tot=',s%p%npart_tot
           close(10)
        else
           write(*,*)'File '//TRIM(filename)//' not found'
           s%p%npart_tot=0
        endif
     else
        s%p%npart_tot=0
     endif

     ! If no particle found, no need to read
     if(s%p%npart_tot>0)then
        ! Call recursive slave routine
        allocate(input_array(1:storage_size(npart_tot)/32))
        input_array=transfer(npart_tot,input_array)
        call r_input_part_ascii(pst,input_array,2)
        deallocate(input_array)
     endif

  endif

  if(s%r%star)then

     ! Compute total number of stars in file
     if(TRIM(s%r%initfile(s%r%levelmin)).NE.' ')then
        filename=TRIM(s%r%initfile(s%r%levelmin))//'/ic_star'
        inquire(file=filename, exist=file_exist)
        if(file_exist)then
           write(*,*)'Opening file '//TRIM(filename)
           open(10,file=filename,form='formatted')
           nstar_tot=0
           do
#if NDIM==1
              read(10,*,end=102)xx1,vv1,mm1,zz1,tt1
              if(ABS(xx1)<s%r%box_size(1)/2.0d0)then
#endif
#if NDIM==2
              read(10,*,end=102)xx1,xx2,vv1,vv2,mm1,zz1,tt1
              if(ABS(xx1)<s%r%box_size(1)/2.0d0.AND.ABS(xx2)<s%r%box_size(2)/2.0d0)then
#endif
#if NDIM==3
              read(10,*,end=102)xx1,xx2,xx3,vv1,vv2,vv3,mm1,zz1,tt1
              if(ABS(xx1)<s%r%box_size(1)/2.0d0.AND.ABS(xx2)<s%r%box_size(2)/2.0d0.AND.ABS(xx3)<s%r%box_size(3)/2.0d0)then
#endif
                 nstar_tot=nstar_tot+1
                 s%g%mass_star_tot=s%g%mass_star_tot+mm1
              endif
           end do
102        continue
           s%star%npart_tot=nstar_tot
           write(*,*)'Found nstar_tot=',s%star%npart_tot
           close(10)
        else
           write(*,*)'File '//TRIM(filename)//' not found'
           s%star%npart_tot=0
        endif
     else
        s%star%npart_tot=0
     endif

     ! If no particle found, no need to read
     if(s%star%npart_tot>0)then
        ! Call recursive slave routine
        allocate(input_array(1:storage_size(nstar_tot)/32))
        input_array=transfer(nstar_tot,input_array)
        call r_input_star_ascii(pst,input_array,2)
        deallocate(input_array)
     endif

  endif

  ! Load sink particles with the appropriate subroutine
  call m_input_sink_ascii(pst)

  if(s%r%core)then

     ! Compute total number of cores in file
     if(TRIM(s%r%initfile(s%r%levelmin)).NE.' ')then
        filename=TRIM(s%r%initfile(s%r%levelmin))//'/ic_core'
        inquire(file=filename, exist=file_exist)
        if(file_exist)then
           write(*,*)'Opening file '//TRIM(filename)
           open(10,file=filename,form='formatted')
           ncore_tot=0
           do
#if NDIM==1
              read(10,*,end=104)xx1,vv1,jj1,mm1,tt1
              if(ABS(xx1)<s%r%box_size(1)/2.0d0)then
#endif
#if NDIM==2
              read(10,*,end=104)xx1,xx2,vv1,vv2,jj1,jj2,mm1,tt1
              if(ABS(xx1)<s%r%box_size(1)/2.0d0.AND.ABS(xx2)<s%r%box_size(2)/2.0d0)then
#endif
#if NDIM==3
              read(10,*,end=104)xx1,xx2,xx3,vv1,vv2,vv3,jj1,jj2,jj3,mm1,tt1
              if(ABS(xx1)<s%r%box_size(1)/2.0d0.AND.ABS(xx2)<s%r%box_size(2)/2.0d0.AND.ABS(xx3)<s%r%box_size(3)/2.0d0)then
#endif
                 ncore_tot=ncore_tot+1
                 s%g%mass_sink_tot=s%g%mass_sink_tot+mm1
              endif
           end do
104        continue
           s%core%npart_tot=ncore_tot
           write(*,*)'Found ncore_tot=',s%core%npart_tot
           close(10)
        else
           write(*,*)'File '//TRIM(filename)//' not found'
           s%core%npart_tot=0
        endif
     else
        s%core%npart_tot=0
     endif

     ! If no particle found, no need to read
     if(s%core%npart_tot>0)then
        ! Call recursive slave routine
        allocate(input_array(1:storage_size(ncore_tot)/32))
        input_array=transfer(ncore_tot,input_array)
        call r_input_core_ascii(pst,input_array,2)
        deallocate(input_array)
     endif

  endif

  end associate

end subroutine m_input_part_ascii
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
recursive subroutine r_input_part_ascii(pst,input_array,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer,dimension(1:input_size)::input_array
  !--------------------------------------------------------------------
  ! This routine is the recursive slave procedure to read and dispatch
  ! particles from a Ramses restart file.
  !--------------------------------------------------------------------
  integer(kind=8)::npart_tot
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_INPUT_PART_ASCII,pst%iUpper+1,input_size,0,input_array)
     call r_input_part_ascii(pst%pLower,input_array,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     npart_tot=transfer(input_array,npart_tot)
     call input_part_ascii(pst%s%mdl,pst%s%r,pst%s%g,pst%s%p,npart_tot)
  endif

end subroutine r_input_part_ascii
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine input_part_ascii(mdl,r,g,p,npart_tot)
  use mdl_module
  use amr_parameters, only: ndim, flen
  use amr_commons, only: run_t,global_t
  use pm_commons, only: part_t
  implicit none
  type(mdl_t)::mdl
  type(run_t)::r
  type(global_t)::g
  type(part_t)::p
  integer(kind=8)::npart_tot
  !------------------------------------------------------------
  ! Allocate particle-based arrays.
  ! Read particles positions and velocities from various files
  ! including gadget, ascii or restart files.
  ! grafic initial conditions are performed after the AMR grid
  ! has been constructed.
  !------------------------------------------------------------
  integer::jpart_loc
  integer::i,ilun,icpu
  integer(kind=8)::indglob
  integer(kind=8)::jpart,npart,nremain
  integer(kind=8),dimension(1:g%ncpu+1)::start_ind
  real(kind=8)::xx1,xx2,xx3,vv1,vv2,vv3,mm1
  character(LEN=flen)::filename

  !--------------------------------------
  ! Compute starting index for each cpu
  !--------------------------------------
  p%npart_tot=npart_tot
  npart=npart_tot/g%ncpu
  nremain=npart_tot-int(npart,kind=8)*g%ncpu
  start_ind(1)=1
  do icpu=1,g%ncpu
     if(icpu.LE.nremain)then
        start_ind(icpu+1)=start_ind(icpu)+npart+1
     else
        start_ind(icpu+1)=start_ind(icpu)+npart
     endif
  end do

  !--------------------------------------
  ! Read ASCII initial conditions file
  !--------------------------------------
  filename=TRIM(r%initfile(r%levelmin))//'/ic_part'
  open(10,file=filename,form='formatted')
  jpart=0
  indglob=0
  jpart_loc=0
  do
#if NDIM==1
     read(10,*,end=100)xx1,vv1,mm1
     if(ABS(xx1)<r%box_size(1)/2.0d0)then
#endif
#if NDIM==2
     read(10,*,end=100)xx1,xx2,vv1,vv2,mm1
     if(ABS(xx1)<r%box_size(1)/2.0d0.AND.ABS(xx2)<r%box_size(2)/2.0d0)then
#endif
#if NDIM==3
     read(10,*,end=100)xx1,xx2,xx3,vv1,vv2,vv3,mm1
     if(ABS(xx1)<r%box_size(1)/2.0d0.AND.ABS(xx2)<r%box_size(2)/2.0d0.AND.ABS(xx3)<r%box_size(3)/2.0d0)then
#endif
        jpart=jpart+1
        indglob=indglob+1
        if(jpart >= start_ind(g%myid) .and. jpart < start_ind(g%myid+1))then
           jpart_loc=jpart_loc+1
           if(jpart_loc>r%npartmax)then
              write(*,*)'Maximum number of particles incorrect'
              write(*,*)'npartmax should be greater than',start_ind(2)
              call mdl_abort(mdl)
           endif
#if NDIM>0
           p%xp(jpart_loc,1)=xx1+r%box_size(1)/2.0
           p%vp(jpart_loc,1)=vv1
#endif
#if NDIM>1
           p%xp(jpart_loc,2)=xx2+r%box_size(2)/2.0
           p%vp(jpart_loc,2)=vv2
#endif
#if NDIM>2
           p%xp(jpart_loc,3)=xx3+r%box_size(3)/2.0
           p%vp(jpart_loc,3)=vv3
#endif
           p%mp(jpart_loc  )=mm1*r%ic_scale_m
           p%idp(jpart_loc )=indglob
           p%levelp(jpart_loc)=r%levelmin
        end if
     endif
  end do
100 continue
  close(10)

  p%npart=jpart_loc

  ! Put all particles in levelmin
  p%headp=p%npart+1
  p%tailp=p%npart
  p%headp(r%levelmin)=1
  p%tailp(r%levelmin)=p%npart
  if(ANY(.not.r%periodic(1:ndim)))then
     p%headp(r%levelmin-1)=1
     p%tailp(r%levelmin-1)=0
  endif

end subroutine input_part_ascii
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
recursive subroutine r_input_star_ascii(pst,input_array,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer,dimension(1:input_size)::input_array
  !--------------------------------------------------------------------
  ! This routine is the recursive slave procedure to read and dispatch
  ! particles from a Ramses restart file.
  !--------------------------------------------------------------------
  integer(kind=8)::npart_tot
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_INPUT_STAR_ASCII,pst%iUpper+1,input_size,0,input_array)
     call r_input_star_ascii(pst%pLower,input_array,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     npart_tot=transfer(input_array,npart_tot)
     call input_star_ascii(pst%s%mdl,pst%s%r,pst%s%g,pst%s%star,npart_tot)
  endif

end subroutine r_input_star_ascii
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine input_star_ascii(mdl,r,g,p,npart_tot)
  use mdl_module
  use amr_parameters, only: ndim, flen
  use amr_commons, only: run_t,global_t
  use pm_commons, only: part_t
  implicit none
  type(mdl_t)::mdl
  type(run_t)::r
  type(global_t)::g
  type(part_t)::p
  integer(kind=8)::npart_tot
  !------------------------------------------------------------
  ! Allocate particle-based arrays.
  ! Read particles positions and velocities from various files
  ! including gadget, ascii or restart files.
  ! grafic initial conditions are performed after the AMR grid
  ! has been constructed.
  !------------------------------------------------------------
  integer::jpart_loc
  integer::i,ilun,icpu
  integer(kind=8)::indglob
  integer(kind=8)::jpart,npart,nremain
  integer(kind=8),dimension(1:g%ncpu+1)::start_ind
  real(kind=8)::xx1,xx2,xx3,vv1,vv2,vv3,mm1,zz1,tt1
  character(LEN=flen)::filename

  !--------------------------------------
  ! Compute starting index for each cpu
  !--------------------------------------
  p%npart_tot=npart_tot
  npart=npart_tot/g%ncpu
  nremain=npart_tot-int(npart,kind=8)*g%ncpu
  start_ind(1)=1
  do icpu=1,g%ncpu
     if(icpu.LE.nremain)then
        start_ind(icpu+1)=start_ind(icpu)+npart+1
     else
        start_ind(icpu+1)=start_ind(icpu)+npart
     endif
  end do

  !--------------------------------------
  ! Read ASCII initial conditions file
  !--------------------------------------
  filename=TRIM(r%initfile(r%levelmin))//'/ic_star'
  open(10,file=filename,form='formatted')
  jpart=0
  indglob=0
  jpart_loc=0
  do
#if NDIM==1
     read(10,*,end=100)xx1,vv1,mm1,zz1,tt1
     if(ABS(xx1)<r%box_size(1)/2.0d0)then
#endif
#if NDIM==2
     read(10,*,end=100)xx1,xx2,vv1,vv2,mm1,zz1,tt1
     if(ABS(xx1)<r%box_size(1)/2.0d0.AND.ABS(xx2)<r%box_size(2)/2.0d0)then
#endif
#if NDIM==3
     read(10,*,end=100)xx1,xx2,xx3,vv1,vv2,vv3,mm1,zz1,tt1
     if(ABS(xx1)<r%box_size(1)/2.0d0.AND.ABS(xx2)<r%box_size(2)/2.0d0.AND.ABS(xx3)<r%box_size(3)/2.0d0)then
#endif
        jpart=jpart+1
        indglob=indglob+1
        if(jpart >= start_ind(g%myid) .and. jpart < start_ind(g%myid+1))then
           jpart_loc=jpart_loc+1
           if(jpart_loc>r%nstarmax)then
              write(*,*)'Maximum number of star particles incorrect'
              write(*,*)'nstarmax should be greater than',start_ind(2)
              call mdl_abort(mdl)
           endif
#if NDIM>0
           p%xp(jpart_loc,1)=xx1+r%box_size(1)/2.0
           p%vp(jpart_loc,1)=vv1
#endif
#if NDIM>1
           p%xp(jpart_loc,2)=xx2+r%box_size(2)/2.0
           p%vp(jpart_loc,2)=vv2
#endif
#if NDIM>2
           p%xp(jpart_loc,3)=xx3+r%box_size(3)/2.0
           p%vp(jpart_loc,3)=vv3
#endif
           p%mp(jpart_loc  )=mm1
           p%zp(jpart_loc  )=zz1
           p%tp(jpart_loc  )=tt1
           p%idp(jpart_loc )=indglob
           p%levelp(jpart_loc)=r%levelmin
        end if
     endif
  end do
100 continue
  close(10)

  p%npart=jpart_loc

  ! Put all particles in levelmin
  p%headp=p%npart+1
  p%tailp=p%npart
  p%headp(r%levelmin)=1
  p%tailp(r%levelmin)=p%npart
  if(ANY(.not.r%periodic(1:ndim)))then
     p%headp(r%levelmin-1)=1
     p%tailp(r%levelmin-1)=0
  endif

end subroutine input_star_ascii
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_input_sink_ascii(pst)
  use ramses_commons, only: pst_t
  use amr_parameters, only: flen
  implicit none
  type(pst_t)::pst
  !--------------------------------------------------------------------
  ! This routine is the master procedure to read and dispatch sink
  ! particles from an ascii file.
  ! --------------------------------------------------------------------
  real(kind=8)::xx1,xx2,xx3,vv1,vv2,vv3,mm1,tt1,jj1,jj2,jj3
  integer(kind=8)::nsink_tot
  character(LEN=flen)::filename
  integer,allocatable,dimension(:)::input_array
  logical::file_exist

  associate(s=>pst%s)

  if(s%r%sink)then
     ! Compute total number of sinks in file
     if(TRIM(s%r%initfile(s%r%levelmin)).NE.' ')then
        filename=TRIM(s%r%initfile(s%r%levelmin))//'/ic_sink'
        inquire(file=filename, exist=file_exist)
        if(file_exist)then
           write(*,*)'Opening file '//TRIM(filename)
           open(10,file=filename,form='formatted')
           nsink_tot=0
           do
#if NDIM==1
              read(10,*,end=103)xx1,vv1,jj1,mm1,tt1
              if(ABS(xx1)<s%r%box_size(1)/2.0d0)then
#endif
#if NDIM==2
              read(10,*,end=103)xx1,xx2,vv1,vv2,jj1,jj2,mm1,tt1
              if(ABS(xx1)<s%r%box_size(1)/2.0d0.AND.ABS(xx2)<s%r%box_size(2)/2.0d0)then
#endif
#if NDIM==3
              read(10,*,end=103)xx1,xx2,xx3,vv1,vv2,vv3,jj1,jj2,jj3,mm1,tt1
              if(ABS(xx1)<s%r%box_size(1)/2.0d0.AND.ABS(xx2)<s%r%box_size(2)/2.0d0.AND.ABS(xx3)<s%r%box_size(3)/2.0d0)then
#endif
                 nsink_tot=nsink_tot+1
                 s%g%mass_sink_tot=s%g%mass_sink_tot+mm1
              endif
           end do
103        continue
           s%sink%npart_tot=nsink_tot
           write(*,*)'Found nsink_tot=',s%sink%npart_tot
           close(10)
        else
           write(*,*)'File '//TRIM(filename)//' not found'
           s%sink%npart_tot=0
        endif
     else
        s%sink%npart_tot=0
     endif

     ! If no particle found, no need to read
     if(s%sink%npart_tot>0)then
        ! Call recursive slave routine
        allocate(input_array(1:storage_size(nsink_tot)/32))
        input_array=transfer(nsink_tot,input_array)
        call r_input_sink_ascii(pst,input_array,2)
        deallocate(input_array)
     endif
  endif

  end associate

end subroutine m_input_sink_ascii
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
recursive subroutine r_input_sink_ascii(pst,input_array,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer,dimension(1:input_size)::input_array
  !--------------------------------------------------------------------
  ! This routine is the recursive slave procedure to read and dispatch
  ! particles from a Ramses restart file.
  !--------------------------------------------------------------------
  integer(kind=8)::npart_tot
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_INPUT_SINK_ASCII,pst%iUpper+1,input_size,0,input_array)
     call r_input_sink_ascii(pst%pLower,input_array,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     npart_tot=transfer(input_array,npart_tot)
     call input_sink_ascii(pst%s%mdl,pst%s%r,pst%s%g,pst%s%sink,npart_tot)
  endif

end subroutine r_input_sink_ascii
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine input_sink_ascii(mdl,r,g,p,npart_tot)
  use mdl_module
  use amr_parameters, only: ndim, flen
  use amr_commons, only: run_t,global_t
  use pm_commons, only: part_t
  implicit none
  type(mdl_t)::mdl
  type(run_t)::r
  type(global_t)::g
  type(part_t)::p
  integer(kind=8)::npart_tot
  !------------------------------------------------------------
  ! Allocate particle-based arrays.
  ! Read particles positions and velocities from various files
  ! including gadget, ascii or restart files.
  ! grafic initial conditions are performed after the AMR grid
  ! has been constructed.
  !------------------------------------------------------------
  integer::jpart_loc
  integer::i,ilun,icpu
  integer(kind=8)::indglob
  integer(kind=8)::jpart,npart,nremain
  integer(kind=8),dimension(1:g%ncpu+1)::start_ind
  real(kind=8)::xx1,xx2,xx3,vv1,vv2,vv3,mm1,zz1,tt1,jj1,jj2,jj3
  character(LEN=flen)::filename

  !--------------------------------------
  ! Compute starting index for each cpu
  !--------------------------------------
  p%npart_tot=npart_tot
  npart=npart_tot/g%ncpu
  nremain=npart_tot-int(npart,kind=8)*g%ncpu
  start_ind(1)=1
  do icpu=1,g%ncpu
     if(icpu.LE.nremain)then
        start_ind(icpu+1)=start_ind(icpu)+npart+1
     else
        start_ind(icpu+1)=start_ind(icpu)+npart
     endif
  end do

  !--------------------------------------
  ! Read ASCII initial conditions file
  !--------------------------------------
  filename=TRIM(r%initfile(r%levelmin))//'/ic_sink'
  open(10,file=filename,form='formatted')
  jpart=0
  indglob=0
  jpart_loc=0
  do
#if NDIM==1
     read(10,*,end=100)xx1,vv1,jj1,mm1,tt1
     if(ABS(xx1)<r%box_size(1)/2.0d0)then
#endif
#if NDIM==2
     read(10,*,end=100)xx1,xx2,vv1,vv2,jj1,jj2,mm1,tt1
     if(ABS(xx1)<r%box_size(1)/2.0d0.AND.ABS(xx2)<r%box_size(2)/2.0d0)then
#endif
#if NDIM==3
     read(10,*,end=100)xx1,xx2,xx3,vv1,vv2,vv3,jj1,jj2,jj3,mm1,tt1
     if(ABS(xx1)<r%box_size(1)/2.0d0.AND.ABS(xx2)<r%box_size(2)/2.0d0.AND.ABS(xx3)<r%box_size(3)/2.0d0)then
#endif
        jpart=jpart+1
        indglob=indglob+1
        if(jpart >= start_ind(g%myid) .and. jpart < start_ind(g%myid+1))then
           jpart_loc=jpart_loc+1
           if(jpart_loc>r%nsinkmax)then
              write(*,*)'Maximum number of sink particles incorrect'
              write(*,*)'nsinkmax should be greater than',start_ind(2)
              call mdl_abort(mdl)
           endif
#if NDIM>0
           p%xp(jpart_loc,1)=xx1+r%box_size(1)/2.0
           p%vp(jpart_loc,1)=vv1
           p%jp(jpart_loc,1)=jj1
#endif
#if NDIM>1
           p%xp(jpart_loc,2)=xx2+r%box_size(2)/2.0
           p%vp(jpart_loc,2)=vv2
           p%jp(jpart_loc,2)=jj2
#endif
#if NDIM>2
           p%xp(jpart_loc,3)=xx3+r%box_size(3)/2.0
           p%vp(jpart_loc,3)=vv3
           p%jp(jpart_loc,3)=jj3
#endif
           p%mp(jpart_loc  )=mm1
           p%tp(jpart_loc  )=tt1
           p%idp(jpart_loc )=indglob
           p%levelp(jpart_loc)=r%levelmin
        end if
     endif
  end do
100 continue
  close(10)

  p%npart=jpart_loc

  ! Put all particles in levelmin
  p%headp=p%npart+1
  p%tailp=p%npart
  p%headp(r%levelmin)=1
  p%tailp(r%levelmin)=p%npart
  if(ANY(.not.r%periodic(1:ndim)))then
     p%headp(r%levelmin-1)=1
     p%tailp(r%levelmin-1)=0
  endif

  ! Fields the file does not give: zero acceleration, and the merging
  ! time and id of a new sink
  if(allocated(p%fp))p%fp(1:p%npart,1:ndim)=0d0
  if(allocated(p%tm))p%tm(1:p%npart)=-1000d0
  if(allocated(p%idm))p%idm(1:p%npart)=0

end subroutine input_sink_ascii
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
recursive subroutine r_input_core_ascii(pst,input_array,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer,dimension(1:input_size)::input_array
  !--------------------------------------------------------------------
  ! This routine is the recursive slave procedure to read and dispatch
  ! cores from an ascii file.
  !--------------------------------------------------------------------
  integer(kind=8)::npart_tot
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_INPUT_CORE_ASCII,pst%iUpper+1,input_size,0,input_array)
     call r_input_core_ascii(pst%pLower,input_array,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     npart_tot=transfer(input_array,npart_tot)
     call input_core_ascii(pst%s%mdl,pst%s%r,pst%s%g,pst%s%core,npart_tot)
  endif

end subroutine r_input_core_ascii
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine input_core_ascii(mdl,r,g,p,npart_tot)
  use mdl_module
  use amr_parameters, only: ndim, flen
  use amr_commons, only: run_t,global_t
  use pm_commons, only: part_t
  implicit none
  type(mdl_t)::mdl
  type(run_t)::r
  type(global_t)::g
  type(part_t)::p
  integer(kind=8)::npart_tot
  !------------------------------------------------------------
  ! Allocate particle-based arrays.
  ! Read particles positions and velocities from various files
  ! including gadget, ascii or restart files.
  ! grafic initial conditions are performed after the AMR grid
  ! has been constructed.
  !------------------------------------------------------------
  integer::jpart_loc
  integer::i,ilun,icpu
  integer(kind=8)::indglob
  integer(kind=8)::jpart,npart,nremain
  integer(kind=8),dimension(1:g%ncpu+1)::start_ind
  real(kind=8)::xx1,xx2,xx3,vv1,vv2,vv3,mm1,zz1,tt1,jj1,jj2,jj3
  character(LEN=flen)::filename

  !--------------------------------------
  ! Compute starting index for each cpu
  !--------------------------------------
  p%npart_tot=npart_tot
  npart=npart_tot/g%ncpu
  nremain=npart_tot-int(npart,kind=8)*g%ncpu
  start_ind(1)=1
  do icpu=1,g%ncpu
     if(icpu.LE.nremain)then
        start_ind(icpu+1)=start_ind(icpu)+npart+1
     else
        start_ind(icpu+1)=start_ind(icpu)+npart
     endif
  end do

  !--------------------------------------
  ! Read ASCII initial conditions file
  !--------------------------------------
  filename=TRIM(r%initfile(r%levelmin))//'/ic_core'
  open(10,file=filename,form='formatted')
  jpart=0
  indglob=0
  jpart_loc=0
  do
#if NDIM==1
     read(10,*,end=100)xx1,vv1,jj1,mm1,tt1
     if(ABS(xx1)<r%box_size(1)/2.0d0)then
#endif
#if NDIM==2
     read(10,*,end=100)xx1,xx2,vv1,vv2,jj1,jj2,mm1,tt1
     if(ABS(xx1)<r%box_size(1)/2.0d0.AND.ABS(xx2)<r%box_size(2)/2.0d0)then
#endif
#if NDIM==3
     read(10,*,end=100)xx1,xx2,xx3,vv1,vv2,vv3,jj1,jj2,jj3,mm1,tt1
     if(ABS(xx1)<r%box_size(1)/2.0d0.AND.ABS(xx2)<r%box_size(2)/2.0d0.AND.ABS(xx3)<r%box_size(3)/2.0d0)then
#endif
        jpart=jpart+1
        indglob=indglob+1
        if(jpart >= start_ind(g%myid) .and. jpart < start_ind(g%myid+1))then
           jpart_loc=jpart_loc+1
           if(jpart_loc>r%ncoremax)then
              write(*,*)'Maximum number of core particles incorrect'
              write(*,*)'ncoremax should be greater than',start_ind(2)
              call mdl_abort(mdl)
           endif
#if NDIM>0
           p%xp(jpart_loc,1)=xx1+r%box_size(1)/2.0
           p%vp(jpart_loc,1)=vv1
           p%jp(jpart_loc,1)=jj1
#endif
#if NDIM>1
           p%xp(jpart_loc,2)=xx2+r%box_size(2)/2.0
           p%vp(jpart_loc,2)=vv2
           p%jp(jpart_loc,2)=jj2
#endif
#if NDIM>2
           p%xp(jpart_loc,3)=xx3+r%box_size(3)/2.0
           p%vp(jpart_loc,3)=vv3
           p%jp(jpart_loc,3)=jj3
#endif
           p%mp(jpart_loc  )=mm1
           p%tp(jpart_loc  )=tt1
           p%idp(jpart_loc )=indglob
           p%levelp(jpart_loc)=r%levelmin
        end if
     endif
  end do
100 continue
  close(10)

  p%npart=jpart_loc

  ! Put all particles in levelmin
  p%headp=p%npart+1
  p%tailp=p%npart
  p%headp(r%levelmin)=1
  p%tailp(r%levelmin)=p%npart
  if(ANY(.not.r%periodic(1:ndim)))then
     p%headp(r%levelmin-1)=1
     p%tailp(r%levelmin-1)=0
  endif

  ! Fields the file does not give: zero acceleration and metallicity, the
  ! internal energy of gas at the sound speed floor, and accretion counted
  ! from step 0
  if(allocated(p%fp))p%fp(1:p%npart,1:ndim)=0d0
  if(allocated(p%zp))p%zp(1:p%npart)=0d0
  if(allocated(p%up))p%up(1:p%npart)=r%smallc**2/(r%gamma-1d0)
  if(allocated(p%nacc))p%nacc(1:p%npart)=0

end subroutine input_core_ascii
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
end module input_part_ascii_module
