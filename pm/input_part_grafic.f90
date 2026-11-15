module input_part_grafic_module
contains
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_input_part_grafic(pst)
  use ramses_commons, only: pst_t
  implicit none
  type(pst_t)::pst
  !--------------------------------------------------------------------
  ! This routine is the master procedure to read and dispatch particles
  ! from a Ramses restart file.
  !--------------------------------------------------------------------
  integer,allocatable,dimension(:)::input_array

  associate(s=>pst%s)

  if(s%r%verbose)write(*,*)'Entering input_part_grafic'
  if(TRIM(s%r%filetype).NE.'grafic')return
  if(s%r%nrestart>0)return

  if(s%r%part)then
      ! Compute total number of particles in file
      if(TRIM(s%r%initfile(s%r%levelmin)).NE.' ')then
         s%p%npart_tot=s%g%n1(s%r%levelmin)*s%g%n2(s%r%levelmin)*s%g%n3(s%r%levelmin)
         write(*,'(A,I0)')' Found npart_tot= ',s%p%npart_tot
      else
         s%p%npart_tot=0
      endif

      ! If no particle found, no need to read
      if(s%p%npart_tot==0)then
         return
      endif

      ! Call recursive slave routine
      allocate(input_array(1:storage_size(s%p%npart_tot)/32))
      input_array=transfer(s%p%npart_tot,input_array)
      call r_input_part_grafic(pst,input_array,storage_size(s%p%npart_tot)/32)
      deallocate(input_array)
  endif

  if(s%r%trac)then
      ! Compute total number of particles in file
      if(TRIM(s%r%initfile(s%r%levelmin)).NE.' ')then
         s%trac%npart_tot=s%g%n1(s%r%levelmin)*s%g%n2(s%r%levelmin)*s%g%n3(s%r%levelmin)
         write(*,'(A,I0)')' Found ntrac_tot= ',s%trac%npart_tot
      else
         s%trac%npart_tot=0
      endif

      ! If no particle found, no need to read
      if(s%trac%npart_tot==0)then
         return
      endif

      ! Call recursive slave routine
      allocate(input_array(1:storage_size(s%trac%npart_tot)/32))
      input_array=transfer(s%trac%npart_tot,input_array)
      call r_input_trac_grafic(pst,input_array,storage_size(s%trac%npart_tot)/32)
      deallocate(input_array)
  endif

  if(s%r%dust)then
     ! Compute total number of particles in file
     if(TRIM(s%r%initfile(s%r%levelmin)).NE.' ')then
        s%dust%npart_tot=s%g%n1(s%r%levelmin)*s%g%n2(s%r%levelmin)*s%g%n3(s%r%levelmin)
        write(*,'(A,I0)')' Found ndust_tot= ',s%dust%npart_tot
     else
        s%dust%npart_tot=0
     endif

     ! If no particle found, no need to read
     if(s%dust%npart_tot==0)then
        return
     endif

     ! Call recursive slave routine
     allocate(input_array(1:storage_size(s%dust%npart_tot)/32))
     input_array=transfer(s%dust%npart_tot,input_array)
     call r_input_dust_grafic(pst,input_array,storage_size(s%dust%npart_tot)/32)
     deallocate(input_array)
  endif

  end associate

end subroutine m_input_part_grafic
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
recursive subroutine r_input_part_grafic(pst,input_array,input_size)
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
     rID = mdl_send_request(pst%s%mdl,MDL_INPUT_PART_GRAFIC,pst%iUpper+1,input_size,0,input_array)
     call r_input_part_grafic(pst%pLower,input_array,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     npart_tot=transfer(input_array,npart_tot)
     call input_part_grafic(pst%s%mdl,pst%s%r,pst%s%g,pst%s%p,npart_tot)
  endif

end subroutine r_input_part_grafic
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine input_part_grafic(mdl,r,g,p,npart_tot)
  use mdl_module, only: mdl_t, mdl_abort
  use amr_parameters, only: ndim, twotondim, flen
  use amr_commons, only: run_t, global_t
  use pm_commons, only: part_t
  implicit none
  type(mdl_t)::mdl
  type(run_t)::r
  type(global_t)::g
  type(part_t)::p
  integer(kind=8)::npart_tot
  !----------------------------------------------------
  ! Reading initial conditions from single GRAFIC file
  !----------------------------------------------------
  integer::icpu,ipart,idim
  integer::i1,i2,i3,i1_min,i1_max,i2_min,i2_max,i3_min,i3_max
  integer::plane_size
  real(kind=8)::dx,xx1,xx2,xx3
  real(kind=8)::dispmax=0.0
  integer,dimension(1:g%ncpu)::npart_loc
  integer(kind=8)::ipart_grafic
  integer(kind=8),dimension(1:g%ncpu+1)::start_ind
  real(kind=4),dimension(:,:),allocatable::init_plane,init_plane_x,init_plane_m
  character(LEN=flen)::filename,filename_x,filename_m
  character(LEN=5)::nchar
  logical::ok,error,keep_part,read_pos=.false.,read_mass=.false.

  !-------------------------------------------
  ! Mesh size at levelmin in normalised units
  !-------------------------------------------
  dx=0.5d0**r%levelmin

  !--------------------------------------
  ! Compute starting index for each cpu
  !--------------------------------------
  p%npart_tot=npart_tot
  do icpu=1,g%ncpu+1
     start_ind(icpu)=((icpu-1)*npart_tot)/g%ncpu
     if(icpu>1)then
        npart_loc(icpu-1)=start_ind(icpu)-start_ind(icpu-1)
     endif
  end do
  p%npart=npart_loc(g%myid)

  ! Check that local number of particles does not exceed maximum
  if(p%npart > r%npartmax)then
     write(*,*)'ERROR: CPU ',g%myid,' has too many particles: ',p%npart,' > ',r%npartmax
     call mdl_abort(mdl)
  endif

  !--------------------------------------
  ! Initialize particles in planes
  !--------------------------------------
  plane_size=g%n1(r%levelmin)*g%n2(r%levelmin)
  i3_min=start_ind(g%myid)/plane_size
  i3_max=(start_ind(g%myid+1)-1)/plane_size

  ipart=1
  ipart_grafic=i3_min*plane_size

  ! Initialize positions, masses and ids
  do i3=i3_min,i3_max
     do i2=0,g%n2(r%levelmin)-1
        do i1=0,g%n1(r%levelmin)-1
           keep_part=(ipart_grafic>=start_ind(g%myid).AND.ipart<=p%npart)
           if(keep_part)then
              xx1=(dble(i1)+0.5)/dble(g%n1(r%levelmin))
              xx2=(dble(i2)+0.5)/dble(g%n2(r%levelmin))
              xx3=(dble(i3)+0.5)/dble(g%n3(r%levelmin))
              if(ndim>0)p%xp(ipart,1)=xx1
              if(ndim>1)p%xp(ipart,2)=xx2
              if(ndim>2)p%xp(ipart,3)=xx3
              p%mp(ipart)=0.5d0**(3*r%levelmin)*(1.0d0-g%omega_b/g%omega_m)
              p%idp(ipart)=ipart_grafic+1
              ipart=ipart+1
           endif
           ipart_grafic=ipart_grafic+1
        end do
     end do
  end do

  !--------------------------------------
  ! Allocate temporary arrays
  !--------------------------------------
  allocate(init_plane(1:g%n1(r%levelmin),1:g%n2(r%levelmin)))

  filename_x=TRIM(r%initfile(r%levelmin))//'/ic_poscx'
  INQUIRE(file=filename_x,exist=ok)
  read_pos=.false.
  if(ok)then
     read_pos=.true.
     allocate(init_plane_x(1:g%n1(r%levelmin),1:g%n2(r%levelmin)))
  endif

  filename_m=TRIM(r%initfile(r%levelmin))//'/ic_massc'
  INQUIRE(file=filename_m,exist=ok)
  read_mass=.false.
  if(ok)then
     read_mass=.true.
     allocate(init_plane_m(1:g%n1(r%levelmin),1:g%n2(r%levelmin)))
  endif

  !--------------------------------------
  ! Loop over displacements dimensions
  !--------------------------------------
  do idim=1,ndim

     ! Read dark matter Zeldovich initial displacement field
     if(idim==1)filename=TRIM(r%initfile(r%levelmin))//'/ic_velcx'
     if(idim==2)filename=TRIM(r%initfile(r%levelmin))//'/ic_velcy'
     if(idim==3)filename=TRIM(r%initfile(r%levelmin))//'/ic_velcz'
     if(g%myid==1)write(*,*)'Reading '//TRIM(filename)
     open(10,file=filename,form='unformatted')
     rewind 10
     read(10) ! skip first line
     do i3=0,i3_min-1
        read(10) ! skip unnecessary planes
     end do

     ! If present, read higher order dark matter initial displacement field
     if(idim==1)filename_x=TRIM(r%initfile(r%levelmin))//'/ic_poscx'
     if(idim==2)filename_x=TRIM(r%initfile(r%levelmin))//'/ic_poscy'
     if(idim==3)filename_x=TRIM(r%initfile(r%levelmin))//'/ic_poscz'
     if(read_pos)then
        if(g%myid==1)write(*,*)'Reading '//TRIM(filename_x)
        open(11,file=filename_x,form='unformatted')
        rewind 11
        read(11) ! skip first line
        do i3=0,i3_min-1
           read(11) ! skip unnecessary planes
        end do
     else
        if(g%myid==1)write(*,*)'Missing '//TRIM(filename_x)
     end if

     ! Read useful planes
     ipart=1
     ipart_grafic=i3_min*plane_size
     do i3=i3_min,i3_max
        read(10)((init_plane(i1,i2),i1=1,g%n1(r%levelmin)),i2=1,g%n2(r%levelmin))
        if(read_pos)then
           read(11)((init_plane_x(i1,i2),i1=1,g%n1(r%levelmin)),i2=1,g%n2(r%levelmin))
        endif
        ! Rescale initial displacement field to code units
        init_plane = init_plane*g%dfact(r%levelmin)*dx/g%dxini(r%levelmin)/g%vfact(r%levelmin)
        if(read_pos)then
           init_plane_x = init_plane_x/g%boxlen_ini
        endif
        do i2=1,g%n2(r%levelmin)
           do i1=1,g%n1(r%levelmin)
              keep_part=(ipart_grafic>=start_ind(g%myid).AND.ipart<=p%npart)
              if(keep_part)then
                 p%vp(ipart,idim)=init_plane(i1,i2)
                 if(.not. read_pos)then
                    dispmax=max(dispmax,abs(init_plane(i1,i2)/dx))
                 else
                    p%xp(ipart,idim)=p%xp(ipart,idim)+init_plane_x(i1,i2)
                    dispmax=max(dispmax,abs(init_plane_x(i1,i2)/dx))
                 endif
                 ipart=ipart+1
              endif
              ipart_grafic=ipart_grafic+1
           end do
        end do
     end do
     ! End loop over planes
     close(10)
     if(read_pos)close(11)

  end do
  ! End loop over dimensions

  ! Deallocate temporary array
  deallocate(init_plane)
  if(read_pos)deallocate(init_plane_x)

  ! Move particle according to Zeldovich approximation
  if(.not. read_pos)then
     do ipart=1,p%npart
        p%xp(ipart,1:ndim)=p%xp(ipart,1:ndim)+p%vp(ipart,1:ndim)
     enddo
  endif

  ! Scale displacement to velocity
  do ipart=1,p%npart
     p%vp(ipart,1:ndim)=g%vfact(1)*p%vp(ipart,1:ndim)
  end do

  ! Periodic box
  do ipart=1,p%npart
     do idim=1,ndim
        if(r%periodic(idim))then
           if(p%xp(ipart,idim)< 0.0d0           )p%xp(ipart,idim)=p%xp(ipart,idim)+r%box_size(idim)
           if(p%xp(ipart,idim)>=r%box_size(idim))p%xp(ipart,idim)=p%xp(ipart,idim)-r%box_size(idim)
        end if
     end do
  end do

  ! If mass file present, read dark matter particle masses
  if(read_mass)then
     if(g%myid==1)write(*,*)'Reading '//TRIM(filename_m)
     open(10,file=filename_m,form='unformatted')
     rewind 10
     read(10) ! skip first line
     do i3=0,i3_min-1
        read(10) ! skip unnecessary planes
     end do
     ipart=1
     ipart_grafic=i3_min*plane_size
     ! Loop over planes
     do i3=i3_min,i3_max
        read(10)((init_plane_m(i1,i2),i1=1,g%n1(r%levelmin)),i2=1,g%n2(r%levelmin))
        do i2=1,g%n2(r%levelmin)
           do i1=1,g%n1(r%levelmin)
              keep_part=(ipart_grafic>=start_ind(g%myid).AND.ipart<=p%npart)
              if(keep_part)then
                 p%mp(ipart)=0.5d0**(3*r%levelmin)*init_plane_m(i1,i2)
                 ipart=ipart+1
              endif
              ipart_grafic=ipart_grafic+1
           end do
        end do
     end do
     ! End loop over planes
     close(10)
     deallocate(init_plane_m)
  else
     if(g%myid==1)write(*,*)'Missing '//TRIM(filename_m)
  endif

  ! Compute particle initial level
  do ipart=1,p%npart
     p%levelp(ipart)=r%levelmin
  end do

  ! Put all particles in levelmin
  p%headp=p%npart+1
  p%tailp=p%npart
  p%headp(r%levelmin)=1
  p%tailp(r%levelmin)=p%npart
  if(ANY(.not.r%periodic(1:ndim)))then
     p%headp(r%levelmin-1)=1
     p%tailp(r%levelmin-1)=0
  endif

end subroutine input_part_grafic
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
recursive subroutine r_input_trac_grafic(pst,input_array,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer,dimension(1:input_size)::input_array
  !--------------------------------------------------------------------
  ! This routine is the recursive slave procedure to read and dispatch
  ! particles from a grafic file.
  !--------------------------------------------------------------------
  integer(kind=8)::npart_tot
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_INPUT_TRAC_GRAFIC,pst%iUpper+1,input_size,0,input_array)
     call r_input_trac_grafic(pst%pLower,input_array,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     npart_tot=transfer(input_array,npart_tot)
     call input_trac_grafic(pst%s%mdl,pst%s%r,pst%s%g,pst%s%trac,npart_tot)
  endif

end subroutine r_input_trac_grafic
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine input_trac_grafic(mdl,r,g,p,npart_tot)
   use mdl_module, only: mdl_t, mdl_abort
   use amr_parameters, only: ndim, twotondim, flen
   use amr_commons, only: run_t, global_t
   use pm_commons, only: part_t
   implicit none
   type(mdl_t)::mdl
   type(run_t)::r
   type(global_t)::g
   type(part_t)::p
   integer(kind=8)::npart_tot
   !----------------------------------------------------
   ! Reading initial conditions from single GRAFIC file
   !----------------------------------------------------
   integer::icpu,ipart,idim,ipercell,n_per_cell
   integer::i1,i2,i3,i1_min,i1_max,i2_min,i2_max,i3_min,i3_max
   integer::plane_size
   real(kind=8)::dx,xx1,xx2,xx3
   real(kind=8)::dispmax=0.0
   integer,dimension(1:g%ncpu)::npart_loc
   integer(kind=8)::ipart_grafic
   integer(kind=8),dimension(1:g%ncpu+1)::start_ind
   real(kind=4),dimension(:,:),allocatable::init_plane,init_plane_x
   real(kind=8),dimension(:),allocatable::tdx,tdy,tdz
   character(LEN=flen)::filename,filename_x
   character(LEN=5)::nchar
   logical::ok,error,keep_part,read_pos=.false.

   !-------------------------------------------
   ! Mesh size at levelmin in normalised units
   !-------------------------------------------
   dx=0.5d0**r%levelmin

   !--------------------------------------
   ! Compute starting index for each cpu
   !--------------------------------------
   p%npart_tot=npart_tot
   do icpu=1,g%ncpu+1
      start_ind(icpu)=((icpu-1)*npart_tot)/g%ncpu
      if(icpu>1)then
         npart_loc(icpu-1)=start_ind(icpu)-start_ind(icpu-1)
      endif
   end do
   ! Number of tracers per cell
   n_per_cell = max(1,r%ntrac_per_cell)
   p%npart=n_per_cell*npart_loc(g%myid)

   ! Check that local number of tracer particles does not exceed maximum
   if(p%npart > r%ntracmax)then
      write(*,*)'ERROR: CPU ',g%myid,' has too many tracer particles: ',p%npart,' > ',r%ntracmax
      call mdl_abort(mdl)
   endif

   !--------------------------------------
   ! Initialize particles in planes
   !--------------------------------------
   plane_size=g%n1(r%levelmin)*g%n2(r%levelmin)
   i3_min=start_ind(g%myid)/plane_size
   i3_max=(start_ind(g%myid+1)-1)/plane_size

   ! Precompute subcell offsets
   allocate(tdx(1:n_per_cell),tdy(1:n_per_cell),tdz(1:n_per_cell))
  if(r%part_subcell_positions)then
     call part_subcell_positions(n_per_cell,tdx,tdy,tdz)
  else
     tdx=0.d0 ; tdy=0.d0 ; tdz=0.d0
  endif

   ipart=1
   ipart_grafic=i3_min*plane_size

   ! Initialize positions, masses and ids
   do i3=i3_min,i3_max
      do i2=0,g%n2(r%levelmin)-1
         do i1=0,g%n1(r%levelmin)-1
           do ipercell=1,n_per_cell
              keep_part=(ipart_grafic>=start_ind(g%myid).AND.ipart<=p%npart)
              if(keep_part)then
                 xx1=(dble(i1)+0.5)/dble(g%n1(r%levelmin))
                 xx2=(dble(i2)+0.5)/dble(g%n2(r%levelmin))
                 xx3=(dble(i3)+0.5)/dble(g%n3(r%levelmin))
                 if(ndim>0)p%xp(ipart,1)=xx1+tdx(ipercell)*dx
                 if(ndim>1)p%xp(ipart,2)=xx2+tdy(ipercell)*dx
                 if(ndim>2)p%xp(ipart,3)=xx3+tdz(ipercell)*dx
                 p%mp(ipart)=0.5d0**(3*r%levelmin)/n_per_cell
                 p%idp(ipart)=ipart_grafic*n_per_cell + (ipercell-1) + 1
                 ipart=ipart+1
              endif
           end do
           ipart_grafic=ipart_grafic+1
         end do
      end do
   end do

   !--------------------------------------
   ! Allocate temporary arrays
   !--------------------------------------
   allocate(init_plane(1:g%n1(r%levelmin),1:g%n2(r%levelmin)))
   filename_x=TRIM(r%initfile(r%levelmin))//'/ic_poscx'
   INQUIRE(file=filename_x,exist=ok)
   read_pos=.false.
   if(ok)then
      read_pos=.true.
      allocate(init_plane_x(1:g%n1(r%levelmin),1:g%n2(r%levelmin)))
   else
      if(g%myid==1)write(*,*)'File '//TRIM(filename_x)//' not found.'
   endif

   !--------------------------------------
   ! Loop over spatial dimensions
   !--------------------------------------
   do idim=1,ndim

      ! Read particle initial velocities
      if(idim==1)filename=TRIM(r%initfile(r%levelmin))//'/ic_velcx'
      if(idim==2)filename=TRIM(r%initfile(r%levelmin))//'/ic_velcy'
      if(idim==3)filename=TRIM(r%initfile(r%levelmin))//'/ic_velcz'

      if(g%myid==1)write(*,*)'Reading file '//TRIM(filename)
      open(10,file=filename,form='unformatted')
      rewind 10
      read(10) ! skip first line
      do i3=0,i3_min-1
         read(10) ! skip unnecessary planes
      end do

      ! If present, read initial position displacement field
      if(read_pos)then
         if(idim==1)filename_x=TRIM(r%initfile(r%levelmin))//'/ic_poscx'
         if(idim==2)filename_x=TRIM(r%initfile(r%levelmin))//'/ic_poscy'
         if(idim==3)filename_x=TRIM(r%initfile(r%levelmin))//'/ic_poscz'

         if(g%myid==1)write(*,*)'Reading file '//TRIM(filename_x)
         open(11,file=filename_x,form='unformatted')
         rewind 11
         read(11) ! skip first line
         do i3=0,i3_min-1
            read(11) ! skip unnecessary planes
         end do
      end if

      ! Read useful planes
     ipart=1
     ipart_grafic=i3_min*plane_size
      do i3=i3_min,i3_max
         read(10)((init_plane(i1,i2),i1=1,g%n1(r%levelmin)),i2=1,g%n2(r%levelmin))
         if(read_pos)then
            read(11)((init_plane_x(i1,i2),i1=1,g%n1(r%levelmin)),i2=1,g%n2(r%levelmin))
         endif
         ! Rescale initial displacement field to code units
         init_plane = init_plane*g%dfact(r%levelmin)*dx/g%dxini(r%levelmin)/g%vfact(r%levelmin)
         if(read_pos)then
            init_plane_x = init_plane_x/g%boxlen_ini
         endif
        do i2=1,g%n2(r%levelmin)
           do i1=1,g%n1(r%levelmin)
              do ipercell=1,n_per_cell
                 keep_part=(ipart_grafic>=start_ind(g%myid).AND.ipart<=p%npart)
                 if(keep_part)then
                    p%vp(ipart,idim)=init_plane(i1,i2)
                    if(.not. read_pos)then
                       dispmax=max(dispmax,abs(init_plane(i1,i2)/dx))
                    else
                       p%xp(ipart,idim)=p%xp(ipart,idim)+init_plane_x(i1,i2)
                       dispmax=max(dispmax,abs(init_plane_x(i1,i2)/dx))
                    endif
                    ipart=ipart+1
                 endif
              end do
              ipart_grafic=ipart_grafic+1
           end do
        end do
      end do
      ! End loop over planes
      close(10)
      if(read_pos)close(11)

   end do
   ! End loop over dimensions

  ! Deallocate temporary array
  deallocate(init_plane)
  if(read_pos)deallocate(init_plane_x)
  if(allocated(tdx))deallocate(tdx,tdy,tdz)

  ! Periodic box
  do ipart=1,p%npart
     do idim=1,ndim
        if(r%periodic(idim))then
           if(p%xp(ipart,idim)< 0.0d0           )p%xp(ipart,idim)=p%xp(ipart,idim)+r%box_size(idim)
           if(p%xp(ipart,idim)>=r%box_size(idim))p%xp(ipart,idim)=p%xp(ipart,idim)-r%box_size(idim)
        end if
     end do
  end do

   ! Compute particle initial level
   do ipart=1,p%npart
      p%levelp(ipart)=r%levelmin
   end do

   ! Put all particles in levelmin
   p%headp=p%npart+1
   p%tailp=p%npart
   p%headp(r%levelmin)=1
   p%tailp(r%levelmin)=p%npart

end subroutine input_trac_grafic
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
recursive subroutine r_input_dust_grafic(pst,input_array,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer,dimension(1:input_size)::input_array
  !--------------------------------------------------------------------
  ! This routine is the recursive slave procedure to read and dispatch
  ! dust particles from a grafic file.
  !--------------------------------------------------------------------
  integer(kind=8)::npart_tot
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_INPUT_DUST_GRAFIC,pst%iUpper+1,input_size,0,input_array)
     call r_input_dust_grafic(pst%pLower,input_array,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     npart_tot=transfer(input_array,npart_tot)
     call input_dust_grafic(pst%s%mdl,pst%s%r,pst%s%g,pst%s%dust,npart_tot)
  endif

end subroutine r_input_dust_grafic
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine input_dust_grafic(mdl,r,g,p,npart_tot)
  use mdl_module, only: mdl_t, mdl_abort
  use amr_parameters, only: ndim, twotondim, flen
  use amr_commons, only: run_t, global_t
  use pm_commons, only: part_t
  implicit none
  type(mdl_t)::mdl
  type(run_t)::r
  type(global_t)::g
  type(part_t)::p
  integer(kind=8)::npart_tot
  !----------------------------------------------------
  ! Reading initial conditions from single GRAFIC file
  !----------------------------------------------------
  integer::icpu,ipart,idim,ipercell,n_per_cell
  integer::i1,i2,i3,i1_min,i1_max,i2_min,i2_max,i3_min,i3_max
  integer::plane_size
  real(kind=8)::dx,xx1,xx2,xx3
  real(kind=8)::dispmax=0.0
  integer,dimension(1:g%ncpu)::npart_loc
  integer(kind=8)::ipart_grafic
  integer(kind=8),dimension(1:g%ncpu+1)::start_ind
  real(kind=4),dimension(:,:),allocatable::init_plane,init_plane_x
  real(kind=8),dimension(:),allocatable::tdx,tdy,tdz
  character(LEN=flen)::filename,filename_x
  character(LEN=5)::nchar
  logical::ok,error,keep_part,read_pos=.false.

  !-------------------------------------------
  ! Mesh size at levelmin in normalised units
  !-------------------------------------------
  dx=0.5d0**r%levelmin

  !--------------------------------------
  ! Compute starting index for each cpu
  !--------------------------------------
  p%npart_tot=npart_tot
  do icpu=1,g%ncpu+1
     start_ind(icpu)=((icpu-1)*npart_tot)/g%ncpu
     if(icpu>1)then
        npart_loc(icpu-1)=start_ind(icpu)-start_ind(icpu-1)
     endif
  end do
  ! Number of dust per cell
  n_per_cell = MAX(1,r%ndust_per_cell)
  p%npart=n_per_cell*npart_loc(g%myid)

  ! Check that local number of dust particles does not exceed maximum
  if(p%npart > r%ndustmax)then
     write(*,*)'ERROR: CPU ',g%myid,' has too many dust particles: ',p%npart,' > ',r%ndustmax
     call mdl_abort(mdl)
  endif

  !--------------------------------------
  ! Initialize particles in planes
  !--------------------------------------
  plane_size=g%n1(r%levelmin)*g%n2(r%levelmin)
  i3_min=start_ind(g%myid)/plane_size
  i3_max=(start_ind(g%myid+1)-1)/plane_size

  ! Precompute subcell offsets
  allocate(tdx(1:n_per_cell),tdy(1:n_per_cell),tdz(1:n_per_cell))
  if(r%part_subcell_positions)then
     call part_subcell_positions(n_per_cell,tdx,tdy,tdz)
  else
     tdx=0.d0 ; tdy=0.d0 ; tdz=0.d0
  endif

  ipart=1
  ipart_grafic=i3_min*plane_size

  ! Initialize positions, masses and ids
  do i3=i3_min,i3_max
     do i2=0,g%n2(r%levelmin)-1
        do i1=0,g%n1(r%levelmin)-1
           do ipercell=1,n_per_cell
              keep_part=(ipart_grafic>=start_ind(g%myid).AND.ipart<=p%npart)
              if(keep_part)then
                 xx1=(dble(i1)+0.5)/dble(g%n1(r%levelmin))
                 xx2=(dble(i2)+0.5)/dble(g%n2(r%levelmin))
                 xx3=(dble(i3)+0.5)/dble(g%n3(r%levelmin))
                 if(ndim>0)p%xp(ipart,1)=xx1+tdx(ipercell)*dx
                 if(ndim>1)p%xp(ipart,2)=xx2+tdy(ipercell)*dx
                 if(ndim>2)p%xp(ipart,3)=xx3+tdz(ipercell)*dx
                 p%mp(ipart)=r%dust_to_gas_mass_ratio * 0.5d0**(3*r%levelmin)/n_per_cell
                 p%idp(ipart)=ipart_grafic*n_per_cell + (ipercell-1) + 1
                 p%size(ipart)=r%grain_size_parameter
                 p%charge(ipart)=r%grain_charge_parameter
                 ipart=ipart+1
              endif
           end do
           ipart_grafic=ipart_grafic+1
        end do
     end do
  end do

  !--------------------------------------
  ! Allocate temporary arrays
  !--------------------------------------
  allocate(init_plane(1:g%n1(r%levelmin),1:g%n2(r%levelmin)))
  filename_x=TRIM(r%initfile(r%levelmin))//'/ic_poscx'
  INQUIRE(file=filename_x,exist=ok)
  read_pos=.false.
  if(ok)then
     read_pos=.true.
     allocate(init_plane_x(1:g%n1(r%levelmin),1:g%n2(r%levelmin)))
  else
     if(g%myid==1)write(*,*)'File '//TRIM(filename_x)//' not found.'
  endif

  !--------------------------------------
  ! Loop over spatial dimensions
  !--------------------------------------
  do idim=1,ndim

     ! Read particle initial velocities
     if(idim==1)filename=TRIM(r%initfile(r%levelmin))//'/ic_velcx'
     if(idim==2)filename=TRIM(r%initfile(r%levelmin))//'/ic_velcy'
     if(idim==3)filename=TRIM(r%initfile(r%levelmin))//'/ic_velcz'

     if(g%myid==1)write(*,*)'Reading file '//TRIM(filename)
     open(10,file=filename,form='unformatted')
     rewind 10
     read(10) ! skip first line
     do i3=0,i3_min-1
        read(10) ! skip unnecessary planes
     end do

     ! If present, read initial position displacement field
     if(read_pos)then
        if(idim==1)filename_x=TRIM(r%initfile(r%levelmin))//'/ic_poscx'
        if(idim==2)filename_x=TRIM(r%initfile(r%levelmin))//'/ic_poscy'
        if(idim==3)filename_x=TRIM(r%initfile(r%levelmin))//'/ic_poscz'

        if(g%myid==1)write(*,*)'Reading file '//TRIM(filename_x)
        open(11,file=filename_x,form='unformatted')
        rewind 11
        read(11) ! skip first line
        do i3=0,i3_min-1
           read(11) ! skip unnecessary planes
        end do
     end if

     ! Read useful planes
     ipart=1
     ipart_grafic=i3_min*plane_size
     do i3=i3_min,i3_max
        read(10)((init_plane(i1,i2),i1=1,g%n1(r%levelmin)),i2=1,g%n2(r%levelmin))
        if(read_pos)then
           read(11)((init_plane_x(i1,i2),i1=1,g%n1(r%levelmin)),i2=1,g%n2(r%levelmin))
        endif
        ! Rescale initial displacement field to code units
        init_plane = init_plane*g%dfact(r%levelmin)*dx/g%dxini(r%levelmin)/g%vfact(r%levelmin)
        if(read_pos)then
           init_plane_x = init_plane_x/g%boxlen_ini
        endif
        do i2=1,g%n2(r%levelmin)
           do i1=1,g%n1(r%levelmin)
              do ipercell=1,n_per_cell
                 keep_part=(ipart_grafic>=start_ind(g%myid).AND.ipart<=p%npart)
                 if(keep_part)then
                    p%vp(ipart,idim)=init_plane(i1,i2)
                    if(.not. read_pos)then
                       dispmax=max(dispmax,abs(init_plane(i1,i2)/dx))
                    else
                       p%xp(ipart,idim)=p%xp(ipart,idim)+init_plane_x(i1,i2)
                       dispmax=max(dispmax,abs(init_plane_x(i1,i2)/dx))
                    endif
                    ipart=ipart+1
                 endif
              end do
              ipart_grafic=ipart_grafic+1
           end do
        end do
     end do
     ! End loop over planes
     close(10)
     if(read_pos)close(11)

  end do
  ! End loop over dimensions

  ! Deallocate temporary array
  deallocate(init_plane)
  if(read_pos)deallocate(init_plane_x)
  if(allocated(tdx))deallocate(tdx,tdy,tdz)

  ! Periodic box
  do ipart=1,p%npart
     do idim=1,ndim
        if(r%periodic(idim))then
           if(p%xp(ipart,idim)< 0.0d0           )p%xp(ipart,idim)=p%xp(ipart,idim)+r%box_size(idim)
           if(p%xp(ipart,idim)>=r%box_size(idim))p%xp(ipart,idim)=p%xp(ipart,idim)-r%box_size(idim)
        endif
     end do
  end do

  ! Compute particle initial level
  do ipart=1,p%npart
     p%levelp(ipart)=r%levelmin
  end do

  ! Put all particles in levelmin
  p%headp=p%npart+1
  p%tailp=p%npart
  p%headp(r%levelmin)=1
  p%tailp(r%levelmin)=p%npart

end subroutine input_dust_grafic

! Evenly distribute offsets within a cell for nt sub-particles
subroutine part_subcell_positions(nt,trx,try,trz)
  integer::nt,ipart,i
  real(kind=8)::phi3,avx,avy,avz
  real(kind=8),dimension(1:nt)::trx,try,trz
  real(kind=8),dimension(1:3),save::phi3a

  phi3 = 1.2207440846057596d0
  do i=1,3
    phi3a(i)=phi3**(-i)
  end do
  do ipart = 1,nt
    trx(ipart) = mod(0.5d0+ ipart*phi3a(1),1.0d0)
    try(ipart) = mod(0.5d0+ ipart*phi3a(2),1.0d0)
    trz(ipart) = mod(0.5d0+ ipart*phi3a(3),1.0d0)
  end do

  avx=sum(trx)/nt
  avy=sum(try)/nt
  avz=sum(trz)/nt
  do ipart=1,nt
    trx(ipart)=mod(trx(ipart)-avx+0.5d0,1.0d0)-0.5d0
    try(ipart)=mod(try(ipart)-avy+0.5d0,1.0d0)-0.5d0
    trz(ipart)=mod(trz(ipart)-avz+0.5d0,1.0d0)-0.5d0
  end do
end subroutine part_subcell_positions
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
end module input_part_grafic_module
