module input_part_ramses_module

  type :: out_input_star_t
     real(kind=8)::mass
  end type out_input_star_t

contains
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_input_part_ramses(pst)
  use mdl_module
  use amr_parameters, only: ndim, flen
  use ramses_commons, only: pst_t
  use output_amr_module, only: input_header
  use pm_parameters, only: PART_TYPE, STAR_TYPE, SINK_TYPE, TREE_TYPE
  use pm_parameters, only: CORE_TYPE
  implicit none
  type(pst_t)::pst
  !--------------------------------------------------------------------
  ! This routine is the master procedure to read and dispatch particles
  ! from a Ramses output file.
  !--------------------------------------------------------------------
  integer::icpu,ilun,ncpu_file
  integer(kind=8)::npart_tot_file,npart_tot_check
  character(LEN=5)::ncharcpu
  character(LEN=flen)::file_head,file_part
  integer,allocatable,dimension(:)::npart_file
  type(out_input_star_t)::output

  if(pst%s%r%verbose)write(*,*)'Entering input_part_ramses'

  if(pst%s%r%part)then

     ! Read particle files header
     file_head=TRIM(pst%s%r%initfile(pst%s%r%levelmin))//'/part_header.txt'
     call input_header(pst%s%r,pst%s%g,file_head,npart_tot_file,ncpu_file)
     write(*,'(" Output folder has ",I12," DM particles")')npart_tot_file

     ! Allocate local array
     allocate(npart_file(0:ncpu_file))

     ! Read number of particles in each file
     npart_tot_check=0
     npart_file(0)=PART_TYPE
     do icpu=1,ncpu_file
        call title(icpu,ncharcpu)
        file_part=TRIM(pst%s%r%initfile(pst%s%r%levelmin))//'/part.'//TRIM(ncharcpu)
        ilun=10
        open(unit=ilun,file=TRIM(file_part),access="stream",action="read",form='unformatted')
        read(ilun,POS=5)npart_file(icpu)
        npart_tot_check=npart_tot_check+npart_file(icpu)
        close(ilun)
     end do
     if(npart_tot_check.NE.npart_tot_file)then
        write(*,*)' Input file corrupted'
        call mdl_abort(pst%s%mdl)
     endif

     ! Call recursive slave routine
     call r_input_part_ramses(pst,npart_file,ncpu_file+1,output,2)
     write(*,*)'Total mass in dark matter=',output%mass

     ! Deallocate local array
     deallocate(npart_file)

  endif

  if(pst%s%r%star)then

     ! Read star particle files header
     file_head=TRIM(pst%s%r%initfile(pst%s%r%levelmin))//'/star_header.txt'
     call input_header(pst%s%r,pst%s%g,file_head,npart_tot_file,ncpu_file)
     write(*,'(" Output folder has ",I12," star particles")')npart_tot_file

     ! Allocate local array
     allocate(npart_file(0:ncpu_file))

     ! Read number of particles in each file
     npart_tot_check=0
     npart_file(0)=STAR_TYPE
     do icpu=1,ncpu_file
        call title(icpu,ncharcpu)
        file_part=TRIM(pst%s%r%initfile(pst%s%r%levelmin))//'/star.'//TRIM(ncharcpu)
        ilun=10
        open(unit=ilun,file=TRIM(file_part),access="stream",action="read",form='unformatted')
        read(ilun,POS=5)npart_file(icpu)
        npart_tot_check=npart_tot_check+npart_file(icpu)
        close(ilun)
     end do
     if(npart_tot_check.NE.npart_tot_file)then
        write(*,*)' Input file corrupted'
        call mdl_abort(pst%s%mdl)
     endif

     ! Call recursive slave routine
     call r_input_part_ramses(pst,npart_file,ncpu_file+1,output,2)
     write(*,*)'Total mass in stars=',output%mass
     pst%s%g%mass_star_tot=output%mass

     ! Deallocate local array
     deallocate(npart_file)

  endif

  if(pst%s%r%sink)then

     ! Read star particle files header
     file_head=TRIM(pst%s%r%initfile(pst%s%r%levelmin))//'/sink_header.txt'
     call input_header(pst%s%r,pst%s%g,file_head,npart_tot_file,ncpu_file)
     write(*,'(" Output folder has ",I12," sink particles")')npart_tot_file

     ! Allocate local array
     allocate(npart_file(0:ncpu_file))

     ! Read number of particles in each file
     npart_tot_check=0
     npart_file(0)=SINK_TYPE
     do icpu=1,ncpu_file
        call title(icpu,ncharcpu)
        file_part=TRIM(pst%s%r%initfile(pst%s%r%levelmin))//'/sink.'//TRIM(ncharcpu)
        ilun=10
        open(unit=ilun,file=TRIM(file_part),access="stream",action="read",form='unformatted')
        read(ilun,POS=5)npart_file(icpu)
        npart_tot_check=npart_tot_check+npart_file(icpu)
        close(ilun)
     end do
     if(npart_tot_check.NE.npart_tot_file)then
        write(*,*)' Input file corrupted'
        call mdl_abort(pst%s%mdl)
     endif

     ! Call recursive slave routine
     call r_input_part_ramses(pst,npart_file,ncpu_file+1,output,2)
     write(*,*)'Total mass in sinks=',output%mass
     pst%s%g%mass_sink_tot=output%mass

     ! Deallocate local array
     deallocate(npart_file)

  endif

  if(pst%s%r%core)then

     ! Read core particle files header
     file_head=TRIM(pst%s%r%initfile(pst%s%r%levelmin))//'/core_header.txt'
     call input_header(pst%s%r,pst%s%g,file_head,npart_tot_file,ncpu_file)
     write(*,'(" Output folder has ",I12," core particles")')npart_tot_file

     ! Allocate local array
     allocate(npart_file(0:ncpu_file))

     ! Read number of particles in each file
     npart_tot_check=0
     npart_file(0)=CORE_TYPE
     do icpu=1,ncpu_file
        call title(icpu,ncharcpu)
        file_part=TRIM(pst%s%r%initfile(pst%s%r%levelmin))//'/core.'//TRIM(ncharcpu)
        ilun=10
        open(unit=ilun,file=TRIM(file_part),access="stream",action="read",form='unformatted')
        read(ilun,POS=5)npart_file(icpu)
        npart_tot_check=npart_tot_check+npart_file(icpu)
        close(ilun)
     end do
     if(npart_tot_check.NE.npart_tot_file)then
        write(*,*)' Input file corrupted'
        call mdl_abort(pst%s%mdl)
     endif

     ! Call recursive slave routine
     call r_input_part_ramses(pst,npart_file,ncpu_file+1,output,2)
     write(*,*)'Total mass in cores=',output%mass
     pst%s%g%mass_sink_tot=pst%s%g%mass_sink_tot+output%mass

     ! Deallocate local array
     deallocate(npart_file)

  endif

end subroutine m_input_part_ramses
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
recursive subroutine r_input_part_ramses(pst,input_array,input_size,output,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use pm_parameters, only: PART_TYPE, STAR_TYPE, SINK_TYPE, TREE_TYPE
  use pm_parameters, only: CORE_TYPE
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer,dimension(1:input_size)::input_array
  type(out_input_star_t)::output, next_output

  integer::rID
  integer::p_type

  !--------------------------------------------------------------------
  ! This routine is the recursive slave procedure to read and dispatch
  ! particles from a Ramses output file.
  !--------------------------------------------------------------------

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_INPUT_PART_RAMSES,pst%iUpper+1,input_size,output_size,input_array)
     call r_input_part_ramses(pst%pLower,input_array,input_size,output,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_output)
     output%mass=output%mass+next_output%mass
  else
     p_type=input_array(1)
     if(p_type==PART_TYPE)then
        call input_part_ramses(pst%s%r,pst%s%g,pst%s%p,input_size-1,input_array(2:input_size),output%mass)
     endif
     if(p_type==STAR_TYPE)then
        call input_part_ramses(pst%s%r,pst%s%g,pst%s%star,input_size-1,input_array(2:input_size),output%mass)
     endif
     if(p_type==SINK_TYPE)then
        call input_part_ramses(pst%s%r,pst%s%g,pst%s%sink,input_size-1,input_array(2:input_size),output%mass)
     endif
     if(p_type==CORE_TYPE)then
        call input_part_ramses(pst%s%r,pst%s%g,pst%s%core,input_size-1,input_array(2:input_size),output%mass)
     endif
  endif

end subroutine r_input_part_ramses
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine input_part_ramses(r,g,p,ncpu_file,npart_file,mpart_loc)
  use amr_parameters, only: ndim, i8b, flen
  use amr_commons, only: run_t, global_t
  use pm_parameters, only: PART_TYPE, STAR_TYPE, SINK_TYPE
  use pm_parameters, only: CORE_TYPE
  use pm_commons, only: part_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(part_t)::p
  real(kind=8)::mpart_loc
  integer::ncpu_file
  integer,dimension(1:ncpu_file)::npart_file
  !------------------------------------------------------------
  ! Read particles positions and velocities from a Ramses
  ! output file and allocate particle-based arrays.
  !------------------------------------------------------------
  integer::ipart,ipart_old
  integer::i,idim,icpu,ileft,iright,nrest
  integer(kind=8)::ipos,iskip
  integer::istart,iend
  integer(kind=8)::nleft,nright,npart_tot
  integer(kind=8),dimension(0:ncpu_file)::ncum_file

  real(kind=4),allocatable,dimension(:)::xsp
  integer,allocatable,dimension(:)::isp
  integer(i8b),allocatable,dimension(:)::isp8

  character(LEN=5)::ncharcpu
  character(LEN=10)::prefix
  character(LEN=flen)::file_part

  !-------------------------------------
  ! Compute local particle number
  !-------------------------------------
  ncum_file(0)=0
  do icpu=1,ncpu_file
     ncum_file(icpu)=ncum_file(icpu-1)+npart_file(icpu)
  end do
  npart_tot=ncum_file(ncpu_file)
  p%npart_tot=npart_tot

  p%npart=npart_tot/g%ncpu
  nleft=(g%myid-1)*p%npart
  nright=g%myid*p%npart
  nrest=npart_tot-p%npart*g%ncpu
  if(g%myid.LE.nrest)then
     p%npart=p%npart+1
     nleft=nleft+(g%myid-1)
     nright=nright+g%myid
  else
     nleft=nleft+nrest
     nright=nright+nrest
  endif

  ! Compute interval of file to open for current process
  ileft=0
  iright=-1
  if(nright.GT.nleft)then
     do icpu=1,ncpu_file
        if(icpu>1)then
           if(ncum_file(icpu).GT.nleft.AND.ncum_file(icpu-1).LT.nright)then
              if(ileft==0)ileft=icpu
              iright=MAX(icpu,iright)
           endif
        else
           if(ncum_file(icpu).GT.nleft)then
              if(ileft==0)ileft=icpu
              iright=MAX(icpu,iright)
           endif
        endif
     end do
  endif

  ! Determine prefix based on particle type
  if(p%type==PART_TYPE)prefix='part'
  if(p%type==STAR_TYPE)prefix='star'
  if(p%type==SINK_TYPE)prefix='sink'
  if(p%type==CORE_TYPE)prefix='core'

  ! Loop over relevant files
  ipart=0
  ipart_old=0
  mpart_loc=0.0d0
  do icpu=ileft,iright
     if(icpu>1)then
        istart=MAX(nleft-ncum_file(icpu-1),0)+1
        iend=MIN(nright-ncum_file(icpu-1),npart_file(icpu))
     else
        istart=nleft+1
        iend=MIN(nright,npart_file(icpu))
     endif

     ! Open the particle file
     call title(icpu,ncharcpu)
     file_part=TRIM(r%initfile(r%levelmin))//'/'//TRIM(prefix)//'.'//TRIM(ncharcpu)
     open(unit=10,file=TRIM(file_part),access="stream",action="read",form='unformatted')

     iskip=9

     allocate(xsp(istart:iend))

     ! Read positions
     do idim=1,ndim
        ipos=iskip+4_8*(istart-1)
        read(10,POS=ipos)xsp
        ipart=ipart_old
        do i=istart,iend
           ipart=ipart+1
           p%xp(ipart,idim)=xsp(i)
        end do
        iskip=iskip+4_8*npart_file(icpu)
     end do

     ! Read velocities
     do idim=1,ndim
        ipos=iskip+4_8*(istart-1)
        read(10,POS=ipos)xsp
        ipart=ipart_old
        do i=istart,iend
           ipart=ipart+1
           p%vp(ipart,idim)=xsp(i)
        end do
        iskip=iskip+4_8*npart_file(icpu)
     end do

     ! Read masses
     ipos=iskip+4_8*(istart-1)
     read(10,POS=ipos)xsp
     ipart=ipart_old
     do i=istart,iend
        ipart=ipart+1
        p%mp(ipart)=xsp(i)
        mpart_loc=mpart_loc+xsp(i)
     end do
     iskip=iskip+4_8*npart_file(icpu)

     ! Skip potential (diagnostic only, not read back)
     if(allocated(p%phip))then
        iskip=iskip+4_8*npart_file(icpu)
     endif

     ! Read metallicity
     if(allocated(p%zp))then
        ipos=iskip+4_8*(istart-1)
        read(10,POS=ipos)xsp
        ipart=ipart_old
        do i=istart,iend
           ipart=ipart+1
           p%zp(ipart)=xsp(i)
        end do
        iskip=iskip+4_8*npart_file(icpu)
     endif

     ! Read accelerations
     if(allocated(p%fp))then
        do idim=1,ndim
           ipos=iskip+4_8*(istart-1)
           read(10,POS=ipos)xsp
           ipart=ipart_old
           do i=istart,iend
              ipart=ipart+1
              p%fp(ipart,idim)=xsp(i)
           end do
           iskip=iskip+4_8*npart_file(icpu)
        end do
     endif

     ! Read Angular momenta
     if(allocated(p%jp))then
        do idim=1,ndim
           ipos=iskip+4_8*(istart-1)
           read(10,POS=ipos)xsp
           ipart=ipart_old
           do i=istart,iend
              ipart=ipart+1
              p%jp(ipart,idim)=xsp(i)
           end do
           iskip=iskip+4_8*npart_file(icpu)
        end do
     endif

     ! Read birth time
     if(allocated(p%tp))then
        ipos=iskip+4_8*(istart-1)
        read(10,POS=ipos)xsp
        ipart=ipart_old
        do i=istart,iend
           ipart=ipart+1
           p%tp(ipart)=xsp(i)
        end do
        iskip=iskip+4_8*npart_file(icpu)
     endif

     ! Read merging time
     if(allocated(p%tm))then
        ipos=iskip+4_8*(istart-1)
        read(10,POS=ipos)xsp
        ipart=ipart_old
        do i=istart,iend
           ipart=ipart+1
           p%tm(ipart)=xsp(i)
        end do
        iskip=iskip+4_8*npart_file(icpu)
     endif

     ! Read specific internal energy
     if(allocated(p%up))then
        ipos=iskip+4_8*(istart-1)
        read(10,POS=ipos)xsp
        ipart=ipart_old
        do i=istart,iend
           ipart=ipart+1
           p%up(ipart)=xsp(i)
        end do
        iskip=iskip+4_8*npart_file(icpu)
     endif

     allocate(isp(istart:iend))

     ! Read level
     ipos=iskip+4_8*(istart-1)
     read(10,POS=ipos)isp
     ipart=ipart_old
     do i=istart,iend
        ipart=ipart+1
        p%levelp(ipart)=isp(i)
     end do
     iskip=iskip+4_8*npart_file(icpu)

     ! Read last accretion step
     if(allocated(p%nacc))then
        ipos=iskip+4_8*(istart-1)
        read(10,POS=ipos)isp
        ipart=ipart_old
        do i=istart,iend
           ipart=ipart+1
           p%nacc(ipart)=isp(i)
        end do
        iskip=iskip+4_8*npart_file(icpu)
     endif

     deallocate(isp)

     allocate(isp8(istart:iend))

     ! Read identity
#ifndef LONGINT
     ipos=iskip+4_8*(istart-1)
#else
     ipos=iskip+8_8*(istart-1)
#endif
     read(10,POS=ipos)isp8
     ipart=ipart_old
     do i=istart,iend
        ipart=ipart+1
        p%idp(ipart)=isp8(i)
     end do
#ifndef LONGINT
     iskip=iskip+4_8*npart_file(icpu)
#else
     iskip=iskip+8_8*npart_file(icpu)
#endif

     ! Read merging identity
     if(allocated(p%idm))then
#ifndef LONGINT
        ipos=iskip+4_8*(istart-1)
#else
        ipos=iskip+8_8*(istart-1)
#endif
        read(10,POS=ipos)isp8
        ipart=ipart_old
        do i=istart,iend
           ipart=ipart+1
           p%idm(ipart)=isp8(i)
        end do
#ifndef LONGINT
        iskip=iskip+4_8*npart_file(icpu)
#else
        iskip=iskip+8_8*npart_file(icpu)
#endif
     endif

     ! Read tracking identity
     if(allocated(p%idt))then
#ifndef LONGINT
        ipos=iskip+4_8*(istart-1)
#else
        ipos=iskip+8_8*(istart-1)
#endif
        read(10,POS=ipos)isp8
        ipart=ipart_old
        do i=istart,iend
           ipart=ipart+1
           p%idt(ipart)=isp8(i)
        end do
#ifndef LONGINT
        iskip=iskip+4_8*npart_file(icpu)
#else
        iskip=iskip+8_8*npart_file(icpu)
#endif
     endif

     deallocate(isp8)

     ! Read size
     if(allocated(p%size))then
        ipos=iskip+4_8*(istart-1)
        read(10,POS=ipos)xsp
        ipart=ipart_old
        do i=istart,iend
           ipart=ipart+1
           p%size(ipart)=xsp(i)
        end do
        iskip=iskip+4_8*npart_file(icpu)
     endif

     ! Read charge
     if(allocated(p%charge))then
        ipos=iskip+4_8*(istart-1)
        read(10,POS=ipos)xsp
        ipart=ipart_old
        do i=istart,iend
           ipart=ipart+1
           p%charge(ipart)=xsp(i)
        end do
        iskip=iskip+4_8*npart_file(icpu)
     endif

     deallocate(xsp)

     ! Close the particle file
     close(10)
     ipart_old=ipart

  end do
  ! End loop over files

  ! Periodic box
  do ipart=1,p%npart
     do idim=1,ndim
        if(r%periodic(idim))then
           if(p%xp(ipart,idim)< 0.0d0           )p%xp(ipart,idim)=p%xp(ipart,idim)+r%box_size(idim)
           if(p%xp(ipart,idim)>=r%box_size(idim))p%xp(ipart,idim)=p%xp(ipart,idim)-r%box_size(idim)
        endif
     end do
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

end subroutine input_part_ramses
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
end module input_part_ramses_module
