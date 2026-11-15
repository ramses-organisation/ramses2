module ramses_commons
  use amr_commons, only: run_t,global_t,mesh_t,hydro_params_t
  use pm_commons, only: part_t
  use turb_commons, only: turb_t
  use mdl_module, only: mdl_t
  use clfind_commons, only: clump_t
  use cooling_module, only: cooling_t
  use coolrates_module, only: neq_cooling_t
  use SED_module, only: sed_table_t

  type ramses_t

     type(run_t)::r
     type(global_t)::g
     type(mesh_t),pointer::m => null()
     type(mesh_t),pointer::m_mg => null()
     type(part_t)::p
     type(part_t)::star
     type(part_t)::sink
     type(part_t)::core
     type(part_t)::tree
     type(part_t)::trac
     type(part_t)::dust
     type(part_t)::gas
     type(clump_t),pointer::c => null()
     type(turb_t)::turb
     type(cooling_t)::cool
     type(neq_cooling_t)::tables
     type(sed_table_t)::SED
     type(hydro_params_t)::h_params
     type(mdl_t),pointer::mdl => null()

  end type ramses_t

  type pst_t

     type(ramses_t),pointer::s => null()
     type(pst_t),pointer::pLower => null()
     integer::iUpper = -1
     integer::nLower = 0
     integer::nUpper = 0

  end type pst_t

contains
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine open_file(s,filename,nskip,ilun)
  use hydro_parameters, only: nvar
  use rt_parameters, only: nrtvar
  use cr_parameters, only: ncruvar
  use amr_parameters, only: ndim,flen,twotondim
#ifndef WITHOUTMPI
  use mpi
#endif
  type(ramses_t)::s
  character(len=flen),intent(in)::filename
  integer(kind=8),dimension(s%r%levelmin:s%r%nlevelmax),intent(out)::nskip
  integer,intent(out)::ilun
#ifndef WITHOUTMPI
  integer::info,reqsend,reqrecv,tag1=100,tag2=200
  integer,dimension(MPI_STATUS_SIZE)::reqstatus
#endif
  character(LEN=flen)::fileloc
  character(LEN=5)::nchar
  integer,dimension(1:s%r%nfile+1)::istart
  integer::i,ifile,ncpufile,nremain,ilevel,ierr
  integer(kind=8)::iskip=0
  integer(kind=8),dimension(s%r%levelmin:s%r%nlevelmax)::noct
  logical::file_exist

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

  ncpufile=g%ncpu/r%nfile
  nremain=g%ncpu-r%nfile*ncpufile
  istart(1)=1
  do i=1,r%nfile
     if(i.LE.nremain)then
        istart(i+1)=istart(i)+ncpufile+1
     else
        istart(i+1)=istart(i)+ncpufile
     endif
  end do

  ifile=1
  do i=1,r%nfile
     if(g%myid.GE.istart(i).AND.g%myid.LT.istart(i+1))then
        ifile=i
        exit
     end if
  end do

#ifndef WITHOUTMPI
  if(g%myid.EQ.istart(ifile+1)-1)then
     noct=m%noct(s%r%levelmin:s%r%nlevelmax)
     if(g%myid.GT.istart(ifile))then
     call MPI_ISEND(noct,r%nlevelmax-r%levelmin+1,MPI_INTEGER8,g%myid-2,tag1,MPI_COMM_WORLD,reqsend,info)
     call MPI_WAIT(reqsend,reqstatus,info)
     endif
  elseif(g%myid.GT.istart(ifile))then
     call MPI_IRECV(noct,r%nlevelmax-r%levelmin+1,MPI_INTEGER8,g%myid,tag1,MPI_COMM_WORLD,reqrecv,info)
     call MPI_WAIT(reqrecv,reqstatus,info)
     noct=noct+m%noct(s%r%levelmin:s%r%nlevelmax)
     call MPI_ISEND(noct,r%nlevelmax-r%levelmin+1,MPI_INTEGER8,g%myid-2,tag1,MPI_COMM_WORLD,reqsend,info)
     call MPI_WAIT(reqsend,reqstatus,info)
  elseif(g%myid.EQ.istart(ifile))then
     call MPI_IRECV(noct,r%nlevelmax-r%levelmin+1,MPI_INTEGER8,g%myid,tag1,MPI_COMM_WORLD,reqrecv,info)
     call MPI_WAIT(reqrecv,reqstatus,info)
     noct=noct+m%noct(s%r%levelmin:s%r%nlevelmax)
  else
     write(*,*)'Problem in open_file'
  endif
#else
  noct=m%noct(s%r%levelmin:s%r%nlevelmax)
#endif

  if(g%myid.EQ.istart(ifile))then

     ilun=10+istart(ifile)
     call title(ifile,nchar)
     fileloc=TRIM(filename)//TRIM(nchar)
     inquire(file=fileloc, exist=file_exist)
     if (file_exist) then
        open(unit=ilun,file=fileloc,iostat=ierr)
        close(ilun,status="delete")
     end if

     if(index(filename,'clump').NE.0)then
        open(unit=ilun,file=fileloc,form='formatted')
        write(ilun,'(240A)')'     index       parent     halo      ncell      npart'//&
             & '    rho_min            rho_max            rho_sad            rho_ave    '//&
             & '        mpatch             mass               r200               rmax               c200   '//&
             & '            pos_x              pos_y              pos_z      '//&
             & '        vel_x              vel_y              vel_z'
     else
        open(unit=ilun,file=fileloc,access="stream",action="write",form='unformatted')
        write(ilun)ndim
#ifdef MHD
        if(index(filename,'hydro').NE.0)write(ilun)nvar+3
#else
        if(index(filename,'hydro').NE.0)write(ilun)nvar
#endif
        if(index(filename,'grav').NE.0)write(ilun)ndim+1
        if(index(filename,'peak').NE.0)write(ilun)3
        if(index(filename,'rt').NE.0)write(ilun)nrtvar
        if(index(filename,'cr').NE.0)write(ilun)ncruvar
        write(ilun)r%levelmin
        write(ilun)r%nlevelmax
        do ilevel=r%levelmin,r%nlevelmax
           if(noct(ilevel)>2147483647_8)then
              write(*,*)'Too many octs at level=',ilevel
              write(*,*)'Increase nfile=',r%nfile
              stop
           endif
           write(ilun)int(noct(ilevel),kind=4)
        end do
     endif

     if(index(filename,'amr').NE.0)iskip=13+(r%nlevelmax-r%levelmin+1)*4
     if(index(filename,'hydro').NE.0)iskip=17+(r%nlevelmax-r%levelmin+1)*4
     if(index(filename,'grav').NE.0)iskip=17+(r%nlevelmax-r%levelmin+1)*4
     if(index(filename,'peak').NE.0)iskip=17+(r%nlevelmax-r%levelmin+1)*4
     if(index(filename,'rt').NE.0)iskip=17+(r%nlevelmax-r%levelmin+1)*4
     if(index(filename,'cr').NE.0)iskip=17+(r%nlevelmax-r%levelmin+1)*4

     do ilevel=r%levelmin,r%nlevelmax
        nskip(ilevel)=iskip
        if(index(filename,'amr').NE.0)iskip=iskip+(4*ndim+4)*noct(ilevel)
#ifdef MHD
        if(index(filename,'hydro').NE.0)iskip=iskip+(4*twotondim*(nvar+3))*noct(ilevel)
#else
        if(index(filename,'hydro').NE.0)iskip=iskip+(4*twotondim*nvar)*noct(ilevel)
#endif
        if(index(filename,'grav').NE.0)iskip=iskip+(4*twotondim*(ndim+1))*noct(ilevel)
        if(index(filename,'peak').NE.0)iskip=iskip+(4*twotondim*3)*noct(ilevel)
        if(index(filename,'rt').NE.0)iskip=iskip+(4*twotondim*nrtvar)*noct(ilevel)
        if(index(filename,'cr').NE.0)iskip=iskip+(4*twotondim*ncruvar)*noct(ilevel)
     end do

  elseif(g%myid.GT.istart(ifile))then

#ifndef WITHOUTMPI
     call MPI_IRECV(nskip,r%nlevelmax-r%levelmin+1,MPI_INTEGER8,g%myid-2,tag2,MPI_COMM_WORLD,reqrecv,info)
     call MPI_WAIT(reqrecv,reqstatus,info)
#endif
     ilun=10+istart(ifile)
     call title(ifile,nchar)
     fileloc=TRIM(filename)//TRIM(nchar)
     inquire(file=fileloc, exist=file_exist)
     if (.not. file_exist) then
        write(*,*)'Problem in open_file'
        stop
     end if

     if(index(filename,'clump').NE.0)then
        open(unit=ilun,file=fileloc,form='formatted',access='append')
     else
        open(unit=ilun,file=fileloc,access="stream",action="write",form='unformatted')
     endif

  end if

  end associate

end subroutine open_file
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine close_file(s,filename,nskip,ilun)
  use hydro_parameters, only: nvar
  use rt_parameters, only: nrtvar
  use cr_parameters, only: ncruvar
  use amr_parameters, only: ndim,flen,twotondim
#ifndef WITHOUTMPI
  use mpi
#endif
  type(ramses_t)::s
  character(len=flen),intent(in)::filename
  integer(kind=8),dimension(s%r%levelmin:s%r%nlevelmax),intent(inout)::nskip
  integer,intent(in)::ilun
#ifndef WITHOUTMPI
  integer::info,reqsend,tag2=200
  integer,dimension(MPI_STATUS_SIZE)::reqstatus
#endif
  integer,dimension(1:s%r%nfile+1)::istart
  integer::ncpufile,nremain
  integer::i,ifile,ilevel
  integer(kind=8)::iskip

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

  ncpufile=g%ncpu/r%nfile
  nremain=g%ncpu-r%nfile*ncpufile
  istart(1)=1
  do i=1,r%nfile
     if(i.LE.nremain)then
        istart(i+1)=istart(i)+ncpufile+1
     else
        istart(i+1)=istart(i)+ncpufile
     endif
  end do

  ifile=1
  do i=1,r%nfile
     if(g%myid.GE.istart(i).AND.g%myid.LT.istart(i+1))then
        ifile=i
        exit
     end if
  end do

  close(ilun)

  if(g%myid.LT.istart(ifile+1)-1)then

     do ilevel=r%levelmin,r%nlevelmax
        iskip=nskip(ilevel)
        if(index(filename,'amr').NE.0)iskip=iskip+(4*ndim+4)*int(m%noct(ilevel),kind=8)
#ifdef MHD
        if(index(filename,'hydro').NE.0)iskip=iskip+(4*twotondim*(nvar+3))*int(m%noct(ilevel),kind=8)
#else
        if(index(filename,'hydro').NE.0)iskip=iskip+(4*twotondim*nvar)*int(m%noct(ilevel),kind=8)
#endif
        if(index(filename,'grav').NE.0)iskip=iskip+(4*twotondim*(ndim+1))*int(m%noct(ilevel),kind=8)
        if(index(filename,'peak').NE.0)iskip=iskip+(4*twotondim*3)*int(m%noct(ilevel),kind=8)
        if(index(filename,'rt').NE.0)iskip=iskip+(4*twotondim*nrtvar)*int(m%noct(ilevel),kind=8)
        if(index(filename,'cr').NE.0)iskip=iskip+(4*twotondim*ncruvar)*int(m%noct(ilevel),kind=8)
        nskip(ilevel)=iskip
     end do

#ifndef WITHOUTMPI
     call MPI_ISEND(nskip,r%nlevelmax-r%levelmin+1,MPI_INTEGER8,g%myid,tag2,MPI_COMM_WORLD,reqsend,info)
     call MPI_WAIT(reqsend,reqstatus,info)
#endif
  end if

  end associate

end subroutine close_file
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine open_part_file(s,p,filename,nskip,ilun)
  use amr_parameters, only: ndim,flen,twotondim,i8b
  use pm_commons, only: part_t
#ifndef WITHOUTMPI
  use mpi
#endif
  type(ramses_t)::s
  type(part_t)::p
  character(len=flen),intent(in)::filename
  integer(kind=8), dimension(1:p%nvaralloc+1),intent(out)::nskip
  integer, intent(out)::ilun
#ifndef WITHOUTMPI
  integer::info,reqsend,reqrecv,tag1=100,tag2=200
  integer,dimension(MPI_STATUS_SIZE)::reqstatus
#endif
  character(LEN=flen)::fileloc
  character(LEN=5)::nchar
  integer,dimension(1:s%r%nfile+1)::istart
  integer::i,idim,ivar,ifile,ncpufile,nremain,ierr
  integer(kind=8)::npart
  logical::file_exist

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

  ncpufile=g%ncpu/r%nfile
  nremain=g%ncpu-r%nfile*ncpufile
  istart(1)=1
  do i=1,r%nfile
     if(i.LE.nremain)then
        istart(i+1)=istart(i)+ncpufile+1
     else
        istart(i+1)=istart(i)+ncpufile
     endif
  end do

  ifile=1
  do i=1,r%nfile
     if(g%myid.GE.istart(i).AND.g%myid.LT.istart(i+1))then
        ifile=i
        exit
     end if
  end do

#ifndef WITHOUTMPI
  if(g%myid.EQ.istart(ifile+1)-1)then
     npart=p%npart
     if(g%myid.GT.istart(ifile))then
     call MPI_ISEND(npart,1,MPI_INTEGER8,g%myid-2,tag1,MPI_COMM_WORLD,reqsend,info)
     call MPI_WAIT(reqsend,reqstatus,info)
     endif
  elseif(g%myid.GT.istart(ifile))then
     call MPI_IRECV(npart,1,MPI_INTEGER8,g%myid,tag1,MPI_COMM_WORLD,reqrecv,info)
     call MPI_WAIT(reqrecv,reqstatus,info)
     npart=npart+p%npart
     call MPI_ISEND(npart,1,MPI_INTEGER8,g%myid-2,tag1,MPI_COMM_WORLD,reqsend,info)
     call MPI_WAIT(reqsend,reqstatus,info)
  elseif(g%myid.EQ.istart(ifile))then
     call MPI_IRECV(npart,1,MPI_INTEGER8,g%myid,tag1,MPI_COMM_WORLD,reqrecv,info)
     call MPI_WAIT(reqrecv,reqstatus,info)
     npart=npart+p%npart
  else
     write(*,*)'Problem in open_part_file'
  endif
#else
  npart=p%npart
#endif

  if(g%myid.EQ.istart(ifile))then

     ilun=10+istart(ifile)
     call title(ifile,nchar)
     fileloc=TRIM(filename)//TRIM(nchar)
     inquire(file=fileloc, exist=file_exist)
     if (file_exist) then
        open(unit=ilun,file=fileloc,iostat=ierr)
        close(ilun,status="delete")
     end if
     open(unit=ilun,file=fileloc,access="stream",action="write",form='unformatted')

     write(ilun)ndim
     if(npart>2147483647_8)then
        write(*,*)'Too many particles in output_part'
        write(*,*)'Increase nfile=',r%nfile
        stop
     endif
     write(ilun)int(npart,kind=4)
     nskip(1)=9

     ! Positions
     do ivar=1,ndim
        nskip(ivar+1)=nskip(ivar)+4*npart
     end do
     ! Velocities
     do ivar=ndim+1,2*ndim
        nskip(ivar+1)=nskip(ivar)+4*npart
     end do
     ! Masses
     ivar=2*ndim+1
     nskip(ivar+1)=nskip(ivar)+4*npart
     ! Potentials
     if(allocated(p%phip))then
        ivar=ivar+1
        nskip(ivar+1)=nskip(ivar)+4*npart
     endif
     ! Metallicities
     if(allocated(p%zp))then
        ivar=ivar+1
        nskip(ivar+1)=nskip(ivar)+4*npart
     endif
     ! Acceleration
     if(allocated(p%fp))then
        do idim=1,ndim
           ivar=ivar+1
           nskip(ivar+1)=nskip(ivar)+4*npart
        end do
     endif
     ! Angular momentum
     if(allocated(p%jp))then
        do idim=1,ndim
           ivar=ivar+1
           nskip(ivar+1)=nskip(ivar)+4*npart
        end do
     endif
     ! Birth times
     if(allocated(p%tp))then
        ivar=ivar+1
        nskip(ivar+1)=nskip(ivar)+4*npart
     endif
     ! Merging times
     if(allocated(p%tm))then
        ivar=ivar+1
        nskip(ivar+1)=nskip(ivar)+4*npart
     endif
     ! Specific internal energy
     if(allocated(p%up))then
        ivar=ivar+1
        nskip(ivar+1)=nskip(ivar)+4*npart
     endif
     ! Levels
     ivar=ivar+1
     nskip(ivar+1)=nskip(ivar)+4*npart
     ! Last accretion steps
     if(allocated(p%nacc))then
        ivar=ivar+1
        nskip(ivar+1)=nskip(ivar)+4*npart
     endif
     ! Identities
     ivar=ivar+1
     nskip(ivar+1)=nskip(ivar)+i8b*npart
     ! Merging identities
     if(allocated(p%idm))then
        ivar=ivar+1
        nskip(ivar+1)=nskip(ivar)+i8b*npart
     endif
     ! Tracking identities
     if(allocated(p%idt))then
        ivar=ivar+1
        nskip(ivar+1)=nskip(ivar)+i8b*npart
     endif
     if(allocated(p%size))then
        ivar=ivar+1
        nskip(ivar+1)=nskip(ivar)+4*npart
     endif
     if(allocated(p%charge))then
        ivar=ivar+1
        nskip(ivar+1)=nskip(ivar)+4*npart
     endif

  elseif(g%myid.GT.istart(ifile))then

#ifndef WITHOUTMPI
     call MPI_IRECV(nskip,p%nvaralloc+1,MPI_INTEGER8,g%myid-2,tag2,MPI_COMM_WORLD,reqrecv,info)
     call MPI_WAIT(reqrecv,reqstatus,info)
#endif
     ilun=10+istart(ifile)
     call title(ifile,nchar)
     fileloc=TRIM(filename)//TRIM(nchar)
     inquire(file=fileloc, exist=file_exist)
     if (.not. file_exist) then
        write(*,*)'Problem in open_part_file'
        stop
     end if
     open(unit=ilun,file=fileloc,access="stream",action="write",form='unformatted')

  end if

  end associate

end subroutine open_part_file
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine close_part_file(s,p,filename,nskip,ilun)
  use amr_parameters, only: ndim,flen,twotondim,i8b
  use pm_commons, only: part_t
#ifndef WITHOUTMPI
  use mpi
#endif
  type(ramses_t)::s
  type(part_t)::p
  character(len=flen),intent(in)::filename
  integer(kind=8), dimension(1:p%nvaralloc+1),intent(inout)::nskip
  integer, intent(in)::ilun
#ifndef WITHOUTMPI
  integer::info,reqsend,tag2=200
  integer,dimension(MPI_STATUS_SIZE)::reqstatus
#endif
  integer,dimension(1:s%r%nfile+1)::istart
  integer::ncpufile,nremain
  integer::i,ivar,idim,ifile

  associate(r=>s%r,g=>s%g)

  ncpufile=g%ncpu/r%nfile
  nremain=g%ncpu-r%nfile*ncpufile
  istart(1)=1
  do i=1,r%nfile
     if(i.LE.nremain)then
        istart(i+1)=istart(i)+ncpufile+1
     else
        istart(i+1)=istart(i)+ncpufile
     endif
  end do

  ifile=1
  do i=1,r%nfile
     if(g%myid.GE.istart(i).AND.g%myid.LT.istart(i+1))then
        ifile=i
        exit
     end if
  end do

  close(ilun)

  if(g%myid.LT.istart(ifile+1)-1)then

     ! Positions
     do ivar=1,ndim
        nskip(ivar)=nskip(ivar)+4*p%npart
     end do
     ! Velocities
     do ivar=ndim+1,2*ndim
        nskip(ivar)=nskip(ivar)+4*p%npart
     end do
     ! Masses
     ivar=2*ndim+1
     nskip(ivar)=nskip(ivar)+4*p%npart
     ! Potentials
     if(allocated(p%phip))then
        ivar=ivar+1
        nskip(ivar)=nskip(ivar)+4*p%npart
     endif
     ! Metallicities
     if(allocated(p%zp))then
        ivar=ivar+1
        nskip(ivar)=nskip(ivar)+4*p%npart
     endif
     ! Accelerations
     if(allocated(p%fp))then
        do idim=1,ndim
           ivar=ivar+1
           nskip(ivar)=nskip(ivar)+4*p%npart
        end do
     endif
     ! Angular momenta
     if(allocated(p%jp))then
        do idim=1,ndim
           ivar=ivar+1
           nskip(ivar)=nskip(ivar)+4*p%npart
        end do
     endif
     ! Birth times
     if(allocated(p%tp))then
        ivar=ivar+1
        nskip(ivar)=nskip(ivar)+4*p%npart
     endif
     ! Merging times
     if(allocated(p%tm))then
        ivar=ivar+1
        nskip(ivar)=nskip(ivar)+4*p%npart
     endif
     ! Specific internal energy
     if(allocated(p%up))then
        ivar=ivar+1
        nskip(ivar)=nskip(ivar)+4*p%npart
     endif
     ! Levels
     ivar=ivar+1
     nskip(ivar)=nskip(ivar)+4*p%npart
     ! Last accretion steps
     if(allocated(p%nacc))then
        ivar=ivar+1
        nskip(ivar)=nskip(ivar)+4*p%npart
     endif
     ! Identities
     ivar=ivar+1
     nskip(ivar)=nskip(ivar)+i8b*p%npart
     ! Merging identities
     if(allocated(p%idm))then
        ivar=ivar+1
        nskip(ivar)=nskip(ivar)+i8b*p%npart
     endif
     ! Tracking identities
     if(allocated(p%idt))then
        ivar=ivar+1
        nskip(ivar)=nskip(ivar)+i8b*p%npart
     endif
     if(allocated(p%size))then
        ivar=ivar+1
        nskip(ivar)=nskip(ivar)+4*p%npart
     endif
     if(allocated(p%charge))then
        ivar=ivar+1
        nskip(ivar)=nskip(ivar)+4*p%npart
     endif

#ifndef WITHOUTMPI
     call MPI_ISEND(nskip,p%nvaralloc+1,MPI_INTEGER8,g%myid,tag2,MPI_COMM_WORLD,reqsend,info)
     call MPI_WAIT(reqsend,reqstatus,info)
#endif
  end if

  end associate

end subroutine close_part_file
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
end module ramses_commons
!##############################################################
!##############################################################
!##############################################################
!##############################################################
module clean_stop_module
contains
recursive subroutine r_clean_stop(pst)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_CLEAN_STOP,pst%iUpper+1)
     call r_clean_stop(pst%pLower)
     call mdl_get_reply(pst%s%mdl,rID,0)
  endif

end subroutine r_clean_stop
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine m_clean_abort(pst)
  ! Stop all ranks cleanly from the master, while the workers wait in mdl_wait
  use ramses_commons, only: pst_t
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  type(pst_t)::pst
#ifndef WITHOUTMPI
  integer::info
#endif

  call r_clean_stop(pst)
#ifndef WITHOUTMPI
  call MPI_FINALIZE(info)
#endif
  stop

end subroutine m_clean_abort
end module clean_stop_module
