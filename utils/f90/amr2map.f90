program amr2map
  !--------------------------------------------------------------------------
  ! This code projects RAMSES data onto a map.
  ! This is based in the MINI-RAMSES prototype.
  ! R. Teyssier, Princeton, February 2nd 2023
  !--------------------------------------------------------------------------
  implicit none

  integer,parameter::flen=200

  integer::ndim,twotondim,nvar
  integer::n,i,j,k,type=0,dopeak=0,domax=0,dograv=0,dort=0,backup=0
  integer::ivar,lmax=0
  integer::ilevel,idim,jdim,kdim,ldim,icell
  integer::nlevel
  integer::ind,icpu,ncpu_read

  integer::nx_sample=0,ny_sample=0
  integer::imin,imax,jmin,jmax,kmin,kmax
  integer::nstride,nx_full,ny_full,lmin
  integer::ix,iy,iz,ndom,impi,bit_length,maxdom
  integer::noct_file,noct_tmp
  integer,dimension(1:8)::idom,jdom,kdom,cpu_min,cpu_max
  integer(kind=8)::ipos,iskip_amr,iskip_hydro,noct_skip

  real(KIND=8)::dxline,weight
  real(KIND=8)::xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1
  real(KIND=8)::xxmin,xxmax,yymin,yymax,zzmin,zzmax
  real(KIND=8)::dx,dy,dz,xx,yy,zz
  real(KIND=8)::rho,map,ekin,bx,by,bz,emag,pres
  real(KIND=8)::metmax=0d0

  character(LEN=1)::proj='z'
  character(LEN=5)::nchar,ncharcpu
  character(LEN=flen)::nomfich,repository,outfich,mapfiletype='bin'
  character(LEN=flen)::file_amr,file_hydro
  logical::ok,ok_part,ok_cell,do_max,do_grav,do_peak,do_rt
  logical::backup_file=.false.
  logical::check_ramses_exist

  real(KIND=4),dimension(:,:),allocatable::filemap
  integer,dimension(:),allocatable::cpu_list

  integer,dimension(:),allocatable::ckey,cart_key
  integer::refined_int
  logical,dimension(:),allocatable::refined
  real(KIND=8),dimension(:,:),allocatable::uold
  real(KIND=4),dimension(:,:),allocatable::qold

  type level
     integer::ilevel
     real(KIND=8),dimension(:,:),pointer::map
     real(KIND=8),dimension(:,:),pointer::rho
     integer::imin
     integer::imax
     integer::jmin
     integer::jmax
  end type level

  type params
     integer::ndim
     integer::ncpu
     integer::nvar
     integer::levelmin
     integer::nlevelmax
     real(kind=8)::boxlen
     real(kind=8)::t
     real(kind=8)::gamma
     integer::nhilbert
     integer(kind=8),dimension(:,:,:),allocatable::bound_key
  end type params

  type(level),dimension(1:100)::mapgrid
  type(params)::p

  write(*,*)'Starting amr2map'

  ! Read amr2map parameters
  call read_params

  ! Check that all files exist
  if(.NOT. check_ramses_exist(repository))then
     write(*,*)'Repository '//TRIM(repository)//' incomplete.'
     write(*,*)'Stopping.'
     stop
  endif
  if(index(repository,'output')==0)backup_file=.true.

  ! Read RAMSES params
  call read_ramses_params
  write(*,'(A,1PE12.5)')' time= ',p%t
  ndim=p%ndim
  twotondim=2**ndim
  nvar=p%nvar

  !------------------------------
  ! Set up geometry parameters
  !------------------------------
  if(lmax==0)then
     lmax=p%nlevelmax
  endif
  write(*,'(A,I0)')' Working resolution = ',2**lmax
  zzmax=1.0
  zzmin=0.0
  if(p%ndim>2)then
     if (proj=='x')then
        idim=2
        jdim=3
        kdim=1
        xxmin=ymin ; xxmax=ymax
        yymin=zmin ; yymax=zmax
        zzmin=xmin ; zzmax=xmax
     else if (proj=='y') then
        idim=1
        jdim=3
        kdim=2
        xxmin=xmin ; xxmax=xmax
        yymin=zmin ; yymax=zmax
        zzmin=ymin ; zzmax=ymax
     else
        idim=1
        jdim=2
        kdim=3
        xxmin=xmin ; xxmax=xmax
        yymin=ymin ; yymax=ymax
        zzmin=zmin ; zzmax=zmax
     end if
  else
     idim=1
     jdim=2
     xxmin=xmin ; xxmax=xmax
     yymin=ymin ; yymax=ymax
  end if

  !-----------------------------
  ! Allocate map hierarchy
  !-----------------------------
  do ilevel=p%levelmin,lmax
     nx_full=2**ilevel
     ny_full=2**ilevel
     imin=int(xxmin*dble(nx_full))+1
     imax=int(xxmax*dble(nx_full))+1
     jmin=int(yymin*dble(ny_full))+1
     jmax=int(yymax*dble(ny_full))+1
     allocate(mapgrid(ilevel)%map(imin:imax,jmin:jmax))
     allocate(mapgrid(ilevel)%rho(imin:imax,jmin:jmax))
     mapgrid(ilevel)%map(:,:)=0.0
     mapgrid(ilevel)%rho(:,:)=0.0
     mapgrid(ilevel)%imin=imin
     mapgrid(ilevel)%imax=imax
     mapgrid(ilevel)%jmin=jmin
     mapgrid(ilevel)%jmax=jmax
  end do

  !-----------------------------------------------
  ! Compute projected variables
  !----------------------------------------------
  allocate(cpu_list(1:p%ncpu))
  allocate(ckey(1:ndim),cart_key(1:ndim))
  allocate(refined(1:twotondim))
  if(backup_file)then
     allocate(uold(1:twotondim,1:nvar))
  else
     allocate(qold(1:twotondim,1:nvar))
  endif

  ! Loop over levels
  do ilevel=p%levelmin,lmax

     dx=0.5**ilevel
     zmax=1.0/dx
     dxline=1
     if(ndim==3)dxline=dx

!     write(*,*)'ilevel=',ilevel,dx,zmax

     ! Get cpu list within bounding box
!     call cmp_cpu_list(xmin,xmax,ymin,ymax,zmin,zmax,ilevel,ncpu_read,cpu_list)
     ncpu_read=p%ncpu
     cpu_list(1)=1

     ! Loop over processor files
     do k=1,ncpu_read
        icpu=k !cpu_list(k)
        call title(icpu,ncharcpu)

        ! Prepare reading the AMR file
        file_amr=TRIM(repository)//'/amr.'//TRIM(ncharcpu)

        noct_skip=0
        open(unit=10,file=file_amr,access="stream",action="read",form='unformatted')
        do i=p%levelmin,ilevel-1
           ipos=13+4*(i-p%levelmin)
           read(10,POS=ipos)noct_tmp
           noct_skip=noct_skip+noct_tmp
        end do
        ipos=13+4*(ilevel-p%levelmin)
        read(10,POS=ipos)noct_file
        iskip_amr=13+4*(p%nlevelmax-p%levelmin+1)+(4*ndim+4)*noct_skip

        ! Prepare reading the HYDRO file
        if(do_grav)then
           file_hydro=TRIM(repository)//'/grav.'//TRIM(ncharcpu)
        else if(do_peak)then
           file_hydro=TRIM(repository)//'/peak_grid.'//TRIM(ncharcpu)
        else if(do_rt)then
           file_hydro=TRIM(repository)//'/rt.'//TRIM(ncharcpu)
        else
           file_hydro=TRIM(repository)//'/hydro.'//TRIM(ncharcpu)
        endif
        open(unit=11,file=file_hydro,access="stream",action="read",form='unformatted')
        if(backup_file)then
           iskip_hydro=17+4*(p%nlevelmax-p%levelmin+1)+(8*twotondim*nvar)*noct_skip
        else
           iskip_hydro=17+4*(p%nlevelmax-p%levelmin+1)+(4*twotondim*nvar)*noct_skip
        endif
        ! Loop over useful octs in file
        do i=1,noct_file

           ! Read values from files
           ipos=iskip_amr+(4*ndim+4)*(i-1)
           read(10,POS=ipos)ckey
           ipos=iskip_amr+(4*ndim+4)*(i-1)+4*ndim
           read(10,POS=ipos)refined_int
           do ind=1,twotondim
              refined(ind)=btest(refined_int,ind-1)
           end do
           if(backup_file)then
              ipos=iskip_hydro+(8*twotondim*nvar)*(i-1)
              read(11,POS=ipos)uold
           else
              ipos=iskip_hydro+(4*twotondim*nvar)*int((i-1),kind=8)
              read(11,POS=ipos)qold
           endif
           ! Loop over 2**ndim cells
           do ind=1,twotondim

              ! Compute Cartesian key
              do ldim=1,ndim
                 nstride=2**(ldim-1)
                 cart_key(ldim)=2*ckey(ldim)+MOD((ind-1)/nstride,2)
              end do

              ok_cell=.not.refined(ind).OR.ilevel.eq.lmax

              if(backup_file)then
                 ! Get variable to map out
                 rho = uold(ind,1)
                 select case (type)
                 case (-1) ! Cpu map
                    map = icpu
                 case (0) ! Refinement map
                    map = ilevel
                 case (15) ! Temperature
                    ekin = 0.5*uold(ind,2)**2/rho
                    ekin = ekin+0.5*uold(ind,3)**2/rho
                    ekin = ekin+0.5*uold(ind,4)**2/rho
                    pres = (p%gamma-1)*(uold(ind,5)-ekin)
                    if(do_max)then
                       map = pres/rho
                    else
                       map = pres
                    endif
                 case (16) ! Magnetic pressure
                    bx = 0.5*(uold(ind,nvar-5)+uold(ind,nvar-2))
                    by = 0.5*(uold(ind,nvar-4)+uold(ind,nvar-1))
                    bz = 0.5*(uold(ind,nvar-3)+uold(ind,nvar  ))
                    emag = 0.5*(bx**2+by**2+bz**2)
                    if(do_max)then
                       map = emag
                    else
                       map = rho*emag
                    endif
                 case default ! Passive scalar
                    if(do_max)then
                       map = uold(ind,type)
                    else
                       map = rho*uold(ind,type)
                    endif
                    metmax=max(metmax,uold(ind,type))
                 end select
              else
                 ! Get variable to map out
                 rho = qold(ind,1)
                 select case (type)
                 case (-1) ! Cpu map
                    map = icpu
                 case (0) ! Refinement map
                    map = ilevel
                    metmax = max(metmax,map)
                 case (15) ! Temperature
                    pres = qold(ind,5)
                    metmax = max(metmax,pres/rho)
                    if(do_max)then
                       map = pres/rho
                    else
                       map = pres
                    endif
                 case (16) ! Magnetic pressure
                    pres = qold(ind,6)**2+qold(ind,7)**2+qold(ind,8)**2
                    metmax = max(metmax,pres)
                    if(do_max)then
                       map = pres
                    else
                       map = rho*pres
                    endif
                 case (66) ! H2 fraction
                    metmax = max(metmax,real(qold(ind,type),kind=8))
                    if(do_max)then
                       map = qold(ind,type)
                    else
                       map = rho*(1.-qold(ind,6)-qold(ind,7))
                    endif
                 case default ! Passive scalar
                    metmax = max(metmax,real(qold(ind,type),kind=8))
                    if(do_max)then
                       map = qold(ind,type)
                    else
                       map = rho*qold(ind,type)
                    endif
                 end select
              endif

              if(ok_cell)then
                 ix=cart_key(idim)+1
                 iy=cart_key(jdim)+1
                 if(ndim==3)then
                    zz=(cart_key(kdim)+0.5)*dx
                 endif

                 if(ndim==3)then
                    weight=(min(zz+dx/2.,zzmax)-max(zz-dx/2.,zzmin))/dx
                    weight=min(1.0d0,max(weight,0.0d0))
                 else
                    weight=1.0
                 endif

                 if(    ix>=mapgrid(ilevel)%imin.and.&
                      & iy>=mapgrid(ilevel)%jmin.and.&
                      & ix<=mapgrid(ilevel)%imax.and.&
                      & iy<=mapgrid(ilevel)%jmax)then

                    if(do_max)then
                       if(weight>0d0)then
                          mapgrid(ilevel)%map(ix,iy)=max(mapgrid(ilevel)%map(ix,iy),map)
                          mapgrid(ilevel)%rho(ix,iy)=max(mapgrid(ilevel)%rho(ix,iy),rho)
                       endif
                    else
                       mapgrid(ilevel)%map(ix,iy)=mapgrid(ilevel)%map(ix,iy)+map*dxline*weight/(zzmax-zzmin)
                       mapgrid(ilevel)%rho(ix,iy)=mapgrid(ilevel)%rho(ix,iy)+rho*dxline*weight/(zzmax-zzmin)
                    endif
                 endif
              end if

           end do
           ! End loop over cells

        end do
        ! End loop over octs

        ! Close file
        close(10)
        close(11)

     end do
     ! End loop over cpu

  end do
  ! End loop over levels

  write(*,*)'Data read and projected.'
  select case (type)
  case (0) ! Refinement map
     write(*,'(A,I0)')" Max level of refinement = ",nint(metmax)
  case (15) ! Temperature
     write(*,'(A,1PE12.5)')" Max sound speed = ",sqrt(metmax)
  case (16) ! Magnetic pressure
     write(*,'(A,1PE12.5)')" Max field strength = ",sqrt(metmax)
  case default ! Passive scalar
     write(*,'(A,I0,A,1PE12.5)')" Max variable ",type," = ",metmax
  end select

  nx_full=2**lmax
  ny_full=2**lmax
  imin=int(xxmin*dble(nx_full))+1
  imax=int(xxmax*dble(nx_full))
  jmin=int(yymin*dble(ny_full))+1
  jmax=int(yymax*dble(ny_full))

  do ix=imin,imax
     xmin=((ix-0.5)/2**lmax)
     do iy=jmin,jmax
        ymin=((iy-0.5)/2**lmax)
        do ilevel=p%levelmin,lmax-1
           ndom=2**ilevel
           i=int(xmin*ndom)+1
           j=int(ymin*ndom)+1
           if(do_max) then
              mapgrid(lmax)%map(ix,iy)=max(mapgrid(lmax)%map(ix,iy), &
                   & mapgrid(ilevel)%map(i,j))
              mapgrid(lmax)%rho(ix,iy)=max(mapgrid(lmax)%rho(ix,iy), &
                   & mapgrid(ilevel)%rho(i,j))
           else
              mapgrid(lmax)%map(ix,iy)=mapgrid(lmax)%map(ix,iy) + &
                   & mapgrid(ilevel)%map(i,j)
              mapgrid(lmax)%rho(ix,iy)=mapgrid(lmax)%rho(ix,iy) + &
                   & mapgrid(ilevel)%rho(i,j)
           endif
        end do
     end do
  end do

  if(.not. do_max)then
     write(*,'(A,1PE12.5)')' Norm rho= ',sum(mapgrid(lmax)%rho(imin:imax,jmin:jmax))/(imax-imin+1)/(jmax-jmin+1)
  endif

  ! Output file
  nomfich=TRIM(outfich)
  write(*,*)'Writing map in '//TRIM(nomfich)

  if (mapfiletype=='bin')then
     write(*,*)'as unformatted binary file'
     open(unit=20,file=nomfich,form='unformatted')
     if(nx_sample==0)then
        write(20)p%t, (xxmax-xxmin)*p%boxlen, (yymax-yymin)*p%boxlen, (zzmax-zzmin)*p%boxlen
        write(20)imax-imin+1,jmax-jmin+1
        allocate(filemap(imax-imin+1,jmax-jmin+1))
        if(do_max)then
           filemap=mapgrid(lmax)%map(imin:imax,jmin:jmax)
        else
           filemap=mapgrid(lmax)%map(imin:imax,jmin:jmax)/mapgrid(lmax)%rho(imin:imax,jmin:jmax)
        endif
        write(20)filemap
     else
        if(ny_sample==0)ny_sample=nx_sample
        write(20)p%t, xxmax-xxmin, yymax-yymin, zzmax-zzmin
        write(20)nx_sample+1,ny_sample+1
        allocate(filemap(0:nx_sample,0:ny_sample))
        do i=0,nx_sample
           ix=int(dble(i)/dble(nx_sample)*dble(imax-imin+1))+imin
           ix=min(ix,imax)
           do j=0,ny_sample
              iy=int(dble(j)/dble(ny_sample)*dble(jmax-jmin+1))+jmin
              iy=min(iy,jmax)
              if(do_max)then
                 filemap(i,j)=mapgrid(lmax)%map(ix,iy)
              else
                 filemap(i,j)=mapgrid(lmax)%map(ix,iy)/mapgrid(lmax)%rho(ix,iy)
              endif
           end do
        end do
        write(20)filemap
     endif
     close(20)
  endif
  if (mapfiletype=='ascii')then
     write(*,*)'formatted text file'
     open(unit=20,file=nomfich,form='formatted')
     if(nx_sample==0)then
        do j=jmin,jmax
           do i=imin,imax
              xx=xxmin+(dble(i-imin)+0.5)/dble(imax-imin+1)*(xxmax-xxmin)
              yy=yymin+(dble(j-jmin)+0.5)/dble(jmax-jmin+1)*(yymax-yymin)
              if(do_max)then
                 write(20,*)xx,yy,mapgrid(lmax)%map(i,j)
              else
                 write(20,*)xx,yy,mapgrid(lmax)%map(i,j)/mapgrid(lmax)%rho(i,j)
              endif
           end do
           write(20,*) " "
        end do
     else
        if(ny_sample==0)ny_sample=nx_sample
        allocate(filemap(0:nx_sample,0:ny_sample))
        do i=0,nx_sample
           ix=int(dble(i)/dble(nx_sample)*dble(imax-imin+1))+imin
           ix=min(ix,imax)
           do j=0,ny_sample
              iy=int(dble(j)/dble(ny_sample)*dble(jmax-jmin+1))+jmin
              iy=min(iy,jmax)
              if(do_max)then
                 filemap(i,j)=mapgrid(lmax)%map(ix,iy)
              else
                 filemap(i,j)=mapgrid(lmax)%map(ix,iy)/mapgrid(lmax)%rho(ix,iy)
           endif
           end do
        end do
        do j=0,ny_sample
           do i=0,nx_sample
              xx=xxmin+dble(i)/dble(nx_sample)*(xxmax-xxmin)
              yy=yymin+dble(j)/dble(ny_sample)*(yymax-yymin)
              write(20,*)xx,yy,filemap(i,j)
           end do
           write(20,*) " "
        end do
     endif
     close(20)
  endif

contains

  subroutine read_params

    implicit none

    integer       :: i,n
    integer       :: iargc
    character(len=4)   :: opt
    character(len=flen) :: arg
    LOGICAL       :: bad, ok
    real(kind=8)  :: xcen=0.5,ycen=0.5,zcen=0.5,rad=0

    n = iargc()
    if (n < 4) then
       print *, 'usage: amr2map   -inp  input_dir'
       print *, '                 -out  output_file'
       print *, '                 [-dir axis] '
       print *, '                 [-xce xcen] '
       print *, '                 [-yce ycen] '
       print *, '                 [-zce zcen] '
       print *, '                 [-rad rad] '
       print *, '                 [-xmi xmin] '
       print *, '                 [-xma xmax] '
       print *, '                 [-ymi ymin] '
       print *, '                 [-yma ymax] '
       print *, '                 [-zmi zmin] '
       print *, '                 [-zma zmax] '
       print *, '                 [-lma lmax] '
       print *, '                 [-typ type] '
       print *, '                 [-fil filetype] '
       print *, '                 [-max maxi] '
       print *, '                 [-gra read_grav] '
       print *, '                 [-pk  read_pk] '
       print *, '                 [-rt  read_rt] '
       print *, '                 [-nx  nx] '
       print *, '                 [-ny  ny] '
       print *, 'ex: amr2map -inp output_00001 -out map.dat'// &
            &   ' -dir z -xmi 0.1 -xma 0.7 -lma 12'
       print *, ' '
       print *, ' type :-1 = cpu number'
       print *, ' type : 0 = ref. level (default)'
       print *, ' type : 1 = gas density'
       print *, ' type : 2 = X velocity'
       print *, ' type : 3 = Y velocity'
       print *, ' type : 4 = Z velocity'
       print *, ' type : 5 = gas pressure'
       print *, ' type : 6 = gas metallicity'
       print *, ' '
       print *, ' maxi : 0 = average along line of sight (default)'
       print *, ' maxi : 1 = maximum along line of sight'
       print *, '        (always 1 for type -1 and 0)'
       stop
    end if

    do i = 1,n,2
       call getarg(i,opt)
       if (i == n) then
          print '("option ",a2," has no argument")', opt
          stop 2
       end if
       call getarg(i+1,arg)
       select case (opt)
       case ('-inp')
          repository = trim(arg)
       case ('-out')
          outfich = trim(arg)
       case ('-dir')
          proj = trim(arg)
       case ('-xce')
          read (arg,*) xcen
       case ('-yce')
          read (arg,*) ycen
       case ('-zce')
          read (arg,*) zcen
       case ('-rad')
          read (arg,*) rad
       case ('-xmi')
          read (arg,*) xmin
       case ('-xma')
          read (arg,*) xmax
       case ('-ymi')
          read (arg,*) ymin
       case ('-yma')
          read (arg,*) ymax
       case ('-zmi')
          read (arg,*) zmin
       case ('-zma')
          read (arg,*) zmax
       case ('-lma')
          read (arg,*) lmax
       case ('-nx')
          read (arg,*) nx_sample
       case ('-ny')
          read (arg,*) ny_sample
       case ('-typ')
          read (arg,*) type
       case ('-fil')
          read (arg,*) mapfiletype
       case ('-max')
          read (arg,*) domax
       case ('-gra')
          read (arg,*) dograv
       case ('-pk')
          read (arg,*) dopeak
       case ('-rt')
          read (arg,*) dort
       case default
          print '("unknown option ",a2," ignored")', opt
       end select
    end do

    do_max=.false.
    if(domax==1)do_max=.true.
    ! Refinement level and cpu number are never averaged
    if(type==0.or.type==-1)do_max=.true.
    do_grav=.false.
    if(dograv==1)do_grav=.true.
    do_peak=.false.
    if(dopeak==1)do_peak=.true.
    do_rt=.false.
    if(dort==1)do_rt=.true.

    if(rad>0)then
       xmin=xcen-rad
       xmax=xcen+rad
       ymin=ycen-rad
       ymax=ycen+rad
       zmin=zcen-rad
       zmax=zcen+rad
    endif

    return

  end subroutine read_params
  !================================================================
  !================================================================
  !================================================================
  !================================================================
  subroutine read_ramses_params
    !-----------------------------------------------
    ! Read RAMSES parameters file
    !-----------------------------------------------
    integer::ilun,ilevel,noutput,skip
    integer::nfile,ncpu

    nomfich=TRIM(repository)//'/params.bin'

    ilun=10
    open(unit=ilun,file=nomfich,access="stream",action="read",form='unformatted')
    read(ilun,POS=1)nfile
    read(ilun,POS=5)ncpu
    read(ilun,POS=9)p%ndim
    read(ilun,POS=13)p%levelmin
    read(ilun,POS=17)p%nlevelmax
    read(ilun,POS=21)p%boxlen
    read(ilun,POS=29)noutput
    skip=4*(11+4*noutput)+1
    read(ilun,POS=skip)p%t
    skip=skip+4*(2+4*p%nlevelmax+2+2*17)
    read(ilun,POS=skip)p%gamma
    skip=skip+8+48
    read(ilun,POS=skip)p%nhilbert
    allocate(p%bound_key(1:p%nhilbert,0:p%ncpu,p%levelmin:p%nlevelmax))
    skip=skip+4
    read(ilun,POS=skip)p%bound_key
    close(ilun)

    if(backup_file)then
       p%ncpu=ncpu
    else
       p%ncpu=nfile
    endif
    if(do_grav)then
       file_hydro=TRIM(repository)//'/grav.00001'
    else if(do_peak)then
       file_hydro=TRIM(repository)//'/peak_grid.00001'
    else if(do_rt)then
       file_hydro=TRIM(repository)//'/rt.00001'
    else
       file_hydro=TRIM(repository)//'/hydro.00001'
    endif
    open(unit=ilun,file=file_hydro,access="stream",action="read",form='unformatted')
    read(ilun,POS=5)p%nvar
    close(ilun)

  end subroutine read_ramses_params
  !================================================================
  !================================================================
  !================================================================
  !================================================================
end program amr2map

function check_ramses_exist(repository)
  !-----------------------------------------------
  ! Check that RAMSES files are there
  !-----------------------------------------------
  logical::check_ramses_exist
  character(len=*)::repository
  integer::ipos
  character(LEN=5)::char
  character(LEN=len(repository)+16)::nomfich
  character(LEN=len(repository)+16)::nomfich_grav
  logical::ok,ok_hydro,ok_grav

  check_ramses_exist=.true.
  ipos=INDEX(repository,'output_')
  nomfich_grav=TRIM(repository)//'/grav.00001'
  inquire(file=nomfich_grav, exist=ok_grav) ! verify input file
  nomfich=TRIM(repository)//'/hydro.00001'
  inquire(file=nomfich, exist=ok_hydro) ! verify input file
  if ( .not. ok_grav .and. .not. ok_hydro ) then
     print *,TRIM(nomfich_grav)//' not found.'
     print *,TRIM(nomfich)//' not found.'
     check_ramses_exist=.false.
  endif
  nomfich=TRIM(repository)//'/amr.00001'
  inquire(file=nomfich, exist=ok) ! verify input file
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found.'
     check_ramses_exist=.false.
  endif
  nomfich=TRIM(repository)//'/params.bin'
  inquire(file=nomfich, exist=ok) ! verify input file
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found.'
     check_ramses_exist=.false.
  endif

end function check_ramses_exist
!================================================================
!================================================================
!================================================================
!================================================================

!=======================================================================
subroutine title(n,nchar)
!=======================================================================
  implicit none
  integer::n
  character*5::nchar

  character*1::nchar1
  character*2::nchar2
  character*3::nchar3
  character*4::nchar4
  character*5::nchar5

  if(n.ge.10000)then
     write(nchar5,'(i5)') n
     nchar = nchar5
  elseif(n.ge.1000)then
     write(nchar4,'(i4)') n
     nchar = '0'//nchar4
  elseif(n.ge.100)then
     write(nchar3,'(i3)') n
     nchar = '00'//nchar3
  elseif(n.ge.10)then
     write(nchar2,'(i2)') n
     nchar = '000'//nchar2
  else
     write(nchar1,'(i1)') n
     nchar = '0000'//nchar1
  endif


end subroutine title

!!$subroutine cmp_cpu_list(xmin,xmax,ymin,ymax,zmin,zmax,ilevel,ncpu_read,cpu_list)
!!$  use amr_parameters, only: ndim,nhilbert,flen
!!$  implicit none
!!$  real(dp)::xmin,xmax,ymin,ymax,zmin,zmax
!!$  integer::ilevel,ncpu_read
!!$  integer,dimension(1:1)::cpu_list
!!$  !----------------------------------------------
!!$  ! Set up Hilbert key for optimal file reading
!!$  !----------------------------------------------
!!$  logical,dimension(1:ncpu)::cpu_read
!!$  integer(kind=4),dimension(1:nvector),save::dummy_state
!!$  integer(kind=8),dimension(1:nvector,1:ndim)::ix
!!$  integer(kind=8),dimension(1:nvector,1:nhilbert)::hk
!!$  integer(kind=8),dimension(1:nhilbert)::one_key,order_min,order_max
!!$  integer(kind=8),dimension(1:nhilbert,1:8)::bounding_min,bounding_max
!!$  integer,dimension(1:8)::idom,jdom,kdom,cpu_min,cpu_max
!!$  real(dp)::dmax,dx,maxdom
!!$  integer::i,j,ilev,lmin,imin,imax,jmin,jmax,kmin,kmax,bit_length,ndom,impi
!!$  logical::ok1,ok2
!!$
!!$  ncpu_read=1 !ncpu
!!$  do i=1,ncpu
!!$     cpu_list(i)=i
!!$  end do
!!$  return
!!$
!!$  dmax=max(xmax-xmin,ymax-ymin,zmax-zmin)
!!$  do ilev=1,nlevelmax
!!$     dx=0.5d0**ilev
!!$     if(dx.lt.dmax)exit
!!$  end do
!!$  lmin=MIN(ilev,ilevel)
!!$  bit_length=lmin-1
!!$  maxdom=2**bit_length
!!$  imin=0; imax=0; jmin=0; jmax=0; kmin=0; kmax=0
!!$  if(bit_length>0)then
!!$     imin=int(xmin*dble(maxdom))
!!$     imax=imin+1
!!$     jmin=int(ymin*dble(maxdom))
!!$     jmax=jmin+1
!!$     kmin=int(zmin*dble(maxdom))
!!$     kmax=kmin+1
!!$  endif
!!$
!!$  ndom=1
!!$  if(bit_length>0)ndom=8
!!$  idom(1)=imin; idom(2)=imax
!!$  idom(3)=imin; idom(4)=imax
!!$  idom(5)=imin; idom(6)=imax
!!$  idom(7)=imin; idom(8)=imax
!!$  jdom(1)=jmin; jdom(2)=jmin
!!$  jdom(3)=jmax; jdom(4)=jmax
!!$  jdom(5)=jmin; jdom(6)=jmin
!!$  jdom(7)=jmax; jdom(8)=jmax
!!$  kdom(1)=kmin; kdom(2)=kmin
!!$  kdom(3)=kmin; kdom(4)=kmin
!!$  kdom(5)=kmax; kdom(6)=kmax
!!$  kdom(7)=kmax; kdom(8)=kmax
!!$
!!$  one_key=0
!!$  one_key(1)=1
!!$  do i=1,ndom
!!$     if(bit_length>0)then
!!$        ix(1,1)=idom(i)
!!$        if(ndim>1)ix(1,2)=jdom(i)
!!$        if(ndim>2)ix(1,3)=kdom(i)
!!$        call hilbert_key(ix,hk,dummy_state,0,bit_length,1)
!!$        order_min(1:nhilbert)=hk(1,1:nhilbert)
!!$     else
!!$        order_min(1:nhilbert)=0
!!$     endif
!!$     bounding_min(1:nhilbert,i)=order_min
!!$     bounding_max(1:nhilbert,i)=order_min+one_key
!!$  end do
!!$
!!$  if(ilevel>bit_length+1)then
!!$     do ilev=bit_length+1,ilevel
!!$        do i=1,ndom
!!$           bounding_min(1:nhilbert,i)=refine_key(bounding_min(1:nhilbert,i),ilevel)
!!$        end do
!!$     end do
!!$  endif
!!$
!!$  cpu_min=0; cpu_max=0
!!$  do impi=1,ncpu
!!$     do i=1,ndom
!!$        ok1=ge_keys(bounding_min(1:nhilbert,i),bound_key_level(1:nhilbert,impi-1,ilevel))
!!$        ok2=gt_keys(bound_key_level(1:nhilbert,impi,ilevel),bounding_min(1:nhilbert,i))
!!$        if(ok1.and.ok2)then
!!$           cpu_min(i)=impi
!!$        endif
!!$        ok1=gt_keys(bounding_max(1:nhilbert,i),bound_key_level(1:nhilbert,impi-1,ilevel))
!!$        ok2=ge_keys(bound_key_level(1:nhilbert,impi,ilevel),bounding_max(1:nhilbert,i))
!!$        if(ok1.and.ok2)then
!!$           cpu_max(i)=impi
!!$        endif
!!$     end do
!!$  end do
!!$
!!$  ncpu_read=0
!!$  do i=1,ndom
!!$     do j=cpu_min(i),cpu_max(i)
!!$        if(.not. cpu_read(j))then
!!$           ncpu_read=ncpu_read+1
!!$           cpu_list(ncpu_read)=j
!!$           cpu_read(j)=.true.
!!$        endif
!!$     enddo
!!$  enddo
!!$
!!$end subroutine cmp_cpu_list
!================================================================
!================================================================
!================================================================
!================================================================
