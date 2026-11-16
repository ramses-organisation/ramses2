module input_hydro_grafic_module
contains
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_input_hydro_grafic(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::ilevel

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_INPUT_HYDRO_GRAFIC,pst%iUpper+1,input_size,0,ilevel)
     call r_input_hydro_grafic(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call input_hydro_grafic(pst%s%mdl,pst%s%r,pst%s%g,pst%s%m,ilevel)
  endif

end subroutine r_input_hydro_grafic
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine input_hydro_grafic(mdl,r,g,m,ilevel)
  use mdl_module
  use amr_parameters, only: ndim,twotondim,nvector,flen
  use hydro_parameters, only: nvar
  use amr_commons, only: run_t,global_t,mesh_t
  implicit none
  type(mdl_t)::mdl
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer::ilevel
  !--------------------------------------
  ! Compute initial conditions from files
  ! with the grafic format.
  !--------------------------------------
  integer::igrid,ilun
  integer::ind,idim,ivar
  integer::i1,i2,i3,i1_min,i1_max,i2_min,i2_max,i3_min,i3_max
  integer::buf_count

  real(kind=8)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(kind=8)::dx,rr,vx,vy=0,vz=0,ek,ei,pp,xx1,xx2,xx3,dx_loc
#ifdef MHD
  real(kind=8)::bx,by,bz,bdefault
#endif

  real(kind=8),allocatable,dimension(:,:,:)::init_array
  real(kind=4),allocatable,dimension(:,:)::init_plane

  logical::error,ok_file,ok_velb
  character(LEN=flen)::filename
  character(LEN=5)::ncharvar

  if(m%noct(ilevel)==0)return

  ! Conversion factor from user units to cgs units
  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Mesh size at level ilevel in normalised units
  dx=0.5D0**ilevel

  ! Mesh size at level ilevel in code units
  dx_loc=r%boxlen*dx

  !-------------------------------------------------------------------------
  ! First step: compute level boundaries in terms of initial condition array
  !-------------------------------------------------------------------------
  i1_min=g%n1(ilevel); i1_max=1
  i2_min=g%n2(ilevel); i2_max=1
  i3_min=g%n3(ilevel); i3_max=1
  do igrid=m%head(ilevel),m%tail(ilevel)
     do ind=1,twotondim
        ! Coordinates in normalised units (between 0 and 1)
        xx1=(2*m%grid(igrid)%ckey(1)+MOD((ind-1)  ,2)+0.5)*dx - m%skip(1)/r%boxlen
#if NDIM>1
        xx2=(2*m%grid(igrid)%ckey(2)+MOD((ind-1)/2,2)+0.5)*dx - m%skip(2)/r%boxlen
#endif
#if NDIM>2
        xx3=(2*m%grid(igrid)%ckey(3)+MOD((ind-1)/4,2)+0.5)*dx - m%skip(3)/r%boxlen
#endif
        ! Scale to integer coordinates in the frame of the file
        xx1=(xx1*(g%dxini(ilevel)/dx)-g%xoff1(ilevel))/g%dxini(ilevel)
#if NDIM>1
        xx2=(xx2*(g%dxini(ilevel)/dx)-g%xoff2(ilevel))/g%dxini(ilevel)
#endif
#if NDIM>2
        xx3=(xx3*(g%dxini(ilevel)/dx)-g%xoff3(ilevel))/g%dxini(ilevel)
#endif
        ! Compute min and max
        i1_min=MIN(i1_min,int(xx1)+1); i1_max=MAX(i1_max,int(xx1)+1)
#if NDIM>1
        i2_min=MIN(i2_min,int(xx2)+1); i2_max=MAX(i2_max,int(xx2)+1)
#endif
#if NDIM>2
        i3_min=MIN(i3_min,int(xx3)+1); i3_max=MAX(i3_max,int(xx3)+1)
#endif
     end do
  end do
  error=.false.
  if(i1_min<1.or.i1_max>g%n1(ilevel))error=.true.
  if(i2_min<1.or.i2_max>g%n2(ilevel))error=.true.
  if(i3_min<1.or.i3_max>g%n3(ilevel))error=.true.
  if(error) then
     write(*,*)'Some grid are outside initial conditions sub-volume'
     write(*,*)'for ilevel=',ilevel
     write(*,*)'and processor=',g%myid
     write(*,*)i1_min,i1_max
     write(*,*)i2_min,i2_max
     write(*,*)i3_min,i3_max
     write(*,*)g%n1(ilevel),g%n2(ilevel),g%n3(ilevel)
     call mdl_abort(mdl)
  end if

  !------------------------------------------
  ! Second step: read initial condition files
  !------------------------------------------
  ! Allocate initial conditions array
  allocate(init_array(i1_min:i1_max,i2_min:i2_max,i3_min:i3_max))
  allocate(init_plane(1:g%n1(ilevel),1:g%n2(ilevel)))

  ! Loop over input variables
  do ivar=1,nvar
     if(r%cosmo)then
        ! Read baryons initial overdensity and displacement at a=aexp
        if(ivar==1)filename=TRIM(r%initfile(ilevel))//'/ic_deltab'
        if(ivar==2)then
           filename=TRIM(r%initfile(ilevel))//'/ic_velbx'
           INQUIRE(file=filename,exist=ok_velb)
           if(.not. ok_velb)then
              filename=TRIM(r%initfile(ilevel))//'/ic_velcx'
           endif
        endif
        if(ivar==3)then
           filename=TRIM(r%initfile(ilevel))//'/ic_velby'
           INQUIRE(file=filename,exist=ok_velb)
           if(.not. ok_velb)then
              filename=TRIM(r%initfile(ilevel))//'/ic_velcy'
           endif
        endif
        if(ivar==4)then
           filename=TRIM(r%initfile(ilevel))//'/ic_velbz'
           INQUIRE(file=filename,exist=ok_velb)
           if(.not. ok_velb)then
              filename=TRIM(r%initfile(ilevel))//'/ic_velcz'
           endif
        endif
        if(ivar==5)filename=TRIM(r%initfile(ilevel))//'/ic_tempb'
     else
        ! Read primitive variables
        if(ivar==1)filename=TRIM(r%initfile(ilevel))//'/ic_d'
        if(ivar==2)filename=TRIM(r%initfile(ilevel))//'/ic_u'
        if(ivar==3)filename=TRIM(r%initfile(ilevel))//'/ic_v'
        if(ivar==4)filename=TRIM(r%initfile(ilevel))//'/ic_w'
        if(ivar==5)filename=TRIM(r%initfile(ilevel))//'/ic_p'
     endif
     call title(ivar,ncharvar)
     if(ivar>5)then
        call title(ivar-5,ncharvar)
        filename=TRIM(r%initfile(ilevel))//'/ic_pvar_'//TRIM(ncharvar)
     endif

     INQUIRE(file=filename,exist=ok_file)
     if(ok_file)then
        ! Reading the existing file
        if(g%myid==1)write(*,*)"Reading "//TRIM(filename)
        open(10,file=filename,form='unformatted')
        rewind 10
        read(10) ! skip first line
        do i3=1,i3_min-1
           read(10)
        end do
        do i3=i3_min,i3_max
           read(10) ((init_plane(i1,i2),i1=1,g%n1(ilevel)),i2=1,g%n2(ilevel))
           init_array(i1_min:i1_max,i2_min:i2_max,i3) = init_plane(i1_min:i1_max,i2_min:i2_max)
        end do
        close(10)
     else
        ! If file doesn't exist, initialize variable to default value
        ! In most cases, this is zero (you can change that if necessary)
        if(g%myid==1)write(*,*)"Missing "//TRIM(filename)
        init_array=0d0
        if(r%metal.and.ivar==r%imetal)then
           init_array=r%z_ave*0.02
        endif
     endif

     ! For cosmo runs, rescale initial conditions to code units
     if(r%cosmo)then
        if(.not. r%cooling .and. .not. r%neq_chem)then
           ! Compute approximate average temperature in K
           g%T2_start=1.356d-2/g%aexp**2
           if(g%myid==1)write(*,*)"Approximate temperature assumed using an adiabatic evolution, T/mu (K) = ",g%T2_start
        endif
        if(ivar==1)init_array=(1.0+g%dfact(ilevel)*init_array)*g%omega_b/g%omega_m
        if(ivar==2)init_array=g%dfact(ilevel)*g%vfact(1)*dx_loc/g%dxini(ilevel)*init_array/g%vfact(ilevel)
        if(ivar==3)init_array=g%dfact(ilevel)*g%vfact(1)*dx_loc/g%dxini(ilevel)*init_array/g%vfact(ilevel)
        if(ivar==4)init_array=g%dfact(ilevel)*g%vfact(1)*dx_loc/g%dxini(ilevel)*init_array/g%vfact(ilevel)
        if(ivar==5)init_array=(1.0+init_array)*g%T2_start/scale_T2
     endif

     ! Loop over cells
     do igrid=m%head(ilevel),m%tail(ilevel)
        do ind=1,twotondim
           ! Coordinates in normalised units (between 0 and 1)
           xx1=(2*m%grid(igrid)%ckey(1)+MOD((ind-1)  ,2)+0.5)*dx - m%skip(1)/r%boxlen
#if NDIM>1
           xx2=(2*m%grid(igrid)%ckey(2)+MOD((ind-1)/2,2)+0.5)*dx - m%skip(2)/r%boxlen
#endif
#if NDIM>2
           xx3=(2*m%grid(igrid)%ckey(3)+MOD((ind-1)/4,2)+0.5)*dx - m%skip(3)/r%boxlen
#endif
           ! Scale to integer coordinates in the frame of the file
           xx1=(xx1*(g%dxini(ilevel)/dx)-g%xoff1(ilevel))/g%dxini(ilevel)
#if NDIM>1
           xx2=(xx2*(g%dxini(ilevel)/dx)-g%xoff2(ilevel))/g%dxini(ilevel)
#endif
#if NDIM>2
           xx3=(xx3*(g%dxini(ilevel)/dx)-g%xoff3(ilevel))/g%dxini(ilevel)
#endif
           ! Compute integer coordinates in the frame of the file
           i1=int(xx1)+1
           i2=1; i3=1
#if NDIM>1
           i2=int(xx2)+1
#endif
#if NDIM>2
           i3=int(xx3)+1
#endif
#ifdef HYDRO
           ! Scatter to corresponding primitive variable
           m%uold(ind,ivar,igrid)=init_array(i1,i2,i3)
#endif
        end do
     end do
     ! End loop over cells

  end do
  ! End loop over input variables

  !------------------------------------------
  ! Read magnetic field initial condition files (MHD)
  !------------------------------------------
#ifdef MHD
  do ivar=1,6
     ! Map filename per face component
     if(ivar==1) filename=TRIM(r%initfile(ilevel))//'/ic_bxleft'
     if(ivar==2) filename=TRIM(r%initfile(ilevel))//'/ic_byleft'
     if(ivar==3) filename=TRIM(r%initfile(ilevel))//'/ic_bzleft'
     if(ivar==4) filename=TRIM(r%initfile(ilevel))//'/ic_bxright'
     if(ivar==5) filename=TRIM(r%initfile(ilevel))//'/ic_byright'
     if(ivar==6) filename=TRIM(r%initfile(ilevel))//'/ic_bzright'

     ! Set default background B by component group
     if(ivar==1 .or. ivar==4) bdefault=r%A_ave
     if(ivar==2 .or. ivar==5) bdefault=r%B_ave
     if(ivar==3 .or. ivar==6) bdefault=r%C_ave

     INQUIRE(file=filename,exist=ok_file)
     if(ok_file)then
        if(g%myid==1)write(*,*)"Reading "//TRIM(filename)
        open(10,file=filename,form='unformatted')
        rewind 10
        read(10)
        do i3=1,i3_min-1
           read(10)
        end do
        do i3=i3_min,i3_max
           read(10) ((init_plane(i1,i2),i1=1,g%n1(ilevel)),i2=1,g%n2(ilevel))
           init_array(i1_min:i1_max,i2_min:i2_max,i3) = init_plane(i1_min:i1_max,i2_min:i2_max)
        end do
        close(10)
     else
        if(g%myid==1)write(*,*)"Missing "//TRIM(filename)
        init_array = bdefault
     endif

     do igrid=m%head(ilevel),m%tail(ilevel)
        do ind=1,twotondim
           xx1=(2*m%grid(igrid)%ckey(1)+MOD((ind-1)  ,2)+0.5)*dx - m%skip(1)/r%boxlen
#if NDIM>1
           xx2=(2*m%grid(igrid)%ckey(2)+MOD((ind-1)/2,2)+0.5)*dx - m%skip(2)/r%boxlen
#endif
#if NDIM>2
           xx3=(2*m%grid(igrid)%ckey(3)+MOD((ind-1)/4,2)+0.5)*dx - m%skip(3)/r%boxlen
#endif
           xx1=(xx1*(g%dxini(ilevel)/dx)-g%xoff1(ilevel))/g%dxini(ilevel)
#if NDIM>1
           xx2=(xx2*(g%dxini(ilevel)/dx)-g%xoff2(ilevel))/g%dxini(ilevel)
#endif
#if NDIM>2
           xx3=(xx3*(g%dxini(ilevel)/dx)-g%xoff3(ilevel))/g%dxini(ilevel)
#endif
           i1=int(xx1)+1
           i2=1; i3=1
#if NDIM>1
           i2=int(xx2)+1
#endif
#if NDIM>2
           i3=int(xx3)+1
#endif
           m%bold(ind,ivar,igrid)=init_array(i1,i2,i3)
        end do
     end do
  end do
#endif

  ! Deallocate initial conditions array
  deallocate(init_array)
  deallocate(init_plane)

  !----------------------------------------------------------------
  ! For cosmology runs: compute pressure, prevent negative density
  !----------------------------------------------------------------
  if(r%cosmo)then
     ! Loop over grids
     do igrid=m%head(ilevel),m%tail(ilevel)
        ! Loop over cells
        do ind=1,twotondim
#ifdef HYDRO
           ! Prevent negative densities
           rr=max(dble(m%uold(ind,1,igrid)),0.1*g%omega_b/g%omega_m)
           m%uold(ind,1,igrid)=rr
           ! Compute pressure from temperature and density
           m%uold(ind,5,igrid)=m%uold(ind,1,igrid)*m%uold(ind,5,igrid)
#endif
        end do
        ! End loop over cells
     end do
     ! End loop over grids
  end if

  !-------------------------------------
  ! If required, compute initial entropy
  !-------------------------------------
  if(r%entropy)then
     ! Loop over grids
     do igrid=m%head(ilevel),m%tail(ilevel)
        ! Loop over cells
        do ind=1,twotondim
#ifdef HYDRO
           ! Compute entropy from pressure and density
           m%uold(ind,r%ientropy,igrid)=m%uold(ind,5,igrid)/m%uold(ind,1,igrid)**r%gamma
#endif
        end do
        ! End loop over cells
     end do
     ! End loop over grids
  end if

  !-----------------------------------------
  ! If required, compute refinement map
  !-----------------------------------------
#ifdef GRAV
#ifdef HYDRO
  if(r%ivar_refine>0)then
     ! Loop over grids
     do igrid=m%head(ilevel),m%tail(ilevel)
        ! Loop over cells
        do ind=1,twotondim
           ! Compute initial refinement map for zoom-in simulations
           ! only if next level file exists
           if(r%initfile(ilevel+1) .ne.' ')then
              if(m%uold(ind,r%ivar_refine,igrid)>r%var_cut_refine)then
                 m%nref(ind,igrid)=r%m_refine(ilevel)*1.1d0
              else
                 m%nref(ind,igrid)=0.0d0
              endif
           else
              m%nref(ind,igrid)=0.0d0
           endif
        end do
        ! End loop over cells
     end do
     ! End loop over grids
  end if
#endif
#endif

  !---------------------------------------------------
  ! Third step: compute initial conservative variables
  !---------------------------------------------------
  ! Loop over grids
  do igrid=m%head(ilevel),m%tail(ilevel)
     ! Loop over cells
     do ind=1,twotondim
#ifdef HYDRO
        ! Compute total energy density
        rr=m%uold(ind,1,igrid)
        vx=m%uold(ind,2,igrid)
        vy=m%uold(ind,3,igrid)
        vz=m%uold(ind,4,igrid)
        pp=m%uold(ind,5,igrid)
        ek=0.5d0*rr*(vx**2+vy**2+vz**2)
        ei=pp/(r%gamma-1.0)
#ifdef MHD
        ! Add magnetic energy using face-averaged, cell-centered B
        bx=0.5d0*(m%bold(ind,1,igrid)+m%bold(ind,4,igrid))
        by=0.5d0*(m%bold(ind,2,igrid)+m%bold(ind,5,igrid))
        bz=0.5d0*(m%bold(ind,3,igrid)+m%bold(ind,6,igrid))
        m%uold(ind,5,igrid)=ei+ek+0.5d0*(bx*bx+by*by+bz*bz)
#else
        m%uold(ind,5,igrid)=ei+ek
#endif
        ! Compute momentum density
        do idim=1,3
           rr=m%uold(ind,1,igrid)
           m%uold(ind,idim+1,igrid)=rr*m%uold(ind,idim+1,igrid)
        end do
#if NVAR>5
        ! Compute passive scalar density
        do ivar=6,nvar
           rr=m%uold(ind,1,igrid)
           m%uold(ind,ivar,igrid)=rr*m%uold(ind,ivar,igrid)
        enddo
#endif
#endif
     end do
     ! End loop over cells
  end do
  ! End loop over grids

end subroutine input_hydro_grafic
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_input_refmap_grafic(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::ilevel
  !--------------------------------------------------------------------
  ! This routine reads the refinbement map or mask from a grafic file
  ! for zoom-in simulation when there is no corresponding hydro variable.
  !--------------------------------------------------------------------
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_INPUT_REFMAP_GRAFIC,pst%iUpper+1,input_size,0,ilevel)
     call r_input_refmap_grafic(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call input_refmap_grafic(pst%s%mdl,pst%s%r,pst%s%g,pst%s%m,ilevel)
  endif

end subroutine r_input_refmap_grafic
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine input_refmap_grafic(mdl,r,g,m,ilevel)
  use mdl_module
  use amr_parameters, only: ndim,twotondim,nvector,flen
  use hydro_parameters, only: nvar
  use amr_commons, only: run_t,global_t,mesh_t
  implicit none
  type(mdl_t)::mdl
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer::ilevel
  !--------------------------------------
  ! Read refinement map or mask  from a
  ! file with the grafic format.
  !--------------------------------------
  integer::igrid,ilun
  integer::ind,idim
  integer::i1,i2,i3,i1_min,i1_max,i2_min,i2_max,i3_min,i3_max
  integer::buf_count

  real(kind=8)::dx,xx1,xx2,xx3

  real(kind=8),allocatable,dimension(:,:,:)::init_array
  real(kind=4),allocatable,dimension(:,:)::init_plane

  logical::error,ok_file
  character(LEN=flen)::filename

  ! Report refinement mask status (before the early return: rank 1 may hold no oct)
  filename=TRIM(r%initfile(ilevel))//'/ic_refmap'
  INQUIRE(file=filename,exist=ok_file)
  ok_file = ok_file .and. (r%initfile(ilevel+1) .ne.' ')
  if(g%myid==1)then
     if(ok_file)write(*,*)"Reading "//TRIM(filename)
     ! If a finer level has no mask, this level is not zoomed
     if(.not.ok_file.and.r%initfile(ilevel+1).ne.' ')write(*,*)"Missing "//TRIM(filename)
  endif

  if(m%noct(ilevel)==0)return

#if NDIM==3

  ! Mesh size at level ilevel in normalised units
  dx=0.5D0**ilevel

  !-------------------------------------------------------------------------
  ! First step: compute level boundaries in terms of initial condition array
  !-------------------------------------------------------------------------
  i1_min=g%n1(ilevel)+1; i1_max=0
  i2_min=g%n2(ilevel)+1; i2_max=0
  i3_min=g%n3(ilevel)+1; i3_max=0
  do igrid=m%head(ilevel),m%tail(ilevel)
     do ind=1,twotondim
        ! Coordinates in normalised units (between 0 and 1)
        xx1=(2*m%grid(igrid)%ckey(1)+MOD((ind-1)  ,2)+0.5)*dx
        xx2=(2*m%grid(igrid)%ckey(2)+MOD((ind-1)/2,2)+0.5)*dx
        xx3=(2*m%grid(igrid)%ckey(3)+MOD((ind-1)/4,2)+0.5)*dx
        ! Scale to integer coordinates in the frame of the file
        xx1=(xx1*(g%dxini(ilevel)/dx)-g%xoff1(ilevel))/g%dxini(ilevel)
        xx2=(xx2*(g%dxini(ilevel)/dx)-g%xoff2(ilevel))/g%dxini(ilevel)
        xx3=(xx3*(g%dxini(ilevel)/dx)-g%xoff3(ilevel))/g%dxini(ilevel)
        ! Compute min and max
        i1_min=MIN(i1_min,int(xx1)+1); i1_max=MAX(i1_max,int(xx1)+1)
        i2_min=MIN(i2_min,int(xx2)+1); i2_max=MAX(i2_max,int(xx2)+1)
        i3_min=MIN(i3_min,int(xx3)+1); i3_max=MAX(i3_max,int(xx3)+1)
     end do
  end do
  error=.false.
  if(i1_min<1.or.i1_max>g%n1(ilevel))error=.true.
  if(i2_min<1.or.i2_max>g%n2(ilevel))error=.true.
  if(i3_min<1.or.i3_max>g%n3(ilevel))error=.true.
  if(error) then
     write(*,*)'Some grid are outside initial conditions sub-volume'
     write(*,*)'for ilevel=',ilevel
     write(*,*)'and processor=',g%myid
     write(*,*)i1_min,i1_max
     write(*,*)i2_min,i2_max
     write(*,*)i3_min,i3_max
     write(*,*)g%n1(ilevel),g%n2(ilevel),g%n3(ilevel)
     call mdl_abort(mdl)
  end if

  !------------------------------------------
  ! Second step: read initial condition files
  !------------------------------------------
  ! Allocate initial conditions array
  allocate(init_array(i1_min:i1_max,i2_min:i2_max,i3_min:i3_max))
  allocate(init_plane(1:g%n1(ilevel),1:g%n2(ilevel)))

  ! Read refinement mask
  if(ok_file)then
     ! Reading the existing file
     open(10,file=filename,form='unformatted')
     rewind 10
     read(10) ! skip first line
     do i3=1,i3_min-1
        read(10)
     end do
     do i3=i3_min,i3_max
        read(10) ((init_plane(i1,i2),i1=1,g%n1(ilevel)),i2=1,g%n2(ilevel))
        init_array(i1_min:i1_max,i2_min:i2_max,i3) = init_plane(i1_min:i1_max,i2_min:i2_max)
     end do
     close(10)
  else
     init_array(i1_min:i1_max,i2_min:i2_max,i3_min:i3_max) = 0.0d0
  endif

  ! Loop over cells
  do igrid=m%head(ilevel),m%tail(ilevel)
     do ind=1,twotondim
        ! Coordinates in normalised units (between 0 and 1)
        xx1=(2*m%grid(igrid)%ckey(1)+MOD((ind-1)  ,2)+0.5)*dx
        xx2=(2*m%grid(igrid)%ckey(2)+MOD((ind-1)/2,2)+0.5)*dx
        xx3=(2*m%grid(igrid)%ckey(3)+MOD((ind-1)/4,2)+0.5)*dx
        ! Scale to integer coordinates in the frame of the file
        xx1=(xx1*(g%dxini(ilevel)/dx)-g%xoff1(ilevel))/g%dxini(ilevel)
        xx2=(xx2*(g%dxini(ilevel)/dx)-g%xoff2(ilevel))/g%dxini(ilevel)
        xx3=(xx3*(g%dxini(ilevel)/dx)-g%xoff3(ilevel))/g%dxini(ilevel)
        ! Compute integer coordinates in the frame of the file
        i1=int(xx1)+1
        i2=int(xx2)+1
        i3=int(xx3)+1
#ifdef GRAV
        ! Scatter to corresponding refinement variable
        if(r%initfile(ilevel+1) .ne.' ')then
           m%nref(ind,igrid)=init_array(i1,i2,i3)*r%m_refine(ilevel)*1.1d0
        else
           m%nref(ind,igrid)=0d0
        endif
#endif
     end do
  end do
  ! End loop over cells

  ! Deallocate initial conditions array
  deallocate(init_array)
  deallocate(init_plane)
#endif
end subroutine input_refmap_grafic
!################################################################
!################################################################
!################################################################
!################################################################
end module input_hydro_grafic_module
