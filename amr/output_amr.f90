module output_amr_module
contains
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_dump_all(pst,write_bkp_file)
  use amr_parameters, only: ndim,flen
  use ramses_commons, only: pst_t
  use output_hydro_module, only: r_output_hydro,file_descriptor_hydro
  use output_poisson_module, only: r_output_poisson,in_output_poisson_t,file_descriptor_poisson
  use output_part_module, only: r_output_part
  use mdl_module, only: mdl_mkdir, mdl_wtime
  use cooling_module, only: output_cool
  use output_rt_module, only: r_output_rt,file_descriptor_rt,output_rtinfo
  use output_cr_module, only: r_output_cr,file_descriptor_cr,output_crinfo
  use turb_commons, only: write_turb_fields
  implicit none
  type(pst_t)::pst
  logical::write_bkp_file

  ! Local variables
  integer::i,dummy(1)
#ifdef NOSYSTEM
  integer::ierr
#endif
  double precision::ttend, ttstart=0.0
  character(LEN=5)::nchar
  character(LEN=flen)::filename,filedir,filecmd
  integer,dimension(1:flen/4)::input_array
  type(in_output_poisson_t)::in_output_poisson

  associate(r=>pst%s%r,g=>pst%s%g,m=>pst%s%m,p=>pst%s%p,star=>pst%s%star, &
   & sink=>pst%s%sink,tree=>pst%s%tree,trac=>pst%s%trac,dust=>pst%s%dust,mdl=>pst%s%mdl)

  if(g%nstep_coarse==g%nstep_coarse_old.and.g%nstep_coarse>0)return
  if(g%nstep_coarse==0.and.r%nrestart>0)return

  if(r%verbose)then
     write(*,*)'Entering dump_all'
     if(write_bkp_file)then
        write(*,*)'Writing restart files'
     else
        write(*,*)'Writing output files'
     endif
  endif
  ! For 1D serial runs, output data to screen
  if(ndim==1.and.g%ncpu==1)call write_screen(r,m)

  if(write_bkp_file)then
     ! Increment backup file counter
     if(r%bkp_modulo>0)g%ifbkp=MOD(g%ifbkp-1,r%bkp_modulo)+1
     call title(g%ifbkp,nchar)
     g%ifbkp=g%ifbkp+1
  else
     ! Increment output file counter
     call title(g%ifout,nchar)
     g%ifout=g%ifout+1
     if(g%t>=r%tout(g%iout).or.g%aexp>=r%aout(g%iout))g%iout=g%iout+1
     g%output_done=.true.
  endif

  ! For 2D and 3D runs, output data to files
  if(ndim>1.or.g%ncpu>1)then
     if(write_bkp_file)then
        filedir='backup_'//TRIM(nchar)//'/'
     else
        filedir='output_'//TRIM(nchar)//'/'
     endif
     call mdl_mkdir(mdl,filedir)
     ! filecmd='mkdir -p '//TRIM(filedir)
! #ifdef NOSYSTEM
!      call PXFMKDIR(TRIM(filedir),LEN(TRIM(filedir)),O'755',ierr)
! #else
!      call system(filecmd)
! #endif

     !-----------------------
     ! Only master process
     !-----------------------
     if(write_bkp_file)then
        write(*,*)'Writing backup file to disk in '//TRIM(filedir)
     else
        write(*,*)'Writing output file to disk in '//TRIM(filedir)
     endif
     ttstart = mdl_wtime(mdl)

     if(r%verbose)write(*,*)'Writing header files'
     if(r%part)then
        filename=TRIM(filedir)//'part_header.txt'
        call output_header(r,g,p   ,filename)
     endif
     if(r%star)then
        filename=TRIM(filedir)//'star_header.txt'
        call output_header(r,g,star,filename)
     endif
     if(r%sink)then
        filename=TRIM(filedir)//'sink_header.txt'
        call output_header(r,g,sink,filename)
     endif
     if(r%tree)then
        filename=TRIM(filedir)//'tree_header.txt'
        call output_header(r,g,tree,filename)
     endif
     if(r%trac)then
        filename=TRIM(filedir)//'trac_header.txt'
        call output_header(r,g,trac,filename)
     endif
     if(r%dust)then
        filename=TRIM(filedir)//'dust_header.txt'
        call output_header(r,g,dust,filename)
     endif
     if(r%hydro)then
        filename=TRIM(filedir)//'hydro_header.txt'
        call file_descriptor_hydro(r,filename,write_bkp_file)
     end if
     if(r%rt)then
        filename=TRIM(filedir)//'rt_header.txt'
        call file_descriptor_rt(r,filename,write_bkp_file)
        filename=TRIM(filedir)//'rt_info.txt'
        call output_rtinfo(r,g,filename)
     end if
     if(r%cr)then
        filename=TRIM(filedir)//'cr_header.txt'
        call file_descriptor_cr(r,filename)
        filename=TRIM(filedir)//'cr_info.txt'
        call output_crinfo(r,g,filename)
     end if
     if(r%poisson)then
        filename=TRIM(filedir)//'grav_header.txt'
        call file_descriptor_poisson(r,filename,write_bkp_file)
     end if
     if(r%cooling)then
        filename=TRIM(filedir)//'cooling.bin'
        call output_cool(pst%s%cool,filename)
     end if
     if(r%turb)then
        call write_turb_fields(pst%s%r,pst%s%turb,filedir)
     endif
     filename=TRIM(filedir)//'info.txt'
     call output_info(r,g,filename)
     filename=TRIM(filedir)//'namelist.txt'
     call output_namelist(filename)
     filename=TRIM(filedir)//'compilation.txt'
     call output_compil(filename)
     filename=TRIM(filedir)//'params.bin'
     call output_params(r,g,m,filename)
     filename=TRIM(filedir)//'timer.txt'
     call m_output_timer(.true.,filename)

     !-----------------------
     ! All slave processes
     !-----------------------

     ! Output AMR data
     filename=TRIM(filedir)//'amr.'
     input_array=transfer(filename,input_array)
     ! Only write AMR data if output_amr is true OR if this is a backup
     if(r%output_amr .or. write_bkp_file)then
        if(r%verbose)write(*,*)'Writing AMR files'
        call r_output_amr(pst,input_array,flen/4,dummy,0)
     else
        if(r%verbose)write(*,*)'Skipping AMR files'
     endif

     ! Output HYDRO data
     if(r%hydro)then
        filename=TRIM(filedir)//'hydro.'
        input_array=transfer(filename,input_array)
        ! Only write hydro data if output_hydro is true OR if this is a backup
        if(r%output_hydro .or. write_bkp_file)then
           if(r%verbose)write(*,*)'Writing hydro files'
           call r_output_hydro(pst,input_array,flen/4,dummy,0)
        else
           if(r%verbose)write(*,*)'Skipping hydro files'
        endif
     end if

     ! Output GRAV data
     if(r%poisson)then
        filename=TRIM(filedir)//'grav.'
        input_array=transfer(filename,input_array)
        in_output_poisson%filename=TRIM(filedir)//'grav.'
        ! Only write gravity data if output_grav is true OR if this is a backup
        if(r%output_grav .or. write_bkp_file)then
           if(r%verbose)write(*,*)'Writing gravity files'
           call r_output_poisson(pst,in_output_poisson,storage_size(in_output_poisson)/32)
        else
           if(r%verbose)write(*,*)'Skipping gravity files'
        endif
     end if

     ! Output PART data
     if(r%pic)then
        filename=TRIM(filedir) ! Note that suffix will be added later
        input_array=transfer(filename,input_array)
        if(r%verbose)write(*,*)'Writing particle files'
        call r_output_part(pst,input_array,flen/4,dummy,0)
     end if

     ! Output RT data
     if(r%rt)then
        filename=TRIM(filedir)//'rt.'
        input_array=transfer(filename,input_array)
        if(r%verbose)write(*,*)'Writing RT files'
        call r_output_rt(pst,input_array,flen/4,dummy,0)
     end if

     ! Output CR data
     if(r%cr)then
        filename=TRIM(filedir)//'cr.'
        input_array=transfer(filename,input_array)
        if(r%verbose)write(*,*)'Writing CR files'
        call r_output_cr(pst,input_array,flen/4,dummy,0)
     end if

     ttend = mdl_wtime(mdl)
     print '(A,F0.7)',' Time elapsed writing to disk: ',ttend-ttstart

  end if

  end associate

end subroutine m_dump_all
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine output_namelist(filename)
  use amr_parameters, only: flen
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
  character(LEN=flen)::filename
  ! Copy namelist file to output directory
  character::nml_char
  integer::ierr

  open(10,file=namelist_file,access="stream",action="read")
  open(11,file=filename,access="stream",action="write")
  do
     read(10,iostat=ierr)nml_char
     if(ierr.NE.0)exit
     write(11)nml_char
  end do
  close(10)
  close(11)

end subroutine output_namelist
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine output_compil(filename)
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
  real(dp)::real_size
  character(LEN=flen)::filename
  character(LEN=30)::mystring
  ! Copy compilation details to output directory
  OPEN(UNIT=11, FILE=filename, FORM='formatted')
  write(11,'(" compile date = ",A)')TRIM(builddate)
  write(11,'(" compile command = ",A)')TRIM(buildcommand)
  write(11,'(" remote repo = ",A)')TRIM(gitrepo)
  write(11,'(" local branch = ",A)')TRIM(gitbranch)
  write(11,'(" last commit = ",A)')TRIM(githash)
  write(11,'(" PATCH = ",A)')TRIM(patchdir)
  write(mystring,*)nhilbert
  write(11,'(" NHILBERT = ",A)')adjustl(mystring)
  write(mystring,*)nvector
  write(11,'(" NVECTOR = ",A)')adjustl(mystring)
  write(mystring,*)ndim
  write(11,'(" NDIM = ",A)')adjustl(mystring)
  if(sizeof(real_size)==8)then
     write(11,'(" NPRE = 8")')
  else
     write(11,'(" NPRE = 4")')
  endif
  write(mystring,*)nvar
  write(11,'(" NVAR = ",A)')adjustl(mystring)
  write(mystring,*)nener
  write(11,'(" NENER = ",A)')adjustl(mystring)
#ifdef HYDRO
  write(11,'(" HYDRO = 1")')
#else
  write(11,'(" HYDRO = 0")')
#endif
#ifdef GRAV
  write(11,'(" GRAV = 1")')
#else
  write(11,'(" GRAV = 0")')
#endif
#ifdef MHD
  write(11,'(" MHD = 1")')
#else
  write(11,'(" MHD = 0")')
#endif
#ifdef WITHOUTMPI
  write(11,'(" MPI = 0")')
#else
  write(11,'(" MPI = 1")')
#endif
#define COEUR 1
#define COSMO 2
#define MERGER 3
#if UNITS==COSMO
  write(11,'(" UNITS = COSMO")')
#elif UNITS==COEUR
  write(11,'(" UNITS = COEUR")')
#elif UNITS==MERGER
  write(11,'(" UNITS = MERGER")')
#else
  write(11,'(" UNITS = ")')
#endif
#define INSTA 2
#define DOUBLEMACH 3
#define LOOP 4
#define OT 5
#define PONO 6
#if INIT==COEUR
  write(11,'(" INIT = COEUR")')
#elif INIT==INSTA
  write(11,'(" INIT = INSTA")')
#elif INIT==DOUBLEMACH
  write(11,'(" INIT = DOUBLEMACH")')
#elif INIT==LOOP
  write(11,'(" INIT = LOOP")')
#elif INIT==OT
  write(11,'(" INIT = OT")')
#elif INIT==PONO
  write(11,'(" INIT = PONO")')
#else
  write(11,'(" INIT = ")')
#endif

  CLOSE(11)
end subroutine output_compil
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine output_params(r,g,m,filename)
  use amr_parameters, only: ndim,nhilbert,flen
  use amr_commons, only: run_t,global_t,mesh_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  character(LEN=flen)::filename

  ! Local variables
  integer::ilun,ilevel
  character(LEN=flen)::fileloc

  if(r%verbose)write(*,*)'Entering output_params'

  !-----------------------------------
  ! Output run parameters in file
  !-----------------------------------
  ilun=10
  fileloc=TRIM(filename)
  open(unit=ilun,file=fileloc,access="stream"&
       & ,action="write",form='unformatted')
  ! Write grid variables
  write(ilun)r%nfile
  write(ilun)g%ncpu
  write(ilun)ndim
  write(ilun)r%levelmin
  write(ilun)r%nlevelmax
  write(ilun)r%boxlen
  ! Write time variables
  write(ilun)r%noutput,g%iout,g%ifout,g%ifbkp
  write(ilun)r%tout(1:r%noutput)
  write(ilun)r%aout(1:r%noutput)
  write(ilun)g%t
  write(ilun)g%dtold(1:r%nlevelmax)
  write(ilun)g%dtnew(1:r%nlevelmax)
  write(ilun)g%nstep,g%nstep_coarse
  ! Write various constants
  write(ilun)g%const,g%mass_tot_0,g%rho_tot
  write(ilun)g%omega_m,g%omega_l,g%omega_k,g%omega_b,g%h0,g%aexp_ini,g%boxlen_ini
  write(ilun)g%aexp,g%texp,g%hexp
  write(ilun)g%aexp_old,g%epot_tot_int,g%epot_tot_old
  write(ilun)r%mass_sph
  write(ilun)r%gamma
  ! Write RNG latest seed
  write(ilun)r%seed
  ! Write cpu boundaries
  write(ilun)nhilbert
  do ilevel=r%levelmin,r%nlevelmax
     write(ilun)m%domain(ilevel)%b(1:nhilbert,0:g%ncpu)
  end do
  close(ilun)

end subroutine output_params
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine input_params(mdl,r,g,filename,nfile_file,ncpu_file,levelmin_file,nlevelmax_file)
  use mdl_module
  use amr_parameters, only: ndim, nhilbert, flen
  use amr_commons, only: run_t, global_t
  implicit none
  type(mdl_t)::mdl
  type(run_t)::r
  type(global_t)::g
  character(LEN=flen)::filename
  integer::nfile_file,ncpu_file,levelmin_file,nlevelmax_file,i
  !-----------------------------------
  ! Read run parameters from file.
  ! Note that ncpu, levelmin and nlevelmax
  ! are allowed to vary at restart.
  !-----------------------------------
  integer::ilun
  integer::ndim_file,noutput_file
  integer::noutput_min,nlevelmax_min
  real(kind=8)::mass_sph_file,gamma_file
  real(kind=8),allocatable,dimension(:)::tout,aout,dtold,dtnew
  character(LEN=flen)::fileloc

  if(r%verbose)write(*,*)'Entering input_params'

  ilun=10
  fileloc=TRIM(filename)
  open(unit=ilun,file=fileloc,access="stream"&
       & ,action="read",form='unformatted')
  ! Read grid variables
  read(ilun)nfile_file
  read(ilun)ncpu_file
  read(ilun)ndim_file
  read(ilun)levelmin_file
  read(ilun)nlevelmax_file
  ! Overwrite boxlen with value from file
  read(ilun)r%boxlen
  ! Read time variables
  read(ilun)noutput_file,g%iout,g%ifout,g%ifbkp
  noutput_min=MIN(r%noutput,noutput_file)
  allocate(tout(1:noutput_file))
  allocate(aout(1:noutput_file))
  read(ilun)tout
  read(ilun)aout
  if(noutput_min>0)then
     r%tout(1:noutput_min)=tout(1:noutput_min)
     r%aout(1:noutput_min)=aout(1:noutput_min)
  endif
  deallocate(tout,aout)
  read(ilun)g%t
  nlevelmax_min=MIN(r%nlevelmax,nlevelmax_file)
  allocate(dtold(1:nlevelmax_file))
  allocate(dtnew(1:nlevelmax_file))
  read(ilun)dtold(1:nlevelmax_file)
  read(ilun)dtnew(1:nlevelmax_file)
  if(nlevelmax_min>0)then
     g%dtold(1:nlevelmax_min)=dtold(1:nlevelmax_min)
     g%dtnew(1:nlevelmax_min)=dtnew(1:nlevelmax_min)
  endif
  deallocate(dtold,dtnew)
  read(ilun)g%nstep,g%nstep_coarse
  ! Read various constants
  read(ilun)g%const,g%mass_tot_0,g%rho_tot
  read(ilun)g%omega_m,g%omega_l,g%omega_k,g%omega_b,g%h0,g%aexp_ini,g%boxlen_ini
  g%omega_k=1-g%omega_m-g%omega_l
  read(ilun)g%aexp,g%texp,g%hexp
  read(ilun)g%aexp_old,g%epot_tot_int,g%epot_tot_old
  read(ilun)mass_sph_file
  read(ilun)gamma_file
  ! Read RNG seed to restart the stream
  read(ilun)r%seed
  close(ilun)
  ! For cosmo runs only, as mass_sph is not set in the namelist
  if(r%cosmo)r%mass_sph=mass_sph_file
  ! Check dimensions
  if(ndim.NE.ndim_file)then
     if(g%myid==1)then
        write(*,*)'Incorrect number of space dimensions in restart file'
     endif
     call mdl_abort(mdl)
  endif
  ! Compute movie frame number if applicable
  if(r%imovout>0) then
     do i=2,r%imovout
        if(r%aendmov>0)then
           if(g%aexp>r%aendmov*dble(i-1)/dble(r%imovout).and.g%aexp<r%aendmov*dble(i)/dble(r%imovout)) then
              r%imov=i
           endif
        else
           if(g%t>r%tendmov*dble(i-1)/dble(r%imovout).and.g%t<r%tendmov*dble(i)/dble(r%imovout)) then
              r%imov=i
           endif
        endif
     enddo
  endif

end subroutine input_params
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
recursive subroutine r_output_amr(pst,input_array,input_size,output_array,output_size)
  use mdl_module
  use amr_parameters, only: flen
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array

  character(LEN=flen)::filename
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_OUTPUT_AMR,pst%iUpper+1,input_size,output_size,input_array)
     call r_output_amr(pst%pLower,input_array,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size)
  else
     filename=transfer(input_array,filename)
     if(index(filename,'output')==0)then
        call backup_amr(pst%s%r,pst%s%g,pst%s%m,pst%s%mdl,filename)
     else
        call output_amr(pst%s,filename)
     endif
  endif

end subroutine r_output_amr
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine output_amr(s,filename)
  use amr_parameters, only: ndim,flen,twotondim
  use ramses_commons, only: ramses_t,open_file,close_file
  use mdl_module
  implicit none
  type(ramses_t)::s
  character(LEN=flen)::filename
  !-----------------------------------
  ! Output amr grid to file
  !-----------------------------------
  integer::ilun,ilevel,igrid,ind,refined_int
  integer(kind=8),dimension(s%r%levelmin:s%r%nlevelmax)::nskip

  associate(r=>s%r,g=>s%g,m=>s%m)

  call open_file(s,filename,nskip,ilun)

  do ilevel=r%levelmin,r%nlevelmax
     write(ilun,POS=nskip(ilevel))
     do igrid=m%head(ilevel),m%tail(ilevel)
        write(ilun)m%grid(igrid)%ckey
        refined_int=0
        do ind=1,twotondim
           if(m%grid(igrid)%refined(ind))refined_int=ibset(refined_int,ind-1)
        end do
        write(ilun)refined_int
     end do
  end do

  call close_file(s,filename,nskip,ilun)

  end associate

end subroutine output_amr
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine backup_amr(r,g,m,mdl,filename)
  use amr_parameters, only: ndim, flen, twotondim
  use amr_commons, only: run_t, global_t, mesh_t
  use mdl_module
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(mdl_t)::mdl
  character(LEN=flen)::filename
  !-----------------------------------
  ! Output amr grid to restart file
  !-----------------------------------
  integer::ilun,ierr,ilevel,igrid,ind,refined_int
  character(LEN=flen)::fileloc
  character(LEN=5)::nchar
  logical::file_exist

  ilun=10+mdl_core(mdl)
  call title(g%myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)
  inquire(file=fileloc, exist=file_exist)
  if (file_exist) then
     open(unit=ilun,file=fileloc,iostat=ierr)
     close(ilun,status="delete")
  end if
  open(unit=ilun,file=fileloc,access="stream",action="write",form='unformatted')
  write(ilun)ndim
  write(ilun)r%levelmin
  write(ilun)r%nlevelmax
  do ilevel=r%levelmin,r%nlevelmax
     write(ilun)m%noct(ilevel)
  end do
  do ilevel=r%levelmin,r%nlevelmax
     do igrid=m%head(ilevel),m%tail(ilevel)
        write(ilun)m%grid(igrid)%ckey
        refined_int=0
        do ind=1,twotondim
           if(m%grid(igrid)%refined(ind))refined_int=ibset(refined_int,ind-1)
        end do
        write(ilun)refined_int
     end do
  end do
  close(ilun)

end subroutine backup_amr
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine output_info(r,g,filename)
  use amr_parameters, only: ndim, flen
  use amr_commons, only: run_t, global_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  character(LEN=flen)::filename

  integer::ilun
  real(kind=8)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  character(LEN=flen)::fileloc

  if(r%verbose)write(*,*)'Entering output_info'

  ilun=11

  ! Conversion factor from user units to cgs units
  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Open file
  fileloc=TRIM(filename)
  open(unit=ilun,file=fileloc,form='formatted')

  ! Write run parameters
  write(ilun,'("nfile       =",I11)')r%nfile
  write(ilun,'("ncpu        =",I11)')g%ncpu
  write(ilun,'("ndim        =",I11)')ndim
  write(ilun,'("levelmin    =",I11)')r%levelmin
  write(ilun,'("levelmax    =",I11)')r%nlevelmax
  write(ilun,'("ngridmax    =",I11)')r%ngridmax
  write(ilun,'("nstep_coarse=",I11)')g%nstep_coarse
  write(ilun,*)

  ! Write physical parameters
  write(ilun,'("boxlen      =",E23.15)')r%boxlen
  write(ilun,'("time        =",E23.15)')g%t
  write(ilun,'("texp        =",E23.15)')g%texp
  write(ilun,'("aexp        =",E23.15)')g%aexp
  write(ilun,'("H0          =",E23.15)')g%h0
  write(ilun,'("omega_m     =",E23.15)')g%omega_m
  write(ilun,'("omega_l     =",E23.15)')g%omega_l
  write(ilun,'("omega_k     =",E23.15)')g%omega_k
  write(ilun,'("omega_b     =",E23.15)')g%omega_b
  write(ilun,'("gamma       =",E23.15)')r%gamma
  write(ilun,'("unit_l      =",E23.15)')scale_l
  write(ilun,'("unit_d      =",E23.15)')scale_d
  write(ilun,'("unit_t      =",E23.15)')scale_t
  write(ilun,*)

  close(ilun)

end subroutine output_info
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine output_header(r,g,p,filename)
  use amr_parameters, only: flen
  use amr_commons, only: run_t,global_t
  use pm_commons, only: part_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(part_t)::p
  character(LEN=flen)::filename

  ! Local variables
  integer::ilun
  character(LEN=flen)::fileloc

  if(r%verbose)write(*,*)'Entering output_header'

  ilun=10!+g%myid

  ! Open file
  fileloc=TRIM(filename)
  open(unit=ilun,file=fileloc,form='formatted')

  ! Write header information
  write(ilun,*)'Total number of particles'
  write(ilun,*)p%npart_tot

  write(ilun,*)'Total number of files'
  if(index(filename,'output')==0)then
     write(ilun,*)g%ncpu
  else
     write(ilun,*)r%nfile
  endif

  ! Keep track of what particle fields are present
  write(ilun,*)'Particle fields'
  write(ilun,'(a)',advance='no')'pos vel mass '
  if(allocated(p%phip))then
     write(ilun,'(a)',advance='no')'potential '
  endif
  if(allocated(p%zp))then
     write(ilun,'(a)',advance='no')'metallicity '
  endif
  if(allocated(p%fp))then
     write(ilun,'(a)',advance='no')'accel '
  endif
  if(allocated(p%jp))then
     write(ilun,'(a)',advance='no')'angmom '
  end if
  if(allocated(p%tp))then
     write(ilun,'(a)',advance='no')'birth_date '
  endif
  if(allocated(p%tm))then
     write(ilun,'(a)',advance='no')'merging_date '
  endif
  write(ilun,'(a)',advance='no')'level birth_id '
  if(allocated(p%idm))then
     write(ilun,'(a)',advance='no')'merging_id '
  endif
  if(allocated(p%idt))then
     write(ilun,'(a)',advance='no')'tracking_id '
  endif
  if(allocated(p%size))then
     write(ilun,'(a)',advance='no')'size '
  endif
  if(allocated(p%charge))then
     write(ilun,'(a)',advance='no')'charge '
  endif
  close(ilun)

end subroutine output_header
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine input_header(r,g,filename,npart_tot_file,ncpu_file)
  use amr_parameters, only: i8b,flen
  use amr_commons, only: run_t,global_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  character(LEN=flen)::filename
  integer(kind=8)::npart_tot_file
  integer::ncpu_file

  integer::ilun
  character(LEN=flen)::fileloc

  if(r%verbose)write(*,*)'Entering input_header'

  ilun=10!+g%myid

  ! Write header information
  fileloc=TRIM(filename)
  open(unit=ilun,file=fileloc,form='formatted')
  read(ilun,*)
  read(ilun,*)npart_tot_file
  read(ilun,*)
  read(ilun,*)ncpu_file
  close(ilun)

end subroutine input_header
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
#ifdef GADGET
subroutine savegadget(filename)
  use amr_parameters, only: flen
  use amr_commons
  use hydro_commons
  use pm_commons
  use gadgetreadfilemod
  implicit none
  character(LEN=flen)::filename
  TYPE (gadgetheadertype) :: header
  real,allocatable,dimension(:,:)::pos, vel
  integer(i8b),allocatable,dimension(:)::ids
  integer::i, idim, ipart
  real:: gadgetvfact
  integer::info
  integer(kind=8)::npart_tot
  integer(kind=4)::npart_loc
  real, parameter:: RHOcrit = 2.7755d11

  npart_tot=int(npart,kind=8)

  allocate(pos(ndim, npart), vel(ndim, npart), ids(npart))
  gadgetvfact = 100.0 * boxlen_ini / aexp / SQRT(aexp)

  header%npart = 0
  header%npart(2) = npart
  header%mass = 0
  header%mass(2) = omega_m*RHOcrit*(boxlen_ini)**3/npart_tot/1.d10
  header%time = aexp
  header%redshift = 1.d0/aexp-1.d0
  header%flag_sfr = 0
  header%nparttotal = 0
  header%nparttotal(2) = MOD(npart_tot,4294967296)
  header%flag_cooling = 0
  header%numfiles = ncpu
  header%boxsize = boxlen_ini
  header%omega0 = omega_m
  header%omegalambda = omega_l
  header%hubbleparam = h0/100.0
  header%flag_stellarage = 0
  header%flag_metals = 0
  header%totalhighword = 0
  header%totalhighword(2) = npart_tot/4294967296
  header%flag_entropy_instead_u = 0
  header%flag_doubleprecision = 0
  header%flag_ic_info = 0
  header%lpt_scalingfactor = 0
  header%unused = ' '

  do idim=1,ndim
     ipart=0
     do i=1,npartmax
        if(levelp(i)>0)then
           ipart=ipart+1
           if (ipart .gt. npart) then
                write(*,*) myid, "Ipart=",ipart, "exceeds", npart
                stop
           endif
           pos(idim, ipart)=xp(i,idim) * boxlen_ini
           vel(idim, ipart)=vp(i,idim) * gadgetvfact
           if (idim.eq.1) ids(ipart) = idp(i)
        end if
     end do
  end do

  call gadgetwritefile(filename, myid-1, header, pos, vel, ids)
  deallocate(pos, vel, ids)

end subroutine savegadget
#endif
end module output_amr_module
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
