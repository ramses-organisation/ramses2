! Core clump catalogue output.
module core_clump_output_module

contains
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine output_core_clump_catalogue(pst)
  use amr_parameters, only: flen
  use mdl_module, only: mdl_mkdir
  use ramses_commons, only: pst_t
  use output_clump_module, only: file_descriptor_clump, output_peak_header
  implicit none
  type(pst_t)::pst
  !------------------------------------------------------------------
  ! Write the core clump catalogue of the output just written.
  !------------------------------------------------------------------
  character(LEN=5)::nchar
  character(LEN=flen)::filename,filedir,fileprefix
  integer,dimension(1:flen/4)::input_array
  integer::dummy(1)

#ifdef GRAV
  associate(r=>pst%s%r,g=>pst%s%g,mdl=>pst%s%mdl,p=>pst%s%p,star=>pst%s%star, &
            sink=>pst%s%sink,tree=>pst%s%tree,c=>pst%s%c)

  if(.not.g%output_done)return
  if(.not.(r%output_core_clump.or.r%output_peak_grid.or.r%output_peak_part.or. &
       r%output_peak_star.or.r%output_peak_sink.or.r%output_peak_tree.or. &
       r%output_peak_core))return

  ! Directory of the output just written
  call title(g%ifout-1,nchar)
  filedir='output_'//TRIM(nchar)//'/'
  call mdl_mkdir(mdl,filedir)

  if(r%output_peak_grid)then
     filename=TRIM(filedir)//'core_peak_grid_header.txt'
     call file_descriptor_clump(r,filename)
  endif
  if(r%output_peak_part.and.r%part)then
     filename=TRIM(filedir)//'core_peak_part_header.txt'
     call output_peak_header(r,g,p,filename)
  endif
  if(r%output_peak_star.and.r%star)then
     filename=TRIM(filedir)//'core_peak_star_header.txt'
     call output_peak_header(r,g,star,filename)
  endif
  if(r%output_peak_sink.and.r%sink)then
     filename=TRIM(filedir)//'core_peak_sink_header.txt'
     call output_peak_header(r,g,sink,filename)
  endif
  if(r%output_peak_tree.and.r%tree)then
     filename=TRIM(filedir)//'core_peak_tree_header.txt'
     call output_peak_header(r,g,tree,filename)
  endif
  if(r%output_peak_core.and.r%core)then
     filename=TRIM(filedir)//'core_peak_core_header.txt'
     call output_peak_header(r,g,pst%s%core,filename)
  endif

  fileprefix=' '
  fileprefix=TRIM(filedir)//'core_'
  input_array=transfer(fileprefix,input_array)
  call r_output_core_clump(pst,input_array,flen/4,dummy,0)

  end associate
#endif

end subroutine output_core_clump_catalogue
!###########################################################
!###########################################################
!###########################################################
!###########################################################
recursive subroutine r_output_core_clump(pst,input_array,input_size,output_array,output_size)
  use mdl_module
  use amr_parameters, only: flen
  use ramses_commons, only: pst_t
  use mdl_parameters
  use output_clump_module, only: output_clump_field, output_peak_part
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array

  character(LEN=flen)::filename,fileloc
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_OUTPUT_CORE_CLUMP,pst%iUpper+1,input_size,output_size,input_array)
     call r_output_core_clump(pst%pLower,input_array,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size)
  else
     filename=transfer(input_array,filename)
     ! Write core clump properties
     if(pst%s%r%output_core_clump)then
        call output_core_clump_properties(pst%s,filename)
     endif
     ! Write density, peak id and halo id for AMR cells and particles
     if(pst%s%r%output_peak_grid)then
        fileloc=TRIM(filename)//'peak_grid.'
        call output_clump_field(pst%s,fileloc)
     endif
     if(pst%s%r%output_peak_part.and.pst%s%r%part)then
        fileloc=TRIM(filename)//'peak_part.'
        call output_peak_part(pst%s,pst%s%p,fileloc)
     endif
     if(pst%s%r%output_peak_star.and.pst%s%r%star)then
        fileloc=TRIM(filename)//'peak_star.'
        call output_peak_part(pst%s,pst%s%star,fileloc)
     endif
     if(pst%s%r%output_peak_sink.and.pst%s%r%sink)then
        fileloc=TRIM(filename)//'peak_sink.'
        call output_peak_part(pst%s,pst%s%sink,fileloc)
     endif
     if(pst%s%r%output_peak_tree.and.pst%s%r%tree)then
        fileloc=TRIM(filename)//'peak_tree.'
        call output_peak_part(pst%s,pst%s%tree,fileloc)
     endif
     if(pst%s%r%output_peak_core.and.pst%s%r%core)then
        fileloc=TRIM(filename)//'peak_core.'
        call output_peak_part(pst%s,pst%s%core,fileloc)
     endif
  endif

end subroutine r_output_core_clump
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine output_core_clump_properties(s,filename)
  use amr_parameters, only: flen
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t)::s
  character(LEN=flen)::filename
  !------------------------------------------------------------------
  ! Properties of the relevant core clumps, one line per peak.
  !------------------------------------------------------------------
  character(LEN=flen)::fileloc
  integer::ilun,j,hid,gid,bound_flag,contracting_flag
  real(kind=8)::mass,r200b,rmax,concentration,rho_ave

  associate(c=>s%c)

  fileloc=TRIM(filename)//'clump.'
  call open_core_clump_file(s,fileloc,ilun)

  do j=1,c%npeak
     ! Global peak id
     gid=int(j+c%npeak_cum(s%g%myid-1))
     ! Without halo merging, each clump is its own halo
     if(c%saddle_threshold>0)then
        hid=c%ind_halo(j)
     else
        hid=c%new_peak(j)
     endif
     if (    c%relevance(j) > c%relevance_threshold.AND. &
          & c%clump_mass(j) > c%mass_threshold.AND. &
          & c%n_cells(j) > 0)then

        mass=c%clump_mass(j)
        r200b=0d0; rmax=0d0; concentration=0d0
        rho_ave=c%clump_mass(j)/c%clump_vol(j)

        bound_flag=0
        if(c%bound(j))bound_flag=1
        contracting_flag=0
        if(c%contracting(j))contracting_flag=1

        write(ilun,'(I10,3(1X,I10),16(1X,1PE18.9E2),2(1X,I3))')&
              gid&
             ,c%new_peak(j)&
             ,hid&
             ,c%n_cells(j)&
             ,c%min_dens(j)&
             ,c%max_dens(j)&
             ,c%tidal_dens(j)&
             ,rho_ave&
             ,c%clump_mass(j)&
             ,mass&
             ,r200b&
             ,rmax&
             ,concentration&
             ,c%peak_pos(j,1),c%peak_pos(j,2),c%peak_pos(j,3)&
             ,c%peak_vel(j,1),c%peak_vel(j,2),c%peak_vel(j,3)&
             ,c%virial_ratio(j)&
             ,bound_flag&
             ,contracting_flag
     endif
  end do

  call close_core_clump_file(s,ilun)

  end associate

end subroutine output_core_clump_properties
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine open_core_clump_file(s,filename,ilun)
  use amr_parameters, only: flen
  use ramses_commons, only: ramses_t
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  type(ramses_t)::s
  character(LEN=flen),intent(in)::filename
  integer,intent(out)::ilun
  !------------------------------------------------------------------
  ! Open the core clump file of this rank; the first rank of each file
  ! writes the column names.
  !------------------------------------------------------------------
  character(LEN=flen)::fileloc
  character(LEN=5)::nchar
  integer,dimension(1:s%r%nfile+1)::istart
  integer::i,ifile,ncpufile,nremain,ierr,dummy,tag
  logical::file_exist
#ifndef WITHOUTMPI
  integer::info,reqrecv
  integer,dimension(MPI_STATUS_SIZE)::reqstatus
#endif

  associate(r=>s%r,g=>s%g)

  tag=901
  dummy=1
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

  ilun=10+istart(ifile)
  call title(ifile,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)

  if(g%myid.EQ.istart(ifile))then
     inquire(file=fileloc, exist=file_exist)
     if (file_exist) then
        open(unit=ilun,file=fileloc,iostat=ierr)
        close(ilun,status='delete')
     end if
     open(unit=ilun,file=fileloc,form='formatted')
     write(ilun,'(A)')'     index     parent       halo      ncell'//&
          '    rho_min            rho_max            rho_sad            rho_ave    '//&
          '        mpatch             mass               r200               rmax               c200   '//&
          '            pos_x              pos_y              pos_z      '//&
          '        vel_x              vel_y              vel_z      '//&
          ' virial_ratio bound contracting'
  elseif(g%myid.GT.istart(ifile))then
#ifndef WITHOUTMPI
     call MPI_IRECV(dummy,1,MPI_INTEGER,g%myid-2,tag,MPI_COMM_WORLD,reqrecv,info)
     call MPI_WAIT(reqrecv,reqstatus,info)
#endif
     inquire(file=fileloc, exist=file_exist)
     if (.not. file_exist) then
        write(*,*)'Error: core clump file not found in open_core_clump_file'
        stop
     end if
     open(unit=ilun,file=fileloc,form='formatted',access='append')
  endif

  end associate

end subroutine open_core_clump_file
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine close_core_clump_file(s,ilun)
  use ramses_commons, only: ramses_t
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  type(ramses_t)::s
  integer,intent(in)::ilun
  !------------------------------------------------------------------
  ! Close the core clump file and let the next rank of the file append.
  !------------------------------------------------------------------
  integer,dimension(1:s%r%nfile+1)::istart
  integer::i,ifile,ncpufile,nremain,dummy,tag
#ifndef WITHOUTMPI
  integer::info,reqsend
  integer,dimension(MPI_STATUS_SIZE)::reqstatus
#endif

  associate(r=>s%r,g=>s%g)

  tag=901
  dummy=1
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
#ifndef WITHOUTMPI
     call MPI_ISEND(dummy,1,MPI_INTEGER,g%myid,tag,MPI_COMM_WORLD,reqsend,info)
     call MPI_WAIT(reqsend,reqstatus,info)
#endif
  endif

  end associate

end subroutine close_core_clump_file

end module core_clump_output_module
