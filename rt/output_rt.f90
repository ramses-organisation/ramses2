module output_rt_module
contains
!###################################################
!###################################################
!###################################################
!###################################################
recursive subroutine r_output_rt(pst,input_array,input_size,output_array,output_size)
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
     rID = mdl_send_request(pst%s%mdl,MDL_OUTPUT_RT,pst%iUpper+1,input_size,output_size,input_array)
     call r_output_rt(pst%pLower,input_array,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size)
  else
     filename=transfer(input_array,filename)
     if(index(filename,'output')==0)then
        call backup_rt(pst%s%r,pst%s%g,pst%s%m,pst%s%mdl,filename)
     else
        call output_rt(pst%s,filename)
     endif
  endif

end subroutine r_output_rt
!###################################################
!###################################################
!###################################################
!###################################################
subroutine output_rt(s,filename)
  use amr_parameters, only: ndim,twotondim, flen, nbuf_io
  use rt_parameters, only: nrtvar, nrtgrp
  use ramses_commons, only: ramses_t, open_file, close_file
  use posix_io, only: pio_t, pio_open, pio_write, pio_close
  use iso_c_binding, only: c_loc
  implicit none
  type(ramses_t)::s
  character(LEN=flen)::filename
  !-----------------------------------
  ! Output rt data in file
  !-----------------------------------
  integer::ilevel,igrid,ilun,igrp,idim,ind,i,nchunk
  integer(kind=8),dimension(s%r%levelmin:s%r%nlevelmax)::nskip
  integer(kind=8)::offset,nbyte
  real(kind=4),allocatable,dimension(:,:,:),target::qbuf
  type(pio_t)::f
  real(kind=8),dimension(1:twotondim,1:nrtvar)::qold
  real(kind=8),dimension(1:twotondim,1:nrtvar)::rtuold
  logical::overflow_reported=.false.

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)
#ifdef DO_RT
  ! Buffer octs to write in large chunks
  nchunk=max(1,nbuf_io/(4*twotondim*nrtvar))
  allocate(qbuf(1:twotondim,1:nrtvar,1:nchunk))
  call open_file(s,filename,nskip,ilun)
  call pio_open(f,ilun)
  overflow_reported=.false.
  do ilevel=r%levelmin,r%nlevelmax
     offset=nskip(ilevel)-1
     i=0
     do igrid=m%head(ilevel),m%tail(ilevel)
        rtuold=m%rtuold(:,:,igrid)
        ! The photon density is floored; flux components keep their sign
        do ind=1,twotondim
            do igrp = 1, nrtgrp
              qold(ind,1+(igrp-1)*(ndim+1)) = max(rtuold(ind,1+(igrp-1)*(ndim+1))*g%rt_c(ilevel),tiny(0d0))
              do idim = 1, ndim
                qold(ind,1+idim+(igrp-1)*(ndim+1)) = rtuold(ind,1+idim+(igrp-1)*(ndim+1))
              end do
            end do
        end do
        i=i+1
        qbuf(:,:,i)=real(qold,kind=4)
        if(maxval(abs(dble(qbuf(:,:,i)))).gt.1d31 .and. s%g%myid==1 .and. .not. overflow_reported) then
            print*,'The RT variables have very high values and are overflowing in the outputs'
            overflow_reported=.true.
        endif
        if(i==nchunk.or.igrid==m%tail(ilevel))then
           nbyte=int(4*twotondim*nrtvar,kind=8)*i
           call pio_write(f,c_loc(qbuf),nbyte,offset)
           offset=offset+nbyte
           i=0
        endif
     end do
  end do
  call pio_close(f)
  call close_file(s,filename,nskip,ilun)
  deallocate(qbuf)
#endif
  end associate

end subroutine output_rt
!###################################################
!###################################################
!###################################################
!###################################################
subroutine backup_rt(r,g,m,mdl,filename)
  use amr_parameters, only: ndim,flen,dp,twotondim,nbuf_io
  use rt_parameters, only: nrtvar
  use amr_commons, only: run_t,global_t,mesh_t
  use posix_io, only: pio_t,pio_open,pio_write,pio_close
  use iso_c_binding, only: c_loc
  use mdl_module
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(mdl_t)::mdl
  character(LEN=flen)::filename

  integer::ilevel,igrid,ilun,ierr,nchunk,i
  character(LEN=5)::nchar
  character(LEN=flen)::fileloc
  logical::file_exist
  integer(kind=8)::offset,nbyte
  real(dp),allocatable,dimension(:,:,:),target::dbuf
  type(pio_t)::f
  ! Buffer octs to write in large chunks
  nchunk=max(1,nbuf_io/((storage_size(1.0_dp)/8)*twotondim*nrtvar))
  allocate(dbuf(1:twotondim,1:nrtvar,1:nchunk))
  ilun=10+mdl_core(mdl)
  call title(g%myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)
  inquire(file=fileloc,exist=file_exist)
  if (file_exist) then
     open(unit=ilun,file=fileloc,iostat=ierr)
     close(ilun,status="delete")
  end if
  open(unit=ilun,file=fileloc,access="stream",action="write",form='unformatted')
  write(ilun)ndim
  write(ilun)nrtvar
  write(ilun)r%levelmin
  write(ilun)r%nlevelmax
  do ilevel=r%levelmin,r%nlevelmax
     write(ilun)m%noct(ilevel)
  enddo
  ! Bulk data starts right after the header
  inquire(unit=ilun,pos=offset)
  offset=offset-1
  call pio_open(f,ilun)
  do ilevel=r%levelmin,r%nlevelmax
#ifdef DO_RT
     i=0
     do igrid=m%head(ilevel),m%tail(ilevel)
        i=i+1
        dbuf(:,:,i)=m%rtuold(:,:,igrid)
        if(i==nchunk.or.igrid==m%tail(ilevel))then
           nbyte=int((storage_size(1.0_dp)/8)*twotondim*nrtvar,kind=8)*i
           call pio_write(f,c_loc(dbuf),nbyte,offset)
           offset=offset+nbyte
           i=0
        endif
     end do
#endif
  enddo
  call pio_close(f)
  close(ilun)
  deallocate(dbuf)
end subroutine backup_rt
!########################################################################
!########################################################################
!########################################################################
!########################################################################
subroutine file_descriptor_rt(r,filename,write_bkp_file)
  use amr_parameters, only: ndim,flen
  use rt_parameters, only: nrtgrp
  use amr_commons, only: run_t
  implicit none
  type(run_t)::r
  character(LEN=flen)::filename
  logical::write_bkp_file
  character(len=1), dimension(1:3), parameter :: dim_keys = ["x", "y", "z"]

  character(LEN=flen)::fileloc
  integer::igrp,idim,ilun

  if(r%verbose)write(*,*)'Entering file_descriptor_rt'

  ilun=11

  ! Open file
  fileloc=TRIM(filename)
  open(unit=ilun,file=fileloc,form='formatted')

  ! Write variable names in backup file
  write(ilun,'("nvar        =",I11)')nrtgrp*(1+ndim)
  do igrp = 1, nrtgrp
     write(ilun,'("variable #",I2,": photon_flux_", i0.2)')1+(igrp-1)*(ndim+1), igrp
     do idim = 1, ndim
        write(ilun,'("variable #",I2,": photon_flux_", i0.2, "_", a)')1+idim+(igrp-1)*(ndim+1), igrp, dim_keys(idim)
     end do
  end do
  close(ilun)

end subroutine file_descriptor_rt

!########################################################################
!########################################################################
!########################################################################
!########################################################################
subroutine output_rtinfo(r, g, filename)

! Output rt information into info_rt_XXXXX.txt
!------------------------------------------------------------------------
  use amr_parameters, only: flen
  use amr_commons, only: run_t, global_t
  use hydro_parameters, only: nion
  use rt_parameters, only: nrtvar, nrtgrp
  implicit none
  type(run_t)::r
  type(global_t)::g
  character(LEN=flen)::filename, fileloc
  integer :: ilun
  real(kind=8)::scale_np,scale_fp
!------------------------------------------------------------------------
  if (r%verbose) write(*,*)'Entering output_rtinfo'

  ilun=11

#ifdef DO_RT
  ! Conversion factor from user units to cgs units
  call rt_units(r,g,scale_np,scale_fp)
#endif

  ! Open file
  fileloc=TRIM(filename)
  open(unit=ilun,file=fileloc,form='formatted')

  ! Write run parameters
  write(ilun,'("nrtvar       = ", I11)') nrtvar
  write(ilun,'("nrtgrp       = ", I11)') nrtgrp
  write(ilun,'("nion         = ", I11)') nion
  write(ilun,'("iion         = ", I11)') r%iions

  ! Write cooling parameters
  write(ilun,'("x_h          = ", E23.15)') r%x_h
  write(ilun,'("y_he         = ", E23.15)') r%y_he

  ! Write physical parameters
  write(ilun,'("unit_np      = ", E23.15)') scale_np
  write(ilun,'("unit_fp      = ", E23.15)') scale_fp
  write(ilun,'("rt_c_fraction= ", 100(E15.7))') r%rt_c_fraction(r%levelmin:r%nlevelmax)
  write(ilun,*)

  ! Write photon group properties
  call write_group_props(r, .false., ilun)

  close(ilun)

end subroutine output_rtInfo

!########################################################################
!########################################################################
!########################################################################
!########################################################################
subroutine output_rt_groups(r, filename)

! Output the photon group properties, for a restart
!------------------------------------------------------------------------
  use amr_parameters, only: flen
  use amr_commons, only: run_t
  use rt_parameters, only: nrtgrp
  implicit none
  type(run_t)::r
  character(LEN=flen)::filename
  integer :: ilun
!------------------------------------------------------------------------
  ilun=11
  open(unit=ilun,file=TRIM(filename),access="stream",action="write",form='unformatted')
  write(ilun)nrtgrp,size(r%group_csn)
  write(ilun)r%group_egy,r%group_csn,r%group_cse
  close(ilun)

end subroutine output_rt_groups

!########################################################################
!########################################################################
!########################################################################
!########################################################################
subroutine input_rt_groups(r, g, filename)

! Read back the photon group properties written by output_rt_groups
!------------------------------------------------------------------------
  use amr_parameters, only: flen
  use amr_commons, only: run_t, global_t
  use rt_parameters, only: nrtgrp
  implicit none
  type(run_t)::r
  type(global_t)::g
  character(LEN=flen)::filename
  integer :: ilun, nrtgrp_file, ncs_file
  logical :: file_exist
!------------------------------------------------------------------------
  inquire(file=filename, exist=file_exist)
  if(.not.file_exist)then
     if(g%myid==1)write(*,*)'Missing '//TRIM(filename)
     return
  endif
  if(g%myid==1)write(*,*)'Reading '//TRIM(filename)
  ilun=11
  open(unit=ilun,file=TRIM(filename),access="stream",action="read",form='unformatted')
  read(ilun)nrtgrp_file,ncs_file
  if(nrtgrp_file/=nrtgrp.or.ncs_file/=size(r%group_csn))then
     write(*,*)'Error: '//TRIM(filename)//' does not match the photon groups of this run'
     stop
  endif
  read(ilun)r%group_egy,r%group_csn,r%group_cse
  close(ilun)

end subroutine input_rt_groups

!########################################################################
!########################################################################
!########################################################################
!########################################################################
subroutine write_group_props(r, update, lun)

! Write photon group properties to file or std output.
! lun => File identifier (use 6 for std. output)
!------------------------------------------------------------------------
  use amr_commons, only: run_t
  use rt_parameters, only: nrtgrp
#ifdef DO_RTZ
  use rtz_module, only: elements, n_elements
#endif
  implicit none
  type(run_t)::r
  logical :: update
  integer :: ip, lun
  integer :: iE
!------------------------------------------------------------------------
  if (.not. update) then
     write(lun,*) 'Photon group properties=------------------------------ '
  else
     write(lun,*) 'Photon properties have been changed to--------------- '
  end if
  write(lun, 901) r%group_L0(:)
  write(lun, 902) r%group_L1(:)
#ifndef DO_RTZ
  write(lun, 903) r%spec2group(:)
#endif
  do ip = 1, nrtgrp
     write(lun, 907) ip
     write(lun, 904) r%group_egy(ip)
#ifdef DO_RTZ
     do iE=1,n_elements
        if (elements(iE)%atomic_number.gt.0) then
           write(lun, 905) elements(iE)%element_name, r%group_csn(ip,iE,1:elements(iE)%n_ions+elements(iE)%n_mol)
           write(lun, 905) elements(iE)%element_name, r%group_cse(ip,iE,1:elements(iE)%n_ions+elements(iE)%n_mol)
        end if
     end do
#else
     write(lun, 905) r%group_csn(ip,:)
     write(lun, 906) r%group_cse(ip,:)
#endif
  end do
  write (lun,*) '=-----------------------------------------------------'

901 format ('  groupL0  [eV]  = ', 20f12.3)
902 format ('  groupL1  [eV]  = ', 20f12.3)
903 format ('  spec2group     = ', 20I12)
904 format ('  egy      [eV]  = ', 1pe12.3)
#ifdef DO_RTZ
905 format ('  csn    [cm^2]  = ', A20, 27(1pe12.3))
906 format ('  cse    [cm^2]  = ', A20, 27(1pe12.3))
#else
905 format ('  csn    [cm^2]  = ', 20(1pe12.3))
906 format ('  cse    [cm^2]  = ', 20(1pe12.3))
#endif
907 format ('  --=Group', I2)

end subroutine write_group_props

!########################################################################
!########################################################################
!########################################################################
!########################################################################
subroutine output_photon_emission_stats(r,g)

! Output and reset  statistics on stellar emission of radiation.
! star rt feedback statistics
!------------------------------------------------------------------------
  use amr_commons, only: run_t, global_t
  use constants, only: yr2sec
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  type(run_t)::r
  type(global_t)::g
  real(kind=8) :: step_nPhot_all, step_nStar_all, step_mStar_all
  real(kind=8) :: scale_l, scale_t, scale_d, scale_v, scale_nh, scale_T2
#ifndef WITHOUTMPI
  integer :: ierr
#endif
!------------------------------------------------------------------------
  call units(r, g, scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2)

  step_nPhot_all = 0d0 ; step_nStar_all = 0d0 ; step_mStar_all = 0d0
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(g%step_nPhot,           step_nPhot_all,  1,      &
       MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(g%step_nStar,           step_nStar_all,  1,      &
          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(g%step_mStar,           step_mStar_all,  1,      &
          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  g%step_nPhot  = step_nPhot_all
  g%step_nStar  = step_nStar_all
  g%step_mStar  = step_mStar_all
#endif
  g%tot_nPhot = g%tot_nPhot + g%step_nPhot
  if (g%myid .eq. 1) then
      if(r%cosmo) then
          write(*, 113) g%step_nPhot, g%tot_nPhot                    &
                  , g%step_nStar/g%dtnew(r%levelmin)                 &
                  , g%step_mStar/g%dtnew(r%levelmin)                 &
                  , g%aexp, g%texp
      else
          write(*, 114) g%step_nPhot, g%tot_nPhot                    &
                  , g%step_nStar/g%dtnew(r%levelmin)                 &
                  , g%step_mStar/g%dtnew(r%levelmin)                 &
                  , g%t*scale_t/yr2sec
      endif
  endif
  g%step_nPhot = 0d0 ; g%step_nStar = 0d0 ; g%step_mStar = 0d0
113 format(' Stellar radiation(phot/step, phot/tot, *, */Msun, aexp, texp)= '  &
                                                             , 10(1pe11.3))
114 format(' Stellar radiation(phot/step, phot/tot, *, */Msun, t_yr)= '  &
                                                             , 10(1pe11.3))
end subroutine output_photon_emission_stats

end module output_rt_module
