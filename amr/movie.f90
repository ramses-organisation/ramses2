module movie_module
contains
!=======================================================================
!=======================================================================
!=======================================================================
!=======================================================================
subroutine m_output_frame(pst)
  use amr_parameters, only: ndim, nvector, twotondim, flen
  use hydro_parameters, only: nvar
  use rt_parameters, only: nrtgrp
  use ramses_commons, only: pst_t
  use output_amr_module, only: output_info
  use mdl_module, only: mdl_mkdir
  implicit none
  type(pst_t)::pst

  ! Master driver for movie frame output. Runs on rank 1 only (by call
  ! structure: amr_step -> m_output_frame). For each projection it asks
  ! the MDL service r_output_frame to compute a per-CPU map for each
  ! enabled variable and writes the reduced result to disk.

  character(len=1)::temp_string
  character(len=5)::istep_str,dummy
  character(len=flen)::moviedir,infofile
  character(len=flen),dimension(0:nvar+2+nrtgrp)::moviefiles
  integer::ilun,kk,ll,ind_proj,input_size,output_size
  integer,dimension(:),allocatable::input_array,output_array
  real(kind=8)::delx,dely,delz,timer
  real(kind=8),dimension(:),allocatable::data_frame,dens
  real(kind=4),dimension(:),allocatable::data_single
  logical::is_mean,is_sum,is_min,is_max

  associate(r=>pst%s%r,g=>pst%s%g,m=>pst%s%m,mdl=>pst%s%mdl)

  do ind_proj=1,LEN(trim(r%proj_axis))

#if NDIM > 1
     if(r%imov<1)r%imov=1
     if(r%imov>r%imovout)return

     is_mean = (r%method_frame(ind_proj)(1:4).eq.'mean')
     is_sum  = (trim(r%method_frame(ind_proj)).eq.'sum')
     is_min  = (trim(r%method_frame(ind_proj)).eq.'min')
     is_max  = (trim(r%method_frame(ind_proj)).eq.'max')

     write(*,*)'Computing and dumping movie frame'

     call title(r%imov, istep_str)
     write(temp_string,'(I1)') ind_proj
     moviedir = 'movie'//trim(temp_string)//'/'
     write(*,*) "Writing frame ", istep_str
     if(.not.g%withoutmkdir) call mdl_mkdir(mdl,moviedir)

     infofile = trim(moviedir)//'info_'//trim(istep_str)//'.txt'
     call output_info(r,g,infofile)

     moviefiles(0) = trim(moviedir)//'colden_'//trim(istep_str)//'.map'
     moviefiles(1) = trim(moviedir)//'dens_'//trim(istep_str)//'.map'
     moviefiles(2) = trim(moviedir)//'vx_'//trim(istep_str)//'.map'
     moviefiles(3) = trim(moviedir)//'vy_'//trim(istep_str)//'.map'
     moviefiles(4) = trim(moviedir)//'vz_'//trim(istep_str)//'.map'
     moviefiles(5) = trim(moviedir)//'temp_'//trim(istep_str)//'.map'
#if NVAR>5
     do ll=6,NVAR
        write(dummy,'(I3.1)') ll
        moviefiles(ll) = trim(moviedir)//'var'//trim(adjustl(dummy))//'_'//trim(istep_str)//'.map'
     end do
#endif
     do ll=nvar+1,nvar+nrtgrp
        write(dummy,'(I3.1)') ll-nvar
        moviefiles(ll) = trim(moviedir)//'Fp'//trim(adjustl(dummy))//'_'//trim(istep_str)//'.map'
     end do
     moviefiles(NVAR+nrtgrp+1) = trim(moviedir)//'dm_'//trim(istep_str)//'.map'
     moviefiles(NVAR+nrtgrp+2) = trim(moviedir)//'stars_'//trim(istep_str)//'.map'

     allocate(data_single(1:r%nw_frame*r%nh_frame))
     allocate(data_frame (1:r%nw_frame*r%nh_frame))
     allocate(dens       (1:r%nw_frame*r%nh_frame))
     input_size=2
     output_size=2*r%nw_frame*r%nh_frame
     allocate(input_array(1:input_size))
     allocate(output_array(1:output_size))
     input_array(1)=ind_proj

     ! Frame extent on screen (for writing in the file header). The same
     ! axis-swap logic lives inside output_frame; keep it consistent here.
     if(r%cosmo) then
        timer = g%aexp
     else
        timer = g%t
     endif
     if(r%proj_axis(ind_proj:ind_proj).eq.'x')then
        delx=r%deltay_frame(ind_proj*2-1)+r%deltay_frame(ind_proj*2)/g%aexp
        dely=r%deltaz_frame(ind_proj*2-1)+r%deltaz_frame(ind_proj*2)/g%aexp
        delz=r%deltax_frame(ind_proj*2-1)+r%deltax_frame(ind_proj*2)/g%aexp
     elseif(r%proj_axis(ind_proj:ind_proj).eq.'y')then
        delx=r%deltax_frame(ind_proj*2-1)+r%deltax_frame(ind_proj*2)/g%aexp
        dely=r%deltaz_frame(ind_proj*2-1)+r%deltaz_frame(ind_proj*2)/g%aexp
        delz=r%deltay_frame(ind_proj*2-1)+r%deltay_frame(ind_proj*2)/g%aexp
     else
        delx=r%deltax_frame(ind_proj*2-1)+r%deltax_frame(ind_proj*2)/g%aexp
        dely=r%deltay_frame(ind_proj*2-1)+r%deltay_frame(ind_proj*2)/g%aexp
        delz=r%deltaz_frame(ind_proj*2-1)+r%deltaz_frame(ind_proj*2)/g%aexp
     endif

     if(r%hydro) then

        ! Weight map (used as denominator for mean_* methods). For sum/min/max
        ! we skip this; for mean_* the leaf decides what to put in the map
        ! based on r%method_frame(ind_proj).
        if(is_mean)then
           input_array(2)=0
           call r_output_frame(pst,input_array,input_size,output_array,output_size)
           dens=transfer(output_array,dens)
        else
           dens=1d0
        endif

        do kk=1,NVAR+nrtgrp
           if(r%movie_vars(kk).eq.1)then
              input_array(2)=kk
              call r_output_frame(pst,input_array,input_size,output_array,output_size)
              data_frame=transfer(output_array,data_frame)
              if(is_mean)then
                 where(dens>0d0)
                    data_frame=data_frame/dens
                 elsewhere
                    data_frame=0d0
                 end where
              else if(is_min)then
                 where(data_frame.ge.1d-3*huge(0.0))data_frame=0d0
              else if(is_max)then
                 where(data_frame.le.-1d-3*huge(0.0))data_frame=0d0
              endif
              call write_map(moviefiles(kk),data_frame,data_single,r%nw_frame,r%nh_frame,timer,delx,dely,delz)
           end if
        end do

        ! Optionally write the weight map itself (useful for diagnostics).
        if(is_mean .and. r%movie_vars(0).eq.1)then
           call write_map(moviefiles(0),dens,data_single,r%nw_frame,r%nh_frame,timer,delx,dely,delz)
        endif

     endif

     ! Dark matter projection - runs even without hydro (DMO).
     if(r%movie_vars(NVAR+nrtgrp+1).eq.1) then
        input_array(2) = NVAR+nrtgrp+1
        call r_output_frame(pst,input_array,input_size,output_array,output_size)
        data_frame = transfer(output_array,data_frame)
        call write_map(moviefiles(NVAR+nrtgrp+1),data_frame,data_single,r%nw_frame,r%nh_frame,timer,delx,dely,delz)
     end if

     ! Star projection.
     if(r%movie_vars(NVAR+nrtgrp+2).eq.1) then
        input_array(2) = NVAR+nrtgrp+2
        call r_output_frame(pst,input_array,input_size,output_array,output_size)
        data_frame = transfer(output_array,data_frame)
        call write_map(moviefiles(NVAR+nrtgrp+2),data_frame,data_single,r%nw_frame,r%nh_frame,timer,delx,dely,delz)
     end if

     deallocate(data_single)
     deallocate(data_frame)
     deallocate(dens)
     deallocate(input_array)
     deallocate(output_array)

     ! Advance counter only after the last projection.
     if(ind_proj.eq.len(trim(r%proj_axis))) then
        r%imov=r%imov+1
        if(r%aendmov>0)then
           do while(((r%aendmov-r%astartmov)*dble(r%imov)/dble(r%imovout)+r%astartmov<g%aexp).and.(r%imov.lt.r%imovout))
              r%imov=r%imov+1
           end do
        end if
        if(r%tendmov>0)then
           do while(((r%tendmov-r%tstartmov)*dble(r%imov)/dble(r%imovout)+r%tstartmov<g%t).and.(r%imov.lt.r%imovout))
              r%imov=r%imov+1
           end do
        end if
     endif

#endif

  enddo

  end associate

end subroutine m_output_frame
!=======================================================================
!=======================================================================
!=======================================================================
!=======================================================================
subroutine write_map(filename,frame,single,nw,nh,timer,delx,dely,delz)
  use amr_parameters, only: flen
  implicit none
  character(len=flen),intent(in)::filename
  real(kind=8),intent(in),dimension(:)::frame
  real(kind=4),intent(inout),dimension(:)::single
  integer,intent(in)::nw,nh
  real(kind=8),intent(in)::timer,delx,dely,delz
  integer::ilun
  ilun=10
  open(ilun,file=trim(filename),form='unformatted')
  single=real(frame,kind=4)
  rewind(ilun)
  write(ilun)timer,delx,dely,delz
  write(ilun)nw,nh
  write(ilun)single
  close(ilun)
end subroutine write_map
!=======================================================================
!=======================================================================
!=======================================================================
!=======================================================================
recursive subroutine r_output_frame(pst,input_array,input_size,output_array,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  use hilbert
  use hydro_parameters, only: nvar
  use rt_parameters, only: nrtgrp
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array

  integer,dimension(:),allocatable::next_output_array
  integer::ind_proj,ind_var
  real(kind=8),dimension(:),allocatable::map,next_map
  integer::rID
  character(len=10)::method

  ind_proj=input_array(1)
  ind_var =input_array(2)
  method  = pst%s%r%method_frame(ind_proj)

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_OUTPUT_FRAME,pst%iUpper+1,input_size,output_size,input_array)
     call r_output_frame(pst%pLower,input_array,input_size,output_array,output_size)
     allocate(next_output_array(1:output_size))
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_output_array)
     allocate(map(1:pst%s%r%nw_frame*pst%s%r%nh_frame))
     allocate(next_map(1:pst%s%r%nw_frame*pst%s%r%nh_frame))
     map=transfer(output_array,map)
     next_map=transfer(next_output_array,next_map)
     ! Hydro variables (ind_var in 0..nvar+nrtgrp) follow the method's reduction;
     ! particle (DM/star) maps are always summed.
     if(ind_var.ge.0 .and. ind_var.le.nvar+nrtgrp)then
        select case(trim(method))
        case('min')
           map=min(map,next_map)
        case('max')
           map=max(map,next_map)
        case default
           map=map+next_map
        end select
     else
        map=map+next_map
     endif
     output_array=transfer(map,output_array)
     deallocate(next_output_array)
     deallocate(map)
     deallocate(next_map)
  else
     allocate(map(1:pst%s%r%nw_frame*pst%s%r%nh_frame))
     ! Initialise the local map according to the reduction operator.
     select case(trim(method))
     case('min')
        map=1d-3*huge(0d0)
     case('max')
        map=-1d-3*huge(0d0)
     case default
        map=0d0
     end select
     call output_frame(pst%s,ind_proj,ind_var,pst%s%r%nw_frame*pst%s%r%nh_frame,map)
     output_array=transfer(map,output_array)
     deallocate(map)
  endif

end subroutine r_output_frame
!=======================================================================
!=======================================================================
!=======================================================================
!=======================================================================
subroutine output_frame(s,ind_proj,ind_var,map_size,map)
  use amr_parameters, only:ndim, nvector, twotondim
  use hydro_parameters, only: nvar, nener, ie
  use rt_parameters, only: nrtgrp
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t)::s
  integer::ind_proj,map_size,ind_var
  real(kind=8),dimension(1:map_size)::map
  real(kind=8),parameter::pi=3.14159265358979323846d0

  ! Local variables
  integer::nlevelmax_frame,nstride
  integer::ind,ind_map
  integer::imin,imax,jmin,jmax,ii,jj
  real(kind=8)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(kind=8)::xcen,ycen,zcen
  real(kind=8)::xleft_frame,xright_frame,yleft_frame,yright_frame,zleft_frame,zright_frame
  real(kind=8)::xleft,xright,yleft,yright,zleft,zright
  real(kind=8)::xxleft,xxright,yyleft,yyright,xxcen,yycen
  real(kind=8)::delx,dely,delz
  real(kind=8)::dx_frame,dy_frame,dx
  real(kind=8)::dx_cell,dy_cell,dz_cell,dvol,weight,uvar,e
  logical::ok
  real(kind=8),dimension(1:ndim)::xx
  real(kind=8)::temp,ekk
  integer::igrid,idim,ilevel,ip
  integer::box_xidx,box_yidx,box_zidx
  real(kind=8)::theta_cam,phi_cam,dist_cam,fov_camera
  real(kind=8)::alpha,beta,pers_corr,timer,timer_end
  real(kind=8)::base_dist,tend_theta_eff,tend_phi_eff,span
  real(kind=8)::xtmp,ytmp,ztmp
  real(kind=8)::xcentre,ycentre,zcentre,dx_proj
  real(kind=8)::xpc,ypc,xmap,ymap,zmap
  logical::is_sphere,is_square,is_cube
  logical::is_mean,is_mean_mass,is_mean_dens,is_mean_vol,is_sum,is_min,is_max
  logical::perspective
  character(len=6)::shader
  character(len=10)::method
#if NENER>0
  integer::irad
#endif
  associate(r=>s%r,g=>s%g,m=>s%m)

  if(r%levelmax_frame==0)then
     nlevelmax_frame=r%nlevelmax
  else if (r%levelmax_frame.gt.r%nlevelmax)then
     nlevelmax_frame=r%nlevelmax
  else
     nlevelmax_frame=r%levelmax_frame
  endif

  ! Method / shader booleans
  shader = r%shader_frame(ind_proj)
  method = r%method_frame(ind_proj)
  is_square    = (trim(shader).eq.'square')
  is_sphere    = (trim(shader).eq.'sphere')
  is_cube      = (trim(shader).eq.'cube')   ! Treated as square here (full cube shader not ported)
  is_mean_mass = (trim(method).eq.'mean_mass')
  is_mean_dens = (trim(method).eq.'mean_dens')
  is_mean_vol  = (trim(method).eq.'mean_vol')
  is_sum       = (trim(method).eq.'sum')
  is_min       = (trim(method).eq.'min')
  is_max       = (trim(method).eq.'max')
  is_mean      = (method(1:4).eq.'mean')
  if(.not.(is_mean.or.is_sum.or.is_min.or.is_max)) is_mean_mass = .true.
  perspective  = r%perspective_camera(ind_proj)

  ! Conversion factor from user units to cgs units
  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Camera evolution timer (cosmological vs non-cosmological).
  if(r%cosmo)then
     timer     = g%aexp
     timer_end = r%aendmov
  else
     timer     = g%t
     timer_end = r%tendmov
  endif

  ! Camera evolution — formulas mirror ramses/amr/movie.f90:237-259 exactly so
  ! the two codes produce identical camera tracks for matching namelists.
  ! Note: in the cosmological branch, ramses gates phi on tstart_theta_camera
  ! (not tstart_phi_camera); the non-cosmological branch uses tstart_phi_camera.
  ! dist_cam is always gated on tstart_theta_camera / tend_theta_camera in both
  ! modes — there is no separate tstart_dist_camera.

  ! Default base distance to boxlen when unset.
  base_dist = r%dist_camera(ind_proj)
  if(base_dist.le.0d0) base_dist = r%boxlen

  ! Default tend_*_camera to the global movie endpoint when <=0.
  tend_theta_eff = r%tend_theta_camera(ind_proj)
  if(tend_theta_eff.le.0d0) tend_theta_eff = timer_end
  tend_phi_eff   = r%tend_phi_camera(ind_proj)
  if(tend_phi_eff  .le.0d0) tend_phi_eff   = timer_end

  ! Span = movie duration in the active timer.
  if(r%cosmo)then
     span = r%aendmov - r%astartmov
  else
     span = r%tendmov - r%tstartmov
  endif

  if(span.le.0d0)then
     theta_cam = r%theta_camera(ind_proj)*pi/180d0
     phi_cam   = r%phi_camera  (ind_proj)*pi/180d0
     dist_cam  = base_dist
  else
     theta_cam = r%theta_camera(ind_proj)*pi/180d0 &
          & + min(max(timer-r%tstart_theta_camera(ind_proj),0d0),tend_theta_eff) &
          &   *r%dtheta_camera(ind_proj)*pi/180d0/span
     if(r%cosmo)then
        ! ramses cosmo branch: phi gated on tstart_THETA_camera (see ramses movie.f90:246)
        phi_cam = r%phi_camera(ind_proj)*pi/180d0 &
             & + min(max(timer-r%tstart_theta_camera(ind_proj),0d0),tend_phi_eff) &
             &   *r%dphi_camera(ind_proj)*pi/180d0/span
     else
        phi_cam = r%phi_camera(ind_proj)*pi/180d0 &
             & + min(max(timer-r%tstart_phi_camera(ind_proj),0d0),tend_phi_eff) &
             &   *r%dphi_camera(ind_proj)*pi/180d0/span
     endif
     dist_cam = base_dist &
          & + min(max(timer-r%tstart_theta_camera(ind_proj),0d0),tend_theta_eff) &
          &   *r%ddist_camera(ind_proj)/span
  endif

  ! Default focal_camera — matches ramses movie.f90:258. Mutates r%focal_camera
  ! so the perspective routines (project_point / project_particles) downstream
  ! use the same defaulted value.
  if(r%focal_camera(ind_proj).le.0d0 .or. r%focal_camera(ind_proj).gt.base_dist) then
     r%focal_camera(ind_proj) = dist_cam
  endif

  ! Field-of-view from raw deltax_frame including the cosmological /aexp term,
  ! independent of proj_axis (matches ramses movie.f90:234,259).
  fov_camera = atan( &
       & (r%deltax_frame(ind_proj*2-1)+r%deltax_frame(ind_proj*2)/g%aexp)/2d0 &
       & /max(r%focal_camera(ind_proj),1d-30))

  ! Compute frame centre (with optional cosmological polynomial drift in aexp).
  ! After the axis-swap, (xcen,ycen,zcen) are screen-x, screen-y, line-of-sight centres
  ! expressed in world coordinates.
  if(r%proj_axis(ind_proj:ind_proj).eq.'x')then
     xcen=r%ycentre_frame(ind_proj*4-3)+r%ycentre_frame(ind_proj*4-2)*g%aexp+r%ycentre_frame(ind_proj*4-1)*g%aexp**2+r%ycentre_frame(ind_proj*4)*g%aexp**3
     ycen=r%zcentre_frame(ind_proj*4-3)+r%zcentre_frame(ind_proj*4-2)*g%aexp+r%zcentre_frame(ind_proj*4-1)*g%aexp**2+r%zcentre_frame(ind_proj*4)*g%aexp**3
     zcen=r%xcentre_frame(ind_proj*4-3)+r%xcentre_frame(ind_proj*4-2)*g%aexp+r%xcentre_frame(ind_proj*4-1)*g%aexp**2+r%xcentre_frame(ind_proj*4)*g%aexp**3
     box_xidx=2; box_yidx=3; box_zidx=1
  elseif(r%proj_axis(ind_proj:ind_proj).eq.'y')then
     xcen=r%xcentre_frame(ind_proj*4-3)+r%xcentre_frame(ind_proj*4-2)*g%aexp+r%xcentre_frame(ind_proj*4-1)*g%aexp**2+r%xcentre_frame(ind_proj*4)*g%aexp**3
     ycen=r%zcentre_frame(ind_proj*4-3)+r%zcentre_frame(ind_proj*4-2)*g%aexp+r%zcentre_frame(ind_proj*4-1)*g%aexp**2+r%zcentre_frame(ind_proj*4)*g%aexp**3
     zcen=r%ycentre_frame(ind_proj*4-3)+r%ycentre_frame(ind_proj*4-2)*g%aexp+r%ycentre_frame(ind_proj*4-1)*g%aexp**2+r%ycentre_frame(ind_proj*4)*g%aexp**3
     box_xidx=1; box_yidx=3; box_zidx=2
  else
     xcen=r%xcentre_frame(ind_proj*4-3)+r%xcentre_frame(ind_proj*4-2)*g%aexp+r%xcentre_frame(ind_proj*4-1)*g%aexp**2+r%xcentre_frame(ind_proj*4)*g%aexp**3
     ycen=r%ycentre_frame(ind_proj*4-3)+r%ycentre_frame(ind_proj*4-2)*g%aexp+r%ycentre_frame(ind_proj*4-1)*g%aexp**2+r%ycentre_frame(ind_proj*4)*g%aexp**3
     zcen=r%zcentre_frame(ind_proj*4-3)+r%zcentre_frame(ind_proj*4-2)*g%aexp+r%zcentre_frame(ind_proj*4-1)*g%aexp**2+r%zcentre_frame(ind_proj*4)*g%aexp**3
     box_xidx=1; box_yidx=2; box_zidx=3
  endif

  ! Compute frame size (screen extents).
  if(r%proj_axis(ind_proj:ind_proj).eq.'x')then
     delx=r%deltay_frame(ind_proj*2-1)+r%deltay_frame(ind_proj*2)/g%aexp
     dely=r%deltaz_frame(ind_proj*2-1)+r%deltaz_frame(ind_proj*2)/g%aexp
     delz=r%deltax_frame(ind_proj*2-1)+r%deltax_frame(ind_proj*2)/g%aexp
  elseif(r%proj_axis(ind_proj:ind_proj).eq.'y')then
     delx=r%deltax_frame(ind_proj*2-1)+r%deltax_frame(ind_proj*2)/g%aexp
     dely=r%deltaz_frame(ind_proj*2-1)+r%deltaz_frame(ind_proj*2)/g%aexp
     delz=r%deltay_frame(ind_proj*2-1)+r%deltay_frame(ind_proj*2)/g%aexp
  else
     delx=r%deltax_frame(ind_proj*2-1)+r%deltax_frame(ind_proj*2)/g%aexp
     dely=r%deltay_frame(ind_proj*2-1)+r%deltay_frame(ind_proj*2)/g%aexp
     delz=r%deltaz_frame(ind_proj*2-1)+r%deltaz_frame(ind_proj*2)/g%aexp
  endif

  dx_frame=delx/dble(r%nw_frame)
  dy_frame=dely/dble(r%nh_frame)

  ! Clamp centre against the correct world-axis bound for each role.
  xcen=min(max(xcen,delx/2.0d0),r%box_size(box_xidx)-delx/2.0d0)
  ycen=min(max(ycen,dely/2.0d0),r%box_size(box_yidx)-dely/2.0d0)
  zcen=min(max(zcen,delz/2.0d0),r%box_size(box_zidx)-delz/2.0d0)

  xleft_frame =xcen-delx/2.0d0
  xright_frame=xcen+delx/2.0d0
  yleft_frame =ycen-dely/2.0d0
  yright_frame=ycen+dely/2.0d0
  zleft_frame =zcen-delz/2.0d0
  zright_frame=zcen+delz/2.0d0

  ! ---------- HYDRO CELL PROJECTION ----------
  ! Skip the AMR walk entirely when the requested variable is particle-only.
  if(ind_var.le.nvar+nrtgrp .and. r%hydro) then

  do ilevel=r%levelmin,nlevelmax_frame

     dx=r%boxlen/2**ilevel
     dx_proj=(dx/2.0d0)*r%smooth_frame(ind_proj)

     do igrid=m%head(ilevel),m%tail(ilevel)

        do ind=1,twotondim

           ! Cell centre in world coordinates.
           do idim=1,ndim
              nstride=2**(idim-1)
              xx(idim)=(2*m%grid(igrid)%ckey(idim)+MOD((ind-1)/nstride,2)+0.5d0)*dx-m%skip(idim)
           end do

           ok=(.NOT.m%grid(igrid)%refined(ind)).or.(ilevel==nlevelmax_frame)

#ifdef HYDRO
           ! varmin/varmax filter on ivar_frame (matches ramses' behaviour)
           if(ok .and. r%ivar_frame.ge.1 .and. r%ivar_frame.le.nvar)then
              uvar = m%uold(ind,r%ivar_frame,igrid)
              if(r%ivar_frame.gt.1) uvar = uvar/max(real(m%uold(ind,1,igrid),kind=8),r%smallr)
              ok = ok .and. (uvar.ge.r%varmin_frame(ind_proj))
              ok = ok .and. (uvar.le.r%varmax_frame(ind_proj))
           endif
#endif

           if(.not.ok) cycle

           ! Project cell centre into screen coords with optional camera rotation.
           call project_point(xx,xcen,ycen,zcen,theta_cam,phi_cam, &
                & dist_cam,r%focal_camera(ind_proj),fov_camera,perspective, &
                & r%proj_axis(ind_proj:ind_proj),xcentre,ycentre,zcentre,pers_corr,ok)
           if(.not.ok) cycle
           dx_proj = (dx/2.0d0)*r%smooth_frame(ind_proj)*pers_corr

#if NDIM>2
           xleft  = xcentre-dx_proj
           xright = xcentre+dx_proj
           yleft  = ycentre-dx_proj
           yright = ycentre+dx_proj
           zleft  = zcentre-dx/2.0d0
           zright = zcentre+dx/2.0d0
           if(    xright.lt.xleft_frame.or.xleft.ge.xright_frame.or.&
                & yright.lt.yleft_frame.or.yleft.ge.yright_frame.or.&
                & zright.lt.zleft_frame.or.zleft.ge.zright_frame)cycle
#else
           xleft  = xcentre-dx_proj
           xright = xcentre+dx_proj
           yleft  = ycentre-dx_proj
           yright = ycentre+dx_proj
           if(    xright.lt.xleft_frame.or.xleft.ge.xright_frame.or.&
                & yright.lt.yleft_frame.or.yleft.ge.yright_frame)cycle
#endif

           ! Pixel index range that may overlap this cell.
           if(xleft>xleft_frame)then
              imin=min(int((xleft-xleft_frame)/dx_frame)+1,r%nw_frame)
           else
              imin=1
           endif
           imax=min(int((xright-xleft_frame)/dx_frame)+1,r%nw_frame)
           if(yleft>yleft_frame)then
              jmin=min(int((yleft-yleft_frame)/dy_frame)+1,r%nh_frame)
           else
              jmin=1
           endif
           jmax=min(int((yright-yleft_frame)/dy_frame)+1,r%nh_frame)

#if NDIM>2
           dz_cell=min(zright_frame,zright)-max(zleft_frame,zleft)
#endif
           do ii=imin,imax
              xxleft =xleft_frame+dble(ii-1)*dx_frame
              xxright=xxleft+dx_frame
              xxcen  =0.5d0*(xxleft+xxright)
              dx_cell=min(xxright,xright)-max(xxleft,xleft)
              do jj=jmin,jmax
                 yyleft =yleft_frame+dble(jj-1)*dy_frame
                 yyright=yyleft+dy_frame
                 yycen  =0.5d0*(yyleft+yyright)
                 dy_cell=min(yyright,yright)-max(yyleft,yleft)
                 ! Pixel-cell intersection volume in world units.
                 dvol=dx_cell*dy_cell
#if NDIM>2
                 dvol=dvol*dz_cell
#endif
                 ! Shader gating: square uses the box intersection, sphere
                 ! requires the pixel centre to fall within the cell's disk.
                 if(is_sphere)then
                    xpc=xxcen-xcentre
                    ypc=yycen-ycentre
                    if((xpc*xpc+ypc*ypc).gt.(dx_proj*dx_proj)) cycle
                 endif

                 ind_map=ii+(jj-1)*r%nw_frame

#ifdef HYDRO
                 ! Pick the projection weight according to the method.
                 if(is_mean_mass)then
                    weight = dvol*max(dble(m%uold(ind,1,igrid)),r%smallr)
                 elseif(is_mean_dens)then
                    weight = dvol*max(dble(m%uold(ind,1,igrid)),r%smallr)
                 elseif(is_mean_vol)then
                    weight = dvol
                 elseif(is_sum)then
                    weight = 1d0
                 else
                    weight = dvol*max(dble(m%uold(ind,1,igrid)),r%smallr)
                 endif

                 ! Compute uvar for the requested variable.
                 if(ind_var==0)then
                    ! Weight map only - just accumulate the weight.
                    if(is_min)then
                       map(ind_map)=min(map(ind_map),weight)
                    else if(is_max)then
                       map(ind_map)=max(map(ind_map),weight)
                    else
                       map(ind_map)=map(ind_map)+weight
                    endif
                    cycle
                 else if(ind_var==1)then
                    uvar = max(dble(m%uold(ind,1,igrid)),r%smallr)
                 else if(ind_var.ge.2 .and. ind_var.le.ndim+1)then
                    uvar = m%uold(ind,ind_var,igrid)/max(dble(m%uold(ind,1,igrid)),r%smallr)
                 else if(ind_var==5)then
                    ekk=0d0
                    do idim=1,3
                       ekk=ekk+0.5d0*m%uold(ind,idim+1,igrid)**2/max(dble(m%uold(ind,1,igrid)),r%smallr)
                    enddo
#if NENER>0
                    do irad=1,nener
                      ekk = ekk+m%uold(ind,5+irad,igrid)
                    enddo
#endif
#ifdef MHD
                    ekk=ekk+0.125d0*sum((m%bold(ind,1:3,igrid)+m%bold(ind,4:6,igrid))**2)
#endif
                    temp=(r%gamma-1d0)*(m%uold(ind,5,igrid)-ekk)
                    temp=max(temp/max(dble(m%uold(ind,1,igrid)),r%smallr),r%smallc**2)*scale_T2
                    uvar = temp
#ifdef DO_RT
                 else if(ind_var.ge.nvar+1 .and. ind_var.le.nvar+nrtgrp)then
                    uvar = m%rtuold(ind,1+(ind_var-nvar-1)*(ndim+1),igrid)*g%rt_c(ilevel)
#endif
                 else
                    uvar = m%uold(ind,ind_var,igrid)
                 endif

                 if(is_min)then
                    map(ind_map)=min(map(ind_map),uvar)
                 else if(is_max)then
                    map(ind_map)=max(map(ind_map),uvar)
                 else if(is_sum)then
                    map(ind_map)=map(ind_map)+uvar
                 else
                    map(ind_map)=map(ind_map)+weight*uvar
                 endif
#endif
              end do
           end do

        end do
        ! End loop over cells

     end do
     ! End loop over grids
  end do
  ! End loop over levels

  endif ! hydro section

  ! ---------- PARTICLE PROJECTION ----------
#if NDIM>2
  ! Dark matter: ind_var == NVAR+nrtgrp+1
  if(ind_var.eq.NVAR+nrtgrp+1)then
     call project_particles(s%p, .true., r%zoom_only_frame(ind_proj), r%mass_cut_refine, &
          & ind_proj,r,xcen,ycen,zcen,theta_cam,phi_cam,dist_cam, &
          & r%focal_camera(ind_proj),fov_camera,perspective, &
          & xleft_frame,xright_frame,yleft_frame,yright_frame, &
          & zleft_frame,zright_frame,dx_frame,dy_frame,r%nw_frame,r%nh_frame, &
          & map,map_size)
  endif

  ! Stars: ind_var == NVAR+nrtgrp+2
  if(ind_var.eq.NVAR+nrtgrp+2)then
     call project_particles(s%star, .false., .false., -1d0, &
          & ind_proj,r,xcen,ycen,zcen,theta_cam,phi_cam,dist_cam, &
          & r%focal_camera(ind_proj),fov_camera,perspective, &
          & xleft_frame,xright_frame,yleft_frame,yright_frame, &
          & zleft_frame,zright_frame,dx_frame,dy_frame,r%nw_frame,r%nh_frame, &
          & map,map_size)
  endif
#endif

  end associate

end subroutine output_frame
!=======================================================================
!=======================================================================
!=======================================================================
!=======================================================================
subroutine project_point(xx,xcen,ycen,zcen,theta_cam,phi_cam,dist_cam, &
     & focal_cam,fov_cam,perspective,proj_ax,xcentre,ycentre,zcentre,pers_corr,ok)
  ! Rotate xx around (xcen,ycen,zcen) by (theta,phi) and project to screen
  ! coordinates matching ramses' convention. proj_ax is one of 'x','y','z'.
  use amr_parameters, only: ndim
  implicit none
  real(kind=8),dimension(1:ndim),intent(in)::xx
  real(kind=8),intent(in)::xcen,ycen,zcen,theta_cam,phi_cam,dist_cam,focal_cam,fov_cam
  logical,intent(in)::perspective
  character(len=1),intent(in)::proj_ax
  real(kind=8),intent(out)::xcentre,ycentre,zcentre,pers_corr
  logical,intent(inout)::ok

  real(kind=8)::x1,x2,x3,xtmp,ytmp,ztmp,alpha,beta

  pers_corr=1d0

#if NDIM>2
  x1 = xx(1)-xcen
  x2 = xx(2)-ycen
  x3 = xx(3)-zcen
  ! XY rotation by theta_cam
  xtmp = cos(theta_cam)*x1+sin(theta_cam)*x2
  ytmp = cos(theta_cam)*x2-sin(theta_cam)*x1
  x1 = xtmp
  x2 = ytmp
  ! YZ rotation by phi_cam
  ytmp = cos(phi_cam)*x2+sin(phi_cam)*x3
  ztmp = cos(phi_cam)*x3-sin(phi_cam)*x2
  x2 = ytmp
  x3 = ztmp

  if(proj_ax.eq.'x')then
     if(dist_cam-x1.lt.0d0) then; ok=.false.; return; endif
     if(perspective)then
        alpha = atan(x2/(dist_cam-x1))
        beta  = atan(x3/(dist_cam-x1))
        if(abs(alpha)/2d0.gt.fov_cam) then; ok=.false.; return; endif
        if(abs(beta)/2d0.gt.fov_cam)  then; ok=.false.; return; endif
        pers_corr = focal_cam/(dist_cam-x1)
        x2 = x2*pers_corr
        x3 = x3*pers_corr
     endif
     xcentre = x2+ycen
     ycentre = x3+zcen
     zcentre = x1+xcen
  elseif(proj_ax.eq.'y')then
     if(dist_cam-x2.lt.0d0) then; ok=.false.; return; endif
     if(perspective)then
        alpha = atan(x1/(dist_cam-x2))
        beta  = atan(x3/(dist_cam-x2))
        if(abs(alpha)/2d0.gt.fov_cam) then; ok=.false.; return; endif
        if(abs(beta)/2d0.gt.fov_cam)  then; ok=.false.; return; endif
        pers_corr = focal_cam/(dist_cam-x2)
        x1 = x1*pers_corr
        x3 = x3*pers_corr
     endif
     xcentre = x1+xcen
     ycentre = x3+zcen
     zcentre = x2+ycen
  else
     if(dist_cam-x3.lt.0d0) then; ok=.false.; return; endif
     if(perspective)then
        alpha = atan(x1/(dist_cam-x3))
        beta  = atan(x2/(dist_cam-x3))
        if(abs(alpha)/2d0.gt.fov_cam) then; ok=.false.; return; endif
        if(abs(beta)/2d0.gt.fov_cam)  then; ok=.false.; return; endif
        pers_corr = focal_cam/(dist_cam-x3)
        x1 = x1*pers_corr
        x2 = x2*pers_corr
     endif
     xcentre = x1+xcen
     ycentre = x2+ycen
     zcentre = x3+zcen
  endif
#elif NDIM>1
  x1 = xx(1)-xcen
  x2 = xx(2)-ycen
  xtmp = cos(theta_cam)*x1+sin(theta_cam)*x2
  ytmp = cos(theta_cam)*x2-sin(theta_cam)*x1
  xcentre = xtmp+xcen
  ycentre = ytmp+ycen
  zcentre = 0d0
#endif
end subroutine project_point
!=======================================================================
!=======================================================================
!=======================================================================
!=======================================================================
subroutine project_particles(p,is_dm,zoom_only_proj,mass_cut, &
     & ind_proj,r,xcen,ycen,zcen,theta_cam,phi_cam,dist_cam, &
     & focal_cam,fov_cam,perspective, &
     & xleft_frame,xright_frame,yleft_frame,yright_frame, &
     & zleft_frame,zright_frame,dx_frame,dy_frame,nw,nh,map,map_size)
  use amr_parameters, only: ndim
  use amr_commons, only: run_t
  use pm_commons, only: part_t
  implicit none
  type(part_t)::p
  type(run_t)::r
  logical,intent(in)::is_dm,zoom_only_proj,perspective
  real(kind=8),intent(in)::mass_cut
  integer,intent(in)::ind_proj,nw,nh,map_size
  real(kind=8),intent(in)::xcen,ycen,zcen,theta_cam,phi_cam,dist_cam,focal_cam,fov_cam
  real(kind=8),intent(in)::xleft_frame,xright_frame,yleft_frame,yright_frame
  real(kind=8),intent(in)::zleft_frame,zright_frame,dx_frame,dy_frame
  real(kind=8),dimension(1:map_size),intent(inout)::map

  integer::ip,ii,jj,ind_map
  real(kind=8)::xmap,ymap,zmap,pers_corr
  real(kind=8),dimension(1:ndim)::xx
  logical::ok

  if(p%npart.le.0) return
  if(.not.allocated(p%xp)) return
  if(.not.allocated(p%mp)) return

  do ip=1,p%npart
     xx(1)=p%xp(ip,1)
#if NDIM>1
     xx(2)=p%xp(ip,2)
#endif
#if NDIM>2
     xx(3)=p%xp(ip,3)
#endif
     ok=.true.
     call project_point(xx,xcen,ycen,zcen,theta_cam,phi_cam,dist_cam, &
          & focal_cam,fov_cam,perspective,r%proj_axis(ind_proj:ind_proj), &
          & xmap,ymap,zmap,pers_corr,ok)
     if(.not.ok) cycle

     if(xmap.lt.xleft_frame .or. xmap.ge.xright_frame) cycle
     if(ymap.lt.yleft_frame .or. ymap.ge.yright_frame) cycle
#if NDIM>2
     if(zmap.lt.zleft_frame .or. zmap.ge.zright_frame) cycle
#endif

     ii = max(1,min(int((xmap-xleft_frame)/dx_frame)+1, nw))
     jj = max(1,min(int((ymap-yleft_frame)/dy_frame)+1, nh))
     ind_map = ii+(jj-1)*nw

     if(is_dm .and. mass_cut.gt.0d0 .and. zoom_only_proj)then
        if(p%mp(ip).ge.mass_cut) cycle
     endif

     map(ind_map) = map(ind_map) + p%mp(ip)
  end do
end subroutine project_particles
!=======================================================================
!=======================================================================
!=======================================================================
!=======================================================================
subroutine set_movie_vars(r)
  use amr_parameters, only: ndim
  use amr_commons, only: run_t
  use hydro_parameters, only: nvar
  use rt_parameters, only: nrtgrp
  ! Map textual movie variable names to the integer index used by the
  ! variable dispatch in output_frame.
  type(run_t)::r
  integer::tmp,ll
  character(LEN=5)::dummy

  if(ANY(r%movie_vars_txt=='dens '))r%movie_vars(1)=1
  if(ANY(r%movie_vars_txt=='vx   '))r%movie_vars(2)=1
  if(ANY(r%movie_vars_txt=='vy   '))r%movie_vars(3)=1
  if(ANY(r%movie_vars_txt=='vz   '))r%movie_vars(4)=1
  if(ANY(r%movie_vars_txt=='temp '))r%movie_vars(5)=1
#if NVAR>5
  do ll=6,nvar
     write(dummy,'(I3.1)') ll
     if(ANY(r%movie_vars_txt=='var'//trim(adjustl(dummy))//' '))r%movie_vars(ll)=1
  end do
#endif
  do ll=nvar+1,nvar+nrtgrp
     tmp=ll-nvar
     write(dummy,'(I3.1)') tmp
     if(ANY(r%movie_vars_txt=='Fp'//trim(adjustl(dummy))//' '))r%movie_vars(ll)=1
  end do
  if(ANY(r%movie_vars_txt=='dm   '))r%movie_vars(NVAR+nrtgrp+1)=1
  if(ANY(r%movie_vars_txt=='stars'))r%movie_vars(NVAR+nrtgrp+2)=1

end subroutine set_movie_vars
!=======================================================================
!=======================================================================
!=======================================================================
!=======================================================================
end module movie_module
