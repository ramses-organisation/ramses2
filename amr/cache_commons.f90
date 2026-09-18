module cache_commons
  use amr_parameters, only: ndim, twotondim, nbin
  use hydro_parameters, only: nvar
  use rt_parameters, only: nrtvar, nrtgrp
  use cr_parameters, only: ncruvar
  use call_back

  ! Communication-related taghs
  integer::flush_tag=1000,msg_tag=100,request_tag=10
  integer::flush_tag_clump=2000,msg_tag_clump=200,request_tag_clump=20

  ! Software cache parameters
  integer,parameter::operation_initflag=1,operation_upload=2,operation_godunov=3,operation_smooth=4
  integer,parameter::operation_hydro=5,operation_refine=6,operation_derefine=7,operation_loadbalance=8
  integer,parameter::operation_phi=9,operation_rho=10,operation_multipole=11,operation_cg=12
  integer,parameter::operation_build_mg=13,operation_restrict_mask=14,operation_mg=15
  integer,parameter::operation_restrict_res=16,operation_scan=17,operation_interpol=18
  integer,parameter::operation_split=19,operation_kick=20

  integer,parameter::domain_decompos_amr=1,domain_decompos_mg=2

  ! Message size
  integer,parameter::ntilemax=16  ! Fetch message buffer size
  integer,parameter::nflushmax=128  ! Flush message buffer size
  integer,parameter::nbuffermax=128  ! Initial number of send buffers

  ! Combiner rules
  integer,parameter::COMBINER_EXIST=1
  integer,parameter::COMBINER_CREATE=2

  ! Message element type
  type msg_int4
     integer(kind=4),dimension(1:twotondim)::int4
  end type msg_int4
  type msg_realdp
     integer(kind=4),dimension(1:twotondim)::int4
#ifdef HYDRO
     real(kind=8),dimension(1:twotondim,1:nvar)::realdp
#ifdef TRCFLX
     real(kind=8),dimension(1:twotondim,1:2*ndim+1)::realdp_mflux
#endif
#endif
#ifdef MHD
     real(kind=8),dimension(1:twotondim,1:6)::realdp_mhd
#endif
#ifdef DO_RT
     real(kind=8),dimension(1:twotondim,1:nrtvar)::realdp_rt
#endif
#ifdef DO_CR
     real(kind=8),dimension(1:twotondim,1:ncruvar)::realdp_cr
#endif
  end type msg_realdp
  type msg_small_realdp
     real(kind=8),dimension(1:twotondim)::realdp
  end type msg_small_realdp
  type msg_int4_small_realdp
     integer(kind=4),dimension(1:twotondim)::flg
     integer(kind=4),dimension(1:twotondim)::ref
     real(kind=8),dimension(1:twotondim)::realdp
  end type msg_int4_small_realdp
  type msg_twin_realdp
     real(kind=8),dimension(1:twotondim)::realdp_phi
     real(kind=8),dimension(1:twotondim)::realdp_dis
  end type msg_twin_realdp
  type msg_three_realdp
     real(kind=8),dimension(1:twotondim)::realdp_phi
     real(kind=8),dimension(1:twotondim)::realdp_phi_old
     real(kind=8),dimension(1:twotondim)::realdp_dis
  end type msg_three_realdp
  type msg_nvar_realdp
     real(kind=8),dimension(1:twotondim,1:nvar)::realdp_hydro
  end type msg_nvar_realdp
  ! Explicit mflux message types (keep mflux transfers self-contained)
  type msg_hydro_mflux
     real(kind=8),dimension(1:twotondim,1:nvar)::realdp_hydro
#ifdef TRCFLX
     real(kind=8),dimension(1:twotondim,1:2*ndim+1)::realdp_mflux
#endif
  end type msg_hydro_mflux
  type msg_upload_hydro_mflux_mhd
     integer(kind=4),dimension(1:twotondim)::int4
     real(kind=8),dimension(1:twotondim,1:nvar)::realdp_hydro
#ifdef TRCFLX
     real(kind=8),dimension(1:twotondim,1:2*ndim+1)::realdp_mflux
#endif
#ifdef MHD
     real(kind=8),dimension(1:twotondim,1:6)::realdp_mhd
#endif
  end type msg_upload_hydro_mflux_mhd
  type msg_large_realdp
     integer(kind=4),dimension(1:twotondim)::int4
#ifdef HYDRO
     real(kind=8),dimension(1:twotondim,1:nvar)::realdp_hydro
#endif
#ifdef DO_RT
     real(kind=8),dimension(1:twotondim,1:nrtvar)::realdp_rt
#endif
#ifdef DO_CR
     real(kind=8),dimension(1:twotondim,1:ncruvar)::realdp_cr
#endif
#ifdef MHD
     real(kind=8),dimension(1:twotondim,1:6)::realdp_mhd
#endif
#ifdef GRAV
     real(kind=8),dimension(1:twotondim,1:ndim+2)::realdp_poisson
#endif
#ifdef TRCFLX
     real(kind=8),dimension(1:twotondim,1:2*ndim+1)::realdp_mflux
#endif
  end type msg_large_realdp
  type msg_rt_emissivity_realdp
     real(kind=8),dimension(1:twotondim,1:nrtgrp)::realdp
  end type msg_rt_emissivity_realdp
  type msg_saddle_clump
     integer(kind=8)::nbor
     real(kind=8)::dens
  end type msg_saddle_clump
  type msg_merge_clump
     integer(kind=8)::npeak
     real(kind=8)::mdens
  end type msg_merge_clump
  type msg_prop_clump
     integer::ncell
     integer(kind=8)::ind
     real(kind=8)::dens
     real(kind=8)::mass
     real(kind=8)::vol
     real(kind=8)::rad
     real(kind=8),dimension(1:ndim)::pos
     real(kind=8),dimension(1:ndim)::vel
  end type msg_prop_clump
  type msg_halo_clump
     integer::ihalo
     real(kind=8)::mhalo
  end type msg_halo_clump
  type msg_mbin_clump
     integer::npart
     integer::lev
     integer(kind=8)::ind
     real(kind=8)::rad
     real(kind=8),dimension(1:ndim)::pos
     real(kind=8),dimension(1:nbin)::mbin
  end type msg_mbin_clump
  type msg_unbind_clump
     integer::lev
     integer(kind=8)::ind
     real(kind=8)::rad
     real(kind=8)::mass
     real(kind=8),dimension(1:ndim)::pos
     real(kind=8),dimension(1:ndim)::vel
     real(kind=8),dimension(1:nbin)::mbin
  end type msg_unbind_clump
  type msg_sink_clump
     integer::lev
     real(kind=8),dimension(1:ndim)::pos
     real(kind=8),dimension(1:ndim)::vel
     real(kind=8),dimension(1:ndim)::acc
  end type msg_sink_clump
  type msg_tree_clump
     integer::id
     integer::lev
     real(kind=8)::rad
     real(kind=8)::mass
     real(kind=8),dimension(1:ndim)::pos
     real(kind=8),dimension(1:ndim)::vel
     real(kind=8),dimension(1:ndim)::acc
  end type msg_tree_clump
  type msg_tree_minid
     integer::id
     integer(kind=8)::ind
     real(kind=8)::mass
  end type msg_tree_minid

  ! Cache call back functions
  type(cache_f)       ::pack_fetch
  type(cache_unpack_f)::unpack_fetch
  type(cache_f)       ::pack_flush
  type(cache_unpack_f)::unpack_flush
  type(cache_init_f)  ::init_flush
  type(cache_bound_f) ::init_bound

  type(cache_f_clump)       ::pack_fetch_clump
  type(cache_unpack_f_clump)::unpack_fetch_clump
  type(cache_f_clump)       ::pack_flush_clump
  type(cache_unpack_f_clump)::unpack_flush_clump
  type(cache_init_f_clump)  ::init_flush_clump

end module cache_commons
