module hydro_parameters
  use amr_parameters, only: ndim

  ! Number of independant variables
#ifndef NENER
  integer,parameter::nener=0
#else
  integer,parameter::nener=NENER
#endif
#ifndef NVAR
  integer,parameter::nvar=5+nener
#else
  integer,parameter::nvar=NVAR
#endif

#ifdef MHD
  integer,parameter::nprim=NVAR+3
  integer,parameter::ie=8
#else
  integer,parameter::nprim=NVAR
  integer,parameter::ie=5
#endif

#ifdef NION
  integer,parameter::nion=NION  ! # of ionization fractions species
#else
  integer,parameter::nion=1
#endif

  integer,parameter::solver_llf=1
  integer,parameter::solver_hll=2
  integer,parameter::solver_hllc=3
  integer,parameter::solver_hlld=4
  integer,parameter::solver_roe=5
  integer,parameter::solver_upwind=6
#if defined(_METAL) || defined(_CUDA)
  integer,parameter::solver_uct_hlld=7
#endif
#ifdef _CUDA
  integer,parameter::solver_uct_llf=8
  integer,parameter::solver_uct_hll=9
#endif

  integer,parameter::solver2d_llf=1
  integer,parameter::solver2d_hllf=2
  integer,parameter::solver2d_hlla=3
  integer,parameter::solver2d_hlld=4
  integer,parameter::solver2d_roe=5
  integer,parameter::solver2d_upwind=6

end module hydro_parameters

module const

  ! Some useful constant
  real(kind=8),parameter ::bigreal = 1.0d+30
  real(kind=8),parameter ::zero = 0.0
  real(kind=8),parameter ::one = 1.0
  real(kind=8),parameter ::two = 2.0
  real(kind=8),parameter ::three = 3.0
  real(kind=8),parameter ::four = 4.0
  real(kind=8),parameter ::two3rd = 0.6666666666666667
  real(kind=8),parameter ::half = 0.5
  real(kind=8),parameter ::third = 0.33333333333333333
  real(kind=8),parameter ::forth = 0.25
  real(kind=8),parameter ::sixth = 0.16666666666666667

end module const
