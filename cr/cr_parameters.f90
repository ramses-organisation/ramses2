module cr_parameters
  use amr_parameters, only: ndim

#ifdef NCRGRP
  integer,parameter::ncrgrp=NCRGRP       ! # of CR groups (set in Makefile)
#else
  integer,parameter::ncrgrp=1
#endif
  integer,parameter::ncruvar=ncrgrp*ndim    ! # of CR flux variables, stored in cruold and crunew
  real(kind=8),parameter::smallEcr=1d-50 !  Minimum CR density

end module cr_parameters
