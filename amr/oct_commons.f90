module oct_commons
  use amr_parameters

  ! Oct object
  type oct
     integer(kind=8),dimension(1:nhilbert)::hkey
     integer(kind=4),dimension(1:ndim)::ckey
     logical(kind=1),dimension(1:twotondim)::refined
     integer(kind=4)::lev
  end type oct

end module oct_commons
