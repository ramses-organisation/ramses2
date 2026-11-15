module init_cr_module

contains
!###############################################
!###############################################
!###############################################
!###############################################
recursive subroutine r_init_cr(pst)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_INIT_CR,pst%iUpper+1)
     call r_init_cr(pst%pLower)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call init_cr(pst%s%r,pst%s%m)
  endif

end subroutine r_init_cr
!###############################################
!###############################################
!###############################################
!###############################################
subroutine init_cr(r,m)
  use amr_parameters, only: nsubgrid
  use amr_commons, only: run_t,mesh_t
  implicit none
  type(run_t)::r
  type(mesh_t)::m

  ! Initialise workspace for hydro kernels
  associate(h => m%cr_w)
    call h%kernel%init(2*nsubgrid)
  end associate

end subroutine init_cr
!###############################################
!###############################################
!###############################################
!###############################################
end module init_cr_module
