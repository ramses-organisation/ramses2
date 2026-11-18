!################################################################
!################################################################
!################################################################
!################################################################
module color_module
  !--------------------------------------------------------------
  ! ANSI colors for the standard output.
  ! Colors are used only if stdout is a terminal. This can be
  ! overridden with the environment variable RAMSES_COLOR set to
  ! always or never. A non-empty NO_COLOR disables colors.
  ! When colors are off, all functions return an empty string,
  ! so that the output is unchanged.
  !--------------------------------------------------------------
  implicit none
  private
  public :: init_color, c_bold, c_dim, c_under, c_cyan, c_green, c_yellow, c_reset, c_mark

  logical :: use_color = .false.

contains
  !--------------------------------------------------------------
  subroutine init_color()
    use iso_c_binding, only: c_int
    interface
       integer(c_int) function isatty(fd) bind(c,name='isatty')
         import :: c_int
         integer(c_int), value :: fd
       end function isatty
    end interface
    character(len=16)::val
    integer::ilen,istat

    use_color = (isatty(1_c_int)==1)

    call get_environment_variable('RAMSES_COLOR',val,ilen,istat)
    if(istat==0)then
       if(trim(val)=='always')use_color=.true.
       if(trim(val)=='never') use_color=.false.
    endif

    call get_environment_variable('NO_COLOR',val,ilen,istat)
    if(istat==0.and.ilen>0)use_color=.false.

  end subroutine init_color
  !--------------------------------------------------------------
  function esc(code) result(s)
    character(len=*), intent(in) :: code
    character(len=:), allocatable :: s
    if(use_color)then
       s = achar(27)//'['//code//'m'
    else
       s = ''
    endif
  end function esc
  !--------------------------------------------------------------
  function c_bold() result(s)
    character(len=:), allocatable :: s
    s = esc('1')
  end function c_bold
  function c_dim() result(s)
    character(len=:), allocatable :: s
    s = esc('2')
  end function c_dim
  function c_under() result(s)
    character(len=:), allocatable :: s
    s = esc('4')
  end function c_under
  function c_cyan() result(s)
    character(len=:), allocatable :: s
    s = esc('36')
  end function c_cyan
  function c_green() result(s)
    character(len=:), allocatable :: s
    s = esc('32')
  end function c_green
  function c_yellow() result(s)
    character(len=:), allocatable :: s
    s = esc('33')
  end function c_yellow
  function c_reset() result(s)
    character(len=:), allocatable :: s
    s = esc('0')
  end function c_reset
  !--------------------------------------------------------------
  ! Milestone marker: UTF-8 "four pointed star" (U+2726) and a blank
  function c_mark() result(s)
    character(len=:), allocatable :: s
    if(use_color)then
       s = achar(226)//achar(156)//achar(166)//' '
    else
       s = ''
    endif
  end function c_mark
end module color_module
