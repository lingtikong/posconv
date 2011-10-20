!******************************************************************
!***  Module of prec, to adapt double precision to all machine. ***
!******************************************************************
module prec
implicit none
   integer, parameter :: q = SELECTED_REAL_KIND(10)
   !
contains
   !
   function is_numeric(string)
   implicit none
      character(len=*), intent(in) :: string
      logical :: is_numeric
      real    :: x
      integer :: e
      read(string,*,iostat=e) x
      is_numeric = e == 0
   end function is_numeric
end module

module iounits
implicit none
   integer, parameter :: ioin  = 8
   integer, parameter :: ioout = 9
   integer, parameter :: iotmp = 10
end module
