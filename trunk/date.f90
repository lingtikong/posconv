! Module to call date_and_time
module datetime
implicit none
   !----------------------------------------------------------
   character (len=10) :: datestr, timestr
   integer            :: datevalue(8)
   !----------------------------------------------------------
   !
contains
   !
   subroutine date_time( )
   implicit none
   call date_and_time( date=datestr, time=timestr, values=datevalue )
   return
   end subroutine date_time
end module
