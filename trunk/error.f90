!***********************************************************************
!***            Subroutine for warning and error                     ***
!***********************************************************************
subroutine warn( sub, info, ierr)
implicit none
   !--------------------------------------------------------------------
   character(len=10) :: sub
   character(len=80) :: info
   integer           :: ierr
   !--------------------------------------------------------------------
   write(*, '(/, 10X, 60("W") )')
   write(*, '( 10X, "WWW", 14X, "WARNING from ", A12, 15X, "WWW")') sub
   write(*, '( 10X, "WWW", 54X, "WWW")')
   write(*, '( 10X, "WWW", 2X, A50, 2X, "WWW" )'   ) trim( info(1:50) )
   if ( trim( info(51:) ).ne.'' ) &
   &  write(*, '( 10X, "WWW", 2X, A30, 22X, "WWW" )'  ) trim( info(51:) )
   write(*, '( 10X, "WWW", 2X, "Error Code: ", I4, 36X, "WWW" )' ) ierr
   write(*, '( 10X, "WWW", 54X, "WWW", /, 60("W") )')
   !--------------------------------------------------------------------
return
end subroutine
!
subroutine error( sub, info, ierr )
implicit none
   !--------------------------------------------------------------------
   character(len=10) :: sub
   character(len=80) :: info
   integer           :: ierr
   !--------------------------------------------------------------------
   if ( ierr.eq.0 ) return
   !
   write(*, '(/, 10X, 60("E") )')
   !
   write(*, '( 10X, "EEE", 16X, "ERROR from ", A12, 15X, "EEE")') sub
   write(*, '( 10X, "EEE", 54X, "EEE")')
   write(*, '( 10X, "EEE", 2X, A50, 2X, "EEE" )'   ) trim( info(1:50) )
   if ( trim( info(51:) ).ne.'' ) &
   &  write(*, '( 10X, "EEE", 2X, A30, 22X, "EEE" )'  ) trim( info(51:) )
   write(*, '( 10X, "EEE", 2X, "Error Code: ", I4, 36X, "EEE" )' ) ierr
   write(*, '( 10X, "EEE", 54X, "EEE", /, 10X, 60("E") )')
   !
   if ( ierr.gt.0 ) stop
   !---------------------------------------------------------------------
return
end subroutine
