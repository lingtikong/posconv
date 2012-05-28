!*******************************************************************************
!*** Subroutines to read/write PWSCF position card.                          ***
!*******************************************************************************
subroutine readpwscf
use prec
use cell
use iounits
implicit none
   !----------------------------------------------------------------------------
   info = 'Sorry, I can write but cannot read this kind of coordination file!'
   call error( subname, info, 2 )
   !
end subroutine
!
subroutine writepwscf
use prec
use cell
use iounits
implicit none
   !----------------------------------------------------------------------------
   integer :: i, j, ibrav = 14
   !----------------------------------------------------------------------------
   !
   if ( cartesian ) call car2dir
   !
   write( ioout, 90 )
   !
   do i = 1, natom
      j = attyp(i)
      write( ioout, 100 ) Ename(j), atpos(:, i)
   enddo
   !
   write( ioout, 101 )
   write( ioout, 200 ) "&system"
   write( ioout, 300 ) "ibrav     = ", ibrav
   write( ioout, 400 ) "celldm(1) = ", a * 1.889726D0
   write( ioout, 400 ) "celldm(2) = ", b/a
   write( ioout, 400 ) "celldm(3) = ", c/a
   write( ioout, 400 ) "celldm(4) = ", cosbc
   write( ioout, 400 ) "celldm(5) = ", cosac
   write( ioout, 400 ) "celldm(6) = ", cosab
   write( ioout, 500 ) "nat       = ", natom
   write( ioout, 500 ) "ntyp      = ", ntype
   write( ioout, 600 )
   !
 90 format( "ATOMIC_POSITIONS  {crystal}" )
100 format( A2, 2X, 3(F20.16, 2X) )
101 format( // )
200 format( A  )
300 format( 2X, A, I3 )
400 format( 2X, A, F16.12 )
500 format( 2X, A, I5 )
600 format( /, 2X, "Note: ibrav might be wrong!" )
   !
return
end subroutine
