!*******************************************************************************
!*** Subroutines to read/write siesta STRUCT_IN format position file         ***
!*******************************************************************************
subroutine readsiesta
use prec
use cell
use iounits
implicit none
   !-----------------------------------------------------------------
   integer            :: ioerr, i, j
   real(q)            :: radum(3), vol
   character (len=256):: input
   character (len=20 ):: strtmp
   integer, external  :: typescreen
   !-----------------------------------------------------------------
   subname = 'readsiesta'
   !
   alat = 1.D0
   ! read lattice vectors
   read( ioin, *, iostat=ioerr ) (axis(:,i), i=1,3)
   ! Read total number of atoms
   read( ioin, *, iostat=ioerr ) natom
   call error( subname, info, ioerr )
   if ( natom.lt.1 ) then
      info = 'Too few atoms read from: '//trim(infile)
      call error( subname, info, abs(natom-1) )
   endif
   allocate( atpos(3, natom), atrel(3,natom), attyp(natom) )
   atrel = 1
   !
   title = "Siesta STRUCT_IN info from"//trim(infile)
   !
   ntype = 0
   ! Read atomic positions
   do i = 1, natom
      read( ioin, '(A)', iostat=ioerr ) input
      call error( subname, info, ioerr )
      read(input, *, iostat=ioerr ) ETmp, idum, atpos(:,i)
      attyp(i) = typescreen( ETmp )
   enddo
   !
   if ( ntype.lt.1 ) then
      info = 'Too few elements identified from:'//trim(infile)
      call error( subname, info, abs(ntype-1) )
   endif
   !
   ! To assign all other cell related variables
   allocate( ntm(ntype), EName(ntype) )
   forall( i=1:ntype ) ntm(i) = count(attyp.eq.i)
   EName = Eread(1:ntype)
   !
   call volume( axis, vol )
   cartesian = .false.
   call axis2abc()
   !
return
end subroutine
!
subroutine writesiesta
use prec
use cell
use iounits
implicit none
   !----------------------------------------------------------------------------
   integer :: i, j
   !----------------------------------------------------------------------------
   !
   if ( cartesian ) call dir2car
   !
   write( ioout, 100 ) (axis(i,:), i=1, 3)
   write( ioout, 200 ) natom
   !
   do i = 1, natom
      write( ioout, 300 ) attyp(i), i, atpos(:, i)
   enddo
   !
100 format( 9(1x,F20.15) )
200 format( I8  )
300 format( I3, I6, 3(1X, F20.12) )
   !
return
end subroutine
