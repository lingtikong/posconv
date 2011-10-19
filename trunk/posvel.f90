!*******************************************************************************
!*** Subroutines to read/write the position file in format of ReaxFF save.   ***
!*******************************************************************************
subroutine readvel
use prec
use cell
use iounits
implicit none
   !----------------------------------------------------------------------------
   integer :: ioerr, i, j, iatom
   integer, external :: typescreen
   !----------------------------------------------------------------------------
   subname = 'readvel'
   !
   ntype = 0
   !
   ! Read lattice parameters:
   read( ioin, *, iostat=ioerr)
   call error( subname, info, ioerr )
   read( ioin, *, iostat=ioerr) latt(1:3)
   read( ioin, *, iostat=ioerr) latt(4:6)
   ! Read total number of atoms
   read( ioin, *, iostat=ioerr ) natom
   call error( subname, info, ioerr )
   !
   if ( natom.lt.1 ) then
      info = 'Too few atoms read from: '//trim(infile)
      call error( subname, info, abs(natom-1) )
   endif
   allocate( atpos(3, natom), atrel(3,natom), attyp(natom) )
   atrel = 1
   !
   title = 'Atomic position read from:'//trim(infile)
   !
   ! Read atomic positions
   do i = 1, natom
      read( ioin, *, iostat=ioerr ) atpos(:,i), ETmp
      call error( subname, info, ioerr )
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
   forall ( i=1:ntype ) ntm(i) = count( attyp.eq.i )
   EName = Eread(1:ntype)
   !
   call abc2axis
   cartesian = .true.
   !
return
end subroutine
!
subroutine writevel
use prec
use cell
use iounits
implicit none
   !----------------------------------------------------------------------------
   integer :: i
   !----------------------------------------------------------------------------
   !
   if ( .not.cartesian ) call dir2car
   write( ioout, 100 )
   write( ioout, 200 ) latt(1:3)
   write( ioout, 200 ) latt(4:6)
   write( ioout, 300 ) natom
   !
   do i = 1, natom
      write( ioout, 400 ) atpos(:, i), EName( attyp(i) )
   enddo
   !
100 format( 1X, "Lattice parameters:" )
200 format( 3(1X,F14.8)  )
300 format( I5, 1X, "Atom coordinates (Angstrom):" )
400 format( 3(1X, E23.15), 1X, A2 )
   !
return
end subroutine
