!*******************************************************************************
!*** Subroutines to read/write XYZ format position file                      ***
!*******************************************************************************
subroutine readxyz
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
   subname = 'readxyz'
   !
   ntype = 0
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
   ! Read title of the coordinationfile
   read( ioin, '(A)', iostat=ioerr ) title
   call error( subname, info, ioerr )
   !
   ! Read atomic positions
   do i = 1, natom
      read( ioin, '(A)', iostat=ioerr ) input
      call error( subname, info, ioerr )
      idum = 0
      read(input, *, iostat=ioerr ) ETmp, atpos(:,i), strtmp, idum, radum
      if ( ioerr.eq.0.and.(idum.ge.1.and.idum.le.3) ) axis(idum, :) = radum
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
   alat  = 1.D0
   !
   call volume( axis, vol )
   if ( vol.lt.1.D-6 ) then
      !
      info  = 'Error encoutered while reading standard input!'
      write(*, '(/, 10x, "Please input the lattice constant of the cell:", $ )' )
      read (*, *, iostat=ioerr) alat
      call error( subname, info, ioerr )
      !
      write(*, '( 10x, "Please input the value for vector A: ", $ )' )
      read (*, *, iostat=ioerr) axis(1,:)
      call error( subname, info, ioerr )
      !
      write(*, '( 10x, "Please input the value for vector B: ", $ )' )
      read (*, *, iostat=ioerr) axis(2,:)
      call error( subname, info, ioerr )
      !
      write(*, '( 10x, "Please input the value for vector C: ", $ )' )
      read (*, *, iostat=ioerr) axis(3,:)
      call error( subname, info, ioerr )
      !
   endif
   cartesian = .true.
   call axis2abc()
   !
return
end subroutine
!
subroutine writexyz
use prec
use cell
use iounits
implicit none
   !----------------------------------------------------------------------------
   integer :: i, j
   character (len=512):: newtitle
   !----------------------------------------------------------------------------
   !
   if ( .not.cartesian ) call dir2car()
   !
   write( newtitle, 010)
   do i = 1, 3
      write( newtitle, 020) trim(newtitle), axis(i,:)*alat
   enddo
   write( newtitle, 030) trim(newtitle), trim(title)
   !
   write( ioout, 100 ) natom
   write( ioout, 200 ) trim(newtitle)
   !
   do i = 1, min(4,natom)
      write( ioout, 300 ) EName( attyp(i) ), atpos(:, i)
      select case ( i )
      case ( 1:3 )
         write( ioout, 320 ) i, axis(i,:)*alat
      case default
         write( ioout, * )
      end select
   enddo
   do i = 5, natom
      write( ioout, 350 ) EName( attyp(i) ), atpos(:, i)
   enddo
   !
   if (natom < 4) then
     write(*,'(/,10x,"Lattice vector info:")')
     do i = 1, 3
       write(*, 320) i, axis(i,:)*alat
     enddo
   endif
   !
010 format( "lattice='")
020 format( A, 3(1X, F15.10) )
030 format( A, "' title='", A, "'")
100 format( I8 )
200 format( A  )
300 format( A2, 3(1X, F20.15), $ )
320 format( " crystal_vector ", I2, 3(1x, F20.15) ) 
350 format( A2, 3(1X, F20.15)    )
   !
return
end subroutine
