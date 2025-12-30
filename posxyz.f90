!*******************************************************************************
!*** Subroutines to read/write XYZ format position file                      ***
!*******************************************************************************
subroutine readxyz
use prec
use cell
use iounits
implicit none
   !-----------------------------------------------------------------
   integer            :: ioerr, i, j, il, ir, flag = 0
   real(q)            :: radum(3), vol
   character (len=512):: input
   character (len=20 ):: strtmp
   character (len=2  ):: delimiter = "'"
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
   ! analyze the title for extended xyz format
   ! convert all letters in titile as lower-case
   do i = 1, len(title)
      if (title(i:i) >= 'A' .and. title(i:i) <= 'Z') then
         title(i:i) = char(ichar(title(i:i)) - ichar('A') + ichar('a'))
      endif
   enddo
   !
   i = index(title, 'lattice')
   if (i > 0) then
      il = index(title(i+1:), delimiter)
      if (il == 0) then
         delimiter = '"'
         il = index(title(i+1:), delimiter)
      endif
      il = il + i + 1
      ir = index(title(il:), delimiter) + il - 2
      read(title(il:ir), *, iostat = ioerr) axis
      if (ioerr == 0) then
         flag = 1
      else
         axis = 0.
      endif
   endif
   !
   ! Read atomic positions
   do i = 1, natom
      read( ioin, '(A)', iostat=ioerr ) input
      call error( subname, info, ioerr )
      idum = 0
      read(input, *, iostat=ioerr ) ETmp, atpos(:,i), strtmp, idum, radum
      if ( flag.eq.0.and.ioerr.eq.0.and.(idum.ge.1.and.idum.le.3) ) axis(idum, :) = radum
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
   character (len=512):: newtitle ! extended xyz format
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
030 format( A, "' pbc='T T T' properties=species:S:1:pos:R:3 title='", A, "'")
100 format( I8 )
200 format( A  )
300 format( A2, 3(1X, F20.15), $ )
320 format( " crystal_vector ", I2, 3(1x, F20.15) ) 
350 format( A2, 3(1X, F20.15)    )
   !
return
end subroutine
