!******************************************************************************
!***    Subroutine to read VASP format position file                        ***
!******************************************************************************
subroutine readvasp
use prec
use cell
use iounits
implicit none
   !----------------------------------------------------------------------------
   integer             :: ioerr, i, j, iatom, istr, iend
   real(q)             :: radum(3), vol
   character (len=256) :: input, elements
   character (len=20 ) :: strtmp, flag
   character (len=1 )  :: cdum, cadum(3)
   logical             :: seldyn = .false., f_eread = .false.
   !----------------------------------------------------------------------------
   subname = 'readvasp'
   !
   seldyn    = .false.
   ! Read title
   read( ioin, '(A)', iostat=ioerr ) title
   call error( subname, info, ioerr )
   ! Read lattice constant
   read( ioin, *, iostat=ioerr ) alat
   call error( subname, info, ioerr )
   ! Read cell vectors
   read( ioin, *, iostat=ioerr ) axis(1, :)
   call error( subname, info, ioerr )
   read( ioin, *, iostat=ioerr ) axis(2, :)
   call error( subname, info, ioerr )
   read( ioin, *, iostat=ioerr ) axis(3, :)
   call error( subname, info, ioerr )
   ! Read number of atoms
   read( ioin, '(A)', iostat=ioerr ) input
   call error( subname, info, ioerr )
   if ( .not.is_numeric(input) ) then
      elements = input  ! vasp 5 writes elements info before defining # of atoms
      f_eread  = .true.
      read( ioin, '(A)', iostat=ioerr ) input
      call error( subname, info, ioerr )
   endif
   ! If NMax is changed in module 'cell', the number of '0' here should
   ! also be changed.
   input = trim(input)//" 0 0 0 0 0 0 0 0 0 0"
   read(input, *, iostat=ioerr ) iadum
   call error( subname, info, ioerr )
   natom = sum( iadum )
   ntype = count( iadum.gt.0 )
   !
   if ( natom.lt.1 ) then
      info = 'Total number of atoms read from:'//trim(infile)// &
           & ' is less than 1!'
      call error( subname, info, (natom-1) )
   endif
   if ( ntype.lt.1 ) then
      info = 'Total type of atoms read from:'//trim(infile)// &
           & ' is less than 1!'
      call error( subname, info, (ntype-1) )
   endif
   !
   ! Initialize allocatable arrays
   allocate( atpos(3, natom), atrel(3, natom), attyp(natom), ntm( ntype ), EName( ntype ) )
   atrel = 1
   read( input, *, iostat=ioerr ) ntm
   call error( subname, info, ioerr )
   !
   ! Continue to read POSCAR
   ! To read the tag for selective dynamics and/or coordination type
   read( ioin, '(A)', iostat=ioerr ) input
   call error( subname, info, ioerr )
   read(input, *, iostat=ioerr ) cdum
   call error( subname, info, ioerr )
   if ( cdum.eq.'s'.or.cdum.eq.'S' ) then
      seldyn = .true.
      read( ioin, '(A)', iostat=ioerr ) input
      call error( subname, info, ioerr )
      read(input, *, iostat=ioerr ) cdum
      call error( subname, info, ioerr )
   endif
   !
   cartesian = ( (cdum.eq.'C').or.(cdum.eq.'c').or.(cdum.eq.'K').or.(cdum.eq.'k') )
   !
   j = 1
   ! To read atomic positions
   if ( seldyn ) then
      do i = 1, natom
         read( ioin, *, iostat=ioerr) atpos(:,i), cadum
         where ( cadum.eq.'F'.or.cadum.eq.'f' ) atrel(:,i) = 0
      enddo
   else
      do i = 1, natom
         read( ioin, *, iostat=ioerr) atpos(:,i)
      enddo
   endif
   do i = 1, ntype
      istr = sum(ntm(1:i-1))+1
      iend = sum(ntm(1:i))
      attyp(istr:iend) = i
   enddo
   !
   ! To assign all other cell related variables
   if ( alat.lt.0.D0 ) then
      call volume( axis, vol )
      alat = ( -alat / vol ) ** ( 1.D0/3.D0 )
   endif
   !
   call axis2abc()
   !
   ! Ask the user to input the element(s) name
   info = 'Wrong input while the system tries to get the element names!'
   if (f_eread) then
      read(elements, *, iostat=ioerr) EName
   else
      write(*, '(/, 10x, "Please input the element name for the ", I2, " kind(s) of atoms" )' ) ntype
      write(*, '(   10x, "in the sequence of POSCAR : ", $ )' )
      read (*, *, iostat=ioerr ) EName
   endif
   call error( subname, info, ioerr )
   ERead(1:ntype) = EName
   !
return
end subroutine
!
subroutine writevasp
use prec
use cell
use iounits
implicit none
   !----------------------------------------------------------------------------
   character(len=40  ) :: fmtstr1 = '(??(A6,2X) )', fmtstr2 = '(??(I6, 1X) )'
   integer             :: i, j
   logical             :: seldyn
   !----------------------------------------------------------------------------
   if ( cartesian ) call car2dir
   call reduceaxis
   !
   write( fmtstr1(2:3), 100 ) ntype
   write( fmtstr2(2:3), 100 ) ntype
   write( ioout, 110 ) trim(title)
   write( ioout, 120 ) alat
   write( ioout, 150 ) axis(1, :)
   write( ioout, 150 ) axis(2, :)
   write( ioout, 150 ) axis(3, :)
   write( ioout, fmtstr1 ) EName
   write( ioout, fmtstr2 ) ntm
   !
   seldyn = sum(atrel).lt.(3*natom)
   if ( seldyn ) write( ioout, 160 )
   write( ioout, 180 )
   !
   if ( seldyn ) then
      do j = 1, ntype
      do i = 1, natom
         if ( attyp(i).eq.typeID(j) ) write( ioout, 200 ) atpos(:, i), atrel(:,i)
      enddo
      enddo
   else
      do j = 1, ntype
      do i = 1, natom
         if ( attyp(i).eq.typeID(j) ) write( ioout, 250 ) atpos(:, i)
      enddo
      enddo
   endif
   !
100 format( I2 )
110 format( A  )
120 format( F20.14 )
150 format( 3(F20.14,2X) )
160 format( "Selective dynamics" )
180 format( "direct" )
200 format( 3( F20.14, 2X ), 3L2 )
250 format( 3( F20.14, 2X ) )
return
end subroutine
