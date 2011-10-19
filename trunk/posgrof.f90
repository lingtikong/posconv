!*******************************************************************************
!*** Subroutines to read/write the position file in groF format.             ***
!*******************************************************************************
subroutine readgroF
use prec
use cell
use iounits
implicit none
   !----------------------------------------------------------------------------
   integer            :: ioerr, i, j
   character (len=256):: input
   integer, external  :: typescreen
   !----------------------------------------------------------------------------
   subname = 'readgroF'
   !
   ntype = 0
   !
   ! Read total number of atoms
   read(ioin, *, iostat=ioerr) natom
   call error( subname, info, ioerr )
   !
   if ( natom.lt.1 ) then
      info = 'Too few atoms read from: '//trim(infile)
      call error( subname, info, abs(natom-1) )
   endif
   allocate( atpos(3, natom), atrel(3,natom), attyp(natom) )
   !
   do i = 1, natom
      read(ioin, '(A)', iostat=ioerr ) input
      call error( subname, info, ioerr )
      input = trim(input)//" 1 1 1"
      read(input, *, iostat=ioerr) idum, atpos(:,i), ETmp, atrel(:,i)
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
   !
   input = info
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
   write(*, '(/, 10x, "As identified, there are ", I2, " kinds of atoms in the cell," )' ) ntype
   write(*, '(   10x, "please input their respective element name(s) :", $ )' )
   read (*, *, iostat=ioerr) EName
   call error( subname, info, ioerr )
   !
   info = input
   cartesian = .false.
   call axis2abc
   !
return
end subroutine
!
subroutine writegroF
use prec
use cell
use iounits
implicit none
   !----------------------------------------------------------------------------
   integer :: i, j, ibrav = 14
   !----------------------------------------------------------------------------
   !
   if ( cartesian ) call car2dir
   write( ioout, 100 ) natom, natom
   do i = 1, natom
      write( ioout, 200 ) i, atpos(:, i), attyp(i), atrel(:,i)
   enddo
   !
   write(*, 300 )
   write(*, 310 )
   write(*, 320 )
   write(*, 330 ) a * 1.0D-10
   write(*, 340 ) b * 1.0D-10
   write(*, 350 ) c * 1.0D-10
   write(*, 360 )
   write(*, 370 ) CHAR(34), trim(outfile), CHAR(34)
   write(*, 380 )
   !
   write(*, 400 )
   write(*, 320 )
   write(*, 410 )
   write(*, 420 ) ibrav 
   write(*, 430 ) natom
   write(*, 440 ) alat
   write(*, 360 )
   write(*, 450 ) CHAR(34), trim(outfile), CHAR(34)
   write(*, 460 ) "/"
   write(*, 360 ) "Cell Parameters"
   write(*, 470 ) axis(1, :)
   write(*, 470 ) axis(2, :)
   write(*, 470 ) axis(3, :)
   write(*, 480 )
   write(*, 380 )
   !
   write(*, 500 )
   write(*, 320 )
   write(*, 410 )
   write(*, 430 ) natom
   write(*, 440 ) alat
   write(*, 510 ) char(34), trim(outfile), char(34)
   write(*, 460 ) "/"
   write(*, 360 ) "Cell Parameters"
   write(*, 470 ) axis(1, :)
   write(*, 470 ) axis(2, :)
   write(*, 470 ) axis(3, :)
   write(*, 380 )
   !
100 format( I6, 2X, I6 )
200 format( I6, 2X, 3(F20.16,2X), 4(I2, 2X) )
300 format( /, 80("=") )
310 format( "Remember to modify the following information for groF:" )
320 format( 60("-") )
330 format( "setenv boxlx  ", E20.10 )
340 format( "setenv boxly  ", E20.10 )
350 format( "setenv boxlz  ", E20.10 )
360 format( A )
370 format( "setenv cnfile ", A1, A, A1 )
380 format( 80("=") )
400 format( /, "Remember to modify the following information for pceam:" )
410 format( 5X, "&INPUT" )
420 format( 5X, "ibrav  =", I3 )
430 format( 5X, "natom  =", I6 )
440 format( 5X, "alat   =", F12.6 )
450 format( 5X, "CONFIG =", A1, A, A1 )
460 format( 5X, A )
470 format( 3(F20.16,2X) )
480 format( /, "ibrav should be correct according to your lattice." )
500 format( /, "Remember to modify the following information for MSS:" )
510 format( 5X, "poscar=", 3A )
   !
return
end subroutine
