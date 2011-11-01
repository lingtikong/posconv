!*******************************************************************************
!*** Subroutine to write the position file in format of MS Modeling:         ***
!*** Car and Res.                                                            ***
!*******************************************************************************
subroutine readcar
use prec
use cell
use iounits
implicit none
   !-----------------------------------------------------------------
   integer            :: ioerr, i
   real(q)            :: radum(3), rdum
   logical            :: lpbc
   character (len=256):: input
   character (len=20 ):: ctmp, stmp
   integer, external  :: typescreen
   !-----------------------------------------------------------------
   subname = 'readcar'
   !
   ntype = 0
   natom = 0
   !
   do while ( .true. )
      read( ioin, '(A)', iostat=ioerr ) input
      if ( ioerr.ne.0 ) exit
      read( input, *, iostat=ioerr ) stmp, radum, ctmp, idum, ctmp, ctmp, rdum
      if ( ioerr.eq.0 ) natom = natom + 1
   enddo
   if ( natom.lt.1 ) then
      info = 'Too few atoms read from: '//trim(infile)
      call error( subname, info, abs(natom-1) )
   endif
   allocate( atpos(3, natom), atrel(3,natom), attyp(natom) )
   atrel = 1
   !
   rewind( ioin )
   info = 'Error encountered while reading file'//trim(infile)
   !
   read( ioin, *, iostat=ioerr ) 
   read( ioin, '(A)', iostat=ioerr ) stmp
   call error( subname, info, ioerr )
   stmp = stmp(5:6)
   lpbc = stmp.eq.'ON'.or.stmp.eq.'on'.or.stmp.eq.'On'.or.stmp.eq.'oN'
   !
   read( ioin, '(A)', iostat=ioerr ) title
   read( ioin, *, iostat= ioerr )
   if ( lpbc ) then
      read( ioin, *, iostat=ioerr ) stmp, latt
   else
      latt(1:3) = 200.D0
      latt(4:6) = 90.D0
   endif
   call error( subname, info, ioerr )
   !
   ! Read atomic positions
   do i = 1, natom
      read( ioin, 100, iostat=ioerr ) stmp, atpos(:,i), stmp, idum, stmp, ETmp, rdum
      attyp(i) = typescreen( ETmp )
   enddo
100 format( A5, 3(1X, F14.9), 1X, A4, 1X, I1, 6X, A2, 6X, A2, 1X, F6.3 )
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
   cartesian = .true.
   call abc2axis()
   !
return
end subroutine

subroutine writecar
use prec
use cell
use iounits
use datetime
implicit none
   !----------------------------------------------------------------------------
   character(len=10  ) :: fmtres = '(A?,I?)', atsn
   integer             :: i, nd
   !----------------------------------------------------------------------------
   if ( .not.cartesian ) call dir2car
   write( ioout, 100 )
   write( ioout, 200 )
   write( ioout, 300 ) trim(title)
   write( ioout, 400 ) datestr(1:4), datestr(5:6), datestr(7:8), timestr(1:2), timestr(3:4), timestr(5:6)
   write( ioout, 500 ) latt
   !
   do i = 1, natom
      write( fmtres(3:3), 550 ) len_trim(EName(attyp(i)))
      nd = int( log10(dble(i)) ) + 1
      write( fmtres(6:6), 550 ) nd
      write( atsn, fmtres ) EName(attyp(i)), i
      write( ioout, 800 ) atsn, atpos(:,i), 1, EName(attyp(i)), 0.D0
   enddo
   write( ioout, 600 )
   !
100 format( "!BIOSYM archive 3" )
200 format( "PBC=ON" )
300 format( A )
400 format( "!Converted ", A4, "-", A2, "-", A2, 2X, A2,":",A2,":",A2 )
500 format( "PBC", 6F12.6, 1X, "(P1)" )
550 format( I1 )
600 format( "end", /, "end" )
800 format( A5, 3(1X, F14.9), 1X, "XXXX", 1X, I1, 6X, "xx", 6X, A2, 1X, F6.3 )
   !
return
end subroutine
!
subroutine writeres
use prec
use cell
use iounits
use datetime
implicit none
   !----------------------------------------------------------------------------
   character(len=10  ) :: fmtres = '(A?,I?)', atsn
   integer             :: i, j, k, ibrav = 14, nd
   !----------------------------------------------------------------------------
   !
   if ( cartesian ) call car2dir
   write( ioout, 100 ) trim(title)
   write( ioout, 200 ) 0.D0, latt
   write( ioout, 300 )
   write( ioout, 400 )
   do i = 1, ntype
      write( ioout, 500 ) EName(i)
   enddo
   write( ioout, * )
   do j = 1, ntype
      k = 0
      write( fmtres(3:3), 600 ) len_trim(EName(j))
      do i = 1, natom
         if ( attyp(i).ne.j ) cycle
         !
         k = k + 1
         nd = int( log10(dble(k)) ) + 1
         write( fmtres(6:6), 600 ) nd
         write( atsn, fmtres ) EName(j), k
         !
         write( ioout, 700 ) trim(atsn), j, atpos(:,i), 1.D0, 0.D0
         !
      enddo
      !
   enddo
   write( ioout, 800 )
   !
100 format( "TITL", 2X, A )
200 format( "CELL", 7F12.6 )
300 format( "LATT -1" )
400 format( "SFAC", $ )
500 format( A3, $     )
600 format( I1        )
700 format( A5, 1X, I3, 5(1X,F9.5) )
800 format( "END"     )
900 format( 10X, "Atomic configuration is written to: ", A )
   !
return
end subroutine
