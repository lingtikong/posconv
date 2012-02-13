!*******************************************************************************
!*** Subroutines to read/write XYZ format position file                      ***
!*******************************************************************************
subroutine read_abinit_xyz
use prec
use cell
use iounits
implicit none
   !-----------------------------------------------------------------
   integer            :: ioerr, i, j
   real(q)            :: radum(3), vol, factor
   character (len=256):: input, lunit
   character (len=20 ):: strtmp
   integer, external  :: typescreen
   !-----------------------------------------------------------------
   subname = 'read_abinit_xyz'
   !
   ntype = 0
   ! Read total number of atoms
   read( ioin, *, iostat=ioerr ) natom, lunit
   call error( subname, info, ioerr )
   if ( natom.lt.1 ) then
      info = 'Too few atoms read from: '//trim(infile)
      call error( subname, info, abs(natom-1) )
   endif
   factor = 1.D0
   if ((lunit.eq."atomic").or.(lunit.eq."atomicd0").or.(lunit.eq."Bohr") ) factor = 0.529177249D0

   allocate( atpos(3, natom), atrel(3,natom), attyp(natom) )
   atrel = 1
   !
   ! Read title of the coordinationfile
   read( ioin, '(A)', iostat=ioerr ) title
   call error( subname, info, ioerr )
   !
   read( title, *, iostat=ioerr) strtmp, a, b, c
   if ( b < 1.D0 ) b = 1.D0 ! Just in case this is the surface mode, where b is ignored for abinit
   !
   alpha = 0.D0
   beta  = 0.D0
   gamma = 0.D0
   !
   call abc2axis()
   !
   ! Read atomic positions
   do i = 1, natom
      read( ioin, '(A)', iostat=ioerr ) ETmp, radum
      call error( subname, info, ioerr )
      attyp(i)   = typescreen( ETmp )
      atpos(i,:) = radum * factor
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
   cartesian = .true.
   call axis2abc()
   !
return
end subroutine
!
subroutine write_abinit_xyz
use prec
use cell
use iounits
implicit none
   !----------------------------------------------------------------------------
   integer :: i, j
   !----------------------------------------------------------------------------
   !
   if ( .not.cartesian ) call dir2car
   !
   write( ioout, 100 ) natom
   write( ioout, 200 ) a, b, c
   !
   do i = 1, natom
      write( ioout, 300 ) EName( attyp(i) ), atpos(:, i)
   enddo
   !
100 format( I8, " angstroemd0" )
200 format( "periodic",3(1x,E24.17) )
300 format( A2, 3(1X, E24.17) )
   !
return
end subroutine
