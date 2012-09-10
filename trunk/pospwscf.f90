!*******************************************************************************
!*** Subroutines to read/write PWSCF position card.                          ***
!*******************************************************************************
subroutine readpwscf
use prec
use cell
use iounits
implicit none
   !----------------------------------------------------------------------------
   character (len=512):: input
   !
   integer :: ibrav, i, j, ptyp, ioerr
   real(q) :: celldm(6), vol
   real(q) :: u, v, w, tx, ty, tz
   !----------------------------------------------------------------------------
   real(q), parameter :: bohr2ang = 0.5291772108D0
   integer, external  :: typescreen
   !----------------------------------------------------------------------------
   subname = 'readpwscf'
   info = 'Error while reading stdin.'
   !
   ntype = 0
   ! ask for total number of atoms
   write(*,'(/,10x,"Please input the total number of atoms: ", $)')
   read(*, *, iostat=ioerr ) natom
   call error( subname, info, ioerr )
   if ( natom.lt.1 ) then
      info = 'Too few atoms read from input.'
      call error( subname, info, abs(natom-1) )
   endif
   allocate( atpos(3, natom), atrel(3,natom), attyp(natom) )
   atrel = 1
   !
   ! To get lattice and lattice constant info
   write(*,'(10x,"Please input the lattice type (ibrav) : ", $)')
   read(*, *, iostat=ioerr ) ibrav
   call error( subname, info, ioerr )
   !
   write(*,'(10x,"Please input the alat of the lattice, in bohr: ", $)')
   read(*, *, iostat=ioerr ) alat
   call error( subname, info, ioerr )
   alat = alat * bohr2ang
   !
   do i = 1, 3
   do j = 1, 3
      axis(i,j) = 0.D0
   enddo
   enddo
   !
   ! To get lattice basis vector info
   select case (ibrav)
   case (0) ! Arbitray lattice
      !
      write(*,'(10x,"Please input the basis vectors of your lattice: ")')
      do i = 1, 3
         read(*, *, iostat=ioerr) axis(i,:)
         call error( subname, info, ioerr )
      enddo
      !
   case (1) ! SC
      !
      axis(1,1) = 1.D0
      axis(2,2) = 1.D0
      axis(3,3) = 1.D0
      !
   case (2) ! FCC
      !
      axis(1,1) = -0.5D0
      axis(1,3) =  0.5D0
      axis(2,2) =  0.5D0
      axis(2,3) =  0.5D0
      axis(3,1) = -0.5D0
      axis(3,2) =  0.5D0
      !
   case (3) ! BCC
      !
      axis(1,1) =  0.5D0
      axis(1,2) =  0.5D0
      axis(1,3) =  0.5D0
      axis(2,1) = -0.5D0
      axis(2,2) =  0.5D0
      axis(2,3) =  0.5D0
      axis(3,1) = -0.5D0
      axis(3,2) = -0.5D0
      axis(3,3) =  0.5D0
      !
   case (4) ! Hexagonal and Trigonal P
      !
      write(*,'(10x,"Please input the celldm(3) (c/a) of the lattice: ", $)')
      read(*, *, iostat=ioerr) celldm(3)
      call error( subname, info, ioerr )
      !
      axis(1,1) =  1.0D0
      axis(2,1) = -0.5D0
      axis(2,2) = sqrt(3.D0)*0.5D0
      axis(3,3) = celldm(3)
      !
   case (5) ! Trigonal R, 3 fold axis c
      !
      write(*,'(10x,"Please input the celldm(4) (cosAlpha) of the lattice: ", $)')
      read(*, *, iostat=ioerr) celldm(4)
      call error( subname, info, ioerr )
      !
      tx = sqrt((1.D0-celldm(4))*0.5D0)
      ty = sqrt((1.D0-celldm(4))/6.0D0)
      tz = sqrt((1.D0+2.D0*celldm(4))/3.0D0)
      !
      axis(1,1) =  tx
      axis(1,2) = -ty
      axis(1,3) =  tz
      axis(2,2) =  ty + ty
      axis(2,3) =  tz
      axis(3,1) = -tx
      axis(3,2) = -ty
      axis(3,3) =  tz
      !
   case (-5) ! Trigonal R, 3 fold axis <111>
      !
      write(*,'(10x,"Please input the celldm(4) (cosAlpha) of the lattice: ", $)')
      read(*, *, iostat=ioerr) celldm(4)
      call error( subname, info, ioerr )
      !
      tx = sqrt((1.D0-celldm(4))*0.5D0)
      ty = sqrt((1.D0-celldm(4))/6.0D0)
      tz = sqrt((1.D0+2.D0*celldm(4))/3.0D0)
      u = (tz - 2.D0*sqrt(2.D0)*ty) / sqrt(3.D0)
      v = (tz + sqrt(2.D0)*ty) / sqrt(3.D0)
      !
      axis(1,1) =  u
      axis(1,2) =  v
      axis(1,3) =  v
      axis(2,1) =  v
      axis(2,2) =  u
      axis(2,3) =  v
      axis(3,1) =  v
      axis(3,2) =  v
      axis(3,3) =  u
      !
   case (6) ! Tetragonal P (st)
      !
      write(*,'(10x,"Please input the celldm(3) (c/a) of the lattice: ", $)')
      read(*, *, iostat=ioerr) celldm(3)
      call error( subname, info, ioerr )
      !
      axis(1,1) = 1.D0 / sqrt(3.D0)
      axis(2,2) = 1.D0 / sqrt(3.D0)
      axis(3,3) = celldm(3) / sqrt(3.D0)
      !
   case (7) ! Tetragonal I (bct)
      !
      write(*,'(10x,"Please input the celldm(3) (c/a) of the lattice: ", $)')
      read(*, *, iostat=ioerr) celldm(3)
      call error( subname, info, ioerr )
      !
      axis(1,1) =  0.5D0
      axis(1,2) = -0.5D0
      axis(1,3) =  celldm(3)*0.5D0
      axis(2,1) =  0.5D0
      axis(2,2) =  0.5D0
      axis(2,3) =  celldm(3)*0.5D0
      axis(3,1) = -0.5D0
      axis(3,2) = -0.5D0
      axis(3,3) =  celldm(3)*0.5D0
      !
   case (8) ! Orthorhombic P
      !
      write(*,'(10x,"Please input the celldm(2) (b/a) of the lattice: ", $)')
      read(*, *, iostat=ioerr) celldm(2)
      call error( subname, info, ioerr )
      write(*,'(10x,"Please input the celldm(3) (c/a) of the lattice: ", $)')
      read(*, *, iostat=ioerr) celldm(3)
      call error( subname, info, ioerr )
      !
      axis(1,1) = 1.D0
      axis(2,2) = celldm(2)
      axis(3,3) = celldm(3)
      !
   case (9) ! Orthorhombic base-centered (bco)
      !
      write(*,'(10x,"Please input the celldm(2) (b/a) of the lattice: ", $)')
      read(*, *, iostat=ioerr) celldm(2)
      call error( subname, info, ioerr )
      write(*,'(10x,"Please input the celldm(3) (c/a) of the lattice: ", $)')
      read(*, *, iostat=ioerr) celldm(3)
      call error( subname, info, ioerr )
      !
      axis(1,1) = 0.5D0
      axis(1,2) = celldm(2)*0.5D0
      axis(2,1) = -0.5D0
      axis(2,2) = celldm(2)*0.5D0
      axis(3,3) = celldm(3)
      !
   case (-9) ! Orthorhombic base-centered (bco), alternative
      !
      write(*,'(10x,"Please input the celldm(2) (b/a) of the lattice: ", $)')
      read(*, *, iostat=ioerr) celldm(2)
      call error( subname, info, ioerr )
      write(*,'(10x,"Please input the celldm(3) (c/a) of the lattice: ", $)')
      read(*, *, iostat=ioerr) celldm(3)
      call error( subname, info, ioerr )
      !
      axis(1,1) = 0.5D0
      axis(1,2) = -celldm(2)*0.5D0
      axis(2,1) = 0.5D0
      axis(2,2) = celldm(2)*0.5D0
      axis(3,3) = celldm(3)
      !
   case (10) ! Orthorhombic face-centered (fco)
      !
      write(*,'(10x,"Please input the celldm(2) (b/a) of the lattice: ", $)')
      read(*, *, iostat=ioerr) celldm(2)
      call error( subname, info, ioerr )
      write(*,'(10x,"Please input the celldm(3) (c/a) of the lattice: ", $)')
      read(*, *, iostat=ioerr) celldm(3)
      call error( subname, info, ioerr )
      !
      axis(1,1) = 0.5D0
      axis(1,3) = celldm(3)*0.5D0
      axis(2,1) = 0.5D0
      axis(2,2) = celldm(2)*0.5D0
      axis(3,2) = celldm(2)*0.5D0
      axis(3,3) = celldm(3)*0.5D0
      !
   case (11) ! Orthorhombic body-centered (bco)
      !
      write(*,'(10x,"Please input the celldm(2) (b/a) of the lattice: ", $)')
      read(*, *, iostat=ioerr) celldm(2)
      call error( subname, info, ioerr )
      write(*,'(10x,"Please input the celldm(3) (c/a) of the lattice: ", $)')
      read(*, *, iostat=ioerr) celldm(3)
      call error( subname, info, ioerr )
      !
      axis(1,1) = 0.5D0
      axis(1,2) = celldm(2)*0.5D0
      axis(1,3) = celldm(3)*0.5D0
      axis(2,1) = -0.5D0
      axis(2,2) = celldm(2)*0.5D0
      axis(2,3) = celldm(3)*0.5D0
      axis(3,1) = -0.5D0
      axis(3,2) = -celldm(2)*0.5D0
      axis(3,3) =  celldm(3)*0.5D0
      !
   case (12) ! Monoclinic P, unique axis c
      !
      write(*,'(10x,"Please input the celldm(2) (b/a) of the lattice: ", $)')
      read(*, *, iostat=ioerr) celldm(2)
      call error( subname, info, ioerr )
      write(*,'(10x,"Please input the celldm(3) (c/a) of the lattice: ", $)')
      read(*, *, iostat=ioerr) celldm(3)
      call error( subname, info, ioerr )
      write(*,'(10x,"Please input the celldm(4) (cosAB) of the lattice: ", $)')
      read(*, *, iostat=ioerr) celldm(4)
      call error( subname, info, ioerr )
      !
      axis(1,1) = 1.D0
      axis(2,1) = celldm(2)*celldm(4)
      axis(2,2) = celldm(2)*sqrt(1.D0-celldm(4)*celldm(4))
      axis(3,3) = celldm(3)
      !
   case (-12) ! Monoclinic P, unique axis b
      !
      write(*,'(10x,"Please input the celldm(2) (b/a) of the lattice: ", $)')
      read(*, *, iostat=ioerr) celldm(2)
      call error( subname, info, ioerr )
      write(*,'(10x,"Please input the celldm(3) (c/a) of the lattice: ", $)')
      read(*, *, iostat=ioerr) celldm(3)
      call error( subname, info, ioerr )
      write(*,'(10x,"Please input the celldm(5) (cosAC) of the lattice: ", $)')
      read(*, *, iostat=ioerr) celldm(5)
      call error( subname, info, ioerr )
      !
      axis(1,1) = 1.D0
      axis(2,2) = celldm(2)
      axis(3,1) = celldm(3)*sqrt(1.D0-celldm(5)*celldm(5))
      axis(3,3) = celldm(3)*celldm(5)
      !
   case (13) ! Monoclinic base-centered
      !
      write(*,'(10x,"Please input the celldm(2) (b/a) of the lattice: ", $)')
      read(*, *, iostat=ioerr) celldm(2)
      call error( subname, info, ioerr )
      write(*,'(10x,"Please input the celldm(3) (c/a) of the lattice: ", $)')
      read(*, *, iostat=ioerr) celldm(3)
      call error( subname, info, ioerr )
      write(*,'(10x,"Please input the celldm(4) (cosAB) of the lattice: ", $)')
      read(*, *, iostat=ioerr) celldm(4)
      call error( subname, info, ioerr )
      !
      axis(1,1) = 0.5D0
      axis(1,3) = -0.5D0*celldm(3)
      axis(2,1) = celldm(2)*celldm(4)
      axis(2,2) = celldm(2)*sqrt(1.D0-celldm(4)*celldm(4))
      axis(3,1) = 0.5D0
      axis(3,3) = celldm(3)*0.5D0
      !
   case (14) ! Triclinic
      !
      write(*,'(10x,"Please input the celldm(2) (b/a) of the lattice: ", $)')
      read(*, *, iostat=ioerr) celldm(2)
      call error( subname, info, ioerr )
      write(*,'(10x,"Please input the celldm(3) (c/a) of the lattice: ", $)')
      read(*, *, iostat=ioerr) celldm(3)
      call error( subname, info, ioerr )
      write(*,'(10x,"Please input the celldm(4) (cosBC) of the lattice: ", $)')
      read(*, *, iostat=ioerr) celldm(4)
      call error( subname, info, ioerr )
      write(*,'(10x,"Please input the celldm(5) (cosAC) of the lattice: ", $)')
      read(*, *, iostat=ioerr) celldm(5)
      call error( subname, info, ioerr )
      write(*,'(10x,"Please input the celldm(6) (cosAB) of the lattice: ", $)')
      read(*, *, iostat=ioerr) celldm(6)
      call error( subname, info, ioerr )
      !
      axis(1,1) = 1.D0
      axis(2,1) = celldm(2)*celldm(6)
      axis(2,2) = celldm(2)*sqrt(1.D0-celldm(6)*celldm(6))
      axis(3,1) = celldm(3)*celldm(5)
      axis(3,2) = celldm(3)*(celldm(4)-celldm(5)*celldm(6))/sqrt(1.D0-celldm(6)*celldm(6))
      axis(3,3) = celldm(3)*sqrt(1.D0 + 2.D0*celldm(4)*celldm(5)*celldm(6) - celldm(4)*celldm(4) &
      &         - celldm(5)*celldm(5)-celldm(6)*celldm(6))/sqrt(1.D0-celldm(6)*celldm(6))
   case default
      info = 'Sorry, Unknown ibrav!'
      call error( subname, info, 2 )
   end select
   !
   ! Define title of the coordinationfile
   title = 'PWSCF type from'//infile
   !
   ! Read first line, analyse position format
   read( ioin, '(A)', iostat=ioerr ) input
   j = len_trim(input)
   if ( j < 17 ) then
      info = 'Error while reading '//trim(infile)
      call error(subname, info, 16-j)
   endif
   do i = 17, j
      if (ichar(input(i:i)).lt.ichar("a")) input(i:i) = char( ichar(input(i:i)) - ichar('A') + ichar('a'))
   enddo
   !
   ptyp = 1
   if (matches('bohr', input(17:)) ) then
      ptyp = 2
   else if (matches('crystal', input(17:)) )  then
      ptyp = 3
   else if (matches('angstrom', input(17:)) )  then
      ptyp = 4
   endif
   !
   info = 'Error while reading atomic position info from '//trim(infile)
   ! Read atomic positions
   do i = 1, natom
      read( ioin, '(A)', iostat=ioerr ) input
      call error( subname, info, ioerr )
      read(input, *, iostat=ioerr ) ETmp, atpos(:,i), atrel(:,i)
      if ( ioerr.ne.0 ) then
         read(input, *, iostat=ioerr ) ETmp, atpos(:,i)
         call error( subname, info, ioerr )
         atrel(:,i) = 1
      endif

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
   !
   select case (ptyp)
   case ( 1 )
      do i = 1, natom
      do j = 1, 3
         atpos(j, i) = atpos(j, i) * alat
      enddo
      enddo
      cartesian = .true.
   case ( 2 )
      do i = 1, natom
      do j = 1, 3
         atpos(j, i) = atpos(j, i) * bohr2ang
      enddo
      enddo
      cartesian = .true.
   case ( 3 )
      cartesian = .false.
   case ( 4 )
      cartesian = .true.
   end select
   !
   call axis2abc()
   !
contains
   !
   logical function matches(str1, str2)
   implicit none
      !------------------------------------------------------
      character (len=*), intent(in) :: str1, str2
      integer :: len1, len2, l
      !------------------------------------------------------
      len1 = len_trim(str1)
      len2 = len_trim(str2)

      do l = 1, len2-len1+1
         if (str1(1:len1) == str2(l:l+len1-1) ) then
            matches = .true.
            return
         endif
      enddo
      matches = .false.

      return
   end function matches
   !---------------------------------------------------------
end subroutine
!
subroutine writepwscf
use prec
use cell
use iounits
implicit none
   !----------------------------------------------------------------------------
   integer :: i, j, ibrav = 14
   real(q), parameter :: bohr2ang = 0.5291772108D0
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
   write( ioout, 400 ) "celldm(1) = ", a / bohr2ang
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
