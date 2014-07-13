!*******************************************************************************
!*** Subroutines to read/write atom format dump file of LAMMPS               ***
!*******************************************************************************
subroutine read_lmp_atom
use prec
use cell
use iounits
implicit none
   !-----------------------------------------------------------------
   integer            :: ioerr, i, j, jp, jd
   real(q)            :: radum(3), vol, xlo, xhi, ylo, yhi, zlo, zhi
   real(q)            :: xy, xz, yz
   character (len=256):: input
   character (len=20 ):: strtmp
   integer, external  :: typescreen
   !-----------------------------------------------------------------
   subname = 'read_lmp_atom'
   !
   ! Read total number of atoms
   read( ioin, '(A)', iostat=ioerr ) input
   read( ioin, *, iostat=ioerr ) idum
   call error( subname, info, ioerr )
   write(title,'("LAMMPS dump atom at timestep: ", I12)') idum
   read( ioin, '(A)', iostat=ioerr ) input
   read( ioin, *, iostat=ioerr ) natom
   call error( subname, info, ioerr )
   !
   if ( natom.lt.1 ) then
      info = 'Too few atoms read from: '//trim(infile)
      call error( subname, info, abs(natom-1) )
   endif
   allocate( atpos(3, natom), atrel(3,natom), attyp(natom) )
   atrel = 1
   xy = 0.D0
   xz = 0.D0
   yz = 0.D0
   ! Read box dimensions
   read( ioin, '(A)', iostat=ioerr ) input
   read( ioin, '(A)', iostat=ioerr ) input
   read( input, *, iostat=ioerr ) xlo, xhi, xy
   if (ioerr.ne.0) read( input, *, iostat=ioerr ) xlo, xhi

   read( ioin, '(A)', iostat=ioerr ) input
   read( input, *, iostat=ioerr ) ylo, yhi, xz
   if (ioerr.ne.0) read( input, *, iostat=ioerr ) ylo, yhi

   read( ioin, '(A)', iostat=ioerr ) input
   read( input, *, iostat=ioerr ) zlo, zhi, yz
   if (ioerr.ne.0) read( input, *, iostat=ioerr ) zlo, zhi
   call error( subname, info, ioerr )
   !
   xlo = xlo - MIN(0., MIN(xy, MIN(xz, xy+xz)))
   xhi = xhi - MAX(0., MAX(xy, MAX(xz, xy+xz)))

   ylo = ylo - MIN(0.,yz)
   yhi = yhi - MAX(0.,yz)
   !
   axis = 0.D0
   axis(1,1) = xhi - xlo
   axis(2,1) = xy
   axis(2,2) = yhi - ylo
   axis(3,1) = xz
   axis(3,2) = yz
   axis(3,3) = zhi - zlo
   !
   ! Read atomic positions
   read( ioin, '(A)', iostat=ioerr )
   do i = 1, natom
      read( ioin, '(A)', iostat=ioerr ) input
      call error( subname, info, ioerr )
      read(input, *, iostat=ioerr ) j, Etmp, radum
      read(Etmp, *, iostat=ioerr) jp
      jd = typescreen(Etmp)
      typeID(jd) = jp
      attyp(j)   = jp
      atpos(:,j) = radum
   enddo
   !
   if ( ntype.lt.1 ) then
      info = 'Too few elements identified from:'//trim(infile)
      call error( subname, info, abs(ntype-1) )
   endif
   !
   ! To assign all other cell related variables
   allocate( ntm(ntype), EName(ntype) )
   do i = 1, ntype
      ntm(i) = count(attyp.eq.i)
   enddo
   Ename = Eread(1:ntype)
   alat  = 1.D0
   !
   cartesian = .false.
   call axis2abc()
   !
return
end subroutine
!
subroutine write_lmp_atom
use prec
use cell
use iounits
implicit none
   !----------------------------------------------------------------------------
   integer :: i, j
   real(q)            :: xlo, xhi, ylo, yhi, zlo, zhi
   real(q)            :: xy, xz, yz
   !----------------------------------------------------------------------------
   !
   if ( cartesian ) call car2dir
!   if ( abs(alpha-90.D0).gt.1.D0.or.abs(beta-90.D0).gt.1.D0.or.abs(gamma-90.D0).gt.1.D0 ) then
!      subname = 'write_lmp_atom'
!      info = 'Box is not orthogonal for LAMMPS atom style dump!'
!      call warn(subname, info, 0)
!   endif
   !
   write( ioout, 100 )
   write( ioout, 110 ) 0
   write( ioout, 200 ) 
   write( ioout, 210 ) natom
   xlo = 0.
   ylo = 0.
   zlo = 0.
   xhi = axis(1,1) * alat
   yhi = axis(2,2) * alat
   zhi = axis(3,3) * alat
   xy  = axis(2,1) * alat
   xz  = axis(3,1) * alat
   yz  = axis(3,2) * alat
   if ( (xy*xy + xz*xz + yz*yz) > 1.e-10 ) then
      xlo = xlo + min(0., xy, xz, xy+xz)
      xhi = xhi + max(0., xy, xz, xy+xz)
      ylo = ylo + min(0., yz)
      yhi = yhi + max(0., yz)
      write( ioout, 301 )
      write( ioout, 320 ) xlo, xhi, xy
      write( ioout, 320 ) ylo, yhi, xz
      write( ioout, 320 ) zlo, zhi, yz
   else 
      write( ioout, 300 )
      write( ioout, 310 ) 0.D0, a
      write( ioout, 310 ) 0.D0, b
      write( ioout, 310 ) 0.D0, c
   endif
   write( ioout, 400 )
   !
   do i = 1, natom
      write( ioout, 410 ) i, attyp(i), atpos(:,i)
   enddo
   !
100 format("ITEM: TIMESTEP" )
110 format(I15)
200 format("ITEM: NUMBER OF ATOMS" )
210 format(I10)
300 format("ITEM: BOX BOUNDS xx yy zz" )
301 format("ITEM: BOX BOUNDS xy xz yz xx yy zz" )
310 format(2(1x,F20.10))
320 format(3(1x,F20.10))
400 format("ITEM: ATOMS id type xs ys zs" )
410 format(I10,1x,I3,3(1x,F20.10))
   !
return
end subroutine
