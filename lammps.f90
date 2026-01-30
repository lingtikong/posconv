! To read LAMMPS full style data file, and write LAMMPS full and/or atomic style datafile
!
module lmp_full
use prec
implicit none
   !----------------------------------------------------------------------------
   integer :: atoms = 0, bonds=0, angles=0, dihedrals=0
   integer :: atomtypes=0, bondtypes=0, angletypes=0, dihedraltypes=0
   real(q) :: xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz
   !
   real(q), allocatable :: Masses(:), Velocities(:,:)
   integer, allocatable :: MolID(:), images(:,:)
   integer, allocatable :: Bondlist(:,:), Anglelist(:,:), Dihedrallist(:,:)
   character(len=256)   :: lmptitle
   character(len=512),allocatable :: PairCoeffs(:), BondCoeffs(:)
   character(len=512),allocatable :: AngleCoeffs(:), DihedralCoeffs(:)
   !
   logical :: lmp_full_read = .false.
   !----------------------------------------------------------------------------
end module

! To read LAMMPS full style data file
subroutine readlmpfull()
use cell
use iounits
use lmp_full
implicit none
   !----------------------------------------------------------------------------
   real(q)             :: dqdum(3), chgdum
   integer             :: ioerr, rderr, i, aid, mid, tid, idumarray(5)
   character (len=512) :: oneline
   character (len=50 ) :: strtmp, keyword
   !-----------------------------------------------------------------
   subname   = 'readlmpfull'
   atoms     = 0
   bonds     = 0
   angles    = 0
   dihedrals = 0
   atomtypes = 0
   bondtypes = 0
   angletypes = 0
   dihedraltypes = 0
   xlo = 0.D0
   xhi = 0.D0
   ylo = 0.D0
   yhi = 0.D0
   zlo = 0.D0
   zhi = 0.D0
   xy  = 0.D0
   xz  = 0.D0
   yz  = 0.D0
   !
   read(ioin,'(A)', iostat=ioerr) lmptitle
   read(ioin,'(A)', iostat=ioerr) oneline
   do while (ioerr.eq.0)
      i = index(oneline, '#')
      if (i.gt.0) then
         keyword = oneline(1:i-1)
      else
         keyword = trim(oneline)
      endif

      select case ( trim(keyword) )
      case ( 'Masses' )
         if (atomtypes.lt.1) then
            info = 'Masses defined before number of atom types!'
            call error( subname, info, 1)
         endif
         info = 'Error while reading mass info!'
         read(ioin, '(A)', iostat=ioerr) oneline
         if (ioerr.ne.0) call error( subname, info, 1)
         if (allocated(Masses)) deallocate(Masses)
         allocate( Masses(atomtypes) )
         do i = 1, atomtypes
            read(ioin, *, iostat=rderr) idum, dqdum(1)
            if (rderr.ne.0) call error( subname, info, 1)
            if (idum.gt.atomtypes) then
               info = 'Type id greater than atom_types!'
               call error( subname, info, idum-atomtypes)
            endif
            Masses(idum) = dqdum(1)
         enddo

      case ( 'Pair Coeffs' )
         if (atomtypes.lt.1) then
            info = 'Pair Coeffs defined before number of atom types!'
            call error( subname, info, 1)
         endif
         info = 'Error while reading Pair Coeffs info!'
         read(ioin, '(A)', iostat=ioerr) oneline
         if (ioerr.ne.0) call error( subname, info, 1)
         if (allocated(PairCoeffs)) deallocate(PairCoeffs)
         allocate( PairCoeffs(atomtypes) )
         do i = 1, atomtypes
            read(ioin, '(A)', iostat=rderr) PairCoeffs(i)
            if (rderr.ne.0) call error( subname, info, 1)
         enddo

      case ( 'Bond Coeffs' )
         if (bondtypes.lt.1) then
            info = 'Bond Coeffs defined before number of bond types!'
            call error( subname, info, 1)
         endif
         info = 'Error while reading Bond Coeffs info!'
         read(ioin, '(A)', iostat=ioerr) oneline
         if (ioerr.ne.0) call error( subname, info, 1)
         if (allocated(BondCoeffs)) deallocate(BondCoeffs)
         allocate( BondCoeffs(bondtypes) )
         do i = 1, bondtypes
            read(ioin, '(A)', iostat=rderr) BondCoeffs(i)
            if (rderr.ne.0) call error( subname, info, 1)
         enddo

      case ( 'Angle Coeffs' )
         if (angletypes.lt.1) then
            info = 'Angle Coeffs defined before number of angle types!'
            call error( subname, info, 1)
         endif
         info = 'Error while reading Angle Coeffs info!'
         read(ioin, '(A)', iostat=ioerr) oneline
         if (ioerr.ne.0) call error( subname, info, 1)
         if (allocated(AngleCoeffs)) deallocate(AngleCoeffs)
         allocate( AngleCoeffs(angletypes) )
         do i = 1, angletypes
            read(ioin, '(A)', iostat=rderr) AngleCoeffs(i)
            if (rderr.ne.0) call error( subname, info, 1)
         enddo

      case ( 'Dihedral Coeffs' )
        if (dihedraltypes.lt.1) then
            info = 'Dihedral Coeffs defined before number of dihedral types!'
            call error( subname, info, 1)
         endif
         info = 'Error while reading Dihedral Coeffs info!'
         read(ioin, '(A)', iostat=ioerr) oneline
         if (ioerr.ne.0) call error( subname, info, 1)
         if (allocated(DihedralCoeffs)) deallocate(DihedralCoeffs)
         allocate( DihedralCoeffs(dihedraltypes) )
         do i = 1, dihedraltypes
            read(ioin, '(A)', iostat=rderr) DihedralCoeffs(i)
            if (rderr.ne.0) call error( subname, info, 1)
         enddo

      case ( 'Atoms' )
         if (atoms.lt.1) then
            info = 'Atomic position defined before number of atoms!'
            call error( subname, info, 1)
         endif
         if ( allocated( atpos ) ) deallocate(atpos)
         if ( allocated( atrel ) ) deallocate(atrel)
         if ( allocated( attyp ) ) deallocate(attyp)
         if ( allocated( atchg ) ) deallocate(atchg)
         if ( allocated( MolID ) ) deallocate(MolID)
         if ( allocated(images ) ) deallocate(images)
         allocate( atpos(3,atoms), attyp(atoms), atrel(3,atoms), atchg(atoms), MolID(atoms), images(3,atoms) ) 
         atrel = 1
         info = 'Error while reading atomic positions!'
         read(ioin, '(A)', iostat=ioerr) oneline
         do i = 1, atoms
            read(ioin, '(A)', iostat=rderr) oneline
            if (rderr.ne.0) call error( subname, info, 1)
            oneline = trim(oneline)//" 0 0 0"
            read(oneline,*,iostat=rderr) aid, mid, tid, chgdum, dqdum, idumarray(1:3)
            if (rderr.eq.0) then
               atpos(:,aid) = dqdum
               attyp(aid)   = tid
               atchg(aid)   = chgdum
               MolID(aid)   = mid
               images(:,aid)= idumarray(1:3)
            else
               call error( subname, info, 1)
            endif
         enddo

      case ( 'Velocities' )
         if (atoms.lt.1) then
            info = 'Velocities defined before number of atoms!'
            call error( subname, info, 1)
         endif
         if ( allocated( Velocities ) ) deallocate(Velocities)
         allocate( Velocities(3,atoms) )
         info = 'Error while reading velocities!'
         read(ioin, '(A)', iostat=ioerr) oneline
         do i = 1, atoms
            read(ioin, *, iostat=rderr) aid, dqdum
            if (rderr.eq.0) then
               Velocities(:,aid) = dqdum
            else
               call error( subname, info, 1)
            endif
         enddo

      case ( 'Bonds' )
         if (bonds.lt.1) then
            info = 'Bond list defined before number of bonds!'
            call error( subname, info, 1)
         endif
         if ( allocated( Bondlist ) ) deallocate(Bondlist)
         allocate( Bondlist(3,bonds) )
         info = 'Error while reading bond list!'
         read(ioin, '(A)', iostat=ioerr) oneline
         do i = 1, bonds
            read(ioin, *, iostat=rderr) idum, idumarray(1:3)
            if (rderr.eq.0) then
               Bondlist(:,idum) = idumarray(1:3)
            else
               call error( subname, info, 1)
            endif
         enddo

      case ( 'Angles' )
         if (angles.lt.1) then
            info = 'Angle list defined before number of angles!'
            call error( subname, info, 1)
         endif
         if ( allocated( Anglelist ) ) deallocate(Anglelist)
         allocate( Anglelist(4,angles) )
         info = 'Error while reading bond list!'
         read(ioin, '(A)', iostat=ioerr) oneline
         do i = 1, angles
            read(ioin, *, iostat=rderr) idum, idumarray(1:4)
            if (rderr.eq.0) then
               Anglelist(:,idum)   = idumarray(1:4)
            else
               call error( subname, info, 1)
            endif
         enddo

      case ( 'Dihedrals' )
         if (dihedrals.lt.1) then
            info = 'Dihedral list defined before number of dihedrals!'
            call error( subname, info, 1)
         endif
         if ( allocated( Dihedrallist ) ) deallocate(Dihedrallist)
         allocate( Dihedrallist(5,dihedrals) )
         info = 'Error while reading bond list!'
         read(ioin, '(A)', iostat=ioerr) oneline
         do i = 1, dihedrals
            read(ioin, *, iostat=rderr) idum, idumarray(1:5)
            if (rderr.eq.0) then
               Dihedrallist(:,idum) = idumarray(1:5)
            else
               call error( subname, info, 1)
            endif
         enddo

      case default
         read(oneline, *, iostat=rderr) dqdum, strtmp
         if (rderr.eq.0.and.trim(strtmp).eq.'xy') then
            xy = dqdum(1)
            xz = dqdum(2)
            yz = dqdum(3)
         else
            read(oneline,*,iostat=rderr) dqdum(1:2), strtmp
            if (rderr.eq.0) then
               select case ( trim(strtmp) )
               case ( 'xlo' )
                  xlo = dqdum(1)
                  xhi = dqdum(2)
               case ( 'ylo' )
                  ylo = dqdum(1)
                  yhi = dqdum(2)
               case ( 'zlo' )
                  zlo = dqdum(1)
                  zhi = dqdum(2)
               end select
            else
               read(oneline,*,iostat=rderr) idum, strtmp
               if (rderr.eq.0) then
                  select case ( trim(strtmp) )
                  case ( 'atoms' )
                     atoms = idum
                  case ( 'bonds' )
                     bonds = idum
                  case ( 'angles' )
                     angles = idum
                  case ( 'dihedrals' )
                     dihedrals = idum
                  case ( 'atom' )
                     atomtypes = idum
                  case ( 'bond' )
                     bondtypes = idum
                  case ( 'angle' )
                     angletypes = idum
                  case ( 'dihedral' )
                     dihedraltypes = idum
                  end select
               endif
            endif
         endif
      end select
      read(ioin,'(A)',iostat=ioerr) oneline
   enddo
   !
   if ( (xhi.le.xlo).or.(yhi.le.ylo).or.(zhi.le.zlo) ) then
      info = 'Box size info is incorrect!'
      call error( subname, info, 2)
   endif
   if ( atoms.lt.1.or.atomtypes.lt.1 ) then
      info = 'No atom identified in file!'
      call error( subname, info, 3)
   endif
   !
   title = trim(lmptitle)
   ntype = atomtypes
   natom = atoms
   !
   axis(1,:) = (/ xhi-xlo, 0.D0, 0.D0 /)
   axis(2,:) = (/ xy, yhi-ylo, 0.D0 /)
   axis(3,:) = (/ xz, yz, zhi-zlo /)
   alat      = 1.D0
   !
   ! To assign all other cell related variables
   allocate( ntm(ntype), EName(ntype) )
   forall( i=1:ntype ) ntm(i) = count(attyp.eq.i)
   do i = 1, ntype
      write(EName(i),'(I2)') i
   enddo
   Eread(1:ntype) = EName(1:ntype)
   !
   cartesian = .true.
   call axis2abc()
   !
   lmp_full_read = .true.
   !
return
end subroutine

!*******************************************************************************
!*** Subroutines to write LAMMPS atomic format position file                 ***
!*******************************************************************************
subroutine writelmp_atomic
use prec
use cell
use iounits
use lmp_full
implicit none
   !----------------------------------------------------------------------------
   integer             :: i, istr, iend, ichg, ioerr
   character (len=100) :: fmtstr
   !----------------------------------------------------------------------------
   call abc2axis()
   if ( .not.cartesian ) call dir2car
   !
   ichg = 0
   if ( allocated(atchg) ) then
      write(*, '(/,10x,"Charge info exists, do you want to output as charge style? (y/n)[n]: ", $)')
      read(*, '(A)', iostat = ioerr) fmtstr
      if ( ioerr.eq.0.and.(fmtstr.eq.'y'.or.fmtstr.eq.'Y') ) ichg = 1
   endif

   write( ioout, 500 ) trim( title )
   fmtstr = '(I??,2x,A)'
   write(fmtstr(3:4),'(I2)') int(log10(dble(natom)))+1
   write( ioout, fmtstr ) natom, "atoms"
   write( ioout, fmtstr ) ntype, "atom types"
   write( ioout, * )
   !
   if ( lmp_full_read ) then
      write( ioout, 520 ) xlo, xhi, "xlo xhi"
      write( ioout, 520 ) ylo, yhi, "ylo yhi"
      write( ioout, 520 ) zlo, zhi, "zlo zhi"
      if (abs(xy)+abs(xz)+abs(yz).gt.1.D-8) then
         write( ioout, 525 ) xy, xz, yz
      endif
   else
      write( ioout, 520 ) 0.D0, axis(1,1)*alat, "xlo xhi"
      write( ioout, 520 ) 0.D0, axis(2,2)*alat, "ylo yhi"
      write( ioout, 520 ) 0.D0, axis(3,3)*alat, "zlo zhi"
      if (abs(cosac)+abs(cosbc)+abs(cosab).gt.1.D-8) then
         write( ioout, 525 ) axis(2,1)*alat, axis(3,1:2)*alat
      endif
   endif
   write( ioout, 530 )
   !
   if (ichg.eq.0) then
      write( fmtstr,'("(I",I2.2,",1x,I2,3(1x,F15.8))")') int(log10(dble(natom)+0.1)+1)
      do i = 1, natom
         write( ioout, fmtstr ) i, attyp(i), atpos(:,i)
      enddo
   else
      write( fmtstr,'("(I",I2.2,",1x,I2,1x,F8.4,3(1x,F15.8))")') int(log10(dble(natom)+0.1)+1)
      do i = 1, natom
         write( ioout, fmtstr ) i, attyp(i), atchg(i), atpos(:,i)
      enddo
   endif
   !
500 format(A,/)
510 format(I10,2x,A,/)
520 format(2(F15.8,1x),A)
525 format(3(F15.8,1x),"xy xz yz")
530 format(/,"Atoms",/)
   !
return
end subroutine

! To write LAMMPS full style data file; it requires to read a full style datafile
! as input. Thus enable one to use modify some properties of the atomic positions.
subroutine writelmp_full()
use cell
use iounits
use lmp_full
implicit none
   !----------------------------------------------------------------------------
   integer             :: i, istr, iend
   character (len=100) :: fmtstr
   !----------------------------------------------------------------------------
   subname = 'WriteLMPFull'
   if ( .not.lmp_full_read ) then
      info = 'To write LAMMPS full format, you need to read in LAMMPS full format first!'
      call error( subname, info, 1)
   endif
   if ( .not.cartesian ) call dir2car
   !
   fmtstr = '(I??,2x,A)'
   write(fmtstr(3:4),'(I2)') int(log10(dble(max(atoms,bonds,angles,dihedrals))))+1
   write( ioout, 500 ) trim( lmptitle )
   write( ioout, fmtstr ) atoms, "atoms"
   write( ioout, fmtstr ) bonds, "bonds"
   write( ioout, fmtstr ) angles, "angles"
   write( ioout, fmtstr ) dihedrals, "dihedrals"
   write( ioout, * )
   write(fmtstr(3:4),'(I2)') int(log10(dble(max(atomtypes,bondtypes,angletypes,dihedraltypes))))+1
   write( ioout, fmtstr ) atomtypes, "atom types"
   write( ioout, fmtstr ) bondtypes, "bond types"
   write( ioout, fmtstr ) angletypes, "angle types"
   write( ioout, fmtstr ) dihedraltypes, "dihedral types"
   write( ioout, * )
   write( ioout, 520 ) xlo, xhi, "xlo xhi"
   write( ioout, 520 ) ylo, yhi, "ylo yhi"
   write( ioout, 520 ) zlo, zhi, "zlo zhi"
   if (abs(xy)+abs(xz)+abs(yz).gt.0.D0) write( ioout, 525 ) xy, xz, yz
   !
   if ( allocated(Masses) ) then
      write( ioout, 530 ) "Masses"
      do i = 1, atomtypes
         write( ioout, 540 ) i, Masses(i)
      enddo
   endif
   if ( allocated(PairCoeffs) ) then
      write( ioout, 530 ) "Pair Coeffs"
      do i = 1, atomtypes
         write( ioout, 550 ) trim(PairCoeffs(i))
      enddo
   endif
   if ( allocated(BondCoeffs) ) then
      write( ioout, 530 ) "Bond Coeffs"
      do i = 1, bondtypes
         write( ioout, 550 ) trim(BondCoeffs(i))
      enddo
   endif
   if ( allocated(AngleCoeffs) ) then
      write( ioout, 530 ) "Angle Coeffs"
      do i = 1, angletypes
         write( ioout, 550 ) trim(AngleCoeffs(i))
      enddo
   endif
   if ( allocated(DihedralCoeffs) ) then
      write( ioout, 530 ) "Dihedral Coeffs"
      do i = 1, dihedraltypes
         write( ioout, 550 ) trim(DihedralCoeffs(i))
      enddo
   endif
   fmtstr = '(I??,1x,I??,1x,I?,1x,F10.5,3(1x,F19.12),3(1x,I??))'
   write(fmtstr(3:4),'(I2)') int(log10(dble(atoms)))+1
   write(fmtstr(10:11),'(I2)') int(log10(dble(maxval(MolID))))+1
   write(fmtstr(17:17),'(I1)') int(log10(dble(atomtypes))) + 1
   write(fmtstr(47:48),'(I2)') int(log10(dble(maxval(abs(images))))) + 2 ! Could be negative
   !
   write( ioout, 530 ) "Atoms"
   do i = 1, atoms
      write( ioout, fmtstr ) i, MolID(i), attyp(i), atchg(i), atpos(:,i), images(:,i)
   enddo
   !
   if ( allocated(Velocities) ) then
      write( ioout, 530 ) "Velocities"
      fmtstr = '(I??,3(1x,F20.12))'
      write(fmtstr(3:4),'(I2)') int(log10(dble(atoms)))+1
      do i = 1, atoms
         write( ioout, fmtstr ) i, Velocities(:,i)
      enddo
   endif
   if ( allocated(Bondlist) ) then
      write( ioout, 530 ) "Bonds"
      fmtstr = '(I??,1x,I?,2(1x,I??))'
      write(fmtstr(3:4),'(I2)') int(log10(dble(bonds)))+1
      write(fmtstr(10:10),'(I1)') int(log10(dble(bondtypes)))+1
      write(fmtstr(18:19),'(I2)') int(log10(dble(atoms)))+1
      do i = 1, bonds
         write( ioout, fmtstr ) i, Bondlist(:,i)
      enddo
   endif
   if ( allocated(Anglelist) ) then
      write( ioout, 530 ) "Angles"
      fmtstr = '(I??,1x,I?,3(1x,I??))'
      write(fmtstr(3:4),'(I2)') int(log10(dble(angles)))+1
      write(fmtstr(10:10),'(I1)') int(log10(dble(angletypes)))+1
      write(fmtstr(18:19),'(I2)') int(log10(dble(atoms)))+1
      do i = 1, angles
         write( ioout, fmtstr ) i, Anglelist(:,i)
      enddo
   endif
   if ( allocated(Dihedrallist) ) then
      write( ioout, 530 ) "Dihedrals"
      fmtstr = '(I??,1x,I?,4(1x,I??))'
      write(fmtstr(3:4),'(I2)')   int(log10(dble(dihedrals)))+1
      write(fmtstr(10:10),'(I1)') int(log10(dble(dihedraltypes)))+1
      write(fmtstr(18:19),'(I2)') int(log10(dble(atoms)))+1
      do i = 1, dihedrals
         write( ioout, fmtstr ) i, Dihedrallist(:,i)
      enddo
   endif

500 format(A,/)
510 format(I10,2x,A)
520 format(2(F15.8,1x),A)
525 format(3(F15.8,1x),"xy xz yz")
530 format(/,A,/)
540 format(I3,1x,F15.6)
550 format(A)
   !
return
end subroutine
