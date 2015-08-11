subroutine readpos
use cell
use iounits
implicit none
   !----------------------------------------------------------------------------
   integer :: ioerr, i
   !----------------------------------------------------------------------------
   axis    = 0.D0
   Eread   = ''
   !
   do i = 1, NMax
      typeID(i) = i
   enddo
   !
   if ( allocated( atpos  ) ) deallocate( atpos   )
   if ( allocated( atrel  ) ) deallocate( atrel   )
   if ( allocated( attyp  ) ) deallocate( attyp   )
   if ( allocated( ntm    ) ) deallocate( ntm     )
   if ( allocated( EName  ) ) deallocate( EName   )
   !
   open( ioin, file = infile, status = 'old', action = 'read', iostat = ioerr )
   info = 'Error encountered while reading file:'//trim(infile)
   !
   select case ( postyp )
   case ( 1 ) ! VASP
      call readvasp
   case ( 2 ) ! PWSCF POSITION CARD
      call readpwscf
   case ( 3 ) ! xyz file
      call readxyz
   case ( 4 ) ! groF and/or MSS and/or PCEAM format
      call readgroF
   case ( 5 ) ! BGF file
      call readbgf
   case ( 6 ) ! ReaxFF save file, moldyn.vel
      call readvel
   case ( 7 ) ! MS Car file
      call readcar
   case ( 8 ) ! LAMMPS full
      call readlmpfull
   case ( 9 ) ! LAMMPS dump atom
      call read_lmp_atom
   case ( 10 ) ! SIESTA STRUCT_IN
      call readsiesta
   case ( 11 ) ! Abinit/BigDFT xyz
      call read_abinit_xyz
   case ( 12 ) ! MS RES
      call readres
   case ( 13 ) ! ARTn
      call readart
   case default
      info = 'Un-supported format!'
      call error( subname, info, 1 )
   end select
   close( ioin )
   !
   if ( cartesian ) call car2dir()
   !-----------------------------------------------------------------
return
end subroutine
