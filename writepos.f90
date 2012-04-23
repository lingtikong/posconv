subroutine writepos
use cell
use iounits
implicit none
   !----------------------------------------------------------------------------
   if ( natom.lt.1 ) return
   if ( len_trim(outfile).eq.0 ) outfile = 'atomcfg'
   !
   open( ioout, file = outfile, status ='unknown', action='write' )
   !
   select case( outtyp )
   case ( 1 ) ! VASP POSCAR
      call writevasp
   case ( 2 ) ! PWSCF ATOMIC_POSITIONS CARD
      call writepwscf
   case ( 3 ) ! xyz
      call writexyz
   case ( 4 ) ! groF
      call writegroF
   case ( 5 ) ! BGF
      call writebgf
   case ( 6 ) ! ReaxFF save file
      call writevel
   case ( 7 ) ! Insight II ( MS Modeling Car )
      call writecar
   case ( 8 ) ! MS Modeling Res
      call writeres
   case ( 9 ) ! LAMMPS atomic style
      call writelmp_atomic
   case ( 10 ) ! LAMMPS full style
      call writelmp_full
   case ( 11 )
      call write_lmp_atom
   case ( 12 )
      call writesiesta ! SIESTA STRUCT_IN
   case ( 13 )         ! Abinit/BigDFT xyz
      call write_abinit_xyz
   case ( 14 )
      call latgen
   case default
      info = "Un-support format!"
      call error( subname, info, 1 )
   end select
   !
   close( ioout  )
   write( *, 100 ) trim(outfile)
100 format( 10X, "Atomic configuration is written to: ", A )
   !
return
END subroutine
