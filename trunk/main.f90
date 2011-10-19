!*******************************************************************************
!***                          POSCONVERT                                     ***
!*** Program to convert between different configuraion file formats.         ***
!***                                                                         ***
!*** L.T. Kong,   Sep. 2006                                                  ***
!*** konglt@gmail.com                                                        ***
!*******************************************************************************
program main
use prec
use cell
use datetime
implicit none
   !----------------------------------------------------------------------------
   character (len=40) :: input
   character (len=1 ) :: direction(3) = (/ "X", "Y", "Z" /)
   integer            :: ioerr, iref, nchar, idir, istr, iend, i, j
   real(q)            :: reflect, factor, dpdum, postmp(3)
   logical            :: fexst, first = .true.
   !----------------------------------------------------------------------
   subname = 'main'
   info    = 'Wrong input!'
   call date_time
   !
   write(*, '(//, 40x, "POSCONVERT" )' )
   write(*, '( 35x, A4,"-",A2,"-",A2,2X,A2,":",A2,":",A2 )' ) datestr(1:4), &
   &    datestr(5:6), datestr(7:8), timestr(1:2), timestr(3:4), timestr(5:6)
   write(*, '(/, 10x, "Program to convert atomic configuration files"  )' )
   write(*, '(//, 10x, 23("="), 2x, "Input Configuration", 3x, 23("=") )' )
   write(*, '( 10x, "Please select input configuration format:"        )' )
   write(*, '( 15x, "1. VASP POSCAR;"                                  )' )
   write(*, '( 15x, "2. PWSCF POSITION CARD;"                          )' )
   write(*, '( 15x, "3. xyz file;"                                     )' )
   write(*, '( 15x, "4. groF and/or MSS format;"                       )' )
   write(*, '( 15x, "5. BGF (ReaxFF);"                                 )' )
   write(*, '( 15x, "6. ReaxFF save file;"                             )' )
   write(*, '( 15x, "7. MS Car file;"                                  )' )
   write(*, '( 15x, "8. LAMMPS full;"                                  )' )
   write(*, '( 15x, "9. LAMMPS dump atom;"                             )' )
   write(*, '( 15x, "0. Exit."                                         )' )
   write(*, '( 10x, "Your choice [3]:", $                              )' )
   read (*, '(A)', iostat=ioerr ) input
   Ninpostypes = 9  ! change its value if new type to read is enabled.
   call error( subname, info, ioerr )
   !
   if ( input.eq.'' ) then
      postyp = 3
   else
      read(input, *, iostat=ioerr ) postyp
      call error( subname, info, ioerr )
      if ( postyp.lt.1.or.postyp.gt.Ninpostypes ) stop
   endif
   !
   infile = defaultIns( postyp )
   !
   if ( iargc().ge.1 ) then
      call getarg(1, infile)
      write(*,'(/,10x,"Atomic configuration will be read from: ",A,/,10x,70("="))') trim(infile)
   else
      write(*, '(/,10x, "Please input the configuration filename [", A, "]: ", $    )' ) trim(infile)
      read (*, '(A)', iostat=ioerr ) input
      call error( subname, info, ioerr )
      if ( input.ne.'' ) then
         read( input, *, iostat=ioerr ) infile
         call error( subname, info, ioerr )
      endif
      write(*, '( 10x, 30("="), 2x, "Thanks", 2x, 30("=")                 )' )
   endif
   !
   inquire( file=infile, exist =fexst )
   if ( .not.fexst ) then
      info = "File:"//trim(infile)//" not found!"
      call error(subname, info, 1)
   endif
   ang2rad = atan(1.D0) / 45.D0
   rad2ang = 1.D0/ang2rad
   !
   call readpos
   !
   call DisplayTypeRead
   !
   call OperatePos
   !
   write(*, '(//, 10x, 23("="), 2x, "Output Configuration", 2x, 23("=") )' )
   write(*, '( 10x, "Please select output configuration format:"        )' )
   write(*, '( 15x, "1. VASP POSCAR;"                                  )' )
   write(*, '( 15x, "2. PWSCF POSITION CARD;"                          )' )
   write(*, '( 15x, "3. xyz file;"                                     )' )
   write(*, '( 15x, "4. groF and/or MSS format;"                       )' )
   write(*, '( 15x, "5. BGF (ReaxFF);"                                 )' )
   write(*, '( 15X, "6. ReaxFF save file;"                             )' )
   write(*, '( 15X, "7. Car (MS);"                                     )' )
   write(*, '( 15X, "8. Res (MS);"                                     )' )
   write(*, '( 15X, "9. LAMMPS atomic;"                                )' )
   write(*, '( 14X, "10. LAMMPS full;"                                 )' )
   write(*, '( 14X, "11. LAMMPS dump atom;"                            )' )
   write(*, '( 15x, "0. Exit."                                         )' )
   write(*, '( 10x, "Your choice [3]:", $                              )' )
   read (*, '(A)', iostat=ioerr ) input
   Noutpostypes = 11 ! change if new format is enabled
   call error( subname, info, ioerr )
   !
   if ( input.eq.'' ) then
      outtyp = 3
   else
      read(input, *, iostat=ioerr ) outtyp
      call error( subname, info, ioerr )
      if ( outtyp.lt.1.or.outtyp.gt.Noutpostypes ) then
         write( *, '( 10x,  70("=") )' )
         stop
      endif
   endif
   !
   outfile = defaultOut( outtyp )
   !
   write(*, '(/,10x, "Please input the output configuration filename[", A, "]: ", $    )' ) trim(outfile)
   read (*, '(A)', iostat=ioerr ) input
   call error( subname, info, ioerr )
   if ( input.ne.'' ) then
      read( input, *, iostat=ioerr ) outfile
      call error( subname, info, ioerr )
   endif
   !
   nchar = len_trim(outfile)
   !
   select case ( outtyp )
   case ( 3 )
      if ( nchar.lt.4 ) then
         outfile = trim(outfile)//".xyz"
      else
         input = outfile((nchar-3):nchar)
         do i = 2, 4
            if (input(i:i).ge.'a'.and.input(i:i).le.'z') input(i:i) = char(ichar(input(i:i)) -32)
         enddo
         if ( input.ne.'.XYZ' ) outfile = trim(outfile)//".xyz"
      endif
   case ( 7 )
      if ( nchar.lt.4 ) then
         outfile = trim(outfile)//".car"
      else
         input = outfile((nchar-3):nchar)
         do i = 2, 4
            if (input(i:i).ge.'a'.and.input(i:i).le.'z') input(i:i) = char(ichar(input(i:i)) -32)
         enddo
         if ( input.ne.'.CAR' ) outfile = trim(outfile)//".car"
      endif
   case ( 8 )
      if ( nchar.lt.4 ) then
         outfile = trim(outfile)//".res"
      else
         input = outfile((nchar-3):nchar)
         do i = 2, 4
            if (input(i:i).ge.'a'.and.input(i:i).le.'z') input(i:i) = char(ichar(input(i:i)) -32)
         enddo
         if ( input.ne.'.RES' ) outfile = trim(outfile)//".res"
      endif
   end select
   write(*, '( 10x, 30("="), 2x, "Thanks", 2x, 30("=")                 )' )
   !
   call writepos
   !
   write(*, '(/, 10x, "Mission completed!" )' )
   !
stop
end program
