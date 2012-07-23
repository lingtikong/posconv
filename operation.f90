!*******************************************************************************
!*** To make some special operation on the positions.                        ***
!*******************************************************************************
subroutine OperatePos
use prec
use cell
implicit none
   !----------------------------------------------------------------------------
   character (len=40) :: input
   character (len=1 ) :: direction(3) = (/ "X", "Y", "Z" /)
   integer            :: ioerr, iref, nchar, idir, istr, iend, i, j, extend(3), nExt
   integer            :: ii, jj, kk, ip, ntmp
   real(q)            :: reflect, factor, dpdum, postmp(3)
   real(q)            :: rotx(3,3), roty(3,3), rotz(3,3), ralph, rbeta, rgamm
   real(q)            :: rcosa, rsina, rcosb, rsinb, rcosg, rsing
   logical            :: first = .true.
   integer            :: newID(NMax)
   integer, allocatable :: oneDint(:), twoDint(:,:)
   real(q), allocatable :: oneDdbl(:), twoDdbl(:,:)
   !----------------------------------------------------------------------------
   ! Special operation on atomic positions
   do while ( .true. )
      write(*, '(//, 10x, 15("="), 2x, "Special operation on atomic positions", 2x, 14("=") )' )
      write(*, '( 15x, "1. Reflect around y-z plane; | 10. Shift along x direction; " )' )
      write(*, '( 15x, "2. Reflect around x-z plane; | 11. Shift along y direction; " )' )
      write(*, '( 15x, "3. Reflect around x-y plane; | 12. Shift along z direction; " )' )
      write(*, '( 15x, "4. Rescale lattice constant; | 13. Sort along x direction;  " )' )
      write(*, '( 15x, "5. Rescale x axis;           | 14. Sort along y direction;  " )' )
      write(*, '( 15x, "6. Rescale y axis;           | 15. Sort along z direction;  " )' )
      write(*, '( 15x, "7. Rescale z axis;           | 16. Realign the molecular;   " )' )
      write(*, '( 15x, "8. Calculate nominal radius; | 17. Rotate the lattice;   " )' )
      write(*, '( 15x, "9. Apply pbc condition;      | 18. Calculate neighbor list;")')
      write(*, '( 14x,"20. Reset atomic types;       | 19. Shift the origin of cell;")')
      write(*, '( 14x,"21. Extend original cell;     |  ",$)' )
      !
      if ( first ) then
         write(*, '("0. As is." )' )
         first = .false.
      else
          write(*, '("0. Done." )' )
      endif
      write(*, '( 10x, 70("-") )' )
      write(*, '( 10x, "Your choice [0]: ", $ )' )
      read (*, '(A)', iostat=ioerr ) input
      if ( ioerr.ne.0.or.input.eq.'' ) then
         iref = 0
      else
         read(input, *, iostat=ioerr ) iref
         if ( ioerr.ne.0 ) iref = 0
      endif
      !
      select case ( iref )
      case ( 1:3 ) !To reflect the coordination by taking the x, y, and/or z plane as a mirror
         !
         if ( cartesian ) call car2dir()
         write(*, '(10x, "Please input the position of the reflect mirror [0.0]: ", $ )' )
         read (*, '(A)', iostat=ioerr ) input
         call error( subname, info, ioerr )
         if ( input.eq.'' ) then
            reflect = 0.D0
         else
            read(input, *, iostat=ioerr ) reflect
            call error( subname, info, ioerr )
         endif
         reflect = 1.D0 + reflect
         !
         atpos(iref, :) = mod( atpos(iref,:), 1.D0 )
         where( atpos(iref,:).lt.0.D0 ) atpos(iref,:) = atpos(iref,:) + 1.D0
         atpos(iref, :) = reflect - atpos(iref,:)
         atpos(iref, :) = mod( atpos(iref,:), 1.D0 )
         where( atpos(iref,:).lt.0.D0 ) atpos(iref,:) = atpos(iref,:) + 1.D0
         !
      case ( 4 ) ! Rescale the lattice constanst; if positive, take the new value; if negative, multiply on the current value.
         !
         if ( cartesian ) call car2dir( )
         write(*, '(10x, "Current lattice constant of the cell is : ", F20.12 )' ) alat
         write(*, '(10x, "Please input the new lattice constant   : ", $ )' )
         read (*, '(A)', iostat=ioerr ) input
         call error( subname, info, ioerr )
         if ( input.ne.'' ) then
            read(input, *, iostat=ioerr ) factor
            call error( subname, info, ioerr )
            if ( factor.lt.0.D0 ) then
               factor = - factor
               alat   = alat * factor
            else
               factor = factor / alat
               alat   = factor * alat
            endif
            write( *, '(10x, "Do you want to keep some atoms at the origional coordinations?(y/n)[n]:", $ )' )
            read ( *, '(A)', iostat=ioerr ) input
            call error( subname, info, ioerr )
            if ( input.eq.'Y'.or.input.eq.'y' ) then
               do while ( .true. )
                  write( *, '( 10x, "Please input the index range to keep, 0 0  to exit: ", $ )' )
                  read ( *, *, iostat= ioerr ) istr, iend
                  call error( subname, info, ioerr )
                  if ( istr.lt.1.or.iend.lt.1.or.istr.gt.natom.or.iend.gt.natom.or.istr.gt.iend ) exit
                  atpos(:,istr:iend) = atpos(:,istr:iend) / factor
               enddo
            endif
            call axis2abc()
         endif
         !
      case ( 5:7 ) ! Rescale the x, y, or z axis. If positive, take the new value; if negative, multiply on the current value.
         !
         if ( cartesian ) call car2dir()
         idir = iref - 4
         write(*, '(10x, "The length of current ", A1, " axis is: ", F20.12 )' ) direction(idir), latt(idir)
         write(*, '(10x, "Please input its new length    : ", $ )' )
         read (*, '(A)', iostat=ioerr ) input
         call error( subname, info, ioerr )
         if ( input.ne.'' ) then
            read(input, *, iostat=ioerr ) factor
            call error( subname, info, ioerr )
            if ( factor.lt.0.D0 ) then
               factor     = - factor
               latt(idir) = latt(idir) * factor
            else
               factor     = factor / latt(idir)
               latt(idir) = factor * latt(idir)
            endif
            write( *, '(10x, "Do you want to keep some atoms at the origional coordinations?(y/n)[n]:", $ )' )
            read ( *, '(A)', iostat=ioerr ) input
            call error( subname, info, ioerr )
            if ( input.eq.'Y'.or.input.eq.'y' ) then
               do while ( .true. )
                  write( *, '( 10x, "Please input the index range to keep, 0 0  to exit: ", $ )' )
                  read ( *, *, iostat= ioerr ) istr, iend
                  call error( subname, info, ioerr )
                  if ( istr.lt.1.or.iend.lt.1.or.istr.gt.natom.or.iend.gt.natom.or.istr.gt.iend ) exit
                  atpos(idir,istr:iend) = atpos(idir,istr:iend) / factor
               enddo
            endif
            !
            call abc2axis
         endif
         !
      case ( 8 ) ! To calcualte the nominal radius of the current configuration.
         call calnominal
         write(*, '(10x, 70("=") )' )
         stop
      case ( 9 ) ! apply periodic boundary condition on the coordinates
         if ( cartesian ) call car2dir()
         atpos = mod(atpos, 1.D0)
         where ( atpos.lt.0.D0 ) atpos = atpos + 1.D0
      case ( 10:12 ) ! shift along certain direction
         idir = iref - 9
         write(*,'(/,10x,"Please indicate the way to shift :")')
         write(*,'(  15x,"1. shift in cartesian by certain amount;")')
         write(*,'(  15x,"2. shift in fractional by certain amount;")')
         write(*,'(  15x,"0. Exit.")')
         write(*,'(  10x,"Your choice[0]: ",$)')
         read (*,'(A)',iostat=ioerr) input
         if (ioerr.eq.0.and.input.ne.'') then
            read(input,*,iostat=ioerr) idum
            if (ioerr.ne.0) idum = 0
            if ( idum.ne.0 ) then
               write(*, '(10x, "Please input the value that will be shifted to zero ",$)')
               write(*, '(10x, "in ", A1, " direction: ", $ )' ) direction(idir)
               read (*, '(A)', iostat=ioerr ) input
               call error( subname, info, ioerr )
               dpdum = 0.D0
               if ( input.ne.'' ) read(input, *, iostat=ioerr ) dpdum

               if ( idum.eq.1 ) then
                  if (.not.cartesian) call dir2car()
                  atpos(idir,:) = atpos(idir,:) - dpdum
               elseif ( idum.eq.2 ) then
                  if (cartesian) call car2dir()
                  atpos(idir,:) = atpos(idir,:) - dpdum
               endif
            endif
         endif
      case ( 13:15 ) ! sort the coordination file along certain direction
         write(*, '(10x, "You want to sort by (1) ascending or (2) descending ? [1]:", $ )' )
         read (*, '(A)', iostat=ioerr) input
         idum = 1
         if ( ioerr.eq.0 ) then
            read( input, *, iostat=ioerr ) idum
            if ( ioerr.ne.0 ) idum = 1
         endif
         idum = min(2, max(1,idum))
         !
         idir = iref - 12
         do i = 1, natom-1
            do j = i+1, natom
               if ( atpos(idir,j).gt.atpos(idir,i) ) then
                  if ( idum.eq.2 ) then
                     postmp     = atpos(:,i)
                     atpos(:,i) = atpos(:,j)
                     atpos(:,j) = postmp
                     istr       = attyp(i)
                     attyp(i)   = attyp(j)
                     attyp(j)   = istr
                  endif
               elseif ( atpos(idir,j).lt.atpos(idir,i) ) then
                  if ( idum.eq.1 ) then
                     postmp     = atpos(:,i)
                     atpos(:,i) = atpos(:,j)
                     atpos(:,j) = postmp
                     istr       = attyp(i)
                     attyp(i)   = attyp(j)
                     attyp(j)   = istr
                  endif
               endif
            enddo
         enddo
      case ( 16 ) ! Rotate the molecular
         call realign_molecular
      case ( 17 ) ! rotate the lattice
         write(*,'(/10x,"Please input the degrees to rotate around x, y and z axis: ", $)')
         read (*, *, iostat=ioerr) ralph, rbeta, rgamm
         if (ioerr.ne.0) then
            ralph = 0.D0
            rbeta = 0.D0
            rgamm = 0.D0
         endif
         ralph = ralph * ang2rad
         rbeta = rbeta * ang2rad
         rgamm = rgamm * ang2rad
         rcosa = cos(ralph)
         rsina = sin(ralph)
         rcosb = cos(rbeta)
         rsinb = sin(rbeta)
         rcosg = cos(rgamm)
         rsing = sin(rgamm)
         rotx = reshape((/1.D0,0.D0,0.D0,0.D0,rcosa,rsina,0.D0,-rsina,rcosa/),(/3,3/))
         roty = reshape((/rcosb,0.D0,-rsinb,0.D0,1.D0,0.D0,rsinb,0.D0,rcosb/),(/3,3/))
         rotz = reshape((/rcosg,rsing,0.D0,-rsing,rcosg,0.D0,0.D0,0.D0,1.D0/),(/3,3/))
         !
         if (cartesian) call car2dir()
         axis = matmul(rotx,axis)
         axis = matmul(roty,axis)
         axis = matmul(rotz,axis)
         call axis2abc()
         call dir2car()
         call abc2axis()
         !
      case ( 18 ) 
         call neighcal()
      case ( 19 )
         if ( .not.cartesian ) call dir2car
         write(*, '(10x, "Please indicate the atom index or coordinate to shift to origin: ",$)')
         read (*, '(A)', iostat=ioerr ) input
         if ( ioerr.eq.0 ) then
            read(input, *,iostat=ioerr) postmp
            if (ioerr.ne.0) then
               read(input,*,iostat=ioerr) istr
               if (ioerr.ne.0.or.istr.le.0.or.istr.gt.natom) then
                  postmp = 0.D0
               else
                  postmp = atpos(:,istr)
               endif
            endif
            forall (i=1:natom) atpos(:,i) = atpos(:,i) - postmp
            write(*,'(/,10x,"Original position of the new origin: ",3F12.6)') postmp
         endif
      case ( 20 ) ! Reset atomic types
         write(*,'(/,10x,"*** Caution: this might not always work as expected ***")');
         write(*,'(10x,"Atomic types current assigned:",/,12x,"Assigned  type IDs: ", $)')
         do i = 1, ntype
            write(*, '(I4,$)') typeID(i)
         enddo
         write(*, '(/,12x,"Read type name/IDs: ", $)')
         do i = 1, ntype
            write(*, '(A4,$)') trim(Eread(i))
         enddo
         write(*,'(/,10x, "Your desired new IDs: ", $)');
         read(*, *) newID(1:ntype)
         do i = 1, natom
            do j = 1, ntype
               if ( attyp(i).eq.typeID(j) ) exit
            enddo
            attyp(i) = newID(j)
         enddo
         typeID(1:ntype) = newID(1:ntype)

         write(*,'(/,10x,"Atomic types are now assigned as:",/,12x,"Assigned  type IDs: ", $)')
         do i = 1, ntype
            write(*, '(I4,$)') typeID(i)
         enddo
         write(*, '(/,12x,"Read type name/IDs: ", $)')
         do i = 1, ntype
            write(*, '(A4,$)') trim(Eread(i))
         enddo
         write(*,*)
      case ( 21 ) ! Extend the original cell
         call DisplayCellInfo
         write(*,'(/,10x,"Please input the number of extensions in x, y, and z [1 1 1]:", $)')
         read (*, '(A)', iostat=ioerr ) input
         if (ioerr.eq.0) then
            read(input, *, iostat=ioerr) extend
            if (ioerr.ne.0) cycle
            do i = 1, 3
               extend(i) = MAX(1,extend(i))
            enddo
            nExt = extend(1)*extend(2)*extend(3)
            if (nExt.le.1) cycle

            if (cartesian) call car2dir()
            do i = 1, 3
            do j = 1, 3
               axis(i,j) = axis(i,j)*real(extend(i))
            enddo
            enddo
            ntm = ntm * nExt
            ntmp = natom * nExt

            allocate(twoDdbl(3,natom), twoDint(3,natom), oneDint(natom))
            twoDdbl = atpos
            twoDint = atrel
            oneDint = attyp
            deallocate(atpos, atrel, attyp)
            allocate(atpos(3,ntmp), atrel(3,ntmp), attyp(ntmp))
            if (allocated(atchg)) then
               allocate(oneDdbl(natom))
               oneDdbl = atchg
               deallocate(atchg)
               allocate(atchg(ntmp))
            endif
            ip = 0
            do ii = 1, extend(1)
            do jj = 1, extend(2)
            do kk = 1, extend(3)
               do i = 1, natom
                  ip = ip + 1
                  atpos(:,ip) = twoDdbl(:,i) + real((/ ii-1, jj-1, kk-1 /))
                  atrel(:,ip) = twoDint(:,i)
                  attyp(ip) = oneDint(i)
                  if (allocated(atchg)) atchg(ip) = oneDdbl(i)
               enddo
            enddo
            enddo
            enddo
            deallocate(oneDint, twoDint, twoDdbl)
            if (allocated(oneDdbl)) deallocate(oneDdbl)
            natom = ntmp
            call axis2abc()
            call DisplayCellInfo
         endif
      case default
         exit
      end select
      write(*, '(10x, 70("=") )' )
   end do
   write(*, '(10x, 70("=") )' )
   !
return
end subroutine
