!*******************************************************************************
!*** Subroutines to read/write BGF format position file                      ***
!*******************************************************************************
subroutine readbgf
use prec
use cell
use iounits
implicit none
   !----------------------------------------------------------------------------
   integer            :: ioerr, i, j, iatom
   real(q)            :: vol
   character (len=256):: input
   character (len=80 ):: strtmp, flag, strarray(6)
   logical            :: lpread = .false.
   integer, external  :: typescreen, typeindex
   !----------------------------------------------------------------------------
   subname = 'readbgf'
   !
   natom = 0
   ntype = 0
   latt  = 0.D0
   !
   do while ( .true. )
      read( ioin,  '(A)', iostat=ioerr ) input
      if ( ioerr.ne.0 ) exit
      read(input, *, iostat=ioerr ) flag
      if ( ioerr.ne.0 ) exit
      !
      select case ( trim(flag) )
      case ( "DESCRP" )
         read( input, '(7X, A)', iostat=ioerr ) title
         call error(subname, info, ioerr )
      case ( "FIXATOMS" )
         read( input, *, iostat=ioerr ) strtmp, i, j
         call error(subname, info, ioerr )
      case ( "CRYSTX" )
         read( input, *, iostat=ioerr ) strtmp, latt
         call error(subname, info, ioerr )
         lpread = .true.
      case ( "HETATM" )
         read( input, 50, iostat=ioerr ) strarray(1:2)
         call error(subname, info, ioerr )
         read(strarray(2),*,iostat=ioerr) ETmp
         call error(subname, info, ioerr )
         idum  = typescreen( ETmp )
         natom = natom + 1
      case ( "END" )
         exit
      end select
   enddo
   !
   if ( natom.lt.1 ) then
      info = 'Too few atoms read from: '//trim(infile)
      call error( subname, info, abs(natom-1) )
   endif
   !
   if ( ntype.lt.1 ) then
      info = 'Too few elements identified from:'//trim(infile)
      call error( subname, info, abs(ntype-1) )
   endif
   !
   rewind( ioin )
   !
   allocate( atpos(3, natom), atrel(3, natom), attyp(natom), ntm(ntype), EName(ntype), atchg(natom) )
   !
   atrel = 1
   EName = Eread(1:ntype)
   !
   iatom = 0
   do while ( .true. )
      read( ioin,  '(A)', iostat=ioerr ) input
      if ( ioerr.ne.0 ) exit
      read(input, *, iostat=ioerr ) flag
      if ( ioerr.ne.0 ) exit
      !
      select case ( trim(flag) )
      case ( "HETATM" )
         iatom = iatom + 1
         read( input, 100, iostat=ioerr) strarray(1), i, strarray(2), atpos(:,iatom), strarray(3:4), atchg(iatom)
         call error(subname, info, ioerr )
         read(strarray(3),*,iostat=ioerr) ETmp
         call error(subname, info, ioerr )
         attyp(iatom)  = typeindex( ETmp )
      case ( "FIXATOMS" )
         read( input, *, iostat=ioerr ) strtmp, i, j
         call error( subname, info, ioerr )
         atrel(:, i:j) = 0
      case ( "END" )
         exit
      end select
      !
 50 format( A61, A5 )
100 format( A7,I5,A18,3F10.5,1x,A5,A6,F8.5 )
   enddo
   !
   ! To assign all other cell related variables
   if ( title.eq.'' ) title = 'Position read from:'//trim(infile)
   forall ( i=1:ntype) ntm(i) = count( attyp.eq.i )
   !
   if ( .not.lpread ) then
      input = info
      info  = 'Error encoutered while reading standard input!'
      write(*, '(/, 10X, "Please input the lattice parameters," )' )
      write(*, '( 10x, "a, b, and c in A :", $ )' )
      read (*, *, iostat=ioerr ) latt(1:3)
      call error(subname, info, ioerr)
      write(*, '( 10x, "angles in degree :", $ )' )
      read (*, *, iostat=ioerr ) latt(4:6)
      call error(subname, info, ioerr)
      info = input
   endif
   !
   call abc2axis
   !
   call volume( axis, vol)
   if ( vol.lt.1.D-6 ) then
      info = 'The read or input lattice parameters might be wrong!'
      call error(subname, info, 1)
   endif
   cartesian = .true.
   !
return
end subroutine
!
subroutine writebgf
use prec
use cell
use iounits
implicit none
   !----------------------------------------------------------------------------
   integer :: i, istr, iend, ioerr
   real(q) :: chgtmp
   character (len=200) :: chgfile, oneline
   logical             :: fexst
   !----------------------------------------------------------------------------
   !
   if ( .not.allocated( atchg ) ) then
      allocate( atchg(natom) )
      atchg = 0.D0
      write(*, 100)
      read(*,'(A)', iostat=ioerr) oneline
      if ( ioerr.eq.0.and.(oneline.eq.'y'.or.oneline.eq.'Y') ) then
         write(*,200)
         read(*,*,iostat=ioerr) chgfile
         if ( ioerr.eq.0 ) then
            inquire( file=chgfile, exist=fexst )
            if ( .not.fexst ) then
               write(*, 250) trim(chgfile)
            else
               open(iotmp, file=chgfile, action='read', iostat=ioerr)
               read(iotmp, *, iostat=ioerr) i, chgtmp
               do while ( ioerr.eq.0 )
                  if ( i.ge.1.and.i.le.natom ) atchg(i) = chgtmp
                  read(iotmp, *, iostat=ioerr) i, chgtmp
               enddo
               close(iotmp)
            endif
         endif
      endif
   endif
100 format(/,10x,"No charge properties assigned yet, would you like to read in the charge on some atoms? (y/n)[n]: ",$)
200 format(  10x,"Please input the file that carries the charge, with two columns each line.",/,&
   &         10x,"The first line is the atom index, and the second line is the charge on that",/,&
   &         10x,"atom. Now please input the file name: ",$)
250 format(/,10x,"File ",A," not found, all charges will be set to be 0!")
   !
   call abc2axis
   if ( .not.cartesian ) call dir2car
   !
   write( ioout, 501 ) 200
   write( ioout, 510 ) trim( title )
   write( ioout, 520 ) ''
   write( ioout, 530 ) 'Reactive force field'
   !
   istr = 0
   iend = 0
   do i = 1, natom
      if ( sum(atrel(:,i)).lt.3 ) then
         if ( istr.eq.0 ) then
            istr = i
            iend = istr
         else
            iend = i
         endif
         if ( i.eq.natom ) write(ioout, 540) istr, iend
      else
         if ( istr.gt.0.and.iend.ge.istr ) then
            write(ioout, 540) istr, iend
            istr = 0
            iend = 0
         endif
      endif
   enddo
   !
   write( ioout, 550 ) latt
   write( ioout, 552 )
   !write( ioout, 560 )
   !
   do i = 1, natom
      write( ioout, 555 ) i, EName( attyp(i) ), atpos(:, i), EName( attyp(i) ), 0, 0, atchg(i)
   enddo
   write( ioout, 570 )
   !
500 format( "BIOGRF", i4 )
501 format( "XTLGRF", i4 )
510 format( "DESCRP", 1x, a )
520 format( "REMARK", 1x, a )
530 format( "FFIELD", 1x, a )
540 format( "FIXATOMS", 2i6 )
550 format( "CRYSTX", 1X, 6f11.5 )
552 format( "FORMAT ATOM (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5,i2,i4,f10.5)" )
555 format( "HETATM", 1x, i5, 1x, a5,1x, 3x,1x,1x,1x,5x,3f10.5,1x,a5,i3,i2,1x,f8.5," 0 0 0." )
560 format( "#CELLS", 23x, "6i5" )
570 format( "END" )
   !
return
end subroutine
