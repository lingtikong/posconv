!*******************************************************************************
!*** Subroutines to read/write XYZ format position file                      ***
!*******************************************************************************
subroutine readart
use prec
use cell
use iounits
implicit none
   !-----------------------------------------------------------------
   integer            :: ioerr, i, j, runid, MaxAtom
   real(q)            :: radum(3), vol, eng
   character (len=512):: input
   character (len=20 ):: strtmp
   integer, external  :: typescreen
   !-----------------------------------------------------------------
   real(q),allocatable:: postmp(:,:)
   integer,allocatable:: typtmp(:)
   !-----------------------------------------------------------------
   subname = 'readart'
   !
   ntype = 0
   ! Read run id
   read( ioin, *, iostat=ioerr ) strtmp, runid
   call error( subname, info, ioerr )
   ! Read energy
   read( ioin, *, iostat=ioerr ) strtmp, eng
   call error( subname, info, ioerr )
   ! write title
   write(title,'("ART run_id:", I5, " Reference_energy:", F15.10)') runid, eng
   do i = 1, 3
   do j = 1, 3
      axis(i,j) = 0.D0
   enddo
   enddo
   ! read box size
   alat = 1.D0
   read( ioin, '(A)', iostat=ioerr ) input
   read(input, *, iostat=ioerr) strtmp, axis(1,1), axis(2,2), axis(3,3), axis(2,1), axis(3,1), axis(3,2)
   if (ioerr.ne.0) read(input, *, iostat=ioerr ) strtmp, axis(1,1), axis(2,2), axis(3,3)
   !
   ! read atoms
   MaxAtom = 2000
   allocate( atpos(3, MaxAtom), attyp(MaxAtom) )
    
   natom = 0
   do while (.true.)
      read( ioin, *, iostat=ioerr) ETmp, radum
      if (ioerr.ne.0) exit

      if (natom >= MaxAtom) then
         allocate(postmp(3, MaxAtom), typtmp(MaxAtom))
         postmp = atpos
         typtmp = attyp
         MaxAtom = MaxAtom + 2000
         deallocate(atpos, attyp)
         allocate( atpos(3, MaxAtom), attyp(MaxAtom) )
         atpos(:,1:natom) = postmp
         attyp(1:natom) = typtmp
         deallocate(postmp, typtmp)
      endif
      natom = natom + 1
      atpos(:,natom) = radum
      attyp(natom) = typescreen(ETmp)
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
   alat  = 1.D0
   !
   allocate( atrel(3,natom) )
   atrel = 1
   !
   cartesian = .true.
   call axis2abc()
   !
return
end subroutine
!
subroutine writeart
use prec
use cell
use iounits
implicit none
   !----------------------------------------------------------------------------
   integer :: i
   !----------------------------------------------------------------------------
   !
   if ( .not.cartesian ) call dir2car
   !
   write( ioout, 100 )
   write( ioout, 200 )
   write( ioout, 350 ) axis(1,1), axis(2,2), axis(3,3), axis(2,1), axis(3,1), axis(3,2)
   !
   do i = 1, natom
      write( ioout, 300 ) attyp(i), atpos(:, i)
   enddo
   !
100 format( "run_id: 0000" )
200 format( "total_energy: 0." )
300 format( I2, 3(1X, F20.15)  )
350 format( "P", 6(1X, F20.15) )
   !
return
end subroutine
