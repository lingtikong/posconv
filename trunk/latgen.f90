!******************************************************************************
!***    Subroutine to write out atomic configurations for latgen            ***
!******************************************************************************
subroutine latgen
use prec
use cell
use iounits
implicit none
   !----------------------------------------------------------------------------
   integer  :: i, j
   character(len=100 ) :: fmtstr1 = '("attyp[",I??,"] = ", I2,";")'
   character(len=100 ) :: fmtstr2 = '("atpos[",I??,"][",I1,"] = ", F16.10,";")'
   !----------------------------------------------------------------------------
   if ( cartesian ) call car2dir
   call reduceaxis
   !
   write( ioout, '("nucell = ", I4,";")') natom
   write( ioout, '("ntype  = ", I2,";")') ntype
   write( ioout, '(/,"// alat  = ", F8.4,";")') alat
   do i = 1, 3
      write( ioout, *)
      do j = 1, 3
         if (axis(i,j)*axis(i,j).gt.1.e-8) then
            write( ioout, '("latvec[",I1,"][",I1,"] = ", F16.10,";")') i-1, j-1, axis(i,j)
         endif
      enddo
   enddo
   write( ioout, * )
   write( ioout, '("atpos = memory->create(atpos,nucell,3,",A,"atpos",A,");")') char(34),char(34)
   write( ioout, '("attyp = memory->create(attyp,nucell,",A,"attyp",A,");")') char(34),char(34)
   !
   write( fmtstr1(12:13),'(I2.2)') int(log10(natom-0.1))+1
   write( fmtstr2(12:13),'(I2.2)') int(log10(natom-0.1))+1
   !
   do i = 1, natom
      write( ioout, *)
      write( ioout, fmtstr1) i-1, attyp(i)
      do j = 1, 3
         write( ioout, fmtstr2) i-1, j-1, atpos(j,i)
      enddo
   enddo
   write( ioout, *)
   write( ioout, '("initialized = 1;",/,"//break;")')
   !
return
end subroutine
