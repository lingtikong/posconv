!-------------------------------------------------------------------------------
! Subroutine to realign the molecular
!-------------------------------------------------------------------------------
subroutine realign_molecular
use prec
use cell
implicit none
   !----------------------------------------------------------------------------
   integer :: i, ia, ib, ic, ioerr
   real(q) :: X1(3), X2(3), X3(3), R12(3), R13(3), Cross(3), costh, r1, r2, r3
   real(q) :: posold(3), posnew(3), mat1(3,3), mat2(3,3), mat3(3,3), sinth
   !----------------------------------------------------------------------------
   character(len=200) :: input
   !---------------------------------------------------------------------------- 
   subname = 'realign_molecular'
   info    = 'Wrong input!'
   !
   if ( .not.cartesian ) call dir2car( )
   write( *, 10 ) natom
   write( *, 20 ) latt
   write( *, 50 )
   read ( *, *, iostat=ioerr ) ia
   call error( subname, info, ioerr )
   if ( ia.lt.1.or.ia.gt.natom ) then
      write( info, '("Index must not be within 1 to ", I6 )' ) natom
      call error( subname, info, 1 )
   endif
   posold = atpos(:,ia)
   X1     = posold
   !
   write( *, 60 )
   read ( *, *, iostat=ioerr ) ib
   call error( subname, info, ioerr )
   if ( ib.lt.1.or.ib.gt.natom ) then
      write( info, '("Index must not be within 1 to ", I6 )' ) natom
      call error( subname, info, 1 )
   endif
   write( *, 70 )
   read ( *, *, iostat=ioerr ) ic
   call error( subname, info, ioerr )
   if ( ic.lt.1.or.ic.gt.natom ) then
      write( info, '("Index must not be within 1 to ", I6 )' ) natom
      call error( subname, info, 1 )
   endif
   X2 = atpos(:, ib)
   X3 = atpos(:, ic)
   !
   R12 = X2 - X1
   R13 = X3 - X1
   !
   call vec_cross( R12, R13, Cross )
   !
   r1 = sqrt( sum(R12*R12) )
   r2 = sqrt( sum(R13*R13) )
   r3 = sqrt( sum(Cross*Cross) )
   !
   sinth = r3 / ( r1 * r2 )
   costh = sum( R12 * R13 ) / ( r1 * r2 )
   !
   mat3(:,1) = (/ 0.D0, -r3, 0.D0 /)
   mat3(:,2) = (/ r1,  0.D0, 0.D0 /)
   mat3(:,3) = (/ r2*costh, 0.D0, r2*sinth /)
   !
   mat2(:,1) = Cross
   mat2(:,2) = R12
   mat2(:,3) = R13
   !
   call matinv( 3, mat2, mat1 )
   mat2 = matmul( mat3, mat1 )
   !
   forall ( i=1:natom) atpos(:,i) = atpos(:,i) - posold
   atpos = matmul( mat2, atpos )
   !
   posnew = latt(1:3) * 0.5D0
   !
   write( *, 80 ) posnew
   read ( *, '(A)', iostat=ioerr ) input
   call error( subname, info, ioerr )
   if ( input.ne.'' ) then
      read( input, *, iostat=ioerr ) posnew
      call error( subname, info, ioerr )
   endif
   !
   forall ( i=1:natom) atpos(:,i) = atpos(:,i) + posnew
   !
   write( *, 90 )
   read ( *, '(A)', iostat=ioerr ) input
   if ( ioerr.eq.0 ) then
      if ( input.eq.'y'.or.input.eq.'Y' ) then
         call car2dir( )
         atpos = mod(atpos, 1.D0)
         where ( atpos.lt.0.D0 ) atpos = atpos + 1.D0
      endif
   endif
   !
10 format( 10x,"Atom index range: 1-", I6        )
20 format( 10x,"Current box size:  ", 6(F7.3,1x) )
50 format( 10x,"Please input the atom index of the rotation center     : ", $ )
60 format( 10x,"Please input the atom index that will sit on the x-axis: ", $ )
70 format( 10x,"Please input the atom index that will sit in x-z plane : ", $ )
80 format( 10x,"Please input the new coordinate of the rotation center[",3F8.3,"]: ", $ )
90 format( 10x,"Do you want to apply periodic boundary condition ?(y/n)[n]: ", $ )
   !
return
end subroutine
