!*******************************************************************************
! Subroutine to calculate the nominal radius of a cluster
! No periodic boundary condition is applied.
!*******************************************************************************
subroutine calnominal
use prec
use cell
implicit none
   !----------------------------------------------------------------------------
   real(q) :: Rsum, Rsum2, poscen(3), Rnom, Rsdv
   real(q) :: RIJ(3), rsq, r
   integer :: i, j
   !----------------------------------------------------------------------------
   if ( .not.cartesian ) call dir2car
   !
   forall(i=1:3) poscen(i) = sum(atpos(i,:))/dble(natom)
   !
   Rsum  = 0.D0
   Rsum2 = 0.D0
   !
   do i = 1, natom
      RIJ = atpos(:,i) - poscen
      rsq = sum( RIJ * RIJ )
      r   = dsqrt( rsq )
      !
      Rsum2 = Rsum2 + rsq
      Rsum  = Rsum  + r
   enddo
   !
   Rnom = Rsum / dble(natom)
   Rsdv = dsqrt( ( Rsum2 - dble(natom)*Rnom*Rnom )/ dble(natom) )
   !
   write(*, '(/, 10x, "The nominal radius of the cluster is: ", F12.8, " +/-", F12.8 )' ) Rnom, Rsdv
   !
return
end subroutine
