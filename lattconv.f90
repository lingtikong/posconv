!******************************************************************************
!*** This subroutine transfer the coordinate from cartesian into fractional ***
!****************************************************************************** 
subroutine car2dir
use prec
use cell, only: axis, atpos, alat, cartesian
implicit none
   !
   real(q)  :: invaxis(3, 3)
   !------------------------------------------------------------------
   if ( .not.cartesian ) return
   call matinv( 3, axis, invaxis )
   invaxis   = invaxis / alat
   atpos     = matmul(transpose(invaxis), atpos)
   cartesian = .false.
   !
return
end subroutine
!
!******************************************************************************
!*** This subroutine transfer the coordinate from fractional into cartesian ***
!******************************************************************************
subroutine dir2car
use prec
use cell, only: axis, atpos, alat, cartesian
implicit none
   !------------------------------------------------------------------
   if (cartesian) return
   atpos     = matmul(transpose(axis), atpos) * alat
   cartesian = .true.
return
end subroutine
!
subroutine vec_cross( vec1, vec2, vec)
use prec
!-----------------------------------------------------------------
! Subroutine to calculate the cross product
! of two dimensional 3 vectors
!
! vec1, vec2 : (in ) vectors, double precision
! vec        : (out) vector,  double precision
!-----------------------------------------------------------------
implicit none
   !--------------------------------------------------------------
   real(q), intent(in) :: vec1(3), vec2(3)
   real(q), intent(out):: vec(3)
   !--------------------------------------------------------------
   vec(1) = vec1(2) * vec2(3) - vec1(3) * vec2(2)
   vec(2) = vec1(3) * vec2(1) - vec1(1) * vec2(3)
   vec(3) = vec1(1) * vec2(2) - vec1(2) * vec2(1)
   !
return
end subroutine
!
!*****************************************************************
!
subroutine volume( axis, vol )
use prec
!-----------------------------------------------------------------
! Subroutine to calculate the volume of a 3x3 matrix
!-----------------------------------------------------------------
implicit none
   !--------------------------------------------------------------
   real(q), intent(in) :: axis(3, 3)
   real(q), intent(out):: vol
   !--------------------------------------------------------------
   real(q)  :: ratmp(3)
   !--------------------------------------------------------------
   call vec_cross( axis(1,:), axis(2,:), ratmp )
   vol  = sum( ratmp * axis(3,:) )
return
end subroutine

subroutine abc2axis( )
use cell
implicit none
   !--------------------------------------------------------------
   cosbc = cos( alpha * ang2rad )
   cosac = cos( beta  * ang2rad )
   cosab = cos( gamma * ang2rad )
   !
   alat = a
   axis = 0.D0
   !axis(1, 1) =  a * sin( gamma * ang2rad )
   !axis(1, 2) =  a * cos( gamma * ang2rad )
   axis(1, 1) = a
   axis(2, 1) = b * cosab
   axis(2, 2) = b * sin( gamma * ang2rad )
   !
   axis(3, 1) = c * cosac
   !axis(3, 2) = c * ( cosbc - cosac* cosab  ) / sin( gamma * ang2rad )
   axis(3, 2) = b*c*(cosbc - cosab*cosac)/axis(2,2)
   axis(3, 3) = sqrt( c*c - axis(3,1) * axis(3,1) - axis(3,2) * axis(3,2) )
   !
   axis       = axis / a
   !-------------------------------------------------------------
return
end subroutine
!
subroutine reduceaxis()
use cell
implicit none
   !--------------------------------------------------------------
   real(q) :: la, rla
   integer :: i, j
   !--------------------------------------------------------------
   la  = sqrt(sum(axis(1,:)*axis(1,:)))
   rla = 1.D0/la
   
   alat = alat * la
   do i = 1, 3
   do j = 1, 3
      axis(i,j) = axis(i,j) * rla
   enddo
   enddo
return
end subroutine
!
subroutine axis2abc( )
use cell
implicit none
   !--------------------------------------------------------------
   a       = sqrt( sum( axis(1,:)*axis(1,:) ) )
   b       = sqrt( sum( axis(2,:)*axis(2,:) ) )
   c       = sqrt( sum( axis(3,:)*axis(3,:) ) )
   cosbc   = sum( axis(2,:)*axis(3,:) ) / ( b*c )
   cosab   = sum( axis(1,:)*axis(2,:) ) / ( a*b )
   cosac   = sum( axis(1,:)*axis(3,:) ) / ( a*c )
   alpha   = acos( cosbc )
   beta    = acos( cosac )
   gamma   = acos( cosab )
   latt(1:3) = latt(1:3) * alat
   latt(4:6) = latt(4:6) * rad2ang
   !--------------------------------------------------------------
return
end subroutine

subroutine matinv(n, MatA, Mat)
use prec
implicit none
   !
   integer, intent(in)  :: n
   real(q), intent(in)  :: MatA(n,n)
   real(q), intent(out) :: Mat(n,n)
   !
   integer :: i,icol,irow,j,k,l,ll,idr,idc,rl,cl
   integer :: indxc(n),indxr(n),ipiv(n);
   real(q) :: big, dum, pivinv

   Mat = MatA
   ipiv = 0

   do i = 1, n
      big = 0.D0
      do j = 1, n
         if ( ipiv(j).ne.1 ) then
            do k=1, n
               if (ipiv(k).eq.0) then
                  if (abs(Mat(j,k)).ge.big) then
                     big = abs(Mat(j,k))
                     irow = j
                     icol = k
                  endif
               else if (ipiv(k) > 1) then

               endif
            enddo
         endif
      enddo
      ipiv(icol) = ipiv(icol) + 1
      if ( irow.ne.icol) then
         do l=1, n
            dum = Mat(irow,l)
            Mat(irow,l) = Mat(icol,l)
            Mat(icol,l) = dum
         enddo
      endif
      indxr(i) = irow
      indxc(i) = icol
      if (Mat(icol,icol).eq.0.) then

      endif
      pivinv = 1./Mat(icol,icol)
      Mat(icol,icol) = 1.D0
      do l=1, n
         Mat(icol,l) = Mat(icol,l) * pivinv;
      enddo
      do ll=1, n
         if (ll.ne.icol) then
            dum = Mat(ll, icol)
            Mat(ll,icol) = 0.
            do l=1, n
               Mat(ll,l) = Mat(ll,l) - Mat(icol,l)*dum
            enddo
         endif
      enddo
   enddo

   do l=n, 0, -1
      rl = indxr(l)
      cl = indxc(l)
      if ( rl.ne.cl ) then
         do k = 1, n
            dum = Mat(k, rl)
            Mat(k,rl) = Mat(k, cl)
            Mat(k,cl) = dum
         enddo
      endif
   enddo
end subroutine
