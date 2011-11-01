integer function typescreen( Ein )
use cell, only: ntype, Eread, NMax
implicit none
   !----------------------------------------------------------------------------
   character ( len=2 ), intent(in) :: Ein
   !----------------------------------------------------------------------------
   integer            :: i, j
   character (len=80) :: errinfo
   character (len=10) :: sname = 'typescreen'
   !----------------------------------------------------------------------------
   j  = 0
   do i = 1, ntype
      if ( Ein.eq.Eread(i) ) then
         j = i
         exit
      endif
   enddo
   !
   if ( j.eq.0 ) then
      ntype = ntype + 1
      j     = ntype
      if ( ntype.gt.NMax ) then
         errinfo = 'Too many kinds of atoms found!'
         call error( sname, errinfo, ntype )
      endif
      Eread(ntype)  = Ein
   endif
   typescreen = j
return
end function typescreen

integer function typeindex( Ein )
use cell, only: ntype, Eread
implicit none
   !--------------------------------------------------------------
   character ( len=2 ), intent(in) :: Ein
   !--------------------------------------------------------------
   integer            :: i, j
   character (len=80) :: errinfo
   character (len=10) :: sname = 'typeindex'
   !--------------------------------------------------------------
   typeindex  = 0
   !
   do i = 1, ntype
      if ( Ein.eq.Eread(i) ) then
         typeindex = i
         exit
      endif
   enddo
   !
   if ( typeindex.eq.0 ) then
      errinfo = 'Error while identifying element!'
      call error( sname, errinfo, 1 )
   endif
   !
return
end function typeindex

subroutine DisplayTypeRead()
use cell, only: ntype, Eread, typeID
implicit none
   !--------------------------------------------------------------
   integer  :: i
   !--------------------------------------------------------------
   write(*, '(/,10x,70("="))')
   write(*, '(10x,"Comparison of atom types read and assigned:",/,10x,70("-"))')
   write(*, '(12x,"Assigned type  ID: ", $)')
   do i = 1, ntype
      write(*, '(I4,$)') typeID(i)
   enddo
   write(*, '(/,12x,"Read type name/ID: ", $)')
   do i = 1, ntype
      write(*, '(A4,$)') trim(Eread(i))
   enddo
   write(*, '(/,10x,70("="),/)')
   !
return
end subroutine
