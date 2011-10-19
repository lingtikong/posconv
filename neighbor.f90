! subroutine to calculate the neighbor list accordance to distance
subroutine neighcal()
use prec
use cell,    only: atpos, attyp, natom, axis, alat, cartesian
use iounits, only: iotmp
implicit none
   !----------------------------------------------------------------------------
   integer              :: NbMax = 40, i, j, JJ, ioerr, nbonds, ibond, idum
   integer              :: pairi, pairj, ip, jp
   logical              :: allpairs
   integer, allocatable :: nblist(:,:), bonds(:,:), nbfind(:)
   real(q)              :: RI(3), RJ(3), XIJ(3,1), RIJ(3), rsq, rcutsq(2)
   equivalence (XIJ, RIJ)
   integer              :: outkin = 2, bondtyp = 1
   character (len=100)  :: input, fmtstr='(I?,??(1x,I?))', outname='nblist.dat'
   integer, save        :: iprev = 0
   !----------------------------------------------------------------------------
   if ( allocated(nblist) ) deallocate( nblist )
   if ( allocated(nbfind) ) deallocate( nbfind )
   allocate( nblist(NbMax,natom), nbfind(natom) )
   nbfind = 0
   nblist = 0
   !
   if (cartesian) call car2dir()
   !
   write(*,'(/,10x,"Please input the lower and upper bounds of bond length, if only")')
   write(*,'(10x,"one number is provided, the lowere bound is set to 0 : ", $)')
   read (*, '(A)', iostat = ioerr ) input
   if ( ioerr.ne.0.or.input.eq.'' ) return
   read(input, *, iostat=ioerr) rcutsq
   if (ioerr.ne.0) then
      rcutsq(1) = 0.D0
      read(input, *, iostat=ioerr) rcutsq(2)
      if (ioerr.ne.0) return
   endif
   rcutsq = rcutsq * rcutsq
   !
   pairi = 0
   pairj = 0
   write(*,'(/,10x,"Please input the type pair to count neighbors, 0 for all types. [0 0]:",$)')
   read(*,'(A)',iostat=ioerr) input
   if (ioerr.eq.0.and.input.ne.'') then
      read(input,*,iostat=ioerr) pairi, pairj
      if (ioerr.ne.0) then
         pairi = 0
         pairj = 0
      endif
   endif
   !
   allpairs = pairi.eq.0.and.pairj.eq.0
   !
   do i = 1, natom-1
      RI = atpos(:,i)
      ip = attyp(i)
      do j = i+1, natom
         jp = attyp(j)

         if ( allpairs.or.(pairi.eq.0.and.(ip.eq.pairj.or.jp.eq.pairj))    &
            .or.(pairj.eq.0.and.(ip.eq.pairi.or.jp.eq.pairi)).or.          &
            (ip.eq.pairi.and.jp.eq.pairj).or.(jp.eq.pairi.and.ip.eq.pairj) ) then
            RIJ = atpos(:,j) - RI
            RIJ = RIJ - NINT(RIJ)
            XIJ = matmul(transpose(axis), XIJ)*alat
            !
            rsq = RIJ(1)*RIJ(1) + RIJ(2)*RIJ(2) + RIJ(3)*RIJ(3)
            if ( rsq.gt.rcutsq(1).and.rsq.le.rcutsq(2) ) then
               nbfind(i) = nbfind(i) + 1
               nbfind(j) = nbfind(j) + 1
               if ( nbfind(i).gt.NbMax.or.nbfind(j).gt.NbMax ) then
                  write(*,'(/,10x,"Too many neighbors found, please increase NbMax and recompile the code!")')
                  return
               endif
               nblist( nbfind(i), i) = j
               nblist( nbfind(j), j) = i
            endif
         endif
      enddo
   enddo
   !
   write(*,'(/,10x,"Please indicate if you want to output (1) neighbor list or (2) bonds list [2]: ", $)')
   read (*, '(A)', iostat=ioerr ) input
   if ( ioerr.eq.0.and.input.ne.'' ) then
      read(input, *, iostat=ioerr) outkin
      if (ioerr.ne.0) outkin = 2
   endif
   write(*,'(10x,"Please input the output file name [nblist.dat]: ", $ )')
   read (*, '(A)', iostat=ioerr ) input
   if ( ioerr.eq.0.and.input.ne.'' ) outname = trim(input)
   open(iotmp, file=outname, status='unknown', action='write')
   !
   if ( outkin.eq.1 ) then
      idum = int(log10(dble(natom))) + 1
      write(fmtstr(3:3),  '(I1)') idum
      write(fmtstr(12:12),'(I1)') idum
      idum = maxval(nbfind)
      write(fmtstr(5:6), '(I2.2)') idum
      !
      write(iotmp,'("# index  neighbors")')
      do i = 1, natom
         write(iotmp, fmtstr) nblist(1:idum,i)
      enddo
   else
      nbonds = count(nblist.gt.0)/2
      if ( nbonds.lt.1 ) then
         write(*, '(/,10x,"No neighbors found!")')
         return
      endif
      if ( allocated(bonds) ) deallocate(bonds)
      allocate( bonds(2,nbonds) )
      ibond = 0
      do i = 1, natom
         do JJ = 1, nbfind(i)
            j  = nblist(JJ,i)
            if (j.gt.i) then
               ibond = ibond + 1
               bonds(:,ibond) = (/i, j/)
            endif
         enddo
      enddo
      !
      write(*,'(10x,"Please assign a bond type to the bonds found [1]: ", $ )')
      read (*, '(A)', iostat=ioerr) input
      if (ioerr.eq.0.and.input.ne.'') then
         read(input, *, iostat=ioerr) bondtyp
         if (ioerr.ne.0) bondtyp = 1
      endif
      ! reset index
      write(*,'(10x,"Previous bond index is ",I10,", reset? (y/n, num)[n]: ", $)') iprev
      read (*,'(A)', iostat=ioerr) input
      if (ioerr.eq.0) then
         if (input.eq.'y'.or.input.eq.'Y') then
            iprev = 0
         else if (input.ne.'n'.and.input.ne.'N') then
            read(input,*,iostat=ioerr) idum
            if (ioerr.eq.0) iprev = idum
         endif
      endif

      write(fmtstr,'("(I",I2.2,",1x,I2,2(1x,I",I2.2,"))")') int(log10(dble(nbonds+iprev)))+1, int(log10(dble(natom)))+1
      write(iotmp,'(/,"Bonds",/)')
      do i = 1, nbonds
         write(iotmp, fmtstr) i+iprev, bondtyp, bonds(:,i)
      enddo
      iprev = iprev + nbonds
   endif
   close(iotmp)
   !
return
end subroutine
