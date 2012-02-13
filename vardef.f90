!*******************************************************************************
!***         Definition of common variables.                                 ***
!*******************************************************************************
module cell
use prec
implicit none
   !----------------------------------------------------------------------------
   integer :: natom = 0, ntype = 0, postyp, outtyp
   real(q) :: a, b, c, alpha, beta, gamma, latt(6), cosbc, cosab, cosac
   real(q) :: alat, axis(3,3), rad2ang, ang2rad
   logical :: cartesian = .false.
   equivalence (latt(1), a), (latt(2), b), (latt(3), c), (latt(4), alpha),     &
   &           (latt(5), beta), (latt(6), gamma)
   !----------------------------------------------------------------------------
   real(q), allocatable :: atpos(:,:), atchg(:)
   integer, allocatable :: atrel(:,:), attyp(:), ntm(:)
   !----------------------------------------------------------------------------
   character ( len=80 ) :: title, infile, outfile, info
   character ( len=20 ) :: subname
   character ( len=2  ) :: ETmp
   !----------------------------------------------------------------------------
   character (len=2), allocatable :: EName(:)
   !----------------------------------------------------------------------------
   ! Used to identify atomic type
   integer, parameter :: NMax = 10
   integer            :: iadum(NMax) = 0, idum = 0, typeID(NMax)
   character (len=2)  :: ERead(NMax)
   !----------------------------------------------------------------------------
   ! Default input/output file names:
   character (len=80) :: defaultIns(11) = (/ 'POSCAR        ', 'atomcfg       ',&
   &                     'atomcfg.xyz   ',  'dynmat.icfg   ', 'fort.90       ',&
   &                     'moldyn.vel    ',  'noname        ', 'data.in       ',&
   &                     'dump.lammpstrj',  'STRUCT_IN     ', 'posinp.xyz    '/),  &
   &                     defaultOut(13)= (/ 'POSCAR        ', 'atomcfg       ',&
   &                     'atomcfg.xyz   ',  'dynmat.icfg   ', 'geo           ',&
   &                     'moldyn.vel    ',  'atomcfg.car   ', 'atomcfg.res   ',&
   &                     'data.pos      ',  'data.out      ', 'dump.lammpstrj',&
   &                     'STRUCT_IN     ',  'posinp.xyz    ' /)
   !----------------------------------------------------------------------------
   integer            :: Ninpostypes, Noutpostypes
   !----------------------------------------------------------------------------
end module
