module ModelParameters
    
use dvode_kinds_module, only: wp => dvode_wp
    
implicit none

logical :: INCLUDE_ORGANICS, INCLUDE_COCONDENSE, INCLUDE_BAT
real(wp) :: W, ac, at, sftc, TmpCurrent
integer ndrop, ninorg, norg, Msft, Cubes
real(wp), allocatable :: InorganicBase(:,:), OrganicProperties(:,:), SoluteProperties(:,:), DropSurfTens(:)

end module ModelParameters
