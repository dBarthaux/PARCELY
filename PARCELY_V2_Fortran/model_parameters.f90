module ModelParameters
    
use dvode_kinds_module, only: dvode_wp
    
implicit none

logical :: INCLUDE_ORGANICS, INCLUDE_COCONDENSE
real(dvode_wp) :: W, ac, at, sftc, TmpCurrent
integer ndrop, ninorg, norg, Msft, Cubes
real(dvode_wp), allocatable :: InorganicBase(:,:), OrganicProperties(:,:), SoluteProperties(:,:), DropSurfTens(:)

end module ModelParameters