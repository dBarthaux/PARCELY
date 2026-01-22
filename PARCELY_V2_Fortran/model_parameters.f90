module ModelParameters
    
use dvode_kinds_module, only: dvode_wp
    
implicit none

logical :: INCLUDE_ORGANICS
real(dvode_wp) :: W, ac, at, sftc
integer ndrop, ninorg, norg, Msft
real(dvode_wp), allocatable :: InorganicBase(:,:), OrganicProperties(:,:), SoluteProperties(:,:), DropSurfTens(:)

end module ModelParameters