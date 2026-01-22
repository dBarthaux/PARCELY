module ODEquation

contains

!==========================================================================================================================

subroutine parcelydydt(me, neq, t, y, ydot)
    
use dvode_module
use dvode_kinds_module,     only: dvode_wp
use ieee_arithmetic,        only: ieee_is_nan
use DropFuncs,              only: DropGrowthRate, DropletSurfaceTension
use EnvEquations,           only: PressureTime, TemperatureTime, SaturationTime
use SolFuncs,               only: SoluteUpdate, Cocondense
use EnvironmentConstants,   only: Na, Mw, pi, rhoL
use ModelParameters,        only: INCLUDE_ORGANICS, SoluteProperties, ndrop, norg, DropSurfTens, OrganicProperties

implicit none

class(dvode_t), intent(inout)       :: me       ! Solver instance
integer                             :: i, neq   !-, Number of equations
real(dvode_wp)                      :: t, y(neq), ydot(neq), Tmp, Prs, RH, WetGrowthSums
real(dvode_wp), dimension(ndrop)    :: WaterMass
real(dvode_wp), allocatable         :: yorg(:,:), OrganicMass(:,:), dOrgdt(:,:)

allocate(yorg(ndrop,norg), OrganicMass(ndrop,norg), dOrgdt(ndrop,norg))

! Last 3 entries are temperature, pressure, and RH
Tmp  = y(neq-2)
Prs  = y(neq-1)
RH   = y(neq)
! ydot is the array of values for dy/dt
ydot = 0.0

if (INCLUDE_ORGANICS) then
    yorg = RESHAPE(y(ndrop+1:neq-3), [ndrop, norg])
    ! Extract the masses of each organic
    do i = 1, norg
        OrganicMass(:,i) = yorg(:,i)
    end do
end if

! Sum of growth rates for the environmental equations
WetGrowthSums = 0.0
! Get water mass of each droplet, kg
WaterMass = (4*pi/3)*rhoL*(y(1:ndrop)**3.0 - SoluteProperties(:,2)**3.0)
! Update droplet surface tension
call DropletSurfaceTension(Tmp, WaterMass, OrganicMass, y(1:ndrop))

if (INCLUDE_ORGANICS) then
    ! Update mixture properties (kappa, density, molar mass)
    call SoluteUpdate(OrganicMass)
    ! Calculate dm/dt of organics
    call Cocondense(OrganicMass, Tmp, WaterMass, y(1:ndrop), dOrgdt)
    ydot(ndrop+1:neq-3) = RESHAPE(dOrgdt, [ndrop*norg])
end if
    
! Loop to calculate dr/dt of water
do i = 1, ndrop
    ydot(i) = DropGrowthRate(Tmp, Prs, RH, y(i), SoluteProperties(i,2), SoluteProperties(i,5), DropSurfTens(i))
    WetGrowthSums = WetGrowthSums + ydot(i)*4*pi*rhoL*y(i)**2
end do
    
! Compute atmospheric variable rates
ydot(neq-2) = TemperatureTime(Tmp, Prs, RH, WetGrowthSums)
ydot(neq-1) = PressureTime(Tmp, Prs, RH)
ydot(neq)   = SaturationTime(Tmp, Prs, RH, WetGrowthSums)

end subroutine parcelydydt

end module ODEquation