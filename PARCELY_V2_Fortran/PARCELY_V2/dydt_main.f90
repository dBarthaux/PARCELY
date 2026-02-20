module ODEquations

contains

!==========================================================================================================================

subroutine parcelydydt(me, neq, t, y, ydot)
    
use dvode_module
use dvode_kinds_module,     only: wp => dvode_wp
use DropFuncs,              only: DropGrowthRate, DropletSurfaceTension
use EnvEquations,           only: PressureTime, TemperatureTime, SaturationTime
use SolFuncs,               only: SoluteUpdate, Cocondense
use EnvironmentConstants,   only: Na, Mw, pi, rhoL
use ModelParameters,        only: INCLUDE_ORGANICS, INCLUDE_COCONDENSE, SoluteProperties, ndrop, norg, DropSurfTens, OrganicProperties

implicit none

class(dvode_t), intent(inout)       :: me       ! Solver instance
integer                             :: i, neq   !-, Number of equations
real(wp)                            :: t, y(neq), ydot(neq), Tmp, Prs, RH, WetGrowthSums
real(wp), dimension(ndrop)          :: WaterMass
real(wp), allocatable               :: yorgcond(:,:), OrganicGas(:), OrganicCond(:,:), dOrgConddt(:,:), dOrgGasdt(:), WaterRad(:), OriginalGas(:)

allocate(yorgcond(ndrop,norg), OrganicGas(norg), OrganicCond(ndrop,norg), dOrgConddt(ndrop,norg), dOrgGasdt(norg), WaterRad(ndrop), OriginalGas(norg))

! Last 3 entries are temperature, pressure, and RH
Tmp  = y(neq-2)
Prs  = y(neq-1)
RH   = y(neq)
! ydot is the array of values for dy/dt
ydot = 0.0_wp

if (INCLUDE_ORGANICS) then
    yorgcond = RESHAPE(y(ndrop+1:neq-3-norg), shape(yorgcond))
    ! Gas-phase organic, molec/cm3
    OrganicGas = y(neq-3-norg+1:neq-3)
    ! Extract the masses of each organic
    do i = 1, norg
        ! Condensed phase organic, molec
        OrganicCond(:,i) = yorgcond(:,i)
    end do
end if

WaterRad = EXP(y(1:ndrop))
! Sum of growth rates for the environmental equations, kg/s
WetGrowthSums = 0.0_wp
! Get water mass of each droplet, kg
WaterMass = (4_wp*pi/3_wp)*rhoL*(WaterRad**3 - SoluteProperties(:,2)**3)
! Update droplet surface tension
call DropletSurfaceTension(Tmp, WaterMass, OrganicCond, WaterRad)

! Comment this section out and rebuild if you want initial organics but no cocondensation
if ((INCLUDE_ORGANICS) .and. (INCLUDE_COCONDENSE)) then
    ! Update mixture properties (kappa, density, molar mass)
    call SoluteUpdate(OrganicCond)
    ! Calculate dm/dt of organics
    call Cocondense(OrganicCond, OrganicGas, WaterMass, WaterRad, Tmp, dOrgConddt, dOrgGasdt)
    ydot(ndrop+1:neq-3-norg) = RESHAPE(dOrgConddt, [ndrop*norg])
    ! ln rate
    ydot(ndrop+1:neq-3-norg) = ydot(ndrop+1:neq-3-norg)
    ydot(neq-3-norg+1:neq-3) = dOrgGasdt
end if

! Calculate dr/dt of water, m/s
do i = 1, ndrop
    ydot(i) = DropGrowthRate(Tmp, Prs, RH, WaterRad(i), SoluteProperties(i,2), SoluteProperties(i,5), DropSurfTens(i))
end do
WetGrowthSums = SUM(ydot(1:ndrop)*4_wp*pi*rhoL*WaterRad**2)
ydot(1:ndrop) = ydot(1:ndrop)/WaterRad

! Compute atmospheric variable rates
ydot(neq-2) = TemperatureTime(Tmp, Prs, RH, WetGrowthSums)
ydot(neq-1) = PressureTime(Tmp, Prs, RH)
ydot(neq)   = SaturationTime(Tmp, Prs, RH, WetGrowthSums)

end subroutine parcelydydt

!==========================================================================================================================

end module ODEquations
