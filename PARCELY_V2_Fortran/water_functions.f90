module WaterFuncs

use dvode_kinds_module, only: dvode_wp
use EnvironmentConstants, only: cpl, cpi, Rv

contains

!==========================================================================================================================

!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Function to calculate the pure component surface tension of water using the        *
!*   equation from Kalova and Mares (2018). https://doi.org/10.1063/1.5081640           *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Dan Barthaux                                                                       *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2025/11/04                                                      *
!*                                                                                      *
!****************************************************************************************

function WaterSurfaceTension(T) result(Sigma)

implicit none

real(dvode_wp), intent(in)  :: T        ! Ambient temperature of air, Kelvin
real(dvode_wp)              :: CritT    ! Critical temperature of water, Kelvin
real(dvode_wp)              :: Tau      ! Dimensionless variable
real(dvode_wp)              :: Sigma    ! Surface tension, J/m2

CritT = 647.15
Tau = 1 - T/CritT
Sigma = 241.322*(Tau**1.26)*(1 - 0.0589*(Tau**0.5) - 0.56917*Tau)*1e-3

end function WaterSurfaceTension

!==========================================================================================================================
    
!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Function to calculate the latent heat of condensation of water from Ambaum (2020). *
!*   https://doi.org/10.1002/qj.3899                                                    *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Dan Barthaux                                                                       *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2025/11/04                                                      *
!*                                                                                      *
!****************************************************************************************
    
function LatentHeatEvap(T) result(LHEvap)

implicit none

real(dvode_wp), intent(in)  :: T        ! Ambient temperature of air, Kelvin
real(dvode_wp)              :: T0       ! Temperature of water at its triple-point, Kelvin
real(dvode_wp)              :: cpv      ! Isobaric mass heat capacity of water vapor, J/kg/K
real(dvode_wp)              :: L0       ! Latent heat of vaporization at triple point, J/kg
real(dvode_wp)              :: LHEvap   ! Latent heat of vaporization, J/kg

T0  = 273.16
cpv = 2040
L0  = 2.501e6
LHEvap = L0 - (cpl - cpv)*(T - T0)

end function LatentHeatEvap

!==========================================================================================================================

!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Function to calculate the latent heat of sublimation of water from Ambaum (2020).  *
!*   https://doi.org/10.1002/qj.3899                                                    *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Dan Barthaux                                                                       *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2025/11/04                                                      *
!*                                                                                      *
!****************************************************************************************

function LatentHeatSub(T) result(LHSub)

implicit none

real(dvode_wp), intent(in)  :: T        ! Ambient temperature of air, Kelvin
real(dvode_wp)              :: T0       ! Temperature of water at its triple-point, Kelvin
real(dvode_wp)              :: cpv      ! Isobaric mass heat capacity of water vapor, J/kg/K
real(dvode_wp)              :: L0       ! Latent heat of vaporization at triple point, J/kg
real(dvode_wp)              :: LHSub    ! Latent heat of sublimation, J/kg

T0  = 273.16
cpv = 1885
L0  = 2.260e6
LHSub = L0 - (cpi - cpv)*(T - T0)

end function LatentHeatSub
    
!==========================================================================================================================

!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Function to calculate the saturation vapor pressure of water from Ambaum (2020).   *
!*   https://doi.org/10.1002/qj.3899                                                    *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Dan Barthaux                                                                       *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2025/11/04                                                      *
!*                                                                                      *
!****************************************************************************************

function SatVapWater(T, Phase) result(SVPWater)

implicit none

real(dvode_wp), intent(in)  :: T            ! Ambient temperature of air, Kelvin
character(3), intent(in)    :: Phase        ! Phase of water, "ice" or "liq" for liquid
real(dvode_wp)              :: T0           ! Temperature of water at its triple-point, Kelvin
real(dvode_wp)              :: es0          ! Saturation vapor pressure at triple point, Pa
real(dvode_wp)              :: cpvl         ! Isobaric mass heat capacity of vapor [liquid SVP], J/kg/K
real(dvode_wp)              :: cpvi         ! Isobaric mass heat capacity of vapor [ice SVP], J/kg/K
real(dvode_wp)              :: L0l          ! Latent heat of vaporization at triple point, J/kg
real(dvode_wp)              :: L0i          ! Latent heat of sublimation at triple point, J/kg
real(dvode_wp)              :: Ll           ! Latent heat of vaporization at T, J/kg
real(dvode_wp)              :: Li           ! Latent heat of sublimation at T, J/kg
real(dvode_wp)              :: SVPWater     ! Saturation vapor pressure of water vapor, Pa

T0 = 273.16
es0 = 611.655
cpvl = 2040
cpvi = 1885
L0l = 2.501e6
L0i = 2.260e6

if (Phase == 'liq') then
    Ll = LatentHeatEvap(T)
    SVPWater = es0*(T0/T)**((cpl - cpvl)/Rv)*EXP(L0l/(Rv*T0) - Ll/(Rv*T))
else
    Li = LatentHeatSub(T)
    SVPWater = es0*(T0/T)**((cpi - cpvi)/Rv)*EXP(L0i/(Rv*T0) - Li/(Rv*T))
end if
    
end function SatVapWater

!==========================================================================================================================
    
end module WaterFuncs