module WaterFuncs

use dvode_kinds_module, only: wp => dvode_wp
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

pure elemental function WaterSurfaceTension(T) result(Sigma)

implicit none

real(wp), intent(in)  :: T        ! Ambient temperature of air, Kelvin
real(wp)              :: CritT    ! Critical temperature of water, Kelvin
real(wp)              :: Tau      ! Dimensionless variable
real(wp)              :: Sigma    ! Surface tension, J/m2

CritT = 647.15_wp
Tau = 1_wp - T/CritT
Sigma = 241.322_wp*(Tau**1.26_wp)*(1_wp - 0.0589_wp*(Tau**0.5_wp) - 0.56917_wp*Tau)*1e-3_wp

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
    
pure elemental function LatentHeatEvap(T) result(LHEvap)

implicit none

real(wp), intent(in)  :: T        ! Ambient temperature of air, Kelvin
real(wp)              :: T0       ! Temperature of water at its triple-point, Kelvin
real(wp)              :: cpv      ! Isobaric mass heat capacity of water vapor, J/kg/K
real(wp)              :: L0       ! Latent heat of vaporization at triple point, J/kg
real(wp)              :: LHEvap   ! Latent heat of vaporization, J/kg

T0  = 273.16_wp
cpv = 2040_wp
L0  = 2.501e6_wp
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

pure elemental function LatentHeatSub(T) result(LHSub)

implicit none

real(wp), intent(in)  :: T        ! Ambient temperature of air, Kelvin
real(wp)              :: T0       ! Temperature of water at its triple-point, Kelvin
real(wp)              :: cpv      ! Isobaric mass heat capacity of water vapor, J/kg/K
real(wp)              :: L0       ! Latent heat of vaporization at triple point, J/kg
real(wp)              :: LHSub    ! Latent heat of sublimation, J/kg

T0  = 273.16_wp
cpv = 1885_wp
L0  = 2.260e6_wp
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

pure elemental function SatVapWater(T, Phase) result(SVPWater)

implicit none

real(wp), intent(in)        :: T            ! Ambient temperature of air, Kelvin
character(3), intent(in)    :: Phase        ! Phase of water, "ice" or "liq" for liquid
real(wp)                    :: T0           ! Temperature of water at its triple-point, Kelvin
real(wp)                    :: es0          ! Saturation vapor pressure at triple point, Pa
real(wp)                    :: cpvl         ! Isobaric mass heat capacity of vapor [liquid SVP], J/kg/K
real(wp)                    :: cpvi         ! Isobaric mass heat capacity of vapor [ice SVP], J/kg/K
real(wp)                    :: L0l          ! Latent heat of vaporization at triple point, J/kg
real(wp)                    :: L0i          ! Latent heat of sublimation at triple point, J/kg
real(wp)                    :: Ll           ! Latent heat of vaporization at T, J/kg
real(wp)                    :: Li           ! Latent heat of sublimation at T, J/kg
real(wp)                    :: SVPWater     ! Saturation vapor pressure of water vapor, Pa

T0 = 273.16_wp
es0 = 611.655_wp
cpvl = 2040_wp
cpvi = 1885_wp
L0l = 2.501e6_wp
L0i = 2.260e6_wp

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