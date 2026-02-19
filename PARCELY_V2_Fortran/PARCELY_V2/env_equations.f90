module EnvEquations

use dvode_kinds_module, only: wp => dvode_wp
use EnvironmentConstants, only: g, Rd, cp, Rv, eta
use ModelParameters, only: W, Cubes

contains

!==========================================================================================================================

!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Function to calculate the change in ambient pressure over time taken from          *
!*   Rothenberg and Wang (2017) https://doi.org/10.5194/gmd-10-1817-2017                *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Dan Barthaux                                                                       *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2025/11/04                                                      *
!*                                                                                      *
!****************************************************************************************

pure elemental function PressureTime(Tmp, Prs, RH) result (dPdt)

use AirFuncs, only: VirtualTemp

implicit none

real(wp), intent(in)  :: Tmp      ! Temperature, K
real(wp), intent(in)  :: Prs      ! Pressure, Pa
real(wp), intent(in)  :: RH       ! Relative humidity, -
real(wp)              :: Tv       ! Virtual temperature, K
real(wp)              :: dPdt     ! Change in pressure over time, Pa/s

Tv = VirtualTemp(Tmp, Prs, RH)
dPdt = -g*Prs*W/(Rd*Tv)

end function PressureTime

!==========================================================================================================================

!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Function to calculate the change in ambient temperature over time taken from       *
!*   Rothenberg and Wang (2017) https://doi.org/10.5194/gmd-10-1817-2017                *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Dan Barthaux                                                                       *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2025/11/04                                                      *
!*                                                                                      *
!****************************************************************************************

pure elemental function TemperatureTime(Tmp, Prs, RH, WetGrowthSums) result(dTdt)

use WaterFuncs, only: LatentHeatEvap
use AirFuncs, only: DensityAir

implicit none

real(wp), intent(in)  :: Tmp              ! Temperature, K
real(wp), intent(in)  :: Prs              ! Pressure, Pa
real(wp), intent(in)  :: RH               ! Relative humidity, -
real(wp), intent(in)  :: WetGrowthSums    ! Sum of condensational growth rates (water vapor) of the droplets, dm/dt [kg/s]
real(wp)              :: L                ! Latent heat of evaporation (of water vapor), J/kg
real(wp)              :: rhoA             ! Density of dry air, kg/m3
real(wp)              :: dTdt             ! Change in temperature over time, K/s
real(wp)              :: cbrl             ! Cubes in real format

cbrl = REAL(Cubes)
L = LatentHeatEvap(Tmp)
rhoA = DensityAir(Tmp, Prs, RH)
dTdt = -g*W/cp + (L/(cp*rhoA))*WetGrowthSums/(cbrl*1e-6_wp) ! per cubic centimeter factor

end function TemperatureTime

!==========================================================================================================================

!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Function to calculate the change in ambient relative humidity over time taken from *
!*   Rothenberg and Wang (2017) https://doi.org/10.5194/gmd-10-1817-2017                *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Dan Barthaux                                                                       *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2025/11/04                                                      *
!*                                                                                      *
!****************************************************************************************

pure elemental function SaturationTime(Tmp, Prs, RH, WetGrowthSums) result(dRHdt)

use WaterFuncs, only: LatentHeatEvap, SatVapWater
use AirFuncs, only: DensityAir

implicit none

real(wp), intent(in)  :: Tmp              ! Temperature, K
real(wp), intent(in)  :: Prs              ! Pressure, Pa
real(wp), intent(in)  :: RH               ! Relative humidity, -
real(wp), intent(in)  :: WetGrowthSums    ! Sum of condensational growth rates (water vapor) of the droplets, dm/dt [kg/s]
real(wp)              :: L                ! Latent heat of evaporation (of water vapor), J/kg
real(wp)              :: es               ! Saturation vapor pressure of water, Pa
real(wp)              :: rhoA             ! Density of dry air, kg/m3
real(wp)              :: Q1, Q2           ! Adiabatic cooling and condensational growth terms respectively, -
real(wp)              :: dRHdt            ! Change in ambient saturation ratio/relative humidity over time, 1/s
real(wp)              :: cbrl             ! Cubes in real format

es = SatVapWater(Tmp, 'liq')
L = LatentHeatEvap(Tmp)
rhoA = DensityAir(Tmp, Prs, RH)

cbrl = REAL(Cubes)

! For adiabatic cooling
Q1 = (1_wp/Tmp)*((L*g)/(Rv*Tmp*cp)) - g/(Rd*Tmp)
! For condensation onto droplets
Q2 = Prs/(eta*es) + (1_wp/Tmp)*(L**2/(Rv*Tmp*cp))

dRHdt = Q1*W - (Q2/rhoA)*WetGrowthSums/(cbrl*1e-6_wp) ! per cubic centimeter factor

end function SaturationTime

!==========================================================================================================================

end module EnvEquations
