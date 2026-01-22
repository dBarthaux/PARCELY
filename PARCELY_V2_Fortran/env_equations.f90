module EnvEquations

use dvode_kinds_module, only: dvode_wp
use EnvironmentConstants, only: g, Rd, cp, Rv, eta
use ModelParameters, only: W

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

function PressureTime(Tmp, Prs, RH) result (dPdt)

use AirFuncs, only: VirtualTemp

implicit none

real(dvode_wp), intent(in)  :: Tmp      ! Temperature, K
real(dvode_wp), intent(in)  :: Prs      ! Pressure, Pa
real(dvode_wp), intent(in)  :: RH       ! Relative humidity, -
real(dvode_wp)              :: Tv       ! Virtual temperature, K
real(dvode_wp)              :: dPdt     ! Change in pressure over time, Pa/s

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

function TemperatureTime(Tmp, Prs, RH, WetGrowthSums) result(dTdt)

use WaterFuncs, only: LatentHeatEvap
use AirFuncs, only: DensityAir

implicit none

real(dvode_wp), intent(in)  :: Tmp              ! Temperature, K
real(dvode_wp), intent(in)  :: Prs              ! Pressure, Pa
real(dvode_wp), intent(in)  :: RH               ! Relative humidity, -
real(dvode_wp), intent(in)  :: WetGrowthSums    ! Sum of condensational growth rates (water vapor) of the droplets, dm/dt [kg/s]
real(dvode_wp)              :: L                ! Latent heat of evaporation (of water vapor), J/kg
real(dvode_wp)              :: rhoA             ! Density of dry air, kg/m3
real(dvode_wp)              :: dTdt             ! Change in temperature over time, K/s

L = LatentHeatEvap(Tmp)
rhoA = DensityAir(Tmp, Prs, RH)
dTdt = -g*W/cp + (L/(cp*rhoA))*WetGrowthSums/((1e-2)**3) ! per cubic centimeter factor

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

function SaturationTime(Tmp, Prs, RH, WetGrowthSums) result(dRHdt)

use WaterFuncs, only: LatentHeatEvap, SatVapWater
use AirFuncs, only: DensityAir

implicit none

real(dvode_wp), intent(in)  :: Tmp              ! Temperature, K
real(dvode_wp), intent(in)  :: Prs              ! Pressure, Pa
real(dvode_wp), intent(in)  :: RH               ! Relative humidity, -
real(dvode_wp), intent(in)  :: WetGrowthSums    ! Sum of condensational growth rates (water vapor) of the droplets, dm/dt [kg/s]
real(dvode_wp)              :: L                ! Latent heat of evaporation (of water vapor), J/kg
real(dvode_wp)              :: es               ! Saturation vapor pressure of water, Pa
real(dvode_wp)              :: rhoA             ! Density of dry air, kg/m3
real(dvode_wp)              :: Q1, Q2           ! Adiabatic cooling and condensational growth terms respectively, -
real(dvode_wp)              :: dRHdt            ! Change in ambient saturation ratio/relative humidity over time, 1/s

es = SatVapWater(Tmp, 'liq')
L = LatentHeatEvap(Tmp)
rhoA = DensityAir(Tmp, Prs, RH)

! For adiabatic cooling
Q1 = (1/Tmp)*((L*g)/(Rv*Tmp*cp)) - g/(Rd*Tmp)
! For condensation onto droplets
Q2 = Prs/(eta*es) + (1/Tmp)*(L**2/(Rv*Tmp*cp))

dRHdt = Q1*W - (Q2/rhoA)*WetGrowthSums/((1e-2)**3) ! per cubic centimeter factor

end function SaturationTime

!==========================================================================================================================

end module EnvEquations