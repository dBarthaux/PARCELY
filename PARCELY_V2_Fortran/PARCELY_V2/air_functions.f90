module AirFuncs

use EnvironmentConstants, only: Rd, eta
use dvode_kinds_module, only: wp => dvode_wp

contains

!==========================================================================================================================

!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Function to calculate the coefficient of diffusivity of air given a temperature T  *
!*   and pressure P. The equation is from Seinfeld and Pandis Chapter 17, equivalently  *
!*   also from Pruppacher and Klett (1997).                                             *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Dan Barthaux                                                                       *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2025/11/04                                                      *
!*                                                                                      *
!****************************************************************************************

pure elemental function DiffusivityAir(Tmp, Prs) result(Dair)

implicit none

real(wp), intent(in)  :: Tmp      ! Ambient temperature of air, Kelvin
real(wp), intent(in)  :: Prs      ! Ambient pressure, Pa
real(wp)              :: T0       ! Reference temperature, Kelvin
real(wp)              :: P0       ! Reference pressure, Pa
real(wp)              :: Dair     ! Diffusivity of air, m^2/s

T0 = 273.15_wp
P0 = 101325.0_wp

Dair = 1e-4_wp*0.211_wp*((Tmp/T0)**1.94_wp)*(P0/Prs)

end function DiffusivityAir

!==========================================================================================================================

!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Function to calculate the thermal conductivity of air given a temperature T.       *
!*   The equation is from Seinfeld and Pandis Chapter 17, equivalently also from        *
!*    Pruppacher and Klett (1997).                                                      *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Dan Barthaux                                                                       *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2025/11/04                                                      *
!*                                                                                      *
!**************************************************************************************** 
    
pure elemental function ThermConductAir(Tmp) result(Kair)

implicit none

real(wp), intent(in)  :: Tmp      ! Ambient temperature of air, Kelvin
real(wp)              :: T0       ! Reference temperature, Kelvin
real(wp)              :: Kair     ! Thermal conductivity of air, J/m/s/K

T0 = 273.15_wp
Kair = 4.1868e-3_wp*(5.69_wp + 0.017_wp*(Tmp - T0))

end function ThermConductAir

!==========================================================================================================================

!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Function to calculate the density of air given a temperature T, pressure P, and    *
!*   relative humidity RH.                                                              *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Dan Barthaux                                                                       *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2025/11/04                                                      *
!*                                                                                      *
!**************************************************************************************** 

pure elemental function DensityAir(Tmp, Prs, RH) result(Rhoair)

use WaterFuncs, only: SatVapWater

implicit none

real(wp), intent(in)  :: Tmp      ! Ambient temperature of air, Kelvin
real(wp), intent(in)  :: Prs      ! Ambient pressure, Pa
real(wp), intent(in)  :: RH       ! Relative humidity, -
real(wp)              :: es       ! Saturation vapor pressure, Pa
real(wp)              :: e        ! Partial pressure of water vapor, Pa
real(wp)              :: Rhoair   ! Denisty of dry air, kg/m^3

es = SatVapWater(Tmp, 'liq')
e = es*RH
Rhoair = (Prs - e)/(Rd*Tmp)

end function DensityAir  
    
!==========================================================================================================================

!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Function to calculate the virtual air temperature.                                 *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Dan Barthaux                                                                       *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2025/11/04                                                      *
!*                                                                                      *
!**************************************************************************************** 

pure elemental function VirtualTemp(Tmp, Prs, RH) result(Tv)

use WaterFuncs, only: SatVapWater

implicit none

real(wp), intent(in)  :: Tmp      ! Ambient temperature of air, Kelvin
real(wp), intent(in)  :: Prs      ! Ambient pressure, Pa
real(wp), intent(in)  :: RH       ! Relative humidity, -
real(wp)              :: es       ! Saturation vapor pressure, Pa
real(wp)              :: e        ! Partial pressure of water vapor, Pa
real(wp)              :: w        ! Mixing ratio of water vapor, -
real(wp)              :: Tv       ! Virtual temperature, K

es = SatVapWater(Tmp, 'liq')
e = es*RH
w = eta*e/Prs
Tv = Tmp*((1_wp + w/eta)/(1_wp + w))

end function VirtualTemp

!==========================================================================================================================

end module AirFuncs
