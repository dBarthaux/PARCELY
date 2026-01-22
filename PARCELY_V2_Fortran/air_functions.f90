module AirFuncs

use EnvironmentConstants, only: Rd, eta
use dvode_kinds_module, only: dvode_wp

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

function DiffusivityAir(Tmp, Prs) result(Dair)

implicit none

real(dvode_wp), intent(in)  :: Tmp      ! Ambient temperature of air, Kelvin
real(dvode_wp), intent(in)  :: Prs      ! Ambient pressure, Pa
real(dvode_wp)              :: T0       ! Reference temperature, Kelvin
real(dvode_wp)              :: P0       ! Reference pressure, Pa
real(dvode_wp)              :: Dair     ! Diffusivity of air, m^2/s

T0 = 273.15
P0 = 101325.0

Dair = 1e-4*0.211*((Tmp/T0)**1.94)*(P0/Prs)

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
    
function ThermConductAir(Tmp) result(Kair)

implicit none

real(dvode_wp), intent(in)  :: Tmp      ! Ambient temperature of air, Kelvin
real(dvode_wp)              :: T0       ! Reference temperature, Kelvin
real(dvode_wp)              :: Kair     ! Thermal conductivity of air, J/m/s/K

T0 = 273.15
Kair = 4.1868e-3*(5.69 + 0.017*(Tmp - T0))

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

function DensityAir(Tmp, Prs, RH) result(Rhoair)

use WaterFuncs, only: SatVapWater

implicit none

real(dvode_wp), intent(in)  :: Tmp      ! Ambient temperature of air, Kelvin
real(dvode_wp), intent(in)  :: Prs      ! Ambient pressure, Pa
real(dvode_wp), intent(in)  :: RH       ! Relative humidity, -
real(dvode_wp)              :: es       ! Saturation vapor pressure, Pa
real(dvode_wp)              :: e        ! Partial pressure of water vapor, Pa
real(dvode_wp)              :: Rhoair   ! Denisty of dry air, kg/m^3

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

function VirtualTemp(Tmp, Prs, RH) result(Tv)

use WaterFuncs, only: SatVapWater

implicit none

real(dvode_wp), intent(in)  :: Tmp      ! Ambient temperature of air, Kelvin
real(dvode_wp), intent(in)  :: Prs      ! Ambient pressure, Pa
real(dvode_wp), intent(in)  :: RH       ! Relative humidity, -
real(dvode_wp)              :: es       ! Saturation vapor pressure, Pa
real(dvode_wp)              :: e        ! Partial pressure of water vapor, Pa
real(dvode_wp)              :: w        ! Mixing ratio of water vapor, -
real(dvode_wp)              :: Tv       ! Virtual temperature, K

es = SatVapWater(Tmp, 'liq')
e = es*RH
w = eta*e/Prs
Tv = Tmp*((1 + w/eta)/(1 + w))

end function VirtualTemp

!==========================================================================================================================

end module AirFuncs