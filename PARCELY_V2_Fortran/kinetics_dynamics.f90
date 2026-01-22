module KinDynFuncs

use EnvironmentConstants, only: rhoL, g, cp, Rd, Rv, pi
use dvode_kinds_module, only: dvode_wp

contains

!==========================================================================================================================
    
!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Function to calculate the drag coefficient for a droplet using a fitted power-law  *
!*   function to data from Shafrir and Gal-Chen (1971).                                 *
!*   https://doi.org/10.1175/1520-0469(1971)028%3C0741:ANSOCE%3E2.0.CO;2                *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Dan Barthaux                                                                       *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2025/11/04                                                      *
!*                                                                                      *
!**************************************************************************************** 
    
function DragCoefficient(Radius) result(Cd)

implicit none

real(dvode_wp), intent(in)  :: Radius   ! Droplet radius, m
real(dvode_wp)              :: Cd       ! Drag coefficient, unitless

! Coefficient with 95% confidence bounds
! a = 9.86e-10    #(7.956e-10, 1.176e-09)
! b = -2.375      #(-2.394, -2.357)
! c = 1.014       #(0.9607, 1.068)

Cd = 9.86e-10*(Radius**(-2.375)) + 1.014

end function DragCoefficient
    
!==========================================================================================================================
    
!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Function to determine the terminal velocity of a droplet as a function of its      *
!*   radius, using the analytic expression from force balance.                          *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Dan Barthaux                                                                       *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2025/11/04                                                      *
!*                                                                                      *
!**************************************************************************************** 

function TerminalVelocity(Radius, T, P, RH) result(Vt)

use AirFuncs, only: DensityAir

implicit none

real(dvode_wp), intent(in)  :: Radius   ! Droplet radius, m
real(dvode_wp), intent(in)  :: T        ! Ambient temperature of air, Kelvin
real(dvode_wp), intent(in)  :: P        ! Ambient pressure, Pa
real(dvode_wp), intent(in)  :: RH       ! Relative humidity, -
real(dvode_wp)              :: Cd       ! Drag coefficient, -
real(dvode_wp)              :: Vt       ! Terminal velocity of a droplet, m/s
real(dvode_wp)              :: Rhoair   ! Density of dry air, kg/m^3

Rhoair = DensityAir(T, P, RH)
Cd = DragCoefficient(Radius)
Vt = SQRT(8*Radius*g*(rhoL - Rhoair)/(3*Cd*Rhoair))

end function TerminalVelocity

!==========================================================================================================================

end module KinDynFuncs