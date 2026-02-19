module EnvironmentConstants

use dvode_kinds_module, only: wp => dvode_wp

!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Collection of constants used in PARCELY.                                           *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Dan Barthaux                                                                       *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2025/11/04                                                      *
!*                                                                                      *
!**************************************************************************************** 

implicit none

real(wp), parameter :: g         = 9.81                   ! Gravitational acceleration, m/s^2
real(wp), parameter :: rhoL      = 997.0                  ! Density of liquid water, kg/m^3
real(wp), parameter :: Ma        = 28.97e-3               ! Molar mass of dry air, kg/mol
real(wp), parameter :: Mw        = 18.01528e-3            ! Molar Mass of water, kg/mol
real(wp), parameter :: R         = 8.31446261815          ! Universal molar gas constant, J/K/mol
real(wp), parameter :: cp        = 1005.0                 ! Specific heat capacity (constant pressure), J/kg/K
real(wp), parameter :: Rv        = 461.0                  ! Gas constant for water vapour, J/kg/K
real(wp), parameter :: Rd        = 287.0                  ! Gas constant of dry air, J/kg/K
real(wp), parameter :: eta       = Rd/Rv
real(wp), parameter :: kB        = 1.380649e-23           ! Boltzmann constant, J/K
real(wp), parameter :: Na        = 6.02214076e23          ! Avogadro's number, molecules/mol
real(wp), parameter :: R_alt     = 8.20573e-5             ! Ideal gas constant, m3 atm/K/mol
real(wp), parameter :: cpl       = 4220                   ! Isobaric heat capacity of water vapor, J/K
real(wp), parameter :: cpi       = 2097                   ! Isobaric heat capacity of ice, J/K
real(wp), parameter :: pi        = 4*atan(1.0_wp)         ! Pi

end module EnvironmentConstants