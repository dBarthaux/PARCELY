module DropFuncs

use dvode_kinds_module, only: dvode_wp
use EnvironmentConstants, only: Mw, R, rhoL, Rv, pi, Na, Rd, cp
use ModelParameters, only: norg, ndrop, ninorg, OrganicProperties, InorganicBase, Msft, sftc, DropSurfTens, ac, at

contains

!==========================================================================================================================

!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Function to calculate the saturation ratio of water vapor at the surface of a      *
!*   droplet using Kappa-Koehler theory by Petters and Kreidenweis (2007).              *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Dan Barthaux                                                                       *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2025/11/04                                                      *
!*                                                                                      *
!****************************************************************************************

function KappaKoehler(T, SolRad, SolKappa, WetRad, DropSrfc) result(SatDrop)

use WaterFuncs, only: WaterSurfaceTension

implicit none

real(dvode_wp), intent(in)  :: T                ! Ambient temperature of air, Kelvin
real(dvode_wp), intent(in)  :: SolRad           ! Droplet dry radius (solute radius), m
real(dvode_wp), intent(in)  :: SolKappa         ! Droplet kappa (solute hygroscopicity), -
real(dvode_wp), intent(in)  :: WetRad           ! Droplet wet radius (total radius), m
real(dvode_wp), intent(in)  :: DropSrfc         ! Droplet surface tension, J/m^2
real(dvode_wp)              :: SatDrop          ! Droplet water vapor saturation term of the Kappa-Koehler equation
real(dvode_wp)              :: ActivityTerm     ! Activity term of the Kappa-Koehler equation
real(dvode_wp)              :: SurfTensTerm     ! Surface tension term of the Kappa-Koehler equation

ActivityTerm = (WetRad**3 - SolRad**3)/(WetRad**3 - (SolRad**3)*(1 - SolKappa))
SurfTensTerm = EXP((2*DropSrfc*Mw)/(T*R*WetRad*rhoL))

SatDrop = ActivityTerm*SurfTensTerm

end function KappaKoehler

!==========================================================================================================================

!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Function to calculate the condensational (water vapor) growth rate of a droplet.   *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Dan Barthaux                                                                       *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2025/11/04                                                      *
!*                                                                                      *
!****************************************************************************************

function DropGrowthRate(Tmp, Prs, RH, WetRad, SolRad, SolKappa, DropSrfc) result(GrowthRate)

use WaterFuncs, only: SatVapWater, LatentHeatEvap

implicit none

real(dvode_wp), intent(in)  :: Tmp          ! Ambient temperature of air, Kelvin
real(dvode_wp), intent(in)  :: Prs          ! Ambient pressure, Pa
real(dvode_wp), intent(in)  :: RH           ! Relative humidity, -
real(dvode_wp), intent(in)  :: WetRad       ! Droplet wet radius (total radius), m
real(dvode_wp), intent(in)  :: SolRad       ! Droplet dry radius (solute radius), m
real(dvode_wp), intent(in)  :: SolKappa     ! Droplet kappa (solute hygroscopicity), -
real(dvode_wp), intent(in)  :: DropSrfc     ! Droplet surface tension, J/m^2
real(dvode_wp)              :: SatDrop      ! Droplet water vapor saturation term of the Kappa-Koehler equation
real(dvode_wp)              :: Kn           ! Adjusted thermal conductivity of air, J/m/K/s
real(dvode_wp)              :: Dn           ! Adjusted diffusivity of air, m^2/s
real(dvode_wp)              :: Fk           ! Heat conduction term, m s/kg
real(dvode_wp)              :: Fd           ! Molecular diffusion term, m s/kg
real(dvode_wp)              :: GrowthRate   ! Droplet growth rate, molec/s
real(dvode_wp)              :: es           ! Saturation vapor pressure, Pa
real(dvode_wp)              :: L            ! Latent heat of vaporization, J/kg

es = SatVapWater(Tmp, 'liq')
SatDrop = KappaKoehler(Tmp, SolRad, SolKappa, WetRad, DropSrfc)

call KineticEffects(WetRad, Tmp, Prs, RH, Kn, Dn)
L = LatentHeatEvap(Tmp)
Fk = (L/(Rv*Tmp) - 1.0)*L*rhoL/(Kn*Tmp)
Fd = Rv*Tmp*rhoL/(Dn*es)

GrowthRate = (RH - SatDrop)/(WetRad*(Fk + Fd))

end function DropGrowthRate

!==========================================================================================================================

!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Wrapper function to help calculate the equilibrium radius of a droplet using       *
!*   Kappa-Koehler theory by Petters and Kreidenweis (2007).                            *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Dan Barthaux                                                                       *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2025/11/04                                                      *
!*                                                                                      *
!****************************************************************************************

function InvertedKoehler(WetRad, Tmp, RH, SolRad, SolKappa) result(SatDrop)

use WaterFuncs, only: WaterSurfaceTension

implicit none

real(dvode_wp), intent(in)  :: Tmp              ! Ambient temperature of air, Kelvin
real(dvode_wp), intent(in)  :: RH               ! Ambient relative humidity, -
real(dvode_wp), intent(in)  :: SolRad           ! Droplet dry radius (solute radius), m
real(dvode_wp), intent(in)  :: SolKappa         ! Droplet kappa (solute hygroscopicity), -
real(dvode_wp), intent(in)  :: WetRad           ! Droplet wet radius (total radius), m
real(dvode_wp)              :: SatDrop          ! Droplet water vapor saturation term of the Kappa-Koehler equation
real(dvode_wp)              :: SigWater         ! Surface tension of pure water, J/m^2
real(dvode_wp)              :: ActivityTerm     ! Activity term of the Kappa-Koehler equation
real(dvode_wp)              :: SurfTensTerm     ! Surface tension term of the Kappa-Koehler equation

SigWater = WaterSurfaceTension(Tmp)
ActivityTerm = (WetRad**3 - SolRad**3)/(WetRad**3 - (SolRad**3)*(1 - SolKappa))
SurfTensTerm = EXP((2*SigWater*Mw)/(Tmp*R*WetRad*rhoL))

SatDrop = ActivityTerm*SurfTensTerm - RH

end function InvertedKoehler

!==========================================================================================================================

!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Wrapper function to help calculate the critical radius of a droplet using          *
!*   Kappa-Koehler theory by Petters and Kreidenweis (2007).                            *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Dan Barthaux                                                                       *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2025/11/04                                                      *
!*                                                                                      *
!****************************************************************************************

function CriticalKoehler(WetRad, Tmp, SolRad, SolKappa) result(SatDrop)

use WaterFuncs, only: WaterSurfaceTension

implicit none

real(dvode_wp), intent(in)  :: Tmp              ! Ambient temperature of air, Kelvin
real(dvode_wp), intent(in)  :: SolRad           ! Droplet dry radius (solute radius), m
real(dvode_wp), intent(in)  :: SolKappa         ! Droplet kappa (solute hygroscopicity), -
real(dvode_wp), intent(in)  :: WetRad           ! Droplet wet radius (total radius), m
real(dvode_wp)              :: SatDrop          ! Droplet water vapor saturation term of the Kappa-Koehler equation
real(dvode_wp)              :: SigWater         ! Surface tension of pure water, J/m^2
real(dvode_wp)              :: ActivityTerm     ! Activity term of the Kappa-Koehler equation
real(dvode_wp)              :: SurfTensTerm     ! Surface tension term of the Kappa-Koehler equation

SigWater = WaterSurfaceTension(Tmp)
ActivityTerm = (WetRad**3 - SolRad**3)/(WetRad**3 - (SolRad**3)*(1 - SolKappa))
SurfTensTerm = EXP((2*SigWater*Mw)/(Tmp*R*WetRad*rhoL))
! Return a negative value
SatDrop = -(ActivityTerm*SurfTensTerm - 1)*100

end function CriticalKoehler

!==========================================================================================================================

!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Function to calculate the equilibrium radius of a water droplet using Kappa-       *
!*   Koehler theory.                                                                    *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Dan Barthaux                                                                       *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2025/11/04                                                      *
!*                                                                                      *
!****************************************************************************************

subroutine EquilibriumRadius(xequil, Tmp, RH, SolRad, SolKappa)

use root_module
use fmin_module

implicit none

real(dvode_wp), intent(in)      :: Tmp, RH, SolRad, SolKappa
real(dvode_wp), intent(out)     :: xequil
integer                         :: iflag
real(dvode_wp)                  :: xcrit, f

! Get the radius at which the droplet activates
xcrit = fmin(g_wrapper, SolRad*1.0001, 1e-4_dvode_wp, 1e-9_dvode_wp)
! Find the equilibrium radius using the radius of the solute (slightly larger) and the activation radius as bounds
call root_scalar('toms748', f_wrapper, SolRad*1.0001, xcrit, xequil, f, iflag)

contains

! Wrapper functions for the minimization/root-finding
function f_wrapper(r) result(val)
    real(dvode_wp), intent(in) :: r
    real(dvode_wp) :: val
    val = InvertedKoehler(r, Tmp, RH, SolRad, SolKappa)
end function f_wrapper

function g_wrapper(r) result(val)
    real(dvode_wp), intent(in) :: r
    real(dvode_wp) :: val
    val = CriticalKoehler(r, Tmp, SolRad, SolKappa)
end function g_wrapper

end subroutine EquilibriumRadius

!==========================================================================================================================

!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Function to calculate thermal conductivity and diffusivity adjusted for kinetic    *
!*   effects using equations from Seinfeld and Pandis (2006)                            *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Dan Barthaux                                                                       *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2025/11/04                                                      *
!*                                                                                      *
!**************************************************************************************** 

subroutine KineticEffects(Radius, T, P, RH, Kn, Dn)

use AirFuncs, only: DensityAir, ThermConductAir, DiffusivityAir

implicit none

real(dvode_wp), intent(in)      :: Radius   ! Droplet radius, m
real(dvode_wp), intent(in)      :: T        ! Ambient temperature of air, Kelvin
real(dvode_wp), intent(in)      :: P        ! Ambient pressure, Pa
real(dvode_wp), intent(in)      :: RH       ! Relative humidity, -
real(dvode_wp)                  :: Kr       ! Thermal conductivity of air, J/m/K/s
real(dvode_wp)                  :: Dr       ! Diffusivity of air, m^2/s
real(dvode_wp)                  :: Rhoair   ! Density of dry air, kg/m^3
real(dvode_wp)                  :: Vv       ! Vapor jump length, m
real(dvode_wp)                  :: VT       ! Thermal jump length, m
real(dvode_wp), intent(out)     :: Kn       ! Adjusted thermal conductivity of air, J/m/K/s
real(dvode_wp), intent(out)     :: Dn       ! Adjusted diffusivity of air, m^2/s

Kr = ThermConductAir(T)
Dr = DiffusivityAir(T, P)
Rhoair = DensityAir(T, P, RH)
    
Vv = 1.3*(8e-8)
VT = 2.16e-7
        
Kn = Kr/(Radius/(Radius + VT) + (Kr/(at*Radius*Rhoair*cp))*SQRT(2.0*pi/(Rd*T)))
Dn = Dr/(Radius/(Radius + Vv) + (Dr/(ac*Radius))*SQRT(2.0*pi/(Rv*T)))

end subroutine KineticEffects

!==========================================================================================================================

!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Function to calculate the droplet surface tension.                                 *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Dan Barthaux                                                                       *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2025/11/04                                                      *
!*                                                                                      *
!****************************************************************************************

subroutine DropletSurfaceTension(Tmp, WaterMass, OrganicMass, WaterRadius)

use WaterFuncs, only: WaterSurfaceTension

implicit none

real(dvode_wp), intent(in)                          :: Tmp ! Ambient temperature of air, Kelvin
real(dvode_wp), dimension(ndrop), intent(in)        :: WaterMass, WaterRadius
real(dvode_wp), dimension(ndrop, norg), intent(in)  :: OrganicMass
integer                                             :: i
real(dvode_wp)                                      :: delz ! Thickness of the monolayer, m 
real(dvode_wp), dimension(ndrop, norg+1)            :: DropMoles, MoleFracs
real(dvode_wp), dimension(norg+1)                   :: AllSfts
real(dvode_wp), dimension(ndrop, norg)              :: Vorganics
real(dvode_wp), dimension(ndrop)                    :: InorganicMoles, TotalMoles, Vmono, Vwater, Vinorg, OrgTension, SurfaceCover, InorgMoles

! Constant value set in the environment file
if (Msft == 1) then
    DropSurfTens = sftc

! Surface tension of pure water, temperature dependent
else if (Msft == 2) then
    DropSurfTens = WaterSurfaceTension(Tmp)

! Direct mole fraction-weighted mixing ratio approach
else if (Msft == 3) then
    ! Get moles of each component in each droplet
    ! Moles of water
    DropMoles(:,norg+1) = WaterMass/Mw
    ! Moles of inorganic
    InorgMoles = InorganicBase(:,1)/InorganicBase(:,3)
    ! Moles of organics
    do i = 1, norg
        DropMoles(:,i) = OrganicMass(:,i)/OrganicProperties(i,3)
    end do
    ! Sum moles
    TotalMoles = SUM(DropMoles, 2) ! + InorgMoles
    ! Calculate mole fractions
    do i = 1, norg+1
        MoleFracs(:,i) = DropMoles(:,i)/TotalMoles
    end do
    ! Get surface tensions of components
    AllSfts(1:norg) = OrganicProperties(:,6)
    AllSfts(norg+1) = WaterSurfaceTension(Tmp)
    ! Assume inorganic does not contribute to surface tension
    DropSurfTens = MATMUL(MoleFracs, AllSfts)

! Organic-film method
else if (Msft == 4) then
    ! Characteristic monolayer thickness taken from AIOMFAC
    delz = 0.2e-9
    ! Volume of the monolayer, m^3
    Vmono = (4*pi/3)*(WaterRadius**3 - (WaterRadius - delz)**3)
    ! Volumes for organics, water, and inorganics from mass and density
    Vwater  = WaterMass/rhoL
    Vinorg  = InorganicBase(:,1)/InorganicBase(:,4)
    do i = 1, norg
        Vorganics(:,i) = OrganicMass(:,i)/OrganicProperties(i,4)
    end do
    ! Determine if there's enough organics to cover the particle
    do i = 1, ndrop
        ! Get volume weighted-mean surface tension of the organics
        OrgTension(i) = SUM(OrganicProperties(:,6)*Vorganics(i,:)/SUM(Vorganics(i,:)))

        if (SUM(Vorganics(i,:)) <= Vmono(i)) then
            ! If not enough, define the surface coverage fraction
            SurfaceCover(i) = SUM(Vorganics(i,:))/Vmono(i)
        else
            ! If enough, set as 1
            SurfaceCover(i) = 1
        end if
    end do
    ! Final droplet surface tension
    DropSurfTens =  OrgTension*SurfaceCover + (1-SurfaceCover)*WaterSurfaceTension(Tmp)
    
end if

end subroutine DropletSurfaceTension

!==========================================================================================================================

end module DropFuncs