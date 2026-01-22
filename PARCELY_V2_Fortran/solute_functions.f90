module SolFuncs

use dvode_kinds_module, only: dvode_wp
use EnvironmentConstants, only: pi, R, Na, Rd_alt, Mw, rhoL
use ModelParameters, only: norg, ninorg, ndrop, OrganicProperties, InorganicBase, SoluteProperties, DropSurfTens

contains

!==========================================================================================================================

!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Function to create a the size distribution for the cloud droplet nuclei/solutes.   *
!*   Options include mono, normal, and log-normal distributions.                        *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Dan Barthaux                                                                       *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2025/11/04                                                      *
!*                                                                                      *
!****************************************************************************************

function SoluteRadiusDistribution(npops, DistConcs, DistRads, DistStds, DistType) result(DryRad)

use DistFuncs, only: NormDist, LogNormDist

implicit none

integer, intent(in)                             :: npops        ! Number of different distributions
integer, dimension(npops), intent(in)           :: DistConcs    ! Number concentration of each distribution
real(dvode_wp), dimension(npops), intent(in)    :: DistRads     ! Distribution means in um
real(dvode_wp), dimension(npops), intent(in)    :: DistStds     ! Distribution standard deviations in um
integer, intent(in)                             :: DistType     ! Type of distribution used [1: monodispersed, 2:normal, 3:log-normal]
real(dvode_wp), allocatable                     :: DryRad(:)    ! Dry/solute/aerosol radii, m
real(dvode_wp), allocatable                     :: randnum(:)
integer                                         :: i, j

allocate(DryRad(SUM(DistConcs)))
DryRad = 1.0

do i = 1, npops
    
    allocate(randnum(DistConcs(i)))

    ! Mono-disperse distribution
    if (DistType == 1) then
        if (i == 1) then
            DryRad(1:DistConcs(i)) = DryRad(1:DistConcs(i))*DistRads(i)
        else
            j = SUM(DistConcs(1:i-1))
            DryRad(j+1:j+DistConcs(i)) = DryRad(j+1:j+DistConcs(i))*DistRads(i)
        endif

    ! Normal distribution
    elseif (DistType == 2) then
        if (i == 1) then
            call NormDist(DistRads(i), DistStds(i), DistConcs(i), randnum)
            DryRad(1:DistConcs(i)) = randnum
        else
            j = SUM(DistConcs(1:i-1))
            call NormDist(DistRads(i), DistStds(i), DistConcs(i), randnum)
            DryRad(j+1:j+DistConcs(i)) = randnum
        endif

    ! Log-normal distribution
    else
        if (i == 1) then
            call LogNormDist(DistRads(i), DistStds(i), DistConcs(i), randnum)
            DryRad(1:DistConcs(i)) = randnum
        else
            j = SUM(DistConcs(1:i-1))
            call LogNormDist(DistRads(i), DistStds(i), DistConcs(i), randnum)
            DryRad(j+1:j+DistConcs(i)) = randnum
        endif

    endif
    
    deallocate(randnum)

end do

! Convert to meters
DryRad = DryRad*1e-6

end function SoluteRadiusDistribution

!==========================================================================================================================

!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Function to determine the initial solute properties depending on the mixing state  *
!*   of the aerosol population as described in Riemer (2013).                           *
!*   doi:10.5194/acp-13-11423-2013                                                      *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Dan Barthaux                                                                       *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2025/11/04                                                      *
!*                                                                                      *
!****************************************************************************************

subroutine InitialSoluteProperties(npops, Xi, DistConcs, DryRad, InorganicProperties)

implicit none

integer, intent(in)                                     :: npops                ! Number of aerosol populations (different size distributions)
real(dvode_wp), intent(in)                              :: Xi                   ! Mixing state
integer, dimension(npops), intent(in)                   :: DistConcs            ! Number concentration of each distribution
real(dvode_wp), dimension(ndrop), intent(in)            :: DryRad               ! Dry/solute/aerosol radii, m
real(dvode_wp), dimension(ninorg, npops+3), intent(in)  :: InorganicProperties

integer                                                 :: i, n, j, NumPop
real(dvode_wp)                                          :: TotalMoles
real(dvode_wp), allocatable                             :: SolKappa(:), SolDensity(:), SolMolarMass(:), InorganicCounts(:), SolMass(:), SolMoles(:), TotalVolume(:)
real(dvode_wp), allocatable                             :: InorganicHomes(:), InorganicAways(:), InorganicMoles(:), InorganicPopFracs(:)
real(dvode_wp), allocatable                             :: InorganicVolumes(:,:), InorganicVolFracs(:,:), InorganicMasses(:,:), InorganicMoleFracs(:), SolMoleFracs(:,:)
real(dvode_wp), dimension(ninorg)                       :: InorganicKappas, InorganicDensities, InorganicMolarMasses
real(dvode_wp), dimension(ninorg, npops)                :: InorganicFractions

allocate(SolKappa(ndrop), SolDensity(ndrop), SolMolarMass(ndrop), InorganicCounts(ndrop), SolMass(ndrop), SolMoles(ndrop), TotalVolume(ndrop), &
         InorganicHomes(ninorg), InorganicAways(ninorg), InorganicMoles(ninorg), InorganicPopFracs(ninorg), InorganicMoleFracs(ninorg), & 
         InorganicVolumes(ndrop,ninorg), InorganicVolFracs(ndrop,ninorg), InorganicMasses(ndrop,ninorg), SolMoleFracs(ndrop,ninorg))

InorganicMolarMasses    = InorganicProperties(:,1)
InorganicDensities      = InorganicProperties(:,2)
InorganicKappas         = InorganicProperties(:,3)
InorganicFractions      = InorganicProperties(:,4:3+npops)

! Inorganic Section (these are the base solute properties and will always be present, i.e. constant)
n = 1
! For every size distribution
do i = 1, SIZE(DistConcs)
    ! For every inorganic
    do j = 1, ninorg
        ! Get the number of particles of inorganic "j" in distribution "i" using its fraction of total particles in distribution "i"
        NumPop = NINT(InorganicFractions(j,i)*DistConcs(i))
        ! Initialize the base properties with those assigned to inorganic "j"
        SolKappa(n:n+NumPop-1) = InorganicKappas(j)
        SolDensity(n:n+NumPop-1) = InorganicDensities(j)
        SolMolarMass(n:n+NumPop-1) = InorganicMolarMasses(j)
        ! Update indices of pure inorganic "j" particles
        InorganicCounts(n:n+NumPop-1) = j
        ! Update counter
        n = n + NumPop
    end do
end do

! If fully externally mixed, every solute is pure, mass calculation is straightforward and no edits to Kappa/density/molar mass
if (Xi == 0) then
    SolMass = (4*pi/3)*SolDensity*(DryRad**3)

! Otherwise, solute average properties are mole fraction-weighted properties
else
    ! Start with initial mass in moles
    SolMoles = (4*pi/3)*SolDensity*(DryRad**3)/(SolMolarMass)
    ! Get the total number of moles
    TotalMoles = SUM(SolMoles)
    ! Get the moles for each inorganic
    do i = 1, ninorg
        InorganicMoles(i) = SUM(SolMoles, MASK = (InorganicCounts == i))
    end do
    ! Get the total mole fractions for each inorganic
    InorganicMoleFracs = InorganicMoles/TotalMoles
    ! Define the Home and Away terms for each inorganic
    InorganicHomes = 1 + Xi*(InorganicMoleFracs - 1)
    InorganicAways = Xi*InorganicMoleFracs
    ! Assign the mole fraction of each inorganic in every droplet
    do n = 1, ndrop
        do i = 1, ninorg
            if (InorganicCounts(n) == i) then
                SolMoleFracs(n,i) = InorganicHomes(i)
            else
                SolMoleFracs(n,i) = InorganicAways(i)
            end if
        end do
    end do
    ! Solute density and molar mass are mole-fraction weighted averages
    SolDensity = MATMUL(SolMoleFracs, InorganicDensities)
    SolMolarMass = MATMUL(SolMoleFracs, InorganicMolarMasses)
    ! Get the masses for each inorganic
    do n = 1, ndrop
        do i = 1, ninorg
            InorganicMasses(n,i) = (4*pi/3)*(SolMoleFracs(n,i)*InorganicDensities(i))*(DryRad(n)**3)
            ! Get volumes of each inorganic
            InorganicVolumes(n,i) = InorganicMasses(n,i)/InorganicDensities(i)
        end do
    end do
    ! Get the total volume
    TotalVolume = SUM(InorganicVolumes, DIM=2)
    ! Get volume fractions of each inorganic
    do i = 1, ninorg
        InorganicVolFracs(:,i) = InorganicVolumes(:,i)/TotalVolume
    end do
    ! Solute kappa is a volume-weighted mixing rule based on Petters and Kreidenweiss (2007)
    SolKappa = MATMUL(InorganicVolFracs, InorganicKappas)
    SolMass = (4*pi/3)*SolDensity*(DryRad**3)

end if

! SoluteProperties: mass, radius, molar mass, density, kappa
allocate(SoluteProperties(ndrop, 5))
SoluteProperties(:,1) = SolMass
SoluteProperties(:,2) = DryRad
SoluteProperties(:,3) = SolMolarMass
SoluteProperties(:,4) = SolDensity
SoluteProperties(:,5) = SolKappa

! Create inorganics file
open(unit=12, file='inorganics_output.txt', status='replace', action='write')
write(12, '(A)') 'Radius (m)   Density (kg/m3)   Molar Mass (g/mol)   Kappa (-)'
do i = 1,ndrop
    write(12, '(4ES15.6)') DryRad(i), SolDensity(i), SolMolarMass(i), SolKappa(i)
end do
close(12)

end subroutine InitialSoluteProperties

!==========================================================================================================================

!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Function to initialize the solute/inorganic population of aerosols/particles.      *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Dan Barthaux                                                                       *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2025/11/04                                                      *
!*                                                                                      *
!****************************************************************************************

subroutine SolPopInitialize(npops, DistType, DistConcs, DistRads, DistStds, Xi, InorganicProperties)

implicit none

integer, intent(in)                                     :: npops, DistType
real(dvode_wp), intent(in)                              :: Xi                   ! Mixing State
integer, dimension(npops), intent(in)                   :: DistConcs            ! Number concentration of each distribution
real(dvode_wp), dimension(npops), intent(in)            :: DistRads, DistStds
real(dvode_wp), dimension(ninorg, 3+npops), intent(in)  :: InorganicProperties
real(dvode_wp), allocatable                             :: DryRad(:)


! Get dry radii
DryRad = SoluteRadiusDistribution(npops, DistConcs, DistRads, DistStds, DistType)
!Set any radius smaller than 5 nm to 5 nm
WHERE (DryRad < 5e-9)
    DryRad = 5e-9
END WHERE

ndrop = SIZE(DryRad)

! Get solute properties
call InitialSoluteProperties(npops, Xi, DistConcs, DryRad, InorganicProperties)

allocate(InorganicBase(ndrop, 5))
InorganicBase = SoluteProperties

end subroutine SolPopInitialize

!==========================================================================================================================

!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Function to update the solute properties of the droplet population.                *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Dan Barthaux                                                                       *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2025/11/04                                                      *
!*                                                                                      *
!****************************************************************************************

subroutine SoluteUpdate(OrganicMass)

implicit none

real(dvode_wp), dimension(ndrop, norg), intent(in) :: OrganicMass

integer                                 :: i, j
real(dvode_wp), dimension(ndrop)        :: InorganicMoles, InorganicVolumes, InorganicMoleFrac, InorganicVolFrac, TotalMoles, TotalVolume
real(dvode_wp), dimension(ndrop, norg)  :: OrganicMoles, OrganicMoleFrac, OrganicVolumes, OrganicVolFrac
!real(dvode_wp), allocatable :: OrgMassKG(:,:)

!allocate(OrgMassKG(ndrop, norg))

!do i = 1, ndrop
!    OrgMassKG(i,:) = OrganicProperties(:,3)*OrganicMass(i,:)
!end do

! Update solute mass
SoluteProperties(:,1) = InorganicBase(:,1) + SUM(OrganicMass, 2)

! Calculate moles
InorganicMoles = InorganicBase(:,1)/InorganicBase(:,3)
do i = 1, ndrop
    OrganicMoles(i,:) = OrganicMass(i,:)/OrganicProperties(:,3)
end do
! Get the total moles
TotalMoles = InorganicMoles + SUM(OrganicMoles, 2)

! Mole fractions
InorganicMoleFrac = InorganicMoles/TotalMoles
do i = 1, ndrop
    OrganicMoleFrac(i,:) = OrganicMoles(i,:)/TotalMoles(i)
end do
! Solute density and molar mass are mole-fraction weighted averages
SoluteProperties(:,3) = InorganicMoleFrac*InorganicBase(:,3) + MATMUL(OrganicMoleFrac, OrganicProperties(:,3))
SoluteProperties(:,4) = InorganicMoleFrac*InorganicBase(:,4) + MATMUL(OrganicMoleFrac, OrganicProperties(:,4))

! Update solute radius
SoluteProperties(:,2) = ((3/(4*pi))*SoluteProperties(:,1)/SoluteProperties(:,4))**(1.0/3.0)

! Get the volumes for each component
InorganicVolumes = InorganicBase(:,1)/InorganicBase(:,4)
do i = 1, ndrop
    OrganicVolumes(i,:) = OrganicMass(i,:)/OrganicProperties(:,4)
end do
! Get the total volume
TotalVolume = InorganicVolumes + SUM(OrganicVolumes, 2)
! Get volume fractions of each component
InorganicVolFrac = InorganicVolumes/TotalVolume
do i = 1, ndrop
    OrganicVolFrac(i,:) = OrganicVolumes(i,:)/TotalVolume(i)
end do
! Solute kappa is a volume-weighted mixing rule based on Petters and Kreidenweiss (2007)
SoluteProperties(:,5) = InorganicVolFrac*InorganicBase(:,5) + MATMUL(OrganicVolFrac, OrganicProperties(:,5))

end subroutine SoluteUpdate

!==========================================================================================================================

!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Function to calculate the condensation rate of organics onto the droplets.         *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Dan Barthaux                                                                       *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2025/11/04                                                      *
!*                                                                                      *
!****************************************************************************************

subroutine Cocondense(OrganicMass, Tmp, WaterMass, WaterRad, dOrgdt)

implicit none

real(dvode_wp), dimension(ndrop, norg), intent(in)  :: OrganicMass
real(dvode_wp), dimension(ndrop), intent(in)        :: WaterMass, WaterRad
real(dvode_wp), intent(in)                          :: Tmp                  ! Ambient temperature, K
real(dvode_wp), dimension(ndrop, norg), intent(out) :: dOrgdt
integer                                             :: i, j
real(dvode_wp), allocatable                         :: MolVol(:), KelFactor(:), EquiPress(:,:), Csurf(:,:), DropDensity(:), OrganicGas(:), &
                                                       MoleFracs(:,:), C0tadj(:), DiffOrgs(:), Psat(:), Denom(:), Dact(:)

allocate(MolVol(ndrop), KelFactor(ndrop), EquiPress(ndrop, norg), Csurf(ndrop, norg), DropDensity(ndrop), &
        OrganicGas(norg), MoleFracs(ndrop, norg+2), C0tadj(norg), DiffOrgs(norg), Psat(norg), Denom(ndrop), Dact(ndrop))

! Get gaseous concentrations by converting total to kg/m3 and subtracting condensed mass
OrganicGas = OrganicProperties(:,1)*1e-9 - SUM(OrganicMass, 1)*1e6 !kg in a cm3 to kg/m3
do i = 1, norg
    if (OrganicGas(i) < 0) then
        OrganicGas(i) = 0.0
    end if
end do

! Saturation concentration, ug/m3
C0tadj = SatConcentration(OrganicProperties(:,2), Tmp)
! Convert C0 to a pure component saturation vapour pressure [Pa]
Psat = 101325*(C0tadj*Rd_alt*Tmp)/(1e9*OrganicProperties(:,3))
! Get mole fraction of each organic in each droplet (including water!)
MoleFracs = DropMoleFrac(WaterMass)

! Get molar volumes of each organic in each droplet
do i = 1, ndrop
    do j = 1, norg
        MolVol(i) = MolVol(i) + MoleFracs(i,j)*OrganicProperties(j,3)
    end do
    ! Add inorganic
    !MolVol(i) = MolVol(i) + MoleFracs(i,norg+1)*InorganicBase(i,3)
    MolVol(i) = MolVol(i) + (1-SUM(MoleFracs(i,1:norg)))*InorganicBase(i,3)
    ! Add water
    !MolVol(i) = MolVol(i) + MoleFracs(i,norg+2)*Mw
end do
! Divide by droplet average density
MolVol = MolVol/SoluteProperties(:,4)
! Kelvin effect
KelFactor = EXP(2*DropSurfTens*MolVol/(R*Tmp*WaterRad))
! Equilibrium pressure for each organic
do i = 1, norg
    ! Equilibrium pressure in atm
    EquiPress(:,i) = MoleFracs(:,i)*Psat(i)*KelFactor*1 ! 1 = activity coefficient
    ! Equilibrium concentration at particle surface in kg/m^3
    Csurf(:,i) = EquiPress(:,i)*101325*OrganicProperties(i,3)/(R*Tmp)
end do

! Organic diffusivities
DiffOrgs = OrganicDiffusivity()

do i = 1, norg
    ! Adjusted for mass accommodation effects, alpha=0.1 (Sahle et al 2013)
    Dact = DiffOrgs(i)/(WaterRad/(WaterRad + (1.12e-7)/2) + (DiffOrgs(i)/(0.1*WaterRad))*SQRT(2*pi*OrganicProperties(i,3)/(R*Tmp)))
    ! Growth rate per droplet, in kg/s
    dOrgdt(:,i) = Dact*(OrganicGas(i) - Csurf(:,i))/(WaterRad*OrganicProperties(i,4))
end do

contains

!==========================================================================================================================

function SatConcentration(C0, Tmp) result(C0tadj)

implicit none

real(dvode_wp), intent(in)                      :: Tmp
real(dvode_wp), dimension(norg), intent(in)     :: C0
real(dvode_wp)                                  :: Tref
real(dvode_wp), dimension(norg)                 :: C0tadj, DelH

Tref = 298.15
!Delh_vap, enthalpy of vaporization term parameterization, in J/mol from Epstein et al (2009)
DelH = (131 - 11*LOG10(C0))*1e3
! New saturation concentration as a function of temperature
C0tadj = C0*(Tref/Tmp)*EXP(-(DelH/R)*(1/Tmp - 1/Tref))

end function SatConcentration

!==========================================================================================================================


function DropMoleFrac(WaterMass) result(MoleFracs)

implicit none

integer                                         :: i
real(dvode_wp), dimension(ndrop), intent(in)    :: WaterMass
real(dvode_wp), dimension(ndrop)                :: WaterMoles, InorganicMoles, TotalMoles
real(dvode_wp), dimension(ndrop, norg+2)        :: DropMoles, MoleFracs

DropMoles(:,norg+1) = InorganicBase(:,1)/InorganicBase(:,3)
DropMoles(:,norg+2) = WaterMass/Mw

do i = 1, norg
    DropMoles(:,i) = OrganicMass(:,i)/OrganicProperties(i,3)
end do

TotalMoles = SUM(DropMoles, 2)

do i = 1, norg+2
    MoleFracs(:,i) = DropMoles(:,i)/TotalMoles
end do

end function DropMoleFrac

!==========================================================================================================================

function OrganicDiffusivity() result(Dorg)

implicit none

real(dvode_wp), dimension(norg) :: Dorg

! Function to calculate the diffusivity of an organic aerosol using its molar mass, 
! taken from Lim, Ho-Jin, Annmarie G. Carlton, and Barbara J. Turpin (2005)
Dorg = (1.9*((OrganicProperties(:,3)*1e3)**(-2/3)))*1e-4
! m^2/s

end function OrganicDiffusivity

!==========================================================================================================================

end subroutine Cocondense

!==========================================================================================================================

end module SolFuncs