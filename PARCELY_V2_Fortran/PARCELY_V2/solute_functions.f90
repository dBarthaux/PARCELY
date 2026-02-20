module SolFuncs

use dvode_kinds_module, only: wp => dvode_wp
use EnvironmentConstants, only: pi, R, Na, R_alt, Mw, rhoL
use ModelParameters, only: norg, ninorg, ndrop, OrganicProperties, InorganicBase, SoluteProperties, DropSurfTens, Cubes

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

function SoluteRadiusDistribution(npops, DistConcs, DistRads, DistStds, DistType, DistSeed) result(DryRad)

use DistFuncs, only: NormDist, LogNormDist

implicit none

integer, intent(in)                             :: npops        ! Number of different distributions
integer, intent(in)                             :: DistSeed     ! Seed for random number generation
integer, dimension(npops), intent(in)           :: DistConcs    ! Number concentration of each distribution
real(wp), dimension(npops), intent(in)          :: DistRads     ! Distribution means in um
real(wp), dimension(npops), intent(in)          :: DistStds     ! Distribution standard deviations in um
integer, intent(in)                             :: DistType     ! Type of distribution used [1: monodispersed, 2:normal, 3:log-normal]
real(wp), allocatable                           :: DryRad(:)    ! Dry/solute/aerosol radii, m
real(wp), allocatable                           :: randnum(:)
integer                                         :: i, j
integer, allocatable                            :: seed(:)

allocate(DryRad(SUM(DistConcs)), seed(SUM(DistConcs)))
DryRad = 1.0_wp

! If seed is -1, random seed used to get distribution
seed = DistSeed
call random_seed(put = seed)

j = 0

do i = 1, npops
    
    allocate(randnum(DistConcs(i)))

    ! Mono-disperse distribution
    if (DistType == 1) then
        if (i == 1) then
            DryRad(1:DistConcs(i)) = DryRad(1:DistConcs(i))*DistRads(i)
        else
            !j = SUM(DistConcs(1:i-1))
            j = j + DistConcs(i)
            DryRad(j+1:j+DistConcs(i)) = DryRad(j+1:j+DistConcs(i))*DistRads(i)
        endif

    ! Normal distribution
    elseif (DistType == 2) then
        if (i == 1) then
            call NormDist(DistRads(i), DistStds(i), DistConcs(i), randnum)
            DryRad(1:DistConcs(i)) = randnum
        else
            !j = SUM(DistConcs(1:i-1))
            j = j + DistConcs(i)
            call NormDist(DistRads(i), DistStds(i), DistConcs(i), randnum)
            DryRad(j+1:j+DistConcs(i)) = randnum
        endif

    ! Log-normal distribution
    else
        if (i == 1) then
            call LogNormDist(DistRads(i), DistStds(i), DistConcs(i), randnum)
            DryRad(1:DistConcs(i)) = randnum
        else
            !j = SUM(DistConcs(1:i-1))
            j = j + DistConcs(i)
            call LogNormDist(DistRads(i), DistStds(i), DistConcs(i), randnum)
            DryRad(j+1:j+DistConcs(i)) = randnum
        endif

    endif
    
    deallocate(randnum)

end do

! Convert to meters
DryRad = DryRad*1e-6_wp

end function SoluteRadiusDistribution

!==========================================================================================================================

function InitialOrganics(r, Mi, rhoi, j) result(NewState)

implicit none

integer, intent(in)                     :: j
real(wp), intent(in)                    :: r, Mi, rhoi ! Dry radius [m], inorganic molar mass [kg/mol], inorganic density [kg/m^3], 
real(wp), dimension(norg)               :: Mo, rhoo, no, masso ! Organic molar mass [kg/mol], density [kg/m^3], moles [-], mass [kg]
real(wp), dimension(norg)               :: wo ! Organic mass fractions
real(wp)                                :: mt, massi, ni, nt, rnew, scaler, wi
real(wp), parameter                     :: initguess = 1e-18 ! Initial guess for total mass, kg
logical                                 :: Match
integer                                 :: i
real(wp), dimension(2+norg)             :: NewState

! Assign names for ease
Mo = OrganicProperties(:,3)
rhoo = OrganicProperties(:,4)
wo = OrganicProperties(:,6+j)
! Get the initial inorganic mass fraction
wi = 1.0_wp - SUM(wo)

! Initialize values
mt = initguess
Match = .FALSE.
scaler = 1.0_wp
NewState = 0.0_wp

do while (.not. Match)
    ! Inorganic mass
    massi = mt*wi
    ! Inorganic moles
    ni = massi/Mi
    ! Organic mass
    masso = mt*wo
    ! Organic moles
    no = masso/Mo
    ! Total moles
    nt = ni + SUM(no)
    ! Corresponding radius
    rnew = ((3.0_wp/(4.0_wp*pi))*((rhoi*ni/nt + SUM(rhoo*(no/nt)))**(-1.0_wp))*(ni*Mi + SUM(no*Mo)))**(1.0_wp/3.0_wp)

    if (abs(rnew*1.0e8_wp - r*1.0e8_wp) < 1.0e-6_wp) then
        Match = .TRUE.
        NewState(1) = rnew
        NewState(2) = ni
        do i = 1, norg
            NewState(i+2) = no(i)
        end do

    else
        ! Scale new total mass by ratio of corresponding radius with target radius
        scaler = r/rnew
        mt = scaler*mt
    end if
end do

end function InitialOrganics

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

subroutine InitialSoluteProperties(npops, Xi, DistConcs, DryRad, InorganicProperties, OrganicCond, OrganicGas, DistSeed)

use ModelParameters, only: INCLUDE_ORGANICS

implicit none

integer, intent(in)                                     :: npops                ! Number of aerosol populations (different size distributions)
real(wp), intent(in)                                    :: Xi                   ! Mixing state
integer, dimension(npops), intent(in)                   :: DistConcs            ! Number concentration of each distribution
real(wp), dimension(ndrop), intent(in)                  :: DryRad               ! Dry/solute/aerosol radii, m
real(wp), dimension(ninorg, npops+3), intent(in)        :: InorganicProperties
integer, intent(in)                                     :: DistSeed             ! Seed for random number generation
real(wp), dimension(ndrop, norg), intent(inout)         :: OrganicCond
real(wp), dimension(norg), intent(inout)                :: OrganicGas

integer                                                 :: i, n, j, NumPop
real(wp)                                                :: TotalMoles, wi
real(wp), allocatable                                   :: SolKappa(:), SolDensity(:), SolMolarMass(:), InorganicCounts(:), SolMass(:), SolMoles(:), TotalVolume(:)
real(wp), allocatable                                   :: InorganicHomes(:), InorganicAways(:), InorganicMoles(:), InorganicPopFracs(:)
real(wp), allocatable                                   :: InorganicVolumes(:,:), InorganicVolFracs(:,:), InorganicMasses(:,:), InorganicMoleFracs(:), SolMoleFracs(:,:)
real(wp), allocatable                                   :: OrganicAdjustment(:,:), TotalOrganics(:)
real(wp), dimension(ninorg)                             :: InorganicKappas, InorganicDensities, InorganicMolarMasses
real(wp), dimension(ninorg, npops)                      :: InorganicFractions
character(len=20)                                       :: seedstring

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
        ! In cases where the inorganic fractions do not match the total number of particles (e.g., two inorganics both 50% of an odd number)
        if (n+NumPop-1 > ndrop) then
            NumPop = ndrop - (n - 1)
        end if
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
    SolMass = (4_wp*pi/3_wp)*SolDensity*(DryRad**3)

! Otherwise, solute average properties are mole fraction-weighted properties
else
    ! Start with initial mass in moles
    SolMoles = (4_wp*pi/3_wp)*SolDensity*(DryRad**3)/(SolMolarMass)
    ! Get the total number of moles
    TotalMoles = SUM(SolMoles)
    ! Get the moles for each inorganic
    do i = 1, ninorg
        InorganicMoles(i) = SUM(SolMoles, MASK = (InorganicCounts == i))
    end do
    ! Get the total mole fractions for each inorganic
    InorganicMoleFracs = InorganicMoles/TotalMoles
    ! Define the Home and Away terms for each inorganic
    InorganicHomes = 1_wp + Xi*(InorganicMoleFracs - 1_wp)
    InorganicAways = Xi*InorganicMoleFracs
    ! Assign the mole fraction of each inorganic in every droplet
    do i = 1, ninorg
        do n = 1, ndrop
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
    do i = 1, ninorg
        do n = 1, ndrop
            InorganicMasses(n,i) = (4_wp*pi/3_wp)*(SolMoleFracs(n,i)*InorganicDensities(i))*(DryRad(n)**3)
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
    SolMass = (4_wp*pi/3_wp)*SolDensity*(DryRad**3)

end if

write(seedstring, '(I0)') DistSeed
seedstring = trim(seedstring)

! If organics are included and they have a non-zero total initial mass fraction
if (INCLUDE_ORGANICS .AND. SUM(OrganicProperties(:,7:)) > 0.0_wp) then
    
    allocate(OrganicAdjustment(ndrop, 2+norg))

    ! Find the corresponding moles for inorganic and organics for the same dry radius
    j = 1
    n = 1
    do i = 1, ndrop
        OrganicAdjustment(i,:) = InitialOrganics(DryRad(i), SolMolarMass(i), SolDensity(i), j)
        n = n + 1
        if (DistConcs(j) == n-1) then
            j = j + 1
            n = 1
        end if
    end do

    ! Reset the relevant inorganic properties
    allocate(InorganicBase(ndrop, 5))
    InorganicBase(:,1) = OrganicAdjustment(:,2)*SolMolarMass
    InorganicBase(:,2) = ((3_wp/(4_wp*pi))*InorganicBase(:,1)/SolDensity)**(1.0_wp/3.0_wp)
    InorganicBase(:,3) = SolMolarMass
    InorganicBase(:,4) = SolDensity
    InorganicBase(:,5) = SolKappa

    ! Create inorganics file
    open(unit=12, file='../PARCELY_Output/inorganics_output.txt', status='replace', action='write')
    write(12, '(A)') 'Mass (kg)   Radius (m)   Molar Mass (kg/mol)   Density (kg/m3)   Kappa (-)   Seed:'//seedstring
    do i = 1,ndrop
        write(12, '(5ES15.6)') InorganicBase(i,1), InorganicBase(i,2), InorganicBase(i,3), InorganicBase(i,4), InorganicBase(i,5)
    end do
    close(12)

    ! Get the organic mass in molecules
    do n = 1, norg
        OrganicCond(:,n) = OrganicAdjustment(:,2+n)*Na
    end do
        
    ! Replace zeros with low initial value
    where (OrganicCond == 0.0_wp)
        OrganicCond = 1.0e-20_wp
    end where

    ! SoluteProperties: mass (kg), radius (m), molar mass (kg/mol), density (kg/m3), kappa
    allocate(SoluteProperties(ndrop, 5))
    ! Update solute
    call SoluteUpdate(OrganicCond)

    ! Get the corresponding gas-phase organic concentrations, molec/cm3
    allocate(TotalOrganics(norg))
    TotalOrganics = ((OrganicProperties(:,1)*1e-6_wp)/(OrganicProperties(:,3)*1e3_wp))*1e-6_wp*Na
    do n = 1, norg
        ! Make sure the initial condensed mass of each organic does not exceed their respective total mass
        OrganicGas(n) = MAX(TotalOrganics(n) - SUM(OrganicCond(:,n)), 1.0e-20_wp)
    end do

else
    ! SoluteProperties: mass (kg), radius (m), molar mass (kg/mol), density (kg/m3), kappa
    allocate(SoluteProperties(ndrop, 5))
    SoluteProperties(:,1) = SolMass
    SoluteProperties(:,2) = DryRad
    SoluteProperties(:,3) = SolMolarMass
    SoluteProperties(:,4) = SolDensity
    SoluteProperties(:,5) = SolKappa

    ! Create inorganics file
    open(unit=12, file='../PARCELY_Output/inorganics_output.txt', status='replace', action='write')
    write(12, '(A)') 'Mass (kg)   Radius (m)   Molar Mass (kg/mol)   Density (kg/m3)   Kappa (-)   Seed:'//seedstring
    do i = 1,ndrop
        write(12, '(5ES15.6)') SolMass(i), DryRad(i), SolDensity(i), SolMolarMass(i), SolKappa(i)
    end do
    close(12)

    allocate(InorganicBase(ndrop, 5))
    InorganicBase = SoluteProperties
end if

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

subroutine SolPopInitialize(npops, DistType, DistConcs, DistRads, DistStds, Xi, InorganicProperties, OrganicCond, OrganicGas, DistSeed)

implicit none

integer, intent(in)                                     :: npops, DistType, DistSeed
real(wp), intent(in)                                    :: Xi                           ! Mixing State
integer, dimension(npops), intent(in)                   :: DistConcs                    ! Number concentration of each distribution
real(wp), dimension(npops), intent(in)                  :: DistRads, DistStds
real(wp), dimension(ninorg, 3+npops), intent(in)        :: InorganicProperties
real(wp), allocatable                                   :: DryRad(:)
real(wp), dimension(ndrop, norg), intent(inout)         :: OrganicCond
real(wp), dimension(norg), intent(inout)                :: OrganicGas


! Get dry radii
DryRad = SoluteRadiusDistribution(npops, DistConcs, DistRads, DistStds, DistType, DistSeed)
!Set any radius smaller than 5 nm to 5 nm
WHERE (DryRad < 5e-9_wp)
    DryRad = 5e-9_wp
END WHERE

ndrop = SIZE(DryRad)

! Get solute properties
call InitialSoluteProperties(npops, Xi, DistConcs, DryRad, InorganicProperties, OrganicCond, OrganicGas, DistSeed)

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

subroutine SoluteUpdate(OrganicCond)

implicit none

real(wp), dimension(ndrop, norg), intent(in)    :: OrganicCond
real(wp), dimension(ndrop, norg)                :: OrgUg
integer                                         :: i, j
real(wp), dimension(ndrop)                      :: InorganicMoles, InorganicVolumes, InorganicMoleFrac, InorganicVolFrac, TotalMoles, TotalVolume
real(wp), dimension(ndrop, norg)                :: OrganicMoles, OrganicMoleFrac, OrganicVolumes, OrganicVolFrac

! Update solute mass
! Convert organics from molecules to ug
OrgUg = 0.0_wp
do i = 1, norg
    OrgUg(:,i) = 1.0e9_wp*OrganicProperties(i,3)*OrganicCond(:,i)/Na
end do
SoluteProperties(:,1) = InorganicBase(:,1) + SUM(OrgUg, 2)*1.0e-9_wp

! Calculate moles
InorganicMoles = InorganicBase(:,1)/InorganicBase(:,3)
do i = 1, norg
    OrganicMoles(:,i) = OrgUg(:,i)*1e-9_wp/OrganicProperties(i,3)
end do

! Get the total moles
TotalMoles = InorganicMoles + SUM(OrganicMoles, 2)

! Mole fractions
InorganicMoleFrac = InorganicMoles/TotalMoles
do i = 1, ndrop
    OrganicMoleFrac(i,:) = OrganicMoles(i,:)/TotalMoles(i)
end do

! Solute molar mass and density are mole-fraction weighted averages
SoluteProperties(:,3) = InorganicMoleFrac*InorganicBase(:,3) + MATMUL(OrganicMoleFrac, OrganicProperties(:,3))
SoluteProperties(:,4) = InorganicMoleFrac*InorganicBase(:,4) + MATMUL(OrganicMoleFrac, OrganicProperties(:,4))

! Update solute radius
SoluteProperties(:,2) = ((3_wp/(4_wp*pi))*SoluteProperties(:,1)/SoluteProperties(:,4))**(1.0_wp/3.0_wp)

! Get the volumes for each component
InorganicVolumes = InorganicBase(:,1)/InorganicBase(:,4)
do i = 1, ndrop
    OrganicVolumes(i,:) = OrgUg(i,:)*1e-9_wp/OrganicProperties(:,4)
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

subroutine Cocondense(OrganicCond, OrganicGas, WaterMass, WaterRad, Tmp, dOrgConddt, dOrgGasdt)

implicit none

real(wp), dimension(ndrop, norg), intent(in)    :: OrganicCond
real(wp), dimension(norg), intent(in)           :: OrganicGas
real(wp), dimension(ndrop), intent(in)          :: WaterMass, WaterRad
real(wp), intent(in)                            :: Tmp                  ! Ambient temperature, K
real(wp), dimension(ndrop, norg), intent(out)   :: dOrgConddt
real(wp), dimension(norg), intent(out)          :: dOrgGasdt
integer                                         :: i, j
real(wp), allocatable                           :: KelFactor(:,:), EquiPress(:,:), Csurf(:,:), MoleFracs(:,:), C0tadj(:), &
                                                   DiffOrgs(:), Psat(:), DelH(:), ThermVel(:), FreePath(:), Knudsen(:,:), &
                                                   KnudInverse(:,:), Corr1(:,:), Corr2(:,:), Corr3(:,:), Correction(:,:), &
                                                   InorgMolec(:), TotalSolMolec(:), WaterMolec(:), DropMw(:), DropDensity(:)

allocate(KelFactor(ndrop, norg), EquiPress(ndrop, norg), Csurf(ndrop, norg), MoleFracs(ndrop, norg+2), C0tadj(norg), DiffOrgs(norg), &
        Psat(norg), DelH(norg), ThermVel(ndrop), FreePath(norg), Knudsen(ndrop, norg), KnudInverse(ndrop, norg), Corr1(ndrop, norg), &
        Corr2(ndrop, norg), Corr3(ndrop, norg), Correction(ndrop, norg), InorgMolec(ndrop), TotalSolMolec(ndrop), &
        DropMw(ndrop), DropDensity(ndrop))

! Number of molecules in inorganic phase
InorgMolec = Na*InorganicBase(:,1)/InorganicBase(:,3)
! Number of water molecules
WaterMolec = Na*WaterMass/Mw
! Total molecules in particle phase for each droplet
TotalSolMolec = InorgMolec + WaterMolec + SUM(OrganicCond, 2)
! Mole fractions using this definition
MoleFracs = 0.0
do i = 1, norg
    MoleFracs(:,i) = OrganicCond(:,i)/TotalSolMolec
end do
MoleFracs(:,norg+1) = InorgMolec/TotalSolMolec
MoleFracs(:,norg+2) = WaterMolec/TotalSolMolec

! Organic diffusivities
DiffOrgs = 1.9_wp*((OrganicProperties(:,3)*1.0e3_wp)**(-2.0_wp/3.0_wp))
! Mean thermal velocity of each molecule
ThermVel = SQRT((8.0_wp*R*Tmp)/(pi*OrganicProperties(:,3)))
! Mean free path for each molecule, m
FreePath = 3.0_wp*DiffOrgs*1.0e-4_wp/ThermVel
! Calculate Knudsen number
do i = 1, norg
    Knudsen(:,i) = FreePath(i)/SoluteProperties(:,2)
end do
! Calculate the corrections
KnudInverse = 1_wp/Knudsen
Corr1 = (1.33_wp + 0.71_wp*KnudInverse)/(1_wp + KnudInverse)
Corr2 = (4.0_wp*(1.0_wp - 0.1_wp))/(3.0_wp*0.1_wp) !alpha=0.1 (Saleh et al 2013)
Corr3 = 1.0_wp + (Corr1 + Corr2)*Knudsen
Correction = 1.0_wp/Corr3
!Delh_vap, enthalpy of vaporization term parameterization, in J/mol from Epstein et al (2009)
DelH = (131.0_wp - 11.0_wp*LOG10(OrganicProperties(:,2)))*1.0e3_wp
! New saturation concentration as a function of temperature
C0tadj = OrganicProperties(:,2)*(298.15_wp/Tmp)*EXP(-(DelH/R)*(1.0_wp/Tmp - 1.0_wp/298.15_wp))
! Convert C0 to a pure component saturation vapour pressure [atm]
Psat = (C0tadj*R*Tmp)/(101325.0_wp*OrganicProperties(:,3)*1.0e9_wp)

! Equilibrium pressure for each organic
do i = 1, norg
    ! Kelvin effect
    KelFactor(:,i) = EXP(2.0_wp*DropSurfTens*OrganicProperties(i,3)/(R*Tmp*WaterRad*OrganicProperties(i,4)))
    ! Equilibrium pressure in atm
    EquiPress(:,i) = MoleFracs(:,i)*Psat(i)*KelFactor(:,i)*1.0_wp ! 1 = activity coefficient
    ! Equilibrium concentration at particle surface in molec/cm^3
    Csurf(:,i) = EquiPress(:,i)*(Na/(R_alt*1.0e6_wp*Tmp))
    ! Mass transfer
    dOrgConddt(:,i) = 4.0_wp*pi*SoluteProperties(:,2)*1.0e2_wp*DiffOrgs(i)*Correction(:,i)*(OrganicGas(i) - Csurf(:,i))
    dOrgGasdt(i) = -SUM(dOrgConddt(:,i))
end do

end subroutine Cocondense

!==========================================================================================================================

end module SolFuncs
