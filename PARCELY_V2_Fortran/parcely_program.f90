!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   The main program to run PARCELY V2.                                                *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Dan Barthaux                                                                       *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2025/11/04                                                      *
!*                                                                                      *
!****************************************************************************************

program PARCELY_V2

use dvode_module
use dvode_kinds_module,     only: wp => dvode_wp
use EnvironmentConstants,   only: pi, rhoL, Na, Mw
use InputReader,            only: EnvReader, InorganicReader, OrganicReader
use DropFuncs,              only: EquilibriumRadius, KappaKoehler, DropGrowthRate, DropletSurfaceTension
use SolFuncs,               only: SolPopInitialize
use ODEquations,            only: parcelydydt
use ModelParameters,        only: INCLUDE_ORGANICS, Cubes, ndrop, ninorg, norg, DropSurfTens, sftc, OrganicProperties, SoluteProperties

implicit none

type(dvode_t)               :: me
integer, allocatable        :: DistConcs(:), iwork(:)
real(wp), allocatable       :: DistRads(:), DistStds(:), InorganicProperties(:,:), EquilRad(:), y(:), rwork(:), rtol(:), atol(:), OrganicCond(:,:), OrganicGas(:)
real(wp), allocatable       :: OrganicOutput(:,:,:), OrganicGasOutput(:,:)
real(wp), allocatable       :: SaveData(:,:), SaveTimes(:), SaveSolutes(:,:,:), initmasswater(:)
real(wp)                    :: Tmp, RH, Prs, x, t, dt, RunTime, tout, Xi, ntsteps
integer                     :: npops, DistType, i, j, neq, istate, itask, iopt, mf, lrw, liw, itol, counter, Pts, ntint, DistSeed
integer, allocatable        :: seed(:)
character(50)               :: FormatString
character(6), allocatable   :: OrganicNames(:)
character(17+18)            :: filename

! Read in the environmental parameters of the parcel, as well as the aerosol population size information
call EnvReader(Prs, Tmp, RH, RunTime, dt, ntsteps, DistType, npops, DistConcs, DistRads, DistStds, Xi, DistSeed)
! Determine new seed if input random (-1)
allocate(seed(SUM(DistConcs)))
seed = DistSeed
if (DistSeed == 0) then
    call random_seed()
    call random_seed(get = seed)
    DistSeed = seed(1)
end if

! Read in the inorganic/constant-core properties for each particle
call InorganicReader(npops, InorganicProperties)
! Set number of droplets/solutes
ndrop = SUM(DistConcs)
! Set number of distinct inorganic/constant core types
ninorg = SIZE(InorganicProperties, 1)

! If organics included, read in properties of each unique organic
if (INCLUDE_ORGANICS) then
    call OrganicReader(OrganicNames, npops)
    ! Number of unique organics
    norg = SIZE(OrganicProperties, 1)
    allocate(OrganicCond(ndrop, norg))
    allocate(OrganicGas(norg))
    ! Set initial phases of organic mass to something small, molecules
    OrganicCond = 1.0e-20_wp
    OrganicGas = 1.0e-20_wp
    ! Determine number of equations for solver (# particles x # organics + water growth rate equations + # organics gas-phase + environment)
    neq = ndrop*(norg+1) + norg + 3
else
    ! No organics, simple case
    allocate(OrganicProperties(1,1))
    norg = 1
    allocate(OrganicCond(ndrop, 1))
    allocate(OrganicGas(1))
    OrganicCond = 0.0_wp
    OrganicGas = 0.0_wp
    neq = ndrop + 3
end if

! Finalize and create the inorganic/constant-core population
call SolPopInitialize(npops, DistType, DistConcs, DistRads, DistStds, Xi, InorganicProperties, OrganicCond, OrganicGas, DistSeed)

allocate(DropSurfTens(ndrop))
allocate(EquilRad(ndrop))

! Estimate the initial wet-radius (m) of the droplet assuming pure-water surface tension at temperature T
do i = 1, ndrop
    call EquilibriumRadius(x, Tmp, RH, SoluteProperties(i,2), SoluteProperties(i,5))
    EquilRad(i) = x
end do

! Set base surface tension (J/m^2) value for base-case (option 2 "Constant")
allocate(initmasswater(ndrop))
initmasswater = (4_wp*pi/3_wp)*rhoL*(EquilRad**3 - SoluteProperties(:,2)**3)
call DropletSurfaceTension(Tmp, initmasswater, OrganicCond, EquilRad)

! Allocate equation array and tolerances arrays
allocate(y(neq))
allocate(rtol(neq), atol(neq))

! Initial conditions
y(1:ndrop) = LOG(EquilRad)
y(neq-2) = Tmp
y(neq-1) = Prs
y(neq) = RH
! Water growth rate tolerances
rtol(1:ndrop) = 1.0e-7_wp
atol(1:ndrop) = 1.0e-8_wp
! Environmental tolerances (dTdt, dPdt, dRHdt)
rtol(neq-2) = 1e-8_wp
atol(neq-2) = 1e-9_wp
rtol(neq-1) = 1e-8_wp
atol(neq-1) = 1e-9_wp
rtol(neq) = 1e-8_wp
atol(neq) = 1e-9_wp
! Tolerances for organic growth rates if applicable
if (INCLUDE_ORGANICS) then
    rtol(ndrop+1:neq-3-norg) = 1.0e-6_wp
    atol(ndrop+1:neq-3-norg) = 1.0e-8_wp
    rtol(neq-3-norg+1:neq-3) = 1.0e-6_wp
    atol(neq-3-norg+1:neq-3) = 1.0e-8_wp
    ! Initial guess for organic concentration in every droplet
    y(ndrop+1:neq-3-norg) = RESHAPE(OrganicCond, [ndrop*norg])
    ! Amount of organics in the gas-phase, molecules/cm^3
    y(neq-3-norg+1:neq-3) = OrganicGas
end if

! Solver settings
itol = 1
istate = 1
itask = 1
iopt = 0
mf = 22 
lrw = 20 + neq*(5 + 1) + 3*neq + 2*neq**2 + 2
liw = 30 + neq

allocate(rwork(lrw))
allocate(iwork(liw))

iwork = 0

ntint = INT(ntsteps)

! Set starting time, s
t = 0.0_wp

! Allocate arrays to save data
allocate(SaveData(neq+1, ntint+2), SaveTimes(ntint+2), SaveSolutes(ndrop,6,ntint+2))
! Calculate number of iterations
do i = 1, ntint+2
    SaveTimes(i) = (i-1)*RunTime/ntsteps
end do

SaveData(neq+1,1) = t
SaveData(1:neq,1) = y
SaveSolutes(:,1:5,1) = SoluteProperties
SaveSolutes(:,6,1) = DropSurfTens
! Create solver instance
call me%initialize(parcelydydt)

print *, 'Starting integration...'

! Not to start at 1 and overwrite initial conditions
i = 2
counter = 2
do while (tout < RunTime)
    ! Set time point
    tout = min(t + dt, RunTime)
    ! Solve ODE system with DVODE
    call dvode(me, neq, y ,t, tout, itol, rtol, atol, itask, istate, iopt, &
               rwork, lrw, iwork, liw, mf)

    if (tout >= SaveTimes(counter)) then
        ! Save output to array
        SaveData(neq+1,counter) = tout
        SaveData(1:neq,counter) = y
        ! Save solute properties
        SaveSolutes(:,1:5,counter) = SoluteProperties
        ! Save surface tension
        SaveSolutes(:,6,counter) = DropSurfTens
        ! Update counter
        counter = min(counter + 1, ntint+2)
        print *, tout
    end if

    i = i + 1
end do

if (INCLUDE_ORGANICS) then
    allocate(OrganicOutput(ndrop, norg, ntint+2), OrganicGasOutput(norg, ntint+2))
    do i = 1, ntint+2
        OrganicOutput(:,:,i) = RESHAPE(SaveData(ndrop+1:neq-3-norg,i), [ndrop, norg])
        OrganicGasOutput(:,i) = SaveData(neq-3-norg+1:neq-3,i)
    end do
end if

! Open output files for saving data
! Environment
open(unit=10, file='../PARCELY_Output/environment_output.txt', status='replace', action='write')
write(10, '(A)') 'Time (s)   Temperature (K)   Pressure (Pa)   RH (-)'
! Droplet radius
open(unit=11, file='../PARCELY_Output/droplet_output.txt', status='replace', action='write')
write(11, '(A)') 'Time (s)   Droplet Radius (m) 1, 2, 3...'

! Solute properties
if (INCLUDE_ORGANICS) then
    open(unit=12, file='../PARCELY_Output/solmass_output.txt', status='replace', action='write')
    write(12, '(A)') 'Time (s)   Solute Mass (kg) 1, 2, 3...'
    open(unit=13, file='../PARCELY_Output/solradius_output.txt', status='replace', action='write')
    write(13, '(A)') 'Time (s)   Solute Radius (m) 1, 2, 3...'
    open(unit=14, file='../PARCELY_Output/solmolarmass_output.txt', status='replace', action='write')
    write(14, '(A)') 'Time (s)   Solute Molar Mass (kg/mol) 1, 2, 3...'
    open(unit=15, file='../PARCELY_Output/soldensity_output.txt', status='replace', action='write')
    write(15, '(A)') 'Time (s)   Solute Density (kg/m^3) 1, 2, 3...'
    open(unit=16, file='../PARCELY_Output/solkappa_output.txt', status='replace', action='write')
    write(16, '(A)') 'Time (s)   Solute Kappa (-) 1, 2, 3...'
    open(unit=17, file='../PARCELY_Output/organicgas_output.txt', status='replace', action='write')
    write(17, '(A)') 'Time (s)   Organic Gas-Phase (molec/cm^3) 1, 2, 3...'
    ! Organic mass
    do i = 1, norg
        filename = '../PARCELY_Output/' // OrganicNames(i) // '_output.txt'
        open(unit=i+18, file=filename, status='replace', action='write')
        write(i+18, '(A)') 'Time (s)   Organic Mass (molec) 1, 2, 3...'
    end do    
end if

open(unit=18, file='../PARCELY_Output/soltension_output.txt', status='replace', action='write')
write(18, '(A)') 'Time (s)   Solute Surface Tension (J/m^-2) 1, 2, 3...'

write(FormatString, '(A,I0,A)') '(F10.4,', ndrop, 'ES15.6)'

do i = 1, ntint+1
    write(10, '(F10.4, 3ES15.6)') SaveData(neq+1,i), SaveData(neq-2,i), SaveData(neq-1,i), SaveData(neq,i)
    write(11, FormatString) SaveData(neq+1,i), EXP(SaveData(1:ndrop,i))
    write(18, FormatString) SaveData(neq+1,i), SaveSolutes(:,6,i)

    if (INCLUDE_ORGANICS) then
        write(12, FormatString) SaveData(neq+1,i), SaveSolutes(:,1,i)
        write(13, FormatString) SaveData(neq+1,i), SaveSolutes(:,2,i)
        write(14, FormatString) SaveData(neq+1,i), SaveSolutes(:,3,i)
        write(15, FormatString) SaveData(neq+1,i), SaveSolutes(:,4,i)
        write(16, FormatString) SaveData(neq+1,i), SaveSolutes(:,5,i)
        write(17, FormatString) SaveData(neq+1,i), OrganicGasOutput(:,i)

        do j = 1, norg
            write(j+18, FormatString) SaveData(neq+1,i), OrganicOutput(:,j,i)
        end do
    end if
end do

! Close files
if (INCLUDE_ORGANICS) then
    do i = 10, 18+norg
        close(i)
    end do
else
    close(10)
    close(11)
    close(18)
end if

print *, 'Integration complete.'
print *, "Press ENTER to exit."
read(*,*)

end program PARCELY_V2