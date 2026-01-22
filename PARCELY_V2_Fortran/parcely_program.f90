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
use dvode_kinds_module,     only: dvode_wp
use EnvironmentConstants,   only: pi, rhoL, Na, Mw
use InputReader,            only: EnvReader, InorganicReader, OrganicReader
use DropFuncs,              only: EquilibriumRadius, KappaKoehler, DropGrowthRate
use SolFuncs,               only: SolPopInitialize
use ODEquation,             only: parcelydydt
use ModelParameters,        only: INCLUDE_ORGANICS, ndrop, ninorg, norg, DropSurfTens, sftc, OrganicProperties, SoluteProperties

implicit none

type(dvode_t)               :: me
integer, allocatable        :: DistConcs(:), iwork(:)
real(dvode_wp), allocatable :: DistRads(:), DistStds(:), InorganicProperties(:,:), EquilRad(:), y(:), rwork(:), rtol(:), atol(:), OrganicMass(:,:), SaveData(:,:)
real(dvode_wp)              :: Tmp, RH, Prs, x, t, dt, RunTime, tout, Xi
integer                     :: npops, DistType, i, j, neq, istate, itask, iopt, mf, lrw, liw, itol, Pts
character(50)               :: FormatString
character(6), allocatable   :: OrganicNames(:)
character(17)               :: filename

! Read in the environmental parameters of the parcel, as well as the aerosol population size information
call EnvReader(Prs, Tmp, RH, RunTime, dt, DistType, npops, DistConcs, DistRads, DistStds, Xi)
! Read in the inorganic/constant-core properties for each particle
call InorganicReader(npops, InorganicProperties)
! Set number of droplets/solutes
ndrop = SUM(DistConcs)
! Set number of distinct inorganic/constant core types
ninorg = SIZE(InorganicProperties, 1)
! Finalize and create the inorganic/constant-core population
call SolPopInitialize(npops, DistType, DistConcs, DistRads, DistStds, Xi, InorganicProperties)

allocate(DropSurfTens(ndrop))
allocate(EquilRad(ndrop))

! Estimate the initial wet-radius of the droplet assuming pure-water surface tension at temperature T
do i = 1, ndrop
    call EquilibriumRadius(x, Tmp, RH, SoluteProperties(i,2), SoluteProperties(i,5))
    EquilRad(i) = x
end do

! Set base surface tension value for base-case (option 2 "Constant")
DropSurfTens = sftc

! If organics included, read in properties of each unique organic
if (INCLUDE_ORGANICS) then
    call OrganicReader(OrganicNames)
    ! Number of unique organics
    norg = SIZE(OrganicProperties, 1)
    allocate(OrganicMass(ndrop, norg))
    ! Set initial condensed organic mass to 0
    OrganicMass = 0.0
    ! Determine number of equations for solver (# particles x # organics + water growth rate equations)
    neq = ndrop*(norg+1) + 3
else
    ! No organics, simple case
    allocate(OrganicProperties(1,1))
    norg = 1
    allocate(OrganicMass(ndrop, 1))
    OrganicMass = 0.0
    neq = ndrop + 3
end if

! Allocate equation array and tolerances arrays
allocate(y(neq))
allocate(rtol(neq), atol(neq))

! Initial conditions
y(1:ndrop) = EquilRad
y(neq-2) = Tmp
y(neq-1) = Prs
y(neq) = RH
! Water growth rate tolerances
rtol(1:ndrop) = 1.0e-16
atol(1:ndrop) = 1.0e-18
! Environmental tolerances (dTdt, dPdt, dRHdt)
rtol(neq-2) = 1e-6
atol(neq-2) = 1e-8
rtol(neq-1) = 1e-6
atol(neq-1) = 1e-8
rtol(neq) = 1e-6
atol(neq) = 1e-8
! Tolerances for organic growth rates if applicable
if (INCLUDE_ORGANICS) then
    rtol(ndrop+1:neq-3) = 1.0e-20
    atol(ndrop+1:neq-3) = 1.0e-22
    y(ndrop+1:neq-3) = 0.0
end if
! Set starting time
t = 0.0
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

! Calculate number of time instances will be saved with the given timestep size and run time
Pts = RunTime/dt
allocate(SaveData(neq+1, Pts+2))

! Create solver instance
call me%initialize(parcelydydt)

print *, 'Starting integration...'

! Save initial conditions
SaveData(1,neq+1) = t
SaveData(1:neq,1) = y

! Not to start at 1 and overwrite initial conditions
i = 2
do while (t < RunTime)
    ! Set time point
    tout = min(t + dt, RunTime)
    ! Solve ODE system with DVODE
    call dvode(me, neq, y ,t, tout, itol, rtol, atol, itask, istate, iopt, &
               rwork, lrw, iwork, liw, mf)
    
    !if (istate /= 2) then
    !    print *, "You moron."
    !    read(*,*)
    !end if

    !write(*,'(A)',advance='no') t
    
    ! Save output to array
    SaveData(neq+1,i) = t
    SaveData(1:neq,i) = y
    i = i + 1
end do

! Open output files for saving data
open(unit=10, file='environment_output.txt', status='replace', action='write')
write(10, '(A)') 'Time (s)   Temperature (K)   Pressure (Pa)   RH (-)'

open(unit=11, file='droplet_output.txt', status='replace', action='write')
write(11, '(A)') 'Time (s)   Droplet Radius (m) 1, 2, 3...'

if (INCLUDE_ORGANICS) then
    do i = 1, norg
        filename = OrganicNames(i) // '_output.txt'
        open(unit=i+11, file=filename, status='replace', action='write')
        write(i+11, '(A)') 'Time (s)   Organic Mass (kg) 1, 2, 3...'
    end do
end if

write(FormatString, '(A,I0,A)') '(F10.4,', ndrop, 'ES15.6)'

do i = 1, Pts+1
    write(10, '(F10.4, 3ES15.6)') SaveData(neq+1,i), SaveData(neq-2,i), SaveData(neq-1,i), SaveData(neq,i)
    write(11, FormatString) SaveData(neq+1,i), SaveData(1:ndrop,i)

    if (INCLUDE_ORGANICS) then
        do j = 1, norg
            write(j+11, FormatString) SaveData(neq+1,i), SaveData(ndrop*j+j:ndrop*j+j+ndrop-1,i)
        end do
    end if    

end do
close(10)
close(11)
if (INCLUDE_ORGANICS) then
    do i = 12, 11+norg
        close(i)
    end do
end if

print *, 'Integration complete.'
print *, "Press ENTER to exit."
read(*,*)

end program PARCELY_V2