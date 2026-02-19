module InputReader

use dvode_kinds_module, only: wp => dvode_wp
use ModelParameters, only: Cubes, W, ac, at, sftc, Msft, OrganicProperties, INCLUDE_ORGANICS, INCLUDE_COCONDENSE

contains

!==========================================================================================================================

!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Function to read in the environmental/distributions/solver input file.             *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Dan Barthaux                                                                       *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2025/11/04                                                      *
!*                                                                                      *
!****************************************************************************************

subroutine EnvReader(Prs, Tmp, RH, RunTime, dt, ntsteps, DistType, npops, DistConcs, DistRads, DistStds, Xi, DistSeed)

implicit none

character(len=256)              :: line
real(wp)                        :: Prs, Tmp, RH, RunTime, dt, Xi, ntsteps
integer                         :: DistType, npops, ios, DistSeed
integer, allocatable            :: DistConcs(:)
real(wp), allocatable           :: DistRads(:), DistStds(:)

open(unit=10, file="parcely_env_input_file.txt", status="old", action="read")

do
    read(10,'(A)', iostat=ios) line
    if (ios /= 0) exit  ! End of file
    if (line(1:1) == '=' .or. line(1:1) == '-' .or. line(1:1) == 'V' .or. line(1:1) == 'E') cycle
    if (trim(line) == '') cycle

    if (index(line, '[npops]') > 0) then
        read(line, *) npops
        allocate(DistConcs(npops), DistRads(npops), DistStds(npops))
    else if (index(line, '[DistSeed]') > 0) then
        read(line, *) DistSeed
    else if (index(line, '[Cubes]') > 0) then
        read(line, *) Cubes
    else if (index(line, '[Prs]') > 0) then
        read(line, *) Prs
    else if (index(line, '[Tmp]') > 0) then
        read(line, *) Tmp
    else if (index(line, '[RH]') > 0) then
        read(line, *) RH
    else if (index(line, '[W]') > 0) then
        read(line, *) W
    else if (index(line, '[ac]') > 0) then
        read(line, *) ac
    else if (index(line, '[at]') > 0) then
        read(line, *) at
    else if (index(line, '[RunTime]') > 0) then
        read(line, *) RunTime
    else if (index(line, '[dt]') > 0) then
        read(line, *) dt
    else if (index(line, '[ntsteps]') > 0) then
        read(line, *) ntsteps
    else if (index(line, '[DistType]') > 0) then
        read(line, *) DistType
    else if (index(line, '[DistConcs]') > 0) then
        read(line, *) DistConcs
    else if (index(line, '[DistRads]') > 0) then
        read(line, *) DistRads
    else if (index(line, '[DistStds]') > 0) then
        read(line, *) DistStds
    else if (index(line, '[Xi]') > 0) then
        read(line, *) Xi
    else if (index(line, '[Msft]') > 0) then
        read(line, *) Msft
    else if (index(line, '[sftc]') > 0) then
        read(line, *) sftc
    else if (index(line, '[INCLUDE_ORGANICS]') > 0) then
        read(line, *) INCLUDE_ORGANICS
    else if (index(line, '[INCLUDE_COCONDENSE]') > 0) then
        read(line, *) INCLUDE_COCONDENSE
    end if
end do

! Multiply concentration by number of cubic centimeters
DistConcs = DistConcs*Cubes

close(10)

end subroutine EnvReader

!==========================================================================================================================

!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Function to read in the inorganic inputs.                                          *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Dan Barthaux                                                                       *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2025/11/04                                                      *
!*                                                                                      *
!****************************************************************************************

subroutine InorganicReader(npops, InorganicProperties)

implicit none

integer, intent(in)                         :: npops
character(len=256)                          :: line
character(len=6)                            :: name
integer                                     :: ios, i, ninorg, j
real(wp), allocatable, intent(out)          :: InorganicProperties(:,:)

open(unit=11, file="parcely_inorganic_input_file.txt", status="old", action="read")

! Need to count how many inorganics before being able to read correctly
ninorg = 0
do
    read(11,'(A)', iostat=ios) line
    if (ios /= 0) exit
    if (line(1:1) == '*') exit
    if (line(1:1) == '=' .or. line(1:9) == 'Inorganic' .or. line(1:4) == 'Name') cycle
    if (trim(line) == '') cycle
    ninorg = ninorg + 1
end do
close(11)

allocate(InorganicProperties(ninorg, 3+npops))

open(unit=11, file="parcely_inorganic_input_file.txt", status="old", action="read")

i = 1
do
    read(11,'(A)', iostat=ios) line
    if (line(1:1) == '*') exit
    if (line(1:1) == '=' .or. line(1:9) == 'Inorganic' .or. line(1:4) == 'Name') cycle
    if (trim(line) == '') cycle
    
    ! InorganicProperties order: MolarMass (g/mol), Density (kg/m3), Kappa, PopFrac_n, PopFrac_n+1...
    read(line, *, iostat=ios) name, (InorganicProperties(i, j+3), j=1,npops), &
                                   InorganicProperties(i,1), InorganicProperties(i,2), InorganicProperties(i,3)
    i = i + 1
end do

close(11)

! Convert molar mass from g/mol to kg/mol
InorganicProperties(:,1) = InorganicProperties(:,1)*1e-3_wp

end subroutine InorganicReader

!==========================================================================================================================

!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Function to read in the inorganic inputs.                                          *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Dan Barthaux                                                                       *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2025/11/04                                                      *
!*                                                                                      *
!****************************************************************************************

subroutine OrganicReader(OrganicNames, npops)

implicit none

integer, intent(in)                     :: npops
character(len=256)                      :: line
character(len=6)                        :: name
integer                                 :: ios, norg, i, j
character(6), allocatable, intent(out)  :: OrganicNames(:)

open(unit=11, file="parcely_organic_input_file.txt", status="old", action="read")

! Need to count how many organics before being able to read correctly
norg = 0
do
    read(11,'(A)', iostat=ios) line
    if (ios /= 0) exit
    if (line(1:1) == '*') exit
    if (line(1:1) == '=' .or. line(1:7) == 'Organic' .or. line(1:4) == 'Name') cycle
    if (trim(line) == '') cycle
    norg = norg + 1
end do
close(11)

allocate(OrganicProperties(norg, 6+npops))
allocate(OrganicNames(norg))

open(unit=11, file="parcely_organic_input_file.txt", status="old", action="read")

i = 1
do
    read(11,'(A)', iostat=ios) line
    if (line(1:1) == '*') exit
    if (line(1:1) == '=' .or. line(1:7) == 'Organic' .or. line(1:4) == 'Name') cycle
    if (trim(line) == '') cycle
    ! Total concentration (ug/m^3), C0 (ug/m^3), Molar mass (g/mol), Density (kg/m^3), Kappa, Surface tension (mJ/m^2), Initial average particle mass fraction
    read(line, *, iostat=ios) OrganicNames(i), OrganicProperties(i,1), OrganicProperties(i,2), OrganicProperties(i,3), & 
                                    OrganicProperties(i,4), OrganicProperties(i,5), OrganicProperties(i,6), &
                                    (OrganicProperties(i, j+6), j=1,npops)
    i = i + 1
end do

close(11)

! Convert molar mass from g/mol to kg/mol
OrganicProperties(:,3) = OrganicProperties(:,3)*1e-3_wp
! Convert surface tension from mJ/m^2 to J/m^2
OrganicProperties(:,6) = OrganicProperties(:,6)*1e-3_wp

end subroutine OrganicReader

!==========================================================================================================================

end module InputReader