module DistFuncs

use dvode_kinds_module, only: wp => dvode_wp
use EnvironmentConstants, only: pi

contains

!==========================================================================================================================

!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Subroutine to create the randomly sampled normal size distribution.                *
!*                                                                                      *
!*                                                                                      *
!*   -> created:        2025/11/04                                                      *
!*                                                                                      *
!****************************************************************************************

subroutine NormDist(mu, sigma, nn, randnum)

implicit none

real(wp), intent(in)                    :: mu, sigma    ! Mean and standard deviation of the normal distribution in um
integer, intent(in)                     :: nn           ! Number of particles
integer                                 :: i
real(wp)                                :: u1, u2, z1
real(wp), dimension(nn), intent(out)    :: randnum
    
do i = 1, nn
    ! Generate two independent uniform random numbers between 0 and 1
    call RANDOM_NUMBER(u1)
    call RANDOM_NUMBER(u2)

    ! Apply the Box-Muller transform to get a standard normal variate (z1)
    z1 = SQRT(-2.0_wp * LOG(u1)) * COS(2.0_wp * pi * u2)

    ! Scale and shift the standard normal variate to the desired mean and standard deviation
    randnum(i) = mu + sigma * z1
end do

end subroutine NormDist

!==========================================================================================================================

!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Subroutine to create the randomly sampled lognormal size distribution.             *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Dan Barthaux                                                                       *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2025/11/04                                                      *
!*                                                                                      *
!****************************************************************************************

subroutine LogNormDist(mu, sigma, nn, randnum)

implicit none

real(wp), intent(in)                    :: mu, sigma    ! Mean and standard deviation of the normal distribution in um
integer, intent(in)                     :: nn           ! Number of particles
real(wp), dimension(nn), intent(out)    :: randnum
real(wp), dimension(nn)                 :: normvals

! Generate normal random variates with mean=log(mu), std=log(sigma)
call NormDist(log(mu), log(sigma), nn, normvals)

! Exponentiate to get log-normal
randnum = exp(normvals)

end subroutine LogNormDist

!==========================================================================================================================

end module DistFuncs