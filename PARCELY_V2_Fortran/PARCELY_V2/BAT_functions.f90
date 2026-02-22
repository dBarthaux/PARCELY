module BATModel

use dvode_kinds_module, only: wp => dvode_wp
use ModelParameters, only: norg

INTEGER(4), PARAMETER, PUBLIC   :: clen = 45
REAL(wp), PARAMETER, PUBLIC     :: M_water = 18.015_wp

contains

!****************************************************************************************
!*   SUBROUTINE BAT_activity_calc_with_refinement_v1                                    *
!*                                                                                      *
!*   :: Purpose ::                                                                      *
!*   Calculates activity organics and water.                                            *
!*   Calculates the water and organic mass fraction associated with every organic       *
!*   species (every binary mixture).                                                    *
!*   This is performed according to the Duhem-Margules relation, which relates the      *
!*   the molar excess Gibbs energy of mixing to the two mole-fraction-based activity    *
!*   coefficients (of water and organics). See equations 13 and 14 of BAT paper.        *
!*                                                                                      *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *
!*   McGill University, Montreal, Quebec (2021),                                        *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created: 15-04-2021                                                             *
!*   -> latest changes:                                                                 *
!****************************************************************************************

!SUBROUTINE BAT_properties_calculation_v1(org_mole_fraction, O2C, H2C, molarmass_ratio, BAT_functional_group, &
!    & N2C_values_densityOnly, ln_func1, ln_func2, ycalc1, ycalc2, activity_calc1, activity_calc2, mass_fraction1, &
!    & mass_fraction2, Gibbs_RT, dGibbs_RTdx2)

SUBROUTINE BAT_properties_calculation_v1(org_mole_fraction, O2C, H2C, molarmass_ratio, BAT_functional_group, N2C_values_densityOnly, ycalc2)

IMPLICIT NONE

! Inputs
REAL(wp), DIMENSION(norg), INTENT(IN)           :: org_mole_fraction
REAL(wp), DIMENSION(norg), INTENT(IN)           :: O2C
REAL(wp), DIMENSION(norg), INTENT(IN)           :: H2C
REAL(wp), DIMENSION(norg), INTENT(IN)           :: N2C_values_densityOnly
REAL(wp), DIMENSION(norg), INTENT(IN)           :: molarmass_ratio
CHARACTER(4), DIMENSION(norg), INTENT(IN)       :: BAT_functional_group

! Outputs
REAL(wp), DIMENSION(norg), INTENT(OUT)          :: ycalc2

REAL(wp), DIMENSION(norg)                       :: ln_func1, ln_func2, ycalc1

!REAL(wp), DIMENSION(:), INTENT(OUT)      :: ln_func1, & ! ln(act. coeff. water)
!                                                & ln_func2, & ! ln(act. coeff. org)
!                                                & ycalc1, & ! Mole-fraction based activity coefficient of water
!                                                & ycalc2, & ! Mole-fraction based activity coefficient of organic species
!                                                & activity_calc1, & ! Activity of water
!                                                & activity_calc2, & ! Activity of organic species
!                                                & mass_fraction1, & ! Mass fraction of water
!                                                & mass_fraction2, & ! Mass fraction of organic species
!                                                & Gibbs_RT, & ! Molar excess Gibbs energy of mixing G_excess/RT
!                                                & dGibbs_RTdx2 ! Derivative of molar excess Gibbs energy of mixing G_excess/RT with respect to the mole
                                                            ! the mole fraction of organic species

! Local variables
INTEGER                                            :: allocstat, do_calc, loop_i, n, n1, j
REAL(wp)                                           :: tran_lowO2C_fractionOne_phase, tran_lowO2C_sigmoid_bend, tran_lowO2C_sigmoid_shift, &
                                                        tran_midO2C_sigmoid_bend, tran_midO2C_sigmoid_shift
REAL(wp), DIMENSION(SIZE(org_mole_fraction, 1))    :: x2, x2_temp, mole_fraction_mask, not_mole_fraction_mask, &
                                                        phi2, dphi2dx2, sum1, sum2, dgemix2dx_temp, gemix_temp
REAL(wp), DIMENSION(SIZE(O2C, 1))                  :: O2C_temp, molarmass_ratio_temp, Mr_massfrac_final, Mr, densityEst, Onephase_O2C, &
                                                        mid_transition, O2C_1phase_delta, weight_1phase, weight_2phase, &
                                                        O2C_1phase_delta_norm, weight_1phase_norm, rhor, scaledMr
REAL(wp), DIMENSION(SIZE(org_mole_fraction, 1), 2) :: gemix, dgemix2dx
CHARACTER(clen)                                    :: fitpar_name1, fitpar_name2
REAL(wp), DIMENSION(10)                            :: fitpar_lowO2C, fitpar_midO2C, fitpar_highO2C, fitpar_1phase, fitpar_2phase, fitpar
REAL(wp), DIMENSION(2)                             :: coeff

tran_lowO2C_fractionOne_phase = 0.189974476118418_wp
tran_lowO2C_sigmoid_bend = 79.2606902175984_wp
tran_lowO2C_sigmoid_shift = 0.0604293454322489_wp
tran_midO2C_sigmoid_bend = 75.0159268221068_wp
tran_midO2C_sigmoid_shift = 0.000947111285750515_wp

! Shift for equivalent O/C and MW, for activity calculation
Mr_massfrac_final = molarmass_ratio ! Used in final mass fraction calc and not for activity coefficent calc.
!PRINT *, O2C, molarmass_ratio, BAT_functional_group
CALL convert_chemical_structure_to_OH_eqv_v3(O2C, molarmass_ratio, BAT_functional_group, O2C_temp, molarmass_ratio_temp)

! Force org mole fraction to be 1 or less
x2 = org_mole_fraction
WHERE (org_mole_fraction > 1.0_wp) x2 = 1.0_wp
WHERE (org_mole_fraction < 0.0_wp) x2 = 0.0_wp
        
! Replace infinite dilution of zero to small value
WHERE (org_mole_fraction < 1.0E-20_wp) x2 = 1.0E-20_wp

Mr = molarmass_ratio_temp  ! molar1/molar2

!  Properties
!  Estimate density of org. components using translated O/C, and molecular weight ratio (to OH-equiv.) and model by Girolami (1994)
CALL Org_density_Estimate_KGv1((M_water/Mr), O2C_temp, H2C, N2C_values_densityOnly, densityEst)
! Estimate O/C limit of miscibility point
CALL single_phase_O2C_point_KGv3(Mr, Onephase_O2C)

! Get region transition properties (BAT model fit parameters for transition from one O/C region to another O/C region).
mid_transition(1) = Onephase_O2C(1) * 0.75_wp
IF (O2C_temp(1) < mid_transition(1)) THEN !  lower to mid O/C region

    ! Data point trasfer weight
    O2C_1phase_delta(1) = O2C_temp(1) - (Onephase_O2C(1) * tran_lowO2C_fractionOne_phase)
    weight_1phase(1) = 1.0_wp/(1.0_wp + EXP(-tran_lowO2C_sigmoid_bend * &
    & (O2C_1phase_delta(1) - tran_lowO2C_sigmoid_shift))) ! logistic transfer function 1/(1+e^-(75*x))

    ! Normalize to end point so at mid_transition weight 2 is 1.
    O2C_1phase_delta_norm(1) = O2C_temp(1) - (mid_transition(1) * tran_lowO2C_fractionOne_phase)
    weight_1phase_norm(1) = 1.0_wp/(1.0_wp + EXP(-tran_lowO2C_sigmoid_bend * &
        (O2C_1phase_delta_norm(1) - tran_lowO2C_sigmoid_shift))) ! Logistic transfer function 1/(1+e^-(75*x))

    weight_1phase(1) = weight_1phase(1)/weight_1phase_norm(1)
    weight_2phase(1) = 1.0_wp - weight_1phase(1)

    fitpar_name1 = 'midO2C'
    fitpar_name2 = 'lowO2C'

ELSE IF (O2C_temp(1) < Onephase_O2C(1) * 2) THEN ! mid to high O/C region

    O2C_1phase_delta(1) = O2C_temp(1) - Onephase_O2C(1)
    weight_1phase(1) = 1.0_wp/(1.0_wp + EXP(-tran_midO2C_sigmoid_bend * &
        (O2C_1phase_delta(1) - tran_midO2C_sigmoid_shift)))
    weight_2phase(1) = 1.0_wp - weight_1phase(1)

    fitpar_name1 = 'highO2C'
    fitpar_name2 = 'midO2C'

ELSE ! High only region

    fitpar_name1 = 'highO2C'
    fitpar_name2 = 'NaN'

    weight_2phase(1) = 0.0_wp
    weight_1phase(1) = 1.0_wp
        
END IF

! Fit properties data
fitpar_lowO2C = [7.089476E+00_wp, -0.6226781_wp, -7.711860E+00_wp, -1.000000E+02_wp, -3.885941E+01_wp, &
    3.081244E-09_wp, -1.000000E+02_wp, 6.188812E+01_wp, -5.988895E+00_wp, 6.940689E+00_wp]

fitpar_midO2C = [5.872214E+00_wp, -0.9740486_wp, -4.535007E+00_wp, -1.000000E+02_wp, -5.129327E+00_wp, &
    2.109751E+00_wp, -2.809232E+01_wp, -2.367683E+01_wp, -1.219164E+00_wp, 4.742729E+00_wp]

fitpar_highO2C = [5.921550E+00_wp, -1.000000E+02_wp, -2.528295E+00_wp, -1.000000E+02_wp, -3.883017E+00_wp, &
    1.353916E+00_wp, -7.898128E+00_wp, -1.160145E+01_wp, -0.07868187_wp, 3.650860E+00_wp]

! Get fit parameters in into correct phases: phase 1
IF (fitpar_name1 == 'highO2C') THEN
    fitpar_1phase = fitpar_highO2C
ELSE IF (fitpar_name1 == 'midO2C') THEN
    fitpar_1phase = fitpar_midO2C
ELSE IF (fitpar_name1 == 'lowO2C') THEN
    fitpar_1phase = fitpar_lowO2C
END IF

! Get fit parameters in into correct phases: phase 2
IF (fitpar_name2 == 'highO2C') THEN
    fitpar_2phase = fitpar_highO2C
ELSE IF (fitpar_name2 == 'midO2C') THEN
    fitpar_2phase = fitpar_midO2C
ELSE IF (fitpar_name2 == 'lowO2C') THEN
    fitpar_2phase = fitpar_lowO2C
END IF

!! For biphasic line calculations
!IF (BAT_functional_group(1) == 'initial fitting for biphasic transition') THEN
!
!    ! Inital fit to biphasic line
!    fitpar_1phase = [5.885109E+00_wp, -9.849013E-01_wp, -4.731250E+00_wp, -6.227207E+00_wp, -5.201652E+00_wp, &
!        2.320286E+00_wp, -3.082297E+01_wp, -2.584037E+01_wp, -1.237227E+00_wp, 4.069905E+00_wp]
!    weight_2phase(1) = 0.0_wp
!    weight_1phase(1) = 1.0_wp
!
!END IF

gemix = 0.0_wp ! Molar excess Gibbs energy of mixing
dgemix2dx = 0.0_wp ! Derivative of the molar excess Gibbs energy of mixing with respect to the mole fraction of organic species
do_calc = 0

DO loop_i = 1, 2
    IF (loop_i == 1) THEN
        IF (weight_1phase(1) < EPSILON(1.0_wp)) THEN
            do_calc = 0
        ELSE
            fitpar = fitpar_1phase
            do_calc = 1
        END IF
    ELSE
        IF (weight_2phase(1) < EPSILON(1.0_wp)) THEN
            do_calc = 0
        ELSE
            fitpar = fitpar_2phase
            do_calc = 1
        END IF
	END IF

    IF (do_calc == 1) THEN
        n = SIZE(fitpar, 1)
        rhor(1) = 0.997_wp/densityEst(1) ! Assumes water

        ! scaledMr is the scaled molar mass ratio of this mixture's components
        scaledMr(1) = Mr(1) * fitpar(n) * (1.0_wp + O2C_temp(1))**fitpar(n-1)

        ! phi2 is a scaled volume fraction, which is determined from mole fractions, scaledMr
        ! and an estimate of water/organic density ratio (rhor).
				
        phi2 = x2/(x2 + (1.0_wp - x2) * scaledMr(1)/rhor(1)) ! and phi1 = 1 - phi2
        dphi2dx2 = (scaledMr(1)/rhor(1)) * (1.0_wp/(x2 + (1.0_wp - x2) * scaledMr(1)/rhor(1)))**2 ! the derivative of phi2 with respect to x2 (mole fraction)

        n1 = (n-2)/4
        coeff(1:n1) = fitpar(1:n1) * EXP(fitpar(n1+1:2*n1)*O2C_temp(1)) + fitpar(2*n1+1:3*n1) * &
            EXP(fitpar(3*n1+1:4*n1)*Mr(1))

        sum1 = 0.0_wp
        DO j = 1, n1
            sum1 = sum1 + coeff(j) * (1.0_wp - 2.0_wp * phi2)**(j-1)
        END DO

        sum2 = 0.0_wp
        DO j = 2, n1 ! first derivative
            sum2 = sum2 + 2.0_wp * REAL(j-1) * coeff(j) * (1.0_wp - 2.0_wp*phi2)**(j-2)
        END DO
            
        gemix(:,loop_i) = phi2 * (1.0_wp - phi2) * sum1
        dgemix2dx(:,loop_i) = ((1.0_wp - 2.0_wp * phi2) * sum1 + phi2 * (1.0_wp - phi2) * (-sum2)) * dphi2dx2

        END IF
END DO

gemix_temp = gemix(:,1) * weight_1phase(1) + gemix(:,2) * weight_2phase(1)
dgemix2dx_temp = dgemix2dx(:,1) * weight_1phase(1) + dgemix2dx(:,2) * weight_2phase(1)

! Calculate the function value funcx1 (= y1(x2)) at point with w2:
ln_func1 = gemix_temp - x2 * dgemix2dx_temp ! The func value for component 1 = LOG(activity coeff. water)
ln_func2 = gemix_temp + (1.0_wp - x2) * dgemix2dx_temp ! The func value of the component 2 = LOG(activity coeff. of the organic)

WHERE (ln_func1 < -690.7755_wp) ln_func1 = -690.7755_wp
WHERE (ln_func1 > 690.7755_wp) ln_func1 = 690.7755_wp
ycalc1 = EXP(ln_func1) ! activity coefficient water

WHERE (ln_func2 < -690.7755_wp) ln_func2 = -690.7755_wp
WHERE (ln_func2 > 690.7755_wp) ln_func2 = 690.7755_wp
ycalc2 = EXP(ln_func2) ! activity coefficient org

!activity_calc1 = ycalc1 * (1.0_wp - x2) ! activity water
!activity_calc2 = ycalc2 * x2 ! avtivity org

!mass_fraction1 = (1.0_wp - x2) * Mr_massfrac_final(1)/((1.0_wp - x2) * (Mr_massfrac_final(1) - 1.0_wp) + 1.0_wp) ! mass fraction water
!mass_fraction2 = 1.0_wp - mass_fraction1 ! mass fraction org
!Gibbs_RT = gemix_temp ! G_excess/RT
!dGibbs_RTdx2 = dgemix2dx_temp ! d(G_excess/RT)/d (mole fraction org)

END SUBROUTINE BAT_properties_calculation_v1

!-----------------------------------------------------------------------------------------------------------------------------------------------

!****************************************************************************************
!*   SUBROUTINE convert_chemical_structure_to_OH_eqv_v3                                 *
!*                                                                                      *
!*   :: Purpose ::                                                                      *
!*   Translate the properties of functionalized molecules to a hypothetical             *
!*   OH-equivalent molecule of modified O/C and molecular weight using shift fit data.  *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *
!*   McGill University, Montreal, Quebec (2021),                                        *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created: 15-04-2021                                                             *
!*   -> latest changes:                                                                 *
!****************************************************************************************

SUBROUTINE convert_chemical_structure_to_OH_eqv_v3(O2C, molarmass_ratio, shift_method, O2C_eqv, molarmass_ratio_eqv)

IMPLICIT NONE

! Input
REAL(wp), DIMENSION(:), INTENT(IN)       :: O2C, molarmass_ratio ! Elemental O/C ratio and molecular weight of each org. species
CHARACTER(4), DIMENSION(:), INTENT(IN)   :: shift_method ! functional group of each organic species

! Output
REAL(wp), DIMENSION(:), INTENT(OUT)      :: O2C_eqv, molarmass_ratio_eqv ! Translated molecule properties (to OH-equivalent):
                                                                            ! Hypothetical OH-equivalent O/C and molecular weight ratios

! Local variables
REAL(wp), DIMENSION(4)                   :: fit_shift
REAL(wp)                                 :: MW_old, MW_new
INTEGER                                  :: sp_i, calc_shift

DO sp_i = 1, SIZE(O2C, 1)
    calc_shift = 1

    IF  ((shift_method(sp_i) == 'hydroxyl') .OR. (shift_method(sp_i) == 'carboxyl')) THEN
        calc_shift = 0
    ELSE IF ((shift_method(sp_i) == 'hydroperoxideSOA') .OR. (shift_method(sp_i) == 'SOA chemicals')) THEN
        fit_shift = [0.000149021455554774_wp, 0.00473627706154738_wp, 0.869057801919811_wp, 0.564783492434737_wp]
    ELSE IF (shift_method(sp_i) == 'PEG') THEN
        fit_shift = [0.00544768078879267_wp, 3.86433605822883_wp, -0.267168022244528_wp, 0.255486696379870_wp]
    ELSE IF (shift_method(sp_i) == 'ketone') THEN
        fit_shift = [0.00453425820008037_wp, 0.000648450309348991_wp, 0.138144408029760_wp, 0.352454367906330_wp]
    ELSE IF (shift_method(sp_i) == 'hydroperoxide') THEN
        fit_shift = [8.17160348517250E-06_wp, 4.53182593743281E-07_wp, 0.966089559236154_wp, 0.459433193460024_wp]
    ELSE IF (shift_method(sp_i) == 'ether') THEN
        fit_shift = [2.44336043644347E-05_wp, 0.000158316167479487_wp, 0.284974095922167_wp, 0.229338647030993_wp]
    ELSE IF (shift_method(sp_i) == 'ester') THEN
        fit_shift = [-1.29324633442629_wp, 0.00108128314380665_wp, 1.24051435479678_wp, 0.405354156019572_wp]
    ELSE
        !PRINT *, 'no O2C and MW system selected'
        calc_shift = 0
    END IF

    IF (calc_shift == 1) THEN
        O2C_eqv(sp_i) = O2C(sp_i)/(1.0_wp + fit_shift(3)*EXP(-(O2C(sp_i))*fit_shift(1)))
        MW_old = (18.016_wp/molarmass_ratio(sp_i))
        MW_new = MW_old/(1.0_wp + fit_shift(4)*EXP(-(MW_old)*fit_shift(2)))
        molarmass_ratio_eqv(sp_i) = 18.016_wp/MW_new
    ELSE
        O2C_eqv(sp_i) = O2C(sp_i)
        molarmass_ratio_eqv(sp_i) = molarmass_ratio(sp_i)
    END IF
END DO

END SUBROUTINE convert_chemical_structure_to_OH_eqv_v3

!-----------------------------------------------------------------------------------------------------------------------------------------------

!****************************************************************************************
!*   SUBROUTINE Org_density_Estimate_KGv1                                               *
!*                                                                                      *
!*   :: Purpose ::                                                                      *
!*   Subroutine that computes the density of pure liquid organic                        *
!*   compounds at 298.15 K using a simple structure-based method                        *
!*   (G. S. Girolami J. Chem. Educ. 71, 962-964, 1994).                                 *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *
!*   McGill University, Montreal, Quebec (2021),                                        *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created: 15-04-2021                                                             *
!*   -> latest changes:                                                                 *
!****************************************************************************************
 
SUBROUTINE Org_density_Estimate_KGv1(M, O2C, H2C, N2C, densityEst)

IMPLICIT NONE
        
! Inputs
REAL(wp), DIMENSION(:), INTENT(IN)  :: M, O2C, H2C, N2C ! elemental ratios (O/C, H/C, N/C)
                                                            ! and molecular weight of organic species

! Outputs
REAL(wp), DIMENSION(:), INTENT(OUT) :: densityEst       ! density estimate for each organic species

! Local variables
REAL(wp), DIMENSION(SIZE(M ,1))     :: H2Cest, NC, rho1, N2Cest
REAL(wp) :: MC, MO, MH, MN
INTEGER :: i, sizeM

sizeM = SIZE(M,1)

!  DensityEst(M, O2C, H2C) ,  H:C input is optional (when known, otherwise enter a negative value)
! the molar masses of the carbon, oxygen and hydrogen atoms in [g/mol]
MC = 12.0100_wp
MO = 16.0000_wp
MH = 1.00800_wp
MN = 14.0067_wp
        
N2Cest = N2C
H2Cest = H2C
        
! 1) Set the N2C value if not provided from input (if a negative value is given at input)
WHERE (N2C < 0.0_wp) N2Cest = 0.0_wp
        
! 2) Estimate the H2C value if not provided from input (if a negative value is given at input)
WHERE (H2C < 0.0_wp) H2Cest = 2.0_wp - O2C ! estimate H2C assuming an aliphatic compound with H2C = 2 in absence of oxygen-bearing.
                                                    ! functional groups, then correct for oxygen content assuming a-1 slope (Van Krevelen Diagram of typical SOA).

! 3) Compute the approx. number of carbon atoms per organic molecule
NC = M/(MC + H2Cest*MH + O2C*MO + N2Cest*MN)

! 4) compute density estimate based on method by Girolami (1994)
!  here no correction is applied for rings and aromatic compounds (due to limited info at input)
rho1 = M/(5.0_wp*NC*(2.0_wp + H2Cest + O2C*2.0_wp + N2Cest*2.0_wp))

DO i = 1, SIZE(NC, 1)
    densityEst(i) = rho1(i)*(1.000_wp + MINVAL([NC(i) * O2C(i) * 0.100_wp + NC(i) * N2Cest(i) * &
        0.100_wp, 0.300_wp], 1))
! density in [g/cm^3]; here it is scaled assuming that most of the oxygen atoms are able
! to make H-bonds (donor or acceptor).
END DO

END SUBROUTINE Org_density_Estimate_KGv1

!-----------------------------------------------------------------------------------------------------------------------------------------------

!****************************************************************************************
!*   SUBROUTINE single_phase_O2C_point_KGv3                                             *
!*                                                                                      *
!*   :: Purpose ::                                                                      *
!*   Sigmoid function that represents the miscibility limit in O/C vs molar mass space. *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Camilo Serrano, Andreas Zuend, Kyle Gorkowski                                      *
!*   McGill University, Montreal, Quebec (2021),                                        *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created: 15-04-2021                                                             *
!*   -> latest changes:                                                                 *
!****************************************************************************************

SUBROUTINE single_phase_O2C_point_KGv3(molarmass_ratio, O2C_single_phase_cross_point)

IMPLICIT NONE
        
! Inputs
REAL(wp), DIMENSION(:), INTENT(IN)  :: molarmass_ratio
! Outputs
REAL(wp), DIMENSION(:), INTENT(OUT) :: O2C_single_phase_cross_point

O2C_single_phase_cross_point = 0.205_wp/(1.0_wp + EXP(26.6_wp * (molarmass_ratio - 0.12_wp)))**0.843_wp + &
    & 0.225_wp

END SUBROUTINE single_phase_O2C_point_KGv3

!==========================================================================================================================

end module BATModel