!***************************************************************************
! drift_correction.f90
! --------------------
! Copyright (C) 2011-2019, LI-COR Biosciences, Inc.  All Rights Reserved.
! Author: Gerardo Fratini
!
! This file is part of EddyPro®.
!
! NON-COMMERCIAL RESEARCH PURPOSES ONLY - EDDYPRO® is licensed for 
! non-commercial academic and government research purposes only, 
! as provided in the EDDYPRO® End User License Agreement. 
! EDDYPRO® may only be used as provided in the End User License Agreement
! and may not be used or accessed for any commercial purposes.
! You may view a copy of the End User License Agreement in the file
! EULA_NON_COMMERCIAL.rtf.
!
! Commercial companies that are LI-COR flux system customers 
! are encouraged to contact LI-COR directly for our commercial 
! EDDYPRO® End User License Agreement.
!
! EDDYPRO® contains Open Source Components (as defined in the 
! End User License Agreement). The licenses and/or notices for the 
! Open Source Components can be found in the file LIBRARIES-ENGINE.txt.
!
! EddyPro® is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!***************************************************************************
!
! \brief       Corrects concentration biases due to instrumental drifts and
!              based on calibration-check data
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine DriftCorrection(Set, nrow, ncol, locCol, ncol2, nCalibEvents, InitialTimestamp)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    integer, intent(in) :: ncol
    integer, intent(in) :: ncol2
    integer, intent(in) :: nCalibEvents
    type(ColType), intent(in) :: locCol(ncol2)
    type(DateType), intent(in) :: InitialTimestamp
    real(kind = dbl), intent(inout) :: Set(nrow, ncol)
    !> Local variables
    integer :: i
    integer :: nPeriods
    integer, external :: NumOfPeriods
    real(kind = dbl) :: lDrift(GHGNumVar)
    real(kind = dbl) :: MeanAbs(ncol)
    real(kind = dbl) :: TempFact

    real(kind = dbl) :: tmp
    real(kind = dbl) :: tmp2


    !> Calculate best guesses of cell and air temperature and pressure
    call AirAndCellParameters()

    !> Calculate drift for the current period, depending on the chosen method
    lDrift(co2:h2o) = error
    select case (trim(adjustl(DriftCorr%method)))
        case ('linear')
            do i = 1, nCalibEvents
                if (Calib(i-1)%ts <= InitialTimestamp &
                    .and. InitialTimestamp <= Calib(i)%ts) then
                    nPeriods = &
                        NumOfPeriods(Calib(i-1)%ts, InitialTimestamp, DateStep)
                    lDrift(co2:h2o) = Calib(i-1)%offset(co2:h2o) &
                        + (Calib(i)%offset(co2:h2o) - Calib(i-1)%offset(co2:h2o)) &
                        / Calib(i)%numPeriods * nPeriods
                    exit
                end if
            end do
            
        case ('signal_strength')
            !> In case of signal strength proxy, calculate
            !> drift based on signal strength
            !> Temperature dependency (LI-7200 manual REv5, Eq. 3-32)

            if (Ambient%Tcell > 0d0) then
                TempFact = 0.6d0 + 0.4d0 / (1d0 &
                    + DriftCorr%b * dexp(DriftCorr%c * (Ambient%Tcell - 273.15d0)))
            else
                TempFact = 1d0
            end if

            !> Detect relevant drift data (all data for periods between
            !> t1 and t2 are stored in Calib(t1))
            do i = 1, nCalibEvents
                if (Calib(i-1)%ts <= InitialTimestamp &
                    .and. InitialTimestamp <= Calib(i)%ts) then
                    where (refCounts(co2:h2o) /= error &
                        .and. Calib(i)%ri(co2:h2o) /= error &
                        .and. Calib(i)%rf(co2:h2o) /= error &
                        .and. Calib(i)%offset(co2:h2o) /= error)
                        lDrift(co2:h2o) = &
                            (refCounts(co2:h2o) - Calib(i)%ri(co2:h2o) * TempFact) &
                            / (Calib(i)%rf(co2:h2o) - Calib(i)%ri(co2:h2o) * TempFact) &
                            * Calib(i)%offset(co2:h2o)
                    elsewhere
                        lDrift(co2:h2o) = error
                    end where
                    exit
                end if
            end do
    end select


!> Only for debug, eliminate!
tmp  = Set(1, co2)
tmp2 = Set(1, h2o)

    !> This call only to calculate chi_h2o, needed for equivalent pressure
    call MoleFractionsAndMixingRatios()

    !> If chi(h2o) could not be calculated, set it to zero, which in this 
    !> context means not accounting for broadening effects
    if (Stats%chi(h2o) == error) Stats%chi(h2o) = 0d0

    !> Convert to density/press (note the 10^3 to get to Abs/P with P in kPa)
    !> co2
    if (locCol(co2)%measure_type /= 'molar_density') then
        where (Set(:, co2) /= error)
            !> last parenthesis here is for P_ec = P*[1 + 0.15 chi_h2o]
            Set(:, co2) = Set(:, co2) / Ru / Ambient%Tcell &
                / (1d0 + 0.15d0 * Stats%chi(h2o) * 1d-3)
        end where
    else
        where (Set(:, co2) /= error)
            !> last parenthesis here is for P_ec = P*[1 + 0.15 chi_h2o]
            Set(:, co2) = Set(:, co2) &
                / (Ambient%Pcell / 1d3 * (1d0 + 0.15d0 * Stats%chi(h2o) * 1d-3))
        end where
    end if

    !> h2o
    if (locCol(h2o)%measure_type /= 'molar_density') then
        where (Set(:, h2o) /= error)
            Set(:, h2o) = Set(:, h2o) / Ru / Ambient%Tcell * 1d3
        end where
    else
        where (Set(:, h2o) /= error)
            Set(:, h2o) = Set(:, h2o) / (Ambient%Pcell / 1d3)
        end where
    end if

    !> Convert densities/press to absorbances/press
    call PolyVal(DriftCorr%inv_cal(0:6, co2), 6, Set(:, co2), size(Set, 1), Set(:, co2))
    call PolyVal(DriftCorr%inv_cal(0:6, h2o), 6, Set(:, h2o), size(Set, 1), Set(:, h2o))

    !> Calculate mean absorptances
    call AverageNoError(Set, size(Set, 1), size(Set, 2), MeanAbs, error)

    !> Remove absorptance/press drift if detected for current period
    !> Note: it can be demonstrated (see Fratini et al. 2014, BG, Eqs. 10-11), that:
    !> abs_theor = (mean_abs_meas - bias) + d_abs_meas / (1 - bias * P [kPa])
    !> where all terms are intended as normalized by pressure.
    if (lDrift(co2) /= error) then
        where (Set(:, co2) /= error .and. MeanAbs(co2) /= error)
            Set(:, co2) = (MeanAbs(co2) - lDrift(co2)) + &
                (Set(:, co2) - MeanAbs(co2)) &
                / (1d0 - lDrift(co2) * Ambient%Pcell * 1d-3)
        end where
    end if
    if (lDrift(h2o) /= error) then
        where (Set(:, h2o) /= error .and. MeanAbs(h2o) /= error)
            Set(:, h2o) = (MeanAbs(h2o) - lDrift(h2o)) + &
                (Set(:, h2o) - MeanAbs(h2o)) &
                / (1d0 - lDrift(h2o) * Ambient%Pcell * 1d-3)
        end where
    end if

    !> Convert absorptances/press back to density/press
    call PolyVal(DriftCorr%dir_cal(0:6, co2), 6, Set(:, co2), size(Set, 1), Set(:, co2))
    call PolyVal(DriftCorr%dir_cal(0:6, h2o), 6, Set(:, h2o), size(Set, 1), Set(:, h2o))

    !> Convert density/press back to concentration or density
    !> co2
    if (locCol(co2)%measure_type /= 'molar_density') then
        where (Set(:, co2) /= error)
            Set(:, co2) = Set(:, co2) * Ru * Ambient%Tcell &
                * (1d0 + 0.15d0 * Stats%chi(h2o) * 1d-3)
        end where
    else
        where (Set(:, co2) /= error)
            Set(:, co2) = Set(:, co2) * &
                (Ambient%Pcell / 1d3 * (1d0 + 0.15d0 * Stats%chi(h2o) * 1d-3))
        end where
    end if
    !> h2o
    if (locCol(h2o)%measure_type /= 'molar_density') then
        where (Set(:, h2o) /= error)
            Set(:, h2o) = Set(:, h2o) * Ru * Ambient%Tcell / 1d3
        end where
    else
        where (Set(:, h2o) /= error)
            Set(:, h2o) = Set(:, h2o) * (Ambient%Pcell / 1d3)
        end where
    end if

!> ONLY FOR DEBUG, ELIMINATE!
write(987,*) InitialTimestamp, &
    (tmp - Set(1, co2)) * Ru * Ambient%Tcell / Ambient%Pcell * 1d3

end subroutine DriftCorrection

!***************************************************************************
!
! \brief       Calculate mean reference counts from raw data if available
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ReferenceCounts(Set, nrow, ncol)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    integer, intent(in) :: ncol
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    !> local variables
    integer :: i
    real(kind = dbl) :: meanSet(ncol)


    !> Calculate mean values from raw data
    call AverageNoError(Set, nrow,  ncol, meanSet, error)

    !> Extract mean counts from average, if available
    refCounts = error
    do i = 1, NumCol
        select case (Col(i)%var)
            case('co2_ref')
                refCounts(co2) = meanSet(i)
            case('h2o_ref')
                refCounts(h2o) = meanSet(i)
        end select
    end do
end subroutine ReferenceCounts
