!***************************************************************************
! drift_correction.f90
! --------------------
! Copyright (C) 2013-2014, LI-COR Biosciences
!
! This file is part of EddyPro (TM).
!
! EddyPro (TM) is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! EddyPro (TM) is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with EddyPro (TM).  If not, see <http://www.gnu.org/licenses/>.
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
    integer :: nTotPeriods
    integer, external :: NumOfPeriods
    real(kind = dbl) :: lDrift(GHGNumVar)
    real(kind = dbl), allocatable :: MeanAbs(:)
    real(kind = dbl), allocatable :: TempFact

    real(kind = dbl) :: tmp
    real(kind = dbl) :: tmp2


    !> To be refined: find best estimate of cell temperature
    if (Stats%Mean(tc) <= 0d0) Stats%Mean(tc) = Stats%Mean(ts)


    !> Calculate drift for the current period, depending on the method chosen
    lDrift(co2:h2o) = error
    select case (trim(adjustl(DriftCorr%method)))
        case ('linear')
            do i = 1, nCalibEvents
                if (Calib(i-1)%ts <= InitialTimestamp .and. InitialTimestamp <= Calib(i)%ts) then
                    nTotPeriods = NumOfPeriods(Calib(i-1)%ts, Calib(i)%ts, DateStep)
                    nPeriods    = NumOfPeriods(Calib(i-1)%ts, InitialTimestamp, DateStep)
                    lDrift(co2:h2o) = Calib(i)%offset(co2:h2o) &
                        + (Calib(i)%offset(co2:h2o) - Calib(i-1)%offset(co2:h2o)) / nTotPeriods * nPeriods
                    exit
                end if
            end do

        case ('signal_strength')
            !> In case of signal strength proxy, calculate drift based on signal strength
            !> Temperature dependency (LI-7200 manual REv5, Eq. 3-32)
            if (Stats%mean(tc) > 0d0) then
                TempFact = 0.6d0 + 0.4d0 / (1d0 + DriftCorr%b * dexp(DriftCorr%c * (Stats%Mean(tc) - 273.16d0)))
            else
                TempFact = 1d0
            end if

            !> Detect relevant drift data (all data for periods between t1 and t2 are stored in Calib(t1))
            do i = 1, nCalibEvents
                if (Calib(i-1)%ts <= InitialTimestamp .and. InitialTimestamp <= Calib(i)%ts) then
                    where (refCounts(co2:h2o) /= error .and. Calib(i)%ri(co2:h2o) /= error &
                        .and. Calib(i)%rf(co2:h2o) /= error .and. Calib(i)%offset(co2:h2o) /= error)
                        lDrift(co2:h2o) = (refCounts(co2:h2o) - Calib(i)%ri(co2:h2o) * TempFact) &
                            / (Calib(i)%rf(co2:h2o) - Calib(i)%ri(co2:h2o) * TempFact) * Calib(i)%offset(co2:h2o)
                    elsewhere
                        lDrift(co2:h2o) = error
                    end where
                    exit
                end if
            end do
    end select


    tmp  = Set(1, co2)
    tmp2 = Set(1, h2o)

    !> This call only to calculate chi_h2o, needed for equivalent pressure
    call MoleFractionsAndMixingRatios()

    !> Convert to density/press (note the 10^3 to get to Abs/P with P in kPa)
    !> co2
    if (locCol(co2)%measure_type /= 'molar density') then
        where (Set(:, co2) /= error)
            Set(:, co2) = Set(:, co2) / Ru / Stats%Mean(tc) / (1d0 + 0.15d0 * Stats%chi(h2o) * 1d-3)    !> last parenthesis for P_ec = P*[1 + 0.15 chi_h2o]
        end where
    else
        where (Set(:, co2) /= error)
            Set(:, co2) = Set(:, co2) / (Stats%Mean(pc) / 1d3 * (1d0 + 0.15d0 * Stats%chi(h2o) * 1d-3)) !> last parenthesis for P_ec = P*[1 + 0.15 chi_h2o]
        end where
    end if
    !> h2o
    if (locCol(h2o)%measure_type /= 'molar density') then
        where (Set(:, h2o) /= error)
            Set(:, h2o) = Set(:, h2o) / Ru / Stats%Mean(tc) * 1d3
        end where
    else
        where (Set(:, h2o) /= error)
            Set(:, h2o) = Set(:, h2o) / (Stats%Mean(pc) / 1d3)
        end where
    end if

    !> Convert densities/press to absorbances/press
    call PolyVal(DriftCorr%inv_cal(0:6, co2), 6, Set(:, co2), size(Set, 1), Set(:, co2))
    call PolyVal(DriftCorr%inv_cal(0:6, h2o), 6, Set(:, h2o), size(Set, 1), Set(:, h2o))

    !> Calculate mean absorptances
    if (.not. allocated (MeanAbs)) allocate(MeanAbs(size(Set, 2)))
    call AverageNoError(Set, size(Set, 1), size(Set, 2), MeanAbs, error)

    !> Remove absorptance/press drift if detected for current period
    !> Note: it can be demonstrated (see Fratini et al. 2013, BGD, Eqs. 10-11), that:
    !> abs_theor = (mean_abs_meas - bias) + d_abs_meas / (1 - bias * P [kPa])
    !> where all terms are intended as normalized by pressure.
    if (lDrift(co2) /= error) then
        where (Set(:, co2) /= error .and. MeanAbs(co2) /= error)
            Set(:, co2) = (MeanAbs(co2) - lDrift(co2)) + &
                (Set(:, co2) - MeanAbs(co2)) / (1d0 - lDrift(co2) * Stats%Mean(pc) * 1d-3)
        end where
    end if
    if (lDrift(h2o) /= error) then
        where (Set(:, h2o) /= error .and. MeanAbs(h2o) /= error)
            Set(:, h2o) = (MeanAbs(h2o) - lDrift(h2o)) + &
                (Set(:, h2o) - MeanAbs(h2o)) / (1d0 - lDrift(h2o) * Stats%Mean(pc) * 1d-3)
        end where
    end if
    if (allocated(MeanAbs)) deallocate(MeanAbs)

    !> Convert absorptances/press back to density/press
    call PolyVal(DriftCorr%dir_cal(0:6, co2), 6, Set(:, co2), size(Set, 1), Set(:, co2))
    call PolyVal(DriftCorr%dir_cal(0:6, h2o), 6, Set(:, h2o), size(Set, 1), Set(:, h2o))

    !> Convert density/press back to concentration or density
    !> co2
    if (locCol(co2)%measure_type /= 'molar density') then
        where (Set(:, co2) /= error)
            Set(:, co2) = Set(:, co2) * Ru * Stats%Mean(tc) * (1d0 + 0.15d0 * Stats%chi(h2o) * 1d-3)
        end where
    else
        where (Set(:, co2) /= error)
            Set(:, co2) = Set(:, co2) * (Stats%Mean(pc) / 1d3 * (1d0 + 0.15d0 * Stats%chi(h2o) * 1d-3))
        end where
    end if
    !> h2o
    if (locCol(h2o)%measure_type /= 'molar density') then
        where (Set(:, h2o) /= error)
            Set(:, h2o) = Set(:, h2o) * Ru * Stats%Mean(tc) / 1d3
        end where
    else
        where (Set(:, h2o) /= error)
            Set(:, h2o) = Set(:, h2o) * (Stats%Mean(pc) / 1d3)
        end where
    end if

    write(124,*) InitialTimestamp, tmp - Set(1, co2), tmp2 - Set(1, h2o)
end subroutine DriftCorrection
