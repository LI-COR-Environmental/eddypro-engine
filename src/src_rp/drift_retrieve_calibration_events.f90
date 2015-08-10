!***************************************************************************
! drift_retrieve_calibration_events.f90
! -------------------------------------
! Copyright (C) 2011-2015, LI-COR Biosciences
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
!*******************************************************************************
!
! \brief       Reads external dynamic (time-varying) metadata file and retrieve
!              cleaning and calibration events information
!              alternative metadata file
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!*******************************************************************************
subroutine driftRetrieveCalibrationEvents(nCalibEvents)
    use m_rp_global_var
    implicit none
    !> In/out variables
    integer, intent(out) :: nCalibEvents
    !> local variables
    integer, parameter :: cal_date  = 1
    integer, parameter :: cal_time  = 2
    integer, parameter :: cal_t     = 3
    integer, parameter :: cleaning  = 4
    integer, parameter :: co2_ref   = 9
    integer, parameter :: h2o_ref   = 10
    integer, parameter :: ch4_ref   = 11
    integer, parameter :: gas4_ref  = 12
    integer :: open_status
    integer :: read_status
    integer :: cnt
    integer :: i
    integer :: j
    integer :: mdcol(GHGNumVar * 2)
    integer :: var_num
    integer :: sepa
    real(kind = dbl) :: biased_abs
    real(kind = dbl) :: unbiased_abs
    character(32) :: text_vars(128)
    character(LongInstringLen) :: dataline
    logical :: co2_bias_is_negative
    logical :: h2o_bias_is_negative


    !> Open dynamic metadata file
    write(*,'(a)') ' Looking for calibration information &
        &in dynamic metadata file:  '
    write(*,'(a)') '  ' // AuxFile%DynMD(1:len_trim(AuxFile%DynMD))
    open(udf, file = AuxFile%DynMD, status = 'old', iostat = open_status)

    Calib = CalibType('', '', tsNull, .false., nint(error), &
        error, error, error, error, error, error)

    nCalibEvents = 0
    mdcol = 0
    if (open_status == 0) then
        !> Look in the header for CO2/H2O calibration information
        read(udf, '(a)', iostat = read_status) dataline
        cnt = 0
        do
            sepa = index(dataline, ',')
            if (sepa == 0) sepa = len_trim(dataline) + 1
            if (len_trim(dataline) == 0) exit
            cnt = cnt + 1
            !> Date and time
            if (dataline(1:sepa - 1) &
                == trim(adjustl(StdDynMDVars(dynmd_date)))) &
                mdcol(cal_date) = cnt
            if (dataline(1:sepa - 1) &
                == trim(adjustl(StdDynMDVars(dynmd_time)))) &
                mdcol(cal_time) = cnt
            !> calibration data columns
            if (dataline(1:sepa - 1) == 'co2_offset')  mdcol(co2)   = cnt
            if (dataline(1:sepa - 1) == 'h2o_offset')  mdcol(h2o)   = cnt
            if (dataline(1:sepa - 1) == 'ch4_offset')  mdcol(ch4)   = cnt
            if (dataline(1:sepa - 1) == 'gas4_offset') mdcol(gas4)  = cnt
            if (dataline(1:sepa - 1) == 'co2_ref')  mdcol(co2_ref)  = cnt
            if (dataline(1:sepa - 1) == 'h2o_ref')  mdcol(h2o_ref)  = cnt
            if (dataline(1:sepa - 1) == 'ch4_ref')  mdcol(ch4_ref)  = cnt
            if (dataline(1:sepa - 1) == 'gas4_ref') mdcol(gas4_ref) = cnt
            if (dataline(1:sepa - 1) == 'calibration_temperature') &
                mdcol(cal_t) = cnt
            if (dataline(1:sepa - 1) == 'cleaning') mdcol(cleaning) = cnt

            dataline = dataline(sepa + 1: len_trim(dataline))
        end do

        if (sum(mdcol(co2:gas4)) == 0) then
            close(udf)
            return
        end if

        !> Reads file and store available calibration data
        i = 0
        do
            read(udf, '(a)', iostat = read_status) dataline
            if (read_status /= 0) exit
            i = i + 1

            !> For each data line, store data in a temporary array as text
            text_vars = 'none'
            var_num = 0
            do
                sepa = index(dataline, separator)
                if (sepa == 0) sepa = len_trim(dataline) + 1
                if (len_trim(dataline) == 0) exit
                var_num = var_num + 1
                text_vars(var_num) = dataline(1:sepa - 1)
                dataline = dataline(sepa + 1: len_trim(dataline))
            end do

            !> Associate stored data to relevant variables, for current dataline
            if (mdcol(cleaning) /= 0 .and. text_vars(mdcol(cleaning)) == '1') &
                Calib(i)%cleaning = .true.
            if (mdcol(cal_date) /= 0) read(text_vars(mdcol(cal_date)), *) &
                Calib(i)%date
            if (mdcol(cal_time) /= 0) read(text_vars(mdcol(cal_time)), *) &
                Calib(i)%time
            if (mdcol(cal_date) /= 0 .and. mdcol(cal_time) /= 0) &
                call DateTimeToDateType(Calib(i)%date, &
                    Calib(i)%time, Calib(i)%ts)

            if (mdcol(cal_t)    /= 0) &
                read(text_vars(mdcol(cal_t)), *) Calib(i)%Tcell
            do j = co2, gas4
                if (mdcol(j) /= 0) &
                read(text_vars(mdcol(j)), *) Calib(i)%offset(j)
            end do
            do j = co2_ref, gas4_ref
                if (mdcol(j) /= 0) &
                read(text_vars(mdcol(j)), *) Calib(i)%ref(j - 4)
            end do
        end do
        nCalibEvents = i
    end if
    close(udf)


    !> Convert offsets into absorptance offsets, considering the error as a
    !> span error thus, evaluating the abs_offset on the cal curve starting
    !> from the reference concentration indicated by the user
    do i = 1, nCalibEvents
        co2_bias_is_negative = .false.
        h2o_bias_is_negative = .false.

        !> co2
        if (Calib(i)%ref(co2) + Calib(i)%offset(co2) < 0d0) &
            co2_bias_is_negative = .true.

        if (mdcol(co2) /= 0 .and. mdcol(co2_ref) /= 0) &
            call PolyVal(DriftCorr%inv_cal(0:6, co2), 6, &
                dabs(Calib(i)%ref(co2) + Calib(i)%offset(co2)) &
                / Calib(i)%Tcell / Ru, &
                1, biased_abs)
        if (co2_bias_is_negative) biased_abs = - biased_abs

        if (mdcol(co2) /= 0 .and. mdcol(co2_ref) /= 0) &
            call PolyVal(DriftCorr%inv_cal(0:6, co2), 6, &
                Calib(i)%ref(co2) / Calib(i)%Tcell / Ru, &
                1, unbiased_abs)

        Calib(i)%offset(co2) = biased_abs - unbiased_abs

        !> h2o
        if (Calib(i)%ref(h2o) + Calib(i)%offset(h2o) < 0d0) &
            h2o_bias_is_negative = .true.

        if (mdcol(h2o) /= 0 .and. mdcol(h2o_ref) /= 0) &
            call PolyVal(DriftCorr%inv_cal(0:6, h2o), &
            6, dabs(Calib(i)%ref(h2o) + Calib(i)%offset(h2o)) &
            / Calib(i)%Tcell / Ru * 1d3, &
            1, biased_abs)
        if (h2o_bias_is_negative) biased_abs = - biased_abs

        if (mdcol(h2o) /= 0 .and. mdcol(h2o_ref) /= 0) &
            call PolyVal(DriftCorr%inv_cal(0:6, h2o), &
            6, Calib(i)%ref(h2o) / Calib(i)%Tcell / Ru * 1d3, &
            1, unbiased_abs)

        Calib(i)%offset(h2o) = biased_abs - unbiased_abs
    end do
    write(*,'(a)') ' Done.'
end subroutine driftRetrieveCalibrationEvents
