!***************************************************************************
! write_outfiles.f90
! ------------------
! Copyright (C) 2007-2011, Eco2s team, Gerardo Fratini
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
!***************************************************************************
!
! \brief       Write all results on (temporary) output files
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WriteIcosOutputRp(StDiff, DtDiff, STFlg, DTFlg)
    use m_rp_global_var
    implicit none
    !> in/out variables
    type(QCType), intent(in) :: StDiff
    type(QCType), intent(in) :: DtDiff
    integer, intent(in) :: STFlg(GHGNumVar)
    integer, intent(in) :: DTFlg(GHGNumVar)
    !> local variables
    integer :: var
    integer :: gas
    integer :: gas1
    integer :: gas2
    integer :: j
    integer :: i
    integer :: indx
    real(kind = dbl), allocatable :: bAggrOut(:)
    character(16000) :: dataline
    character(14) :: tsIso
    include '../src_common/interfaces.inc'

    !> write ICOS output file (csv) 
    call clearstr(dataline)

    !> Start/end imestamps
    tsIso = Stats%start_date(1:4) // Stats%start_date(6:7) // Stats%start_date(9:10) &
                // Stats%start_time(1:2) // Stats%start_time(4:5)
    call AddDatum(dataline, trim(adjustl(tsIso)), separator)
    tsIso = Stats%date(1:4) // Stats%date(6:7) // Stats%date(9:10) &
                // Stats%time(1:2) // Stats%time(4:5)
    call AddDatum(dataline, trim(adjustl(tsIso)), separator)

    !> Potential Radiations
    indx = DateTimeToHalfHourNumber(Stats%date, Stats%time)
    call AddFloatDatumToDataline(PotRad(indx), dataline, EddyProProj%err_label)

    !> Daytime
    if (Stats%daytime) then
        call AddDatum(dataline, '0', separator)
    else
        call AddDatum(dataline, '1', separator)
    endif

    !> Number of records
    !> Number of records teoretically available for current Averaging Interval
    call AddIntDatumToDataline(MaxPeriodNumRecords, dataline, EddyProProj%err_label)
    !> Number of records actually available for current Averaging Interval given length of actual files
    call AddIntDatumToDataline(Essentials%n_in, dataline, EddyProProj%err_label)
    !> Number of records actually available after custom flags filtering
    call AddIntDatumToDataline(Essentials%n_after_custom_flags, dataline, EddyProProj%err_label)
    !> Number of records actually available after wind direction filtering
    call AddIntDatumToDataline(Essentials%n_after_wdf, dataline, EddyProProj%err_label)
    !> Number of valid records for anemometric data
    call AddIntDatumToDataline(Essentials%n(w), dataline, EddyProProj%err_label)
    !> Number of valid records for IRGA data  (N_in – M_diag_IRGA)
    do var = ts, gas4
        call AddIntDatumToDataline(Essentials%n(var), dataline, EddyProProj%err_label)
        end do
    !> Number of valid records available for each main covariance (w/u, w/ts, w/co2, w/h2o, w/ch4, w/gas4)
    call AddIntDatumToDataline(Essentials%n_wcov(u), dataline, EddyProProj%err_label)
    do var = ts, gas4
        call AddIntDatumToDataline(Essentials%n_wcov(var), dataline, EddyProProj%err_label)
        end do

    !> Fluxes
    !> Fluxes level 3 (final fluxes) 
    call AddFloatDatumToDataline(Flux3%tau, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Flux3%H, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Flux3%LE, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Flux3%co2, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Flux3%h2o, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Flux3%ch4, dataline, EddyProProj%err_label, gain=1d3, offset=0d0)
    call AddFloatDatumToDataline(Flux3%gas4, dataline, EddyProProj%err_label, gain=1d3, offset=0d0)

    !> Flux random uncertainties
    if (Essentials%rand_uncer(u) == aflx_error) then
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    else
        call AddFloatDatumToDataline(Essentials%rand_uncer(u), dataline, EddyProProj%err_label)
        end if

    if (Essentials%rand_uncer(ts) == aflx_error) then
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    else
        call AddFloatDatumToDataline(Essentials%rand_uncer(ts), dataline, EddyProProj%err_label)
        end if

    if (Essentials%rand_uncer_LE == aflx_error) then
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    else
        call AddFloatDatumToDataline(Essentials%rand_uncer_LE, dataline, EddyProProj%err_label)
        end if

    if (Essentials%rand_uncer(co2) == aflx_error) then
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    else
        call AddFloatDatumToDataline(Essentials%rand_uncer(co2), dataline, EddyProProj%err_label)
        end if

    if (Essentials%rand_uncer(h2o) == aflx_error) then
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    else
        call AddFloatDatumToDataline(Essentials%rand_uncer(h2o), dataline, EddyProProj%err_label)
        end if

    if (Essentials%rand_uncer(ch4) == aflx_error) then
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    else
        call AddFloatDatumToDataline(Essentials%rand_uncer(ch4) * 1d3, dataline, EddyProProj%err_label)
        end if

    if (Essentials%rand_uncer(gas4) == aflx_error) then
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    else
        call AddFloatDatumToDataline(Essentials%rand_uncer(gas4) * 1d3, dataline, EddyProProj%err_label)
        end if

    !> Additional flux terms (single-point calculation)
    !> Storage fluxes
    call AddFloatDatumToDataline(Stor%H, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Stor%LE, dataline, EddyProProj%err_label)
    do gas = co2, h2o
        call AddFloatDatumToDataline(Stor%of(gas), dataline, EddyProProj%err_label)
        end do
    do gas = ch4, gas4
        call AddFloatDatumToDataline(Stor%of(gas), dataline, EddyProProj%err_label, gain=1d3, offset=0d0)
    end do
    !> Advection fluxes
    do gas = co2, gas4
        if (Stats5%Mean(w) /= error .and. Stats%d(gas) >= 0d0) then
            if (Stats5%Mean(w) /= error .and. Stats%d(gas) /= error) then
                if (gas == co2) then
                    call AddFloatDatumToDataline(Stats5%Mean(w) * Stats%d(gas), &
                        dataline, EddyProProj%err_label, gain=1d3, offset=0d0)
                else if (gas == h2o) then
                    call AddFloatDatumToDataline(Stats5%Mean(w) * Stats%d(gas), dataline, EddyProProj%err_label)
                else if (gas == ch4 .or. gas == gas4) then
                    call AddFloatDatumToDataline(Stats5%Mean(w) * Stats%d(gas), &
                        dataline, EddyProProj%err_label, gain=1d6, offset=0d0)
                end if
            else
                call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            end if
        else
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if
    end do

    !> Turbulence and micromet
    !> Unrotated and rotated wind components
    call AddFloatDatumToDataline(Stats4%Mean(u), dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Stats4%Mean(v), dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Stats4%Mean(w), dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Stats5%Mean(u), dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Stats5%Mean(v), dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Stats5%Mean(w), dataline, EddyProProj%err_label)
    !> wind speed, wind direction, u*, stability, bowen ratio
    call AddFloatDatumToDataline(Ambient%WS, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Ambient%MWS, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Stats4%wind_dir, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Stats4%wind_dir_stdev, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Ambient%us, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Stats%TKE, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Ambient%L, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Ambient%zL, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Ambient%bowen, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Ambient%Ts, dataline, EddyProProj%err_label)

    !> Termodynamics 
    !> Temperature, pressure, RH, VPD, e, es, etc.
    call AddFloatDatumToDataline(Stats7%Mean(ts), dataline, EddyProProj%err_label, gain=1d0, offset=-273.15d0)
    call AddFloatDatumToDataline(Ambient%Ta, dataline, EddyProProj%err_label, gain=1d0, offset=-273.15d0)
    call AddFloatDatumToDataline(Stats%Pr, dataline, EddyProProj%err_label, gain=1d-3, offset=0d0)
    call AddFloatDatumToDataline(Stats%RH, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Ambient%Va, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(RHO%a, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Ambient%RhoCp, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(RHO%w, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Ambient%e, dataline, EddyProProj%err_label, gain=1d-2, offset=0d0)
    call AddFloatDatumToDataline(Ambient%es, dataline, EddyProProj%err_label, gain=1d-2, offset=0d0)
    call AddFloatDatumToDataline(Ambient%Q, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Ambient%VPD, dataline, EddyProProj%err_label, gain=1d-2, offset=0d0)
    call AddFloatDatumToDataline(Ambient%Td, dataline, EddyProProj%err_label, gain=1d0, offset=-273.15d0)
    !> Dry air properties
    call AddFloatDatumToDataline(Ambient%p_d, dataline, EddyProProj%err_label, gain=1d-3, offset=0d0)
    call AddFloatDatumToDataline(RHO%d, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Ambient%Vd, dataline, EddyProProj%err_label)
    !> Specific heat of evaporation
    call AddFloatDatumToDataline(Ambient%lambda, dataline, EddyProProj%err_label)
    !> Wet to dry air density ratio
    call AddFloatDatumToDataline(Ambient%sigma, dataline, EddyProProj%err_label)
    !> Water Use Efficiency
    !>!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*

    !> Gases
    !> Concentrations, densities and "nature" of the raw data 
    !> (mixing ratio, mole fraction, molar density)
    !> Gas concentrations, densities and timelags
    do gas = co2, gas4
        if (E2Col(gas)%present) then
            select case (E2Col(gas)%measure_type)
                case('mixing_ratio')
                    call AddDatum(dataline, '0', separator)
                case('mole_fraction')
                    call AddDatum(dataline, '1', separator)
                case('molar_density')
                    call AddDatum(dataline, '2', separator)
            end select
            call AddFloatDatumToDataline(Stats%d(gas), dataline, EddyProProj%err_label)
            if (gas == ch4 .or. gas == gas4) then
                call AddFloatDatumToDataline(Stats%r(gas), dataline, EddyProProj%err_label, gain=1d3, offset=0d0)
                call AddFloatDatumToDataline(Stats%chi(gas), dataline, EddyProProj%err_label, gain=1d3, offset=0d0)
            else
                call AddFloatDatumToDataline(Stats%r(gas), dataline, EddyProProj%err_label)
                call AddFloatDatumToDataline(Stats%chi(gas), dataline, EddyProProj%err_label)
            end if
        else
            do i = 1, 4
                call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            end do
        end if
    end do
    !> Timelags (calculated, used, min/max/nominal) for all gases
    !> Gas timelags
    do gas = co2, gas4
        if (E2Col(gas)%present) then
            call AddFloatDatumToDataline(Essentials%actual_timelag(gas), dataline, EddyProProj%err_label)
                    call AddFloatDatumToDataline(Essentials%used_timelag(gas), dataline, EddyProProj%err_label)
                    call AddFloatDatumToDataline(E2Col(gas)%def_tl, dataline, EddyProProj%err_label)
                    call AddFloatDatumToDataline(E2Col(gas)%min_tl, dataline, EddyProProj%err_label)
                    call AddFloatDatumToDataline(E2Col(gas)%max_tl, dataline, EddyProProj%err_label)
                else
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if
    end do

!> Basic stats
    !> 25-50-75%
    do var = u, gas4
        if (var == ts) then
            call AddFloatDatumToDataline(Stats6%Median(var), dataline, &
                EddyProProj%err_label, gain=1d0, offset=-273.15d0)
        else
            call AddFloatDatumToDataline(Stats6%Median(var), dataline, EddyProProj%err_label)
                end if
    end do
    do var = u, gas4
        if (var == ts) then
            call AddFloatDatumToDataline(Stats6%Q1(var), dataline, &
                EddyProProj%err_label, gain=1d0, offset=-273.15d0)
        else
            call AddFloatDatumToDataline(Stats6%Q1(var), dataline, EddyProProj%err_label)
                end if
    end do
    do var = u, gas4
        if (var == ts) then
            call AddFloatDatumToDataline(Stats6%Q3(var), dataline, &
                EddyProProj%err_label, gain=1d0, offset=-273.15d0)
        else
            call AddFloatDatumToDataline(Stats6%Q3(var), dataline, EddyProProj%err_label)
                end if
    end do
    !> Standard deviation
    do var = u, gas4
        call AddFloatDatumToDataline(sqrt(Stats7%Cov(var, var)), dataline, EddyProProj%err_label)
        end do
    !> Skwenesses
    do var = u, gas4
        call AddFloatDatumToDataline(Stats7%Skw(var), dataline, EddyProProj%err_label)
        end do
    !> Kurtosis
    do var = u, gas4
        call AddFloatDatumToDataline(Stats7%Kur(var), dataline, EddyProProj%err_label)
        end do
    !> w-covariances 
    call AddFloatDatumToDataline(Stats7%Cov(u, w), dataline, EddyProProj%err_label)
    do var = ts, gas4
        call AddFloatDatumToDataline(Stats7%Cov(w, var), dataline, EddyProProj%err_label)
        end do
    !> Gases covariance matrix
    do gas1 = co2, ch4
        do gas2 = gas1 + 1, gas4 
            call AddFloatDatumToDataline(Stats7%Cov(gas1, gas2), dataline, EddyProProj%err_label)
                end do
    end do

!> Footprint
    call AddFloatDatumToDataline(Foot%peak, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Foot%offset, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Foot%x10, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Foot%x30, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Foot%x50, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Foot%x70, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Foot%x80, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Foot%x90, dataline, EddyProProj%err_label)

!> Intermediate results
    !> Fluxes level 0 (uncorrected fluxes)
    call AddFloatDatumToDataline(Essentials%L, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Essentials%zL, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Flux0%tau, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Flux0%H, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Flux0%LE, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Flux0%co2, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Flux0%h2o, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Flux0%ch4, dataline, EddyProProj%err_label, gain=1d3, offset=0d0)
    call AddFloatDatumToDataline(Flux0%gas4, dataline, EddyProProj%err_label, gain=1d3, offset=0d0)
    !> Fluxes level 1
    call AddFloatDatumToDataline(Flux1%tau, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Flux1%H, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Flux1%LE, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Flux1%co2, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Flux1%h2o, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Flux1%ch4, dataline, EddyProProj%err_label, gain=1d3, offset=0d0)
    call AddFloatDatumToDataline(Flux1%gas4, dataline, EddyProProj%err_label, gain=1d3, offset=0d0)
    !> Fluxes level 2
    call AddFloatDatumToDataline(Flux2%tau, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Flux2%H, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Flux2%LE, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Flux2%co2, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Flux2%h2o, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Flux2%ch4, dataline, EddyProProj%err_label, gain=1d3, offset=0d0)
    call AddFloatDatumToDataline(Flux2%gas4, dataline, EddyProProj%err_label, gain=1d3, offset=0d0)

    !> Tin and Tout                 ******************************************** Add

    !> Temperature, pressure and molar volume 
    !> in the cell of closed-paths, for all gases
    !> Cell parameters              ******************************************** Mke it gas specific like molar volume
    call AddFloatDatumToDataline(Ambient%Tcell, dataline, EddyProProj%err_label, gain=1d0, offset=-273.15d0)
    call AddFloatDatumToDataline(Ambient%Pcell, dataline, EddyProProj%err_label, gain=1d-3, offset=0d0)

    !> Molar volume
    do gas = co2, gas4
        call AddFloatDatumToDataline(E2Col(gas)%Va, dataline, EddyProProj%err_label)
        end do
    !> Evapotranspiration and sensible heat fluxes in the cell of 
    !> closed-paths (for WPL), with timelags of other gases
    call AddFloatDatumToDataline(Flux3%E_co2, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Flux3%E_ch4, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Flux3%E_gas4, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Flux3%Hi_co2, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Flux3%Hi_h2o, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Flux3%Hi_ch4, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Flux3%Hi_gas4, dataline, EddyProProj%err_label)
    !> Burba Terms 
    if (RPsetup%bu_corr /= 'none') then 
        call AddFloatDatumToDataline(Burba%h_bot, dataline, EddyProProj%err_label)
            call AddFloatDatumToDataline(Burba%h_top, dataline, EddyProProj%err_label)
            call AddFloatDatumToDataline(Burba%h_spar, dataline, EddyProProj%err_label)
        else
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    end if
    !> LI-7700 multipliers
    if (E2Col(ch4)%Instr%model(1:len_trim(E2Col(ch4)%Instr%model) - 2) &
        == 'li7700') then
        call AddFloatDatumToDataline(Mul7700%A, dataline, EddyProProj%err_label)
            call AddFloatDatumToDataline(Mul7700%B, dataline, EddyProProj%err_label)
            call AddFloatDatumToDataline(Mul7700%C, dataline, EddyProProj%err_label)
        else
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    end if
    !> WPL Terms                    ********************************************(Individual: H, LE, Pressure)
    !>!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*
    !> Spectral correction factors
    if (E2Col(u)%present) then
        call AddFloatDatumToDataline(BPCF%of(w_u), dataline, EddyProProj%err_label)
        else
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    end if
    if (E2Col(ts)%present) then
        call AddFloatDatumToDataline(BPCF%of(w_ts), dataline, EddyProProj%err_label)
        else
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    end if
    if (E2Col(h2o)%present) then
        call AddFloatDatumToDataline(BPCF%of(w_h2o), dataline, EddyProProj%err_label)
        else
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    end if
    if (E2Col(co2)%present) then
        call AddFloatDatumToDataline(BPCF%of(w_co2), dataline, EddyProProj%err_label)
        else
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    end if
    if (E2Col(h2o)%present) then
        call AddFloatDatumToDataline(BPCF%of(w_h2o), dataline, EddyProProj%err_label)
        else
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    end if
    if (E2Col(ch4)%present) then
        call AddFloatDatumToDataline(BPCF%of(w_ch4), dataline, EddyProProj%err_label)
        else
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    end if
    if (E2Col(gas4)%present) then
        call AddFloatDatumToDataline(BPCF%of(w_gas4), dataline, EddyProProj%err_label)
        else
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    end if

    !> Increasingly filtered w/T covariances (for spectral assessment)
    call AddFloatDatumToDataline(Essentials%degH(NumDegH + 1), dataline, EddyProProj%err_label)
    do j = 1, NumDegH
        call AddFloatDatumToDataline(Essentials%degH(j), dataline, EddyProProj%err_label)
    end do

!> QC details
    !>> Number or records eliminated based on custom flags
    call AddIntDatumToDataline(Essentials%m_custom_flags, dataline, EddyProProj%err_label)
    !>> Number or records eliminated based on wind direction filter
    call AddIntDatumToDataline(Essentials%m_wdf, dataline, EddyProProj%err_label)
    !> Summary of data values eliminated based on diagnostics
    !>> Number or records whose anemometric data was eliminated based on Anemometer diagnostics
    call AddIntDatumToDataline(Essentials%m_diag_anem, dataline, EddyProProj%err_label)
    !>> Number or records whose IRGA data was eliminated based on IRGA diagnostics
    do gas = co2, gas4
        call AddIntDatumToDataline(Essentials%m_diag_irga(gas), dataline, EddyProProj%err_label)
        end do
    !>> Number of values eliminated by the Spike test
    do var = u, gas4
        call AddIntDatumToDataline(Essentials%m_despiking(var), dataline, EddyProProj%err_label)
        end do
    !>> Number of values eliminated by the Absolute Limits test
    do j = u, gas4
        call AddIntDatumToDataline(Essentials%al_s(j), dataline, EddyProProj%err_label)
        end do
    !> VM97 Stats used to calculate flags
    !>> Spikes
    do j = u, gas4
        call AddIntDatumToDataline(Essentials%e2spikes(j), dataline, EddyProProj%err_label)
        end do
    !>> Amplitude resolution
    do j = u, gas4
        call AddFloatDatumToDataline(Essentials%ar_s(j), dataline, EddyProProj%err_label)
        end do
    !>> Dropouts central
    do j = u, gas4
        call AddFloatDatumToDataline(Essentials%do_s_ctr(j), dataline, EddyProProj%err_label)
        end do
    !>> Dropouts extremes
    do j = u, gas4
        call AddFloatDatumToDataline(Essentials%do_s_ext(j), dataline, EddyProProj%err_label)
        end do
    !>> Higher moments Skewness
    do j = u, gas4
        call AddFloatDatumToDataline(Essentials%sk_s_skw(j), dataline, EddyProProj%err_label)
        end do                          
    !>> Higher moments Kurtosis     
    do j = u, gas4
        call AddFloatDatumToDataline(Essentials%sk_s_kur(j), dataline, EddyProProj%err_label)
        end do
    !>> AoA
    call AddFloatDatumToDataline(Essentials%aa_s, dataline, EddyProProj%err_label)
    !>> Non-steady wind
    call AddFloatDatumToDataline(Essentials%ns_s_rnv(1), dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Essentials%ns_s_rnv(2), dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Essentials%ns_s_rns, dataline, EddyProProj%err_label)
    !> VM97 flags
    call AddDatum(dataline, '8'//CharHF%sr(2:9), separator)
    call AddDatum(dataline, '8'//CharHF%ar(2:9), separator)
    call AddDatum(dataline, '8'//CharHF%do(2:9), separator)
    call AddDatum(dataline, '8'//CharHF%al(2:9), separator)
    call AddDatum(dataline, '8'//CharHF%sk(2:9), separator)
    call AddDatum(dataline, '8'//CharSF%sk(2:9), separator)
    call AddDatum(dataline, '8'//CharHF%ds(2:9), separator)
    call AddDatum(dataline, '8'//CharSF%ds(2:9), separator)
    call AddDatum(dataline, '8'//CharHF%tl(6:9), separator)
    call AddDatum(dataline, '8'//CharSF%tl(6:9), separator)
    call AddDatum(dataline, '8'//CharHF%aa(9:9), separator)
    call AddDatum(dataline, '8'//CharHF%ns(9:9), separator)

    !> Quality test results
    !> Foken stats used to calculate flags
    call AddIntDatumToDataline(STDiff%w_u, dataline, EddyProProj%err_label)
    call AddIntDatumToDataline(STDiff%w_ts, dataline, EddyProProj%err_label)
    call AddIntDatumToDataline(StDiff%w_co2, dataline, EddyProProj%err_label)
    call AddIntDatumToDataline(StDiff%w_h2o, dataline, EddyProProj%err_label)
    call AddIntDatumToDataline(StDiff%w_ch4, dataline, EddyProProj%err_label)
    call AddIntDatumToDataline(StDiff%w_gas4, dataline, EddyProProj%err_label)
    call AddIntDatumToDataline(DtDiff%u, dataline, EddyProProj%err_label)
    call AddIntDatumToDataline(DtDiff%w, dataline, EddyProProj%err_label)
    call AddIntDatumToDataline(DtDiff%ts, dataline, EddyProProj%err_label)
    !> Partial Foken flags
    call AddIntDatumToDataline(STFlg(w_u), dataline, EddyProProj%err_label)
    call AddIntDatumToDataline(STFlg(w_ts), dataline, EddyProProj%err_label)
    call AddIntDatumToDataline(STFlg(w_co2), dataline, EddyProProj%err_label)
    call AddIntDatumToDataline(STFlg(w_h2o), dataline, EddyProProj%err_label)
    call AddIntDatumToDataline(STFlg(w_ch4), dataline, EddyProProj%err_label)
    call AddIntDatumToDataline(STFlg(w_gas4), dataline, EddyProProj%err_label)
    call AddIntDatumToDataline(DTFlg(u), dataline, EddyProProj%err_label)
    call AddIntDatumToDataline(DTFlg(w), dataline, EddyProProj%err_label)
    call AddIntDatumToDataline(DTFlg(ts), dataline, EddyProProj%err_label)
    !> Final Foken flags
    call AddIntDatumToDataline(QCFlag%tau, dataline, EddyProProj%err_label)
    call AddIntDatumToDataline(QCFlag%H, dataline, EddyProProj%err_label)
    call AddIntDatumToDataline(QCFlag%h2o, dataline, EddyProProj%err_label)
    call AddIntDatumToDataline(QCFlag%co2, dataline, EddyProProj%err_label)
    call AddIntDatumToDataline(QCFlag%h2o, dataline, EddyProProj%err_label)
    call AddIntDatumToDataline(QCFlag%ch4, dataline, EddyProProj%err_label)
    call AddIntDatumToDataline(QCFlag%gas4, dataline, EddyProProj%err_label)
    !> Number of calculated spikes per variables       ************************* Evaluate if needed (what if no spike removal selected?)
    call AddIntDatumToDataline(Essentials%e2spikes(u), dataline, EddyProProj%err_label)
    call AddIntDatumToDataline(Essentials%e2spikes(v), dataline, EddyProProj%err_label)
    call AddIntDatumToDataline(Essentials%e2spikes(w), dataline, EddyProProj%err_label)
    call AddIntDatumToDataline(Essentials%e2spikes(ts), dataline, EddyProProj%err_label)
    do var = co2, gas4
        call AddIntDatumToDataline(Essentials%e2spikes(var), dataline, EddyProProj%err_label)
        end do

    !> LI-7x00 diagnostics breakdown
    if (Diag7200%present) then
        call AddIntDatumToDataline(Diag7200%head_detect, dataline, EddyProProj%err_label)
            call AddIntDatumToDataline(Diag7200%t_out, dataline, EddyProProj%err_label)
            call AddIntDatumToDataline(Diag7200%t_in, dataline, EddyProProj%err_label)
            call AddIntDatumToDataline(Diag7200%aux_in, dataline, EddyProProj%err_label)
            call AddIntDatumToDataline(Diag7200%delta_p, dataline, EddyProProj%err_label)
            call AddIntDatumToDataline(Diag7200%chopper, dataline, EddyProProj%err_label)
            call AddIntDatumToDataline(Diag7200%detector, dataline, EddyProProj%err_label)
            call AddIntDatumToDataline(Diag7200%pll, dataline, EddyProProj%err_label)
            call AddIntDatumToDataline(Diag7200%sync, dataline, EddyProProj%err_label)
        else
        do i = 1, 9
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end do
    end if
    if (Diag7500%present) then
        call AddIntDatumToDataline(Diag7500%chopper, dataline, EddyProProj%err_label)
            call AddIntDatumToDataline(Diag7500%detector, dataline, EddyProProj%err_label)
            call AddIntDatumToDataline(Diag7500%pll, dataline, EddyProProj%err_label)
            call AddIntDatumToDataline(Diag7500%sync, dataline, EddyProProj%err_label)
        else
        do i = 1, 4
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end do
    end if
    if (Diag7700%present) then
        call AddIntDatumToDataline(Diag7700%not_ready, dataline, EddyProProj%err_label)
            call AddIntDatumToDataline(Diag7700%no_signal, dataline, EddyProProj%err_label)
            call AddIntDatumToDataline(Diag7700%re_unlocked, dataline, EddyProProj%err_label)
            call AddIntDatumToDataline(Diag7700%bad_temp, dataline, EddyProProj%err_label)
            call AddIntDatumToDataline(Diag7700%laser_temp_unregulated, dataline, EddyProProj%err_label)
            call AddIntDatumToDataline(Diag7700%block_temp_unregulated, dataline, EddyProProj%err_label)
            call AddIntDatumToDataline(Diag7700%motor_spinning, dataline, EddyProProj%err_label)
            call AddIntDatumToDataline(Diag7700%pump_on, dataline, EddyProProj%err_label)
            call AddIntDatumToDataline(Diag7700%top_heater_on, dataline, EddyProProj%err_label)
            call AddIntDatumToDataline(Diag7700%bottom_heater_on, dataline, EddyProProj%err_label)
            call AddIntDatumToDataline(Diag7700%calibrating, dataline, EddyProProj%err_label)
            call AddIntDatumToDataline(Diag7700%motor_failure, dataline, EddyProProj%err_label)
            call AddIntDatumToDataline(Diag7700%bad_aux_tc1, dataline, EddyProProj%err_label)
            call AddIntDatumToDataline(Diag7700%bad_aux_tc2, dataline, EddyProProj%err_label)
            call AddIntDatumToDataline(Diag7700%bad_aux_tc3, dataline, EddyProProj%err_label)
            call AddIntDatumToDataline(Diag7700%box_connected, dataline, EddyProProj%err_label)
        else
        do i = 1, 16
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end do
    end if
    !> AGC/RSSI                         **************************************** May need to adapt header to whether it's AGC or RSSI for 7200/7500
    if(CompareSwVer(E2Col(co2)%instr%sw_ver, SwVerFromString('6.0.0'))) then
        call AddIntDatumToDataline(nint(Essentials%AGC72), dataline, EddyProProj%err_label)
    else
        call AddIntDatumToDataline(-nint(Essentials%AGC72), dataline, EddyProProj%err_label)
    end if
    !> LI-7500
    if(CompareSwVer(E2Col(co2)%instr%sw_ver, SwVerFromString('6.0.0'))) then
        call AddIntDatumToDataline(nint(Essentials%AGC75), dataline, EddyProProj%err_label)
    else
        call AddIntDatumToDataline(-nint(Essentials%AGC75), dataline, EddyProProj%err_label)
    end if
    !> LI-7700
    call AddIntDatumToDataline(nint(Essentials%RSSI77), dataline, EddyProProj%err_label)

!> Processing settings
    !> Whether w-boost calibration was applied
    if (RPsetup%calib_wboost) then
        call AddIntDatumToDataline(1, dataline, EddyProProj%err_label)
    else
        call AddIntDatumToDataline(0, dataline, EddyProProj%err_label)
    end if
    !> Whether AoA calibration was applied
    select case(trim(adjustl(RPsetup%calib_aoa)))
        case('automatic')
            call AddIntDatumToDataline(-1, dataline, EddyProProj%err_label)
        case('none')
            call AddIntDatumToDataline(0, dataline, EddyProProj%err_label)
        case('nakai_06')
            call AddIntDatumToDataline(1, dataline, EddyProProj%err_label)
        case('nakai_12')
            call AddIntDatumToDataline(2, dataline, EddyProProj%err_label)
    end select
    !> Tilt compensation method
    select case(trim(adjustl(Meth%rot)))
        case('none')
            call AddIntDatumToDataline(0, dataline, EddyProProj%err_label)
        case('double_rotation')
            call AddIntDatumToDataline(1, dataline, EddyProProj%err_label)
        case('triple_rotation')
            call AddIntDatumToDataline(2, dataline, EddyProProj%err_label)
        case('planar_fit')
            call AddIntDatumToDataline(3, dataline, EddyProProj%err_label)
        case('planar_fit_no_bias')
            call AddIntDatumToDataline(4, dataline, EddyProProj%err_label)
    end select
    !> Rotation angles
    call AddFloatDatumToDataline(Essentials%yaw, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Essentials%pitch, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Essentials%roll, dataline, EddyProProj%err_label)
    !> Detrending method and time constant
    select case(trim(adjustl(Meth%det)))
        case('ba')
            call AddIntDatumToDataline(0, dataline, EddyProProj%err_label)
        case('ld')
            call AddIntDatumToDataline(1, dataline, EddyProProj%err_label)
        case('rm')
            call AddIntDatumToDataline(2, dataline, EddyProProj%err_label)
        case('ew')
            call AddIntDatumToDataline(3, dataline, EddyProProj%err_label)
    end select
    call AddIntDatumToDataline(RPsetup%Tconst, dataline, EddyProProj%err_label)
    !> Time lag detection method
    select case(trim(adjustl(Meth%tlag)))
        case('none')
            call AddIntDatumToDataline(0, dataline, EddyProProj%err_label)
        case('constant')
            call AddIntDatumToDataline(1, dataline, EddyProProj%err_label)
        case('maxcov&default')
            call AddIntDatumToDataline(2, dataline, EddyProProj%err_label)
        case('maxcov')
            call AddIntDatumToDataline(3, dataline, EddyProProj%err_label)
        case('tlag_opt')
            call AddIntDatumToDataline(4, dataline, EddyProProj%err_label)
    end select
    !> WPL terms
    if (EddyProProj%wpl) then
        call AddIntDatumToDataline(1, dataline, EddyProProj%err_label)
    else
        call AddIntDatumToDataline(0, dataline, EddyProProj%err_label)
    end if
    !> Burba terms
    if (trim(adjustl(RPSetup%bu_corr)) == 'yes') then
        if (RPSetup%bu_multi) then
            call AddIntDatumToDataline(2, dataline, EddyProProj%err_label)
        else
            call AddIntDatumToDataline(1, dataline, EddyProProj%err_label)
        end if
    else
        call AddIntDatumToDataline(0, dataline, EddyProProj%err_label)
    end if
    !> Spectral correction method
    select case(trim(adjustl(EddyProProj%hf_meth)))
        case('none')
            call AddIntDatumToDataline(0, dataline, EddyProProj%err_label)
        case('moncrieff_97')
            call AddIntDatumToDataline(1, dataline, EddyProProj%err_label)
        case('horst_97')
            call AddIntDatumToDataline(2, dataline, EddyProProj%err_label)
        case('ibrom_07')
            call AddIntDatumToDataline(3, dataline, EddyProProj%err_label)
        case('fratini_12')
            call AddIntDatumToDataline(4, dataline, EddyProProj%err_label)
        case('massman_00')
            call AddIntDatumToDataline(5, dataline, EddyProProj%err_label)
    end select
    !> Footprint model
    select case(trim(adjustl(foot_model_used)))
        case('none')
            call AddIntDatumToDataline(0, dataline, EddyProProj%err_label)
        case('kljun_04')
            call AddIntDatumToDataline(1, dataline, EddyProProj%err_label)
        case('kormann_meixner_01')
            call AddIntDatumToDataline(2, dataline, EddyProProj%err_label)
        case('hsieh_00')
            call AddIntDatumToDataline(3, dataline, EddyProProj%err_label)
    end select

!> Metadata
    !> Data logger software version
    call AddIntDatumToDataline(Metadata%logger_swver%major, dataline, EddyProProj%err_label)
    call AddIntDatumToDataline(Metadata%logger_swver%minor, dataline, EddyProProj%err_label)
    call AddIntDatumToDataline(Metadata%logger_swver%revision, dataline, EddyProProj%err_label)
    !> Site location and features
    call AddFloatDatumToDataline(Metadata%lat, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Metadata%lon, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Metadata%alt, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Metadata%canopy_height, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Metadata%d, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Metadata%z0, dataline, EddyProProj%err_label)
    !> Data acquisition settings
    call AddIntDatumToDataline(nint(Metadata%file_length), dataline, EddyProProj%err_label)
    call AddIntDatumToDataline(nint(Metadata%ac_freq), dataline, EddyProProj%err_label)
    !> Flux averaging interval
    call AddIntDatumToDataline(RPsetup%avrg_len, dataline, EddyProProj%err_label)
    !> master anemometer
    call AddCharDatumToDataline(E2Col(u)%instr%firm, dataline, EddyProProj%err_label)
    call AddCharDatumToDataline(E2Col(u)%Instr%model, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(E2Col(u)%Instr%height, dataline, EddyProProj%err_label)
    call AddCharDatumToDataline(E2Col(u)%Instr%wformat, dataline, EddyProProj%err_label)
    call AddCharDatumToDataline(E2Col(u)%Instr%wref, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(E2Col(u)%Instr%north_offset, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(E2Col(u)%Instr%hpath_length, dataline, EddyProProj%err_label, gain=1d2, offset=0d0)
    call AddFloatDatumToDataline(E2Col(u)%Instr%vpath_length, dataline, EddyProProj%err_label, gain=1d2, offset=0d0)
    call AddFloatDatumToDataline(E2Col(u)%Instr%tau, dataline, EddyProProj%err_label)
    !> gas analysers details
    do gas = co2, gas4
        call AddCharDatumToDataline(E2Col(gas)%Instr%firm, dataline, EddyProProj%err_label)
        call AddCharDatumToDataline(E2Col(gas)%Instr%model, dataline, EddyProProj%err_label)
        call AddFloatDatumToDataline(E2Col(gas)%Instr%nsep, dataline, EddyProProj%err_label, gain=1d2, offset=0d0)
        call AddFloatDatumToDataline(E2Col(gas)%Instr%esep, dataline, EddyProProj%err_label, gain=1d2, offset=0d0)
        call AddFloatDatumToDataline(E2Col(gas)%Instr%vsep, dataline, EddyProProj%err_label, gain=1d2, offset=0d0)
        call AddFloatDatumToDataline(E2Col(gas)%Instr%tube_l, dataline, EddyProProj%err_label, gain=1d2, offset=0d0)
        call AddFloatDatumToDataline(E2Col(gas)%Instr%tube_d, dataline, EddyProProj%err_label, gain=1d3, offset=0d0)
        call AddFloatDatumToDataline(E2Col(gas)%Instr%tube_f, dataline, EddyProProj%err_label, gain=6d4, offset=0d0)
        if (gas == h2o) then
            call AddFloatDatumToDataline(E2Col(gas)%Instr%kw, dataline, EddyProProj%err_label)
            call AddFloatDatumToDataline(E2Col(gas)%Instr%ko, dataline, EddyProProj%err_label)
        end if
        call AddFloatDatumToDataline(E2Col(gas)%Instr%hpath_length, dataline, EddyProProj%err_label, gain=1d2, offset=0d0)
        call AddFloatDatumToDataline(E2Col(gas)%Instr%vpath_length, dataline, EddyProProj%err_label, gain=1d2, offset=0d0)
        call AddFloatDatumToDataline(E2Col(gas)%Instr%tau, dataline, EddyProProj%err_label)
    end do

    !> Number and mean values of custom variables
    call AddIntDatumToDataline(NumUserVar, dataline, EddyProProj%err_label)
    if (NumUserVar > 0) then
        do var = 1, NumUserVar
            call AddFloatDatumToDataline(UserStats%Mean(var), dataline, EddyProProj%err_label)
        end do
    end if

    !> All aggregated biomet values in FLUXNET units
    call AddIntDatumToDataline(nbVars, dataline, EddyProProj%err_label)
    if (nbVars > 0) then
        if (.not. allocated(bAggrOut)) allocate(bAggrOut(size(bAggr)))
        if (EddyProProj%icos_standardize_biomet) then
            bAggrOut = bAggrFluxnet
        else
            bAggrOut = bAggr
        end if

        do i = 1, nbVars
            call AddFloatDatumToDataline(bAggrOut(i), dataline, EddyProProj%err_label)
        end do

        if (allocated(bAggrOut)) deallocate(bAggrOut)
    end if
    write(uicos, '(a)') dataline(1:len_trim(dataline) - 1)
end subroutine WriteIcosOutputRp