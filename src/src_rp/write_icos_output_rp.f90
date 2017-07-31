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
!    integer :: prof
    character(16000) :: dataline
    character(DatumLen) :: datum
    character(14) :: iso_basic
    include '../src_common/interfaces.inc'

    !> write ICOS output file (csv) 
    call clearstr(dataline)

    !> Timestamp
    iso_basic = Stats%date(1:4) // Stats%date(6:7) // Stats%date(9:10) &
                // Stats%time(1:2) // Stats%time(4:5) // '00'
    call AddDatum(dataline, trim(adjustl(iso_basic)), separator)

!> Potential Radiations
    indx = DateTimeToHalfHourNumber(Stats%date, Stats%time)
    call WriteDatumFloat(PotRad(indx), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)

!> Daytime
    if (Stats%daytime) then
        call AddDatum(dataline, '1', separator)
    else
        call AddDatum(dataline, '0', separator)
    endif

!> Number of records
    !> Number of records teoretically available for current Averaging Interval
    call WriteDatumInt(MaxPeriodNumRecords, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    !> Number of records actually available for current Averaging Interval given length of actual files
    call WriteDatumInt(Essentials%n_in, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    !> Number of records actually available after custom flags filtering
    call WriteDatumInt(Essentials%n_after_custom_flags, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    !> Number of records actually available after wind direction filtering
    call WriteDatumInt(Essentials%n_after_wdf, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    !> Number of valid records for anemometric data
    call WriteDatumInt(Essentials%n(w), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    !> Number of valid records for IRGA data  (N_in – M_diag_IRGA)
    do var = ts, gas4
        call WriteDatumInt(Essentials%n(var), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end do
    !> Number of valid records available for each main covariance (w/u, w/ts, w/co2, w/h2o, w/ch4, w/gas4)
    call WriteDatumInt(Essentials%n_wcov(u), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    do var = ts, gas4
        call WriteDatumInt(Essentials%n_wcov(var), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end do

!> Fluxes
    !> Fluxes level 3 (final fluxes) 
    call WriteDatumFloat(Flux3%tau, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Flux3%H, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Flux3%LE, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Flux3%co2, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Flux3%h2o, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Flux3%ch4, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Flux3%gas4, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)

!> Flux random uncertainties
    if (Essentials%rand_uncer(u) == aflx_error) then
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    else
        call WriteDatumFloat(Essentials%rand_uncer(u), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end if

    if (Essentials%rand_uncer(ts) == aflx_error) then
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    else
        call WriteDatumFloat(Essentials%rand_uncer(ts), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end if

    if (Essentials%rand_uncer_LE == aflx_error) then
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    else
        call WriteDatumFloat(Essentials%rand_uncer_LE, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end if

    if (Essentials%rand_uncer(co2) == aflx_error) then
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    else
        call WriteDatumFloat(Essentials%rand_uncer(co2), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end if

    if (Essentials%rand_uncer(h2o) == aflx_error) then
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    else
        call WriteDatumFloat(Essentials%rand_uncer(h2o), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end if

    if (Essentials%rand_uncer(ch4) == aflx_error) then
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    else
        call WriteDatumFloat(Essentials%rand_uncer(ch4), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end if

    if (Essentials%rand_uncer(gas4) == aflx_error) then
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    else
        call WriteDatumFloat(Essentials%rand_uncer(gas4), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end if

!> Additional flux terms (single-point calculation)
    !> Storage fluxes
    call WriteDatumFloat(Stor%H, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Stor%LE, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    do gas = co2, gas4
        call WriteDatumFloat(Stor%of(gas), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end do
    !> Advection fluxes
    do gas = co2, gas4
        if (Stats5%Mean(w) /= error .and. Stats%d(gas) >= 0d0) then
            if (gas /= h2o) then
                call WriteDatumFloat(Stats5%Mean(w) * Stats%d(gas) * 1d3, datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            else
                call WriteDatumFloat(Stats5%Mean(w) * Stats%d(gas), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end if
        else
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if
    end do

!> Turbulence and micromet
    !> Unrotated and rotated wind components
    call WriteDatumFloat(Stats4%Mean(u), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Stats4%Mean(v), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Stats4%Mean(w), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Stats5%Mean(u), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Stats5%Mean(v), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Stats5%Mean(w), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    !> wind speed, wind direction, u*, stability, bowen ratio
    call WriteDatumFloat(Ambient%WS, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Ambient%MWS, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Stats4%wind_dir, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Ambient%us, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Stats%TKE, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Ambient%L, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Ambient%zL, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Ambient%bowen, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Ambient%Ts, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)

!> Termodynamics 
    !> Temperature, pressure, RH, VPD, e, es, etc.
    call WriteDatumFloat(Stats7%Mean(ts), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Ambient%Ta, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Stats%Pr, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Stats%RH, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Ambient%Va, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(RHO%a, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Ambient%RhoCp, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(RHO%w, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Ambient%e, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Ambient%es, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Ambient%Q, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Ambient%VPD, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Ambient%Td, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    !> Dry air properties
    call WriteDatumFloat(Ambient%p_d, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(RHO%d, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Ambient%Vd, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    !> Specific heat of evaporation
    call WriteDatumFloat(Ambient%lambda, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    !> Wet to dry air density ratio
    call WriteDatumFloat(Ambient%sigma, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    !> Water USe Efficiency
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
            call WriteDatumFloat(Stats%d(gas), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumFloat(Stats%r(gas), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumFloat(Stats%chi(gas), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        else
            do i = 1, 4
                call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            end do
        end if
    end do
    !> Timelags (calculated, used, min/max/nominal) for all gases
    !> Gas timelags
    do gas = co2, gas4
        call WriteDatumFloat(Essentials%actual_timelag(gas), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Essentials%used_timelag(gas), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(E2Col(gas)%def_tl, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(E2Col(gas)%min_tl, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(E2Col(gas)%max_tl, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end do

!> Basic stats
    !> Mean values
    do var = u, gas4
        call WriteDatumFloat(Stats6%Mean(var), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end do
    !> 25-50-75%
    do var = u, gas4
        call WriteDatumFloat(Stats6%Median(var), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end do
    do var = u, gas4
        call WriteDatumFloat(Stats6%Q1(var), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end do
    do var = u, gas4
        call WriteDatumFloat(Stats6%Q3(var), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end do
    !> Variances 
    do var = u, gas4
        call WriteDatumFloat(Stats7%Cov(var, var), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end do
    !> Skwenesses
    do var = u, gas4
        call WriteDatumFloat(Stats7%Skw(var), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end do
    !> Kurtosis
    do var = u, gas4
        call WriteDatumFloat(Stats7%Kur(var), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end do
    !> w-covariances 
    call WriteDatumFloat(Stats7%Cov(u, w), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    do var = ts, gas4
        call WriteDatumFloat(Stats7%Cov(w, var), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end do
    !> Gases covariance matrix
    do gas1 = co2, ch4
        do gas2 = gas1 + 1, gas4 
            call WriteDatumFloat(Stats7%Cov(gas1, gas2), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        end do
    end do

!> Footprint
    call WriteDatumFloat(Foot%peak, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Foot%offset, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Foot%x10, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Foot%x30, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Foot%x50, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Foot%x70, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Foot%x90, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)

!> Intermediate results
    !> Fluxes level 0 (uncorrected fluxes)
    call WriteDatumFloat(Essentials%L, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Essentials%zL, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Flux0%tau, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Flux0%H, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Flux0%LE, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Flux0%co2, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Flux0%h2o, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Flux0%ch4, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Flux0%gas4, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    !> Fluxes level 1
    call WriteDatumFloat(Flux1%tau, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Flux1%H, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Flux1%LE, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Flux1%co2, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Flux1%h2o, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Flux1%ch4, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Flux1%gas4, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    !> Fluxes level 2
    call WriteDatumFloat(Flux2%tau, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Flux2%H, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Flux2%LE, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Flux2%co2, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Flux2%h2o, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Flux2%ch4, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Flux2%gas4, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    !> Tin and Tout                 ********************************************

    !> Temperature, pressure and molar volume 
    !> in the cell of closed-paths, for all gases
    !> Cell parameters              **************************************** (make it gas specific like molar volume)
    call WriteDatumFloat(Ambient%Tcell, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Ambient%Pcell, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    !> Molar volume
    do gas = co2, gas4
        call WriteDatumFloat(E2Col(gas)%Va, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end do
    !> Evapotranspiration and sensible heat fluxes in the cell of 
    !> closed-paths (for WPL), with timelags of other gases
    call WriteDatumFloat(Flux3%E_co2, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Flux3%E_ch4, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Flux3%E_gas4, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Flux3%Hi_co2, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Flux3%Hi_h2o, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Flux3%Hi_ch4, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Flux3%Hi_gas4, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    !> Burba Terms 
    call WriteDatumFloat(Burba%h_bot, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Burba%h_top, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Burba%h_spar, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    !> LI-7700 multipliers
    call WriteDatumFloat(Mul7700%A, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Mul7700%B, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Mul7700%C, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    !> WPL Terms                    ********************************************(Individual: H, LE, Pressure)
    !>!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*
    !> Spectral correction factors
    call WriteDatumFloat(BPCF%of(w_u), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(BPCF%of(w_ts), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(BPCF%of(w_h2o), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(BPCF%of(w_co2), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(BPCF%of(w_h2o), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(BPCF%of(w_ch4), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(BPCF%of(w_gas4), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    !> Increasingly filtered w/T covariances (for spectral assessment)
    call WriteDatumFloat(Essentials%degH(NumDegH + 1), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    do j = 1, NumDegH
        call WriteDatumFloat(Essentials%degH(j), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end do

!> QC details
    !>> Number or records eliminated based on custom flags
    call WriteDatumInt(Essentials%m_custom_flags, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    !>> Number or records eliminated based on wind direction filter
    call WriteDatumInt(Essentials%m_wdf, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    !> Summary of data values eliminated based on diagnostics
    !>> Number or records whose anemometric data was eliminated based on Anemometer diagnostics
    call WriteDatumInt(Essentials%m_diag_anem, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    !>> Number or records whose IRGA data was eliminated based on IRGA diagnostics
    do gas = co2, gas4
        call WriteDatumInt(Essentials%m_diag_irga(gas), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end do
    !>> Number of values eliminated by the Spike test
    do var = u, gas4
        call WriteDatumInt(Essentials%m_despiking(var), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end do
    !>> Number of values eliminated by the Absolute Limits test
    do j = u, gas4
        call WriteDatumInt(Essentials%al_s(j), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end do
    !> VM97 Stats used to calculate flags
    !>> Spikes
    do j = u, gas4
        call WriteDatumInt(Essentials%e2spikes(j), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end do
    !>> Amplitude resolution
    do j = u, gas4
        call WriteDatumFloat(Essentials%ar_s(j), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end do
    !>> Dropouts central
    do j = u, gas4
        call WriteDatumFloat(Essentials%do_s_ctr(j), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end do
    !>> Dropouts extremes
    do j = u, gas4
        call WriteDatumFloat(Essentials%do_s_ext(j), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end do
    !>> Higher moments Skewness
    do j = u, gas4
        call WriteDatumFloat(Essentials%sk_s_skw(j), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end do                          
    !>> Higher moments Kurtosis     
    do j = u, gas4
        call WriteDatumFloat(Essentials%sk_s_kur(j), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end do
    !>> AoA
    call WriteDatumFloat(Essentials%aa_s, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    !>> Non-steady wind
    call WriteDatumFloat(Essentials%ns_s_rnv(1), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Essentials%ns_s_rnv(2), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Essentials%ns_s_rns, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
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
    call WriteDatumInt(STDiff%w_u, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(STDiff%w_ts, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(StDiff%w_co2, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(StDiff%w_h2o, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(StDiff%w_ch4, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(StDiff%w_gas4, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(DtDiff%u, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(DtDiff%w, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(DtDiff%ts, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    !> Partial Foken flags
    call WriteDatumInt(STFlg(w_u), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(STFlg(w_ts), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(STFlg(w_co2), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(STFlg(w_h2o), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(STFlg(w_ch4), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(STFlg(w_gas4), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(DTFlg(u), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(DTFlg(w), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(DTFlg(ts), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    !> Final Foken flags
    call WriteDatumInt(QCFlag%tau, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(QCFlag%H, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(QCFlag%h2o, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(QCFlag%co2, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(QCFlag%h2o, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(QCFlag%ch4, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(QCFlag%gas4, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    !> Number of calculated spikes per variables       *************************(Eliminare)
    call WriteDatumInt(Essentials%e2spikes(u), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(Essentials%e2spikes(v), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(Essentials%e2spikes(w), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(Essentials%e2spikes(ts), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    do var = co2, gas4
        call WriteDatumInt(Essentials%e2spikes(var), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end do

    !> LI-7x00 diagnostics breakdown
    if (Diag7200%present) then
        call WriteDatumInt(Diag7200%head_detect, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7200%t_out, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7200%t_in, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7200%aux_in, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7200%delta_p, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7200%chopper, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7200%detector, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7200%pll, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7200%sync, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    else
        do i = 1, 9
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end do
    end if
    if (Diag7500%present) then
        call WriteDatumInt(Diag7500%chopper, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7500%detector, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7500%pll, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7500%sync, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    else
        do i = 1, 4
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end do
    end if
    if (Diag7700%present) then
        call WriteDatumInt(Diag7700%not_ready, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7700%no_signal, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7700%re_unlocked, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7700%bad_temp, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7700%laser_temp_unregulated, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7700%block_temp_unregulated, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7700%motor_spinning, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7700%pump_on, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7700%top_heater_on, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7700%bottom_heater_on, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7700%calibrating, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7700%motor_failure, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7700%bad_aux_tc1, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7700%bad_aux_tc2, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7700%bad_aux_tc3, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7700%box_connected, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    else
        do i = 1, 16
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end do
    end if
    !> AGC/RSSI                     **************************************** (may need to adapt header to whether it's AGC or RSSI for 7200/7500) 
    if(CompareSwVer(E2Col(co2)%instr%sw_ver, SwVerFromString('6.0.0'))) then
        call WriteDatumInt(nint(Essentials%AGC72), datum, EddyProProj%err_label)
    else
        call WriteDatumInt(-nint(Essentials%AGC72), datum, EddyProProj%err_label)
    end if
    call AddDatum(dataline, datum, separator)
    !> LI-7500
    if(CompareSwVer(E2Col(co2)%instr%sw_ver, SwVerFromString('6.0.0'))) then
        call WriteDatumInt(nint(Essentials%AGC75), datum, EddyProProj%err_label)
    else
        call WriteDatumInt(-nint(Essentials%AGC75), datum, EddyProProj%err_label)
    end if
    call AddDatum(dataline, datum, separator)
    !> LI-7700
    call WriteDatumInt(nint(Essentials%RSSI77), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)

!> Processing settings
    !> Whether w-boost calibration was applied
    if (RPsetup%calib_wboost) then
        call WriteDatumInt(1, datum, EddyProProj%err_label)
    else
        call WriteDatumInt(0, datum, EddyProProj%err_label)
    end if
    call AddDatum(dataline, datum, separator)
    !> Whether AoA calibration was applied
    select case(trim(adjustl(RPsetup%calib_aoa)))
        case('automatic')
            call WriteDatumInt(-1, datum, EddyProProj%err_label)
        case('none')
            call WriteDatumInt(0, datum, EddyProProj%err_label)
        case('nakai_06')
            call WriteDatumInt(1, datum, EddyProProj%err_label)
        case('nakai_12')
            call WriteDatumInt(2, datum, EddyProProj%err_label)
    end select
    call AddDatum(dataline, datum, separator)
    !> Tilt compensation method
    select case(trim(adjustl(Meth%rot)))
        case('none')
            call WriteDatumInt(0, datum, EddyProProj%err_label)
        case('double_rotation')
            call WriteDatumInt(1, datum, EddyProProj%err_label)
        case('triple_rotation')
            call WriteDatumInt(2, datum, EddyProProj%err_label)
        case('planar_fit')
            call WriteDatumInt(3, datum, EddyProProj%err_label)
        case('planar_fit_no_bias')
            call WriteDatumInt(4, datum, EddyProProj%err_label)
    end select
    call AddDatum(dataline, datum, separator)
    !> Rotation angles
    call WriteDatumFloat(Essentials%yaw, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Essentials%pitch, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Essentials%roll, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    !> Detrending method and time constant
    select case(trim(adjustl(Meth%det)))
        case('ba')
            call WriteDatumInt(0, datum, EddyProProj%err_label)
        case('ld')
            call WriteDatumInt(1, datum, EddyProProj%err_label)
        case('rm')
            call WriteDatumInt(2, datum, EddyProProj%err_label)
        case('ew')
            call WriteDatumInt(3, datum, EddyProProj%err_label)
    end select
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(RPsetup%Tconst, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    !> Time lag detection method
    select case(trim(adjustl(Meth%tlag)))
        case('none')
            call WriteDatumInt(0, datum, EddyProProj%err_label)
        case('constant')
            call WriteDatumInt(1, datum, EddyProProj%err_label)
        case('maxcov&default')
            call WriteDatumInt(2, datum, EddyProProj%err_label)
        case('maxcov')
            call WriteDatumInt(3, datum, EddyProProj%err_label)
        case('tlag_opt')
            call WriteDatumInt(4, datum, EddyProProj%err_label)
    end select
    call AddDatum(dataline, datum, separator)
    !> WPL terms
    if (EddyProProj%wpl) then
        call WriteDatumInt(1, datum, EddyProProj%err_label)
    else
        call WriteDatumInt(0, datum, EddyProProj%err_label)
    end if
    call AddDatum(dataline, datum, separator)
    !> Burba terms
    if (trim(adjustl(RPSetup%bu_corr)) == 'yes') then
        if (RPSetup%bu_multi) then
            call WriteDatumInt(2, datum, EddyProProj%err_label)
        else
            call WriteDatumInt(1, datum, EddyProProj%err_label)
        end if
    else
        call WriteDatumInt(0, datum, EddyProProj%err_label)
    end if
    call AddDatum(dataline, datum, separator)
    !> Spectral correction method
    select case(trim(adjustl(EddyProProj%hf_meth)))
        case('none')
            call WriteDatumInt(0, datum, EddyProProj%err_label)
        case('moncrieff_97')
            call WriteDatumInt(1, datum, EddyProProj%err_label)
        case('horst_97')
            call WriteDatumInt(2, datum, EddyProProj%err_label)
        case('ibrom_07')
            call WriteDatumInt(3, datum, EddyProProj%err_label)
        case('fratini_12')
            call WriteDatumInt(4, datum, EddyProProj%err_label)
        case('massman_00')
            call WriteDatumInt(5, datum, EddyProProj%err_label)
    end select
    call AddDatum(dataline, datum, separator)
    !> Footprint model
    select case(trim(adjustl(foot_model_used)))
        case('none')
        call AddDatum(dataline, '0', separator)
        case('kljun_04')
        call AddDatum(dataline, '1', separator)
        case('kormann_meixner_01')
        call AddDatum(dataline, '2', separator)
        case('hsieh_00')
        call AddDatum(dataline, '3', separator)
    end select

!> Metadata
    !> Data logger software version
    call WriteDatumInt(Metadata%logger_swver%major, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(Metadata%logger_swver%minor, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(Metadata%logger_swver%revision, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    !> Site location and features
    call WriteDatumFloat(Metadata%lat, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Metadata%lon, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Metadata%alt, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Metadata%canopy_height, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Metadata%d, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Metadata%z0, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    !> Data acquisition settings
    call WriteDatumInt(nint(Metadata%file_length), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(nint(Metadata%ac_freq), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    !> Flux averaging interval
    call WriteDatumInt(RPsetup%avrg_len, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    !> master anemometer
    write(datum, *) E2Col(u)%instr%firm(1:len_trim(E2Col(u)%Instr%firm))
    call AddDatum(dataline, datum, separator)
    write(datum, *) E2Col(u)%Instr%model(1:len_trim(E2Col(u)%Instr%model))
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(E2Col(u)%Instr%height, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    write(datum, *) E2Col(u)%Instr%wformat(1:len_trim(E2Col(u)%Instr%wformat))
    call AddDatum(dataline, datum, separator)
    write(datum, *) E2Col(u)%Instr%wref(1:len_trim(E2Col(u)%Instr%wref))
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(E2Col(u)%Instr%north_offset, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(E2Col(u)%Instr%hpath_length, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(E2Col(u)%Instr%vpath_length, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(E2Col(u)%Instr%tau, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    !> gas analysers details
    do gas = co2, gas4
        write(datum, *) E2Col(gas)%Instr%firm(1:len_trim(E2Col(gas)%Instr%firm))
        call AddDatum(dataline, datum, separator)
        write(datum, *) E2Col(gas)%Instr%model(1:len_trim(E2Col(gas)%Instr%model))
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(E2Col(gas)%Instr%nsep, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(E2Col(gas)%Instr%esep, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(E2Col(gas)%Instr%vsep, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(E2Col(gas)%Instr%tube_l, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(E2Col(gas)%Instr%tube_d, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(E2Col(gas)%Instr%tube_f, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(E2Col(gas)%Instr%kw, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(E2Col(gas)%Instr%ko, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(E2Col(gas)%Instr%hpath_length, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(E2Col(gas)%Instr%vpath_length, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(E2Col(gas)%Instr%tau, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end do

    !> Number and mean values of custom variables
    call WriteDatumInt(NumUserVar, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    if (NumUserVar > 0) then
        do var = 1, NumUserVar
            call WriteDatumFloat(UserStats%Mean(var), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        end do
    end if

    !> All aggregated biomet values in FLUXNET units
    call WriteDatumInt(nbVars, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    do i = 1, nbVars
        call WriteDatumFloat(bAggrFluxnet(i), datum, '-9999.')
        call AddDatum(dataline, datum, separator)
    end do
    write(uicos, '(a)') dataline(1:len_trim(dataline) - 1)


end subroutine WriteIcosOutputRp