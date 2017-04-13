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
subroutine WriteIcosOutputRp(init_string, PeriodRecords, PeriodActualRecords, &
    StDiff, DtDiff)
    use m_rp_global_var
    implicit none
    !> in/out variables
    character(*), intent(in) :: init_string
    integer, intent(in) :: PeriodRecords
    integer, intent(in) :: PeriodActualRecords
    type(QCType), intent(in) :: StDiff
    type(QCType), intent(in) :: DtDiff
    !> local variables
    integer :: var
    integer :: gas
    integer :: j
    integer :: i
!    integer :: prof
    character(LongOutstringLen) :: dataline
    character(DatumLen) :: datum
    character(64) :: tmp_init_string
    character(14) :: iso_basic
    include '../src_common/interfaces.inc'

    !> write Essentials output file (csv) for communication
    !> with Fluxes
    if (EddyProProj%out_fluxnet_eddy) then
        call clearstr(dataline)

                                       !**************************************** (Look at whether to limit to u:gas4 everywhere instead of u:pe somewhere

    !> Timestamp
        tmp_init_string = &
            init_string(index(init_string, ',') +1: &
                        index(init_string, ',', .true.) - 1)
        iso_basic = tmp_init_string(1:4) // tmp_init_string(6:7) &
            // tmp_init_string(9:10) // tmp_init_string(12:13)  &
            // tmp_init_string(15:16) // '00'
        call AddDatum(dataline, trim(adjustl(iso_basic)), separator)

    !> Daytime
        if (Stats%daytime) then
            call AddDatum(dataline, '1', separator)
        else
            call AddDatum(dataline, '0', separator)
        endif

    !> Number of records
        !> Number of records teoretically available for current Averaging Interval
        write(datum, *) PeriodRecords
        call AddDatum(dataline, datum, separator)
        !> Number of records actually available for current Averaging Interval (N_in)
        write(datum, *) PeriodActualRecords
        call AddDatum(dataline, datum, separator)
        !> Number of valid records for anemometric data (N_in – M_diag_anemometer)
        !>!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*
        !> Number of valid records for IRGA data  (N_in – M_diag_IRGA)
        !>!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*
        !> Number of valid records available for each main covariance (w/u, w/ts, w/co2, w/h2o, w/ch4, w/gas4)
        !>!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*

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
        write(datum, *) Essentials%rand_uncer(u)
        call AddDatum(dataline, datum, separator)
        write(datum, *) Essentials%rand_uncer(ts)
        call AddDatum(dataline, datum, separator)
        write(datum, *) Essentials%rand_uncer_LE
        call AddDatum(dataline, datum, separator)
        write(datum, *) Essentials%rand_uncer(co2)
        call AddDatum(dataline, datum, separator)
        write(datum, *) Essentials%rand_uncer(h2o)
        call AddDatum(dataline, datum, separator)
        write(datum, *) Essentials%rand_uncer(ch4)
        call AddDatum(dataline, datum, separator)
        write(datum, *) Essentials%rand_uncer(gas4)
        call AddDatum(dataline, datum, separator)

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
        write(datum, *) Stats4%Mean(u)
        call AddDatum(dataline, datum, separator)
        write(datum, *) Stats4%Mean(v)
        call AddDatum(dataline, datum, separator)
        write(datum, *) Stats4%Mean(w)
        call AddDatum(dataline, datum, separator)
        write(datum, *) Stats5%Mean(u)
        call AddDatum(dataline, datum, separator)
        write(datum, *) Stats5%Mean(v)
        call AddDatum(dataline, datum, separator)
        write(datum, *) Stats5%Mean(w)
        call AddDatum(dataline, datum, separator)
        !> wind speed, wind direction, u*, stability, bowen ratio
        write(datum, *) Ambient%WS
        call AddDatum(dataline, datum, separator)
        write(datum, *) Ambient%MWS
        call AddDatum(dataline, datum, separator)
        write(datum, *) Stats4%wind_dir
        call AddDatum(dataline, datum, separator)
        write(datum, *) Ambient%us
        call AddDatum(dataline, datum, separator)
        write(datum, *) Stats7%TKE
        call AddDatum(dataline, datum, separator)
        write(datum, *) Essentials%L
        call AddDatum(dataline, datum, separator)
        write(datum, *) Essentials%zL
        call AddDatum(dataline, datum, separator)
        write(datum, *) Ambient%bowen
        call AddDatum(dataline, datum, separator)
        write(datum, *) Ambient%Ts
        call AddDatum(dataline, datum, separator)

    !> Termodynamics 
        !> Temperature, pressure, RH, VPD, e, es, etc.
        write(datum, *) Stats7%Mean(ts)
        call AddDatum(dataline, datum, separator)
        write(datum, *) Ambient%Ta
        call AddDatum(dataline, datum, separator)
        write(datum, *) Stats%Pr
        call AddDatum(dataline, datum, separator)
        write(datum, *) Stats%RH
        call AddDatum(dataline, datum, separator)
        write(datum, *) Ambient%Va
        call AddDatum(dataline, datum, separator)
        write(datum, *) RHO%a
        call AddDatum(dataline, datum, separator)
        write(datum, *) Ambient%RhoCp
        call AddDatum(dataline, datum, separator)
        write(datum, *) RHO%w
        call AddDatum(dataline, datum, separator)
        write(datum, *) Ambient%e
        call AddDatum(dataline, datum, separator)
        write(datum, *) Ambient%es
        call AddDatum(dataline, datum, separator)
        write(datum, *) Ambient%Q
        call AddDatum(dataline, datum, separator)
        write(datum, *) Ambient%VPD
        call AddDatum(dataline, datum, separator)
        write(datum, *) Ambient%Td
        call AddDatum(dataline, datum, separator)
        !> Dry air properties
        write(datum, *) Ambient%p_d
        call AddDatum(dataline, datum, separator)
        write(datum, *) RHO%d
        call AddDatum(dataline, datum, separator)
        write(datum, *) Ambient%Vd
        call AddDatum(dataline, datum, separator)
        !> Specific heat of evaporation
        write(datum, *) Ambient%lambda
        call AddDatum(dataline, datum, separator)
        !> Wet to dry air density ratio
        write(datum, *) Ambient%sigma
        call AddDatum(dataline, datum, separator)
        !> Water USe Efficiency
        !>!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*

    !> Gases
        !> Concentrations, densities and "nature" of the raw data (mixing ratio, mole fraction, molar density)
        !> Gas concentrations, densities and timelags
        do gas = co2, gas4
            select case (E2Col(gas)%measure_type)
                case('mixing_ratio')
                    call AddDatum(dataline, '0', separator)
                case('mole_fraction')
                    call AddDatum(dataline, '1', separator)
                case('molar_density')
                    call AddDatum(dataline, '2', separator)
            end select
            write(datum, *) Stats%d(gas)
            call AddDatum(dataline, datum, separator)
            write(datum, *) Stats%r(gas)
            call AddDatum(dataline, datum, separator)
            write(datum, *) Stats%chi(gas)
            call AddDatum(dataline, datum, separator)
        end do
        !> Timelags (calculated, used, min/max/nominal) for all gases
        !> Gas timelags
        do gas = co2, gas4
            write(datum, *) Essentials%actual_timelag(gas)
            call AddDatum(dataline, datum, separator)
            write(datum, *) Essentials%used_timelag(gas)
            call AddDatum(dataline, datum, separator)
            write(datum, *) E2Col(gas)%def_tl
            call AddDatum(dataline, datum, separator)
            write(datum, *) E2Col(gas)%min_tl
            call AddDatum(dataline, datum, separator)
            write(datum, *) E2Col(gas)%max_tl
            call AddDatum(dataline, datum, separator)
        end do

    !> Basic stats
        !> Mean values
        do var = u, pe
            write(datum, *) Stats%Mean(var)
            call AddDatum(dataline, datum, separator)
        end do
        !> 25-50-75%
        do var = u, pe
            write(datum, *) Stats%Median(var)
            call AddDatum(dataline, datum, separator)
        end do
        do var = u, pe
            write(datum, *) Stats%Q1(var)
            call AddDatum(dataline, datum, separator)
        end do
        do var = u, pe
            write(datum, *) Stats%Q3(var)
            call AddDatum(dataline, datum, separator)
        end do
        !> Variances 
        do var = u, pe
            write(datum, *) Stats%Cov(var, var)
            call AddDatum(dataline, datum, separator)
        end do
        !> w-covariances 
        do var = u, pe
            if (var == w) cycle
            write(datum, *) Stats%Cov(w, var)
            call AddDatum(dataline, datum, separator)
        end do
        !> Gases covariance matrix
        do gas1 = co2, ch4
            do gas2 = gas1 + 1, gas4 
                write(datum, *) Stats%Cov(gas1, gas2)
                call AddDatum(dataline, datum, separator)
            end do
        end do
        !> Skwenesses
        do var = u, pe
            write(datum, *) Stats%Skw(var)
            call AddDatum(dataline, datum, separator)
        end do
        !> Kurtosis
        do var = u, pe
            write(datum, *) Stats%Kur(var)
            call AddDatum(dataline, datum, separator)
        end do

    !> Intermediate results
        !> Fluxes level 0 (uncorrected fluxes)
        write(datum, *) Flux0%tau
        call AddDatum(dataline, datum, separator)
        write(datum, *) Flux0%H
        call AddDatum(dataline, datum, separator)
        write(datum, *) Flux0%LE
        call AddDatum(dataline, datum, separator)
        write(datum, *) Flux0%co2
        call AddDatum(dataline, datum, separator)
        write(datum, *) Flux0%h2o
        call AddDatum(dataline, datum, separator)
        write(datum, *) Flux0%ch4
        call AddDatum(dataline, datum, separator)
        write(datum, *) Flux0%gas4
        call AddDatum(dataline, datum, separator)
        !> Fluxes level 1
        write(datum, *) Flux1%tau
        call AddDatum(dataline, datum, separator)
        write(datum, *) Flux1%H
        call AddDatum(dataline, datum, separator)
        write(datum, *) Flux1%LE
        call AddDatum(dataline, datum, separator)
        write(datum, *) Flux1%co2
        call AddDatum(dataline, datum, separator)
        write(datum, *) Flux1%h2o
        call AddDatum(dataline, datum, separator)
        write(datum, *) Flux1%ch4
        call AddDatum(dataline, datum, separator)
        write(datum, *) Flux1%gas4
        call AddDatum(dataline, datum, separator)
        !> Fluxes level 2
        write(datum, *) Flux2%tau
        call AddDatum(dataline, datum, separator)
        write(datum, *) Flux2%H
        call AddDatum(dataline, datum, separator)
        write(datum, *) Flux2%LE
        call AddDatum(dataline, datum, separator)
        write(datum, *) Flux2%co2
        call AddDatum(dataline, datum, separator)
        write(datum, *) Flux2%h2o
        call AddDatum(dataline, datum, separator)
        write(datum, *) Flux2%ch4
        call AddDatum(dataline, datum, separator)
        write(datum, *) Flux2%gas4
        call AddDatum(dataline, datum, separator)
        !> Temperature, pressure and molar volume 
        !> in the cell of closed-paths, for all gases
        !> Cell parameters              **************************************** (make it gas specific like molar volume)
        write(datum, *) Ambient%Tcell
        call AddDatum(dataline, datum, separator)
        write(datum, *) Ambient%Pcell
        call AddDatum(dataline, datum, separator)
        !> Molar volume
        do gas = co2, gas4
            write(datum, *) E2Col(gas)%Va
            call AddDatum(dataline, datum, separator)
        end do
        !> Evapotranspiration and sensible heat fluxes in the cell of 
        !> closed-paths (for WPL), with timelags of other gases
        write(datum, *) Flux3%E_co2
        call AddDatum(dataline, datum, separator)
        write(datum, *) Flux3%E_ch4
        call AddDatum(dataline, datum, separator)
        write(datum, *) Flux3%E_gas4
        call AddDatum(dataline, datum, separator)
        write(datum, *) Flux3%Hi_co2
        call AddDatum(dataline, datum, separator)
        write(datum, *) Flux3%Hi_h2o
        call AddDatum(dataline, datum, separator)
        write(datum, *) Flux3%Hi_ch4
        call AddDatum(dataline, datum, separator)
        write(datum, *) Flux3%Hi_gas4
        call AddDatum(dataline, datum, separator)
        !> Burba Terms 
        write(datum, *) Burba%h_bot
        call AddDatum(dataline, datum, separator)
        write(datum, *) Burba%h_top
        call AddDatum(dataline, datum, separator)
        write(datum, *) Burba%h_spar
        call AddDatum(dataline, datum, separator)
        !> LI-7700 multipliers
        write(datum, *) Mul7700%A
        call AddDatum(dataline, datum, separator)
        write(datum, *) Mul7700%B
        call AddDatum(dataline, datum, separator)
        write(datum, *) Mul7700%C
        call AddDatum(dataline, datum, separator)
        !> WPL Terms
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
        write(datum, *) Essentials%degH(NumDegH + 1)
        call AddDatum(dataline, datum, separator)
        do j = 1, NumDegH
            write(datum, *) Essentials%degH(j)
            call AddDatum(dataline, datum, separator)
        end do

    !> QC details
        !> Summary of data values/records eliminated based on diagnostics
        !> Number or records whose IRGA data was eliminated based on IRGA diagnostics (M_diag_IRGA)
        !>!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*
        !> Number or records whose anemometric data was eliminated based on Anemometer diagnostics (M_diag_anemometer)
        !>!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*
        !> Summary of data values/records eliminated based on other filters:
        !> Number or records whose anemometric data was eliminated based on wind direction filter (M_wind_dir)
        !>!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*
        !> Number of values eliminated due to spike test or absolute limits test, by variable (M_spikes_u, M_spikes_v, …, M_abslim_u, M_abslim_v, …) 
        !>!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*
        !> VM97 Stats used to calculate flags
        !>> Spikes
        do j = u, gas4
            call WriteDatumFloat(Essentials%e2spikes(j), datum, EddyProProj%err_label)
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
        !>> Absolute limits             **************************************** (may be done by couting how many vals are outside threshold, per variable)
        !>!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*
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
        !>> Discontinuites
        do j = u, gas4
            do i = 1, 6
                call WriteDatumFloat(Essentials%ds_s_haar_avg(i, j), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
                call WriteDatumFloat(Essentials%ds_s_haar_var(i, j), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end do
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
        !> Foken stats used to calculate flags
        !> Quality test results
        write(datum, *) STDiff%w_u
        call AddDatum(dataline, datum, separator)
        write(datum, *) STDiff%w_ts
        call AddDatum(dataline, datum, separator)
        write(datum, *) StDiff%w_co2
        call AddDatum(dataline, datum, separator)
        write(datum, *) StDiff%w_h2o
        call AddDatum(dataline, datum, separator)
        write(datum, *) StDiff%w_ch4
        call AddDatum(dataline, datum, separator)
        write(datum, *) StDiff%w_gas4
        call AddDatum(dataline, datum, separator)
        write(datum, *) DtDiff%u
        call AddDatum(dataline, datum, separator)
        write(datum, *) DtDiff%w
        call AddDatum(dataline, datum, separator)
        write(datum, *) DtDiff%ts
        call AddDatum(dataline, datum, separator)
        !> Foken flags
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
        !> Number of calculated spikes per variables
        write(datum, *) Essentials%e2spikes(u)
        call AddDatum(dataline, datum, separator)
        write(datum, *) Essentials%e2spikes(v)
        call AddDatum(dataline, datum, separator)
        write(datum, *) Essentials%e2spikes(w)
        call AddDatum(dataline, datum, separator)
        write(datum, *) Essentials%e2spikes(ts)
        call AddDatum(dataline, datum, separator)
        do var = co2, gas4
            write(datum, *) Essentials%e2spikes(var)
            call AddDatum(dataline, datum, separator)
        end do
        !> LI-7x00 diagnostics breakdown
        if (Diag7200%present) then
            write(datum, *) Diag7200%head_detect
            call AddDatum(dataline, datum, separator)
            write(datum, *) Diag7200%t_out
            call AddDatum(dataline, datum, separator)
            write(datum, *) Diag7200%t_in
            call AddDatum(dataline, datum, separator)
            write(datum, *) Diag7200%aux_in
            call AddDatum(dataline, datum, separator)
            write(datum, *) Diag7200%delta_p
            call AddDatum(dataline, datum, separator)
            write(datum, *) Diag7200%chopper
            call AddDatum(dataline, datum, separator)
            write(datum, *) Diag7200%detector
            call AddDatum(dataline, datum, separator)
            write(datum, *) Diag7200%pll
            call AddDatum(dataline, datum, separator)
            write(datum, *) Diag7200%sync
            call AddDatum(dataline, datum, separator)
        else
            do i = 1, 9
                write(datum, *) nint(error)
                call AddDatum(dataline, datum, separator)
            end do
        end if
        if (Diag7500%present) then
            write(datum, *) Diag7500%chopper
            call AddDatum(dataline, datum, separator)
            write(datum, *) Diag7500%detector
            call AddDatum(dataline, datum, separator)
            write(datum, *) Diag7500%pll
            call AddDatum(dataline, datum, separator)
            write(datum, *) Diag7500%sync
            call AddDatum(dataline, datum, separator)
        else
            do i = 1, 4
                write(datum, *) nint(error)
                call AddDatum(dataline, datum, separator)
            end do
        end if
        if (Diag7700%present) then
            write(datum, *) Diag7700%not_ready
            call AddDatum(dataline, datum, separator)
            write(datum, *) Diag7700%no_signal
            call AddDatum(dataline, datum, separator)
            write(datum, *) Diag7700%re_unlocked
            call AddDatum(dataline, datum, separator)
            write(datum, *) Diag7700%bad_temp
            call AddDatum(dataline, datum, separator)
            write(datum, *) Diag7700%laser_temp_unregulated
            call AddDatum(dataline, datum, separator)
            write(datum, *) Diag7700%block_temp_unregulated
            call AddDatum(dataline, datum, separator)
            write(datum, *) Diag7700%motor_spinning
            call AddDatum(dataline, datum, separator)
            write(datum, *) Diag7700%pump_on
            call AddDatum(dataline, datum, separator)
            write(datum, *) Diag7700%top_heater_on
            call AddDatum(dataline, datum, separator)
            write(datum, *) Diag7700%bottom_heater_on
            call AddDatum(dataline, datum, separator)
            write(datum, *) Diag7700%calibrating
            call AddDatum(dataline, datum, separator)
            write(datum, *) Diag7700%motor_failure
            call AddDatum(dataline, datum, separator)
            write(datum, *) Diag7700%bad_aux_tc1
            call AddDatum(dataline, datum, separator)
            write(datum, *) Diag7700%bad_aux_tc2
            call AddDatum(dataline, datum, separator)
            write(datum, *) Diag7700%bad_aux_tc3
            call AddDatum(dataline, datum, separator)
            write(datum, *) Diag7700%box_connected
            call AddDatum(dataline, datum, separator)
        else
            do i = 1, 16
                write(datum, *) error
                call AddDatum(dataline, datum, separator)
            end do
        end if
        !> AGC/RSSI                     **************************************** (may need to adapt header to whether it's AGC or RSSI for 7200/7500) 
        if(CompareSwVer(E2Col(co2)%instr%sw_ver, SwVerFromString('6.0.0'))) then
            write(datum, *)   nint(Essentials%AGC72)
        else
            write(datum, *) - nint(Essentials%AGC72)
        end if
        call AddDatum(dataline, datum, separator)
        !> LI-7500
        if(CompareSwVer(E2Col(co2)%instr%sw_ver, SwVerFromString('6.0.0'))) then
            write(datum, *)   nint(Essentials%AGC75)
        else
            write(datum, *) - nint(Essentials%AGC75)
        end if
        call AddDatum(dataline, datum, separator)
        !> LI-7700
        write(datum, *) nint(Essentials%RSSI77)
        call AddDatum(dataline, datum, separator)

    !> Processing settings
        !> Rotation angles
        write(datum, *) Essentials%yaw
        call AddDatum(dataline, datum, separator)
        write(datum, *) Essentials%pitch
        call AddDatum(dataline, datum, separator)
        write(datum, *) Essentials%roll
        call AddDatum(dataline, datum, separator)
        !> Detrending method and time constant
        write(datum, *) Meth%det
        call AddDatum(dataline, datum, separator)
        write(datum, *) RPsetup%Tconst
        call AddDatum(dataline, datum, separator)

    !> Metadata
        !> Data logger software version
        write(datum, *) Metadata%logger_swver%major
        call AddDatum(dataline, datum, separator)
        write(datum, *) Metadata%logger_swver%minor
        call AddDatum(dataline, datum, separator)
        write(datum, *) Metadata%logger_swver%revision
        call AddDatum(dataline, datum, separator)
        !> Site location and features
        write(datum, *) Metadata%lat
        call AddDatum(dataline, datum, separator)
        write(datum, *) Metadata%lon
        call AddDatum(dataline, datum, separator)
        write(datum, *) Metadata%alt
        call AddDatum(dataline, datum, separator)
        write(datum, *) Metadata%canopy_height
        call AddDatum(dataline, datum, separator)
        write(datum, *) Metadata%d
        call AddDatum(dataline, datum, separator)
        write(datum, *) Metadata%z0
        call AddDatum(dataline, datum, separator)
        !> Data acquisition settings
        write(datum, *) Metadata%file_length
        call AddDatum(dataline, datum, separator)
        write(datum, *) Metadata%ac_freq
        call AddDatum(dataline, datum, separator)
        !> Flux averaging interval
        write(datum, *) RPsetup%avrg_len
        call AddDatum(dataline, datum, separator)
        !> master anemometer
        write(datum, *) E2Col(u)%instr%firm(1:len_trim(E2Col(u)%Instr%firm))
        call AddDatum(dataline, datum, separator)
        write(datum, *) E2Col(u)%Instr%model(1:len_trim(E2Col(u)%Instr%model))
        call AddDatum(dataline, datum, separator)
        write(datum, *) E2Col(u)%Instr%height
        call AddDatum(dataline, datum, separator)
        write(datum, *) E2Col(u)%Instr%wformat(1:len_trim(E2Col(u)%Instr%wformat))
        call AddDatum(dataline, datum, separator)
        write(datum, *) E2Col(u)%Instr%wref(1:len_trim(E2Col(u)%Instr%wref))
        call AddDatum(dataline, datum, separator)
        write(datum, *) E2Col(u)%Instr%north_offset
        call AddDatum(dataline, datum, separator)
        write(datum, *) E2Col(u)%Instr%hpath_length
        call AddDatum(dataline, datum, separator)
        write(datum, *) E2Col(u)%Instr%vpath_length
        call AddDatum(dataline, datum, separator)
        write(datum, *) E2Col(u)%Instr%tau
        call AddDatum(dataline, datum, separator)
        !> gas analysers details
        do gas = co2, gas4
            write(datum, *) E2Col(gas)%Instr%firm(1:len_trim(E2Col(gas)%Instr%firm))
            call AddDatum(dataline, datum, separator)
            write(datum, *) E2Col(gas)%Instr%model(1:len_trim(E2Col(gas)%Instr%model))
            call AddDatum(dataline, datum, separator)
            write(datum, *) E2Col(gas)%Instr%nsep
            call AddDatum(dataline, datum, separator)
            write(datum, *) E2Col(gas)%Instr%esep
            call AddDatum(dataline, datum, separator)
            write(datum, *) E2Col(gas)%Instr%vsep
            call AddDatum(dataline, datum, separator)
            write(datum, *) E2Col(gas)%Instr%tube_l
            call AddDatum(dataline, datum, separator)
            write(datum, *) E2Col(gas)%Instr%tube_d
            call AddDatum(dataline, datum, separator)
            write(datum, *) E2Col(gas)%Instr%tube_f
            call AddDatum(dataline, datum, separator)
            write(datum, *) E2Col(gas)%Instr%kw
            call AddDatum(dataline, datum, separator)
            write(datum, *) E2Col(gas)%Instr%ko
            call AddDatum(dataline, datum, separator)
            write(datum, *) E2Col(gas)%Instr%hpath_length
            call AddDatum(dataline, datum, separator)
            write(datum, *) E2Col(gas)%Instr%vpath_length
            call AddDatum(dataline, datum, separator)
            write(datum, *) E2Col(gas)%Instr%tau
            call AddDatum(dataline, datum, separator)
        end do

    !> Custom variables
        !> Number and mean values of custom variables
        write(datum, *) NumUserVar
        call AddDatum(dataline, datum, separator)
        if (NumUserVar > 0) then
            do var = 1, NumUserVar
                write(datum, *) UserStats%Mean(var)
                call AddDatum(dataline, datum, separator)
            end do
        end if

        write(uicos, '(a)') dataline(1:len_trim(dataline) - 1)
    end if


end subroutine WriteIcosOutputRp