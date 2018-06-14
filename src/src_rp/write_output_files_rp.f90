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
subroutine WriteOutFiles(init_string, PeriodRecords, PeriodActualRecords, &
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
    if (EddyProProj%out_essentials) then
        call clearstr(dataline)
        !> Preliminary file and timestamp information
        call AddDatum(dataline, init_string(1:index(init_string, ',', .true.) - 1), &
            separator)
        write(datum, *) Stats%daytime
        call AddDatum(dataline, datum, separator)
        write(datum, *) PeriodRecords
        call AddDatum(dataline, datum, separator)
        write(datum, *) PeriodActualRecords
        call AddDatum(dataline, datum, separator)
        !> Uncorrected fluxes
        !> Tau, H, LE
        write(datum, *) Flux0%tau
        call AddDatum(dataline, datum, separator)
        write(datum, *) Essentials%rand_uncer(u)
        call AddDatum(dataline, datum, separator)
        write(datum, *) Flux0%H
        call AddDatum(dataline, datum, separator)
        write(datum, *) Essentials%rand_uncer(ts)
        call AddDatum(dataline, datum, separator)

        if(OutVarPresent(h2o)) then
            write(datum, *) Flux0%LE
        else
            write(datum, *) aflx_error
        end if
        call AddDatum(dataline, datum, separator)
        write(datum, *) Essentials%rand_uncer_LE
        call AddDatum(dataline, datum, separator)

        if(OutVarPresent(co2)) then
            write(datum, *) Flux0%co2
        else
            write(datum, *) aflx_error
        end if
        call AddDatum(dataline, datum, separator)
        write(datum, *) Essentials%rand_uncer(co2)
        call AddDatum(dataline, datum, separator)

        if(OutVarPresent(h2o)) then
            write(datum, *) Flux0%h2o
        else
            write(datum, *) aflx_error
        end if
        call AddDatum(dataline, datum, separator)
        write(datum, *) Essentials%rand_uncer(h2o)
        call AddDatum(dataline, datum, separator)

        if(OutVarPresent(ch4)) then
            write(datum, *) Flux0%ch4
        else
            write(datum, *) aflx_error
        end if
        call AddDatum(dataline, datum, separator)
        write(datum, *) Essentials%rand_uncer(ch4)
        call AddDatum(dataline, datum, separator)

        if(OutVarPresent(gas4)) then
            write(datum, *) Flux0%gas4
        else
            write(datum, *) aflx_error
        end if
        call AddDatum(dataline, datum, separator)
        write(datum, *) Essentials%rand_uncer(gas4)
        call AddDatum(dataline, datum, separator)

        !> Storage fluxes
        write(datum, *) Stor%H
        call AddDatum(dataline, datum, separator)
        write(datum, *) Stor%LE
        call AddDatum(dataline, datum, separator)
        write(datum, *) Stor%of(co2)
        call AddDatum(dataline, datum, separator)
        write(datum, *) Stor%of(h2o)
        call AddDatum(dataline, datum, separator)
        write(datum, *) Stor%of(ch4)
        call AddDatum(dataline, datum, separator)
        write(datum, *) Stor%of(gas4)
        call AddDatum(dataline, datum, separator)

        !> Cell "fluxes" (corrected = uncorrected)
        !> Hint, Eint, with timelags of other gases
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

        !> Wind and turbulence
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

        !> Gas concentrations, densities and timelags
        do gas = co2, gas4
            write(datum, *) E2Col(gas)%measure_type
            call AddDatum(dataline, datum, separator)
            write(datum, *) Stats%d(gas)
            call AddDatum(dataline, datum, separator)
            write(datum, *) Stats%r(gas)
            call AddDatum(dataline, datum, separator)
            write(datum, *) Stats%chi(gas)
            call AddDatum(dataline, datum, separator)
        end do

        !> Air properties
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
        !> Others
        write(datum, *) Ambient%lambda
        call AddDatum(dataline, datum, separator)
        write(datum, *) Ambient%sigma
        call AddDatum(dataline, datum, separator)
        !> Cell parameters
        write(datum, *) Ambient%Tcell
        call AddDatum(dataline, datum, separator)
        write(datum, *) Ambient%Pcell
        call AddDatum(dataline, datum, separator)
        !> Cell molar volume for each gas
        do gas = co2, gas4
            write(datum, *) E2Col(gas)%Va
            call AddDatum(dataline, datum, separator)
        end do

        !> LI-7700 multipliers
        write(datum, *) Mul7700%A
        call AddDatum(dataline, datum, separator)
        write(datum, *) Mul7700%B
        call AddDatum(dataline, datum, separator)
        write(datum, *) Mul7700%C
        call AddDatum(dataline, datum, separator)

        !> Burba terms
        write(datum, *) Burba%h_bot
        call AddDatum(dataline, datum, separator)
        write(datum, *) Burba%h_top
        call AddDatum(dataline, datum, separator)
        write(datum, *) Burba%h_spar
        call AddDatum(dataline, datum, separator)

        !> Undegraded and degraded w/T covariances
        write(datum, *) Essentials%degH(NumDegH + 1)
        call AddDatum(dataline, datum, separator)
        do j = 1, NumDegH
            write(datum, *) Essentials%degH(j)
            call AddDatum(dataline, datum, separator)
        end do

        !> Information to be reported, but not used in Flux Computation
        !> Statistics
        !> Variances
        do j = u, pe
            if (j == ti1 .or. j == ti2) cycle
            write(datum, *) Stats7%Cov(j, j)
            call AddDatum(dataline, datum, separator)
        end do
        !> w-covariances
        write(datum, *) Stats7%Cov(w, u)
        call AddDatum(dataline, datum, separator)
        write(datum, *) Stats7%Cov(w, v)
        call AddDatum(dataline, datum, separator)
        write(datum, *) Stats7%Cov(w, ts)
        call AddDatum(dataline, datum, separator)
        write(datum, *) Stats7%Cov(w, co2)
        call AddDatum(dataline, datum, separator)
        write(datum, *) Stats7%Cov(w, h2o)
        call AddDatum(dataline, datum, separator)
        write(datum, *) Stats7%Cov(w, ch4)
        call AddDatum(dataline, datum, separator)
        write(datum, *) Stats7%Cov(w, gas4)
        call AddDatum(dataline, datum, separator)
        write(datum, *) Stats7%Cov(w, tc)
        call AddDatum(dataline, datum, separator)
        write(datum, *) Stats7%Cov(w, pi)
        call AddDatum(dataline, datum, separator)
        write(datum, *) Stats7%Cov(w, te)
        call AddDatum(dataline, datum, separator)
        write(datum, *) Stats7%Cov(w, pe)
        call AddDatum(dataline, datum, separator)
        !> Gas timelags
        do gas = co2, gas4
            write(datum, *) Essentials%used_timelag(gas)
            call AddDatum(dataline, datum, separator)
            write(datum, *) Essentials%def_tlag(gas)
            call AddDatum(dataline, datum, separator)
        end do
        !> Rotation angles
        write(datum, *) Essentials%yaw
        call AddDatum(dataline, datum, separator)
        write(datum, *) Essentials%pitch
        call AddDatum(dataline, datum, separator)
        write(datum, *) Essentials%roll
        call AddDatum(dataline, datum, separator)
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

        !> detrending method
        write(datum, *) Meth%det
        call AddDatum(dataline, datum, separator)
        write(datum, *) RPsetup%Tconst
        call AddDatum(dataline, datum, separator)

        !> Metadata
        write(datum, *) Metadata%logger_swver%major
        call AddDatum(dataline, datum, separator)
        write(datum, *) Metadata%logger_swver%minor
        call AddDatum(dataline, datum, separator)
        write(datum, *) Metadata%logger_swver%revision
        call AddDatum(dataline, datum, separator)

        write(datum, *) Metadata%lat
        call AddDatum(dataline, datum, separator)
        write(datum, *) Metadata%lon
        call AddDatum(dataline, datum, separator)
        write(datum, *) Metadata%alt
        call AddDatum(dataline, datum, separator)
        write(datum, *) Metadata%file_length
        call AddDatum(dataline, datum, separator)
        write(datum, *) RPsetup%avrg_len
        call AddDatum(dataline, datum, separator)
        write(datum, *) Metadata%ac_freq
        call AddDatum(dataline, datum, separator)
        write(datum, *) Metadata%canopy_height
        call AddDatum(dataline, datum, separator)
        write(datum, *) Metadata%d
        call AddDatum(dataline, datum, separator)
        write(datum, *) Metadata%z0
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
        !> gas analysers information
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

        !> Vickers and Mahrt 97 hard flags
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

        !> Spikes for EddyPro variables
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

        !> LI-COR diagnostics
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

        !> AGCs and RSSI
        !> LI-7200
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

        !> Number and mean values of user variables
        write(datum, *) NumUserVar
        call AddDatum(dataline, datum, separator)
        if (NumUserVar > 0) then
            do var = 1, NumUserVar
                write(datum, *) UserStats%Mean(var)
                call AddDatum(dataline, datum, separator)
            end do
        end if
        write(uex, '(a)') dataline(1:len_trim(dataline) - 1)
    end if

    !>***************************************************************
    !>***************************************************************

    !>Write out full output file (main express output)
    if (EddyProProj%out_full) then
        !> Preliminary file and timestamp information
        call clearstr(dataline)
        call AddDatum(dataline, trim(adjustl(init_string)), separator)
        if (Stats%daytime) then
            call AddDatum(dataline, '1', separator)
        else
            call AddDatum(dataline, '0', separator)
        endif
        call WriteDatumInt(PeriodRecords, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(PeriodActualRecords, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)

        !> Corrected fluxes (Level 3)
        !> Tau
        call WriteDatumFloat(Flux3%tau, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(QCFlag%tau, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        if (RUsetup%meth /= 'none') then
            call WriteDatumFloat(Essentials%rand_uncer(u), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if

        !> H
        call WriteDatumFloat(Flux3%H, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(QCFlag%H, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        if (RUsetup%meth /= 'none') then
            call WriteDatumFloat(Essentials%rand_uncer(ts), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if

        !> LE
        if(OutVarPresent(h2o)) then
            call WriteDatumFloat(Flux3%LE, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumInt(QCFlag%h2o, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            if (RUsetup%meth /= 'none') then
                call WriteDatumFloat(Essentials%rand_uncer_LE, datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
             elseif(EddyProProj%fix_out_format) then
                call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            end if
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if

        !> Gases
        if(OutVarPresent(co2)) then
            call WriteDatumFloat(Flux3%co2, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumInt(QCFlag%co2, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            if (RUsetup%meth /= 'none') then
                call WriteDatumFloat(Essentials%rand_uncer(co2), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
             elseif(EddyProProj%fix_out_format) then
                call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            end if
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if

        if(OutVarPresent(h2o)) then
            call WriteDatumFloat(Flux3%h2o, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumInt(QCFlag%h2o, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            if (RUsetup%meth /= 'none') then
                call WriteDatumFloat(Essentials%rand_uncer(h2o), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
             elseif(EddyProProj%fix_out_format) then
                call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            end if
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if

        if(OutVarPresent(ch4)) then
            call WriteDatumFloat(Flux3%ch4, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumInt(QCFlag%ch4, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            if (RUsetup%meth /= 'none') then
                call WriteDatumFloat(Essentials%rand_uncer(ch4), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
             elseif(EddyProProj%fix_out_format) then
                call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            end if
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if

        if(OutVarPresent(gas4)) then
            call WriteDatumFloat(Flux3%gas4, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumInt(QCFlag%gas4, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            if (RUsetup%meth /= 'none') then
                call WriteDatumFloat(Essentials%rand_uncer(gas4), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
             elseif(EddyProProj%fix_out_format) then
                call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            end if
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if

        !> storage fluxes
        call WriteDatumFloat(Stor%H, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        if(OutVarPresent(h2o)) then
            call WriteDatumFloat(Stor%LE, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if
        do gas = co2, gas4
            if(OutVarPresent(gas)) then
                call WriteDatumFloat(Stor%of(gas), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            elseif(EddyProProj%fix_out_format) then
                call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            end if
        end do

        !> vertical advection fluxes
        do gas = co2, gas4
            if(OutVarPresent(gas)) then
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
            elseif(EddyProProj%fix_out_format) then
                call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            end if
        end do

        !> Gas concentrations, densities and timelags
        do gas = co2, gas4
            if (OutVarPresent(gas)) then
                call WriteDatumFloat(Stats%d(gas), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
                call WriteDatumFloat(Stats%chi(gas), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
                call WriteDatumFloat(Stats%r(gas), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
                call WriteDatumFloat(Essentials%used_timelag(gas), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
                if (Essentials%def_tlag(gas)) then
                    call AddDatum(dataline, '1', separator)
                else
                    call AddDatum(dataline, '0', separator)
                endif
            elseif(EddyProProj%fix_out_format) then
                call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
                call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
                call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
                call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
                call AddDatum(dataline, '9', separator)
            end if
        end do

        !> Air properties
        call WriteDatumFloat(Stats7%Mean(ts), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Ambient%Ta, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Stats%Pr, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(RHO%a, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        if (RHO%a /= 0d0 .and. RHO%a /= error) then
            call WriteDatumFloat(Ambient%RhoCp / RHO%a, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        else
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if
        call WriteDatumFloat(Ambient%Va, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        if (Flux3%h2o /= error) then
            call WriteDatumFloat(Flux3%h2o * h2o_to_ET, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        else
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if
        call WriteDatumFloat(RHO%w, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Ambient%e, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Ambient%es, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Ambient%Q, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Stats%RH, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Ambient%VPD, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Ambient%Td, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)

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
        call WriteDatumFloat(Ambient%WS, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Ambient%MWS, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Stats4%wind_dir, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        !> rotation angles
        call WriteDatumFloat(Essentials%yaw, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Essentials%pitch, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Essentials%roll, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)

        !> turbulence
        call WriteDatumFloat(Ambient%us, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Stats7%TKE, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Ambient%L, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Ambient%zL, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Ambient%bowen, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Ambient%Ts, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)

        !> footprint
        if (Meth%foot /= 'none') then
            select case(foot_model_used(1:len_trim(foot_model_used)))
                case('kljun_04')
                call AddDatum(dataline, '0', separator)
                case('kormann_meixner_01')
                call AddDatum(dataline, '1', separator)
                case('hsieh_00')
                call AddDatum(dataline, '2', separator)
            end select
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
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, '9', separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if

        !> Uncorrected fluxes (Level 0)
        !> Tau
        call WriteDatumFloat(Flux0%tau, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(BPCF%of(w_u), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        !> H
        call WriteDatumFloat(Flux0%H, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(BPCF%of(w_ts), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        !> LE
        if(OutVarPresent(h2o)) then
            call WriteDatumFloat(Flux0%LE, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumFloat(BPCF%of(w_h2o), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if
        !> Gases
        if(OutVarPresent(co2)) then
            call WriteDatumFloat(Flux0%co2, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumFloat(BPCF%of(w_co2), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if
        if(OutVarPresent(h2o)) then
            call WriteDatumFloat(Flux0%h2o, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumFloat(BPCF%of(w_h2o), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if
        if(OutVarPresent(ch4)) then
            call WriteDatumFloat(Flux0%ch4, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumFloat(BPCF%of(w_ch4), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if
        if(OutVarPresent(gas4)) then
            call WriteDatumFloat(Flux0%gas4, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumFloat(BPCF%of(w_gas4), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if

        !> Vickers and Mahrt 97 hard flags
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

        !> Spikes for EddyPro variables
        call WriteDatumInt(Essentials%e2spikes(u), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Essentials%e2spikes(v), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Essentials%e2spikes(w), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Essentials%e2spikes(ts), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        do var = co2, gas4
            if(OutVarPresent(var)) then
                call WriteDatumInt(Essentials%e2spikes(var), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            elseif(EddyProProj%fix_out_format) then
                call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            end if
        end do

        !> LI-COR's diagnostic flags
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
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
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
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
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
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if

        !> AGCs and RSSIs for LI-7200 and LI-7500
        if (Diag7200%present) then
            call WriteDatumInt(nint(Essentials%AGC72), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if
        if (Diag7500%present) then
            call WriteDatumInt(nint(Essentials%AGC75), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if

        !> Variances
        do var = u, ts
            call WriteDatumFloat(Stats%Cov(var, var), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        end do
        do gas = co2, gas4
            if(OutVarPresent(gas)) then
                call WriteDatumFloat(Stats%Cov(gas, gas), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            elseif(EddyProProj%fix_out_format) then
                call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            end if
        end do
        !> w-covariances
        call WriteDatumFloat(Stats%Cov(w, ts), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        do gas = co2, gas4
            if(OutVarPresent(gas)) then
                call WriteDatumFloat(Stats%Cov(w, gas), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            elseif(EddyProProj%fix_out_format) then
                call AddDatum(dataline, &
                trim(adjustl(EddyProProj%err_label)), separator)
            end if
        enddo

        !> Mean values of user variables
        if (NumUserVar > 0) then
            do var = 1, NumUserVar
                call WriteDatumFloat(UserStats%Mean(var), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end do
        end if

        write(uflx, '(a)') dataline(1:len_trim(dataline) - 1)
    end if

    !>==========================================================================
    !>==========================================================================
    !> METADATA file
    !> Site-specific information
    if (EddyProProj%out_md) then
        call clearstr(dataline)
        call AddDatum(dataline, init_string, separator)

        !> Site location and characteristics
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
        !> Acquisition setup
        write(datum, *) Metadata%file_length
        call AddDatum(dataline, datum, separator)
        write(datum, *) Metadata%ac_freq
        call AddDatum(dataline, datum, separator)
        !> Master sonic height and north offset
        write(datum, *) E2Col(u)%instr%firm(1:len_trim(E2Col(u)%Instr%firm))
        call AddDatum(dataline, datum, separator)
        write(datum, *) E2Col(u)%instr%model(1:len_trim(E2Col(u)%Instr%model))
        call AddDatum(dataline, datum, separator)
        write(datum, *) E2Col(u)%instr%height
        call AddDatum(dataline, datum, separator)
        write(datum, *) E2Col(u)%instr%wformat
        call AddDatum(dataline, datum, separator)
        write(datum, *) E2Col(u)%instr%wref
        call AddDatum(dataline, datum, separator)
        write(datum, *) E2Col(u)%instr%north_offset
        call AddDatum(dataline, datum, separator)
        write(datum, *) E2Col(u)%instr%hpath_length * 1d2
        call AddDatum(dataline, datum, separator)
        write(datum, *) E2Col(u)%instr%vpath_length  * 1d2
        call AddDatum(dataline, datum, separator)
        write(datum, *) E2Col(u)%instr%tau
        call AddDatum(dataline, datum, separator)
        !> irgas
        do gas = co2, gas4
            if (OutVarPresent(gas)) then
                write(datum, *) E2Col(gas)%instr%firm(1:len_trim(E2Col(gas)%Instr%firm))
                call AddDatum(dataline, datum, separator)
                write(datum, *) E2Col(gas)%instr%model(1:len_trim(E2Col(gas)%Instr%model))
                call AddDatum(dataline, datum, separator)
                write(datum, *) E2Col(gas)%measure_type
                call AddDatum(dataline, datum, separator)
                write(datum, *) E2Col(gas)%instr%nsep * 1d2
                call AddDatum(dataline, datum, separator)
                write(datum, *) E2Col(gas)%instr%esep * 1d2
                call AddDatum(dataline, datum, separator)
                write(datum, *) E2Col(gas)%instr%vsep * 1d2
                call AddDatum(dataline, datum, separator)
                write(datum, *) E2Col(gas)%instr%tube_l * 1d2
                call AddDatum(dataline, datum, separator)
                write(datum, *) E2Col(gas)%instr%tube_d * 1d3
                call AddDatum(dataline, datum, separator)
                write(datum, *) E2Col(gas)%instr%tube_f * 6d4
                call AddDatum(dataline, datum, separator)
                write(datum, *) E2Col(gas)%instr%kw
                call AddDatum(dataline, datum, separator)
                write(datum, *) E2Col(gas)%instr%ko
                call AddDatum(dataline, datum, separator)
                write(datum, *) E2Col(gas)%instr%hpath_length * 1d2
                call AddDatum(dataline, datum, separator)
                write(datum, *) E2Col(gas)%instr%vpath_length * 1d2
                call AddDatum(dataline, datum, separator)
                write(datum, *) E2Col(gas)%instr%tau
                call AddDatum(dataline, datum, separator)
            end if
        end do
        write(umd, '(a)') dataline(1:len_trim(dataline) - 1)
    end if
end subroutine WriteOutFiles
