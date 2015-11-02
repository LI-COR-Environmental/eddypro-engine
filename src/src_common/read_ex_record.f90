!***************************************************************************
! read_ex_record.f90
! ------------------
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
! \brief       Read one record of essentials file. Based on the requested
!              record number, either reads following record (rec_num < 0)
!              or open the file and look for the actual rec_num
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ReadExRecord(FilePath, unt, rec_num, lEx, ValidRecord, EndOfFileReached)
    use m_common_global_var
    !> In/out variables
    character(*), intent(in) :: FilePath
    integer, intent(in) :: rec_num
    logical, intent(out) :: ValidRecord
    logical, intent(out) :: EndOfFileReached
    type (ExType), intent(out) :: lEx
    integer, intent(inout) :: unt
    !> Local variables
    integer :: open_status
    integer :: read_status
    integer :: i
    character(LongOutstringLen) :: dataline


    !> If rec_num > 0,open file and moves to the requested record
    if (rec_num > 0) then
        open(udf, file = trim(adjustl(FilePath)), status = 'old', iostat = open_status)
        if (open_status /= 0) call ExceptionHandler(60)
        unt = udf
        !> Skip header and all records until the requested one
        do i = 1, rec_num
            read(unt, *)
        end do
    end if

    !> Read data line
    ValidRecord = .true.
    EndOfFileReached = .false.
    read(unt, '(a)', iostat = read_status) dataline

    !> Controls on what was read
    if (read_status > 0 .or. index(dataline, 'not_enough_data') /= 0) then
        ValidRecord = .false.
        if (rec_num > 0) close(unt)
        return
    end if
    if (read_status < 0) then
        EndOfFileReached = .true.
        if (rec_num > 0) close(unt)
        return
    end if

    !> Strip file name from dataline
    lEx%fname = dataline(1:index(dataline, separator) - 1)
    dataline = dataline(index(dataline, separator) + 1: len_trim(dataline))
    !> Read timestamp and eliminate if from dataline
    lEx%date = dataline(1:10)
    lEx%time = dataline(12:16)
    dataline = dataline(18: len_trim(dataline))

    !> read rest of results
    read(dataline, *, iostat = read_status) lEx%daytime, lEx%file_records, lEx%used_records, &
        lEx%Flux0%Tau, lEx%rand_uncer(u), lEx%Flux0%H, lEx%rand_uncer(ts), &
        lEx%Flux0%LE, lEx%rand_uncer_LE, lEx%Flux0%co2, lEx%rand_uncer(co2), &
        lEx%Flux0%h2o, lEx%rand_uncer(h2o), lEx%Flux0%ch4, lEx%rand_uncer(ch4), &
        lEx%Flux0%gas4, lEx%rand_uncer(gas4), lEx%Stor%H, lEx%Stor%LE, &
        lEx%Stor%of(co2), lEx%Stor%of(h2o), lEx%Stor%of(ch4), lEx%Stor%of(gas4), &
        lEx%Flux0%E_co2, lEx%Flux0%E_ch4, lEx%Flux0%E_gas4, &
        lEx%Flux0%Hi_co2, lEx%Flux0%Hi_h2o, lEx%Flux0%Hi_ch4, lEx%Flux0%Hi_gas4, &
        lEx%unrot_u, lEx%unrot_v, lEx%unrot_w, lEx%rot_u, lEx%rot_v, lEx%rot_w, &
        lEx%WS, lEx%MWS, lEx%WD, lEx%ustar, lEx%TKE, lEx%L, lEx%zL, lEx%Bowen, lEx%Tstar, &
        lEx%measure_type(co2), lEx%d(co2), lEx%r(co2), lEx%chi(co2), &
        lEx%measure_type(h2o), lEx%d(h2o), lEx%r(h2o), lEx%chi(h2o), &
        lEx%measure_type(ch4), lEx%d(ch4), lEx%r(ch4), lEx%chi(ch4), &
        lEx%measure_type(gas4), lEx%d(gas4), lEx%r(gas4), lEx%chi(gas4), &
        lEx%Ts, lEx%Ta, lEx%Pa, lEx%RH, lEx%Va, lEx%RHO%a, lEx%RhoCp, &
        lEx%RHO%w, lEx%e, lEx%es, lEx%Q, lEx%VPD, lEx%Tdew, &
        lEx%Pd, lEx%RHO%d, lEx%Vd, lEx%lambda, lEx%sigma, &
        lEx%Tcell, lEx%Pcell, lEx%Vcell(co2), lEx%Vcell(h2o), lEx%Vcell(ch4), lEx%Vcell(gas4), &
        lEx%Mul7700%A, lEx%Mul7700%B, lEx%Mul7700%C, &
        lEx%Burba%h_bot, lEx%Burba%h_top, lEx%Burba%h_spar, &
        lEx%degT%cov, lEx%degT%dcov(1:9), &
        lEx%var(u:gas4), lEx%var(tc), lEx%var(pc), lEx%var(te), lEx%var(pe), &
        lEx%cov_w(u), lEx%cov_w(v), lEx%cov_w(ts:gas4), &
        lEx%cov_w(tc), lEx%cov_w(pc), lEx%cov_w(te), lEx%cov_w(pe), &
        lEx%tlag(co2), lEx%def_tlag(co2), lEx%tlag(h2o), lEx%def_tlag(h2o), &
        lEx%tlag(ch4), lEx%def_tlag(ch4), lEx%tlag(gas4), lEx%def_tlag(gas4), &
        lEx%yaw, lEx%pitch, lEx%roll, &
        lEx%st_w_u, lEx%st_w_ts, lEx%st_w_co2, lEx%st_w_h2o, lEx%st_w_ch4, lEx%st_w_gas4, &
        lEx%dt_u, lEx%dt_w, lEx%dt_ts, &
        lEx%det_meth, lEx%det_timec, &
        lEx%logger_swver%major,lEx%logger_swver%minor,lEx%logger_swver%revision, &
        lEx%lat, lEx%lon, lEx%alt, lEx%file_length, &
        lEx%avrg_length, lEx%ac_freq, &
        lEx%canopy_height, lEx%disp_height, lEx%rough_length, &
        lEx%instr(sonic)%firm, lEx%instr(sonic)%model, lEx%instr(sonic)%height, &
        lEx%instr(sonic)%wformat, lEx%instr(sonic)%wref, lEx%instr(sonic)%north_offset, &
        lEx%instr(sonic)%hpath_length, lEx%instr(sonic)%vpath_length, lEx%instr(sonic)%tau, &
        lEx%instr(ico2)%firm, lEx%instr(ico2)%model, lEx%instr(ico2)%nsep, lEx%instr(ico2)%esep, &
        lEx%instr(ico2)%vsep, lEx%instr(ico2)%tube_l, lEx%instr(ico2)%tube_d, &
        lEx%instr(ico2)%tube_f, lEx%instr(ico2)%kw, lEx%instr(ico2)%ko, &
        lEx%instr(ico2)%hpath_length, lEx%instr(ico2)%vpath_length, lEx%instr(ico2)%tau, &
        lEx%instr(ih2o)%firm, lEx%instr(ih2o)%model, lEx%instr(ih2o)%nsep, lEx%instr(ih2o)%esep, &
        lEx%instr(ih2o)%vsep, lEx%instr(ih2o)%tube_l, lEx%instr(ih2o)%tube_d, &
        lEx%instr(ih2o)%tube_f, lEx%instr(ih2o)%kw, lEx%instr(ih2o)%ko, &
        lEx%instr(ih2o)%hpath_length, lEx%instr(ih2o)%vpath_length, lEx%instr(ih2o)%tau, &
        lEx%instr(ich4)%firm, lEx%instr(ich4)%model, lEx%instr(ich4)%nsep, lEx%instr(ich4)%esep, &
        lEx%instr(ich4)%vsep, lEx%instr(ich4)%tube_l, lEx%instr(ich4)%tube_d, &
        lEx%instr(ich4)%tube_f, lEx%instr(ich4)%kw, lEx%instr(ich4)%ko, &
        lEx%instr(ich4)%hpath_length, lEx%instr(ich4)%vpath_length, lEx%instr(ich4)%tau, &
        lEx%instr(igas4)%firm, lEx%instr(igas4)%model, lEx%instr(igas4)%nsep, lEx%instr(igas4)%esep, &
        lEx%instr(igas4)%vsep, lEx%instr(igas4)%tube_l, lEx%instr(igas4)%tube_d, &
        lEx%instr(igas4)%tube_f, lEx%instr(igas4)%kw, lEx%instr(igas4)%ko, &
        lEx%instr(igas4)%hpath_length, lEx%instr(igas4)%vpath_length, lEx%instr(igas4)%tau, &
        lEx%vm_flags(1:12),lEx%spikes(1:GHGNumVar),lEx%licor_flags(1:29), &
        lEx%agc72,lEx%agc75,lEx%rssi77,NumUserVar
    if (read_status /= 0) then
        ValidRecord = .false.
        if (rec_num > 0) close(unt)
        return
    end if

    !> Now read user variables if they exist
    if (NumUserVar > 0) then
        !> Reduce dataline to the user variables
        do ii = 1, 272
            dataline = dataline(index(dataline, ',') + 1: len_trim(dataline))
        end do
        !> Read user variables
        read(dataline, *, iostat = read_status) lEx%user_var(1:NumUserVar)
        if (read_status /= 0) then
            ValidRecord = .false.
            if (rec_num > 0) close(unt)
            return
        end if
    end if

    !> Complete essentials information based on retrieved ones
    call CompleteEssentials(lEx)

    !> Close file only if it wasn't open on entrance
    if (rec_num > 0) close(unt)
end subroutine ReadExRecord

!***************************************************************************
!
! \brief       Complete essentials information, based on those retrieved \n
!              from the file be useful to other programs
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine CompleteEssentials(lEx)
    use m_common_global_var
    implicit none
    !> in/out variables
    type(ExType), intent(inout) :: lEx
    !> local variables
    integer :: gas


    lEx%var_present = .false.
    if (lEx%WS /= error) lEx%var_present(u:w) = .true.
    if (lEx%Ts /= error) lEx%var_present(ts)  = .true.
    if (lEx%Flux0%co2  /= aflx_error) lEx%var_present(co2) = .true.
    if (lEx%Flux0%h2o  /= aflx_error) lEx%var_present(h2o) = .true.
    if (lEx%Flux0%ch4  /= aflx_error) lEx%var_present(ch4) = .true.
    if (lEx%Flux0%gas4 /= aflx_error) lEx%var_present(gas4) = .true.

    lEx%instr(ico2:igas4)%category = 'irga'
    lEx%instr(sonic)%category = 'sonic'
    !> Determine whether gas analysers are open or closed path
    do gas = ico2, igas4
        select case (lEx%instr(gas)%model(1:len_trim(lEx%instr(gas)%model) - 2))
            case ('li7700', 'li7500', 'li7500a', 'li7500rs', 'generic_open_path', &
                'open_path_krypton', 'open_path_lyman')
                lEx%instr(gas)%path_type = 'open'
            case default
                lEx%instr(gas)%path_type = 'closed'
        end select
        if (lEx%instr(gas)%nsep /= error .and. lEx%instr(gas)%esep /= error) then
            lEx%instr(gas)%hsep = dsqrt(lEx%instr(gas)%nsep**2 + lEx%instr(gas)%esep**2)
        elseif (lEx%instr(gas)%nsep /= error) then
            lEx%instr(gas)%hsep = lEx%instr(gas)%nsep
        elseif (lEx%instr(gas)%esep /= error) then
            lEx%instr(gas)%hsep = lEx%instr(gas)%esep
        end if
    end do

    !> Understand software version (AGC (or RSSI) value is negative)
    !> LI-7200
    if (lEx%agc72 < 0) then
        lEx%agc72 =  - lEx%agc72
    else
        co2_new_sw_ver = .true.
    end if
    !> LI-7500
    if (lEx%agc75 < 0) then
        lEx%agc75 =  - lEx%agc75
    else
        co2_new_sw_ver = .true.
    end if
end subroutine CompleteEssentials
