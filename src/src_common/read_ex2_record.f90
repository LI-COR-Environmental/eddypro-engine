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
subroutine ReadEx2Record(FilePath, unt, rec_num, lEx2, ValidRecord, EndOfFileReached)
    use m_common_global_var
    !> In/out variables
    character(*), intent(in) :: FilePath
    integer, intent(in) :: rec_num
    logical, intent(out) :: ValidRecord
    logical, intent(out) :: EndOfFileReached
    type (Ex2Type), intent(out) :: lEx2
    integer, intent(inout) :: unt
    !> Local variables
    integer :: open_status
    integer :: read_status
    integer :: i
    character(16000) :: dataline


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
    ! if (read_status > 0 .or. index(dataline, 'not_enough_data') /= 0) then
    !     ValidRecord = .false.
    !     if (rec_num > 0) close(unt)
    !     return
    ! end if

    if (read_status < 0) then
        EndOfFileReached = .true.
        if (rec_num > 0) close(unt)
        return
    end if

    !> Strip file name from dataline
    ! lEx2%fname = dataline(1:index(dataline, separator) - 1)
    ! dataline = dataline(index(dataline, separator) + 1: len_trim(dataline))
    !> Read timestamp and eliminate if from dataline
    lEx2%timestamp = dataline(1:16)
                                                          !********************* Calculate here lEx2%date and lEx2%time

    dataline = dataline(18: len_trim(dataline))

    ! !> read rest of results
    ! read(dataline, *, iostat = read_status) lEx2%daytime, lEx2%file_records, lEx2%used_records, &
    !     lEx2%Flux0%Tau, lEx2%rand_uncer(u), lEx2%Flux0%H, lEx2%rand_uncer(ts), &
    !     lEx2%Flux0%LE, lEx2%rand_uncer_LE, lEx2%Flux0%co2, lEx2%rand_uncer(co2), &
    !     lEx2%Flux0%h2o, lEx2%rand_uncer(h2o), lEx2%Flux0%ch4, lEx2%rand_uncer(ch4), &
    !     lEx2%Flux0%gas4, lEx2%rand_uncer(gas4), lEx2%Stor%H, lEx2%Stor%LE, &
    !     lEx2%Stor%of(co2), lEx2%Stor%of(h2o), lEx2%Stor%of(ch4), lEx2%Stor%of(gas4), &
    !     lEx2%Flux0%E_co2, lEx2%Flux0%E_ch4, lEx2%Flux0%E_gas4, &
    !     lEx2%Flux0%Hi_co2, lEx2%Flux0%Hi_h2o, lEx2%Flux0%Hi_ch4, lEx2%Flux0%Hi_gas4, &
    !     lEx2%unrot_u, lEx2%unrot_v, lEx2%unrot_w, lEx2%rot_u, lEx2%rot_v, lEx2%rot_w, &
    !     lEx2%WS, lEx2%MWS, lEx2%WD, lEx2%ustar, lEx2%TKE, lEx2%L, lEx2%zL, lEx2%Bowen, lEx2%Tstar, &
    !     lEx2%measure_type(co2), lEx2%d(co2), lEx2%r(co2), lEx2%chi(co2), &
    !     lEx2%measure_type(h2o), lEx2%d(h2o), lEx2%r(h2o), lEx2%chi(h2o), &
    !     lEx2%measure_type(ch4), lEx2%d(ch4), lEx2%r(ch4), lEx2%chi(ch4), &
    !     lEx2%measure_type(gas4), lEx2%d(gas4), lEx2%r(gas4), lEx2%chi(gas4), &
    !     lEx2%Ts, lEx2%Ta, lEx2%Pa, lEx2%RH, lEx2%Va, lEx2%RHO%a, lEx2%RhoCp, &
    !     lEx2%RHO%w, lEx2%e, lEx2%es, lEx2%Q, lEx2%VPD, lEx2%Tdew, &
    !     lEx2%Pd, lEx2%RHO%d, lEx2%Vd, lEx2%lambda, lEx2%sigma, &
    !     lEx2%Tcell, lEx2%Pcell, lEx2%Vcell(co2), lEx2%Vcell(h2o), lEx2%Vcell(ch4), lEx2%Vcell(gas4), &
    !     lEx2%Mul7700%A, lEx2%Mul7700%B, lEx2%Mul7700%C, &
    !     lEx2%Burba%h_bot, lEx2%Burba%h_top, lEx2%Burba%h_spar, &
    !     lEx2%degT%cov, lEx2%degT%dcov(1:9), &
    !     lEx2%var(u:gas4), lEx2%var(tc), lEx2%var(pc), lEx2%var(te), lEx2%var(pe), &
    !     lEx2%cov_w(u), lEx2%cov_w(v), lEx2%cov_w(ts:gas4), &
    !     lEx2%cov_w(tc), lEx2%cov_w(pc), lEx2%cov_w(te), lEx2%cov_w(pe), &
    !     lEx2%tlag(co2), lEx2%def_tlag(co2), lEx2%tlag(h2o), lEx2%def_tlag(h2o), &
    !     lEx2%tlag(ch4), lEx2%def_tlag(ch4), lEx2%tlag(gas4), lEx2%def_tlag(gas4), &
    !     lEx2%yaw, lEx2%pitch, lEx2%roll, &
    !     lEx2%st_w_u, lEx2%st_w_ts, lEx2%st_w_co2, lEx2%st_w_h2o, lEx2%st_w_ch4, lEx2%st_w_gas4, &
    !     lEx2%dt_u, lEx2%dt_w, lEx2%dt_ts, &
    !     lEx2%det_meth, lEx2%det_timec, &
    !     lEx2%logger_swver%major,lEx2%logger_swver%minor,lEx2%logger_swver%revision, &
    !     lEx2%lat, lEx2%lon, lEx2%alt, lEx2%file_length, &
    !     lEx2%avrg_length, lEx2%ac_freq, &
    !     lEx2%canopy_height, lEx2%disp_height, lEx2%rough_length, &
    !     lEx2%instr(sonic)%firm, lEx2%instr(sonic)%model, lEx2%instr(sonic)%height, &
    !     lEx2%instr(sonic)%wformat, lEx2%instr(sonic)%wref, lEx2%instr(sonic)%north_offset, &
    !     lEx2%instr(sonic)%hpath_length, lEx2%instr(sonic)%vpath_length, lEx2%instr(sonic)%tau, &
    !     lEx2%instr(ico2)%firm, lEx2%instr(ico2)%model, lEx2%instr(ico2)%nsep, lEx2%instr(ico2)%esep, &
    !     lEx2%instr(ico2)%vsep, lEx2%instr(ico2)%tube_l, lEx2%instr(ico2)%tube_d, &
    !     lEx2%instr(ico2)%tube_f, lEx2%instr(ico2)%kw, lEx2%instr(ico2)%ko, &
    !     lEx2%instr(ico2)%hpath_length, lEx2%instr(ico2)%vpath_length, lEx2%instr(ico2)%tau, &
    !     lEx2%instr(ih2o)%firm, lEx2%instr(ih2o)%model, lEx2%instr(ih2o)%nsep, lEx2%instr(ih2o)%esep, &
    !     lEx2%instr(ih2o)%vsep, lEx2%instr(ih2o)%tube_l, lEx2%instr(ih2o)%tube_d, &
    !     lEx2%instr(ih2o)%tube_f, lEx2%instr(ih2o)%kw, lEx2%instr(ih2o)%ko, &
    !     lEx2%instr(ih2o)%hpath_length, lEx2%instr(ih2o)%vpath_length, lEx2%instr(ih2o)%tau, &
    !     lEx2%instr(ich4)%firm, lEx2%instr(ich4)%model, lEx2%instr(ich4)%nsep, lEx2%instr(ich4)%esep, &
    !     lEx2%instr(ich4)%vsep, lEx2%instr(ich4)%tube_l, lEx2%instr(ich4)%tube_d, &
    !     lEx2%instr(ich4)%tube_f, lEx2%instr(ich4)%kw, lEx2%instr(ich4)%ko, &
    !     lEx2%instr(ich4)%hpath_length, lEx2%instr(ich4)%vpath_length, lEx2%instr(ich4)%tau, &
    !     lEx2%instr(igas4)%firm, lEx2%instr(igas4)%model, lEx2%instr(igas4)%nsep, lEx2%instr(igas4)%esep, &
    !     lEx2%instr(igas4)%vsep, lEx2%instr(igas4)%tube_l, lEx2%instr(igas4)%tube_d, &
    !     lEx2%instr(igas4)%tube_f, lEx2%instr(igas4)%kw, lEx2%instr(igas4)%ko, &
    !     lEx2%instr(igas4)%hpath_length, lEx2%instr(igas4)%vpath_length, lEx2%instr(igas4)%tau, &
    !     lEx2%vm_flags(1:12),lEx2%spikes(1:GHGNumVar),lEx2%licor_flags(1:29), &
    !     lEx2%agc72,lEx2%agc75,lEx2%rssi77,NumUserVar
    ! if (read_status /= 0) then
    !     ValidRecord = .false.
    !     if (rec_num > 0) close(unt)
    !     return
    ! end if

    ! !> Now read user variables if they exist
    ! if (NumUserVar > 0) then
    !     !> Reduce dataline to the user variables
    !     do ii = 1, 272
    !         dataline = dataline(index(dataline, ',') + 1: len_trim(dataline))
    !     end do
    !     !> Read user variables
    !     read(dataline, *, iostat = read_status) lEx2%user_var(1:NumUserVar)
    !     if (read_status /= 0) then
    !         ValidRecord = .false.
    !         if (rec_num > 0) close(unt)
    !         return
    !     end if
    ! end if

    ! !> Complete essentials information based on retrieved ones
    ! call CompleteEssentials2(lEx2)

    !> Close file only if it wasn't open on entrance
    if (rec_num > 0) close(unt)
end subroutine ReadEx2Record

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
subroutine CompleteEssentials2(lEx2)
    use m_common_global_var
    implicit none
    !> in/out variables
    type(Ex2Type), intent(inout) :: lEx2
    !> local variables
    integer :: gas


    lEx2%var_present = .false.
    if (lEx2%WS /= error) lEx2%var_present(u:w) = .true.
    if (lEx2%Ts /= error) lEx2%var_present(ts)  = .true.
    if (lEx2%Flux0%co2  /= aflx_error) lEx2%var_present(co2) = .true.
    if (lEx2%Flux0%h2o  /= aflx_error) lEx2%var_present(h2o) = .true.
    if (lEx2%Flux0%ch4  /= aflx_error) lEx2%var_present(ch4) = .true.
    if (lEx2%Flux0%gas4 /= aflx_error) lEx2%var_present(gas4) = .true.

    lEx2%instr(ico2:igas4)%category = 'irga'
    lEx2%instr(sonic)%category = 'sonic'
    !> Determine whether gas analysers are open or closed path
    do gas = ico2, igas4
        select case (lEx2%instr(gas)%model(1:len_trim(lEx2%instr(gas)%model) - 2))
            case ('li7700', 'li7500', 'li7500a', 'li7500rs', 'generic_open_path', &
                'open_path_krypton', 'open_path_lyman')
                lEx2%instr(gas)%path_type = 'open'
            case default
                lEx2%instr(gas)%path_type = 'closed'
        end select
        if (lEx2%instr(gas)%nsep /= error .and. lEx2%instr(gas)%esep /= error) then
            lEx2%instr(gas)%hsep = dsqrt(lEx2%instr(gas)%nsep**2 + lEx2%instr(gas)%esep**2)
        elseif (lEx2%instr(gas)%nsep /= error) then
            lEx2%instr(gas)%hsep = lEx2%instr(gas)%nsep
        elseif (lEx2%instr(gas)%esep /= error) then
            lEx2%instr(gas)%hsep = lEx2%instr(gas)%esep
        end if
    end do

    !> Understand software version (AGC (or RSSI) value is negative)
    !> LI-7200
    if (lEx2%agc72 < 0) then
        lEx2%agc72 =  - lEx2%agc72
    else
        co2_new_sw_ver = .true.
    end if
    !> LI-7500
    if (lEx2%agc75 < 0) then
        lEx2%agc75 =  - lEx2%agc75
    else
        co2_new_sw_ver = .true.
    end if
end subroutine CompleteEssentials2
