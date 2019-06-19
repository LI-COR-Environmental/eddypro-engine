!***************************************************************************
! read_ex_record.f90
! ------------------
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
    integer :: flag
    integer :: gas
    integer :: open_status
    integer :: read_status
    integer :: i
    integer :: var
    integer :: ix
    character(9) :: vm97flags(GHGNumVar)
    character(16000) :: dataline
    real(kind = dbl) :: aux(32)
    include 'interfaces_1.inc'

    ! integer, external :: strCharIndex

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
    if (read_status > 0) then
        ValidRecord = .false.
        if (rec_num > 0) close(unt)
        return
    end if
    if (read_status < 0) then
        EndOfFileReached = .true.
        if (rec_num > 0) close(unt)
        return
    end if

    !> Replace error code with -9999
    dataline = replace2(dataline, trim(EddyProProj%err_label), '-9999')

    !> Read timestamps and eliminate them from dataline
    lEx%start_timestamp = dataline(1:12)
    dataline = dataline(14: len_trim(dataline))
    lEx%end_timestamp = dataline(1:12)
    dataline = dataline(14: len_trim(dataline))
    lEx%end_date = lEx%end_timestamp(1:4) // '-' // lEx%end_timestamp(5:6) // '-' // lEx%end_timestamp(7:8) 
    lEx%end_time = lEx%end_timestamp(9:10) // ':' // lEx%end_timestamp(11:12)  

    !> Extract some data
    read(dataline, *, iostat = read_status) lEx%DOY_start, lEx%DOY_end, lEx%fname, lEx%RP, &
        lEx%nighttime_int, lEx%nr_theor, &
        lEx%nr_files, lEx%nr_after_custom_flags, lEx%nr_after_wdf, &
        lEx%nr(u), lEx%nr(ts:gas4), lEx%nr_w(u), lEx%nr_w(ts:gas4), &
        aux(1:8), & !< Skip final fluxes
        lEx%rand_uncer(u), lEx%rand_uncer(ts), &
        lEx%rand_uncer_LE, lEx%rand_uncer_ET, lEx%rand_uncer(co2:gas4), &
        lEx%Stor%H, lEx%Stor%LE, lEx%Stor%ET, lEx%Stor%of(co2:gas4), &
        aux(1:4), & !< Skip advection fluxes
        lEx%unrot_u, lEx%unrot_v, lEx%unrot_w, lEx%rot_u, lEx%rot_v, lEx%rot_w, &
        lEx%WS, lEx%MWS, lEx%WD, lEx%WD_SIGMA, lEx%ustar, lEx%TKE, lEx%L, lEx%zL, lEx%Bowen, lEx%Tstar, &
        lEx%Ts, lEx%Ta, lEx%Pa, lEx%RH, lEx%Va, lEx%RHO%a, lEx%RhoCp, &
        lEx%RHO%w, lEx%e, lEx%es, lEx%Q, lEx%VPD, lEx%Tdew, &
        lEx%Pd, lEx%RHO%d, lEx%Vd, lEx%lambda, lEx%sigma, &
        lEx%measure_type_int(co2), lEx%d(co2), lEx%r(co2), lEx%chi(co2), &
        lEx%measure_type_int(h2o), lEx%d(h2o), lEx%r(h2o), lEx%chi(h2o), &
        lEx%measure_type_int(ch4), lEx%d(ch4), lEx%r(ch4), lEx%chi(ch4), &
        lEx%measure_type_int(gas4), lEx%d(gas4), lEx%r(gas4), lEx%chi(gas4), &
        lEx%act_tlag(co2), lEx%used_tlag(co2), lEx%nom_tlag(co2), lEx%min_tlag(co2), lEx%max_tlag(co2), &
        lEx%act_tlag(h2o), lEx%used_tlag(h2o), lEx%nom_tlag(h2o), lEx%min_tlag(h2o), lEx%max_tlag(h2o),&
        lEx%act_tlag(ch4), lEx%used_tlag(ch4), lEx%nom_tlag(ch4), lEx%min_tlag(ch4), lEx%max_tlag(ch4),&
        lEx%act_tlag(gas4), lEx%used_tlag(gas4), lEx%nom_tlag(gas4), lEx%min_tlag(gas4), lEx%max_tlag(gas4), &
        lEx%stats%median(u:gas4), lEx%stats%Q1(u:gas4), lEx%stats%Q3(u:gas4), &
        (lEx%stats%Cov(var, var), var=u, gas4), lEx%stats%Skw(u:gas4), lEx%stats%Kur(u:gas4), &
        lEx%stats%Cov(w, u), lEx%stats%Cov(w, ts:gas4), lEx%stats%Cov(co2, h2o:gas4), &
        lEx%stats%Cov(h2o, ch4:gas4), lEx%stats%Cov(ch4, gas4), &
        aux(1:8), & !< Skip footprint
        lEx%Flux0%ustar, lEx%Flux0%L, lEx%Flux0%zL, &
        lEx%Flux0%Tau, lEx%Flux0%H, lEx%Flux0%LE, lEx%Flux0%ET, &
        lEx%Flux0%co2, lEx%Flux0%h2o, lEx%Flux0%ch4, lEx%Flux0%gas4, &
        aux(1:16), & !< Skip fluxes level 1 and 2
        lEx%Tcell, lEx%Pcell, lEx%Vcell(co2:gas4), &
        lEx%Flux0%E_co2, lEx%Flux0%E_ch4, lEx%Flux0%E_gas4, &
        lEx%Flux0%Hi_co2, lEx%Flux0%Hi_h2o, lEx%Flux0%Hi_ch4, lEx%Flux0%Hi_gas4, &
        lEx%Burba%h_bot, lEx%Burba%h_top, lEx%Burba%h_spar, &
        lEx%Mul7700%A, lEx%Mul7700%B, lEx%Mul7700%C, &
        aux(1:8), & !< Skip SCFs
        lEx%degT%cov, lEx%degT%dcov(1:9), &
        lEx%spikes(u:gas4)
    ix = strCharIndex(dataline, ',', 258)
    dataline = dataline(ix+1: len_trim(dataline))
    
    !> Copy NREX chunk
    ix = strCharIndex(dataline, ',', 23)
    fluxnetChunks%s(1) = dataline(1: ix-1)
    dataline = dataline(ix+1: len_trim(dataline))

    !> Read out VM flags and Foken QC details
    read(dataline, *, iostat = read_status) vm97flags(u:GHGNumVar), &
        lEx%vm_tlag_hf, lEx%vm_tlag_sf, lEx%vm_aoa_hf, lEx%vm_nshw_hf 
    ix = strCharIndex(dataline, ',', 12)
    dataline = dataline(ix+1: len_trim(dataline))

    !> Rearrage VM flags per test, instead of per variable
    if (vm97flags(u) == '-9999') then
        lEx%vm_flags = '-9999'
    else
        do flag = 1, 8
            lEx%vm_flags(flag)(1:1) = '8'
            lEx%vm_flags(flag)(2:2) = vm97flags(u)(flag + 1: flag + 1)
            lEx%vm_flags(flag)(3:3) = vm97flags(v)(flag + 1: flag + 1)
            lEx%vm_flags(flag)(4:4) = vm97flags(w)(flag + 1: flag + 1)
            lEx%vm_flags(flag)(5:5) = vm97flags(ts)(flag + 1: flag + 1)
            do gas = co2, gas4
                if (vm97flags(gas)(1:1) == '8') then
                    lEx%vm_flags(flag)(gas + 1 : gas + 1) = vm97flags(gas)(flag + 1: flag + 1)
                else
                    lEx%vm_flags(flag)(gas + 1 : gas + 1) = '9'
                end if
            end do
        end do
    end if

    !> Copy KID/ZCD/CORRDIFF/NSR chunk
    ix = strCharIndex(dataline, ',', 30)
    fluxnetChunks%s(2) = dataline(1: ix-1)
    dataline = dataline(ix+1: len_trim(dataline))

    read(dataline, *, iostat = read_status) &
        lEx%TAU_SS, lEx%H_SS, lEx%FC_SS, lEx%FH2O_SS, &
        lEx%FCH4_SS, lEx%FGS4_SS, lEx%U_ITC, lEx%W_ITC, lEx%TS_ITC
    ix = strCharIndex(dataline, ',', 9)
    dataline = dataline(ix+1: len_trim(dataline))

    !> Copy .._TEST
    ix = strCharIndex(dataline, ',', 17)
    fluxnetChunks%s(3) = dataline(1: ix-1)
    dataline = dataline(ix+1: len_trim(dataline))

    !> Read licor IRGA flags
    read(dataline, *, iostat = read_status) lEx%licor_flags(1:29)
    ix = strCharIndex(dataline, ',', 29)
    dataline = dataline(ix+1: len_trim(dataline))

    !> Read AGC/RSSI
    read(dataline, *, iostat = read_status) lEx%agc72,lEx%agc75,lEx%rssi77
    ix = strCharIndex(dataline, ',', 3)
    dataline = dataline(ix+1: len_trim(dataline))

    !> Copy WBOOST_APPLIED thru AXES_ROTATION_METHOD
    ix = strCharIndex(dataline, ',', 3)
    fluxnetChunks%s(4) = dataline(1: ix-1)
    dataline = dataline(ix+1: len_trim(dataline))

    !> Read rotation angles and detrending method/time constant
    read(dataline, *, iostat = read_status) &
        lEx%yaw, lEx%pitch, lEx%roll, lEx%det_meth_int, lEx%det_timec
    ix = strCharIndex(dataline, ',', 5)
    dataline = dataline(ix+1: len_trim(dataline))

    !> Copy TIMELAG_DETECTION_METHOD thru FOOTPRINT_MODEL
    ix = strCharIndex(dataline, ',', 5)
    fluxnetChunks%s(5) = dataline(1: ix-1)
    dataline = dataline(ix+1: len_trim(dataline))

    !> Read out metadata
    read(dataline, *, iostat = read_status) &
        lEx%logger_swver%major,lEx%logger_swver%minor,lEx%logger_swver%revision, &
        lEx%lat, lEx%lon, lEx%alt, &
        lEx%canopy_height, lEx%disp_height, lEx%rough_length, &
        lEx%file_length, lEx%ac_freq, lEx%avrg_length, &
        lEx%instr(sonic)%firm, lEx%instr(sonic)%model, lEx%instr(sonic)%height, &
        lEx%instr(sonic)%wformat, lEx%instr(sonic)%wref, lEx%instr(sonic)%north_offset, &
        lEx%instr(sonic)%hpath_length, lEx%instr(sonic)%vpath_length, lEx%instr(sonic)%tau, &
        lEx%instr(ico2)%firm, lEx%instr(ico2)%model, lEx%instr(ico2)%nsep, lEx%instr(ico2)%esep, &
        lEx%instr(ico2)%vsep, lEx%instr(ico2)%tube_l, lEx%instr(ico2)%tube_d, &
        lEx%instr(ico2)%tube_f, &
        lEx%instr(ico2)%hpath_length, lEx%instr(ico2)%vpath_length, lEx%instr(ico2)%tau, &
        lEx%instr(ih2o)%firm, lEx%instr(ih2o)%model, lEx%instr(ih2o)%nsep, lEx%instr(ih2o)%esep, &
        lEx%instr(ih2o)%vsep, lEx%instr(ih2o)%tube_l, lEx%instr(ih2o)%tube_d, &
        lEx%instr(ih2o)%tube_f, lEx%instr(ih2o)%kw, lEx%instr(ih2o)%ko, &
        lEx%instr(ih2o)%hpath_length, lEx%instr(ih2o)%vpath_length, lEx%instr(ih2o)%tau, &
        lEx%instr(ich4)%firm, lEx%instr(ich4)%model, lEx%instr(ich4)%nsep, lEx%instr(ich4)%esep, &
        lEx%instr(ich4)%vsep, lEx%instr(ich4)%tube_l, lEx%instr(ich4)%tube_d, &
        lEx%instr(ich4)%tube_f, &
        lEx%instr(ich4)%hpath_length, lEx%instr(ich4)%vpath_length, lEx%instr(ich4)%tau, &
        lEx%instr(igas4)%firm, lEx%instr(igas4)%model, lEx%instr(igas4)%nsep, lEx%instr(igas4)%esep, &
        lEx%instr(igas4)%vsep, lEx%instr(igas4)%tube_l, lEx%instr(igas4)%tube_d, &
        lEx%instr(igas4)%tube_f, &
        lEx%instr(igas4)%hpath_length, lEx%instr(igas4)%vpath_length, lEx%instr(igas4)%tau, &
        lEx%ncustom
    ix = strCharIndex(dataline, ',', 68)
    dataline = dataline(ix+1: len_trim(dataline))

    !> Read custom variables
    if (lEx%ncustom > 0) then
        do i = 1, lEx%ncustom
            read(dataline, *, iostat = read_status) lEx%user_var(i)
            ix = strCharIndex(dataline, ',', 1)
            dataline = dataline(ix+1: len_trim(dataline))
        end do
    end if

    !> Put remaining into last chunk
    fluxnetChunks%s(6) = dataline(1: len_trim(dataline))

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
    integer :: igas
    integer :: gas
    integer :: var

    if (lEx%fname == 'not_enough_data') then
        lEx%not_enough_data = .True.
    else
        lEx%not_enough_data = .False.
    end if

    lEx%var_present = .false.
    if (lEx%WS /= error) lEx%var_present(u:w) = .true.
    if (lEx%Ts /= error) lEx%var_present(ts)  = .true.
    if (lEx%Flux0%co2  /= error) lEx%var_present(co2) = .true.
    if (lEx%Flux0%h2o  /= error) lEx%var_present(h2o) = .true.
    if (lEx%Flux0%ch4  /= error) lEx%var_present(ch4) = .true.
    if (lEx%Flux0%gas4 /= error) lEx%var_present(gas4) = .true.

    !> Units adjustments
    if (lEx%Flux0%ch4 /= error) lEx%Flux0%ch4 = lEx%Flux0%ch4 * 1d-3
    if (lEx%Flux0%gas4 /= error) lEx%Flux0%gas4 = lEx%Flux0%gas4 * 1d-3
    if (lEx%rand_uncer(ch4) /= error) lEx%rand_uncer(ch4) = lEx%rand_uncer(ch4) * 1d-3
    if (lEx%rand_uncer(gas4) /= error) lEx%rand_uncer(gas4) = lEx%rand_uncer(gas4) * 1d-3
    if (lEx%Ts /= error) lEx%Ts = lEx%Ts + 273.15d0
    if (lEx%Ta /= error) lEx%Ta = lEx%Ta + 273.15d0
    if (lEx%Tdew /= error) lEx%Tdew = lEx%Tdew + 273.15d0
    if (lEx%Pa /= error) lEx%Pa = lEx%Pa * 1d3
    if (lEx%Pd /= error) lEx%Pd = lEx%Pd * 1d3
    if (lEx%e /= error) lEx%e = lEx%e * 1d2
    if (lEx%es /= error) lEx%es = lEx%es * 1d2
    if (lEx%Pa /= error) lEx%VPD = lEx%VPD * 1d2
    if (lEx%r(ch4) /= error) lEx%r(ch4) = lEx%r(ch4) * 1d-3
    if (lEx%chi(ch4) /= error) lEx%chi(ch4) = lEx%chi(ch4) * 1d-3
    if (lEx%r(gas4) /= error) lEx%r(gas4) = lEx%r(gas4) * 1d-3
    if (lEx%chi(gas4) /= error) lEx%chi(gas4) = lEx%chi(gas4) * 1d-3
    if (lEx%stats%median(ts) /= error) lEx%stats%median(ts) = lEx%stats%median(ts) + 273.15d0
    if (lEx%stats%Q1(ts) /= error) lEx%stats%Q1(ts) = lEx%stats%Q1(ts) + 273.15d0
    if (lEx%stats%Q3(ts) /= error) lEx%stats%Q3(ts) = lEx%stats%Q3(ts) + 273.15d0
    if (lEx%instr(sonic)%hpath_length /= error) lEx%instr(sonic)%hpath_length = lEx%instr(sonic)%hpath_length * 1d-2
    if (lEx%instr(sonic)%vpath_length /= error) lEx%instr(sonic)%vpath_length = lEx%instr(sonic)%vpath_length * 1d-2
    if (lEx%instr(sonic)%nsep /= error) lEx%instr(sonic)%nsep = lEx%instr(sonic)%nsep * 1d-2
    if (lEx%instr(sonic)%esep /= error) lEx%instr(sonic)%esep = lEx%instr(sonic)%esep * 1d-2
    if (lEx%Burba%h_bot == error) lEx%Burba%h_bot = 0d0
    if (lEx%Burba%h_top == error) lEx%Burba%h_top = 0d0
    if (lEx%Burba%h_spar == error) lEx%Burba%h_spar = 0d0    

    !> Variances were actually read as standard deviations
    do var = u, gas4 
        if (lEx%var_present(var)) &
            lEx%stats%Cov(var, var) = lEx%stats%Cov(var, var)**2
    end do

    lEx%instr(sonic)%category = 'sonic'
    lEx%instr(ico2:igas4)%category = 'irga'
    !> Determine whether igas analysers are open or closed path
    do igas = ico2, igas4
        select case (lEx%instr(igas)%model(1:len_trim(lEx%instr(igas)%model) - 2))
        case ('li7700', 'li7500', 'li7500a', 'li7500rs', 'li7500ds', &
            'generic_open_path', 'open_path_krypton', 'open_path_lyman')
                lEx%instr(igas)%path_type = 'open'
            case default
                lEx%instr(igas)%path_type = 'closed'
                if (lEx%instr(igas)%tube_d /= error) &
                    lEx%instr(igas)%tube_d = lEx%instr(igas)%tube_d * 1d-3
                if (lEx%instr(igas)%tube_l /= error) &
                    lEx%instr(igas)%tube_l = lEx%instr(igas)%tube_l * 1d-2
                if (lEx%instr(igas)%tube_f /= error) &
                    lEx%instr(igas)%tube_f = lEx%instr(igas)%tube_f / 6d4
        end select
        if (lEx%instr(igas)%vsep /= error) lEx%instr(igas)%vsep = lEx%instr(igas)%vsep * 1d-2
        if (lEx%instr(igas)%nsep /= error) lEx%instr(igas)%nsep = lEx%instr(igas)%nsep * 1d-2
        if (lEx%instr(igas)%esep /= error) lEx%instr(igas)%esep = lEx%instr(igas)%esep * 1d-2
        if (lEx%instr(igas)%hpath_length /= error) lEx%instr(igas)%hpath_length = lEx%instr(igas)%hpath_length * 1d-2
        if (lEx%instr(igas)%vpath_length /= error) lEx%instr(igas)%vpath_length = lEx%instr(igas)%vpath_length * 1d-2

        if (lEx%instr(igas)%nsep /= error .and. lEx%instr(igas)%esep /= error) then
            lEx%instr(igas)%hsep = dsqrt(lEx%instr(igas)%nsep**2 + lEx%instr(igas)%esep**2)
        elseif (lEx%instr(igas)%nsep /= error) then
            lEx%instr(igas)%hsep = lEx%instr(igas)%nsep
        elseif (lEx%instr(igas)%esep /= error) then
            lEx%instr(igas)%hsep = lEx%instr(igas)%esep
        end if
    end do

    !> Understand software version (AGC (or RSSI) value is negative)
    !> LI-7200
    if (lEx%agc72 < 0 .and. lEx%agc72 /= error) then
        lEx%agc72 =  - lEx%agc72
    else
        co2_new_sw_ver = .true.
    end if
    !> LI-7500
    if (lEx%agc75 < 0 .and. lEx%agc75 /= error) then
        lEx%agc75 =  - lEx%agc75
    else
        co2_new_sw_ver = .true.
    end if

    !> Detrending method from integers to strings
    select case(lEx%det_meth_int)
        case(0)
            lEx%det_meth = 'ba'
        case(1)
            lEx%det_meth = 'ld'
        case(2)
            lEx%det_meth = 'rm'
        case(3)
            lEx%det_meth = 'ew'
    end select

    !> Measurement type from integers to strings
    do gas = co2, gas4
        select case(lEx%measure_type_int(gas))
            case(0)
                lEx%measure_type(gas) = 'mixing_ratio'
            case(1)
                lEx%measure_type(gas) = 'mole_fraction'
            case(2)
                lEx%measure_type(gas) = 'molar_density'
        end select
    end do

    !> Daytime
    lEx%daytime = lEx%nighttime_int == 0

    !> Legacy values to be later replaced with newer (left-hand sides) *********
    lEx%file_records = lEx%nr_files
    lEx%used_records = lEx%nr_after_wdf
    lEx%tlag = lEx%used_tlag
    lEx%def_tlag = lEx%used_tlag == lEx%nom_tlag
    do var = u, gas4
        lEx%var(var) = lEx%stats%Cov(var, var)
    end do
    lEx%cov_w(u) = lEx%stats%cov(w, u)
    lEx%cov_w(ts:gas4) = lEx%stats%cov(w, ts:gas4)
end subroutine CompleteEssentials
