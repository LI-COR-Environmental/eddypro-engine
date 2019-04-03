!***************************************************************************
! init_outfiles_rp.f90
! --------------------
! Copyright (C) 2007-2011, Eco2s team, Gerardo Fratini
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
! \brief       Initializes EddyPro output files
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine InitOutFiles_rp()
    use m_rp_global_var
    use iso_fortran_env
    implicit none
    !> in/out variables
    integer, external :: CreateDir
    !> local variables
    integer :: mkdir_status = 1      ! initializing to false
    integer :: open_status = 1      ! initializing to false
    integer :: dot
    integer :: var
    integer :: gas
    integer :: i
    integer :: j
    character(PathLen) :: Test_Path
    character(64) :: e2sg(E2NumVar)
    character(32) :: usg(NumUserVar)
    character(LongOutstringLen) :: header1
    character(LongOutstringLen) :: header2
    character(LongOutstringLen) :: header3
    character(LongOutstringLen) :: head1_utf8
    character(LongOutstringLen) :: head2_utf8
    character(LongOutstringLen) :: head3_utf8
    character(LongOutstringLen) :: dataline
    logical :: proceed
    include '../src_common/interfaces.inc'


    !> Convenient strings
    e2sg(u)   = 'u_'
    e2sg(v)   = 'v_'
    e2sg(w)   = 'w_'
    e2sg(ts)  = 'ts_'
    e2sg(co2) = 'co2_'
    e2sg(h2o) = 'h2o_'
    e2sg(ch4) = 'ch4_'
    e2sg(gas4) = E2Col(gas4)%label(1:len_trim(E2Col(gas4)%label)) // '_'
    e2sg(tc)  = 'cell_t_'
    e2sg(ti1) = 'inlet_t_'
    e2sg(ti2) = 'outlet_t_'
    e2sg(pi)  = 'cell_p_'
    e2sg(te)  = 'air_t_'
    e2sg(pe)  = 'air_p_'

    call lowercase(e2sg(gas4))
    
    do j = 1, NumUserVar
        usg(j)  = UserCol(j)%label(1:len_trim(UserCol(j)%label)) // '_'
        call lowercase(usg(j))
    end do

    !> Create sub-directory
    !> Stats dir
    proceed = .false.
    do i = 1, 7
        if (RPsetup%out_st(i)) then
            proceed = .true.
            exit
        end if
    end do
    if (proceed) then
        StatsDir = Dir%main_out(1:len_trim(Dir%main_out)) // SubDirStats // slash
        mkdir_status = CreateDir('"' // StatsDir(1:len_trim(StatsDir)) // '"')
    end if
    !> Raw dataset dir
    proceed = .false.
    do i = 1, 7
        if (RPsetup%out_raw(i)) then
            proceed = .true.
            exit
        end if
    end do
    if (proceed) then
        RawDir = Dir%main_out(1:len_trim(Dir%main_out)) // SubDirRaw // slash
        mkdir_status = CreateDir('"' // RawDir(1:len_trim(RawDir)) // '"')
        !> Create subfolders for selected outputs
        if (RPsetup%out_raw(1)) then
            RawSubDir(1) = RawDir(1:len_trim(RawDir)) // 'level_1' // slash
            mkdir_status = CreateDir('"' // RawSubDir(1)(1:len_trim(RawSubDir(1))) // '"')
        end if
        if (RPsetup%out_raw(2)) then
            RawSubDir(2) = RawDir(1:len_trim(RawDir)) // 'level_2' // slash
            mkdir_status = CreateDir('"' // RawSubDir(2)(1:len_trim(RawSubDir(2))) // '"')
        end if
        if (RPsetup%out_raw(3)) then
            RawSubDir(3) = RawDir(1:len_trim(RawDir)) // 'level_3' // slash
            mkdir_status = CreateDir('"' // RawSubDir(3)(1:len_trim(RawSubDir(3))) // '"')
        end if
        if (RPsetup%out_raw(4)) then
            RawSubDir(4) = RawDir(1:len_trim(RawDir)) // 'level_4' // slash
            mkdir_status = CreateDir('"' // RawSubDir(4)(1:len_trim(RawSubDir(4))) // '"')
        end if
        if (RPsetup%out_raw(5)) then
            RawSubDir(5) = RawDir(1:len_trim(RawDir)) // 'level_5' // slash
            mkdir_status = CreateDir('"' // RawSubDir(5)(1:len_trim(RawSubDir(5))) // '"')
        end if
        if (RPsetup%out_raw(6)) then
            RawSubDir(6) = RawDir(1:len_trim(RawDir)) // 'level_6' // slash
            mkdir_status = CreateDir('"' // RawSubDir(6)(1:len_trim(RawSubDir(6))) // '"')
        end if
        if (RPsetup%out_raw(7)) then
            RawSubDir(7) = RawDir(1:len_trim(RawDir)) // 'level_7' // slash
            mkdir_status = CreateDir('"' // RawSubDir(7)(1:len_trim(RawSubDir(7))) // '"')
        end if
    end if

    !> Binned cospectral dir
    if (RPsetup%out_bin_sp) then
        BinCospectraDir = Dir%main_out(1:len_trim(Dir%main_out)) // SubDirBinCospectra // slash
        mkdir_status = CreateDir('"' // BinCospectraDir(1:len_trim(BinCospectraDir)) // '"')
    end if
    !> Binned ogive dir
    if (RPsetup%out_bin_og) then
        BinOgivesDir = Dir%main_out(1:len_trim(Dir%main_out)) // SubDirBinOgives // slash
        mkdir_status = CreateDir('"' // BinOgivesDir(1:len_trim(BinOgivesDir)) // '"')
    end if
    !> Full cospectra dir
    !> (First determine if at least one full (co)spectrum has to be written on output)
    proceed = .false.
    do var = 1, GHGNumVar
        if (RPsetup%out_full_sp(var) .or. RPsetup%out_full_cosp(var)) then
            proceed = .true.
            exit
        end if
    end do
    if (proceed) then
        CospectraDir = Dir%main_out(1:len_trim(Dir%main_out)) // SubDirCospectra // slash
        mkdir_status = CreateDir('"' // CospectraDir(1:len_trim(CospectraDir)) // '"')
    end if

    !> Open full output file and writes header
    if (EddyProProj%out_full) then
        Test_Path = Dir%main_out(1:len_trim(Dir%main_out)) &
                  // EddyProProj%id(1:len_trim(EddyProProj%id)) &
                  // FullOut_FilePadding // Timestamp_FilePadding // CsvExt
        dot = index(Test_Path, CsvExt, .true.) - 1
        FullOut_Path = Test_Path(1:dot) // CsvTmpExt
        open(uflx, file = FullOut_Path, iostat = open_status, encoding = 'utf-8')

        !> Initialize header strings to void
        call Clearstr(header1)
        call Clearstr(header2)
        call Clearstr(header3)
        call Clearstr(head1_utf8)
        call Clearstr(head2_utf8)
        call Clearstr(head3_utf8)

        if (.not. EddyProProj%fix_out_format) then
            !> Initial file and timestamp info
            call AddDatum(header1,'file_info,,,,,,', separator)
            call AddDatum(header2,'filename,date,time,DOY,daytime,file_records,used_records', separator)
            call AddDatum(header3,',[yyyy-mm-dd],[HH:MM],[ddd.ddd],[1=daytime],[#],[#]', separator)

            !> Corrected fluxes (Level 3) and quality flags
            !> Tau
            call AddDatum(header1, 'corrected_fluxes_and_quality_flags,', separator)
            call AddDatum(header2,'Tau,qc_Tau', separator)
            call AddDatum(header3,'[kg+1m-1s-2],[#]', separator)
            if (RUsetup%meth /= 'none') then
                call AddDatum(header1, '', separator)
                call AddDatum(header2,'rand_err_Tau', separator)
                call AddDatum(header3,'[kg+1m-1s-2]', separator)
            end if

            !> H
            call AddDatum(header1, ',', separator)
            call AddDatum(header2, 'H,qc_H', separator)
            call AddDatum(header3, '[W+1m-2],[#]', separator)
            if (RUsetup%meth /= 'none') then
                call AddDatum(header1, '', separator)
                call AddDatum(header2, 'rand_err_H', separator)
                call AddDatum(header3, '[W+1m-2]', separator)
            end if

            !> LE
            if(OutVarPresent(h2o)) then
                call AddDatum(header1, ',', separator)
                call AddDatum(header2, 'LE,qc_LE', separator)
                call AddDatum(header3, '[W+1m-2],[#]', separator)
                if (RUsetup%meth /= 'none') then
                    call AddDatum(header1, '', separator)
                    call AddDatum(header2, 'rand_err_LE', separator)
                    call AddDatum(header3, '[W+1m-2]', separator)
                end if
            end if

            !> Corrected co2 fluxes
            if(OutVarPresent(co2)) then
                call AddDatum(header1, ',', separator)
                call AddDatum(header2, 'co2_flux,qc_co2_flux', separator)
                call AddDatum(header3, '[' // char(181) // 'mol+1s-1m-2],[#]', separator)
                if (RUsetup%meth /= 'none') then
                    call AddDatum(header1, '', separator)
                    call AddDatum(header2, 'rand_err_co2_flux', separator)
                    call AddDatum(header3, '[' // char(181) // 'mol+1s-1m-2]', separator)
                end if
            end if

            !> Corrected h2o fluxes
            if(OutVarPresent(h2o)) then
                call AddDatum(header1, ',', separator)
                call AddDatum(header2,'h2o_flux,qc_h2o_flux', separator)
                call AddDatum(header3,'[mmol+1s-1m-2],[#]', separator)
                if (RUsetup%meth /= 'none') then
                    call AddDatum(header1, '', separator)
                    call AddDatum(header2, 'rand_err_h2o_flux', separator)
                    call AddDatum(header3, '[mmol+1s-1m-2]', separator)
                end if
            end if

            !> Corrected ch4 fluxes
            if(OutVarPresent(ch4)) then
                call AddDatum(header1, ',', separator)
                call AddDatum(header2,'ch4_flux,qc_ch4_flux', separator)
                call AddDatum(header3, '[' // char(181) // 'mol+1s-1m-2],[#]', separator)
                if (RUsetup%meth /= 'none') then
                    call AddDatum(header1, '', separator)
                    call AddDatum(header2, 'rand_err_ch4_flux', separator)
                    call AddDatum(header3, '[' // char(181) // 'mol+1s-1m-2]', separator)
                end if
            end if

            !> Corrected 4th gas fluxes
            if(OutVarPresent(gas4)) then
                call AddDatum(header1, ',', separator)
                call AddDatum(header2, e2sg(gas4)(1:len_trim(e2sg(gas4))) &
                    // 'flux,qc_' // e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'flux', separator)
                call AddDatum(header3, '[' // char(181) // 'mol+1s-1m-2],[#]', separator)
                if (RUsetup%meth /= 'none') then
                    call AddDatum(header1, '', separator)
                    call AddDatum(header2, 'rand_err_' // e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'flux', separator)
                    call AddDatum(header3, '[' // char(181) // 'mol+1s-1m-2]', separator)
                end if
            end if

            !> Storage
            call AddDatum(header1, 'storage_fluxes', separator)
            call AddDatum(header2,'H_strg', separator)
            call AddDatum(header3,'[W+1m-2]', separator)
            if(OutVarPresent(h2o)) call AddDatum(header1, '', separator)
            if(OutVarPresent(h2o)) call AddDatum(header2,'LE_strg', separator)
            if(OutVarPresent(h2o)) call AddDatum(header3,'[W+1m-2]', separator)
            do gas = co2, gas4
                if (gas /= h2o) then
                    if(OutVarPresent(gas)) call AddDatum(header1, '', separator)
                    if(OutVarPresent(gas)) call AddDatum(header2, e2sg(gas)(1:len_trim(e2sg(gas))) // 'strg', separator)
                    if(OutVarPresent(gas)) call AddDatum(header3, '[' // char(181) // 'mol+1s-1m-2]', separator)
                else
                    if(OutVarPresent(gas)) call AddDatum(header1, '', separator)
                    if(OutVarPresent(gas)) call AddDatum(header2, e2sg(gas)(1:len_trim(e2sg(gas))) // 'strg', separator)
                    if(OutVarPresent(gas)) call AddDatum(header3, '[mmol+1s-1m-2]', separator)
                end if
            end do

            !> Advection fluxes
            header1 = header1(1:len_trim(header1)) // 'vertical_advection_fluxes'
            do gas = co2, gas4
                if (gas /= h2o) then
                    if(OutVarPresent(gas)) call AddDatum(header1, '', separator)
                    if(OutVarPresent(gas)) call AddDatum(header2, e2sg(gas)(1:len_trim(e2sg(gas))) // 'v-adv', separator)
                    if(OutVarPresent(gas)) call AddDatum(header3, '[' // char(181) // 'mol+1s-1m-2]', separator)
                else
                    if(OutVarPresent(gas)) call AddDatum(header1, '', separator)
                    if(OutVarPresent(gas)) call AddDatum(header2, e2sg(gas)(1:len_trim(e2sg(gas))) // 'v-adv', separator)
                    if(OutVarPresent(gas)) call AddDatum(header3, '[mmol+1s-1m-2]', separator)
                end if
            end do

            !> Average gas concentrations
            call AddDatum(header1,'gas_densities_concentrations_and_timelags', separator)
            do gas = co2, gas4
                if(OutVarPresent(gas)) call AddDatum(header1, ',,,,', separator)
                if(OutVarPresent(gas)) call AddDatum(header2, e2sg(gas)(1:len_trim(e2sg(gas))) // 'molar_density,' &
                    // e2sg(gas)(1:len_trim(e2sg(gas))) // 'mole_fraction,' &
                    // e2sg(gas)(1:len_trim(e2sg(gas))) // 'mixing_ratio,' &
                    // e2sg(gas)(1:len_trim(e2sg(gas))) // 'time_lag,' &
                    // e2sg(gas)(1:len_trim(e2sg(gas))) // 'def_timelag', separator)
                if (gas /= h2o) then
                    if(OutVarPresent(gas)) call AddDatum(header3, '[mmol+1m-3],[' // char(181) // &
                        'mol+1mol_a-1],[' // char(181) // 'mol+1mol_d-1],[s],[1=default]', separator)
                else
                    if(OutVarPresent(gas)) &
                        call AddDatum(header3, '[mmol+1m-3],[mmol+1mol_a-1],[mmol+1mol_d-1],[s],[1=default]', separator)
                end if
            end do
            !> In Header 1 there is one comma too much, take it away
            header1 = header1(1:len_trim(header1) - 1)

            !> Air properties, wind components and rotation angles
            call AddDatum(header1, 'air_properties,,,,,,,,,,,,,,unrotated_wind,,,rotated_wind&
                          &,,,,,,rotation_angles_for_tilt_correction,,', separator)
            call AddDatum(header2,'sonic_temperature,air_temperature,air_pressure,air_density,air_heat_capacity,air_molar_volume,&
                          &ET,water_vapor_density,e,es,specific_humidity,RH,VPD,Tdew&
                          &,u_unrot,v_unrot,w_unrot,u_rot,v_rot,w_rot,wind_speed,max_wind_speed,wind_dir,yaw,pitch,roll', separator)
            call AddDatum(header3,'[K],[K],[Pa],[kg+1m-3],[J+1kg-1K-1],[m+3mol-1],&
                          &[mm+1hour-1],[kg+1m-3],[Pa],[Pa],[kg+1kg-1],[%],[Pa],[K],&
                          &[m+1s-1],[m+1s-1],[m+1s-1],[m+1s-1],[m+1s-1],[m+1s-1],[m+1s-1],&
                          &[m+1s-1],[deg_from_north],[deg],[deg],[deg]', separator)

            !> Turbulence
            call AddDatum(header1, 'turbulence,,,,,', separator)
            call AddDatum(header2,'u*,TKE,L,(z-d)/L,bowen_ratio,T*', separator)
            call AddDatum(header3,'[m+1s-1],[m+2s-2],[m],[#],[#],[K]', separator)

            !> Footprint, if requested
            if (Meth%foot /= 'none') then
                call AddDatum(header1, 'footprint,,,,,,,', separator)
                call AddDatum(header2,'model,x_peak,x_offset,x_10%,x_30%,x_50%,x_70%,x_90%', separator)
                call AddDatum(header3,'[0=KJ/1=KM/2=HS],[m],[m],[m],[m],[m],[m],[m]', separator)
            end if

            !> uncorrected fluxes
            !> Tau and H
            call AddDatum(header1, 'uncorrected_fluxes,,,', separator)
            call AddDatum(header2,'un_Tau,Tau_scf,un_H,H_scf', separator)
            call AddDatum(header3,'[kg+1m-1s-2],[#],[W+1m-2],[#]', separator)
            !> LE
            if(OutVarPresent(h2o)) call AddDatum(header1, ',', separator)
            if(OutVarPresent(h2o)) call AddDatum(header2,'un_LE,LE_scf', separator)
            if(OutVarPresent(h2o)) call AddDatum(header3,'[W+1m-2],[#]', separator)
            !> Uncorrected gas fluxes (Level 0) and spectral correction factors
            do gas = co2, gas4
                if (gas /= h2o) then
                    if(OutVarPresent(gas)) call AddDatum(header1, ',', separator)
                    if(OutVarPresent(gas)) call AddDatum(header2, 'un_' // e2sg(gas)(1:len_trim(e2sg(gas))) &
                        // 'flux,' // e2sg(gas)(1:len_trim(e2sg(gas))) // 'scf', separator)
                    if(OutVarPresent(gas)) call AddDatum(header3, '[' // char(181) // 'mol+1s-1m-2],[#]', separator)
                else
                    if(OutVarPresent(gas)) call AddDatum(header1, ',', separator)
                    if(OutVarPresent(gas)) call AddDatum(header2, 'un_' // e2sg(gas)(1:len_trim(e2sg(gas))) &
                        // 'flux,' // e2sg(gas)(1:len_trim(e2sg(gas))) // 'scf', separator)
                    if(OutVarPresent(gas)) call AddDatum(header3, '[mmol+1s-1m-2],[#]', separator)
                end if
            end do

            !> Vickers and Mahrt 97 hard and soft flags
            call AddDatum(header1,'statistical_flags,,,,,,,,,,,', separator)
            call AddDatum(header2,'spikes_hf,amplitude_resolution_hf,drop_out_hf,absolute_limits_hf,&
                &skewness_kurtosis_hf,skewness_kurtosis_sf,discontinuities_hf,discontinuities_sf,timelag_hf,&
                &timelag_sf,attack_angle_hf,non_steady_wind_hf', separator)
            call AddDatum(header3,'8u/v/w/ts/co2/h2o/ch4/' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) &
                // ',8u/v/w/ts/co2/h2o/ch4/' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) &
                // ',8u/v/w/ts/co2/h2o/ch4/' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) &
                // ',8u/v/w/ts/co2/h2o/ch4/' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) &
                // ',8u/v/w/ts/co2/h2o/ch4/' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) &
                // ',8u/v/w/ts/co2/h2o/ch4/' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) &
                // ',8u/v/w/ts/co2/h2o/ch4/' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) &
                // ',8u/v/w/ts/co2/h2o/ch4/' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) &
                // ',8co2/h2o/ch4/' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) &
                // ',8co2/h2o/ch4/' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) &
                // ',8aa,8U', separator)

            !> Add spikes for EddyPro variables
            call AddDatum(header1,'spikes,,,', separator)
            call AddDatum(header2,'u_spikes,v_spikes,w_spikes,ts_spikes', separator)
            call AddDatum(header3,'[#],[#],[#],[#]', separator)
            do var = co2, gas4
                if(OutVarPresent(var)) then
                    call AddDatum(header1, '', separator)
                    call AddDatum(header2, e2sg(var)(1:len_trim(e2sg(var))) // 'spikes' , separator)
                    call AddDatum(header3, '[#]', separator)
                end if
            end do

            !> LI-COR's diagnostic flags
            if (Diag7200%present) then
                call AddDatum(header1,'diagnostic_flags_LI-7200,,,,,,,,', separator)
                call AddDatum(header2,'head_detect_LI-7200,t_out_LI-7200,t_in_LI-7200,aux_in_LI-7200,delta_p_LI-7200,&
                    &chopper_LI-7200,detector_LI-7200,pll_LI-7200,sync_LI-7200', separator)
                call AddDatum(header3,'[#_flagged_recs],[#_flagged_recs],[#_flagged_recs],[#_flagged_recs],[#_flagged_recs],&
                    &[#_flagged_recs],[#_flagged_recs],[#_flagged_recs],[#_flagged_recs]', separator)
            end if
            if (Diag7500%present) then
                call AddDatum(header1,'diagnostic_flags_LI-7500,,,', separator)
                call AddDatum(header2,'chopper_LI-7500,detector_LI-7500,pll_LI-7500,sync_LI-7500', separator)
                call AddDatum(header3,'[#_flagged_recs],[#_flagged_recs],[#_flagged_recs],[#_flagged_recs]', separator)
            end if
            if (Diag7700%present) then
                call AddDatum(header1,'diagnostic_flags_LI-7700,,,,,,,,,,,,,,,', separator)
                call AddDatum(header2,'not_ready_LI-7700,no_signal_LI-7700,re_unlocked_LI-7700,bad_temp_LI-7700,&
                    &laser_temp_unregulated_LI-7700,block_temp_unregulated_LI-7700,motor_spinning_LI-7700,&
                    &pump_on_LI-7700,top_heater_on_LI-7700,bottom_heater_on_LI-7700,calibrating_LI-7700,&
                    &motor_failure_LI-7700,bad_aux_tc1_LI-7700,bad_aux_tc2_LI-7700,&
                    &bad_aux_tc3_LI-7700,box_connected_LI-7700', separator)
                call AddDatum(header3,'[#_flagged_recs],[#_flagged_recs],[#_flagged_recs],[#_flagged_recs],[#_flagged_recs],&
                    &[#_flagged_recs],[#_flagged_recs],[#_flagged_recs],[#_flagged_recs],[#_flagged_recs],[#_flagged_recs],&
                    &[#_flagged_recs],[#_flagged_recs],[#_flagged_recs],[#_flagged_recs],[#_flagged_recs]', separator)
            end if

            !> AGCs and RSSIs for LI-7200 and LI-7500
            if (Diag7200%present) then
                if(CompareSwVer(E2Col(co2)%instr%sw_ver, SwVerFromString('6.0.0'))) then
                    call AddDatum(header1,'RSSI_LI-7200', separator)
                    call AddDatum(header2,'mean_value_RSSI_LI-7200', separator)
                    call AddDatum(header3,'[#]', separator)
                else
                    call AddDatum(header1,'AGC_LI-7200', separator)
                    call AddDatum(header2,'mean_value_AGC_LI-7200', separator)
                    call AddDatum(header3,'[#]', separator)
                end if
            end if
            if (Diag7500%present) then
                if(CompareSwVer(E2Col(co2)%instr%sw_ver, SwVerFromString('6.0.0'))) then
                    call AddDatum(header1,'RSSI_LI-7500', separator)
                    call AddDatum(header2,'mean_value_RSSI_LI-7500', separator)
                    call AddDatum(header3,'[#]', separator)
                else
                    call AddDatum(header1,'AGC_LI-7500', separator)
                    call AddDatum(header2,'mean_value_AGC_LI-7500', separator)
                    call AddDatum(header3,'[#]', separator)
                end if
            end if

            !> Variances
            call AddDatum(header1, 'variances,,,', separator)
            call AddDatum(header2, 'u_var,v_var,w_var,ts_var', separator)
            call AddDatum(header3, '[m+2s-2],[m+2s-2],[m+2s-2],[K+2]', separator)
            do gas = co2, gas4
                if(OutVarPresent(gas)) call AddDatum(header1, '', separator)
                if(OutVarPresent(gas)) call AddDatum(header2, e2sg(gas)(1:len_trim(e2sg(gas))) // 'var', separator)
                if(OutVarPresent(gas)) call AddDatum(header3, '--', separator)
            end do
            !> w/ts covariance
            call AddDatum(header1, 'covariances', separator)
            call AddDatum(header2,'w/ts_cov', separator)
            call AddDatum(header3,'[m+1K+1s-1]', separator)
            !> w-gases covariances
            do gas = co2, gas4
                if(OutVarPresent(gas)) call AddDatum(header1, '', separator)
                if(OutVarPresent(gas)) call AddDatum(header2, 'w/' // e2sg(gas)(1:len_trim(e2sg(gas))) // 'cov', separator)
                if(OutVarPresent(gas)) call AddDatum(header3, '--', separator)
            end do

            !> Mean values of user variables
            if (NumUserVar > 0) then
                call AddDatum(header1, 'custom_variables', separator)
                do var = 1, NumUserVar
                    call AddDatum(header1, '', separator)
                    call AddDatum(header2, usg(var)(1:len_trim(usg(var))) // 'mean', separator)
                    call AddDatum(header3,'--', separator)
                end do
            end if
            call latin1_to_utf8(header1, head1_utf8)
            call latin1_to_utf8(header2, head2_utf8)
            call latin1_to_utf8(header3, head3_utf8)

            !> Write on output file
            write(uflx, '(a)') head1_utf8(1:len_trim(head1_utf8) - 1)
            write(uflx, '(a)') head2_utf8(1:len_trim(head2_utf8) - 1)
            write(uflx, '(a)') head3_utf8(1:len_trim(head3_utf8) - 1)

        else
            header1 = 'file_info,,,,,,,corrected_fluxes_and_quality_flags,,,,,,,,,,,,,,,,,,,,,&
                &storage_fluxes,,,,,,vertical_advection_fluxes,,,,&
                &gas_densities_concentrations_and_timelags,,,,,,,,,,,,,,,,,,,,&
                &air_properties,,,,,,,,,,,,,,unrotated_wind,,,rotated_wind,,,,,,&
                &rotation_angles_for_tilt_correction,,,turbulence,,,,,,footprint,,,,,,,,&
                &uncorrected_fluxes_and_spectral_correction_factors_(scf),,,,,,,,,,,,,,&
                &statistical_flags,,,,,,,,,,,,spikes,,,,,,,,&
                &diagnostic_flags_LI-7200,,,,,,,,,&
                &diagnostic_flags_LI-7500,,,,diagnostic_flags_LI-7700,,,,,,,,,,,,,,,,'
                if(CompareSwVer(E2Col(co2)%instr%sw_ver, SwVerFromString('6.0.0'))) then
                    header1 = trim(header1) // 'RSSI_LI-7200,RSSI_LI-7500,variances,,,,,,,,covariances,,,,,'
                else
                    header1 = trim(header1) // 'AGC_LI-7200,AGC_LI-7500,variances,,,,,,,,covariances,,,,,'
                end if
            header2 = 'filename,date,time,DOY,daytime,file_records,used_records,Tau,qc_Tau,rand_err_Tau,&
                &H,qc_H,rand_err_H,LE,qc_LE,rand_err_LE,&
                &co2_flux,qc_co2_flux,rand_err_co2_flux,h2o_flux,qc_h2o_flux,rand_err_h2o_flux,ch4_flux,qc_ch4_flux,&
                &rand_err_ch4_flux,' &
                // e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'flux,qc_' &
                // e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'flux,rand_err_' &
                // e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'flux,H_strg,LE_strg,co2_strg,h2o_strg,ch4_strg,' &
                // e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'strg,co2_v-adv,h2o_v-adv,ch4_v-adv,' &
                // e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'v-adv,co2_molar_density,co2_mole_fraction,&
                &co2_mixing_ratio,co2_time_lag,co2_def_timelag,&
                &h2o_molar_density,h2o_mole_fraction,h2o_mixing_ratio,h2o_time_lag,h2o_def_timelag,&
                &ch4_molar_density,ch4_mole_fraction,ch4_mixing_ratio,ch4_time_lag,ch4_def_timelag,' &
                // e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'molar_density,' &
                // e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'mole_fraction,' &
                // e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'mixing_ratio,' &
                // e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'time_lag,' &
                // e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'def_timelag,&
                &sonic_temperature,air_temperature,air_pressure,air_density,air_heat_capacity,air_molar_volume,&
                &ET,water_vapor_density,e,es,specific_humidity,RH,VPD,Tdew&
                &,u_unrot,v_unrot,w_unrot,u_rot,v_rot,w_rot,wind_speed,max_wind_speed,wind_dir,yaw,pitch,roll,&
                &u*,TKE,L,(z-d)/L,bowen_ratio,T*,model,x_peak,x_offset,x_10%,x_30%,x_50%,x_70%,x_90%,&
                &un_Tau,Tau_scf,un_H,H_scf,un_LE,LE_scf,un_co2_flux,co2_scf,un_h2o_flux,h2o_scf,un_ch4_flux,ch4_scf,un_' &
                // e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'flux,un_' &
                // e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'scf,spikes_hf,amplitude_resolution_hf,drop_out_hf,absolute_limits_hf,&
                &skewness_kurtosis_hf,skewness_kurtosis_sf,discontinuities_hf,discontinuities_sf,timelag_hf,&
                &timelag_sf,attack_angle_hf,non_steady_wind_hf,u_spikes,v_spikes,&
                &w_spikes,ts_spikes,co2_spikes,h2o_spikes,ch4_spikes,' &
                 // e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'spikes,&
                &head_detect_LI-7200,t_out_LI-7200,t_in_LI-7200,aux_in_LI-7200,delta_p_LI-7200,&
                &chopper_LI-7200,detector_LI-7200,pll_LI-7200,sync_LI-7200,&
                &chopper_LI-7500,detector_LI-7500,pll_LI-7500,sync_LI-7500,&
                &not_ready_LI-7700,no_signal_LI-7700,re_unlocked_LI-7700,bad_temp_LI-7700,laser_temp_unregulated_LI-7700,&
                &block_temp_unregulated_LI-7700,motor_spinning_LI-7700,pump_on_LI-7700,top_heater_on_LI-7700,&
                &bottom_heater_on_LI-7700,calibrating_LI-7700,&
                &motor_failure_LI-7700,bad_aux_tc1_LI-7700,bad_aux_tc2_LI-7700,bad_aux_tc3_LI-7700,box_connected_LI-7700,&
                &mean_value_RSSI_LI-7200,mean_value_LI-7500,&
                &u_var,v_var,w_var,ts_var,co2_var,h2o_var,ch4_var,' // e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'var,&
                &w/ts_cov,w/co2_cov,w/h2o_cov,w/ch4_cov,w/' &
                // e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'cov,'
            header3 = ',[yyyy-mm-dd],[HH:MM],[ddd.ddd],[1=daytime],[#],[#],[kg+1m-1s-2],[#],[kg+1m-1s-2],&
                &[W+1m-2],[#],[W+1m-2],[W+1m-2],[#],[W+1m-2],&
                &[' // char(181) // 'mol+1s-1m-2],[#],[' // char(181) // 'mol+1s-1m-2],[mmol+1s-1m-2],[#],[mmol+1s-1m-2],&
                &[' // char(181) // 'mol+1s-1m-2],[#],[' // char(181) // 'mol+1s-1m-2],&
                &[' // char(181) // 'mol+1s-1m-2],[#],[' // char(181) // 'mol+1s-1m-2],&
                &[W+1m-2],[W+1m-2],[' // char(181) // 'mol+1s-1m-2],&
                &[mmol+1s-1m-2],[' // char(181) // 'mol+1s-1m-2],[' // char(181) // 'mol+1s-1m-2],&
                &[' // char(181) // 'mol+1s-1m-2],[mmol+1s-1m-2],[' // char(181) &
                // 'mol+1s-1m-2],[' // char(181) // 'mol+1s-1m-2],&
                &[mmol+1m-3],[' // char(181) // 'mol+1mol_a-1],[' // char(181) // 'mol+1mol_d-1],[s],[1=default],&
                &[mmol+1m-3],[mmol+1mol_a-1],[mmol+1mol_d-1],[s],[1=default],&
                &[mmol+1m-3],[' // char(181) // 'mol+1mol_a-1],[' // char(181) // 'mol+1mol_d-1],[s],[1=default],&
                &[mmol+1m-3],[' // char(181) // 'mol+1mol_a-1],[' // char(181) // 'mol+1mol_d-1],[s],[1=default],&
                &[K],[K],[Pa],[kg+1m-3],[J+1kg-1K-1],[m+3mol-1],[mm],[kg+1m-3],[Pa],[Pa],[kg+1kg-1],[%],[Pa],[K],&
                &[m+1s-1],[m+1s-1],[m+1s-1],[m+1s-1],[m+1s-1],[m+1s-1],[m+1s-1],[m+1s-1],[deg_from_north],[deg],[deg],[deg],&
                &[m+1s-1],[m+2s-2],[m],[#],[#],[K],[0=KJ/1=KM/2=HS],[m],[m],[m],[m],[m],[m],[m],&
                &[kg+1m-1s-2],[#],[W+1m-2],[#],[W+1m-2],[#],[' // char(181) // 'mol+1s-1m-2],[#],[mmol+1s-1m-2],[#],&
                &[' // char(181) // 'mol+1s-1m-2],[#],[' // char(181) // 'mol+1s-1m-2],[#],&
                &8u/v/w/ts/co2/h2o/ch4/' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) &
                // ',8u/v/w/ts/co2/h2o/ch4/' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) &
                // ',8u/v/w/ts/co2/h2o/ch4/' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) &
                // ',8u/v/w/ts/co2/h2o/ch4/' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) &
                // ',8u/v/w/ts/co2/h2o/ch4/' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) &
                // ',8u/v/w/ts/co2/h2o/ch4/' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) &
                // ',8u/v/w/ts/co2/h2o/ch4/' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) &
                // ',8u/v/w/ts/co2/h2o/ch4/' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) &
                // ',8co2/h2o/ch4/' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) &
                // ',8co2/h2o/ch4/' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) &
                // ',8aa,8U,[#],[#],[#],[#],[#],[#],[#],[#],&
                &[#_flagged_recs],[#_flagged_recs],[#_flagged_recs],[#_flagged_recs],[#_flagged_recs],&
                &[#_flagged_recs],[#_flagged_recs],[#_flagged_recs],[#_flagged_recs],&
                &[#_flagged_recs],[#_flagged_recs],[#_flagged_recs],[#_flagged_recs],&
                &[#_flagged_recs],[#_flagged_recs],[#_flagged_recs],[#_flagged_recs],[#_flagged_recs],&
                &[#_flagged_recs],[#_flagged_recs],[#_flagged_recs],[#_flagged_recs],[#_flagged_recs],&
                &[#_flagged_recs],[#_flagged_recs],[#_flagged_recs],[#_flagged_recs],[#_flagged_recs],[#_flagged_recs],&
                &[#],[#],[m+2s-2],[m+2s-2],[m+2s-2],[K+2],--,--,--,--,[m+1s-1K+1],--,--,--,--,'
            !> Mean values of user variables
            if (NumUserVar > 0) then
                call AddDatum(header1, 'custom_variables', separator)
                do var = 1, NumUserVar
                    call AddDatum(header1, '', separator)
                    call AddDatum(header2, usg(var)(1:len_trim(usg(var))) // 'mean', separator)
                    call AddDatum(header3,'--', separator)
                end do
            end if
            call latin1_to_utf8(header1, head1_utf8)
            call latin1_to_utf8(header2, head2_utf8)
            call latin1_to_utf8(header3, head3_utf8)

            !> Write on output file
            write(uflx, '(a)') head1_utf8(1:len_trim(head1_utf8) - 1)
            write(uflx, '(a)') head2_utf8(1:len_trim(head2_utf8) - 1)
            write(uflx, '(a)') head3_utf8(1:len_trim(head3_utf8) - 1)
        end if
    end if

    !>==========================================================================
    !>==========================================================================
    !> Metadata output
    if (EddyProProj%out_md) then
        Test_Path = Dir%main_out(1:len_trim(Dir%main_out)) &
                  // EddyProProj%id(1:len_trim(EddyProProj%id)) &
                  // MetaData_FilePadding // Timestamp_FilePadding // CsvExt
        dot = index(Test_Path, CsvExt, .true.) - 1
        Metadata_Path = Test_Path(1:dot) // CsvTmpExt
        open(umd, file = Metadata_Path, iostat = open_status, encoding = 'utf-8')

        call Clearstr(dataline)
        call AddDatum(dataline,'filename,date,time,DOY,latitude,longitude,altitude,canopy_height,displacement_height,&
            &roughness_length,file_length,acquisition_frequency,&
            &master_sonic_manufacturer,master_sonic_model,master_sonic_height,&
            &master_sonic_wformat,master_sonic_wref,master_sonic_north_offset,&
            &master_sonic_hpath_length,master_sonic_vpath_length,master_sonic_tau', separator)
        if (OutVarPresent(co2)) &
            call AddDatum(dataline,'co2_irga_manufacturer,co2_irga_model,co2_measure_type,co2_irga_northward_separation,&
                &co2_irga_eastward_separation,co2_irga_vertical_separation,&
                &co2_irga_tube_length,co2_irga_tube_diameter,co2_irga_tube_flowrate,&
                &co2_irga_kw,co2_irga_ko,co2_irga_hpath_length,co2_irga_vpath_length,co2_irga_tau', separator)
        if (OutVarPresent(h2o)) &
            call AddDatum(dataline,'h2o_irga_manufacturer,h2o_irga_model,h2o_measure_type,h2o_irga_northward_separation,&
                &h2o_irga_eastward_separation,h2o_irga_vertical_separation,&
                &h2o_irga_tube_length,h2o_irga_tube_diameter,h2o_irga_tube_flowrate,&
                &h2o_irga_kw,h2o_irga_ko,h2o_irga_hpath_length,h2o_irga_vpath_length,h2o_irga_tau', separator)
        if (OutVarPresent(ch4)) &
            call AddDatum(dataline,'ch4_irga_manufacturer,ch4_irga_model,ch4_measure_type,ch4_irga_northward_separation,&
                &ch4_irga_eastward_separation,ch4_irga_vertical_separation,&
                &ch4_irga_tube_length,ch4_irga_tube_diameter,ch4_irga_tube_flowrate,&
                &ch4_irga_kw,ch4_irga_ko,ch4_irga_hpath_length,ch4_irga_vpath_length,ch4_irga_tau', separator)
        if (OutVarPresent(gas4)) &
            call AddDatum(dataline, e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'irga_manufacturer,' &
                // e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'irga_model,' &
                // e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'measure_type,' &
                // e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'irga_northward_separation,' &
                // e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'irga_eastward_separation,' &
                // e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'irga_vertical_separation,' &
                // e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'irga_tube_length,' &
                // e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'irga_tube_diameter,' &
                // e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'irga_tube_flowrate' &
                // e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'irga_kw' &
                // e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'irga_ko' &
                // e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'irga_hpath_length' &
                // e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'irga_vpath_length' &
                // e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'irga_tau', separator)
        write(umd, '(a)') dataline(1:len_trim(dataline) - 1)
    end if

    !>*********************************************************************************************
    !>*********************************************************************************************

    !> Details of stationarity and integral turbulence tests
    if(RPsetup%out_qc_details .and. Meth%qcflag /= 'none') then
        !> Open file
        Test_Path = Dir%main_out(1:len_trim(Dir%main_out)) &
                  // EddyProProj%id(1:len_trim(EddyProProj%id)) &
                  // QCdetails_FilePadding // Timestamp_FilePadding // CsvExt
        dot = index(Test_Path, CsvExt, .true.) - 1
        QCdetails_Path = Test_Path(1:dot) // CsvTmpExt
        open(uqc, file = QCdetails_Path, iostat = open_status, encoding = 'utf-8')

        call Clearstr(header1)
        call Clearstr(header2)
        call Clearstr(header3)
        call Clearstr(head1_utf8)
        call Clearstr(head2_utf8)
        call Clearstr(head3_utf8)
        call AddDatum(header1,'file_info,,,,stationarity test,,', separator)
        call AddDatum(header2,'filename,date,time,DOY,dev(u),dev(w),dev(ts)', separator)
        call AddDatum(header3,',[yyyy-mm-dd],[HH:MM],[ddd.ddd],[%],[%],[%]', separator)
        if (OutVarPresent(co2)) then
            call AddDatum(header1,'', separator)
            call AddDatum(header2,'dev(co2)', separator)
            call AddDatum(header3,'[%]', separator)
        end if
        if (OutVarPresent(h2o)) then
            call AddDatum(header1,'', separator)
            call AddDatum(header2,'dev(h2o)', separator)
            call AddDatum(header3,'[%]', separator)
        end if
        if (OutVarPresent(ch4)) then
            call AddDatum(header1,'', separator)
            call AddDatum(header2,'dev(ch4)', separator)
            call AddDatum(header3,'[%]', separator)
        end if
        if (OutVarPresent(gas4)) then
            call AddDatum(header1,'', separator)
            call AddDatum(header2,'dev(' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1)  // ')', separator)
            call AddDatum(header3,'[%]', separator)
        end if

        call AddDatum(header1,',', separator)
        call AddDatum(header2,'dev(w/u),dev(w/ts)', separator)
        call AddDatum(header3,'[%],[%]', separator)

        if (OutVarPresent(co2)) then
            call AddDatum(header1,'', separator)
            call AddDatum(header2,'dev(w/co2)', separator)
            call AddDatum(header3,'[%]', separator)
        end if
        if (OutVarPresent(h2o)) then
            call AddDatum(header1,'', separator)
            call AddDatum(header2,'dev(w/h2o)', separator)
            call AddDatum(header3,'[%]', separator)
        end if
        if (OutVarPresent(ch4)) then
            call AddDatum(header1,'', separator)
            call AddDatum(header2,'dev(w/ch4)', separator)
            call AddDatum(header3,'[%]', separator)
        end if
        if (OutVarPresent(gas4)) then
            call AddDatum(header1,'', separator)
            call AddDatum(header2,'dev(w/' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1)  // ')', separator)
            call AddDatum(header3,'[%]', separator)
        end if

        call AddDatum(header1,',', separator)
        call AddDatum(header2,'flag(w/u),flag(w/ts)', separator)
        call AddDatum(header3,'[#],[#]', separator)

        if (OutVarPresent(co2)) then
            call AddDatum(header1,'', separator)
            call AddDatum(header2,'flag(w/co2)', separator)
            call AddDatum(header3,'[#]', separator)
        end if
        if (OutVarPresent(h2o)) then
            call AddDatum(header1,'', separator)
            call AddDatum(header2,'flag(w/h2o)', separator)
            call AddDatum(header3,'[#]', separator)
        end if
        if (OutVarPresent(ch4)) then
            call AddDatum(header1,'', separator)
            call AddDatum(header2,'flag(w/ch4)', separator)
            call AddDatum(header3,'[#]', separator)
        end if
        if (OutVarPresent(gas4)) then
            call AddDatum(header1,'', separator)
            call AddDatum(header2,'flag(w/' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1)  // ')', separator)
            call AddDatum(header3,'[#]', separator)
        end if

        call AddDatum(header1,'well-developed_turbulence_test,,,,,', separator)
        call AddDatum(header2,'dev(u),dev(w),dev(ts),flag(u),flag(w),flag(ts)', separator)
        call AddDatum(header3,'[%],[%],[%],[#],[#],[#]', separator)

        call latin1_to_utf8(header1, head1_utf8)
        call latin1_to_utf8(header2, head2_utf8)
        call latin1_to_utf8(header3, head3_utf8)

        !> Write on output file
        write(uqc, '(a)') head1_utf8(1:len_trim(head1_utf8) - 1)
        write(uqc, '(a)') head2_utf8(1:len_trim(head2_utf8) - 1)
        write(uqc, '(a)') head3_utf8(1:len_trim(head3_utf8) - 1)
    end if

    !>*********************************************************************************************
    !>*********************************************************************************************

    !> Statistics files Level 1
    if (RPsetup%out_st(1)) then
        Test_Path = StatsDir(1:len_trim(StatsDir)) &
                  // EddyProProj%id(1:len_trim(EddyProProj%id)) &
                  // Stats1_FilePadding // Timestamp_FilePadding // CsvExt
        dot = index(Test_Path, CsvExt, .true.) - 1
        St1_Path = Test_Path(1:dot) // CsvTmpExt
        open(ust1, file = St1_Path, iostat = open_status, encoding = 'utf-8')

        write(ust1, '(a)') 'first_statistics:_on_raw_data'
        write(ust1, '(a)') 'filename,date,time,DOY,used_records,&
                           &mean(u),mean(v),mean(w),mean(ts),mean(co2),mean(h2o),&
                           &mean(ch4),mean(' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),mean(tc),mean(pc),mean(te),&
                           &mean(pe),WindDirection,&
                           &var(u),var(v),var(w),var(ts),var(co2),var(h2o),&
                           &var(ch4),var(' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),var(tc),var(pc),var(te),var(pe),&
                           &cov(u/v),cov(u/w),cov(u/ts),cov(u/co2),cov(u/h2o),&
                           &cov(u/ch4),cov(u/' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),&
                           &cov(v/w),cov(v/ts),cov(v/co2),cov(v/h2o),cov(v/ch4),&
                           &cov(v/' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),&
                           &cov(w/ts),cov(w/co2),cov(w/h2o),cov(w/ch4),cov(w/' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),&
                           &cov(w/tc),cov(w/pc),cov(w/te),cov(w/pe),&
                           &st_dev(u),st_dev(v),st_dev(w),st_dev(ts),st_dev(co2),st_dev(h2o),&
                           &st_dev(ch4),st_dev(' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),&
                           &st_dev(tc),st_dev(pc),st_dev(te),st_dev(pe),&
                           &skw(u),skw(v),skw(w),skw(ts),skw(co2),skw(h2o),&
                           &skw(ch4),skw(' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),skw(tc),skw(pc),skw(te),skw(pe),&
                           &kur(u),kur(v),kur(w),kur(ts),kur(co2),kur(h2o),&
                           &kur(ch4),kur(' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),kur(tc),kur(pc),kur(te),kur(pe)'
    end if

    !> Statistics files Level 2
    if (RPsetup%out_st(2)) then
        Test_Path = StatsDir(1:len_trim(StatsDir)) &
                  // EddyProProj%id(1:len_trim(EddyProProj%id)) &
                  // Stats2_FilePadding // Timestamp_FilePadding // CsvExt
        dot = index(Test_Path, CsvExt, .true.) - 1
        St2_Path = Test_Path(1:dot) // CsvTmpExt
        open(ust2, file = St2_Path, iostat = open_status, encoding = 'utf-8')

        write(ust2, '(a)') 'second_statistics:_on_raw_data_after_after_despiking'
        write(ust2, '(a)') 'filename,date,time,DOY,used_records,&
                           &mean(u),mean(v),mean(w),mean(ts),mean(co2),mean(h2o),&
                           &mean(ch4),mean(' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),mean(tc),mean(pc),mean(te),&
                           &mean(pe),WindDirection,&
                           &var(u),var(v),var(w),var(ts),var(co2),var(h2o),&
                           &var(ch4),var(' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),var(tc),var(pc),var(te),var(pe),&
                           &cov(u/v),cov(u/w),cov(u/ts),cov(u/co2),cov(u/h2o),&
                           &cov(u/ch4),cov(u/' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),&
                           &cov(v/w),cov(v/ts),cov(v/co2),cov(v/h2o),cov(v/ch4),&
                           &cov(v/' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),&
                           &cov(w/ts),cov(w/co2),cov(w/h2o),cov(w/ch4),cov(w/' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),&
                           &cov(w/tc),cov(w/pc),cov(w/te),cov(w/pe),&
                           &st_dev(u),st_dev(v),st_dev(w),st_dev(ts),st_dev(co2),st_dev(h2o),&
                           &st_dev(ch4),st_dev(' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),&
                           &st_dev(tc),st_dev(pc),st_dev(te),st_dev(pe),&
                           &skw(u),skw(v),skw(w),skw(ts),skw(co2),skw(h2o),&
                           &skw(ch4),skw(' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),skw(tc),skw(pc),skw(te),skw(pe),&
                           &kur(u),kur(v),kur(w),kur(ts),kur(co2),kur(h2o),&
                           &kur(ch4),kur(' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),kur(tc),kur(pc),kur(te),kur(pe)'
    end if

    !> Statistics files Level 3
    if (RPsetup%out_st(3)) then
        Test_Path = StatsDir(1:len_trim(StatsDir)) &
                  // EddyProProj%id(1:len_trim(EddyProProj%id)) &
                  // Stats3_FilePadding // Timestamp_FilePadding // CsvExt
        dot = index(Test_Path, CsvExt, .true.) - 1
        St3_Path = Test_Path(1:dot) // CsvTmpExt
        open(ust3, file = St3_Path, iostat = open_status, encoding = 'utf-8')

        write(ust3, '(a)') 'third_statistics:_on_raw_data_after_after_despiking_and_cross-wind_correction'
        write(ust3, '(a)') 'filename,date,time,DOY,used_records,&
                           &mean(u),mean(v),mean(w),mean(ts),mean(co2),mean(h2o),&
                           &mean(ch4),mean(' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),mean(tc),mean(pc),mean(te),&
                           &mean(pe),WindDirection,&
                           &var(u),var(v),var(w),var(ts),var(co2),var(h2o),&
                           &var(ch4),var(' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),var(tc),var(pc),var(te),var(pe),&
                           &cov(u/v),cov(u/w),cov(u/ts),cov(u/co2),cov(u/h2o),&
                           &cov(u/ch4),cov(u/' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),&
                           &cov(v/w),cov(v/ts),cov(v/co2),cov(v/h2o),cov(v/ch4),&
                           &cov(v/' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),&
                           &cov(w/ts),cov(w/co2),cov(w/h2o),cov(w/ch4),cov(w/' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),&
                           &cov(w/tc),cov(w/pc),cov(w/te),cov(w/pe),&
                           &st_dev(u),st_dev(v),st_dev(w),st_dev(ts),st_dev(co2),st_dev(h2o),&
                           &st_dev(ch4),st_dev(' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),&
                           &st_dev(tc),st_dev(pc),st_dev(te),st_dev(pe),&
                           &skw(u),skw(v),skw(w),skw(ts),skw(co2),skw(h2o),&
                           &skw(ch4),skw(' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),skw(tc),skw(pc),skw(te),skw(pe),&
                           &kur(u),kur(v),kur(w),kur(ts),kur(co2),kur(h2o),&
                           &kur(ch4),kur(' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),kur(tc),kur(pc),kur(te),kur(pe)'
    end if

    !> Statistics files Level 4
    if (RPsetup%out_st(4)) then
        Test_Path = StatsDir(1:len_trim(StatsDir)) &
                  // EddyProProj%id(1:len_trim(EddyProProj%id)) &
                  // Stats4_FilePadding // Timestamp_FilePadding // CsvExt
        dot = index(Test_Path, CsvExt, .true.) - 1
        St4_Path = Test_Path(1:dot) // CsvTmpExt
        open(ust4, file = St4_Path, iostat = open_status, encoding = 'utf-8')

        write(ust4, '(a)') 'forth statistics:_on_raw_data_after_despiking_cross_wind_correction&
                            &_and_angle-of-attack_correction'
        write(ust4, '(a)') 'filename,date,time,DOY,used_records,&
                           &mean(u),mean(v),mean(w),mean(ts),mean(co2),mean(h2o),&
                           &mean(ch4),mean(' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),mean(tc),mean(pc),mean(te),&
                           &mean(pe),WindDirection,&
                           &var(u),var(v),var(w),var(ts),var(co2),var(h2o),&
                           &var(ch4),var(' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),var(tc),var(pc),var(te),var(pe),&
                           &cov(u/v),cov(u/w),cov(u/ts),cov(u/co2),cov(u/h2o),&
                           &cov(u/ch4),cov(u/' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),&
                           &cov(v/w),cov(v/ts),cov(v/co2),cov(v/h2o),cov(v/ch4),&
                           &cov(v/' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),&
                           &cov(w/ts),cov(w/co2),cov(w/h2o),cov(w/ch4),cov(w/' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),&
                           &cov(w/tc),cov(w/pc),cov(w/te),cov(w/pe),&
                           &st_dev(u),st_dev(v),st_dev(w),st_dev(ts),st_dev(co2),st_dev(h2o),&
                           &st_dev(ch4),st_dev(' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),&
                           &st_dev(tc),st_dev(pc),st_dev(te),st_dev(pe),&
                           &skw(u),skw(v),skw(w),skw(ts),skw(co2),skw(h2o),&
                           &skw(ch4),skw(' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),skw(tc),skw(pc),skw(te),skw(pe),&
                           &kur(u),kur(v),kur(w),kur(ts),kur(co2),kur(h2o),&
                           &kur(ch4),kur(' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),kur(tc),kur(pc),kur(te),kur(pe)'
    end if

    !> Statistics files Level 5
    if (RPsetup%out_st(5)) then
        Test_Path = StatsDir(1:len_trim(StatsDir)) &
                  // EddyProProj%id(1:len_trim(EddyProProj%id)) &
                  // Stats5_FilePadding // Timestamp_FilePadding // CsvExt
        dot = index(Test_Path, CsvExt, .true.) - 1
        St5_Path = Test_Path(1:dot) // CsvTmpExt
        open(ust5, file = St5_Path, iostat = open_status, encoding = 'utf-8')

        write(ust5, '(a)') 'fifth_statistics:_on_raw_data_after_despiking_cross_wind_correction&
                            &_angle-of-attack_correction_and_tilt_correction'
        write(ust5, '(a)') 'filename,date,time,DOY,used_records,&
                           &mean(u),mean(v),mean(w),mean(ts),mean(co2),mean(h2o),&
                           &mean(ch4),mean(' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),mean(tc),mean(pc),mean(te),&
                           &mean(pe),WindDirection,&
                           &var(u),var(v),var(w),var(ts),var(co2),var(h2o),&
                           &var(ch4),var(' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),var(tc),var(pc),var(te),var(pe),&
                           &cov(u/v),cov(u/w),cov(u/ts),cov(u/co2),cov(u/h2o),&
                           &cov(u/ch4),cov(u/' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),&
                           &cov(v/w),cov(v/ts),cov(v/co2),cov(v/h2o),cov(v/ch4),&
                           &cov(v/' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),&
                           &cov(w/ts),cov(w/co2),cov(w/h2o),cov(w/ch4),cov(w/' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),&
                           &cov(w/tc),cov(w/pc),cov(w/te),cov(w/pe),&
                           &st_dev(u),st_dev(v),st_dev(w),st_dev(ts),st_dev(co2),st_dev(h2o),&
                           &st_dev(ch4),st_dev(' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),&
                           &st_dev(tc),st_dev(pc),st_dev(te),st_dev(pe),&
                           &skw(u),skw(v),skw(w),skw(ts),skw(co2),skw(h2o),&
                           &skw(ch4),skw(' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),skw(tc),skw(pc),skw(te),skw(pe),&
                           &kur(u),kur(v),kur(w),kur(ts),kur(co2),kur(h2o),&
                           &kur(ch4),kur(' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),kur(tc),kur(pc),kur(te),kur(pe)'
    end if

    !> Statistics files Level 6
    if (RPsetup%out_st(6)) then
        Test_Path = StatsDir(1:len_trim(StatsDir)) &
                  // EddyProProj%id(1:len_trim(EddyProProj%id)) &
                  // Stats6_FilePadding // Timestamp_FilePadding // CsvExt
        dot = index(Test_Path, CsvExt, .true.) - 1
        St6_Path = Test_Path(1:dot) // CsvTmpExt
        open(ust6, file = St6_Path, iostat = open_status, encoding = 'utf-8')

        write(ust6, '(a)') 'sixth statistics:_on_raw_data_after_despiking_cross_wind_correction&
            &_angle-of-attack_correction_tilt_correction_and_time-lag_compensation'
        write(ust6, '(a)') 'filename,date,time,DOY,used_records,&
                           &mean(u),mean(v),mean(w),mean(ts),mean(co2),mean(h2o),&
                           &mean(ch4),mean(' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),mean(tc),mean(pc),mean(te),&
                           &mean(pe),WindDirection,&
                           &var(u),var(v),var(w),var(ts),var(co2),var(h2o),&
                           &var(ch4),var(' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),var(tc),var(pc),var(te),var(pe),&
                           &cov(u/v),cov(u/w),cov(u/ts),cov(u/co2),cov(u/h2o),&
                           &cov(u/ch4),cov(u/' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),&
                           &cov(v/w),cov(v/ts),cov(v/co2),cov(v/h2o),cov(v/ch4),&
                           &cov(v/' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),&
                           &cov(w/ts),cov(w/co2),cov(w/h2o),cov(w/ch4),cov(w/' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),&
                           &cov(w/tc),cov(w/pc),cov(w/te),cov(w/pe),&
                           &st_dev(u),st_dev(v),st_dev(w),st_dev(ts),st_dev(co2),st_dev(h2o),&
                           &st_dev(ch4),st_dev(' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),&
                           &st_dev(tc),st_dev(pc),st_dev(te),st_dev(pe),&
                           &skw(u),skw(v),skw(w),skw(ts),skw(co2),skw(h2o),&
                           &skw(ch4),skw(' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),skw(tc),skw(pc),skw(te),skw(pe),&
                           &kur(u),kur(v),kur(w),kur(ts),kur(co2),kur(h2o),&
                           &kur(ch4),kur(' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),kur(tc),kur(pc),kur(te),kur(pe)'
    end if

    !> Statistics files Level 7
    if (RPsetup%out_st(7)) then
        Test_Path = StatsDir(1:len_trim(StatsDir)) &
                  // EddyProProj%id(1:len_trim(EddyProProj%id)) &
                  // Stats7_FilePadding // Timestamp_FilePadding // CsvExt
        dot = index(Test_Path, CsvExt, .true.) - 1
        St7_Path = Test_Path(1:dot) // CsvTmpExt
        open(ust7, file = St7_Path, iostat = open_status, encoding = 'utf-8')

        write(ust7, '(a)') 'seventh_statistics:seventh_statistics:_on_raw_data_after_despiking_cross_wind_correction&
            &_angle-of-attack_correction_tilt_correction_time-lag_compensation_and_detrending'
        write(ust7, '(a)') 'filename,date,time,DOY,used_records,&
                           &mean(u),mean(v),mean(w),mean(ts),mean(co2),mean(h2o),&
                           &mean(ch4),mean(' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),mean(tc),mean(pc),mean(te),&
                           &mean(pe),WindDirection,&
                           &var(u),var(v),var(w),var(ts),var(co2),var(h2o),&
                           &var(ch4),var(' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),var(tc),var(pc),var(te),var(pe),&
                           &cov(u/v),cov(u/w),cov(u/ts),cov(u/co2),cov(u/h2o),&
                           &cov(u/ch4),cov(u/' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),&
                           &cov(v/w),cov(v/ts),cov(v/co2),cov(v/h2o),cov(v/ch4),&
                           &cov(v/' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),&
                           &cov(w/ts),cov(w/co2),cov(w/h2o),cov(w/ch4),cov(w/' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),&
                           &cov(w/tc),cov(w/pc),cov(w/te),cov(w/pe),&
                           &st_dev(u),st_dev(v),st_dev(w),st_dev(ts),st_dev(co2),st_dev(h2o),&
                           &st_dev(ch4),st_dev(' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),&
                           &st_dev(tc),st_dev(pc),st_dev(te),st_dev(pe),&
                           &skw(u),skw(v),skw(w),skw(ts),skw(co2),skw(h2o),&
                           &skw(ch4),skw(' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),skw(tc),skw(pc),skw(te),skw(pe),&
                           &kur(u),kur(v),kur(w),kur(ts),kur(co2),kur(h2o),&
                           &kur(ch4),kur(' // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // '),kur(tc),kur(pc),kur(te),kur(pe)'
    end if
end subroutine InitOutFiles_rp
