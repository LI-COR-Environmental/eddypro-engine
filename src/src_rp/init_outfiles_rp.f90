!***************************************************************************
! init_outfiles_rp.f90
! --------------------
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
    integer :: today(3), now(3)
    character(8) :: dum_string
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

    do j = 1, NumUserVar
        usg(j)  = UserCol(j)%label(1:len_trim(UserCol(j)%label)) // '_'
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
                    call AddDatum(header2, 'rand_err' // e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'flux', separator)
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
    !> Essentials output, for use in Fluxes
    if (EddyProProj%out_essentials) then
        Test_Path = Dir%main_out(1:len_trim(Dir%main_out)) &
                  // EddyProProj%id(1:len_trim(EddyProProj%id)) &
                  // Essentials_FilePadding // Timestamp_FilePadding // CsvExt
        dot = index(Test_Path, CsvExt, .true.) - 1
        Essentials_Path = Test_Path(1:dot) // CsvTmpExt
        open(uex, file = Essentials_Path, iostat = open_status, encoding = 'utf-8')

        call clearstr(dataline)
        dataline = 'filename,date,time,daytime,file_records,used_records,&
            &Tau,ru_Tau,H,ru_H,LE,ru_LE,co2_flux,ru_co2,h2o_flux,ru_h2o,ch4_flux,ru_ch4,' &
            // e2sg(gas4)(1:len_trim(e2sg(gas4))-1) // '_flux,ru_' // e2sg(gas4)(1:len_trim(e2sg(gas4))-1) // &
            ',H_stor,LE_stor,co2_stor,h2o_stor,ch4_stor,' // e2sg(gas4)(1:len_trim(e2sg(gas4))-1) // '_stor,&
            &E_co2,E_ch4,E_' // e2sg(gas4)(1:len_trim(e2sg(gas4))-1) // ',&
            &H_co2,H_h2o,H_ch4,H_' // e2sg(gas4)(1:len_trim(e2sg(gas4))-1) // ',&
            &unrot_u,unrot_v,unrot_w,rot_u,rot_v,rot_w,wind_speed,max_wind_speed,wind_dir,u*,TKE,L,(z-d)/L,bowen_ratio,T*,&
            &measure_type_co2,d_co2,r_co2,chi_co2,measure_type_h2o,d_h2o,r_h2o,chi_h2o,&
            &measure_type_ch4,d_ch4,r_ch4,chi_ch4,measure_type_' // e2sg(gas4)(1:len_trim(e2sg(gas4))-1) // ',&
            &d_' // e2sg(gas4)(1:len_trim(e2sg(gas4))-1) // ',&
            &r_' // e2sg(gas4)(1:len_trim(e2sg(gas4))-1) // ',&
            &chi_' // e2sg(gas4)(1:len_trim(e2sg(gas4))-1) // ',&
            &t_sonic,t_air,p_air,rh_air,v_air,rho_air,rhocp_air,&
            &rho_w,e,es,q,vpd,td,p_dry,rho_dry,v_dry,lambda,sigma,&
            &t_cell,p_cell,v_cell_co2,v_cell_h2o,v_cell_ch4,v_cell_' // e2sg(gas4)(1:len_trim(e2sg(gas4))-1) // ',&
            &li7700_a,li7700_b,li7700_c,bu_h_bot,bu_h_top,bu_h_spar,&
            &cov_wT,cov_wT_1.626,cov_wT_0.614,cov_wT_0.277,cov_wT_0.133,cov_wT_0.065,cov_wT_0.032&
            &,cov_wT_0.016,cov_wT_0.008,cov_wT_0.004,&
            &var(u),var(v),var(w),var(ts),var(co2),var(h2o),&
            &var(ch4),var(' // e2sg(gas4)(1:len_trim(e2sg(gas4))-1) // '),var(tc),var(pc),var(te),var(pe),&
            &cov(w/u),cov(w/v),cov(w/ts),cov(w/co2),cov(w/h2o),cov(w/ch4),&
            &cov(w/' // e2sg(gas4)(1:len_trim(e2sg(gas4))-1) // '),cov(w/tc),cov(w/pc),cov(w/te),cov(w/pe),&
            &tlag_co2,def_tlag_co2,tlag_h2o,def_tlag_h2o,tlag_ch4,def_tlag_ch4,&
            &tlag_' // e2sg(gas4)(1:len_trim(e2sg(gas4))-1) // ',def_tlag_' // e2sg(gas4)(1:len_trim(e2sg(gas4))-1) // ',&
            &yaw,pitch,roll,&
            &st_w_u,st_w_ts,st_w_co2,st_w_h2o,st_w_ch4,st_w_' // e2sg(gas4)(1:len_trim(e2sg(gas4))-1) // ',&
            &dt_u,dt_w,dt_ts,&
            &detrending_method,dentrending_time_constant,&
            &logger_swver_major,logger_swver_minor,logger_swver_revision,&
            &latitude,longitude,altitude,file_length,&
            &averaging_interval,acquisition_frequency,&
            &canopy_height,displacement_height,roughness_length,&
            &master_sonic_manufacturer,master_sonic_model,master_sonic_height,&
            &master_sonic_wformat,master_sonic_wref,master_sonic_north_offset,&
            &master_sonic_hpath,master_sonic_vpath,master_sonic_tau,&
            &co2_irga_manufacturer,co2_irga_model,&
            &co2_irga_northward_separation,co2_irga_eastward_separation,co2_irga_vertical_separation,&
            &co2_irga_tube_length,co2_irga_tube_diameter,co2_irga_tube_flowrate, &
            &co2_irga_kw,co2_irga_ko,co2_irga_hpath_length,co2_irga_vpath_length,co2_irga_tau, &
            &h2o_irga_manufacturer,h2o_irga_model,&
            &h2o_irga_northward_separation,h2o_irga_eastward_separation,h2o_irga_vertical_separation,&
            &h2o_irga_tube_length,h2o_irga_tube_diameter,h2o_irga_tube_flowrate, &
            &h2o_irga_kw,h2o_irga_ko,h2o_irga_hpath_length,h2o_irga_vpath_length,h2o_irga_tau, &
            &ch4_irga_manufacturer,ch4_irga_model,&
            &ch4_irga_northward_separation,ch4_irga_eastward_separation,ch4_irga_vertical_separation,&
            &ch4_irga_tube_length,ch4_irga_tube_diameter,ch4_irga_tube_flowrate, &
            &ch4_irga_kw,ch4_irga_ko,ch4_irga_hpath_length,ch4_irga_vpath_length,ch4_irga_tau,' &
            // e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'irga_manufacturer,' &
            // e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'irga_model,' &
            // e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'irga_northward_separation,' &
            // e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'irga_eastward_separation,' &
            // e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'irga_vertical_separation,' &
            // e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'irga_tube_length,' &
            // e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'irga_tube_diameter,' &
            // e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'irga_tube_flowrate,' &
            // e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'irga_kw,' &
            // e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'irga_ko,' &
            // e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'irga_hpath_length,' &
            // e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'irga_vpath_length,' &
            // e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'irga_tau,&
            &hf_sr,hf_ar,hf_do,hf_al,hf_sk,sf_sk,hf_ds,sf_ds,hf_tl,sf_tl,hf_aa,sf_ns,&
            &spikes_u,spikes_v,spikes_w,spikes_ts,spikes_co2,spikes_h2o,spikes_ch4,spikes_' &
            // e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1) // ',&
            &li72_head_detect,li72_t_out,li72_t_in,li72_aux_in,li72_delta_p,li72_chopper,li72_detector,li72_pll,li72_sync,&
            &li75_chopper,li75_detector,li75_pll,li75_sync,&
            &li77_not_ready,li77_no_signal,li77_re_unlocked,li77_bad_temp,li77_laser_t_unreg,li77_clock_t_unreg,&
            &li77_motor_spinning,li77_pump_on,li77_top_heater_on,li77_bottom_heater_on,li77_calibrating,li77_motor_failure,&
            &li77_bad_aux_tc1,li77_bad_aux_tc2,li77_bad_aux_tc3,li77_box_connected,&
            &li72_AGC,li75_AGC,li77_RSSI,num_user_var,'
        if (NumUserVar > 0) then
            do i = 1, NumUserVar
                dataline = dataline(1:len_trim(dataline)) &
                    // usg(i)(1:len_trim(usg(i))) // 'mean' // ','
            end do
        end if
        write(uex, '(a)') dataline(1:len_trim(dataline) - 1)
    end if

    !>==========================================================================
    !>==========================================================================
    !> ICOS output
    if (EddyProProj%out_icos) then
        Test_Path = Dir%main_out(1:len_trim(Dir%main_out)) &
                  // EddyProProj%id(1:len_trim(EddyProProj%id)) &
                  // ICOS_FilePadding // Timestamp_FilePadding // CsvExt
        dot = index(Test_Path, CsvExt, .true.) - 1
        ICOS_Path = Test_Path(1:dot) // CsvTmpExt
        open(uicos, file = ICOS_Path, iostat = open_status, encoding = 'utf-8')

        call clearstr(dataline)
        dataline = 'TIMESTAMP,POTENTIAL_RADIATION,DAYTIME,NR_THEOR,&
                  &NR_FILES,NR_AFTER_CUSTOM_FLAGS,NR_AFTER_WDF,NR_WIND,NR_TS,NR_CO2,NR_H20,NR_CH4,NR_GS4,&
                  &NR_W_U,NR_W_TS,NR_W_CO2,NR_W_H2O,NR_W_CH4,NR_W_GS4,&
                  &TAU,H,LE,FC,FH2O,FCH4,FGS4,TAU_RU,H_RU,LE_RU,FC_RU,FH2O_RU,FCH4_RU,FGS4_RU,&
                  &H_STOR,LE_STOR,FC_STOR,FH2O_STOR,FCH4_STOR,FGS4_STOR,&
                  &FC_V_ADV,FH2O_V_ADV,FCH4_V_ADV,FGS4_V_ADV,&
                  &U_UNROT,V_UNROT,W_UNROT,U_ROT,V_ROT,W_ROT,&
                  &WS,MWS,WD,USTAR,TKE,L,ZL,BOWEN,TSTAR,&
                  &T_SONIC,TA,PA,RH,AIR_MV,AIR_DENSITY,AIR_RHO_CP,&
                  &VAPOR_DENSITY,VAPOR_PARTIAL_PRESSURE,VAPOR_PARTIAL_PRESSURE_SAT,SPECIFIC_HUMIDITY,VPD,TDEW,&
                  &DRYAIR_PARTIAL_PRESSURE,DRYAIR_DENSITY,DRYAIR_MV,SPECIFIC_HEAT_EVAP,SIGMA,&
                  &CO2_TYPE,CO2_MOLAR_DENSITY,CO2_MIXING_RATIO,CO2_MOLE_FRACTION,&
                  &H2O_TYPE,H2O_MOLAR_DENSITY,H2O_MIXING_RATIO,H2O_MOLE_FRACTION,&
                  &CH4_TYPE,CH4_MOLAR_DENSITY,CH4_MIXING_RATIO,CH4_MOLE_FRACTION,&
                  &GS4_TYPE,GS4_MOLAR_DENSITY,GS4_MIXING_RATIO,GS4_MOLE_FRACTION,&
                  &CO2_TLAG_ACTUAL,CO2_TLAG_USED,CO2_TLAG_NOMINAL,CO2_TLAG_MIN,CO2_TLAG_MAX,&
                  &H2O_TLAG_ACTUAL,H2O_TLAG_USED,H2O_TLAG_NOMINAL,H2O_TLAG_MIN,H2O_TLAG_MAX,&
                  &CH4_TLAG_ACTUAL,CH4_TLAG_USED,CH4_TLAG_NOMINAL,CH4_TLAG_MIN,CH4_TLAG_MAX,&
                  &GS4_TLAG_ACTUAL,GS4_TLAG_USED,GS4_TLAG_NOMINAL,GS4_TLAG_MIN,GS4_TLAG_MAX,&
                  &U_MEAN,V_MEAN,W_MEAN,TS_MEAN,CO2_MEAN,H2O_MEAN,CH4_MEAN,GS4_MEAN,&
                  &U_MEDIAN,V_MEAN,W_MEDIAN,TS_MEDIAN,CO2_MEDIAN,H2O_MEDIAN,CH4_MEDIAN,GS4_MEDIAN,&
                  &U_Q1,V_Q1,W_Q1,TS_Q1,CO2_Q1,H2O_Q1,CH4_Q1,GS4_Q1,&
                  &U_Q3,V_Q3,W_Q3,TS_Q3,CO2_Q3,H2O_Q3,CH4_Q3,GS4_Q3,&
                  &U_VAR,V_VAR,W_VAR,TS_VAR,CO2_VAR,H2O_VAR,CH4_VAR,GS4_VAR,&
                  &U_SKW,V_SKW,W_SKW,TS_SKW,CO2_SKW,H2O_SKW,CH4_SKW,GS4_SKW,&
                  &U_KUR,V_KUR,W_KUR,TS_KUR,CO2_KUR,H2O_KUR,CH4_KUR,GS4_KUR,&
                  &COV_W_U,COV_W_TS,COV_W_CO2,COV_W_H2O,COV_W_CH4,COV_W_GS4,&
                  &COV_CO2_H2O,COV_CO2_CH4,COV_CO2_GS4,COV_H2O_CH4,COV_H2O_GS4,COV_CH4_GS4,&
                  &FOOTPRINT_PEAK,FOOTPRINT_OFFSET,&
                  &FOOTPRINT_10,FOOTPRINT_30,FOOTPRINT_50,FOOTPRINT_70,FOOTPRINT_90,&
                  &TAU_0,H_0,LE_0,FC_0,FH2O_0,FCH4_0,FGS4_0,&
                  &TAU_1,H_1,LE_1,FC_1,FH2O_1,FCH4_1,FGS4_1,&
                  &TAU_2,H_2,LE_2,FC_2,FH2O_2,FCH4_2,FGS4_2,&
                  &CELLAIR_T,CELLAIR_P,CELLAIR_MV_CO2,CELLAIR_MV_H2O,CELLAIR_MV_CH4,CELLAIR_MV_GS4,&
                  &CELL_E_CO2,CELL_E_CH4,CELL_E_GS4,CELL_H_CO2,CELL_H_H2O,CELL_H_CH4,CELL_H_GS4,&
                  &BU_H_BOT,BU_H_TOP,BU_H_SPAR,LI7700_A,LI7700_B,LI7700_C,&
                  &TAU_SCF,H_SCF,LE_SCF,FC_SCF,FH2O_SCF,FCH4_SCF,FGS4_SCF,&
                  &COV_WT,COV_WT_1.626,COV_WT_0.614,COV_WT_0.277,COV_WT_0.133,COV_WT_0.065,&
                  &COV_WT_0.032,COV_WT_0.016,COV_WT_0.008,COV_WT_0.004,&
                  &M_CUSTOM_FLAGS,M_WDF,M_SONIC_DIAG,&
                  &M_IRGA_DIAG_CO2,M_IRGA_DIAG_H2O,M_IRGA_DIAG_CH4,M_IRGA_DIAG_GS4,&
                  &M_SPIKES_U,M_SPIKES_V,M_SPIKES_W,M_SPIKES_TS,&
                  &M_SPIKES_CO2,M_SPIKES_H2O,M_SPIKES_CH4,M_SPIKES_GS4,&
                  &M_ABSLIM_U,M_ABSLIM_V,M_ABSLIM_W,M_ABSLIM_TS,&
                  &M_ABSLIM_CO2,M_ABSLIM_H2O,M_ABSLIM_CH4,M_ABSLIM_GS4,&
                  &VM97_SPIKES_U,VM97_SPIKES_V,VM97_SPIKES_W,VM97_SPIKES_TS,&
                  &VM97_SPIKES_CO2,VM97_SPIKES_H2O,VM97_SPIKES_CH4,VM97_SPIKES_GS4,&
                  &VM97_AMPRES_U,VM97_AMPRES_V,VM97_AMPRES_W,VM97_AMPRES_TS,&
                  &VM97_AMPRES_CO2,VM97_AMPRES_H2O,VM97_AMPRES_CH4,VM97_AMPRES_GS4,&
                  &VM97_DOUT_C_U,VM97_DOUT_C_V,VM97_DOUT_C_W,VM97_DOUT_C_TS,&
                  &VM97_DOUT_C_CO2,VM97_DOUT_C_H2O,VM97_DOUT_C_CH4,VM97_DOUT_C_GS4,&
                  &VM97_DOUT_X_U,VM97_DOUT_X_V,VM97_DOUT_X_W,VM97_DOUT_X_TS,&
                  &VM97_DOUT_X_CO2,VM97_DOUT_X_H2O,VM97_DOUT_X_CH4,VM97_DOUT_X_GS4,&
                  &VM97_HM_SKW_U,VM97_HM_SKW_V,VM97_HM_SKW_W,VM97_HM_SKW_TS,&
                  &VM97_HM_SKW_CO2,VM97_HM_SKW_H2O,VM97_HM_SKW_CH4,VM97_HM_SKW_GS4,&
                  &VM97_HM_KUR_U,VM97_HM_KUR_V,VM97_HM_KUR_W,VM97_HM_KUR_TS,&
                  &VM97_HM_KUR_CO2,VM97_HM_KUR_H2O,VM97_HM_KUR_CH4,VM97_HM_KUR_GS4,&
                  &VM97_AOA,VM97_NSW_RNV1,VM97_NSW_RNV2,VM97_NSW_RNS,&
                  &VM97_SPIKES_HFLAG,VM97_ABSRES_HFLAG,VM97_DRPOUT_HFLAG,&
                  &VM97_ABSLIM_HFLAG,VM97_HGHMOM_HFLAG,VM97_HGHMOM_SFLAG,&
                  &VM97_DISCON_HFLAG,VM97_DISCON_SFLAG,&
                  &VM97_TIMELAG_HFLAG,VM97_TIMELAG_SFLAG,VM97_AOA_HFLAG,VM97_NSW_HFLAG,&
                  &FK04_ST_W_U,FK04_ST_W_TS,FK04_ST_W_CO2,&
                  &FK04_ST_W_H2O,FK04_ST_W_CH4,FK04_ST_W_GS4,&
                  &FK04_ITC_U,FK04_ITC_W,FK04_ITC_TS,&
                  &FK04_ST_FLAG_W_U,FK04_ST_FLAG_W_TS,FK04_ST_FLAG_W_CO2,&
                  &FK04_ST_FLAG_W_H2O,FK04_ST_FLAG_W_CH4,FK04_ST_FLAG_W_GS4,&
                  &FK04_ITC_FLAG_U,FK04_ITC_FLAG_W,FK04_ITC_FLAG_TS,&
                  &TAU_FK04_FLAG,H_FK04_FLAG,LE_FK04_FLAG,FC_FK04_FLAG,&
                  &FH2O_FK04_FLAG,FCH4_FK04_FLAG,FGS4_FK04_FLAG,&
                  &SPIKES_U,SPIKES_V,SPIKES_W,SPIKES_TS,SPIKES_CO2,&
                  &SPIKES_H2O,SPIKES_CH4,SPIKES_GS4,&
                  &LI7200_HEAD_DETECT,LI7200_T_OUT,LI7200_T_IN,LI7200_AUX_IN,&
                  &LI7200_DELTA_P,LI7200_CHOPPER,LI7200_DETECTOR,LI7200_PLL,LI7200_SYNC,&
                  &LI7500_CHOPPER,LI7500_DETECTOR,LI7500_PLL,LI7500_SYNC,&
                  &LI7700_NOT_READY,LI7700_NO_SIGNAL,LI7700_RE_UNLOCKED,&
                  &LI7700_BAD_TEMP,LI7700_LASER_T_UNREG,LI7700_CLOCK_T_UNREG,&
                  &LI7700_MOTOR_SPINNING,LI7700_PUMP_ON,LI7700_TOP_HEATER_ON,&
                  &LI7700_BOTTOM_HEATER_ON,LI7700_CALIBRATING,LI7700_MOTOR_FAILURE,&
                  &LI7700_BAD_AUX_TC1,LI7700_BAD_AUX_TC2,LI7700_BAD_AUX_TC3,LI7700_BOX_CONNECTED,&
                  &LI7200_AGC_OR_RSSI,LI7500_AGC_OR_RSSI,LI7700_RSSI,&
                  &WBOOST_APPLIED,AOA_METHOD,&
                  &AXES_ROTATION_METHOD,ROT_YAW,ROT_PITCH,ROT_ROLL,&
                  &DETRENDING_METHOD,DENTRENDING_TIME_CONSTANT,&
                  &TIMELAG_DETECTION_METHOD,WPL_APPLIED,BURBA_METHOD,&
                  &SPECTRAL_CORRECTION_METHOD,FOOTPRINT_MODEL,&
                  &LOGGER_SWVER_MAJOR,LOGGER_SWVER_MINOR,LOGGER_SWVER_REVISION,&
                  &LATITUDE,LONGITUDE,ALTITUDE,CANOPY_HEIGHT,DISPLACEMENT_HEIGHT,ROUGHNESS_LENGTH,&
                  &FILE_LENGTH,ACQUISITION_FREQUENCY,AVERAGING_INTERVAL,&
                  &MASTER_SONIC_MANUFACTURER,MASTER_SONIC_MODEL,MASTER_SONIC_HEIGHT,&
                  &MASTER_SONIC_WFORMAT,MASTER_SONIC_WREF,MASTER_SONIC_NORTH_OFFSET,&
                  &MASTER_SONIC_HPATH,MASTER_SONIC_VPATH,MASTER_SONIC_TAU,&
                  &CO2_IRGA_MANUFACTURER,CO2_IRGA_MODEL,&
                  &CO2_IRGA_NORTHWARD_SEPARATION,CO2_IRGA_EASTWARD_SEPARATION,CO2_IRGA_VERTICAL_SEPARATION,&
                  &CO2_IRGA_TUBE_LENGTH,CO2_IRGA_TUBE_DIAMETER,CO2_IRGA_TUBE_FLOWRATE,&
                  &CO2_IRGA_KW,CO2_IRGA_KO,CO2_IRGA_HPATH_LENGTH,CO2_IRGA_VPATH_LENGTH,CO2_IRGA_TAU,&
                  &H2O_IRGA_MANUFACTURER,H2O_IRGA_MODEL,&
                  &H2O_IRGA_NORTHWARD_SEPARATION,H2O_IRGA_EASTWARD_SEPARATION,H2O_IRGA_VERTICAL_SEPARATION,&
                  &H2O_IRGA_TUBE_LENGTH,H2O_IRGA_TUBE_DIAMETER,H2O_IRGA_TUBE_FLOWRATE,&
                  &H2O_IRGA_KW,H2O_IRGA_KO,H2O_IRGA_HPATH_LENGTH,H2O_IRGA_VPATH_LENGTH,H2O_IRGA_TAU,&
                  &CH4_IRGA_MANUFACTURER,CH4_IRGA_MODEL,&
                  &CH4_IRGA_NORTHWARD_SEPARATION,CH4_IRGA_EASTWARD_SEPARATION,CH4_IRGA_VERTICAL_SEPARATION,&
                  &CH4_IRGA_TUBE_LENGTH,CH4_IRGA_TUBE_DIAMETER,CH4_IRGA_TUBE_FLOWRATE,&
                  &CH4_IRGA_KW,CH4_IRGA_KO,CH4_IRGA_HPATH_LENGTH,CH4_IRGA_VPATH_LENGTH,CH4_IRGA_TAU,&
                  &GS4_IRGA_MANUFACTURER,GS4_IRGA_MODEL,&
                  &GS4_IRGA_NORTHWARD_SEPARATION,GS4_IRGA_EASTWARD_SEPARATION,GS4_IRGA_VERTICAL_SEPARATION,&
                  &GS4_IRGA_TUBE_LENGTH,GS4_IRGA_TUBE_DIAMETER,GS4_IRGA_TUBE_FLOWRATE,&
                  &GS4_IRGA_KW,GS4_IRGA_KO,GS4_IRGA_HPATH_LENGTH,GS4_IRGA_VPATH_LENGTH,GS4_IRGA_TAU,'

        !> Add custom variables
        call AddDatum(dataline, 'NUM_CUSTOM_VARS', separator)
        if (NumUserVar > 0) then
            do i = 1, NumUserVar
                dataline = dataline(1:len_trim(dataline)) &
                    // usg(i)(1:len_trim(usg(i))) // 'mean' // ','
            end do
        end if

        !> Add biomet variables
        call AddDatum(dataline, 'NUM_BIOMET_VARS', separator)
        if (nbVars > 0) then
            do i = 1, nbVars
                call AddDatum(dataline, trim(bVars(i)%label), separator)
            end do
        end if

        write(uicos, '(a)') dataline(1:len_trim(dataline) - 1)
    end if

    !>==========================================================================
    !>==========================================================================
    !> FLUXNET EDDY output

    !> Even if it was selected for output, FLUXNET eddy is not created if
    !> FCC will run
    if (EddyProProj%fcc_follows) EddyProProj%out_fluxnet_eddy = .false.

    if (EddyProProj%out_fluxnet_eddy) then
        Test_Path = Dir%main_out(1:len_trim(Dir%main_out)) &
                  // EddyProProj%id(1:len_trim(EddyProProj%id)) &
                  // FLUXNET_EDDY_FilePadding // Timestamp_FilePadding // CsvExt
        dot = index(Test_Path, CsvExt, .true.) - 1
        FLUXNET_EDDY_Path = Test_Path(1:dot) // CsvTmpExt
        open(ufnet_e, file = FLUXNET_EDDY_Path, &
            iostat = open_status, encoding = 'utf-8')

        !> Initialize header strings to void
        call Clearstr(header2)
        call Clearstr(header3)
        call Clearstr(head2_utf8)
        call Clearstr(head3_utf8)

        !> Initial common part
        call AddDatum(header2,'TIMESTAMP', separator)
        call AddDatum(header3,'[yyyymmddHHMMSS]', separator)

        !>======================================================================
        !> Variables located at the EC stations
        !> GASES
        do gas = co2, gas4
            if(OutVarPresent(gas)) then
                select case (gas)
                    case(co2)
                        call AddDatum(header2, 'CO2_1_1_1', separator)
                        call AddDatum(header3, &
                            '[mmolCO2 mol-1]', separator)
                    case(h2o)
                        call AddDatum(header2, 'H2O_1_1_1', separator)
                        call AddDatum(header3, '[mmolH2O mol-1]', separator)
                    case(ch4)
                        call AddDatum(header2, 'CH4_1_1_1', separator)
                        call AddDatum(header3, '[nmolCH4 mol-1]', separator)
                    case(gas4)
                        call AddDatum(header2, &
                            trim(adjustl(e2sg(gas4))) // '1_1_1', separator)
                        call AddDatum(header3, '[nmol' &
                            // trim(adjustl(e2sg(gas4))) //' mol-1]', separator)
                end select
            end if
        end do
        call AddDatum(header2,&
            'TAU_1_1_1,TAU_SSITC_TEST_1_1_1,H_1_1_1,H_SSITC_TEST_1_1_1', separator)
        call AddDatum(header3,'[kg m-1 s-2],[#],[W m-2],[#]', separator)
        if(OutVarPresent(h2o)) then
            call AddDatum(header2,'LE_1_1_1,LE_SSITC_TEST_1_1_1', separator)
            call AddDatum(header3,'[W m-2],[#]', separator)
        end if
        if(OutVarPresent(co2)) then
            call AddDatum(header2,'FC_1_1_1,FC_SSITC_TEST_1_1_1', separator)
            call AddDatum(header3,'[umolCO2 m-2 s-1],[#]', separator)
        end if
        if(OutVarPresent(ch4)) then
            call AddDatum(header2,'FCH4_1_1_1,FCH4_SSITC_TEST_1_1_1', separator)
            call AddDatum(header3,'[nmolCH4 m-2 s-1],[#]', separator)
        end if
        if(OutVarPresent(gas4)) then
            call AddDatum(header2,'F' // trim(adjustl(e2sg(gas4))) &
                // '1_1_1,' // trim(adjustl(e2sg(gas4))) // 'SSITC_TEST_1_1_1', separator)
            call AddDatum(header3,'[nmol'// trim(adjustl(e2sg(gas4))) &
                //' m-2 s-1],[#]', separator)
        end if

        !> MET_WIND and FOOTPRINT
        call AddDatum(header2,'WD_0_0_1,WS_0_0_1,WS_MAX_0_0_1,U_SIGMA_0_0_1,&
            &V_SIGMA_0_0_1,W_SIGMA_0_0_1,USTAR_0_0_1,MO_LENGTH_0_0_1,ZL_0_0_1,&
            &FETCH_MAX_0_0_1,FETCH_70_0_0_1,FETCH_80_0_0_1,FETCH_90_0_0_1', separator)
        call AddDatum(header3,'[Decimal degrees],[m s-1],[m s-1],[m s-1],&
            &[m s-1],[m s-1],[m s-1],[m],[#],[m],[m],[m],[m]', separator)

        !> MET_ATM
        call AddDatum(header2,'PA_0_0_1,RH_0_0_1,TA_0_0_1,VPD_0_0_1,T_SONIC_0_0_1,&
            &T_SONIC_SIGMA_0_0_1', separator)
        call AddDatum(header3,'[kPa],[%],[deg C],[hPa],[deg C],[deg C]', separator)

        call latin1_to_utf8(header2, head2_utf8)
        call latin1_to_utf8(header3, head3_utf8)

        !> Write on output file
        write(ufnet_e, '(a)') head2_utf8(1:len_trim(head2_utf8) - 1)
        write(ufnet_e, '(a)') head3_utf8(1:len_trim(head3_utf8) - 1)
    end if

    !>==========================================================================
    !>==========================================================================
    !> AmeriFlux output
    if (EddyProProj%out_amflux) then
        Test_Path = Dir%main_out(1:len_trim(Dir%main_out)) &
                  // EddyProProj%id(1:len_trim(EddyProProj%id)) &
                  // AmeriFlux_FilePadding // Timestamp_FilePadding // CsvExt
        dot = index(Test_Path, CsvExt, .true.) - 1
        AmeriFlux_Path = Test_Path(1:dot) // CsvTmpExt
        open(uaflx, file = AmeriFlux_Path, iostat = open_status, encoding = 'utf-8')

        write(uaflx, '(a)') 'Sitename:' // Metadata%sitename(1:len_trim(Metadata%sitename))
        write(uaflx, '(a, f12.7, a, f12.7, a, f6.0)') 'Location: Latitude: ', Metadata%lat, &
            ' - Longitude: ', Metadata%lon, ' - Elevation (masl): ', Metadata%alt
        write(uaflx, '(a)') 'Principal investigator: '
        write(uaflx, '(a)') 'Ecosystem type: '
        call idate(today)   ! today(1)=day, (2)=month, (3)=year
        call itime(now)     ! now(1)=hour, (2)=minute, (3)=second
        call clearstr(dataline)
        write(dum_string, '(i4)') today(3)
        dataline(1:5) = dum_string(1:4) // '-'
        write(dum_string, '(i2)') today(2)
        if(dum_string(1:1) == ' ') dum_string(1:1) = '0'
        dataline(6:8) = dum_string(1:2) // '-'
        write(dum_string, '(i2)') today(1)
        if(dum_string(1:1) == ' ') dum_string(1:1) = '0'
        dataline(9:11) = dum_string(1:2) // 'T'
        write(dum_string, '(i2)') now(1)
        if(dum_string(1:1) == ' ') dum_string(1:1) = '0'
        dataline(12:14) = dum_string(1:2) // ':'
        write(dum_string, '(i2)') now(2)
        if(dum_string(1:1) == ' ') dum_string(1:1) = '0'
        dataline(15:16) = dum_string(1:2)
        write(uaflx, '(a)') 'File creation date: ' // dataline(1:len_trim(dataline))
        write(uaflx, '(a)') 'Datapolicy:  -- The AmeriFlux data provided on this site&
            & are freely available and were furnished by individual AmeriFlux scientists who encourage their use.'
        write(uaflx, '(a)') 'Please kindly inform in writing (or e-mail) the appropriate AmeriFlux scientist(s)&
            & of how you intend to use the data and of any publication plans.'
        write(uaflx, '(a)') 'It is also important to contact the AmeriFlux investigator to assure you are downloading&
            & the latest revision of the data and to prevent potential misuse or misinterpretation of the data.'
        write(uaflx, '(a)') 'Please acknowledge the data source as a citation or in the&
            & acknowledgments if no citation is available.'
        write(uaflx, '(a)') 'If the AmeriFlux Principal Investigators (PIs) feel that they should be acknowledged&
            & or offered participation as authors they will let you know.'
        write(uaflx, '(a)') 'And we assume that an agreement on such matters will be reached&
            & before publishing and/or use of the data for publication.'
        write(uaflx, '(a)') 'If your work directly competes with the PIs analysis they may ask that they have the&
            & opportunity to submit a manuscript before you submit one that uses unpublished data.'
        write(uaflx, '(a)') 'In addition when publishing please acknowledge the agency that supported the research. --'
        write(uaflx, '(a)') 'File Origin (4 lines) - To be compiled by AmeriFlux data management group.'
        write(uaflx, '(a)') 'File Origin (4 lines) - To be compiled by AmeriFlux data management group.'
        write(uaflx, '(a)') 'File Origin (4 lines) - To be compiled by AmeriFlux data management group.'
        write(uaflx, '(a)') 'File Origin (4 lines) - To be compiled by AmeriFlux data management group.'
        write(uaflx, '(a)') 'YEAR,GAP,DTIME,DOY,HRMIN,UST,TA,WD,WS,NEE,FC,SFC,H,SH,LE,SLE,FG,TS1,TSdepth1,TS2,TSdepth2,&
            &PREC,RH,PRESS,CO2,VPD,SWC1,SWC2,Rn,PAR,Rg,Rgdif,PARout,RgOut,Rgl,RglOut,&
            &H2O,RE,GPP,CO2top,CO2height,APAR,PARdif,APARpct,ZL'
        write(uaflx, '(a)') 'YEAR,GAP,DTIME,DOY,HRMIN,m/s,deg C,deg,m/s,umol/m2/s,umol/m2/s,&
            &umol/m2/s,W/m2,W/m2,W/m2,W/m2,W/m2,&
            &deg C,cm,deg C,cm,mm,%,kPa,umol/mol,kPa,%,%,W/m2,umol/m2/s,W/m2,W/m2,umol/m2/s,&
            &W/m2,W/m2,W/m2,mmol/mol,umol/m2/s,&
            &umol/m2/s,umol/mol,m,umol/m2/s,umol/m2/s,%,unitless'
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
