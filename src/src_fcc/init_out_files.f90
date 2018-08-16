!***************************************************************************
! init_out_files.f90
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
! \brief       Initializes output files
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine InitOutFiles(lEx)
    use m_fx_global_var
    implicit none
    !> in/out variables
    Type(ExType), intent(in) :: lEx
    !> local variables
    integer :: mkdir_status
    integer :: open_status
    integer :: dot
    integer :: gas
    integer :: i
    character(PathLen) :: Test_Path
    character(64) :: e2sg(E2NumVar)
    character(LongOutstringLen) :: header1
    character(LongOutstringLen) :: header2
    character(LongOutstringLen) :: header3
    character(LongOutstringLen) :: head1_utf8
    character(LongOutstringLen) :: head2_utf8
    character(LongOutstringLen) :: head3_utf8
    integer, external :: CreateDir


    e2sg(u)    = 'u_'
    e2sg(v)    = 'v_'
    e2sg(w)    = 'w_'
    e2sg(ts)   = 'ts_'
    e2sg(co2)  = 'co2_'
    e2sg(h2o)  = 'h2o_'
    e2sg(ch4)  = 'ch4_'
    e2sg(gas4) = g4lab(1:g4l) // '_'

    !> Full output file
    if (EddyProProj%out_full) then
        !> Create output directory if it does not exist
        mkdir_status = CreateDir('"' // Dir%main_out(1:len_trim(Dir%main_out)) // '"')

        !> Open full output file and writes header
        Test_Path = Dir%main_out(1:len_trim(Dir%main_out)) &
                  // EddyProProj%id(1:len_trim(EddyProProj%id)) &
                  // FullOut_FilePadding // Timestamp_FilePadding // CsvExt
        dot = index(Test_Path, CsvExt, .true.) - 1
        FullOut_Path = Test_Path(1:dot) // CsvTmpExt
        open(uflx, file = FullOut_Path, iostat = open_status, encoding = 'utf-8')

        !> Initialize header strings to void
        call clearstr(header1)
        call clearstr(header2)
        call clearstr(header3)
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
            if(fcc_var_present(h2o)) then
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
            if(fcc_var_present(co2)) then
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
            if(fcc_var_present(h2o)) then
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
            if(fcc_var_present(ch4)) then
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
            if(fcc_var_present(gas4)) then
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
            if(fcc_var_present(h2o)) call AddDatum(header1, '', separator)
            if(fcc_var_present(h2o)) call AddDatum(header2,'LE_strg', separator)
            if(fcc_var_present(h2o)) call AddDatum(header3,'[W+1m-2]', separator)
            do gas = co2, gas4
                if (gas /= h2o) then
                    if(fcc_var_present(gas)) call AddDatum(header1, '', separator)
                    if(fcc_var_present(gas)) call AddDatum(header2, e2sg(gas)(1:len_trim(e2sg(gas))) // 'strg', separator)
                    if(fcc_var_present(gas)) call AddDatum(header3, '[' // char(181) // 'mol+1s-1m-2]', separator)
                else
                    if(fcc_var_present(gas)) call AddDatum(header1, '', separator)
                    if(fcc_var_present(gas)) call AddDatum(header2, e2sg(gas)(1:len_trim(e2sg(gas))) // 'strg', separator)
                    if(fcc_var_present(gas)) call AddDatum(header3, '[mmol+1s-1m-2]', separator)
                end if
            end do

            !> Advection fluxes
            header1 = header1(1:len_trim(header1)) // 'vertical_advection_fluxes'
            do gas = co2, n2o
                if (gas /= h2o) then
                    if(fcc_var_present(gas)) call AddDatum(header1, '', separator)
                    if(fcc_var_present(gas)) call AddDatum(header2, e2sg(gas)(1:len_trim(e2sg(gas))) // 'v-adv', separator)
                    if(fcc_var_present(gas)) call AddDatum(header3, '[' // char(181) // 'mol+1s-1m-2]', separator)
                else
                    if(fcc_var_present(gas)) call AddDatum(header1, '', separator)
                    if(fcc_var_present(gas)) call AddDatum(header2, e2sg(gas)(1:len_trim(e2sg(gas))) // 'v-adv', separator)
                    if(fcc_var_present(gas)) call AddDatum(header3, '[mmol+1s-1m-2]', separator)
                end if
            end do

            !> Average gas concentrations
            call AddDatum(header1,'gas_densities_concentrations_and_timelags', separator)
            do gas = co2, gas4
                if(fcc_var_present(gas)) call AddDatum(header1, ',,,,', separator)
                if(fcc_var_present(gas)) call AddDatum(header2, e2sg(gas)(1:len_trim(e2sg(gas))) // 'molar_density,' &
                    // e2sg(gas)(1:len_trim(e2sg(gas))) // 'mole_fraction,' &
                    // e2sg(gas)(1:len_trim(e2sg(gas))) // 'mixing_ratio,' &
                    // e2sg(gas)(1:len_trim(e2sg(gas))) // 'time_lag,' &
                    // e2sg(gas)(1:len_trim(e2sg(gas))) // 'def_timelag', separator)
                if (gas /= h2o) then
                    if(fcc_var_present(gas)) call AddDatum(header3, '[mmol+1m-3],[' // char(181) // &
                        'mol+1mol_a-1],[' // char(181) // 'mol+1mol_d-1],[s],[1=default]', separator)
                else
                    if(fcc_var_present(gas)) &
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
                          &[m+1s-1],[m+1s-1],[m+1s-1],[m+1s-1],[m+1s-1],&
                          &[m+1s-1],[m+1s-1],[m+1s-1],[deg_from_north],[deg],[deg],[deg]', separator)

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
            if(fcc_var_present(h2o)) call AddDatum(header1, ',', separator)
            if(fcc_var_present(h2o)) call AddDatum(header2,'un_LE,LE_scf', separator)
            if(fcc_var_present(h2o)) call AddDatum(header3,'[W+1m-2],[#]', separator)
            !> Uncorrected gas fluxes (Level 0)
            do gas = co2, gas4
                if (gas /= h2o) then
                    if(fcc_var_present(gas)) call AddDatum(header1, ',', separator)
                    if(fcc_var_present(gas)) call AddDatum(header2, 'un_' // e2sg(gas)(1:len_trim(e2sg(gas))) &
                        // 'flux,' // e2sg(gas)(1:len_trim(e2sg(gas))) // 'scf', separator)
                    if(fcc_var_present(gas)) call AddDatum(header3, '[' // char(181) // 'mol+1s-1m-2],[#]', separator)
                else
                    if(fcc_var_present(gas)) call AddDatum(header1, ',', separator)
                    if(fcc_var_present(gas)) call AddDatum(header2, 'un_' // e2sg(gas)(1:len_trim(e2sg(gas))) &
                        // 'flux,' // e2sg(gas)(1:len_trim(e2sg(gas))) // 'scf', separator)
                    if(fcc_var_present(gas)) call AddDatum(header3, '[mmol+1s-1m-2],[#]', separator)
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
            do gas = co2, gas4
                if(fcc_var_present(gas)) then
                    call AddDatum(header1, '', separator)
                    call AddDatum(header2, e2sg(gas)(1:len_trim(e2sg(gas))) // 'spikes' , separator)
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
                if(co2_new_sw_ver) then
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
                if(co2_new_sw_ver) then
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
                if(fcc_var_present(gas)) call AddDatum(header1, '', separator)
                if(fcc_var_present(gas)) call AddDatum(header2, e2sg(gas)(1:len_trim(e2sg(gas))) // 'var', separator)
                if(fcc_var_present(gas)) call AddDatum(header3, '--', separator)
            end do
            !> w/ts covariance
            call AddDatum(header1, 'covariances', separator)
            call AddDatum(header2,'w/ts_cov', separator)
            call AddDatum(header3,'[m+1K+1s-1]', separator)
            !> w-gases covariances
            do gas = co2, gas4
                if(fcc_var_present(gas)) call AddDatum(header1, '', separator)
                if(fcc_var_present(gas)) call AddDatum(header2, 'w/' // e2sg(gas)(1:len_trim(e2sg(gas))) // 'cov', separator)
                if(fcc_var_present(gas)) call AddDatum(header3, '--', separator)
            end do

            !> Mean values of user variables
            if (lEx%ncustom > 0) then
                call AddDatum(header1, 'custom_variables', separator)
                call AddDatum(header2, UserVarHeader(1:len_trim(UserVarHeader)), separator)
                do i = 1, lEx%ncustom
                    call AddDatum(header3, '--', separator)
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
                if(co2_new_sw_ver) then
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
            if (lEx%ncustom > 0) then
                call AddDatum(header1, 'custom_variables', separator)
                call AddDatum(header2, UserVarHeader(1:len_trim(UserVarHeader)), separator)
                do i = 1, lEx%ncustom
                    call AddDatum(header3, '--', separator)
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

    !************************************************************************************************************************************
    !************************************************************************************************************************************
    !> METADATA file
    if (EddyProProj%out_md) then
        !> Create metadata output file name
        !> Open dynamic matadata file
        Test_Path = Dir%main_out(1:len_trim(Dir%main_out)) &
                  // EddyProProj%id(1:len_trim(EddyProProj%id)) &
                  // MetaData_FilePadding // Timestamp_FilePadding // CsvExt
        dot = index(Test_Path, CsvExt, .true.) - 1
        Metadata_Path = Test_Path(1:dot) // CsvTmpExt
        open(umd, file = Metadata_Path, iostat = open_status, encoding = 'utf-8')
        call Clearstr(header1)
        call AddDatum(header1,'filename,date,time,latitude,longitude,altitude,canopy_height,displacement_height,&
            &roughness_length,file_length,acquisition_frequency,&
            &master_sonic_manufacturer,master_sonic_model,master_sonic_height,&
            &master_sonic_wformat,master_sonic_wref,master_sonic_north_offset,&
            &master_sonic_hpath_length,master_sonic_vpath_length,master_sonic_tau', separator)
        if (fcc_var_present(co2)) &
            call AddDatum(header1,'co2_irga_manufacturer,co2_irga_model,co2_measure_type,co2_irga_northward_separation,&
                &co2_irga_eastward_separation,co2_irga_vertical_separation,&
                &co2_irga_tube_length,co2_irga_tube_diameter,co2_irga_tube_flowrate,&
                &co2_irga_kw,co2_irga_ko,co2_irga_hpath_length,co2_irga_vpath_length,co2_irga_tau', separator)
        if (fcc_var_present(h2o)) &
            call AddDatum(header1,'h2o_irga_manufacturer,h2o_irga_model,h2o_measure_type,h2o_irga_northward_separation,&
                &h2o_irga_eastward_separation,h2o_irga_vertical_separation,&
                &h2o_irga_tube_length,h2o_irga_tube_diameter,h2o_irga_tube_flowrate,&
                &h2o_irga_kw,h2o_irga_ko,h2o_irga_hpath_length,h2o_irga_vpath_length,h2o_irga_tau', separator)
        if (fcc_var_present(ch4)) &
            call AddDatum(header1,'ch4_irga_manufacturer,ch4_irga_model,ch4_measure_type,ch4_irga_northward_separation,&
                &ch4_irga_eastward_separation,ch4_irga_vertical_separation,&
                &ch4_irga_tube_length,ch4_irga_tube_diameter,ch4_irga_tube_flowrate,&
                &ch4_irga_kw,ch4_irga_ko,ch4_irga_hpath_length,ch4_irga_vpath_length,ch4_irga_tau', separator)
        if (fcc_var_present(gas4)) &
            call AddDatum(header1, e2sg(gas4)(1:len_trim(e2sg(gas4))) // 'irga_manufacturer,' &
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
        write(umd, '(a)') header1(1:len_trim(header1) - 1)
    end if

    !***************************************************************************
    !***************************************************************************

    if (EddyProProj%out_fluxnet) then
        !> Create output directory if it does not exist
        mkdir_status = CreateDir('"' // Dir%main_out(1:len_trim(Dir%main_out)) // '"')

        Test_Path = Dir%main_out(1:len_trim(Dir%main_out)) &
                  // EddyProProj%id(1:len_trim(EddyProProj%id)) &
                  // FLUXNET_FilePadding // Timestamp_FilePadding // CsvExt
        dot = index(Test_Path, CsvExt, .true.) - 1
        FLUXNET_Path = Test_Path(1:dot) // CsvTmpExt
        open(uflxnt, file = FLUXNET_Path, iostat = open_status, encoding = 'utf-8')

        !> Write on output file
        write(uflxnt, '(a)') trim(fluxnet_header)
    end if

end subroutine InitOutFiles
