!***************************************************************************
! init_outfiles_rp.f90
! --------------------
! Copyright (C) 2019-2019, LI-COR Biosciences
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
! \brief       Initializes Fluxnet output file
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine InitFluxnetFile_rp()
    use m_rp_global_var
    use iso_fortran_env
    implicit none
    !> in/out variables
    integer, external :: CreateDir
    !> local variables
    integer :: open_status = 1      ! initializing to false
    integer :: dot
    integer :: i
    character(PathLen) :: Test_Path
    character(64) :: e2sg(E2NumVar)
    character(32) :: usg(NumUserVar)
    character(LongOutstringLen) :: dataline
    include '../src_common/interfaces.inc'


    Test_Path = Dir%main_out(1:len_trim(Dir%main_out)) &
                // EddyProProj%id(1:len_trim(EddyProProj%id)) &
                // FLUXNET_FilePadding // Timestamp_FilePadding // CsvExt
    dot = index(Test_Path, CsvExt, .true.) - 1
    FLUXNET_Path = Test_Path(1:dot) // CsvTmpExt
    open(uflxnt, file = FLUXNET_Path, iostat = open_status, encoding = 'utf-8')

    call clearstr(dataline)
    dataline = 'TIMESTAMP_START,TIMESTAMP_END,DOY_START,DOY_END,FILENAME_HF,SW_IN_POT,NIGHT,EXPECT_NR,&
                &FILE_NR,CUSTOM_FILTER_NR,WD_FILTER_NR,SONIC_NR,T_SONIC_NR,CO2_NR,H2O_NR,CH4_NR,GS4_NR,&
                &TAU_NR,H_NR,FC_NR,LE_NR,FCH4_NR,FGS4_NR,&
                &TAU,H,LE,ET,FC,FH2O,FCH4,FGS4,TAU_RANDUNC_HF,H_RANDUNC_HF,LE_RANDUNC_HF,ET_RANDUNC_HF,&
                &FC_RANDUNC_HF,FH2O_RANDUNC_HF,FCH4_RANDUNC_HF,FGS4_RANDUNC_HF,&
                &SH_SINGLE,SLE_SINGLE,SET_SINGLE,SC_SINGLE,SH2O_SINGLE,SCH4_SINGLE,SGS4_SINGLE,&
                &FC_VADV,FH2O_VADV,FCH4_VADV,FGS4_VADV,&
                &U_UNROT,V_UNROT,W_UNROT,U,V,W,&
                &WS,WS_MAX,WD,WD_SIGMA,USTAR,TKE,MO_LENGTH,ZL,BOWEN,TSTAR,&
                &T_SONIC,TA_EP,PA_EP,RH_EP,AIR_MV,AIR_DENSITY,AIR_RHO_CP,&
                &VAPOR_DENSITY,VAPOR_PARTIAL_PRESSURE,VAPOR_PARTIAL_PRESSURE_SAT,SPECIFIC_HUMIDITY,VPD_EP,TDEW,&
                &DRYAIR_PARTIAL_PRESSURE,DRYAIR_DENSITY,DRYAIR_MV,SPECIFIC_HEAT_EVAP,VAPOR_DRYAIR_RATIO,&
                &CO2_MEAS_TYPE,CO2_MOLAR_DENSITY,CO2_MIXING_RATIO,CO2,&
                &H2O_MEAS_TYPE,H2O_MOLAR_DENSITY,H2O_MIXING_RATIO,H2O,&
                &CH4_MEAS_TYPE,CH4_MOLAR_DENSITY,CH4_MIXING_RATIO,CH4,&
                &GS4_MEAS_TYPE,GS4_MOLAR_DENSITY,GS4_MIXING_RATIO,GS4,&
                &CO2_TLAG_ACTUAL,CO2_TLAG_USED,CO2_TLAG_NOMINAL,CO2_TLAG_MIN,CO2_TLAG_MAX,&
                &H2O_TLAG_ACTUAL,H2O_TLAG_USED,H2O_TLAG_NOMINAL,H2O_TLAG_MIN,H2O_TLAG_MAX,&
                &CH4_TLAG_ACTUAL,CH4_TLAG_USED,CH4_TLAG_NOMINAL,CH4_TLAG_MIN,CH4_TLAG_MAX,&
                &GS4_TLAG_ACTUAL,GS4_TLAG_USED,GS4_TLAG_NOMINAL,GS4_TLAG_MIN,GS4_TLAG_MAX,&
                &U_MEDIAN,V_MEDIAN,W_MEDIAN,T_SONIC_MEDIAN,&
                &CO2_MEAS_MEDIAN,H2O_MEAS_MEDIAN,CH4_MEAS_MEDIAN,GS4_MEAS_MEDIAN,&
                &U_P25,V_P25,W_P25,T_SONIC_P25,CO2_MEAS_P25,H2O_MEAS_P25,CH4_MEAS_P25,GS4_MEAS_P25,&
                &U_P75,V_P75,W_P75,T_SONIC_P75,CO2_MEAS_P75,H2O_MEAS_P75,CH4_MEAS_P75,GS4_MEAS_P75,&
                &U_SIGMA,V_SIGMA,W_SIGMA,T_SONIC_SIGMA,CO2_MEAS_SIGMA,H2O_MEAS_SIGMA,CH4_MEAS_SIGMA,GS4_MEAS_SIGMA,&
                &U_SKW,V_SKW,W_SKW,T_SONIC_SKW,CO2_MEAS_SKW,H2O_MEAS_SKW,CH4_MEAS_SKW,GS4_MEAS_SKW,&
                &U_KUR,V_KUR,W_KUR,T_SONIC_KUR,CO2_MEAS_KUR,H2O_MEAS_KUR,CH4_MEAS_KUR,GS4_MEAS_KUR,&
                &W_U_COV,W_T_SONIC_COV,W_CO2_MEAS_COV,W_H2O_MEAS_COV,W_CH4_MEAS_COV,W_GS4_MEAS_COV,&
                &CO2_MEAS_H2O_MEAS_COV,CO2_MEAS_CH4_MEAS_COV,CO2_MEAS_GS4_MEAS_COV,&
                &H2O_MEAS_CH4_MEAS_COV,H2O_MEAS_GS4_MEAS_COV,CH4_MEAS_GS4_MEAS_COV,&
                &FETCH_MAX,FETCH_OFFSET,FETCH_10,FETCH_30,FETCH_50,FETCH_70,FETCH_80,FETCH_90,&
                &MO_LENGTH_UNCORR,ZL_UNCORR,&
                &TAU_UNCORR,H_UNCORR,LE_UNCORR,ET_UNCORR,FC_UNCORR,FH2O_UNCORR,FCH4_UNCORR,FGS4_UNCORR,&
                &TAU_STAGE1,H_STAGE1,LE_STAGE1,ET_STAGE1,FC_STAGE1,FH2O_STAGE1,FCH4_STAGE1,FGS4_STAGE1,&
                &TAU_STAGE2,H_STAGE2,LE_STAGE2,ET_STAGE2,FC_STAGE2,FH2O_STAGE2,FCH4_STAGE2,FGS4_STAGE2,&
                &T_CELL,PA_CELL,MV_AIR_CELL_CO2,MV_AIR_CELL_H2O,MV_AIR_CELL_CH4,MV_AIR_CELL_GS4,&
                &FH2O_CELL_CO2,FH2O_CELL_CH4,FH2O_CELL_GS4,H_CELL_CO2,H_CELL_H2O,H_CELL_CH4,H_CELL_GS4,&
                &H_BU_BOT,H_BU_TOP,H_BU_SPAR,SPEC_CORR_LI7700_A,SPEC_CORR_LI7700_B,SPEC_CORR_LI7700_C,&
                &TAU_SCF,H_SCF,LE_SCF,ET_SCF,FC_SCF,FH2O_SCF,FCH4_SCF,FGS4_SCF,&
                &W_T_SONIC_COV_IBROM,W_T_SONIC_COV_IBROM_N1626,W_T_SONIC_COV_IBROM_N0614,&
                &W_T_SONIC_COV_IBROM_N0277,W_T_SONIC_COV_IBROM_N0133,&
                &W_T_SONIC_COV_IBROM_N0065,W_T_SONIC_COV_IBROM_N0032,&
                &W_T_SONIC_COV_IBROM_N0016,W_T_SONIC_COV_IBROM_N0008,W_T_SONIC_COV_IBROM_N0004,&
                &CUSTOM_FILTER_NREX,WD_FILTER_NREX,SONIC_DIAG_NREX,&
                &CO2_DIAG_NREX,H2O_DIAG_NREX,CH4_DIAG_NREX,GS4_DIAG_NREX,&
                &U_SPIKE_NREX,V_SPIKE_NREX,W_SPIKE_NREX,T_SONIC_SPIKE_NREX,&
                &CO2_SPIKE_NREX,H2O_SPIKE_NREX,CH4_SPIKE_NREX,GS4_SPIKE_NREX,&
                &U_ABSLIM_NREX,V_ABSLIM_NREX,W_ABSLIM_NREX,T_SONIC_ABSLIM_NREX,&
                &CO2_ABSLIM_NREX,H2O_ABSLIM_NREX,CH4_ABSLIM_NREX,GS4_ABSLIM_NREX,&
                &U_VM97_TEST,V_VM97_TEST,W_VM97_TEST,T_SONIC_VM97_TEST,&
                &CO2_VM97_TEST,H2O_VM97_TEST,CH4_VM97_TEST,GS4_VM97_TEST,&
                &VM97_TLAG_HF,VM97_TLAG_SF,VM97_AOA_HF,VM97_NSHW_HF,&
                &U_KID,V_KID,W_KID,T_SONIC_KID,CO2_KID,H2O_KID,CH4_KID,GS4_KID,&
                &U_ZCD,V_ZCD,W_ZCD,T_SONIC_ZCD,CO2_ZCD,H2O_ZCD,CH4_ZCD,GS4_ZCD,&
                &TAU_CORRDIFF,H_CORRDIFF,LE_CORRDIFF,ET_CORRDIFF,FC_CORRDIFF,&
                &FH2O_CORRDIFF,FCH4_CORRDIFF,FGS4_CORRDIFF,&
                &TAU_NSR,H_NSR,FC_NSR,FH2O_NSR,FCH4_NSR,FGS4_NSR,&
                &TAU_SS,H_SS,FC_SS,FH2O_SS,FCH4_SS,FGS4_SS,&
                &U_ITC,W_ITC,T_SONIC_ITC,TAU_SS_TEST,H_SS_TEST,FC_SS_TEST,&
                &FH2O_SS_TEST,FCH4_SS_TEST,FGS4_SS_TEST,&
                &U_ITC_TEST,W_ITC_TEST,T_SONIC_ITC_TEST,&
                &TAU_SSITC_TEST,H_SSITC_TEST,LE_SSITC_TEST,ET_SSITC_TEST,FC_SSITC_TEST,&
                &FH2O_SSITC_TEST,FCH4_SSITC_TEST,FGS4_SSITC_TEST,&
                &INST_LI7200_HEAD_DETECT,INST_LI7200_T_OUT,INST_LI7200_T_IN,INST_LI7200_AUX_IN,&
                &INST_LI7200_DELTA_P,INST_LI7200_CHOPPER,INST_LI7200_DETECTOR,INST_LI7200_PLL,INST_LI7200_SYNC,&
                &INST_LI7500_CHOPPER,INST_LI7500_DETECTOR,INST_LI7500_PLL,INST_LI7500_SYNC,&
                &INST_LI7700_NOT_READY,INST_LI7700_NO_SIGNAL,INST_LI7700_RE_UNLOCKED,&
                &INST_LI7700_BAD_TEMP,INST_LI7700_LASER_T_UNREG,INST_LI7700_CLOCK_T_UNREG,&
                &INST_LI7700_MOTOR_SPINNING,INST_LI7700_PUMP_ON,INST_LI7700_TOP_HEATER_ON,&
                &INST_LI7700_BOTTOM_HEATER_ON,INST_LI7700_CALIBRATING,INST_LI7700_MOTOR_FAILURE,&
                &INST_LI7700_BAD_AUX_TC1,INST_LI7700_BAD_AUX_TC2,INST_LI7700_BAD_AUX_TC3,INST_LI7700_BOX_CONNECTED,&
                &INST_LI7200_AGC_OR_RSSI,INST_LI7500_AGC_OR_RSSI,INST_LI7700_RSSI,&
                &WBOOST_APPLIED,AOA_METHOD,&
                &AXES_ROTATION_METHOD,ROT_YAW,ROT_PITCH,ROT_ROLL,&
                &DETRENDING_METHOD,DENTRENDING_TIME_CONSTANT,&
                &TIMELAG_DETECTION_METHOD,WPL_APPLIED,BURBA_METHOD,&
                &SPECTRAL_CORRECTION_METHOD,FOOTPRINT_MODEL,&
                &LOGGER_SWVER_MAJOR,LOGGER_SWVER_MINOR,LOGGER_SWVER_REVISION,&
                &BADM_LOCATION_LAT,BADM_LOCATION_LONG,BADM_LOCATION_ELEV,BADM_HEIGHTC,&
                &DISPLACEMENT_HEIGHT,ROUGHNESS_LENGTH,&
                &FILE_TIME_DURATION,BADM_INST_SAMPLING_INT,BADM_INST_AVERAGING_INT,&
                &MANUFACTURER_SA,BADM_INST_MODEL_SA,BADM_INST_HEIGHT_SA,&
                &BADM_INST_SA_WIND_FORMAT,BADM_INST_SA_GILL_ALIGN,BADM_SA_OFFSET_NORTH,&
                &HPATH_SA,VPATH_SA,RESPONSE_TIME_SA,&
                &MANUFACTURER_GA_CO2,BADM_INST_MODEL_GA_CO2,&
                &BADM_INSTPAIR_NORTHWARD_SEP_GA_CO2,BADM_INSTPAIR_EASTWARD_SEP_GA_CO2,BADM_INSTPAIR_HEIGHT_SEP_GA_CO2,&
                &BADM_INST_GA_CP_TUBE_LENGTH_GA_CO2,BADM_INST_GA_CP_TUBE_IN_DIAM_GA_CO2,BADM_INST_GA_CP_TUBE_FLOW_RATE_GA_CO2,&
                &HPATH_GA_CO2,VPATH_GA_CO2,RESPONSE_TIME_GA_CO2,&
                &MANUFACTURER_GA_H2O,BADM_INST_MODEL_GA_H2O,&
                &BADM_INSTPAIR_NORTHWARD_SEP_GA_H2O,BADM_INSTPAIR_EASTWARD_SEP_GA_H2O,BADM_INSTPAIR_HEIGHT_SEP_GA_H2O,&
                &BADM_INST_GA_CP_TUBE_LENGTH_GA_H2O,BADM_INST_GA_CP_TUBE_IN_DIAM_GA_H2O,BADM_INST_GA_CP_TUBE_FLOW_RATE_GA_H2O,&
                &KRYPTON_HYDRO_KH2O_GA_H2O,KRYPTON_HYDRO_KO2_GA_H2O,&
                &HPATH_GA_H2O,VPATH_GA_H2O,RESPONSE_TIME_GA_H2O,&
                &MANUFACTURER_GA_CH4,BADM_INST_MODEL_GA_CH4,&
                &BADM_INSTPAIR_NORTHWARD_SEP_GA_CH4,BADM_INSTPAIR_EASTWARD_SEP_GA_CH4,BADM_INSTPAIR_HEIGHT_SEP_GA_CH4,&
                &BADM_INST_GA_CP_TUBE_LENGTH_GA_CH4,BADM_INST_GA_CP_TUBE_IN_DIAM_GA_CH4,BADM_INST_GA_CP_TUBE_FLOW_RATE_GA_CH4,&
                &HPATH_GA_CH4,VPATH_GA_CH4,RESPONSE_TIME_GA_CH4,&
                &MANUFACTURER_GA_GS4,BADM_INST_MODEL_GA_GS4,&
                &BADM_INSTPAIR_NORTHWARD_SEP_GA_GS4,BADM_INSTPAIR_EASTWARD_SEP_GA_GS4,BADM_INSTPAIR_HEIGHT_SEP_GA_GS4,&
                &BADM_INST_GA_CP_TUBE_LENGTH_GA_GS4,BADM_INST_GA_CP_TUBE_IN_DIAM_GA_GS4,BADM_INST_GA_CP_TUBE_FLOW_RATE_GA_GS4,&
                &HPATH_GA_GS4,VPATH_GA_GS4,RESPONSE_TIME_GA_GS4,'

            !> If need to reitroduce details of VM, paste this after line:
            !>  "&CO2_ABSLIM_NREX,H2O_ABSLIM_NREX,CH4_ABSLIM_NREX,GS4_ABSLIM_NREX,&"

            !   &U_SPIKE_VM97_NR,V_SPIKE_VM97_NR,W_SPIKE_VM97_NR,T_SONIC_SPIKE_VM97_NR,&
            !   &CO2_SPIKE_VM97_NR,H2O_SPIKE_VM97_NR,CH4_SPIKE_VM97_NR,GS4_SPIKE_VM97_NR,&
            !   &U_AMPRES_VM97,V_AMPRES_VM97,W_AMPRES_VM97,T_SONIC_AMPRES_VM97,&
            !   &CO2_AMPRES_VM97,H2O_AMPRES_VM97,CH4_AMPRES_VM97,GS4_AMPRES_VM97,&
            !   &U_DRPOUT_C_VM97_NR,V_DRPOUT_C_VM97_NR,W_DRPOUT_C_VM97_NR,T_SONIC_DRPOUT_C_VM97_NR,&
            !   &CO2_DRPOUT_C_VM97_NR,H2O_DRPOUT_C_VM97_NR,CH4_DRPOUT_C_VM97_NR,GS4_DRPOUT_C_VM97_NR,&
            !   &U_DRPOUT_X_VM97_NR,V_DRPOUT_X_VM97_NR,W_DRPOUT_X_VM97_NR,T_SONIC_DRPOUT_X_VM97_NR,&
            !   &CO2_DRPOUT_X_VM97_NR,H2O_DRPOUT_X_VM97_NR,CH4_DRPOUT_X_VM97_NR,GS4_DRPOUT_X_VM97_NR,&
            !   &U_SKW_VM97_NR,V_SKW_VM97_NR,W_SKW_VM97_NR,T_SONIC_SKW_VM97_NR,&
            !   &CO2_SKW_VM97_NR,H2O_SKW_VM97_NR,CH4_SKW_VM97_NR,GS4_SKW_VM97_NR,&
            !   &U_KUR_VM97_NR,V_KUR_VM97_NR,W_KUR_VM97_NR,T_SONIC_KUR_VM97_NR,&
            !   &CO2_KUR_VM97_NR,H2O_KUR_VM97_NR,CH4_KUR_VM97_NR,GS4_KUR_VM97_NR,&
            !   &AOA_VM97_NR,WS_SS_ALONG_VM97,WS_SS_CROSS_VM97,WS_SS_VM97,&

            !> If need to reitroduce VM flags for last 3 tests, paste this after:
            !> "&VM97_DISCON_HFLAG,VM97_DISCON_SFLAG,&"

            !   &VM97_TIMELAG_HFLAG,VM97_TIMELAG_SFLAG,VM97_AOA_HFLAG,VM97_NSW_HFLAG,&


    !> Add custom variables
    call AddDatum(dataline, 'NUM_CUSTOM_VARS', separator)
    if (NumUserVar > 0) then
        do i = 1, NumUserVar
            call uppercase(usg(i)) 
            dataline = dataline(1:len_trim(dataline)) &
                // 'CUSTOM_' // usg(i)(1:len_trim(usg(i))) // 'MEAN' // ','
        end do
    end if

    !> Add biomet variables
    call AddDatum(dataline, 'NUM_BIOMET_VARS', separator)

    if (nbVars > 0) then
        do i = 1, nbVars
            if (EddyProProj%fluxnet_standardize_biomet) then
                call AddDatum(dataline, trim(bVars(i)%fluxnet_label), separator)
            else
                call AddDatum(dataline, trim(bVars(i)%label), separator)
            end if
        end do
    end if

    call uppercase(e2sg(gas4)) 
    dataline = replace2(dataline, 'GS4', e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1)) 

    write(uflxnt, '(a)') dataline(1:len_trim(dataline) - 1)

end subroutine