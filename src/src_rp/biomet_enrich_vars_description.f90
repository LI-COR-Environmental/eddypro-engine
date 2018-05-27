!***************************************************************************
! biomet_enrich_vars_description.f90
! ----------------------------------
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
!***************************************************************************
!
! \brief       Complete description of biomet vars with all information that
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
!***************************************************************************
subroutine BiometEnrichVarsDescription()
    use m_rp_global_var
    implicit none
    !> In/out variables
    !> Local variables
    integer :: i
    character(32) :: base_name
    ! character(32) :: qPos
    character(32), external :: positionalQualifier


    do i = 1, nbVars
        ! call BiometInterpretPositionalQualifier(bVars(i))
        base_name = trim(bVars(i)%base_name)
        call uppercase(base_name)
        select case(base_name)
            case('TA', 'T_A', 'T_AIR', 'TAIR')
                bVars(i)%fluxnet_base_name = 'TA'
                bVars(i)%nature = 'TEMPERATURE'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'K'
                bVars(i)%pretty_unit_out = '[K]'
                bVars(i)%fluxnet_unit_out = '[deg C]'
            case('TBC', 'T_BC', 'T_BELOW_CANOPY')
                bVars(i)%fluxnet_base_name = 'T_BC'
                bVars(i)%nature = 'TEMPERATURE'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'K'
                bVars(i)%pretty_unit_out = '[K]'
                bVars(i)%fluxnet_unit_out = '[deg C]'
            case('TR')
                bVars(i)%fluxnet_base_name = 'T_R'
                bVars(i)%nature = 'TEMPERATURE'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'K'
                bVars(i)%pretty_unit_out = '[K]'
                bVars(i)%fluxnet_unit_out = '[deg C]'
            case('TS', 'T_S', 'T_SOIL', 'TSOIL')
                bVars(i)%fluxnet_base_name = 'TS'
                bVars(i)%nature = 'TEMPERATURE'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'K'
                bVars(i)%pretty_unit_out = '[K]'
                bVars(i)%fluxnet_unit_out = '[deg C]'
            case('RH')
                bVars(i)%fluxnet_base_name = 'RH'
                bVars(i)%nature = 'RELATIVE_HUMIDITY'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = '%'
                bVars(i)%pretty_unit_out = '[%]'
                bVars(i)%fluxnet_unit_out = '[%]'
            case('PA', 'P_A', 'PAIR', 'P_AIR')
                bVars(i)%fluxnet_base_name = 'PA'
                bVars(i)%nature = 'PRESSURE'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'PA'
                bVars(i)%pretty_unit_out = '[Pa]'
                bVars(i)%fluxnet_unit_out = '[kPa]'
            case('VPD')
                bVars(i)%fluxnet_base_name = 'VPD'
                bVars(i)%nature = 'PRESSURE'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'PA'
                bVars(i)%pretty_unit_out = '[Pa]'
                bVars(i)%fluxnet_unit_out = '[hPa]'

            case('CO2')
                bVars(i)%fluxnet_base_name = 'CO2'
                bVars(i)%nature = 'CONCENTRATION'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'UMOL+1MOL-1'
                bVars(i)%pretty_unit_out = '[' // char(181) // 'mol/mol]'
                bVars(i)%fluxnet_unit_out = '[umolCO2 mol-1]'
            case('H2O')
                bVars(i)%fluxnet_base_name = 'H2O'
                bVars(i)%nature = 'CONCENTRATION'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'MMOL+1MOL-1'
                bVars(i)%pretty_unit_out = '[mmol/mol]'
                bVars(i)%fluxnet_unit_out = '[mmolH2O mol-1]'
            case('CH4')
                bVars(i)%fluxnet_base_name = 'CH4'
                bVars(i)%nature = 'CONCENTRATION'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'NMOL+1MOL-1'
                bVars(i)%pretty_unit_out = '[nmol/mol]'
                bVars(i)%fluxnet_unit_out = '[nmolCH4 mol-1]'
            case('NO')
                bVars(i)%fluxnet_base_name = 'NO'
                bVars(i)%nature = 'CONCENTRATION'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'NMOL+1MOL-1'
                bVars(i)%pretty_unit_out = '[nmol/mol]'
                bVars(i)%fluxnet_unit_out = '[nmolNO mol-1]'
            case('N2O')
                bVars(i)%fluxnet_base_name = 'N2O'
                bVars(i)%nature = 'CONCENTRATION'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'NMOL+1MOL-1'
                bVars(i)%pretty_unit_out = '[nmol/mol]'
                bVars(i)%fluxnet_unit_out = '[nmolN2O mol-1]'
            case('NO2')
                bVars(i)%fluxnet_base_name = 'NO2'
                bVars(i)%nature = 'CONCENTRATION'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'NMOL+1MOL-1'
                bVars(i)%pretty_unit_out = '[nmol/mol]'
                bVars(i)%fluxnet_unit_out = '[nmolNO2 mol-1]'
            case('O3')
                bVars(i)%fluxnet_base_name = 'O3'
                bVars(i)%nature = 'CONCENTRATION'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'NMOL+1MOL-1'
                bVars(i)%pretty_unit_out = '[nmol/mol]'
                bVars(i)%fluxnet_unit_out = '[nmolO3 mol-1]'
            case('CO')
                bVars(i)%fluxnet_base_name = 'CO'
                bVars(i)%nature = 'CONCENTRATION'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'NMOL+1MOL-1'
                bVars(i)%pretty_unit_out = '[nmol/mol]'
                bVars(i)%fluxnet_unit_out = '[nmolCO mol-1]'
            case('SO2')
                bVars(i)%fluxnet_base_name = 'SO2'
                bVars(i)%nature = 'CONCENTRATION'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'NMOL+1MOL-1'
                bVars(i)%pretty_unit_out = '[nmol/mol]'
                bVars(i)%fluxnet_unit_out = '[nmolSO2 mol-1]'
            case('NH3')
                bVars(i)%fluxnet_base_name = 'NH3'
                bVars(i)%nature = 'CONCENTRATION'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'NMOL+1MOL-1'
                bVars(i)%pretty_unit_out = '[nmol/mol]'
                bVars(i)%fluxnet_unit_out = '[nmolNH3 mol-1]'
            case('NETRAD', 'RN', 'R_N', 'NET_RAD','NET_RADIATION','RNET','R_NET')
                bVars(i)%fluxnet_base_name = 'NETRAD'
                bVars(i)%nature = 'RADIATION'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'W+1M-2'
                bVars(i)%pretty_unit_out = '[W/m^2]'
                bVars(i)%fluxnet_unit_out = '[W m-2]'
            case('R_UVA','RUVA')
                bVars(i)%fluxnet_base_name = 'R_UVA'
                bVars(i)%nature = 'RADIATION'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'W+1M-2'
                bVars(i)%pretty_unit_out = '[W/m^2]'
                bVars(i)%fluxnet_unit_out = '[W m-2]'
            case('R_UVB', 'RUVB')
                bVars(i)%fluxnet_base_name = 'R_UVB'
                bVars(i)%nature = 'RADIATION'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'W+1M-2'
                bVars(i)%pretty_unit_out = '[W/m^2]'
                bVars(i)%fluxnet_unit_out = '[W m-2]'
            case('RG', 'R_G', 'RGLOBAL','R_GLOBAL', 'SWIN','SW_IN')
                bVars(i)%fluxnet_base_name = 'SW_IN'
                bVars(i)%nature = 'RADIATION'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'W+1M-2'
                bVars(i)%pretty_unit_out = '[W/m^2]'
                bVars(i)%fluxnet_unit_out = '[W m-2]'
            case('SWDIR','SW_DIR')
                bVars(i)%fluxnet_base_name = 'SW_DIR'
                bVars(i)%nature = 'RADIATION'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'W+1M-2'
                bVars(i)%pretty_unit_out = '[W/m^2]'
                bVars(i)%fluxnet_unit_out = '[W m-2]'
            case('SWDIF','SW_DIF','RD','R_D')
                bVars(i)%fluxnet_base_name = 'SW_DIF'
                bVars(i)%nature = 'RADIATION'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'W+1M-2'
                bVars(i)%pretty_unit_out = '[W/m^2]'
                bVars(i)%fluxnet_unit_out = '[W m-2]'
            case('SWOUT','SW_OUT','RR', 'R_R')
                bVars(i)%fluxnet_base_name = 'SW_OUT'
                bVars(i)%nature = 'RADIATION'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'W+1M-2'
                bVars(i)%pretty_unit_out = '[W/m^2]'
                bVars(i)%fluxnet_unit_out = '[W m-2]'
            case('SWBC','SWBC_IN','SW_BC','SW_BC_IN')
                bVars(i)%fluxnet_base_name = 'SW_BC_IN'
                bVars(i)%nature = 'RADIATION'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'W+1M-2'
                bVars(i)%pretty_unit_out = '[W/m^2]'
                bVars(i)%fluxnet_unit_out = '[W m-2]'
            case('SWBC_OUT','SW_BC_OUT')
                bVars(i)%fluxnet_base_name = 'SW_BC_OUT'
                bVars(i)%nature = 'RADIATION'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'W+1M-2'
                bVars(i)%pretty_unit_out = '[W/m^2]'
                bVars(i)%fluxnet_unit_out = '[W m-2]'
            case('LWIN','LW_IN')
                bVars(i)%fluxnet_base_name = 'LW_IN'
                bVars(i)%nature = 'RADIATION'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'W+1M-2'
                bVars(i)%pretty_unit_out = '[W/m^2]'
                bVars(i)%fluxnet_unit_out = '[W m-2]'
            case('LWOUT','LW_OUT')
                bVars(i)%fluxnet_base_name = 'LW_OUT'
                bVars(i)%nature = 'RADIATION'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'W+1M-2'
                bVars(i)%pretty_unit_out = '[W/m^2]'
                bVars(i)%fluxnet_unit_out = '[W m-2]'
            case('LWBC_IN','LW_BC_IN')
                bVars(i)%fluxnet_base_name = 'LW_BC_IN'
                bVars(i)%nature = 'RADIATION'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'W+1M-2'
                bVars(i)%pretty_unit_out = '[W/m^2]'
                bVars(i)%fluxnet_unit_out = '[W m-2]'
            case('LWBCOUT','LW_BC_OUT')
                bVars(i)%fluxnet_base_name = 'LW_BC_OUT'
                bVars(i)%nature = 'RADIATION'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'W+1M-2'
                bVars(i)%pretty_unit_out = '[W/m^2]'
                bVars(i)%fluxnet_unit_out = '[W m-2]'

            case('SPEC_RED_IN','SPECRED_IN','SPECREDIN')
                bVars(i)%fluxnet_base_name = 'SPEC_RED_IN'
                bVars(i)%nature = 'PHOTON_FLUX'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'UMOL+1M-2S-1'
                bVars(i)%pretty_unit_out = '[' // char(181) // 'mol/m^2/s]'
                bVars(i)%fluxnet_unit_out = '[umolPhoton m-2s-1]'
            case('SPEC_RED_OUT','SPECRED_OUT','SPECREDOUT')
                bVars(i)%fluxnet_base_name = 'SPEC_RED_OUT'
                bVars(i)%nature = 'PHOTON_FLUX'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'UMOL+1M-2S-1'
                bVars(i)%pretty_unit_out = '[' // char(181) // 'mol/m^2/s]'
                bVars(i)%fluxnet_unit_out = '[umolPhoton m-2s-1]'
            case('SPEC_RED_REFL','SPECRED_REFL','SPECREDREFL')
                bVars(i)%fluxnet_base_name = 'SPEC_RED_REFL'
                bVars(i)%nature = 'INDEX'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = '#'
                bVars(i)%pretty_unit_out = '[#]'
                bVars(i)%fluxnet_unit_out = '[#]'
            case('SPEC_NIR_IN','SPECNIR_IN','SPECNIRIN')
                bVars(i)%fluxnet_base_name = 'SPEC_NIR_IN'
                bVars(i)%nature = 'PHOTON_FLUX'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'UMOL+1M-2S-1'
                bVars(i)%pretty_unit_out = '[' // char(181) // 'mol/m^2/s]'
                bVars(i)%fluxnet_unit_out = '[umolPhoton m-2s-1]'
            case('SPEC_NIR_OUT','SPECNIR_OUT','SPECNIROUT')
                bVars(i)%fluxnet_base_name = 'SPEC_NIR_OUT'
                bVars(i)%nature = 'PHOTON_FLUX'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'UMOL+1M-2S-1'
                bVars(i)%pretty_unit_out = '[' // char(181) // 'mol/m^2/s]'
                bVars(i)%fluxnet_unit_out = '[umolPhoton m-2s-1]'
            case('SPEC_NIR_REFL','SPECNIR_REFL','SPECNIRREFL')
                bVars(i)%fluxnet_base_name = 'SPEC_NIR_REFL'
                bVars(i)%nature = 'INDEX'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = '#'
                bVars(i)%pretty_unit_out = '[#]'
                bVars(i)%fluxnet_unit_out = '[#]'
            case('SPEC_PRI_TGT_IN','SPECPRITGT_IN','SPECPRITGTIN')
                bVars(i)%fluxnet_base_name = 'SPEC_PRI_TGT_IN'
                bVars(i)%nature = 'PHOTON_FLUX'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'UMOL+1M-2S-1'
                bVars(i)%pretty_unit_out = '[' // char(181) // 'mol/m^2/s]'
                bVars(i)%fluxnet_unit_out = '[umolPhoton m-2s-1]'
            case('SPEC_PRI_TGT_OUT','SPECPRITGT_OUT','SPECPRITGTOUT')
                bVars(i)%fluxnet_base_name = 'SPEC_PRI_TGT_OUT'
                bVars(i)%nature = 'PHOTON_FLUX'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'UMOL+1M-2S-1'
                bVars(i)%pretty_unit_out = '[' // char(181) // 'mol/m^2/s]'
                bVars(i)%fluxnet_unit_out = '[umolPhoton m-2s-1]'
            case('SPEC_PRI_TGT_REFL','SPECPRITGT_REFL','SPECPRITGTREFL')
                bVars(i)%fluxnet_base_name = 'SPEC_PRI_TGT_REFL'
                bVars(i)%nature = 'INDEX'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = '#'
                bVars(i)%pretty_unit_out = '[#]'
                bVars(i)%fluxnet_unit_out = '[#]'
            case('SPEC_PRI_REF_IN','SPECPRIREF_IN','SPECPRIREFIN')
                bVars(i)%fluxnet_base_name = 'SPEC_PRI_REF_IN'
                bVars(i)%nature = 'PHOTON_FLUX'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'UMOL+1M-2S-1'
                bVars(i)%pretty_unit_out = '[' // char(181) // 'mol/m^2/s]'
                bVars(i)%fluxnet_unit_out = '[umolPhoton m-2s-1]'
            case('SPEC_PRI_REF_OUT','SPECPRIREF_OUT','SPECPRIREFOUT')
                bVars(i)%fluxnet_base_name = 'SPEC_PRI_REF_OUT'
                bVars(i)%nature = 'PHOTON_FLUX'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'UMOL+1M-2S-1'
                bVars(i)%pretty_unit_out = '[' // char(181) // 'mol/m^2/s]'
                bVars(i)%fluxnet_unit_out = '[umolPhoton m-2s-1]'
            case('SPEC_PRI_REF_REFL','SPECPRIREF_REFL','SPECPRIREFREFL')
                bVars(i)%fluxnet_base_name = 'SPEC_PRI_REF_REFL'
                bVars(i)%nature = 'INDEX'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = '#'
                bVars(i)%pretty_unit_out = '[#]'
                bVars(i)%fluxnet_unit_out = '[#]'

            case('PPFD', 'PPFD_IN')
                bVars(i)%fluxnet_base_name = 'PPFD'
                bVars(i)%nature = 'PHOTON_FLUX'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'UMOL+1M-2S-1'
                bVars(i)%pretty_unit_out = '[' // char(181) // 'mol/m^2/s]'
                bVars(i)%fluxnet_unit_out = '[umolPhoton m-2s-1]'
            case('PPFD_DIR', 'PPFDDIR')
                bVars(i)%fluxnet_base_name = 'PPFD_DIR'
                bVars(i)%nature = 'PHOTON_FLUX'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'UMOL+1M-2S-1'
                bVars(i)%pretty_unit_out = '[' // char(181) // 'mol/m^2/s]'
                bVars(i)%fluxnet_unit_out = '[umolPhoton m-2s-1]'
            case('PPFD_DIF', 'PPFDDIF', 'PPFDD', 'PPFD_D')
                bVars(i)%fluxnet_base_name = 'PPFD_DIF'
                bVars(i)%nature = 'PHOTON_FLUX'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'UMOL+1M-2S-1'
                bVars(i)%pretty_unit_out = '[' // char(181) // 'mol/m^2/s]'
                bVars(i)%fluxnet_unit_out = '[umolPhoton m-2s-1]'
            case('PPFDR', 'PPFD_R', 'PPFD_OUT')
                bVars(i)%fluxnet_base_name = 'PPFD_OUT'
                bVars(i)%nature = 'PHOTON_FLUX'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'UMOL+1M-2S-1'
                bVars(i)%pretty_unit_out = '[' // char(181) // 'mol/m^2/s]'
                bVars(i)%fluxnet_unit_out = '[umolPhoton m-2s-1]'
            case('PPFDBC','PPFD_BC','PPFDBC_IN','PPFD_BC_IN')
                bVars(i)%fluxnet_base_name = 'PPFD_BC_IN'
                bVars(i)%nature = 'PHOTON_FLUX'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'UMOL+1M-2S-1'
                bVars(i)%pretty_unit_out = '[' // char(181) // 'mol/m^2/s]'
                bVars(i)%fluxnet_unit_out = '[umolPhoton m-2s-1]'
            case('PPFDBC_R','PPFD_BC_R','PPFDBC_OUT','PPFD_BC_OUT')
                bVars(i)%fluxnet_base_name = 'PPFD_BC_OUT'
                bVars(i)%nature = 'PHOTON_FLUX'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'UMOL+1M-2S-1'
                bVars(i)%pretty_unit_out = '[' // char(181) // 'mol/m^2/s]'
                bVars(i)%fluxnet_unit_out = '[umolPhoton m-2s-1]'
            case('APAR')
                bVars(i)%fluxnet_base_name = 'APAR'
                bVars(i)%nature = 'PHOTON_FLUX'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'UMOL+1M-2S-1'
                bVars(i)%pretty_unit_out = '[' // char(181) // 'mol/m^2/s]'
                bVars(i)%fluxnet_unit_out = '[umol m-2s-1]'
            case('FAPAR','F_APAR')
                bVars(i)%fluxnet_base_name = 'FAPAR'
                bVars(i)%nature = 'PHOTON_FLUX_FRACTION'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = '%'
                bVars(i)%pretty_unit_out = '[%]'
                bVars(i)%fluxnet_unit_out = '[%]'
            case('FIPAR','F_IPAR')
                bVars(i)%fluxnet_base_name = 'FIPAR'
                bVars(i)%nature = 'PHOTON_FLUX_FRACTION'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = '%'
                bVars(i)%pretty_unit_out = '[%]'
                bVars(i)%fluxnet_unit_out = '[%]'

            case('NDVI')
                bVars(i)%fluxnet_base_name = 'NDVI'
                bVars(i)%nature = 'INDEX'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = '#'
                bVars(i)%pretty_unit_out = '[#]'
                bVars(i)%fluxnet_unit_out = '[#]'
            case('PRI')
                bVars(i)%fluxnet_base_name = 'PRI'
                bVars(i)%nature = 'INDEX'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = '#'
                bVars(i)%pretty_unit_out = '[#]'
                bVars(i)%fluxnet_unit_out = '[#]'

            case('WS', 'WIND_SPEED')
                bVars(i)%fluxnet_base_name = 'WS'
                bVars(i)%nature = 'SPEED'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'M+1S-1'
                bVars(i)%pretty_unit_out = '[m/s]'
                bVars(i)%fluxnet_unit_out = '[m s-1]'
            case('MWS', 'MAX_WIND_SPEED')
                bVars(i)%fluxnet_base_name = 'MWS'
                bVars(i)%nature = 'SPEED'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'M+1S-1'
                bVars(i)%pretty_unit_out = '[m/s]'
                bVars(i)%fluxnet_unit_out = '[m s-1]'
            case('WD','WIND_DIRECTION')
                bVars(i)%fluxnet_base_name = 'WD'
                bVars(i)%nature = 'ANGULAR_DIRECTION'
                bVars(i)%accumul_type = 'ANGULAR_AVERAGING'
                bVars(i)%unit_out = 'DEGREES'
                bVars(i)%pretty_unit_out = '[Degrees_past_North]'
                bVars(i)%fluxnet_unit_out = '[Decimal degrees]'

            case('P', 'P_TOTAL', 'PTOTAL')
                bVars(i)%fluxnet_base_name = 'P'
                bVars(i)%nature = 'PRECIPITATION'
                bVars(i)%accumul_type = 'INTEGRATION'
                bVars(i)%unit_out = 'M'
                bVars(i)%pretty_unit_out = '[m]'
                bVars(i)%fluxnet_unit_out = '[mm]'
            case('PRAIN','P_RAIN')
                bVars(i)%fluxnet_base_name = 'P_RAIN'
                bVars(i)%nature = 'PRECIPITATION'
                bVars(i)%accumul_type = 'INTEGRATION'
                bVars(i)%unit_out = 'M'
                bVars(i)%pretty_unit_out = '[m]'
                bVars(i)%fluxnet_unit_out = '[mm]'
            case('PSNOW','P_SNOW')
                bVars(i)%fluxnet_base_name = 'P_SNOW'
                bVars(i)%nature = 'PRECIPITATION'
                bVars(i)%accumul_type = 'INTEGRATION'
                bVars(i)%unit_out = 'M'
                bVars(i)%pretty_unit_out = '[m]'
                bVars(i)%fluxnet_unit_out = '[mm]'

            case('RUNOFF','RUN_OFF')
                bVars(i)%fluxnet_base_name = 'RUNOFF'
                bVars(i)%nature = 'LENGTH'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'M'
                bVars(i)%pretty_unit_out = '[m]'
                bVars(i)%fluxnet_unit_out = '[mm]'
            case('SNOWD','SNOW_D','D_SNOW','DSNOW')
                bVars(i)%fluxnet_base_name = 'D_SNOW'
                bVars(i)%nature = 'LENGTH'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'M'
                bVars(i)%pretty_unit_out = '[m]'
                bVars(i)%fluxnet_unit_out = '[cm]'
            case('LAI')
                bVars(i)%fluxnet_base_name = 'LAI'
                bVars(i)%nature = 'INDEX'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'M+2M-2'
                bVars(i)%pretty_unit_out = '[m^2/m^2]'
                bVars(i)%fluxnet_unit_out = '[%]'
            case('ALB','ALBEDO')
                bVars(i)%fluxnet_base_name = 'ALBEDO'
                bVars(i)%nature = 'INDEX'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = '%'
                bVars(i)%pretty_unit_out = '[%]'
                bVars(i)%fluxnet_unit_out = '[%]'

            case('SHF', 'G')
                bVars(i)%fluxnet_base_name = 'G'
                bVars(i)%nature = 'HEAT_FLUX'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'W+1M-2'
                bVars(i)%pretty_unit_out = '[W/m^2]'
                bVars(i)%fluxnet_unit_out = '[W m-2]'
            case('SWC')
                bVars(i)%fluxnet_base_name = 'SWC'
                bVars(i)%nature = 'VOLUME_CONCENTRATION'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'M+3M-3'
                bVars(i)%pretty_unit_out = '[m^3/m^3]'
                bVars(i)%fluxnet_unit_out = '[%]'
            case('WATER_TABLE_DEPTH','WTD')
                bVars(i)%fluxnet_base_name = 'WATER_TABLE_DEPTH'
                bVars(i)%nature = 'LENGTH'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'M'
                bVars(i)%pretty_unit_out = '[m]'
                bVars(i)%fluxnet_unit_out = '[cm]'
            case('PBLH','PBL_HEIGHT','HPBL','PBL_H', 'H_PBL')
                bVars(i)%fluxnet_base_name = 'PBLH'
                bVars(i)%nature = 'LENGTH'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'M'
                bVars(i)%pretty_unit_out = '[m]'
                bVars(i)%fluxnet_unit_out = '[m]'

            !> BIOLOGICAL
            case('SAP_DT')
                bVars(i)%fluxnet_base_name = 'SAP_DT'
                bVars(i)%nature = 'DELTA_TEMPERATURE'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'K'
                bVars(i)%pretty_unit_out = '[K]'
                bVars(i)%fluxnet_unit_out = '[deg C]'
            case('SAP_FLOW','SAPFLOW')
                bVars(i)%fluxnet_base_name = 'SAP_FLOW'
                bVars(i)%nature = 'FLOW'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'MMOL+1M-2S-1'
                bVars(i)%pretty_unit_out = '[mmol/m^2/s]'
                bVars(i)%fluxnet_unit_out = '[mmolH2O m-2 s-1]'
            case('STEMFLOW','STEM_FLOW')
                bVars(i)%fluxnet_base_name = 'STEMFLOW'
                bVars(i)%nature = 'LENGTH'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'MM'
                bVars(i)%pretty_unit_out = '[mm]'
                bVars(i)%fluxnet_unit_out = '[mm]'
            case('THROUGHFALL','THROUGH_FALL')
                bVars(i)%fluxnet_base_name = 'THROUGHFALL'
                bVars(i)%nature = 'LENGTH'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'M'
                bVars(i)%pretty_unit_out = '[m]'
                bVars(i)%fluxnet_unit_out = '[mm]'
            case('TC', 'T_C', 'T_CANOPY', 'TCANOPY')
                bVars(i)%fluxnet_base_name = 'T_CANOPY'
                bVars(i)%nature = 'TEMPERATURE'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'K'
                bVars(i)%pretty_unit_out = '[K]'
                bVars(i)%fluxnet_unit_out = '[deg C]'
            case('TBOLE', 'T_BOLE')
                bVars(i)%fluxnet_base_name = 'T_BOLE'
                bVars(i)%nature = 'TEMPERATURE'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'K'
                bVars(i)%pretty_unit_out = '[K]'
                bVars(i)%fluxnet_unit_out = '[deg C]'
            case('DBH')
                bVars(i)%fluxnet_base_name = 'DBH'
                bVars(i)%nature = 'LENGTH'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'M'
                bVars(i)%pretty_unit_out = '[m]'
                bVars(i)%fluxnet_unit_out = '[cm]'
            case('LEAF_WET', 'LEAF_WETNESS')
                bVars(i)%fluxnet_base_name = 'LEAF_WET'
                bVars(i)%nature = 'INDEX'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = '%'
                bVars(i)%pretty_unit_out = '[%]'
                bVars(i)%fluxnet_unit_out = '[%]'

            case default
                bVars(i)%fluxnet_base_name = bVars(i)%base_name
                bVars(i)%nature = 'UNKNOWN'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = trim(bVars(i)%unit_in)
                bVars(i)%pretty_unit_out = '[' // trim(bVars(i)%unit_in) // ']'
                bVars(i)%fluxnet_unit_out = bVars(i)%pretty_unit_out
        end select
    end do

    !> Complete Fluxnet labels adding positional qualifier to Fluxnet base name
    do i = 1, nbVars
        ! qPos = positionalQualifier(bVars(i))
            bVars(i)%fluxnet_label = trim(bVars(i)%fluxnet_base_name) // trim(bVars(i)%pq_string)
    end do
end subroutine BiometEnrichVarsDescription

