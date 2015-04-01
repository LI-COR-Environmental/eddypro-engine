!***************************************************************************
! m_fx_global_var.f90
! -------------------
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
! \brief       Module for global variables in eddypro_fcc
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
module m_fx_global_var
    use m_common_global_var
    implicit none
    save

    !> global variables
    integer :: g4l
    integer :: nstat
    integer :: nfull
    integer, parameter :: ndkf = 60     !< TO BE INCREASED????!!!!!!!
    integer, parameter :: unstable = 1
    integer, parameter :: stable   = 2

    real(kind = dbl) :: float_doy
    real(kind = dbl) :: dkf(ndkf + 1)

    character(11), parameter :: fcc_app = 'EddyPro-FCC'
    character(32) :: g4lab
    character(256) :: UserVarHeader
    character(25), parameter :: SubDirSpecAn = 'eddypro_spectral_analysis'

    logical :: MeanBinSpecAvailable(MaxGasClasses, GHGNumVar)
    logical :: MeanBinCospAvailable(MaxGasClasses, GHGNumVar)
    logical :: MeanStabSpecAvailable(MaxGasClasses, GHGNumVar)
    logical :: MeanStabCospAvailable(MaxGasClasses, GHGNumVar)
    logical :: fcc_var_present(GHGNumVar)

    type(FCCsetupType) :: FCCsetup
    type(FileListType), allocatable :: FullFileList(:)
    type(FileListType), allocatable :: BinnedFileList(:)

    type(MeanSpectraType), allocatable :: MeanBinSpec(:, :)
    type(MeanSpectraType), allocatable :: dMeanBinSpec(:, :)
    type(MeanSpectraType), allocatable :: MeanBinCosp(:, :)

    type(MeanSpectraType) :: MeanStabilityCosp(ndkf, 2)
    type(MassParType) :: MassPar(GHGNumVar, 2)
    type(FluxType) :: Flux1
    type(FluxType) :: Flux2
    type(FluxType) :: Flux3
    type(FCCMetadataType) :: FCCMetadata

    !> tags of the setup ".ini" file for eccoce
    integer, parameter :: Nsn = 109
    integer, parameter :: Nsc = 25
    logical            :: SNTagFound(Nsn)
    logical            :: SCTagFound(Nsc)
    type (Numerical)   :: SNTags(Nsn)
    type (Text)        :: SCTags(Nsc)
    data SNTags(1)%Label   / 'sa_nbins'         / & !< not used
         SNTags(2)%Label   / 'sa_min_smpl'      / &
         SNTags(3)%Label   / 'sa_fmin_co2'      / &
         SNTags(4)%Label   / 'sa_fmax_co2'      / &
         SNTags(5)%Label   / 'sa_fmin_h2o'      / &
         SNTags(6)%Label   / 'sa_fmax_h2o'      / &
         SNTags(7)%Label   / 'sa_fmin_ch4'      / &
         SNTags(8)%Label   / 'sa_fmax_ch4'      / &
         SNTags(9)%Label   / 'sa_fmin_gas4'     / &
         SNTags(10)%Label  / 'sa_fmax_gas4'     / &
         SNTags(11)%Label  / 'sa_min_co2'       / &  !< Depracted from 5.3
         SNTags(12)%Label  / 'sa_min_ch4'       / &  !< Depracted from 5.3
         SNTags(13)%Label  / 'sa_min_gas4'      / &  !< Depracted from 5.3
         SNTags(14)%Label  / 'sa_min_le'        / &  !< Depracted from 5.3
         SNTags(15)%Label  / 'sa_min_h'         / &  !< Depracted from 5.3
         SNTags(16)%Label  / 'sa_hfn_co2_fmin'  / &
         SNTags(17)%Label  / 'sa_hfn_h2o_fmin'  / &
         SNTags(18)%Label  / 'sa_hfn_ch4_fmin'  / &
         SNTags(19)%Label  / 'sa_hfn_gas4_fmin' / &
         SNTags(20)%Label  / 'sa_co2_g1_start'  / &
         SNTags(21)%Label  / 'sa_co2_g1_stop'   / &
         SNTags(22)%Label  / 'sa_co2_g2_start'  / &
         SNTags(23)%Label  / 'sa_co2_g2_stop'   / &
         SNTags(24)%Label  / 'sa_co2_g3_start'  / &
         SNTags(25)%Label  / 'sa_co2_g3_stop'   / &
         SNTags(26)%Label  / 'sa_co2_g4_start'  / &
         SNTags(27)%Label  / 'sa_co2_g4_stop'   / &
         SNTags(28)%Label  / 'sa_co2_g5_start'  / &
         SNTags(29)%Label  / 'sa_co2_g5_stop'   / &
         SNTags(30)%Label  / 'sa_co2_g6_start'  / &
         SNTags(31)%Label  / 'sa_co2_g6_stop'   / &
         SNTags(32)%Label  / 'sa_co2_g7_start'  / &
         SNTags(33)%Label  / 'sa_co2_g7_stop'   / &
         SNTags(34)%Label  / 'sa_co2_g8_start'  / &
         SNTags(35)%Label  / 'sa_co2_g8_stop'   / &
         SNTags(36)%Label  / 'sa_co2_g9_start'  / &
         SNTags(37)%Label  / 'sa_co2_g9_stop'   / &
         SNTags(38)%Label  / 'sa_co2_g10_start' / &
         SNTags(39)%Label  / 'sa_co2_g10_stop'  / &
         SNTags(40)%Label  / 'sa_co2_g11_start' / &
         SNTags(41)%Label  / 'sa_co2_g11_stop'  / &
         SNTags(42)%Label  / 'sa_co2_g12_start' / &
         SNTags(43)%Label  / 'sa_co2_g12_stop'  / &
         SNTags(44)%Label  / 'sa_ch4_g1_start'  / &
         SNTags(45)%Label  / 'sa_ch4_g1_stop'   / &
         SNTags(46)%Label  / 'sa_ch4_g2_start'  / &
         SNTags(47)%Label  / 'sa_ch4_g2_stop'   / &
         SNTags(48)%Label  / 'sa_ch4_g3_start'  / &
         SNTags(49)%Label  / 'sa_ch4_g3_stop'   / &
         SNTags(50)%Label  / 'sa_ch4_g4_start'  / &
         SNTags(51)%Label  / 'sa_ch4_g4_stop'   / &
         SNTags(52)%Label  / 'sa_ch4_g5_start'  / &
         SNTags(53)%Label  / 'sa_ch4_g5_stop'   / &
         SNTags(54)%Label  / 'sa_ch4_g6_start'  / &
         SNTags(55)%Label  / 'sa_ch4_g6_stop'   / &
         SNTags(56)%Label  / 'sa_ch4_g7_start'  / &
         SNTags(57)%Label  / 'sa_ch4_g7_stop'   / &
         SNTags(58)%Label  / 'sa_ch4_g8_start'  / &
         SNTags(59)%Label  / 'sa_ch4_g8_stop'   / &
         SNTags(60)%Label  / 'sa_ch4_g9_start'  / &
         SNTags(61)%Label  / 'sa_ch4_g9_stop'   / &
         SNTags(62)%Label  / 'sa_ch4_g10_start' / &
         SNTags(63)%Label  / 'sa_ch4_g10_stop'  / &
         SNTags(64)%Label  / 'sa_ch4_g11_start' / &
         SNTags(65)%Label  / 'sa_ch4_g11_stop'  / &
         SNTags(66)%Label  / 'sa_ch4_g12_start' / &
         SNTags(67)%Label  / 'sa_ch4_g12_stop'  / &
         SNTags(68)%Label  / 'sa_gas4_g1_start' / &
         SNTags(69)%Label  / 'sa_gas4_g1_stop'  / &
         SNTags(70)%Label  / 'sa_gas4_g2_start' / &
         SNTags(71)%Label  / 'sa_gas4_g2_stop'  / &
         SNTags(72)%Label  / 'sa_gas4_g3_start' / &
         SNTags(73)%Label  / 'sa_gas4_g3_stop'  / &
         SNTags(74)%Label  / 'sa_gas4_g4_start' / &
         SNTags(75)%Label  / 'sa_gas4_g4_stop'  / &
         SNTags(76)%Label  / 'sa_gas4_g5_start' / &
         SNTags(77)%Label  / 'sa_gas4_g5_stop'  / &
         SNTags(78)%Label  / 'sa_gas4_g6_start' / &
         SNTags(79)%Label  / 'sa_gas4_g6_stop'  / &
         SNTags(80)%Label  / 'sa_gas4_g7_start' / &
         SNTags(81)%Label  / 'sa_gas4_g7_stop'  / &
         SNTags(82)%Label  / 'sa_gas4_g8_start' / &
         SNTags(83)%Label  / 'sa_gas4_g8_stop'  / &
         SNTags(84)%Label  / 'sa_gas4_g9_start' / &
         SNTags(85)%Label  / 'sa_gas4_g9_stop'  / &
         SNTags(86)%Label  / 'sa_gas4_g10_start' / &
         SNTags(87)%Label  / 'sa_gas4_g10_stop'  / &
         SNTags(88)%Label  / 'sa_gas4_g11_start' / &
         SNTags(89)%Label  / 'sa_gas4_g11_stop'  / &
         SNTags(90)%Label  / 'sa_gas4_g12_start' / &
         SNTags(91)%Label  / 'sa_gas4_g12_stop'  / &
         SNTags(92)%Label   / 'sa_min_un_ustar' / &
         SNTags(93)%Label   / 'sa_min_un_co2'   / &
         SNTags(94)%Label   / 'sa_min_un_ch4'   / &
         SNTags(95)%Label   / 'sa_min_un_gas4'  / &
         SNTags(96)%Label   / 'sa_min_un_le'    / &
         SNTags(97)%Label   / 'sa_min_un_h'     / &
         SNTags(98)%Label   / 'sa_min_st_ustar' / &
         SNTags(99)%Label   / 'sa_min_st_co2'   / &
         SNTags(100)%Label   / 'sa_min_st_ch4'   / &
         SNTags(101)%Label  / 'sa_min_st_gas4'  / &
         SNTags(102)%Label  / 'sa_min_st_le'    / &
         SNTags(103)%Label  / 'sa_min_st_h'     / &
         SNTags(104)%Label  / 'sa_max_ustar'    / &
         SNTags(105)%Label  / 'sa_max_co2'      / &
         SNTags(106)%Label  / 'sa_max_ch4'      / &
         SNTags(107)%Label  / 'sa_max_gas4'     / &
         SNTags(108)%Label  / 'sa_max_le'       / &
         SNTags(109)%Label  / 'sa_max_h'        /

    data SCTags(1)%Label  / 'sa_start_date'     / &
         SCTags(2)%Label  / 'sa_start_time'     / &       !> not used
         SCTags(3)%Label  / 'sa_end_date'       / &
         SCTags(4)%Label  / 'sa_end_time'       / &       !> not used
         SCTags(5)%Label  / 'start_sa_date'     / &       !> not used
         SCTags(6)%Label  / 'end_sa_date'       / &       !> not used
         SCTags(7)%Label  / 'ex_file'           / &
         SCTags(8)%Label  / 'sa_bin_spectra'    / &
         SCTags(9)%Label  / 'sa_full_spectra'   / &
         SCTags(10)%Label / 'out_path'          / &       !> not used
         SCTags(11)%Label / 'make_dataset'      / &
         SCTags(12)%Label / 'sa_hfn_elim'       / &
         SCTags(13)%Label / 'h_meth'            / &
         SCTags(14)%Label / 'hf_meth'           / &
         SCTags(15)%Label / 'lptf_model'        / &
         SCTags(16)%Label / 'cosp_model'        / &
         SCTags(17)%Label / 'add_sonic_lptf'    / &
         SCTags(18)%Label / 'lf_meth'           / &
         SCTags(19)%Label / 'sa_mode'           / &
         SCTags(20)%Label / 'sa_file'           / &
         SCTags(21)%Label / 'horst_lens'        / &
         SCTags(22)%Label / 'sa_subset'         / &
         SCTags(23)%Label / 'sa_use_vm_flags'   / &
         SCTags(24)%Label / 'sa_use_foken_low'  / &
         SCTags(25)%Label / 'sa_use_foken_mid'  /
end module m_fx_global_var
