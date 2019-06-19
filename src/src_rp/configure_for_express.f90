!***************************************************************************
! configure_for_express.f90
! -------------------------
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
! \brief       Set all entries to predefined values, valid for Express
!              This bypasses all user settings in terms of processing choices
! \author      Gerardo Fratini
! \note        Angle of attack correction is set later in main,
!              because it requires master_sonic
!              which is not know at this stage, when running in embedded mode
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ConfigureForExpress
    use m_rp_global_var
    implicit none


    !> raw data processing methods
    Meth%tlag = 'maxcov&default'
    Meth%det = 'ba'
    Meth%rot = 'double_rotation'
    Meth%qcflag = 'mauder_foken_04'
    Meth%foot = 'kljun_04'
    EddyProProj%hf_meth = 'moncrieff_97'
    EddyProProj%wpl  = .true.
    RUsetup%meth = 'none'
    RPsetup%bu_corr = 'none'
    RPsetup%bu_multi = .false.
    RPsetup%filter_sr = .true.
    RPsetup%filter_al = .true.
    RPsetup%calib_cw = .false.
    RPSetup%despike_vickers97 = .true.
    RPsetup%calib_aoa = 'automatic'

    !> Raw statistical tests
    Test%sr = .true.
    Test%ar = .true.
    Test%do = .true.
    Test%al = .true.
    Test%sk = .true.
    Test%ds = .false.
    Test%tl = .false.
    Test%aa = .false.
    Test%ns = .false.

    !> test parameters
    sr%num_spk = 3
    sr%lim_u   = 3.5d0
    sr%hf_lim  = 1d0
    sr%lim_w   = 5.0d0
    sr%lim_co2 = 3.5d0
    sr%lim_h2o = 3.5d0
    sr%lim_ch4 = 8d0
    sr%lim_gas4 = 8d0
    ar%lim     = 7
    ar%bins    = 100
    ar%hf_lim  = 70
    do%extlim_dw = 10d0
    do%hf1_lim = 10d0
    do%hf2_lim = 6d0
    al%u_max   = 30d0
    al%w_max   = 5d0
    al%t_min   = -40d0
    al%t_max   = 50d0
    al%co2_min = 200d0
    al%co2_max = 900d0
    al%h2o_min = 0d0
    al%h2o_max = 40d0
    al%ch4_min = 1.7d0 * 0.1d0  !< 1.7 ppm is minimum in unpolluted troposphere, 0.1 is safety factor
    al%ch4_max = 1000d0    !< to be better assessed
    al%gas4_min = 0.32d0 * 0.1d0  !< 0.32 ppm is minimum in unpolluted troposphere, 0.1 is safety factor
    al%gas4_max = 1000d0    !< to be better assessed
    sk%hf_skmin = -2d0
    sk%hf_skmax = 2d0
    sk%sf_skmin = -1d0
    sk%sf_skmax = 1d0
    sk%hf_kumin = 1d0
    sk%hf_kumax = 8d0
    sk%sf_kumin = 2d0
    sk%sf_kumax = 5d0
    BurbaPar%l(daytime, bot, 1)     =  0.944d0
    BurbaPar%l(daytime, bot, 2)     =  2.57d0
    BurbaPar%l(daytime, top, 1)     =  1.005d0
    BurbaPar%l(daytime, top, 2)     =  0.24d0
    BurbaPar%l(daytime, spar, 1)    =  1.01d0
    BurbaPar%l(daytime, spar, 2)    =  0.36d0
    BurbaPar%l(nighttime, bot, 1)   =  0.883d0
    BurbaPar%l(nighttime, bot, 2)   =  2.17d0
    BurbaPar%l(nighttime, top, 1)   =  1.008d0
    BurbaPar%l(nighttime, top, 2)   =  -0.41d0
    BurbaPar%l(nighttime, spar, 1)  =  1.01d0
    BurbaPar%l(nighttime, spar, 2)  =  -0.17d0

    RPsetup%offset(u) = 0d0
    RPsetup%offset(v) = 0d0
    RPsetup%offset(w) = 0d0

    !> Output files and other settings
    ! EddyProProj%out_fluxnet  = .false.
    EddyProProj%out_full     = .true.
    EddyProProj%out_md       = .true.
    RPsetup%out_st        = .true.
    RPsetup%out_qc_details = .false.
    RPsetup%out_raw        = .false.
    RPsetup%out_bin_sp     = .false.
    RPsetup%out_bin_og     = .false.
    RPsetup%out_full_sp    = .false.
    RPsetup%out_full_cosp  = .false.
    EddyProProj%out_avrg_cosp = .false.
    EddyProProj%out_avrg_spec = .false.
    EddyProProj%fcc_follows  = .false.
    EddyProProj%make_dataset = .true.

    if (EddyProProj%biomet_data /= 'none') then
        EddyProProj%out_biomet = .true.
    else
        EddyProProj%out_biomet = .false.
    end if

end subroutine ConfigureForExpress
