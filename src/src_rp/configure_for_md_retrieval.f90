!***************************************************************************
! configure_for_md_retrieval.f90
! ------------------------------
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
! \brief       Shut down everything except production of metadata file
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ConfigureForMdRetrieval()
    use m_rp_global_var
    implicit none


    !> raw data processing methods
    Meth%tlag = 'none'
    Meth%det = 'ba'
    Meth%rot = 'none'
    Meth%qcflag = 'none'
    Meth%foot = 'none'
    RUsetup%meth = 'none'
    RPsetup%bu_corr = 'none'
    RPsetup%calib_aoa = 'none'
    RPsetup%bu_multi = .false.
    RPsetup%calib_cw = .false.
    RPsetup%filter_by_raw_flags = .false.
    EddyProProj%use_extmd_file = .false.
    EddyProProj%biomet_data = 'none'
    EddyProProj%wpl = .false.
    EddyProProj%hf_meth = 'none'

    !> Raw statistical tests
    Test%sr = .false.
    Test%ar = .false.
    Test%do = .false.
    Test%al = .false.
    Test%sk = .false.
    Test%ds = .false.
    Test%tl = .false.
    Test%aa = .false.
    Test%ns = .false.
    RPsetup%offset(u) = 0d0
    RPsetup%offset(v) = 0d0
    RPsetup%offset(w) = 0d0

    !> Output files and other settings
    EddyProProj%out_md         = .true.
    EddyProProj%out_fluxnet    = .false.
    EddyProProj%out_full       = .false.
    EddyProProj%out_avrg_cosp  = .false.
    EddyProProj%out_biomet     = .false.
    RPsetup%out_st             = .false.
    RPsetup%filter_sr          = .false.
    RPsetup%filter_al          = .false.
    RPsetup%out_qc_details     = .false.
    RPsetup%out_raw            = .false.
    RPsetup%out_bin_sp         = .false.
    RPsetup%out_bin_og         = .false.
    RPsetup%out_full_sp        = .false.
    RPsetup%out_full_cosp      = .false.
    EddyProProj%fcc_follows    = .false.
    EddyProProj%make_dataset   = .true.
end subroutine ConfigureForMdRetrieval
