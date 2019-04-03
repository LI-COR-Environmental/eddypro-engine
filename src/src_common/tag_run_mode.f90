!***************************************************************************
! tag_run_mode.f90
! ----------------
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
! \brief       Add token to Timestamp_FilePadding to tag run mode
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine TagRunMode()
    use m_common_global_var
    implicit none


    select case(EddyProProj%run_mode)
        case('express')
            Timestamp_FilePadding = trim(Timestamp_FilePadding) // '_exp'
        case('advanced')
            Timestamp_FilePadding = trim(Timestamp_FilePadding) // '_adv'
        case('md_retrieval')
            Timestamp_FilePadding = trim(Timestamp_FilePadding) // '_mdr'
        case default
            Timestamp_FilePadding = trim(Timestamp_FilePadding) // '_nul'
    end select
end subroutine TagRunMode
