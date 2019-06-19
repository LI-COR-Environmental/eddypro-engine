!***************************************************************************
! infer_aoa_method.f90
! --------------------
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
! \brief       Infer most appropriate AoA method based on , \n
!              sonic anemometer model and serial number
! \author      Gerardo Fratini
! \notes
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine InferAoaMethod()
    use m_rp_global_var
    implicit none
    !> in/out variables
    include '../src_common/interfaces_1.inc'

    !> AoA selection based only on sonic model
    select case(MasterSonic%model(1:len_trim(MasterSonic%model) - 2))
        case ('r3_50','r3_100', 'r2')
            RPsetup%calib_aoa = 'nakai_06'
        case ('wm','wmpro')
            RPsetup%calib_aoa = 'nakai_12'
         case default
            RPsetup%calib_aoa = 'none'
    end select
end subroutine InferAoaMethod
