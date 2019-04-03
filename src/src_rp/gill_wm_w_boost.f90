!***************************************************************************
! gill_wm_w_boost.f90
! -------------------
! Copyright (C) 2016-2019, LI-COR Biosciences, Inc.  All Rights Reserved.
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
! \brief       Fix firmware bug for WindMaster and WindMaster Pro with firmware
!              versions of the form 2329.x.y with x < 700.
! \author      Gerardo Fratini
! \notes
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ApplyGillWmWBoost(Set, nrow, ncol)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(inout) :: Set(nrow, ncol)
    include '../src_common/interfaces.inc'

    Set(:, w) = asymmetric_linear_transformation(Set(:, w), size(Set(:, w), 1), &
                                                1.166d0, 0d0, 1.289d0, 0d0)
end subroutine ApplyGillWmWBoost


!***************************************************************************
!
! \brief       "Apply" firmware bug for WindMaster and WindMaster Pro with firmware
!              versions other than [2329.x.y with x < 700]. Needed if one
!              wants to apply the AoA correction of Nakai and Shimoyama 2012.
! \author      Gerardo Fratini
! \notes
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine RemoveGillWmWBoost(Set, nrow, ncol)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(inout) :: Set(nrow, ncol)
    include '../src_common/interfaces.inc'

    Set(:, w) = asymmetric_linear_transformation(Set(:, w), size(Set(:, w), 1), &
                                                1d0/1.166d0, 0d0, 1d0/1.289d0, 0d0)
end subroutine RemoveGillWmWBoost
