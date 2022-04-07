!***************************************************************************
! longest_gap_duration.f90
! ------------------------
! Copyright (C) 2018-2019, LI-COR Biosciences, Inc.  All Rights Reserved.
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
! \brief       Compute length of longest gaps in raw time series. Note that
!              because EddyPro does not handle raw-level timestamps, the
!              only gaps that can be "measured" are those due to elimination 
!              of data points, replaced by error codes. As usual in EddyPro,
!              the assumption is made that all data points are continuous.
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine LongestGapDuration(Set, nrow, ncol)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    !> Local variables
    integer :: var
    integer, external :: LongestVariableGap

    do var = u, GHGNumVar
        if (E2Col(var)%present) then
            Essentials%LGD(var) = LongestVariableGap(Set(:, var), nrow) / Metadata%ac_freq
        else
            Essentials%LGD(var) = error
        end if
    end do
end subroutine LongestGapDuration


integer function LongestVariableGap(arr, nrow)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    real(kind = dbl), intent(in) :: arr(nrow)
    !> Local variables
    integer :: cnt, i

    LongestVariableGap = 0
    i = 1
    do while (i < nrow)
        cnt = 0
        do while (arr(i) == error .and. i < nrow)
            cnt = cnt + 1
            i = i + 1
        end do
        if (cnt > LongestVariableGap) LongestVariableGap = cnt
        i = i + 1
    end do
end function LongestVariableGap