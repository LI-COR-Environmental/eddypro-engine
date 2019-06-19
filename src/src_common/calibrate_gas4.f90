!***************************************************************************
! calibrate_gas4.f90
! ------------------
! Copyright (C) 2007-2011, Eco2s team, Gerardo Fratini
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
! \brief       Calibrate 4th gas if a cal-ref column is available.
!              Note that so far the calibration procedure is fully customized
!              on the needs of a specific O3 analyzer
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine CalibrateGas4(Set, nrow, ncol)
    use m_common_global_var
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(inout) :: Set(nrow, ncol)
    !> local variables
    integer :: j

    do j = 1, NumUserVar
        if (UserCol(j)%var == 'cal-ref') then
            !> Converts mV in ppb and then to ppm (with "/ 1d3")
            Set(:, gas4) = Set(:, gas4) * UserStats%Mean(j) &
            / Stats%Mean(gas4) / 1d3
        end if
    end do
end subroutine CalibrateGas4
