!***************************************************************************
! count_records_and_values.f90
! ----------------------------
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
! \brief       Count either:
!              Available records (any record with at least one non-error value)
!              Available values for a passed variable
!              Available values for a pair of passed variables
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
integer function CountRecordsAndValues(Set, nrow, ncol, var1, var2)
    use m_common_global_var
    implicit none
    !> In/out variables
    integer, intent(in) :: nrow, ncol
    integer, optional, intent(in) :: var1, var2
    real(kind = dbl), intent(in) :: Set(nrow, ncol)


    if (.not. present(var1) .and. .not. present(var2)) then
        !> No var passed, count whole records
        CountRecordsAndValues = count(any(Set(:, 1:ncol) /= error, dim = 2))
        return
    else if (present(var1) .and. .not. present(var2)) then
        CountRecordsAndValues = count(Set(:, var1) /= error)
        return
    else
        CountRecordsAndValues = count(Set(:, var1) /= error .and. Set(:, var2) /= error)
        return
    end if

end function CountRecordsAndValues