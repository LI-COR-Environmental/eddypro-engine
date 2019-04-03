!***************************************************************************
! eliminate_corrupted_variables.f90
! ---------------------------------
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
! \brief       If a variable as more than 30% values set to error, set it
!              as not present.
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine EliminateCorruptedVariables(LocSet, nrow, ncol, skip_period, logout)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(in) :: LocSet(nrow, ncol)
    logical, intent(in)  :: logout
    logical, intent(out) :: skip_period
    !> local variables
    integer :: i
    logical :: mask(nrow)


    if (logout) write(*,'(a)', advance = 'no') '  Verifying time series integrity..'

    do i = 1, ncol
        mask(:) = Locset(:, i) == error
        if (count(mask) > MaxPeriodNumRecords * RPsetup%max_lack/1d2) E2Col(i) = NullCol
    end do

    skip_period = .false.
    if ((.not. E2Col(u)%present) .or. &
        (.not. E2Col(v)%present) .or. &
        (.not. E2Col(w)%present)) skip_period = .true.

    if (logout) write(*,'(a)') ' Done.'
end subroutine EliminateCorruptedVariables
