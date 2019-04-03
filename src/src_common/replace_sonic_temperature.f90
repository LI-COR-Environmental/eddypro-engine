!***************************************************************************
! replace_sonic_temperature.f90
! -----------------------------
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
! \brief       If sonic (or fast) temperature is out-ranged, tries to replace it
!              with any other available fast temperature reading
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ReplaceSonicTemperature(Set, nrow, ncol, UserSet, unrow, uncol)
    use m_common_global_var
    implicit none
    !> In/out variables
    integer, intent(in) :: nrow, ncol
    integer, intent(in) :: unrow, uncol
    real(kind = dbl), intent(in) :: UserSet(unrow,uncol)
    real(kind = dbl), intent(inout) :: Set(nrow, ncol)
    !> Local variables
    integer :: i
    integer :: ord


    !> Look for a suitable temperature reading
    ord = nint(error)
    do i = 1, uncol
        if (UserCol(i)%var == 'ts' .and. &
            (UserStats%Mean(i) > 220d0 .and. UserStats%Mean(i) < 340d0)) ord = i
    end do

    do i = 1, uncol
        if (UserCol(i)%var == 'fast_t' .and. &
            (UserStats%Mean(i) > 220d0 .and. UserStats%Mean(i) < 340d0)) ord = i
    end do

    !> If a suitable temperature was found, replace the current one with the good one
    if (ord /= nint(error)) Set(1:nrow, ts) = UserSet(1:nrow, ord)
end subroutine ReplaceSonicTemperature
