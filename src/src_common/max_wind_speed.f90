!***************************************************************************
! max_wind_speed.f90
! ------------------
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
! \brief       Calculate maximal wind speed in a 3d wind array
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine MaxWindSpeed(Set, nrow, ncol, MaxSpeed)
    use m_common_global_var
    implicit none
    !> In/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl) , intent(in) :: Set(nrow, ncol)
    real(kind = dbl) , intent(out) :: MaxSpeed
    !> Local variables
    integer :: i
    real(kind = dbl)  :: CurrentSpeed

    MaxSpeed = 0d0
    do i = 1, nrow
        CurrentSpeed = dsqrt(Set(i, u)**2 + Set(i, v)**2 + Set(i, w)**2)
        if (Set(i, u) /= error .and. Set(i, v) /= error &
            .and. Set(i, w) /= error .and. &
         CurrentSpeed > MaxSpeed) MaxSpeed = CurrentSpeed
    end do

end subroutine MaxWindSpeed
