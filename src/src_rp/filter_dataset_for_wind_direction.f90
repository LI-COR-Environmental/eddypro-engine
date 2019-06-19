!***************************************************************************
! filter_data_for_wind_direction.f90
! ----------------------------------
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
! \brief       Filters dataset when wind comes from user-selected directions
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo        
!***************************************************************************
subroutine FilterDatasetForWindDirection(Set, nrow, ncol)
    use m_rp_global_var
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(inout) :: Set(nrow, ncol)
    !> lcoal variables
    integer :: i
    integer :: sec
    real(kind = dbl) :: WD


    write(*,'(a)', advance='no') '  Filtering raw data for wind direction filter..'
    Essentials%m_wdf = Essentials%m_wdf + 1
    do i = 1, nrow
        if (any(Set(i, u:w) == error)) cycle

        !> Instantaneous wind direction
        call SingleWindDirection(Set(i, u:w), &
            E2Col(u)%instr%north_offset + magnetic_declination, WD)

        !> Set record to error if wind is coming from an excluded sector
        sec_loop: do sec = 1, RPSetup%wdf_num_secs
            if (WD > RPSetup%wdf_start(sec) .and. WD < RPSetup%wdf_end(sec)) then
                Set (i, :) = error
                Essentials%m_wdf = Essentials%m_wdf + 1
                exit sec_loop
            end if 
        end do sec_loop
    end do
    write(*, '(a)') ' Done.'
    write(*, '(a, i6)') '   Number of records eliminated for wind direction filter: ',  Essentials%m_wdf

end subroutine FilterDatasetForWindDirection