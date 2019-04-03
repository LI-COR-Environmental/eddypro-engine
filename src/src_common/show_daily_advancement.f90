!***************************************************************************
! show_daily_advancement.f90
! --------------------------
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
! \brief	   Prints out date (up to the day) of currently processed period
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine DisplayProgress(progress_type, init_message, tstamp, adv)
    use m_common_global_var
    implicit none
    !> in/out variables
    character(*) :: init_message
    character(*) :: progress_type
    character(*) :: adv
    type (DateType), intent(in) :: tstamp
    !> local variables
    character(10) :: date
    character(5) :: time
    character(15) :: months(12)
    data months(1:12) / 'January', 'February', 'March', &
        'April', 'May', 'June', 'July', 'August', &
        'September', 'October', 'November', 'December' /


    call DateTypeToDateTime(tstamp, date, time)
    select case(trim(adjustl(progress_type)))
    case ('daily')
        write(*, '(a)', advance = adv) init_message // date(9:10) // ' ' // &
            trim(adjustl(months(tstamp%month))) // ' ' // date(1:4) // ' '

    case ('avrg_interval')
        write(*, '(a)', advance = adv) init_message // time(1:5)
    end select

end subroutine DisplayProgress
