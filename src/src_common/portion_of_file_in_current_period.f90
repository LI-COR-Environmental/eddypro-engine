!***************************************************************************
! portion_of_file_in_current_period.f90
! -------------------------------------
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
! \brief       Calculate first and last record of current file that are to be
!              associated to current period
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine PortionOfFileInCurrentPeriod(InitialTimestamp, FinalTimestamp, FileInitialTimestamp, FirstRecord, LastRecord)
    use m_common_global_var
    implicit none
    !> In/out variables
    type(DateType), intent(in) :: InitialTimestamp
    type(DateType), intent(in) :: FinalTimestamp
    type(DateType), intent(in) :: FileInitialTimestamp
    integer, intent(out) :: FirstRecord
    integer, intent(out) :: LastRecord


    !> Note:  LastRecord = -1 means: read to the end of the file.

    !> Initialization
    FirstRecord = nint(error)
    LastRecord = nint(error)

    !> Ideal case: File matches flux period
    if (FileInitialTimestamp == InitialTimestamp .and. &
        FileInitialTimestamp + DatafileDateStep == FinalTimestamp) then
        FirstRecord = 1
        LastRecord = -1
        return
    end if

    !> If beginning of file is same or after InitialTimestamp..
    if (FileInitialTimestamp >= InitialTimestamp) then
        !> ..FirstRecord is 1
        FirstRecord = 1
    else
        !> Function TimeLag calculates the time (in fractions of a day) between the beginning of the file
        !> and the BEGINNING of the current period, then converted to number of records
        FirstRecord = nint(TimeLag(InitialTimestamp, FileInitialTimestamp) * 1440d0 * 60d0 * Metadata%ac_freq) + 1
    end if

    !> If end of file is before (or equal to) FinalTimestamp, LastRecord is -1
    if (FileInitialTimestamp + DatafileDateStep <= FinalTimestamp) then
        LastRecord = -1
    else
        !> Function TimeLag calculates the time (in fractions of a day) between the beginning of the file
        !> and the END of the current period, then converted to number of records
        LastRecord = nint(TimeLag(FinalTimestamp, FileInitialTimestamp) * 1440d0 * 60d0 * Metadata%ac_freq)
    end if

end subroutine PortionOfFileInCurrentPeriod
