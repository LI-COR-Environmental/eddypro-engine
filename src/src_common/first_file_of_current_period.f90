!*******************************************************************************
! first_file_of_current_period.f90
! --------------------------------
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
!*******************************************************************************
!
! \brief       Find name of file (in FileList), closest to current period
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!*******************************************************************************
 subroutine FirstFileOfCurrentPeriod(InitialTimestamp, FinalTimestamp, &
    FileList, NumRawFiles, LatestRawFile, NextRawFile, skip_period)
    use m_common_global_var
    implicit none
    !> In/out variables
    integer, intent(in) :: LatestRawFile
    integer, intent(in) :: NumRawFiles
    type(DateType), intent(in) :: InitialTimestamp
    type(DateType), intent(in) :: FinalTimestamp
    type(FileListType), intent(in) :: FileList(NumRawFiles)
    integer, intent(out) :: NextRawFile
    logical, intent(out) :: skip_period
    !> Local variables
    integer :: i
    logical, external :: FileIsRelevantToCurrentPeriod

    skip_period = .false.


    !> Looks for the first file containing data relevant
    !> to the current period, if any
    NextRawFile = LatestRawFile
    do i = LatestRawFile, NumRawFiles
        if(FileIsRelevantToCurrentPeriod(FileList(i)%name, &
            InitialTimestamp, FinalTimestamp)) then
            NextRawFile = i
            return
        end if

        !> Check if the cycle is moving away from InitialTimestamp.
        !> In that case exit loop. This should keep the loop always very short
        if (FileList(i)%Timestamp >= InitialTimestamp + DateStep) then
            skip_period = .true.
            return
        end if
    end do
end subroutine FirstFileOfCurrentPeriod

!*******************************************************************************
!
! \brief       Detect whether a given raw file is relevant
!              to the given time period
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!*******************************************************************************
logical function FileIsRelevantToCurrentPeriod(Filename, &
    InitialTimestamp, FinalTimestamp)
    use m_common_global_var
    implicit none
    !> In/out variables
    character(*), intent(in) :: Filename
    type(DateType), intent(in) :: InitialTimestamp
    type(DateType), intent(in) :: FinalTimestamp
    !> Local variables
    character(10) :: date
    character(5)  :: time
    type(DateType) :: Timestamp

    !> Retrieve timestamp of beginning of current file
    FileIsRelevantToCurrentPeriod = .false.
    call FilenameToDateTime(Filename, EddyProProj%fname_template, &
        EddyProLog%iso_format, date, time)
    call DateTimeToDateType(date, time, Timestamp)
    if (EddyProLog%tstamp_end) Timestamp = Timestamp - DatafileDateStep
    !> 2. Now check if file contains data relevant to the current period
    if (Timestamp == InitialTimestamp .or. &
        (Timestamp < InitialTimestamp .and. &
        Timestamp + DatafileDateStep > InitialTimestamp) .or. &
        (Timestamp > InitialTimestamp .and. Timestamp < FinalTimestamp)) &
        FileIsRelevantToCurrentPeriod = .true.

end function FileIsRelevantToCurrentPeriod
