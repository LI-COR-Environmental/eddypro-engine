!*******************************************************************************
! first_file_of_current_period.f90
! --------------------------------
! Copyright (C) 2007-2011, Eco2s team, Gerardo Fratini
! Copyright (C) 2011-2014, LI-COR Biosciences
!
! This file is part of EddyPro (TM).
!
! EddyPro (TM) is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! EddyPro (TM) is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with EddyPro (TM).  If not, see <http://www.gnu.org/licenses/>.
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
