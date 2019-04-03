!***************************************************************************
! ts_extract_subperiod_indexes_from_filelist.f90
! ----------------------------------------------
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
! \brief       From a list of chronologically ordered file names,
!              retrieve indexes of a selected subperiod
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine tsExtractSubperiodIndexesFromFilelist(FileList, nrow, &
        StartTimestamp, EndTimestamp, StartIndex, EndIndex)
    use m_common_global_var
    implicit none
    !> In/out variables
    integer, intent(in) :: nrow
    type (FileListType), intent(in) :: FileList(nrow)
    type (DateType), intent(in) :: StartTimestamp
    type (DateType), intent(in) :: EndTimestamp
    integer, intent(out) :: StartIndex
    integer, intent(out) :: EndIndex
    !> Local variables
    integer :: i


    StartIndex = nint(error)
    EndIndex = nint(error)

    !> If there is no overlap, exit setting Indexes to error
    if (StartTimestamp > Filelist(nrow)%timestamp &
        .or. EndTimestamp < Filelist(1)%timestamp) return

    !> Search StartTimestamp
    if (StartTimestamp < Filelist(1)%timestamp)  then
        StartIndex = 1
    else
        do i = 1, nrow
            if (FileList(i)%timestamp >= StartTimestamp) then
                StartIndex = i
                exit
            end if
        end do
    end if

    !> Search EndTimestamp
    if (EndTimestamp > FileList(nrow)%timestamp) then
        EndIndex = nrow
    else
        do i = StartIndex, nrow
            if (FileList(i)%timestamp >= EndTimestamp) then
                EndIndex = i - 1
                exit
            end if
        end do
    end if
end subroutine tsExtractSubperiodIndexesFromFilelist

