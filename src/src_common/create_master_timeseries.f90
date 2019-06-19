!***************************************************************************
! create_master_timeseries.f90
! ----------------------------
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
! \brief       Create timestamp array based on start/end timestamp and step
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine CreateTimeSeries(StartTimestamp, EndTimestamp, &
    Step, RawTimeSeries, nrow, printout)
    use m_common_global_var
    implicit none
    !> In/out variables
    integer, intent(in) :: nrow
    logical, intent(in) :: printout
    type(DateType), intent(in) :: StartTimestamp
    type(DateType), intent(inout) :: EndTimestamp
    type(DateType), intent(in) :: Step
    type(DateType), intent(out) :: RawTimeSeries(nrow)
    !> in/out variables
    integer :: cnt


    if (printout) write(*, '(a)', advance = 'no') ' Creating master time series..'

    !> create master timestamps array
    RawTimeSeries(1) = StartTimestamp
    cnt = 1
    do
        if (RawTimeSeries(cnt) >= EndTimestamp) exit
        cnt = cnt + 1
        if (cnt > nrow) exit
        RawTimeSeries(cnt) = RawTimeSeries(cnt - 1) + Step
    end do

    if (printout) write(*, '(a)') ' Done.'
end subroutine CreateTimeSeries

!***************************************************************************
!
! \brief       Calculate number of periods in time series
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
integer function NumOfPeriods(StartTimestamp, EndTimestamp, Step)
    use m_common_global_var
    implicit none
    !> In/out variables
    type(DateType), intent(in) :: StartTimestamp
    type(DateType), intent(inout) :: EndTimestamp
    type(DateType), intent(in) :: Step
    !> Local variables
    integer :: cnt
    type (DateType) :: Timestamp


    Timestamp = StartTimestamp
    cnt = 1
    do
        if (Timestamp >= EndTimestamp) exit
        cnt = cnt + 1
        Timestamp = Timestamp + Step
    end do
    NumOfPeriods = cnt - 1
end function NumOfPeriods
