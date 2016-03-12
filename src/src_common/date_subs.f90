!***************************************************************************
! date_subs.f90
! -------------
! Copyright (C) 2007-2011, Eco2s team, Gerardo Fratini
! Copyright (C) 2011-2015, LI-COR Biosciences
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
!***************************************************************************
!
! \brief       Collection of subs for handling time stamps
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************

!***************************************************************************
!
! \brief       Weak check of compatibility of provided raw file prototype
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine tsValidateTemplate(Template)
    implicit none
    !> local variables
    character(*), intent(in) :: Template

    !> Weak test
    if ( index(Template, 'yy') == 0 &
    .or. index(Template, 'dd') == 0 &
    .or. index(Template, 'HH') == 0 &
    .or. index(Template, 'MM') == 0) call ExceptionHandler(20)
end subroutine tsValidateTemplate

!***************************************************************************
!
! \brief       Subtracts a date-step to a starting date
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine SubtractDateStep(in_date, in_time, out_date, out_time, Step)
    use m_common_global_var
    implicit none
    !> in/out variables
    character(*), intent(in) :: in_date
    character(*), intent(in) :: in_time
    type(DateType), intent(in) :: Step
    character(*), intent(out) :: out_date
    character(*), intent(out) :: out_time
    !> local variables
    type(DateType) :: LateDate
    type(DateType) :: EarlyDate

    call DateTimeToDateType(in_date, in_time, LateDate)
    EarlyDate = LateDate - Step
    call DateTypeToDateTime(EarlyDate, out_date, out_time)
end subroutine SubtractDateStep

!***************************************************************************
!
! \brief       Adds a date-step to a starting date
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine AddDateStep(in_date, in_time, out_date, out_time, Step)
    use m_common_global_var
    implicit none
    !> in/out variables
    character(*), intent(in) :: in_date
    character(*), intent(in) :: in_time
    type(DateType), intent(in) :: Step
    character(*), intent(out) :: out_date
    character(*), intent(out) :: out_time
    !> local variables
    type(DateType) :: LateDate
    type(DateType) :: EarlyDate

    call DateTimeToDateType(in_date, in_time, EarlyDate)
    LateDate = EarlyDate + Step
    call DateTypeToDateTime(LateDate, out_date, out_time)
end subroutine AddDateStep

!***************************************************************************
!
! \brief       Retrieve date and time from a date string
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine DateStringToDateTime(Datestring, doy_format, date, time)
    use m_common_global_Var
    implicit none
    !> in/out variables
    character(*), intent(in) :: Datestring
    logical, intent(in) :: doy_format
    character(*), intent(out) :: date
    character(*), intent(out) :: time
    !> local variables
    integer :: int_year
    integer :: int_doy
    logical :: isleap

    !> Deterime whether current year is leap
    read(Datestring(1:4), '(i4)') int_year
    isleap = leapyear(int_year)

    time = Datestring(10:11) // ':' // Datestring(12:13)
    if (doy_format) then
        !> If date is in ISO format, date and time are straightforward
        date = Datestring(1:4) // '-' // Datestring(5:6) // '-' // Datestring(7:8)
    else
        !> if date in not in ISO format, the actual 'mm-dd' must be derived from DOY
        if (isleap) then
            read(Datestring(5:7), '(i3)') int_doy
            date = Datestring(1:4) // '-' // DayOfLeapYear(int_doy)(1:5)
        else
            read(Datestring(5:7), '(i3)') int_doy
            date = Datestring(1:4) // '-' // DayOfYear(int_doy)(1:5)
        end if
    end if

end subroutine DateStringToDateTime

!***************************************************************************
!
! \brief       Convert date/time pair to a timestamp
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine DateTimeToDateType(date, time, ld_date)
    use m_common_global_Var
    implicit none
    !> in/out variables
    character(*), intent(in) :: date
    character(*), intent(in) :: time
    Type(DateType), intent(out) :: ld_date
    !> local variables
    integer :: loc_year
    integer :: loc_month
    integer :: loc_day
    integer :: loc_hour
    integer :: loc_min


    read(date(1:4),  '(i4)') loc_year
    read(date(6:7),  '(i2)') loc_month
    read(date(9:10), '(i2)') loc_day
    read(time(1:2),  '(i2)') loc_hour
    read(time(4:5),  '(i2)') loc_min
    ld_date = datetype(loc_year, loc_month, loc_day, loc_hour, loc_min)
end subroutine DateTimeToDateType

!***************************************************************************
!
! \brief       Retrieve timestamp from Filename based on Prototype
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine FilenameToTimestamp(Filename, Prototype, doy_format, Timestamp)
    use m_common_global_var
    implicit none
    !> in/out variables
    character(*), intent(in) :: Filename
    character(*), intent(in) :: Prototype
    type (DateType), intent(out) :: Timestamp
    logical, intent(in) :: doy_format
    !> local variables
    character(13) :: DateString
    character(10) :: date
    character(5) :: time
    integer :: loc_year
    integer :: loc_month
    integer :: loc_day
    integer :: loc_hour
    integer :: loc_min


    !> extract date/time info from file name
    call ParseFileNameWithTemplate(Filename, Prototype, DateString)

    !> If midnights are expressed as 24, change it into 00:00 of day after
    call Change24Into00(DateString, doy_format)

    !> Define timestamp to which filename refers
    call DateStringToDateTime(DateString, doy_format, date, time)

    read(date(1:4),  '(i4)') loc_year
    read(date(6:7),  '(i2)') loc_month
    read(date(9:10), '(i2)') loc_day
    read(time(1:2),  '(i2)') loc_hour
    read(time(4:5),  '(i2)') loc_min
    Timestamp = datetype(loc_year, loc_month, loc_day, loc_hour, loc_min)
end subroutine FilenameToTimestamp

!***************************************************************************
!
! \brief       Retrieve date-time string Filename based on Template
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine FilenameToDateTime(Filename, Template, doy_format, date, time)
    use m_common_global_var
    implicit none
    !> in/out variables
    character(*), intent(in) :: Filename
    character(*), intent(in) :: Template
    logical, intent(in) :: doy_format
    character(*), intent(out) :: date
    character(*), intent(out) :: time
    !> local variables
    character(64) :: DateString


    !> extract date/time info from file name
    call ParseFileNameWithTemplate(Filename, Template, DateString)

    !> If midnights are expressed as 24, change it into 00:00 of day after
    call Change24Into00(trim(adjustl(DateString)), doy_format)

    !> Define timestamp to which filename refers
    call DateStringToDateTime(DateString, doy_format, date, time)
end subroutine FilenameToDateTime

!***************************************************************************
!
! \brief       Convert a timestamp into a date/time pair
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine DateTypeToDateTime(ld_date, date, time)
    use m_common_global_Var
    implicit none
    !> in/out variables
    Type(DateType), intent(in) :: ld_date
    character(*), intent(out) :: date
    character(*), intent(out) :: time
    !> local variables
    character(32) :: datestring
    character(32) :: timestring
    character(32) :: pattern


    pattern = 'yyyymmdd'
    call format_date(ld_date, pattern, datestring)
    pattern = 'HHMM'
    call format_date(ld_date, pattern, timestring)
    date = datestring(1:4) // '-' // datestring(5:6) // '-' // datestring(7:8)
    time = timestring(1:2) // ':' // timestring(3:4)
end subroutine DateTypeToDateTime

!***************************************************************************
!
! \brief       Convert a date string into a timestamp
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine DateStringToTimestamp(Datestring, doy_format, ld_date)
    use m_common_global_Var
    implicit none
    !> in/out variables
    character(*), intent(in) :: Datestring
    logical, intent(in) :: doy_format
    Type(DateType), intent(out) :: ld_date
    !> local variables
    character(10) :: date
    character(5) :: time


    call DateStringToDateTime(Datestring, doy_format, date, time)
    call DateTimeToDateType(date, time, ld_date)
end subroutine DateStringToTimestamp

!***************************************************************************
!
! \brief       For each half hour of the year determines if it is day or \n
!              night-time, based on provided potential radiation
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
logical function IsDaytime(rad, date, time)
    use m_common_global_Var
    implicit none
    !> in/out variables
    real(kind = dbl), intent(in) :: rad(17568)
    character(*), intent(in) :: date
    character(*), intent(in) :: time
    !> local variables
    logical :: isleap
    integer :: int_year
    integer :: int_day
    integer :: int_hour
    integer :: int_min
    integer :: indx


    !> Deterime whether current year is leap
    read(date(1:4), '(i4)') int_year
    isleap = leapyear(int_year)

    indx = 0
    !> Determine closest half hour
    !> First handles the month
    select case(date(6:7))
        case('01')
            indx = 0
        case('02')
            indx = 1488
        case('03')
            if (isleap)       indx = 2880
            if (.not. isleap) indx = 2832
        case('04')
            if (isleap)       indx = 4368
            if (.not. isleap) indx = 4320
        case('05')
            if (isleap)       indx = 5808
            if (.not. isleap) indx = 5760
        case('06')
            if (isleap)       indx = 7296
            if (.not. isleap) indx = 7248
        case('07')
            if (isleap)       indx = 8736
            if (.not. isleap) indx = 8688
        case('08')
            if (isleap)       indx = 10224
            if (.not. isleap) indx = 10176
        case('09')
            if (isleap)       indx = 11712
            if (.not. isleap) indx = 11664
        case('10')
            if (isleap)       indx = 13152
            if (.not. isleap) indx = 13104
        case('11')
            if (isleap)       indx = 14640
            if (.not. isleap) indx = 14592
        case('12')
            if (isleap)       indx = 16080
            if (.not. isleap) indx = 16032
    end select

    !> Handle the days before the current one
    read(date(9:10), '(i2)') int_day
    indx = indx + (int_day - 1) * 48

    !> Handle the hours before current one
    read(time(1:2), '(i2)') int_hour
    indx = indx + int_hour * 2

    !> Handle the minutes
    read(time(4:5), '(i2)') int_min
    select case (int_min)
        case(0:15)
            continue
        case(16:45)
            indx = indx + 1
        case(46:59)
            indx = indx + 2
    end select
    if (indx == 0) indx = 1

    !> Now indx is known, use relevant radiation value to determine daytime
    if (rad(indx) > 10d0) then
        IsDaytime = .true.
    else
        IsDaytime = .false.
    end if
end function IsDaytime

!***************************************************************************
!
! \brief       Calculate decimal DOY from date-time pair
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine DateTimeToDOY(date, time, int_doy, float_doy)
    use m_common_global_Var
    implicit none
    !> in/out variables
    character(*), intent(in) :: date
    character(*), intent(in) :: time
    integer, intent(out) :: int_doy
    real(kind = dbl), intent(out) :: float_doy
    !> local variables
    type(datetype) :: now


    call DateTimeToDateType(date, time, now)
    int_doy = doy(now)
    float_doy = dfloat(int_doy) + 0.02083d0 * (2d0 * now%hour) + 0.02083d0/30d0 * (now%minute)
end subroutine DateTimeToDOY

!***************************************************************************
!
! \brief       Returns whether provided timestamp is within a provided
!              time range (initial + delta)
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine tsRelaxedMatch(tsTest, tsList, nrow, tsRange, side, imatch)
    use m_common_global_Var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    type(DateType), intent(in) :: tsTest
    type(DateType), intent(in) :: tsList(nrow)
    type(DateType), intent(in) :: tsRange
    character(*), intent(in) :: side
    integer, intent(out) :: imatch
    !> Local variables
    integer :: i


    imatch = -1

    select case (trim(adjustl(side)))
        case ('later')
            do i = 1, nrow
                if (tsList(i) <= tsTest .and. tsTest <= tsList(i) + tsRange) then
                    imatch = i
                    return
                end if
            end do
        case ('strictly later')
            do i = 1, nrow
                if (tsList(i) < tsTest .and. tsTest <= tsList(i) + tsRange) then
                    imatch = i
                    return
                end if
            end do
        case ('before')
            do i = 1, nrow
                if (tsList(i) - tsRange <= tsTest.and. tsTest <= tsList(i)) then
                    imatch = i
                    return
                end if
            end do
        case ('strictly before')
            do i = 1, nrow
                if (tsList(i) - tsRange <= tsTest.and. tsTest < tsList(i)) then
                    imatch = i
                    return
                end if
            end do
        case ('either')
            do i = 1, nrow
                if (tsList(i) - tsRange <= tsTest.and. tsTest <= tsList(i) + tsRange) then
                    imatch = i
                    return
                end if
            end do
    end select

end subroutine tsRelaxedMatch


!***************************************************************************
!
! \brief       Infer time step (in seconds) from array of timestamps
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
!***************************************************************************
integer function tsInferTimestep(timestamps, nrow) result(tstep)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    type(DateType) :: timestamps(nrow)
    !> local variables
    integer :: i
    integer :: step
    !integer, allocatable :: Steps(:), tmpSteps(:)
    integer :: maxTimes
    integer :: nsteps

    type StepsType
        integer :: v
        integer :: cnt
    end type StepsType
    type(StepsType), allocatable :: Steps(:), tmpSteps(:)
    type(StepsType), parameter :: nullSteps = StepsType(0, 0)


    !> First, count how many different time steps are in the array
    do i = 2, nrow
        !> Calculate current time step
        step = nint(timelag(timestamps(i), timestamps(i-1)) * 24d0 * 60d0 * 60d0)

        !> Special case of first value
        if (i == 2) then
            allocate(Steps(1))
            Steps = nullSteps
            Steps(1)%v = step
            Steps(1)%cnt = Steps(1)%cnt + 1
        end if

        !> Check if step is new or already contained in Steps
        if (.not. any(Steps%v == step)) then
            !> A bit involved way to extend size of array
            allocate(tmpSteps(size(Steps)+1))
            tmpSteps(1:size(Steps)) = Steps
            deallocate(Steps)
            !move_alloc(tmpSteps, Steps)
            allocate(Steps(size(tmpSteps)))
            Steps=tmpSteps
            deallocate(tmpSteps)
            !> Add new step to array
            Steps(size(Steps):size(Steps))%v = step
            Steps(size(Steps):size(Steps))%cnt = 1
        else
            where(Steps%v == step)
                Steps%cnt = Steps%cnt + 1
            end where
        end if
    end do
    nsteps = size(Steps)

    !> Find the occurrences of the most recurring step
    maxTimes = maxval(Steps%cnt)

    tstep = nint(error)
    do i = 1, maxTimes
        if (Steps(i)%cnt .eq. maxTimes) then
            tstep = Steps(i)%v
            exit
        endif
    end do
    if (allocated(Steps)) deallocate(Steps)
end function tsInferTimestep
