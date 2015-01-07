!***************************************************************************
! biomet_subs.f90
! ---------------
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
!***************************************************************************
!
! \brief       Collection of helper subs for handling biomet files/data
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
!***************************************************************************

!***************************************************************************
!
! \brief       Reads and interprets file header, searching for known variables
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
!***************************************************************************
subroutine RetrieveExtBiometVars(row1, row2, nitems)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nitems
    character(*), intent(inout) :: row1
    character(*), intent(inout) :: row2
    !> local variables
    integer :: i
    integer :: cnt
    integer :: tsCnt
    integer :: var_sep
    integer :: unit_sep
    character(32) :: cap_item1
    character(32) :: item1
    character(32) :: item2

    !> Retrieve variables and units from corresponding strings
    bFileMetadata%tsPattern = ''
    bFileMetadata%tsCols = nint(error)
    cnt = 0
    tsCnt = 0
    do i = 1, nitems
        var_sep = index(row1, bFileMetadata%separator)
        unit_sep = index(row2, bFileMetadata%separator)
        if (var_sep == 0) var_sep = len_trim(row1) + 1
        if (unit_sep == 0) unit_sep = len_trim(row2) + 1
        if (len_trim(row1) == 0) exit
        item1 = row1(1:var_sep - 1)
        cap_item1 = item1
        item2 = row2(1:unit_sep - 1)
        row1 = row1(var_sep + 1: len_trim(row1))
        row2 = row2(unit_sep + 1: len_trim(row2))
        call uppercase(cap_item1)
        if (index(cap_item1, 'TIMESTAMP') /= 0) then
            tsCnt = tsCnt + 1
            bFileMetadata%tsCols(tsCnt) = i
            bFileMetadata%tsPattern = &
                trim(adjustl(bFileMetadata%tsPattern)) &
                // trim(adjustl(item2))
        else
            cnt = cnt + 1
            bVars(cnt)%label = item1
            bVars(cnt)%unit_in = item2
            call ShrinkString(bVars(cnt)%label)
            call uppercase(bVars(cnt)%unit_in)
        end if
    end do
    bFileMetadata%numTsCol = tsCnt
end subroutine RetrieveExtBiometVars

!***************************************************************************
!
! \brief       Retrieve timestamp string from data record given var labels
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
!***************************************************************************
subroutine tsStringFromRec(row, nitems, tsString)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nitems
    character(*), intent(inout) :: row
    character(*), intent(out) :: tsString
    !> local variables
    integer :: i
    integer :: sepa
    character(32) :: item


    tsString = ''
    do i = 1, nitems
        !> Isolate item
        if (len(row) <= 0) exit
        sepa = index(row, bFileMetadata%separator)
        if (sepa == 0) then
            item = row(1:len_trim(row))
        else
            item = row(1:sepa-1)
        end if
        !> If item is part of timestamp add it to tsString
        if (any(bFileMetadata%tsCols == i)) &
            tsString = trim(tsString) // trim(adjustl(item))
        if (sepa == 0) exit

        !> Take away item from row
        row = row(sepa+1 : len_trim(row))
    end do
end subroutine tsStringFromRec

!***************************************************************************
!
! \brief       Retrieve timestamp info from timestamp strings in the file
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine BiometDateTime(pattern, string, date, time, failed)
    use m_rp_global_var
    implicit none
    !> in/out variables
    character(*), intent(in) :: pattern
    character(*), intent(in) :: string
    character(*), intent(out) :: date
    character(*), intent(out) :: time
    logical, intent(out) :: failed
    !> local variables
    integer :: loc_doy
    integer :: start
    integer :: int_year
    integer :: read_status
    logical :: isleap
    logical, external :: is_date, is_time


    date(5:5) = '-'
    date(8:8) = '-'
    time(3:3) = ':'
    failed = .false.

    !> If string is too short exit with error
    if (len_trim(string) < len_trim(pattern)) then
        failed = .true.
        return
    end if

    !> year
    if (index(Pattern, 'yyyy') /= 0) then
        start = index(Pattern, 'yyyy')
        date(1:4) = string(start: start + 3)
    else
        if (index(Pattern, 'yy') /= 0) then
            start = index(Pattern, 'yy')
            if (string(start: start + 1) > '70') then
                date(1:4) = '19' // string(start: start + 1)
            else
                date(1:4) = '20' // string(start: start + 1)
            end if
        else
            date(1:4) = 'xxxx'
        end if
    end if

    if (date(1:4) /= 'xxxx') then
        read(date(1:4), '(i4)', iostat=read_status) int_year
        if (read_status /= 0) then
            failed = .true.
            return
        end if
        isleap = leapyear(int_year)
    else
        isleap = .false.
    end if

    !> month
    if (index(Pattern, 'mm') /= 0) then
        start = index(Pattern, 'mm')
        date(6:7) = string(start : start + 3)
    else
        date(6:7) = 'xx'
    end if

    !> day or DOY
    if (index(Pattern, 'ddd') /= 0) then
        start = index(Pattern, 'ddd')
        call Char2Int(string(start:start + 2), loc_doy, 3)
        if (isleap) then
            date(6:10) = DayOfLeapYear(loc_doy)
        else
            date(6:10) = DayOfYear(loc_doy)
        end if
    elseif (index(Pattern, 'dd') /= 0) then
        start = index(Pattern, 'dd')
        date(9:10) = string(start : start + 1)
    end if

    !> hour
    if (index(Pattern, 'HH') /= 0) then
        start = index(Pattern, 'HH')
        time(1:2) = string(start : start + 1)
    else
        time(1:2) = 'xx'
    end if

    !> minute
    if (index(Pattern, 'MM') /= 0) then
        start = index(Pattern, 'MM')
        time(4:5) = string(start : start + 1)
    else
        time(4:5) = 'xx'
    end if

    !> Check if date and time where correctly read
    failed = .not. (is_date(date) .and. is_time(time))

end subroutine BiometDateTime

!***************************************************************************
!
! \brief       Retrieve time step intrinsic in Biomet file (in minutes)
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
!***************************************************************************
subroutine BiometTimestamp(pattern, string, timestamp, failed)
    use m_rp_global_var
    implicit none
    !> In/out variables
    character(*), intent(in) :: pattern
    character(*), intent(in) :: string
    type(DateType), intent(in) :: timestamp
    logical, intent(out) :: failed
    !> Local variables
    character(10) :: date
    character(5) :: time


    call BiometDateTime(pattern, string, date, time, failed)
    if (failed) return
    call DateTimeToDateType(date, time, timestamp)
end subroutine BiometTimestamp


!***************************************************************************
!
! \brief       Based on variable label, extract name, location, profile
!              replicate and various other properties
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
!***************************************************************************
subroutine BiometEnrichVarsDescription()
    use m_rp_global_var
    implicit none
    !> In/out variables
    !> Local variables
    integer :: i
    character(32) :: item


    do i = 1, nbVars
        call BiometInterpretLabels(bVars(i))
        item = trim(bVars(i)%name)
        call uppercase(item)
        select case(item)
            case('TA', 'TC', 'TBOLE', 'TBC', 'TR', 'TS')
                bVars(i)%nature = 'TEMPERATURE'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'K'
                bVars(i)%pretty_unit_out = '[K]'
            case('RH')
                bVars(i)%nature = 'RELATIVE_HUMIDITY'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = '%'
                bVars(i)%pretty_unit_out = '[%]'
            case('PA')
                bVars(i)%nature = 'PRESSURE'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'PA'
                bVars(i)%pretty_unit_out = '[Pa]'
            case('CO2')
                bVars(i)%nature = 'CONCENTRATION'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'UMOL+1M-2S-1'
                bVars(i)%pretty_unit_out = '[' // char(181) // 'mol/m^2s]'
            case('CH4', 'N2O', 'NO', 'NO2', 'CO', 'SO2', 'O3', 'NH3')
                bVars(i)%nature = 'CONCENTRATION'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'NMOL+1M-2S-1'
                bVars(i)%pretty_unit_out = '[nmol/m^2s]'
            case('H2O')
                bVars(i)%nature = 'CONCENTRATION'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'MMOL+1M-2S-1'
                bVars(i)%pretty_unit_out = '[mmol/m^2s]'
            case('RG', 'RN', 'RD', 'RR', 'R_UVA', 'R_UVB', 'LWIN', 'LWOUT', &
                'SWIN', 'SWOUT', 'SWBC', 'SWDIF')
                bVars(i)%nature = 'RADIATION'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'W+1M-2'
                bVars(i)%pretty_unit_out = '[W/m^2]'
            case('SHF')
                bVars(i)%nature = 'HEAT_FLUX'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'W+1M-2'
                bVars(i)%pretty_unit_out = '[W/m^2]'
            case('PPFD', 'PPFDD', 'PPFDR', 'PPFDBC', 'APAR')
                bVars(i)%nature = 'PHOTON_FLUX'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'UMOL+1M-2S-1'
                bVars(i)%pretty_unit_out = '[' // char(181) // 'mol/m^2s]'
            case('WS', 'MWS')
                bVars(i)%nature = 'SPEED'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'M+1S-1'
                bVars(i)%pretty_unit_out = '[m/s]'
            case('WD')
                bVars(i)%nature = 'ANGULAR_DIRECTION'
                bVars(i)%accumul_type = 'ANGULAR_AVERAGE'
                bVars(i)%unit_out = 'DEGREES'
                bVars(i)%pretty_unit_out = '[deg_past_North]'
            case('P', 'P_RAIN', 'P_SNOW')
                bVars(i)%nature = 'PRECIPITATION'
                bVars(i)%accumul_type = 'INTEGRATION'
                bVars(i)%unit_out = 'M'
                bVars(i)%pretty_unit_out = '[m]'
            case('LAI')
                bVars(i)%nature = 'INDEX'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'M+2M-2'
                bVars(i)%pretty_unit_out = '[m^2/m^2]'
            case('ALB')
                bVars(i)%nature = 'INDEX'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = '#'
                bVars(i)%pretty_unit_out = '[#]'
            case('SAPFLOW', 'STEMFLOW')
                bVars(i)%nature = 'FLOW'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'MMOL+1M-2S-1'
                bVars(i)%pretty_unit_out = '[mmol/m^2s]'
            case('SNOWD')
                bVars(i)%nature = 'LENGTH'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'M'
                bVars(i)%pretty_unit_out = '[m]'
            case('SWC')
                bVars(i)%nature = 'VOLUME_CONCENTRATION'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = 'M+3M-3'
                bVars(i)%pretty_unit_out = '[m^3/m^3]'
            case default
                bVars(i)%nature = 'UNKNOWN'
                bVars(i)%accumul_type = 'AVERAGING'
                bVars(i)%unit_out = trim(bVars(i)%unit_in)
                bVars(i)%pretty_unit_out = '[' // trim(bVars(i)%unit_in) // ']'
        end select
    end do
end subroutine

!***************************************************************************
!
! \brief
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
!***************************************************************************
subroutine BiometAppendLocationSuffix()
    use m_rp_global_var
    implicit none
    !> Local variables
    integer :: i, ii
    integer :: cnt
    integer :: numVarsFound(MaxNumBiometVars)
    character(32) :: label
    character(8) :: loc
    character(32) :: varsFound(MaxNumBiometVars)
    logical, external :: BiometVarHasSuffix


    numVarsFound = 0
    varsFound = ''
    cnt = 0
    do i = 1, nbVars
        label = trim(bVars(i)%label)
        if (.not. BiometVarHasSuffix(label)) then
            !> This if block updates list of variables found
            if (.not. any(varsFound == label)) then
                cnt = cnt + 1
                VarsFound(cnt) = label
                numVarsFound(cnt) = 1
            else
                do ii = 1, cnt
                    if (varsFound(ii) == label) then
                        numVarsFound(cnt) = numVarsFound(cnt) + 1
                        exit
                    end if
                end do
            end if
            write(loc, '(i8)') numVarsFound(cnt)
            bVars(i)%label = trim(bVars(i)%label) // '_' &
                // trim(adjustl(loc)) // '_1_1'
        end if
    end do
end subroutine BiometAppendLocationSuffix


!***************************************************************************
!
! \brief
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
!***************************************************************************
logical function BiometVarHasSuffix(label)
    use m_rp_global_var
    implicit none
    !> In/out variables
    character(*), intent(in) :: label
    !> Local variables
    integer :: i
    integer :: ix, ix2
    integer :: int_num
    integer :: io_status
    character(32) :: token, num

    !> Check if variable name contains 3 underscores and text between them
    !> is an integer number.
    token = trim(label)
    BiometVarHasSuffix = .true.
    do i = 1, 3
        ix = index(token, '_')
        if (ix == 0 .or. ix == len_trim(token)) then
            BiometVarHasSuffix = .false.
            return
        else
            token = token(ix+1:len_trim(token))
            ix2 = index(token, '_')
            if (i < 3) then
                if (ix2 == 0) then
                    BiometVarHasSuffix = .false.
                    return
                else
                    num = token(1:ix2-1)
                end if
            else
                num = token(1:len_trim(token))
            end if
            read(num, *, iostat = io_status) int_num
            if (io_status /= 0) then
                BiometVarHasSuffix = .false.
                return
            end if
        end if
    end do
    return

end function BiometVarHasSuffix


!***************************************************************************
!
! \brief       Order biomet variables by variable label and by profile
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
!***************************************************************************
!subroutine BiometOrderVars()
!    use m_rp_global_var
!    implicit none
!    !> In/out variables
!    !> Local variables
!    integer :: i
!    integer :: last_ord
!    integer :: new_label
!    integer :: matching_tags
!    integer :: new_match
!    character(32) :: unique_tags(size(bVars))
!    character(32) :: tag
!    integer :: ord(size(bVars))
!    type(BiometVarsType) :: boVars(size(bVars))
!
!
!    !> Interpret location information  and store it in bVars attributes
!    do i = 1, size(bVars)
!        call BiometInterpretLabels(bVars(i))
!    end do
!
!    unique_tags = ''
!    new_label = 0
!    ord = 0
!    last_ord = 0
!    do i = 1, size(bVars)
!        tag = trim(bVars(i)%label)
!        matching_tags = 0
!        new_match = 0
!        if (not any(unique_tags, tag)) then
!            new_label = new_label + 1
!            unique_tags(new_label) = tag
!            !> For each new label, scan remaining variables
!            !> to find matching ones
!            do ii = i+1, size(bVars)
!                if (trim(bVars(ii)%label) == tag) then
!                    new_match = new_match + 1
!                    matching_tags(new_match) = ii
!                end if
!            end do
!            !> Now that matching vars have been found, order them by
!            !> repetition, profile and location (in this order)
!            if (any(matching_tags /= 0)) then
!                BiometOrderByLocation(bVars, size(bVars), i, matching_tags, &
!                    last_ord, ord)
!            end if
!        else
!
!        end if

!        if (any(aux, trim(bVars(i)%label))) then
!            else
!                new_label = new_label + 1
!                aux(new_label) = trim(bVars(i)%label)
!        end if
!    end do
!
!    boVars = bVars(ord)
!
!
!    do i = 1, size(boVars)
!        write(*,*) i, trim(boVars(i)%label)
!    end do
!
!    stop
!end subroutine BiometOrderVars


!***************************************************************************
!
! \brief       Parse biomet row and return timestamp and data
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
!***************************************************************************
subroutine BiometParseRow(row, tstamp, vals, ncol, skip_row)
    use m_rp_global_var
    implicit none
    !> In/out variables
    integer, intent(in) :: ncol
    real(kind=dbl), intent(out) :: vals(ncol)
    type(DateType), intent(out) :: tstamp
    logical, intent(out) :: skip_row
    character(*), intent(inout) :: row
    !> Local variables
    integer :: j, jj
    character(64) :: tsString
    character(32) :: item
    character(1)  :: sepa


    skip_row = .false.
    sepa = bFileMetadata%separator

    !> Check that row contains something
    if (len_trim(row) <= len_trim(bFileMetadata%data_label)) then
        skip_row = .true.
        return
    else
        !> Skip label if present
        row = row(len_trim(bFileMetadata%data_label)+1: len_trim(row))
    end if

    tsString = ''
    jj = 0
    do j = 1, nbVars + bFileMetadata%numTsCol
        if (index(row, sepa) /= 0) then
            item = row(1: index(row, sepa)-1)
        else
            item = row(1: len_trim(row))
        end if

        !> If item is empty, something went wrong, so
        !> exit with error
        if (j < nbVars + bFileMetadata%numTsCol .and. len_trim(item) == 0) then
            skip_row = .true.
            return
        end if

        !> Eliminate item from row
        row = row(len_trim(item)+2: len_trim(row))

        !> If row is now empty, something went wrong (row was too short), so
        !> exit with error
        if (j < nbVars + bFileMetadata%numTsCol .and. len_trim(row) == 0) then
            skip_row = .true.
            return
        end if

        !> From item, extract timestamp or biomet value as appropriate
        if (any(bFileMetadata%tsCols(1:bFileMetadata%numTsCol) == j)) then
            tsString = trim(tsString) // trim(item)
        else
            jj = jj + 1
            read(item, *) vals(jj)
        end if
    end do

    !> Retrieve timestamp from timestamp-string
    call FilenameToTimestamp(trim(tsString), &
        bFileMetadata%tsPattern, EddyProLog%iso_format, tstamp)

end subroutine BiometParseRow

!***************************************************************************
!
! \brief       Change the
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine BiometUpdateSelectionOrder()
    use m_rp_global_var
    implicit none
    !> Local variables
    integer :: i, bix
    integer :: origSel(6)


    origSel = bSetup%sel
    do i = 1, bFileMetadata%numTsCol
        do bix = bTa, bRg
            if (bFileMetadata%tsCols(i) < origSel(bix)) &
                bSetup%sel(bix) = bSetup%sel(bix) - 1
        end do
    end do
end subroutine BiometUpdateSelectionOrder

!***************************************************************************
!
! \brief       Take timestamp of biomet record to end of biomet averaging period
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine BiometAdjustTimestamp(timestamp)
    use m_rp_global_var
    implicit none
    !> in/out variables
    type(DateType), intent(inout) :: timestamp


    select case (bFileMetadata%tstamp_ref)
        case ('begin')
            timestamp = timestamp &
                + datetype(0, 0, 0, 0, bFileMetadata%time_step)
        case ('middle')
            timestamp = timestamp &
                + datetype(0, 0, 0, 0, bFileMetadata%time_step/2)
        case ('end')
        return
    end select
end subroutine BiometAdjustTimestamp

!***************************************************************************
!
! \brief       Interpret label and return var, loc, profile, rep
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
!***************************************************************************
subroutine BiometInterpretLabels(bVar)
    use m_rp_global_var
    implicit none
    !> In/out variables
    type(BiometVarsType), intent(inout) :: bVar
    !> Local variables
    character(32) :: s
    integer :: i
    integer :: n
    integer, external :: CountCharInString

    s = bvar%label
    !> Count number of underscores
    n = CountCharInString(s, '_')

    !> If n is less than 3, the variable name was not well defined, so
    !> abort using biomet
    if (n < 3) then
        call ExceptionHandler(73)
        EddyProProj%biomet_data = 'none'
        return
    end if

    !> Variable name
    !> If n is more than 3, var name contains the remaining underscores
    bVar%name = ''
    if (n > 3) then
        do i = 1, n-2
            bVar%name = trim(bVar%name) // s(1:index(s, '_'))
            s = s(index(s, '_') + 1: len_trim(s))
        end do
        bVar%name = bVar%name(1:len_trim(bVar%name)-1)
    else
        bVar%name = s(1:index(s, '_')-1)
        s = s(index(s, '_') + 1 : len_trim(s))
    end if

    !> location
    read(s(1:index(s, '_')-1), '(i3)') bVar%loc
    !> location
    s = s(index(s, '_') + 1 : len_trim(s))
    read(s(1:index(s, '_')-1), '(i3)') bVar%profile
    !> rep
    s = s(index(s, '_') + 1 : len_trim(s))
    read(s(1:len_trim(s)), '(i3)') bVar%rep
end subroutine BiometInterpretLabels

!***************************************************************************
!
! \brief       Convert input units into standard units
! \author      Gerardo Fratini
! \note
!              Radiations (Rg, Rn, Rd, Rr, LWin, LWout, Ruva, Ruvb) \n
!              are not expected to need unit conversion
!              Photons flux densities (PPFD, PPFDd, PPFDr, PPFDbc, APAR) \n
!              are not expected to need unit conversion
!              Alb is not expected to need unit conversion
!              PRI is not expected to need unit conversion
!              SWC is not expected to need unit conversion
!              SHF is not expected to need unit conversion
! \sa
! \bug
! \deprecated
! \test
!***************************************************************************
subroutine BiometStandardUnits()
    use m_rp_global_var
    implicit none
    !> local variables
    integer :: i

    !> Temperatures
    do i = 1, nbVars
        select case(trim(bVars(i)%nature))
            case('TEMPERATURE')
                select case(bVars(i)%unit_in)
                    case('C','�C')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) + 273.16d0
                        end where
                    case('F','�F')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = (bSet(:, i) - 32d0) * 5d0 / 9d0 &
                                + 273.16d0
                        end where
                    case('CK')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 1d-2
                        end where
                    case('CC','C�C')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 1d-2 + 273.16d0
                        end where
                    case('CF','C�F')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = (bSet(:, i) * 1d-2 - 32d0) * 5d0 / 9d0 &
                                + 273.16d0
                        end where
                    case default
                end select

            case('RELATIVE_HUMIDITY')
                select case(bVars(i)%unit_in)
                    case('NUMBER','#','DIMENSIONLESS')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 1d2
                        end where
                    case default
                        continue
                end select

            case('PRESSURE')
                select case(bVars(i)%unit_in)
                    case('HPA')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 1d2
                        end where
                    case('KPA')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 1d3
                        end where
                    case('MMHG', 'TORR')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 133.32d0
                        end where
                    case('PSI')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 6894.6d0
                        end where
                    case('BAR')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 1d5
                        end where
                    case('ATM')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 0.980665d5
                        end where
                    case default
                        continue
                end select

!            !> Precipitation is converted to [m]
!            case('PRECIPITATION')
!                select case(bVars(i)%unit_in)
!                    case('NM')
!                        where (bSet(:, i) /= error)
!                            bSet(:, i) = bSet(:, i) * 1d-6
!                        end where
!                    case('UM')
!                        where (bSet(:, i) /= error)
!                            bSet(:, i) = bSet(:, i) * 1d-3
!                        end where
!                    case('CM')
!                        where (bSet(:, i) /= error)
!                            bSet(:, i) = bSet(:, i) * 10
!                        end where
!                    case('M')
!                        where (bSet(:, i) /= error)
!                            bSet(:, i) = bSet(:, i) * 1d3
!                        end where
!                    case default
!                        continue
!                end select

            !> Lengths
            !> converted to [m]
            case('LENGTH', 'PRECIPITATION')
                select case(bVars(i)%unit_in)
                    case('NM')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 1d-9
                        end where
                    case('UM')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 1d-6
                        end where
                    case('MM')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 1d-3
                        end where
                    case('CM')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 1d-2
                        end where
                    case('KM')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 1d3
                        end where
                    case default
                        continue
                end select

            case('CONCENTRATION')
                select case(trim(adjustl(bVars(i)%label)))
                !> CO2 is converted to PPM
                case ('CO2')
                    select case(bVars(i)%unit_in)
                        case('PPT', 'MMOL/MOL', 'MMOL+1MOL-1', 'MMOLMOL-1')
                            where (bSet(:, i) /= error)
                                bSet(:, i) = bSet(:, i) * 1d3
                            end where
                        case('PPB', 'NMOL/MOL', 'NMOL+1MOL-1', 'NMOLMOL-1')
                            where (bSet(:, i) /= error)
                                bSet(:, i) = bSet(:, i) * 1d-3
                            end where
                        case default
                            continue
                    end select
                !> H2O is converted to PPT
                case ('H2O')
                    select case(bVars(i)%unit_in)
                        case('PPM', 'UMOL/MOL', 'UMOL+1MOL-1', 'UMOLMOL-1')
                            where (bSet(:, i) /= error)
                                bSet(:, i) = bSet(:, i) * 1d-3
                            end where
                        case('PPB', 'NMOL/MOL', 'NMOL+1MOL-1', 'NMOLMOL-1')
                            where (bSet(:, i) /= error)
                                bSet(:, i) = bSet(:, i) * 1d-6
                            end where
                        case default
                            continue
                    end select
                !> Any other gas is converted to PPB
                case default
                    select case(bVars(i)%unit_in)
                        case('PPT', 'MMOL/MOL', 'MMOL+1MOL-1', 'MMOLMOL-1')
                            where (bSet(:, i) /= error)
                                bSet(:, i) = bSet(:, i) * 1d6
                            end where
                        case('PPM', 'UMOL/MOL', 'UMOL+1MOL-1', 'UMOLMOL-1')
                            where (bSet(:, i) /= error)
                                bSet(:, i) = bSet(:, i) * 1d3
                            end where
                        case default
                            continue
                    end select
                end select

            case('SPEED')
                select case(bVars(i)%unit_in)
                    case('CM+1S-1','CM/S','CMS^-1','CMS-1')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 1d-2
                        end where
                    case('MM+1S-1','MM/S','MMS^-1','MMS-1')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 1d-3
                        end where
                    case default
                        continue
                end select

            case('ANGULAR_DIRECTION')
            case('FLOW')
        end select
    end do
end subroutine BiometStandardUnits
!***************************************************************************
!
! \brief       Interpret label and return var, loc, profile, rep
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
!***************************************************************************
!subroutine BiometOrderByLocation(bVar, bVars, n, init, matching_tags, &
!                    last_ord, ord)
!    use m_rp_global_var
!    implicit none
!    !> In/out variables
!    integer, intent(in) :: n, init, last_ord, matching_tags
!    type(BiometVarsType), intent(in) :: bVar
!    integer, intent(inout) :: ord
!    !> Local variables
!    integer :: i

!    do i = init, n
!
!end do
!
!end subroutine BiometOrderByLocation


!***************************************************************************
!
! \brief       Aggregate biomet data according to aggregation method
!              (sum, average, etc..)
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
!***************************************************************************
subroutine BiometAggretate(Set, nrow, ncol, Aggr)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    real(kind = dbl), intent(out) :: Aggr(ncol)
    !> local variables
    integer :: i
    real(kind = dbl) :: vAggr

    !> Aggregate biomet variables over the averaging interval
    do i = 1, size(Set, 2)
        select case(trim(adjustl(bVars(i)%accumul_type)))
            !> Variables that need averaging
            case ('AVERAGING')
                call AverageNoError(bSet(:, i), size(bSet, 1), 1, vAggr, error)
            case ('INTEGRATION')
                call SumNoError(bSet(:, i), size(bSet, 1), 1, vAggr, error)
            case default
        end select
        Aggr(i) = vAggr
    end do
end subroutine BiometAggretate


!***************************************************************************
!
! \brief       Put biomet file names in chronological order in BiometFileList
!              (sum, average, etc..)
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
!***************************************************************************
subroutine BiometFileListInChronologicalOrder(FileList, nrow)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    type (FilelistType), intent(inout) :: FileList(nrow)
    !> local variables
    integer :: i
    integer :: nfl
    integer :: cnt
    integer :: io_status
    integer :: rank(nrow)
    character(64) :: tsString
    character(1024) :: record
    logical :: skip_record
    type(DateType) :: bTimestamp
    type(FileListType) :: TmpFileList(nrow)


    TmpFileList= FileListType('none', 'none', nullTimestamp)

    !> Associate a timestamp to each file, by reading the first valid
    !> record and interpreting the timestamp
    files_loop: do nfl = 1, nrow
        !> Open biomet file
        open(udf, file=FileList(nfl)%path, status = 'old', iostat=io_status)
        !> If a problem occurred while opening file, tag file with
        !> default timestamp so that it can be later eliminated from list
        if (io_status /= 0) then
            FileList(nfl)%timestamp = nullTimestamp
            cycle files_loop
        end if

        !> Retrieve timestamp from first valid record
        recs_loop: do
            read(udf, '(a)', iostat = io_status) record
            if (io_status > 0) cycle recs_loop
            if (io_status < 0) then
                FileList(nfl)%timestamp = nullTimestamp
                cycle files_loop
            end if

            !> Retrieve timestamp info from record
            call tsStringFromRec(record, nbItems, tsString)
            call BiometTimestamp(trim(adjustl(bFileMetadata%tsPattern)), &
                    tsString, bTimestamp, skip_record)
            if (skip_record) cycle recs_loop

            !> Associate timestamp to file in FileList
            FileList(nfl)%timestamp = bTimestamp
            cycle files_loop
        end do recs_loop
    end do files_loop

    !> Sort timestamps in a chronological sequence
    call rank_dates(FileList%timestamp, rank, nrow)

    cnt = 0
    do i = 1, nrow
        if (FileList(rank(i))%timestamp /= nullTimestamp) then
            cnt = cnt + 1
            TmpFileList(cnt) = FileList(rank(i))
        end if
    end do
    FileList = TmpFileList
end subroutine BiometFileListInChronologicalOrder


