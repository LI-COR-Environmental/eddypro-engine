!***************************************************************************
! biomet_subs.f90
! ---------------
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
    character(32), external :: biometBaseName
    logical, external :: BiometValidateVar


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
            !> Variable label
            bVars(cnt)%label = item1
            call ShrinkString(bVars(cnt)%label)
            !> Variable units
            bVars(cnt)%unit_in = item2
            call uppercase(bVars(cnt)%unit_in)


            ! !> Check validity of biomet variable label
            ! if (.not. BiometValidateVar(bVars(cnt))) then
            !     write(*,'(a)')
            !     call ExceptionHandler(73)
            !     EddyProProj%biomet_data = 'none'
            !     return
            ! end if

            !> Retrieve variable base name
            bVars(cnt)%base_name = biometBaseName(bVars(cnt)%label)
        end if
    end do

    !> Validate timestamp pattern
    call tsValidateTemplate(bFileMetadata%tsPattern)

    !> Set ISO timestamp or not
    if (index(bFileMetadata%tsPattern, 'ddd') /= 0) then
        bFileMetadata%tsIso = .false.
    else
        bFileMetadata%tsIso = .true.
    end if

    bFileMetadata%numTsCol = tsCnt

    !> Append suffix if variables have not
    call BiometAppendReplicateSuffix()

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
! \brief       Basic validation of biomet labels
! \author      Gerardo Fratini
! \note        The bWaved machinery isn't active yet
! \sa
! \bug
! \deprecated
! \test
!***************************************************************************
logical function BiometValidateVar(bVar) result(valid)
    use m_rp_global_var
    implicit none
    !> In/out variables
    type(BiometVarsType), intent(in) :: bVar
    !> Local variables
    integer, external :: CountCharInString
    character(32) :: bWaved(10)
    data bWaved(1:1) /'LOGGERTEMP'/


    valid = .true.
    if (any(bWaved == bVar%label) .or. len_trim(bVar%label) == 0) return
    if (CountCharInString(bVar%label, '_') < 3) then
        valid = .false.
        return
    end if
end function BiometValidateVar

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
subroutine BiometAppendReplicateSuffix()
    use m_rp_global_var
    implicit none
    !> Local variables
    integer :: i, ii
    integer :: cnt
    integer :: numVarsFound(nbVars)
    character(32) :: label
    character(8) :: loc
    character(32) :: varsFound(nbVars)
    logical, external :: BiometVarHasSuffix


    numVarsFound = 0
    varsFound = ''
    cnt = 0
    do i = 1, nbVars
        label = trim(bVars(i)%label)
        if (.not. BiometVarHasSuffix(label)) then
            !> This block updates list of variables found
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
            bVars(i)%label = trim(bVars(i)%label) // '_0_0_' // trim(adjustl(loc))
        end if
    end do
end subroutine BiometAppendReplicateSuffix


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
    integer :: nsep
    integer :: i
    integer :: ix, ix2
    integer :: int_num
    integer :: io_status
    character(32) :: token, num
    integer, external :: CountCharInString


    !> Check if variable name contains 3 underscores and text between them
    !> is an integer number.
    token = trim(label)
    BiometVarHasSuffix = .true.

    nsep = CountCharInString(token, '_')
    if (nsep < 3) then
        BiometVarHasSuffix = .false.
        return
    else
        !> Extract substring eligible to be suffix
        do i = 1, nsep - 3
            token = token(index(token, '_')+1: len_trim(token))
        end do
        token = token(index(token, '_'): len_trim(token))
        !> Now we are left with something like _xx_yyy_zzzz. Check that
        !> xx, yyy and zzz are integers
        do i = 1, 3
            ix = index(token, '_')
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
        end do
    end if
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
!        call BiometInterpretPositionalQualifier(bVars(i))
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
    integer :: jj, cnt
    character(64) :: tsString
    character(32) :: item
    character(1)  :: sepa


    skip_row = .false.
    sepa = bFileMetadata%separator
    !> Check that row contains something
    if (len_trim(bFileMetadata%data_label) > 0) then
        if (index(row(:32), trim(bFileMetadata%data_label)) == 0) then
            skip_row = .true.
            return
        else
            !> Skip label if present
            row = row(len_trim(bFileMetadata%data_label)+2: len_trim(row))
        end if
    end if

    tsString = ''
    jj = 0
    cnt = 0
    do
        if (jj > nbVars + bFileMetadata%numTsCol - 1) exit
        select case(index(row, sepa))
            case(0)
                item = row(1: len_trim(row))
            case(1)
                !> Eliminate leading separator
                row = row(2:len_trim(row))
                cycle
            case default
                item = row(1: index(row, sepa)-1)
        end select

        !> If item is empty, cycle
        if (jj < nbVars + bFileMetadata%numTsCol - 1 .and. len_trim(item) == 0) cycle
        jj = jj + 1

        !> Eliminate item from row
        row = row(len_trim(item)+2: len_trim(row))

        !> If row is now empty, something went wrong (row was too short), so
        !> exit with error
        if (jj < nbVars + bFileMetadata%numTsCol .and. len_trim(row) == 0) then
            skip_row = .true.
            return
        end if

        !> From item, extract timestamp or biomet value as appropriate
        if (any(bFileMetadata%tsCols(1:bFileMetadata%numTsCol) == jj)) then
            tsString = trim(tsString) // trim(item)
        else
            cnt = cnt + 1
            read(item, *) vals(cnt)
        end if
    end do

    !> Retrieve timestamp from timestamp-string
    call FilenameToTimestamp(trim(tsString), &
        bFileMetadata%tsPattern, bFileMetadata%tsIso, tstamp)
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
! \brief       Interpret label and return var, hpos, vpos, rep
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
!***************************************************************************
subroutine BiometInterpretPositionalQualifier(bVar)
    use m_rp_global_var
    implicit none
    !> In/out variables
    type(BiometVarsType), intent(inout) :: bVar
    !> Local variables
    character(32) :: s
    character(32), external :: biometBaseName


    s = bVar%label(len_trim(bVar%base_name) + 2: len_trim(bVar%label))
    !> location
    read(s(1:index(s, '_')-1), '(i3)') bVar%hpos
    !> location
    s = s(index(s, '_') + 1 : len_trim(s))
    read(s(1:index(s, '_')-1), '(i3)') bVar%vpos
    !> rep
    s = s(index(s, '_') + 1 : len_trim(s))
    read(s(1:len_trim(s)), '(i3)') bVar%rep
end subroutine BiometInterpretPositionalQualifier

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
                            bSet(:, i) = bSet(:, i) + 273.15d0
                        end where
                    case('F','�F')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = (bSet(:, i) - 32d0) * 5d0 / 9d0 &
                                + 273.15d0
                        end where
                    case('CK')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 1d-2
                        end where
                    case('CC','C�C')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 1d-2 + 273.15d0
                        end where
                    case('CF','C�F')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = (bSet(:, i) * 1d-2 - 32d0) * 5d0 / 9d0 &
                                + 273.15d0
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
subroutine BiometAggregate(Set, nrow, ncol, Aggr)
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
            case ('INTEGRATION')
                call SumNoError(bSet(:, i), size(bSet, 1), 1, vAggr, error)
            case ('AVERAGING')
                call AverageNoError(bSet(:, i), size(bSet, 1), 1, vAggr, error)
            case ('ANGULAR_AVERAGING')
                call AngularAverageNoError(bSet(:, i), size(bSet, 1), 1, vAggr, error)
            case default
        end select
        Aggr(i) = vAggr
    end do
end subroutine BiometAggregate

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
    character(LongInstringLen) :: dataline
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
            read(udf, '(a)', iostat = io_status) dataline
            if (io_status > 0) cycle recs_loop
            if (io_status < 0) then
                FileList(nfl)%timestamp = nullTimestamp
                cycle files_loop
            end if

            call tsStringFromRec(dataline, nbItems, tsString)
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

!***************************************************************************
!
! \brief       Sniff biomet metadata to infer number of biomet vars
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
!***************************************************************************
subroutine biometSniffMetaFile(IniFile, skip_file)
    use m_rp_global_var
    implicit none
    !> In/out variables
    character(*), intent(in) :: IniFile
    logical, intent(out) :: skip_file
    !> Local variables
    integer :: io_status
    integer :: nlines
    integer :: i
    integer :: nvar
    type(text) :: Tags(MaxNLinesIni)


    skip_file = .false.
    open(udf, file = IniFile, status = 'old', iostat = io_status)
    if (io_status /= 0) then
        skip_file = .true.
        return
    else
        !> parse ini file and store all tags found in it
        call StoreIniTags(udf, '', Tags, nlines)
        close(udf)
        nbVars = 0
        do i = 1, nlines
            if (index(Tags(i)%Label, '_variable') /= 0) then
                !> extract var number from label
                read(Tags(i)%Label(8:index(Tags(i)%Label, '_', .true.)-1), *) nvar
                if (nvar > nbVars) nbVars = nvar
            end if
        end do
    end if
end subroutine biometSniffMetaFile

!***************************************************************************
!
! \brief       Initialize data structure needed for reading embedded biomet
!              metadata
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
!***************************************************************************
subroutine biometInitEmbedded()
    use m_rp_global_var
    implicit none
    !> Local variable
    integer :: i
    integer :: ii
    integer :: j
    integer :: nc
    integer :: nn
    integer, parameter :: cinit = 2
    integer, parameter :: ninit = 3
    integer, parameter :: nctags = 7
    integer, parameter :: nntags = 6
    character(32) :: indx
    character(16) :: ctags(nctags)
    data ctags(1:nctags) /'_variable', '_id', '_instrument', '_unit_in', &
        '_unit_out', '_aux1', '_aux2'/
    character(16) :: ntags(nntags)
    data ntags(1:nntags) /'_gain', '_offset', '_aux3', '_aux4', &
        '_aux5', '_aux6'/


    !> Total number of (char and num) tags to be written
    nc = nbVars * nctags + cinit
    nn = nbVars * nntags + ninit

    !> Allocate variable
    if (allocated(BiometCTags)) deallocate(BiometCTags)
    if (allocated(BiometCTagFound)) deallocate(BiometCTagFound)
    if (allocated(BiometNTags)) deallocate(BiometNTags)
    if (allocated(BiometNTagFound)) deallocate(BiometNTagFound)
    allocate(BiometCTags(nc))
    allocate(BiometCTagFound(nc))
    allocate(BiometNTags(nn))
    allocate(BiometNTagFound(nn))


    !> Write fixed values
    BiometCTags(1)%Label = 'biomet_separator'
    BiometCTags(2)%Label = 'biomet_data_label'

    BiometNTags(1)%Label = 'biomet_header_rows'
    BiometNTags(2)%Label = 'biomet_file_duration'
    BiometNTags(3)%Label = 'biomet_data_rate'

    !> Write dynamic values
    do i = 1, nbVars
        !> Character tags
        do j = 1, nctags
            ii = cinit + (i - 1) * nctags + j
            call int2char(i, indx, 0)
            BiometCTags(ii)%Label = 'biomet_' // trim(indx) // trim(ctags(j))
        end do
        !> Numeric tags
        do j = 1, nntags
            ii = ninit + (i - 1) * nntags + j
            call int2char(i, indx, 0)
            BiometNTags(ii)%Label = 'biomet_' // trim(indx) // trim(ntags(j))
        end do
    end do
end subroutine biometInitEmbedded

!***************************************************************************
!
! \brief       Infer standard FLUXNET labels from actual variable labels
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
!***************************************************************************
character(32) function biometBaseName(bLabel) result(base_name)
    use m_rp_global_var
    implicit none
    !> Local variable
    integer :: i
    integer :: n
    character(*), intent(in) :: bLabel
    !> In/out variables
    character(len(bLabel)) :: s

    integer, external :: CountCharInString

    !> Count number of underscores
    n = CountCharInString(bLabel, '_')

    s = bLabel
    base_name = ''
    if (n > 3) then
        do i = 1, n-2
            base_name = trim(base_name) // s(1:index(s, '_'))
            s = s(index(s, '_') + 1: len_trim(s))
        end do
        base_name = base_name(1:len_trim(base_name)-1)
    elseif (n == 3) then
        base_name = s(1:index(s, '_')-1)
    else
        base_name = trim(s)
    end if
end function biometBaseName

!***************************************************************************
!
! \brief       Infer standard FLUXNET labels from actual variable labels
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
!***************************************************************************
character(32) function positionalQualifier(bVar) result(qPositional)
    use m_rp_global_var
    implicit none
    !> Local variable
    type(BiometVarsType), intent(in) :: bVar

    !> In/out variables
    character(8) :: datum
    character(32) :: s

    s = '_'
    call int2char(bVar%hpos, datum, 0)

    s = trim(s) // trim(datum) // '_'
    call int2char(bVar%vpos, datum, 0)
    s = trim(s) // trim(datum) // '_'
    call int2char(bVar%rep, datum, 0)
    qPositional = trim(s) // trim(datum)

end function positionalQualifier
