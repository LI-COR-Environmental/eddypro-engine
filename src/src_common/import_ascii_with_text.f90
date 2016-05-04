!***************************************************************************
! import_ascii_with_text.f90
! --------------------------
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
! \brief       Reads in generic ASCII table that includes text columns \n
!              Selects only hot columns, discarding "ignore" \n
!              and "not_numeric" columns
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ImportAsciiWithText(FirstRecord, LastRecord, LocCol, fRaw, &
    nrow, ncol, N, FileEndReached)
    use m_common_global_var
    implicit none
    !> In/out variables
    integer, intent(in) :: FirstRecord
    integer, intent(in) :: LastRecord
    integer, intent(in) :: nrow
    integer, intent(in) :: ncol
    integer, intent(out) :: N
    type(ColType), intent(inout) :: LocCol(MaxNumCol)
    real(kind = sgl), intent(out) :: fRaw(nrow, ncol)
    logical, intent(out) :: FileEndReached
    !> Local variables
    integer :: i
    integer :: j
    integer :: jj
    integer :: tN
    integer :: intsep
    integer :: io_status
    character(LongInstringLen) :: dataline
    character(DatumLen) :: datum
    character(64) :: lab
    type(ColType) :: TmpCol(MaxNumCol)


    FileEndReached = .false.

    !> Define TmpCol, containing only those with "hot" columns
    TmpCol = NullCol
    jj = 0
    do j = 1, NumCol
        if (LocCol(j)%var /= 'ignore' .and. LocCol(j)%var /= 'not_numeric') then
            jj = jj + 1
            TmpCol(jj) = LocCol(j)
            if (Gas4CalRefCol == j) Gas4CalRefCol = jj
        end if
    end do

    !> Skip header if present
    if (FileInterpreter%header_rows > 0) then
        do i = 1, FileInterpreter%header_rows
            read(unat, '(a)', iostat = io_status) dataline
        end do
    end if

    !> Skip all records until FirstRecord
    if (FirstRecord > 1) then
        do i = 1, FirstRecord - 1
            read(unat, '(a)', iostat = io_status) dataline
        end do
    end if

    !> Read datalines as text and parse them into the array fRaw
    N = 0
    tN = 0
    fRaw = error
    record_loop: do
        tN = tN + 1
        !> Read data line as a string and decide what to do if reading fails
        N = N + 1
        read(unat, '(a)', iostat = io_status) dataline
        if (io_status < 0) then
            N = N - 1
            FileEndReached = .true.
            exit record_loop
        end if
        if (io_status > 0 .or. (io_status == 0 .and. len_trim(dataline) == 0)) then
            N = N - 1
            cycle record_loop
        end if

        !> Record was imported, now if there is a record label, first reads it
        if (len_trim(FileInterpreter%data_label) /= 0 .and. &
            index(FileInterpreter%data_label, 'Not set') == 0) then
            intsep = index(dataline, FileInterpreter%separator)
            if (intsep == 0) intsep = len_trim(dataline) + 1
            if (len_trim(dataline) == 0) exit record_loop
            lab = dataline(1:intsep - 1)
            dataline = dataline(intsep + 1: len_trim(dataline))
            !> If label is different than the expected, cycle
            if (trim(adjustl(lab)) /= trim(adjustl(FileInterpreter%data_label))) then
                N = N - 1
                cycle record_loop
            end if
        end if

        !> Normal exit instruction, if current record is beyond the number of record to be imported
        if (tN > LastRecord - FirstRecord + 1) then
            N = N - 1
            exit record_loop
        end if

        !> Eliminate multiple separators from dataline, but currently only if it's a space
        if (FileInterpreter%separator == '') &
            call StripConsecutiveChar(dataline, FileInterpreter%separator)

        !> Parse variables out of the string
        jj = 0
        il: do j = 1, NumCol
            intsep = index(dataline, FileInterpreter%separator)
            if (intsep == 0) intsep = len_trim(dataline) + 1
            if (len_trim(dataline) == 0) exit
            datum = dataline(1:intsep - 1)
            dataline = dataline(intsep + 1: len_trim(dataline))
            if (LocCol(j)%var /= 'ignore' .and. LocCol(j)%var /= 'not_numeric') then
                jj = jj + 1
                if (EddyProProj%col(E2NumVar + DiagAnem) == j &
                    .and. (index(EddyProProj%master_sonic, 'wm') /= 0 &
                    .or. index(EddyProProj%master_sonic, 'hs') /= 0)) then
                    if (trim(datum) == '0A') datum = '10'
                    if (trim(datum) == '0B') datum = '11'
                end if
                read(datum, *, iostat = io_status) fRaw(N, jj)
                if (io_status /= 0) then
                    N = N - 1
                    cycle record_loop
                end if
            end if
        end do il
    end do record_loop
    LocCol = TmpCol
    close(unat)
end subroutine ImportAsciiWithText
