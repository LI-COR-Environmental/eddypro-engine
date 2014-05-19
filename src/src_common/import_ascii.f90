!***************************************************************************
! import_ascii.f90
! ----------------
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
! \brief       Reads in generic ASCII table that does not include text columns \n
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
subroutine ImportAscii(FirstRecord, LastRecord, LocCol, fRaw, nrow, ncol, N, FileEndReached)
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
    integer :: io_status
    real(kind = sgl) :: TmpfRaw(nrow, NumCol)
    type(ColType) :: TmpCol(MaxNumCol)
    character(62) :: lab


    FileEndReached = .false.
    TmpCol = NullCol
    fRaw = error
    TmpfRaw = error
    N = 0
    tN = 0

    !> Skip header if present
    if (FileInterpreter%header_rows > 0) then
        do i = 1, FileInterpreter%header_rows
            read(unat, *)
        end do
    end if

    !> Skip all records until FirstRecord
    if (FirstRecord > 1) then
        do i = 1, FirstRecord - 1
            read(unat, *, iostat = io_status)
            if (io_status /= 0) return
        end do
    end if

    if (len_trim(FileInterpreter%data_label) == 0 &
        .or. index(FileInterpreter%data_label, 'Not set') /= 0) then
        !> Read datalines directly as real arrays
        record_loop: do
            tN = tN + 1

            !> Normal exit instruction, if current record is beyond LastRecord
            if (tN > LastRecord - FirstRecord + 1) exit record_loop

            N = N + 1
            read(unat, *, iostat = io_status) TmpfRaw(N, 1:NumCol)
            if (io_status < 0 .or. io_status == 5001 .or. io_status == 5008) then
                N = N - 1
                FileEndReached = .true.
                exit record_loop
            end if
            if (io_status == 5010) then
                N = N - 1
                cycle record_loop
            end if
            if (io_status > 0) then
                N = N - 1
                cycle record_loop
            end if
        end do record_loop
    else
        !> Read datalines directly as real arrays and discard first field
        record_loop2: do
            tN = tN + 1

            !> Normal exit instruction, if current record is beyond LastRecord
            if (tN > LastRecord - FirstRecord + 1) exit record_loop2

            N = N + 1
            read(unat, *, iostat = io_status) lab, TmpfRaw(N, 1:NumCol)
            if (io_status < 0 .or. io_status == 5001) exit record_loop2
            if (io_status == 5010) then
                !backspace(unat)
                N = N - 1
                cycle record_loop2
            end if
            if (io_status > 0) then
                N = N - 1
                cycle record_loop2
            end if
            if (trim(adjustl(lab)) /= trim(adjustl(FileInterpreter%data_label))) then
                N = N - 1
                cycle record_loop2
            end if
        end do record_loop2
    end if

    !> Store only hot columns
    jj = 0
    do j = 1, NumCol
        if (LocCol(j)%var /= 'ignore') then
            jj = jj + 1
            TmpCol(jj) = LocCol(j)
            fRaw(1:N, jj) = TmpfRaw(1:N, j)
            if (Gas4CalRefCol == j) Gas4CalRefCol = jj
        end if
    end do
    LocCol = TmpCol
    close(unat)
end subroutine ImportAscii
