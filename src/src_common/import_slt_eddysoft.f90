!***************************************************************************
! import_slt_eddysoft.f90
! -----------------------
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
! \brief       Reads in binary file created by EddyMeas \n
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
subroutine ImportSLTEddySoft(FirstRecord, LastRecord, LocCol, fRaw, nrow, ncol, N, FileEndReached)
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
    ! Local variables
    integer :: io_status
    integer :: i
    integer :: j
    integer :: jj
    integer :: NumAnalog
    integer(kind = 1) :: loc_header(8 + (NumCol - 4) * 2)
    integer(kind = 2) :: IntRec(NumCol)
    real(kind = sgl) :: TmpfRaw(nrow, NumCol)
    type(ColType) :: TmpCol(MaxNumCol)
    logical :: high_res(6) = .false.


    FileEndReached = .false.

    !> Initializations
    NumAnalog = NumCol - 4

    !> Read header
    read(unat, rec = 1, iostat = io_status) (loc_header(j), j = 1, 8 + NumAnalog * 2)
    if (io_status /= 0) then
        N = 0
        return
    end if

    do j = 1, NumAnalog
        if (loc_header(7 + 2 * j) == 0_1 .or. loc_header(7 + 2 * j) / 2_1 * 2_1 == loc_header(7 + 2 * j)) then
            !> if mask byte is even, low resolution is used
            high_res(j) = .false.
        else
            !> if mask byte is odd, high resolution is used
            high_res(j) = .true.
        end if
    end do

    !> read data
    i = 0
    N = 0
    fRaw = error
    IntRec = nint(error)
    record_loop: do
        i = i + 1

        !> Normal exit
        if (N + 1 > LastRecord - FirstRecord + 1) exit record_loop

        read(unat, rec = i + 1, iostat = io_status) (IntRec(j), j = 1, ncol)
        if (io_status /= 0) exit record_loop

        !> Cycle until FirstRecord.
        if (i < FirstRecord) cycle record_loop

        !> Define TmpfRaw
        N = N + 1
        TmpfRaw(N, 1:4) = IntRec(1:4)
        do j = 1, NumAnalog
            if (high_res(j)) then
                TmpfRaw(N, 4 + j) =  (IntRec(4 + j) + 25000.) / 10.
            else
                TmpfRaw(N, 4 + j) =  IntRec(4 + j)
            end if
        end do
    end do record_loop
    N = N - 1

    !> Store only hot columns
    TmpCol = NullCol
    jj = 0
    do j = 1, NumCol
        if (LocCol(j)%var /= 'ignore' .and. LocCol(j)%var /= 'not_numeric') then
            jj = jj + 1
            TmpCol(jj) = LocCol(j)
            fRaw(1:N, jj) = TmpfRaw(1:N, j)
            if (Gas4CalRefCol == j) Gas4CalRefCol = jj
        end if
    end do
    LocCol = TmpCol

    close(unat)
end subroutine ImportSLTEddySoft
