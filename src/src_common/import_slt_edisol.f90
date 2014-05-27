!***************************************************************************
! import_slt_edisol.f90
! ---------------------
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
! \brief       Reads in binary file created by EdiSol \n
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
subroutine ImportSLTEdiSol(rec_len, FirstRecord, LastRecord, LocCol, fRaw, nrow, ncol, N, FileEndReached)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer(kind = 1), intent(in) :: rec_len
    integer, intent(in) :: FirstRecord
    integer, intent(in) :: LastRecord
    integer, intent(in) :: nrow
    integer, intent(in) :: ncol
    integer, intent(out) :: N
    type(ColType), intent(inout) :: LocCol(MaxNumCol)
    real(kind = sgl), intent(out) :: fRaw(nrow, ncol)
    logical, intent(out) :: FileEndReached
    !> local variables
    integer :: io_status
    integer :: i
    integer :: j
    integer :: jj
    integer(kind = 1) :: loc_header(20)
    integer(kind=2) :: IntRec(NumCol)
    real(kind = sgl) :: TmpfRaw(nrow, NumCol)
    type(ColType) :: TmpCol(MaxNumCol)


    FileEndReached = .false.

    !> skip header
    read(unat, rec=1) (loc_header(j), j = 1, rec_len)

    IntRec = nint(error)
    i = 0
    N = 0
    TmpfRaw = error
    record_loop: do
        i = i + 1
        read(unat, rec = i + 1, iostat = io_status) (IntRec(j), j = 1, rec_len / 2_1)
        !> In case of binary files, any problem in reading the file
        !> causes EP to skip it until its end.
        if (io_status /= 0) then
            FileEndReached = .true.
            exit record_loop
        end if

        !> Normal exit
        if (N > LastRecord - FirstRecord) exit record_loop

        !> Cycle until FirstRecord.
        if (i < FirstRecord) cycle record_loop
        N = N + 1
        TmpfRaw(N, 1:NumCol) = IntRec(1:NumCol)
    end do record_loop

    !> Store only hot columns
    fRaw = error
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
end subroutine ImportSLTEdiSol
