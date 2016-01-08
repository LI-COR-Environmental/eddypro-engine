!***************************************************************************
! import_binary.f90
! -----------------
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
! \brief       Reads in generic binary file with fixed-length records \n
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
subroutine ImportBinary(FirstRecord, LastRecord, LocCol, fRaw, nrow, ncol, N, FileEndReached)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: FirstRecord
    integer, intent(in) :: LastRecord
    integer, intent(in) :: nrow
    integer, intent(in) :: ncol
    integer, intent(out) :: N
    type(ColType), intent(inout) :: LocCol(MaxNumCol)
    real(kind = sgl), intent(out) :: fRaw(nrow, ncol)
    logical, intent(out) :: FileEndReached
    !> local variables
    integer(kind = 2) :: IntRec(NumCol)
    integer :: io_status
    integer :: i
    integer :: j
    integer :: jj
    integer :: k
    integer :: rec_num
    integer :: nrec
    integer :: head_nlines
    integer :: LastByteHeader
    integer(kind = 1) :: aux, paux
    integer(kind = 1) :: b(12)
    real(kind = sgl) :: TmpfRaw(nrow, NumCol)
    type(ColType) :: TmpCol(MaxNumCol)


    FileEndReached = .false.

    !> Detect first byte of data (skip ASCII header)
    TmpCol = NullCol
    fRaw = error
    IntRec = nint(error)
    nrec = 0
    aux = 0_1
    head_nlines = 0
    LastByteHeader = 0

    if (Binary%head_nlines > 0) then
        select case(Binary%ascii_head_eol)
            !> Windows case (CR/LF carriage control)
            case ('CR/LF')
            do
                paux = aux
                nrec = nrec + 1
                read(unat, rec = nrec, iostat = io_status) aux
                if (io_status /= 0) exit
                if (aux == 10_1 .and. paux == 13_1) then
                    head_nlines = head_nlines + 1
                    if (head_nlines == Binary%head_nlines) then
                        LastByteHeader = nrec
                        exit
                    end if
                end if
            end do
            !> Unix case (LF carriage control)
            case ('LF')
            do
                nrec = nrec + 1
                read(unat, rec = nrec, iostat = io_status) aux
                if (io_status /= 0) exit
                if (aux == 10_1) then
                    head_nlines = head_nlines + 1
                    if (head_nlines == Binary%head_nlines) then
                        LastByteHeader = nrec
                        exit
                    end if
                end if
            end do
            !> Mac case (CR carriage control)
            case ('CR')
            do
                nrec = nrec + 1
                read(unat, rec = nrec, iostat = io_status) aux
                if (io_status /= 0) exit
                if (aux == 13_1) then
                    head_nlines = head_nlines + 1
                    if (head_nlines == Binary%head_nlines) then
                        LastByteHeader = nrec
                        exit
                    end if
                end if
            end do
        end select
    end if

    !> If present, header has been identified and skipped, now read data
    rec_num = LastByteHeader
    i = 0
    N = 0
    record_loop: do
        i = i + 1
        do j = 1, NumCol
            !> For each variable, read <nbytes> bytes
            do i = 1, Binary%nbytes
                read(unat, rec = rec_num + 2 * j + i - 2, iostat = io_status) b(i)
                !> In case of binary files, any problem in reading the file
                !> causes EP to skip it until its end.
                if (io_status /= 0) then
                    FileEndReached = .true.
                    exit record_loop
                end if
            end do
            !> join bytes to create variable value
            IntRec(j) = 0_2
            if (Binary%little_endian) then
                !> little endian
                do i = 1, Binary%nbytes
                    do k = 0, 7
                        if (btest(b(i), k)) IntRec(j) = IBSET(IntRec(j), k + 8 * (i-1))
                    end do
                end do
            else
                !> big endian
                do i = 1, Binary%nbytes
                    do k = 0, 7
                        if (btest(b(i), k)) IntRec(j) = IBSET(IntRec(j), (Binary%nbytes - i) * 8 + k)
                    end do
                end do
            end if
        end do
        rec_num = rec_num + (NumCol * Binary%nbytes)

        !> Cycle until FirstRecord.
        if (i < FirstRecord) cycle record_loop

        !> Normal exit
        if (N > LastRecord - FirstRecord) exit record_loop

        N = N + 1
        TmpfRaw(N, :) = IntRec(:)
    end do record_loop
    close(udf)

    !> Store only hot columns
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
end subroutine ImportBinary
