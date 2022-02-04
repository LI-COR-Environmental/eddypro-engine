!***************************************************************************
! import_slt_eddysoft.f90
! -----------------------
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
    !MC integer(kind = 1) :: loc_header(8 + (NumCol - 4) * 2)
    !MC integer(kind = 2) :: IntRec(NumCol)
    integer(i1) :: loc_header(8 + (NumCol - 4) * 2)
    integer(i2) :: IntRec(NumCol)
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
        !MC if (loc_header(7 + 2 * j) == 0_1 .or. loc_header(7 + 2 * j) / 2_1 * 2_1 == loc_header(7 + 2 * j)) then
        if (loc_header(7 + 2 * j) == 0_i1 .or. loc_header(7 + 2 * j) / 2_i1 * 2_i1 == loc_header(7 + 2 * j)) then
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
    IntRec = nint(error, i2)
    record_loop: do
        i = i + 1
        read(unat, rec = i + 1, iostat = io_status) (IntRec(j), j = 1, ncol)
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

    !> Store only hot columns
    TmpCol = NullCol
    jj = 0
    do j = 1, NumCol
        if (trim(LocCol(j)%var) /= 'ignore' .and. trim(LocCol(j)%var) /= 'not_numeric') then
            jj = jj + 1
            TmpCol(jj) = LocCol(j)
            fRaw(1:N, jj) = TmpfRaw(1:N, j)
            if (Gas4CalRefCol == j) Gas4CalRefCol = jj
        end if
    end do
    LocCol = TmpCol

    close(unat)
end subroutine ImportSLTEddySoft
