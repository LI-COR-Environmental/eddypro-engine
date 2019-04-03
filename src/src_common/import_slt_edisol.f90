!***************************************************************************
! import_slt_edisol.f90
! ---------------------
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
