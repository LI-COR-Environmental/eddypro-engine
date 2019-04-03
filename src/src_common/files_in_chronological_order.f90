!***************************************************************************
! files_in_chronological_order.f90
! --------------------------------
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
! \brief       Rearranges file name lists in chronological order
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine FilesInChronologicalOrder(FileList, nrow, &
    StartTimestamp, EndTimestamp, indent)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(inout) :: nrow
    character(*), intent(in) :: indent
    type(FileListType), intent(inout) :: FileList(nrow)
    type(DateType), intent(out) :: StartTimestamp
    type(DateType), intent(out) :: EndTimestamp
    !> local variables
    integer :: i
    integer :: mini
    integer :: maxi
    type(FileListType) :: TmpFileList(nrow)
    type(DateType) :: tsList(nrow)
    type(DateType) :: TmpListTimestamp(nrow)
    integer :: rank(nrow)


    write(*,'(a)', advance = 'no') &
        indent // ' Arranging files in chronological order..'
    !> Initialization
    StartTimestamp = datetype(2100, 12, 31, 23, 30)
    EndTimestamp = datetype(1900, 0, 0, 0, 0)

    !> Detect start/end dates and times and store all dates
    do i = 1, nrow
        tsList(i) = FileList(i)%timestamp
        !> update start timestamp
        if (FileList(i)%timestamp < StartTimestamp) &
            StartTimestamp = FileList(i)%timestamp
        !> update end timestamp
        if (FileList(i)%timestamp > EndTimestamp) &
            EndTimestamp = FileList(i)%timestamp
    end do

    !> Sort files in a chronological sequence
    call rank_dates(tsList, rank, nrow)

    do i = 1, nrow
        TmpFileList(i) = FileList(rank(i))
        TmpListTimestamp(i) = tsList(rank(i))
    end do
    FileList = TmpFileList
    tsList = TmpListTimestamp

    !> If applicable, reduces FileList to the user-selected period
    mini = 1
    maxi = nrow
    do i = 1, nrow
        if(tsList(i) >= StartTimestamp) then
            mini = i
            exit
        end if
    end do
    do i = nrow, 1, -1
        if(tsList(i) <= EndTimestamp) then
            maxi = i
            exit
        end if
    end do

    !> Places relevant file names at the beginning of FileList
    do i = mini, maxi
        FileList(i - mini + 1) = FileList(i)
    end do
    nrow = maxi - mini + 1
    write(*,'(a)') ' Done.'
end subroutine FilesInChronologicalOrder

!***************************************************************************
!
! \brief       Ranks dates in a chronological order
! \author      Patched by Gerardo Fratini from original code from Michel Olagnon, \n
!              available in the public domain:
!              www.fortran-2000.com
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine rank_dates(date_array, rank, nrow)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    type(DateType), intent(in) :: date_array(nrow)
    integer, intent(out) :: rank(nrow)
    !> local variables
    integer, allocatable :: JWRKT(:)
    integer :: LMTNA, LMTNC
    integer :: nval, i, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB


    nval = min(size(date_array), size(rank))
    if (nval <= 0) return

    !>  Fill-in the index array, creating ordered couples
    do i = 2, nval, 2
        if (date_array(i-1) < date_array(i)) then
            rank(i-1) = i - 1
            rank(i) = i
        else
            rank(i-1) = i
            rank(i) = i - 1
        end if
    end do
    if (modulo(nval, 2) /= 0) then
        rank(nval) = nval
    end if
    allocate (JWRKT(1:nval))
    LMTNC = 2
    LMTNA = 2

    !> Iteration. Each time, the length of the ordered subsets
    !> is doubled.
    do
        If (LMTNA >= nval) exit
        IWRKF = 0
        LMTNC = 2 * LMTNC
        IWRK = 0
        !> Loop on merges of A and B into C
        do
            IINDA = IWRKF
            IWRKD = IWRKF + 1
            IWRKF = IINDA + LMTNC
            JINDA = IINDA + LMTNA
            if (IWRKF >= NVAL) then
               if (JINDA >= NVAL) exit
               IWRKF = NVAL
            end if
            IINDB = JINDA
            !> Shortcut for the case when the max of A is smaller
            !> than the min of B (no need to do anything)
            if (date_array(rank(JINDA)) <= date_array(rank(JINDA+1))) then
               IWRK = IWRKF
               cycle
            end if
            !> One steps in the C subset, that we create in the final rank array
            do
                if (IWRK >= IWRKF) then
                    !> Make a copy of the rank array for next iteration
                    rank(IWRKD:IWRKF) = JWRKT (IWRKD:IWRKF)
                    exit
                end if
                IWRK = IWRK + 1

                !> We still have unprocessed values in both A and B
                if (IINDA < JINDA) then
                    If (IINDB < IWRKF) then
                         If (date_array(rank(IINDA+1)) > date_array(rank(IINDB+1))) then
                            IINDB = IINDB + 1
                            JWRKT (IWRK) = rank(IINDB)
                        else
                            IINDA = IINDA + 1
                            JWRKT (IWRK) = rank(IINDA)
                        end If
                    else
                        !> Only A still with unprocessed values
                        IINDA = IINDA + 1
                        JWRKT (IWRK) = rank(IINDA)
                    end if
                else
                    !>  Only B still with unprocessed values
                    rank(IWRKD:IINDB) = JWRKT (IWRKD:IINDB)
                    IWRK = IWRKF
                    exit
               end if
            end do
        end do
        LMTNA = 2 * LMTNA
    end do
    !>  Clean up
    deallocate (JWRKT)
    return
end subroutine rank_dates
