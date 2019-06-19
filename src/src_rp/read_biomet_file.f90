!***************************************************************************
! read_biomet_file.f90
! --------------------
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
! \brief       Reads biomet data file and retrieve only necessary data
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ReadBiometFile(BiometFile, skip_file)
    use m_rp_global_var
    implicit none
    !> in/out variables
    character(*), intent(in) :: BiometFile
    logical, intent(out) :: skip_file
    !> local variables
    integer :: i
    integer :: cnt
    integer :: open_status
    integer :: io_status
    integer :: nrow
    integer :: ncol
    character(LongInstringLen) :: dataline
    character(1)  :: sepa
    logical :: skip_row


    sepa = bFileMetadata%separator
    skip_file = .false.

    !> Scan file to retrieve number of rows
    call scanCsvFile(BiometFile, sepa, 0, nrow, ncol, skip_file)
    if (skip_file) return

    open(udf, file = BiometFile, iostat = open_status)
    if (open_status /= 0) then
        skip_file = .true.
        close(udf)
        return
    end if

    !> Skip header
    if (bFileMetadata%nhead > 0) then
        do i = 1, bFileMetadata%nhead
            read(udf, *)
        end do
    end if
    fnbRecs = nrow - bFileMetadata%nhead

    !> Allocate variable holding dataset
    if (allocated(fbSet)) deallocate(fbSet)
    if (allocated(fbTs)) deallocate(fbTs)
    allocate(fbSet(fnbRecs, nbVars))
    allocate(fbTs(fnbRecs))

    !> Read data file
    !> Read file, sort out timestamp and dataset
    fbSet = error
    fbTs = nullTimestamp
    cnt = 0
    do i = 1, fnbRecs
        read(udf, '(a)', iostat = io_status) dataline
        !> Exit instruction
        if (io_status > 0) cycle
        if (io_status < 0) exit
        cnt = cnt + 1
        call BiometParseRow(dataline, fbTs(cnt), fbSet(cnt, :), &
            size(fbSet, 2), skip_row)
        if (skip_row) cycle
    end do
    close(udf)
end subroutine ReadBiometFile
