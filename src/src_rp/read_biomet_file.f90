!***************************************************************************
! read_biomet_file.f90
! --------------------
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
