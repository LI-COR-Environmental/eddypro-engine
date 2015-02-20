!***************************************************************************
! biomet_retrieve_external_data.f90
! ---------------------------------
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
! \brief       Retrieve biomet data for current averaging interval \n
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine BiometRetrieveExternalData(bFileList, bnFiles, bLastFile, &
        bLastRec, tsStart, tsEnd, bDataFound, printout)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: bnFiles
    logical, intent(in) :: printout
    type(DateType), intent(in) :: tsStart
    type(DateType), intent(in) :: tsEnd
    type(FileListType), intent(in) :: bFileList(bnFiles)
    integer, intent(inout) :: bLastFile
    integer, intent(inout) :: bLastRec
    logical, intent(out) :: bDataFound
    !> Local variable
    integer :: i
    integer :: nfl
    integer :: cnt
    integer :: io_status
    real(kind=dbl) :: cSet(size(bSet, 2))
    character(1) :: sepa
    character(LongInstringLen) :: dataline
    logical :: skip_row
    logical :: isopen
    type(DateType) :: cTs
    type(DateType) :: tol


    !> Initialize biomet data to error
    bSet = error
    bAggr = error

    !> exit right away, if biomet files are finished
    if (bLastFile > bnFiles) then
        if (printout) call ExceptionHandler(72)
        return
    end if

    !> Initializations
    sepa = bFileMetadata%separator
    bDataFound = .false.
    tol = Datetype(0, 0, 0, 0, bFileMetadata%tolerance)

    !>Loop on all biomet files
    !>Start from last visited file in bFileList
    cnt = 0
    file_loop: do nfl = bLastFile, bnFiles

        !> Cycle if current bFile does not have timestamp attached to it
        if (bFileList(nfl)%timestamp == nullTimestamp) cycle file_loop

        !> Log
        if (printout) then
            if (nfl == bLastFile) write(*, '(a)') &
                '  Searching biomet data in file: '
            write(*, '(a)') '   ' &
                // trim(adjustl(bFileList(nfl)%path))
        end if

        !> Open current biomet file
        inquire(udf, opened=isopen)
        if (isopen) close(udf)
        open(udf, file = bFileList(nfl)%path, iostat = io_status)
        if (io_status /= 0) cycle file_loop

        !> Skip header if necessary
        if (bFileMetadata%nhead > 0) then
            do i = 1, bFileMetadata%nhead
                read(udf, *, iostat = io_status)
                if (io_status /= 0) cycle file_loop
            end do
        end if

        !> Skip first part of the file
        do i = 1, bLastRec
            read(udf, *, iostat = io_status)
            if (io_status /= 0) cycle file_loop
        end do

        !> Search for first biomet record within current averaging period
        rec_loop: do
            read(udf, '(a)', iostat = io_status) dataline

            !> Exit instruction
            if (io_status > 0) cycle rec_loop
            if (io_status < 0) then
                bLastFile = bLastFile + 1
                bLastRec = 0
                cycle file_loop
            end if

            !> Parse record into cSet and cTs
            call BiometParseRow(dataline, cTs, cSet, size(cSet), skip_row)
            if (skip_row) cycle rec_loop

            call BiometAdjustTimestamp(cTs)

            !> Exit instruction if current biomet timestamp if beyond latest
            !> timestamp relevant to the averaging interval
            if (cTs >= tsEnd + tol) exit file_loop
            bLastRec = bLastRec + 1

            !> if record is relevant, add it to bSet
            if (cTs >= tsStart + tol .and. cTs < tsEnd + tol) then
                bDataFound = .true.
                cnt = cnt + 1
                if (cnt > nbRecs) exit file_loop
                bSet(cnt, :) = cSet
            end if
        end do rec_loop
    end do file_loop

    if (bDataFound) then
        write(LogInteger, '(i3)') cnt
        if (printout) write(*, '(a)') '   ' // trim(adjustl(LogInteger )) &
            // ' biomet record(s) imported.'

        !> Convert data to standard units
        call BiometStandardEddyProUnits()

        !> Calculate mean values of biomet over the averaging interval
        call BiometAggretate(bSet, size(bSet, 1), size(bSet, 2), bAggr)

        !> Convert aggregated values to FLUXNET units
        call BiometStandardFluxnetUnits()
    else
        if (printout) call ExceptionHandler(72)
    end if

    !> Associate values to variables, as selected by user
    do i = bTa, bRg
        if (bSetup%sel(i) > 0) biomet%val(i) = bAggr(bSetup%sel(i))
    end do
    if (printout) write(*,'(a)') '  Done.'

end subroutine BiometRetrieveExternalData
