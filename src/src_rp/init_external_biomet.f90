!***************************************************************************
! init_external_biomet.f90
! ------------------------
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
! \brief       Read biomet files and figure out time step and number of rows
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
!***************************************************************************
subroutine InitExternalBiomet(bFileList, N)
    use m_rp_global_var
    implicit none
    !> In/out variables
    integer, intent(in) :: N
    type (FilelistType), intent(out) :: bFileList(N)
    !> Local variables
    integer :: cnt
    integer :: nRec, fnRec
    integer :: fnbItems
    integer :: lastcnt
    integer :: nfl
    integer :: io_status
    character(LongInstringLen) :: dataline
    character(LongInstringLen) :: dataline2
    character(64) :: tsString
    logical :: failed
    logical :: skip_record
    type(DateType), allocatable :: bTimestamp(:)
    integer, external :: tsInferTimestep
    integer, external :: SplitCount
    integer, external :: countsubstring
    character(len(dataline)), external :: replace
    type(BiometVarsType), allocatable :: lbVars(:)
    logical :: excluded_file(N)


    write(*, '(a)', advance = 'no') ' Interpreting biomet data..'

    !> Retrieve list of biomet files
    if (EddyProProj%biomet_data == 'ext_file') then
        bFileList(1)%path = AuxFile%biomet
        call basename(bFileList(1)%path, bFileList(1)%name, slash)
    elseif (EddyProProj%biomet_data == 'ext_dir') then
        call FileListByExt(Dir%biomet, trim(adjustl(EddyProProj%biomet_tail)), &
            .false., .false., 'none', .false., .false., &
            EddyProProj%biomet_recurse, bFileList, size(bFileList), .false., ' ')
    end if

    !> Loop to retrieve number of rows and cols, so that biomet variables
    !> can be allocated
    nRec = 0
    excluded_file = .false.
    size_loop: do nfl = 1, N

        !> Count number of items and rows in file
        call scanCsvFile(bFileList(nfl)%path, ',', 1, &
            fnRec, fnbItems, failed)

        !> If above failed, pass to next one
        if (failed) then
            excluded_file(nfl) = .true.
            write(*,*)
            write(*, '(a)') '  File: ' // trim(adjustl(bFileList(nfl)%path))
            call ExceptionHandler(2)
            cycle size_loop
        end if

        !> From second file on, if number of items in file is different from
        !> previous one, this is a problem and biomet cannot be used
        if (nfl == 1) then
            nbItems = fnbItems
        else
            if (fnbItems /= nbItems) then
                write(*,*)
                call ExceptionHandler(70)
                EddyProProj%biomet_data = 'none'
                return
            end if
        end if

        !> Update max number of biomet records (the -2 is because we know
        !> external biomet files have a 2-line header)
        nRec = nRec + fnRec - 2
    end do size_loop

    !> Control
    if (nRec < 1) then
        call ExceptionHandler(71)
        EddyProProj%biomet_data = 'none'
        return
    end if

    !> Allocate timestamp variable to store timestamps of biomet data, here only
    !> to the purpose of inferring the time-step
    allocate(bTimestamp(nRec))
    bTimestamp = nullTimestamp

    !> Loop over all biomet files
    lastcnt = 0
    files_loop: do nfl = 1, N

        !> If file was excluded above, cycle
        if (excluded_file(nfl)) cycle files_loop

        !> Open biomet file
        open(udf, file = bFileList(nfl)%path, status = 'old', &
            iostat = io_status)

        !> If above failed, pass to next one
        if (io_status /= 0) then
            call ExceptionHandler(2)
            cycle files_loop
        end if

        !> Read header, retrieve variable names and units
        read(udf, '(a)', iostat = io_status) dataline
        read(udf, '(a)', iostat = io_status) dataline2

        !> If timestamp labels are 'Date' and 'Time', replace with
        !> 'Timestamp_1' and 'Timestamp_2', e.g. Sutron case
        dataline = replace(dataline, 'Date', 'TIMESTAMP_1', len(dataline))
        dataline = replace(dataline, 'Time,', 'TIMESTAMP_2,', len(dataline))

        !> Retrieve number of biomet variables excluding
        !> TIMESTAMP-related items from dataline
        nbVars = SplitCount(dataline, bFileMetadata%separator, &
            'TIMESTAMP', .false.)

        !> Allocate and initialize bVars
        if (allocated(bVars)) deallocate(bVars)
        allocate(bVars(nbVars))
        bVars = nullbVar

        !> Allocate vars for aggregated biomet values
        if (allocated(bAggr)) deallocate(bAggr)
        allocate(bAggr(nbVars))
        if (allocated(bAggrFluxnet)) deallocate(bAggrFluxnet)
        allocate(bAggrFluxnet(nbVars))

        !> Retrieve variables and timestamp prototype from
        !> header (labels and units rows)
        call RetrieveExtBiometVars(dataline, dataline2, nbItems)
        if (EddyProProj%biomet_data == 'none') return

        !> Variables consistency among different biomet files
        if (nfl == 1) then
            allocate(lbVars(nbVars))
            lbVars = bVars
        else
            if (size(lbVars) /= size(bVars) &
                .or. (any(lbVars(:)%label /= bVars(:)%label) &
                .or. any(lbVars(:)%unit_in /= bVars(:)%unit_in))) then
                write(*,'(a)')
                call ExceptionHandler(79)
                EddyProProj%biomet_data = 'none'
                return
            end if
        end if

        !> Start loop on file rows
        cnt = lastcnt
        recs_loop: do
            read(udf, '(a)', iostat = io_status) dataline
            if (io_status > 0) cycle recs_loop
            if (io_status < 0) exit recs_loop

            !> Retrieve timestamp info from rec
            call tsStringFromRec(dataline, nbItems, tsString)
            cnt = cnt + 1

            !> Retrieve timestamp from timestamp string
            call BiometTimestamp(trim(adjustl(bFileMetadata%tsPattern)), &
                tsString, bTimestamp(cnt), skip_record)
            if (skip_record) cycle recs_loop

            end do recs_loop
        close(udf)
        lastcnt = cnt
    end do files_loop
    close(udf)

    !> Assess Timestep
    bFileMetadata%time_step = int(tsInferTimestep(bTimestamp(:nRec), nRec) / 60)
    bFileMetadata%tolerance = bFileMetadata%time_step / 2
    write(*, '(a)')
    write(*, '(a, i6)')    '  Number of variables: ', nbVars
    write(*, '(a, i6)')    '  Number of records:   ', nRec
    write(*, '(a, i6, a)') '  Inferred time-step:  ', &
        bFileMetadata%time_step, ' min'

    !> Determine nbRecs, the maximum number of biomet data available for
    !> each averaging interval
    if (mod(RPsetup%avrg_len, bFileMetadata%time_step) == 0) then
        nbRecs = RPsetup%avrg_len/bFileMetadata%time_step
    else
        nbRecs = RPsetup%avrg_len/bFileMetadata%time_step + 1
    end if

    !> Now that number of vars and of records are known, allocate biomet
    !> dataset and timestamp array
    allocate(bSet(nbRecs, nbVars))
    allocate(bTs(nbRecs))

    !> Fill variables information based on label and other available fields
    call BiometEnrichVarsDescription()

    !> No data label is allowed in external biomet files
    bFileMetadata%data_label = ''

    !> Accounts for timestamp columns
    call BiometUpdateSelectionOrder()

    !> Put biomet files in chronological order, regardless of file names
    call BiometFileListInChronologicalOrder(bFileList, N)

    !> Put biomet vars in order (by variable label and secondary by profile)
    !> NOT DONE FOR THE MOMENT
!    call BiometOrderVars()

    write(*, '(a)') ' Done.'
end subroutine InitExternalBiomet
