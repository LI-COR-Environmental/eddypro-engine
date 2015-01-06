!***************************************************************************
! init_external_biomet.f90
! ------------------------
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
    character(1024) :: record
    character(1024) :: record2
    character(64) :: tsString
    logical :: failed
    logical :: init_bVars
    logical :: skip_record
    type(DateType), allocatable :: bTimestamp(:)
    integer, external :: tsInferTimestep
    integer, external :: SplitCount
    integer, external :: countsubstring
    character(len(record)), external :: replace


    write(*, '(a)', advance = 'no') ' Initializing external biomet usage..'

    !> Retrieve list of biomet files
    if (EddyProProj%biomet_data == 'ext_file') then
        bFileList(1)%path = AuxFile%biomet
        call basename(bFileList(1)%path, bFileList(1)%name, slash)
    elseif (EddyProProj%biomet_data == 'ext_dir') then
        call FileListByExt(Dir%biomet, trim(adjustl(EddyProProj%biomet_tail)), &
            .false., 'none', .false., .false., EddyProProj%biomet_recurse, &
            bFileList, size(bFileList), .false., ' ')
    end if

    !> Loop to retrieve number of rows and cols, so that biomet variables
    !> can be allocated
    nRec = 0
    size_loop: do nfl = 1, N

        !> Count number of items and rows in file
        call scanCsvFile(bFileList(nfl)%path, ',', 1, &
            fnRec, fnbItems, failed)

        !> If above failed, pass to next one
        if (failed) then
            write(*,*)
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
    init_bVars = .true.
    files_loop: do nfl = 1, N

        !> Open biomet file
        open(udf, file = bFileList(nfl)%path, status = 'old', &
            iostat = io_status)

        !> If above failed, pass to next one
        if (io_status /= 0) then
            call ExceptionHandler(2)
            cycle files_loop
        end if

        !> Only ones, read header, retrieve variable names and units
        if (init_bVars) then

            !> Read header lines
            read(udf, '(a)', iostat = io_status) record
            read(udf, '(a)', iostat = io_status) record2

            !> If timestamp labels are 'Date' and 'Time', replace with
            !> 'Timestamp_1' and 'Timestamp_2', e.g. Sutron case
            record = replace(record, 'Date', 'TIMESTAMP_1', len(record))
            record = replace(record, 'Time,', 'TIMESTAMP_2,', len(record))

            !> Retrieve number of biomet variables excluding
            !> TIMESTAMP-related items from record
            nbVars = SplitCount(record, bFileMetadata%separator, 'TIMESTAMP', .false.)

            !> Allocate and initialize bVars
            if (allocated(bVars)) deallocate(bVars)
            allocate(bVars(nbVars))
            bVars = nullbVar

            if (allocated(bAggr)) deallocate(bAggr)
            allocate(bAggr(nbVars))

            !> Retrieve variables and timestamp prototype from
            !> header (labels and units rows)
            call RetrieveExtBiometVars(record, record2, nbItems)

            init_bVars = .false.
        else
            !> Skip header for files other than the first
            read(udf, '(a)', iostat = io_status)
            read(udf, '(a)', iostat = io_status)
        end if

        !> Start loop on file rows
        cnt = lastcnt
        recs_loop: do
            read(udf, '(a)', iostat = io_status) record
            if (io_status > 0) cycle recs_loop
            if (io_status < 0) exit recs_loop

            !> Retrieve timestamp info from rec
            call tsStringFromRec(record, nbItems, tsString)
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
    write(*, '(a, i6)') '  Number of biomet variables: ', nbVars
    write(*, '(a, i6)') '  Number of biomet records:   ', nRec
    write(*, '(a, i6, a)') '  Inferred biomet time step:  ', &
        bFileMetadata%time_step, 'min'

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

    !> No data label is allowed in external biomet files
    call BiometUpdateSelectionOrder()

    !> Put biomet files in chronological order, regardless of file names
    call BiometFileListInChronologicalOrder(bFileList, N)

    !> Put biomet vars in order (by variable label and secondary by profile)
    !> NOT DONE FOR THE MOMENT
!    call BiometOrderVars()

    write(*, '(a)') ' done'
end subroutine InitExternalBiomet
