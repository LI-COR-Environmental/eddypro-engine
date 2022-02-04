!***************************************************************************
! import_current_period.f90
! -------------------------
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
! \brief       Create raw and biomet datasets for current period. Stop importing
!              data if a lag in adjacent files is found (of course only if
!              more than 1 file is needed)
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ImportCurrentPeriod(InitialTimestamp, FinalTimestamp, FileList, &
    NumFiles, FirstFile, LocBypassCol, MaxNumFileRecords, MetaIsNeeded, &
    BiometIsNeeded, logout, Raw, nrow, ncol, N, &
    bDataFound, skip_period, NextFile, LocCol, printout)

    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    integer, intent(in) :: NumFiles
    integer, intent(in) :: FirstFile
    integer, intent(in) :: MaxNumFileRecords
    logical, intent(in) :: BiometIsNeeded
    logical, intent(in) :: logout
    type(ColType), intent(inout) :: LocCol(MaxNumCol)
    type(ColType), intent(in) :: LocBypassCol(MaxNumCol)
    type(FileListType), intent(in) :: FileList(NumFiles)
    type(DateType), intent(in) :: InitialTimestamp
    type(DateType), intent(in) :: FinalTimestamp
    integer, intent(out) :: N
    integer, intent(out) :: NextFile
    real(kind = sgl), intent(out) :: Raw(nrow, ncol)
    logical, intent(in) :: printout
    logical, intent(out) :: bDataFound
    logical, intent(out) :: skip_period
    logical, intent(inout) :: MetaIsNeeded
    !> local variables
    integer :: pN
    integer :: faulty_col
    integer :: CurrentFile
    integer :: FirstRecord
    integer :: LastRecord
    real(kind = sgl), allocatable :: fRaw(:, :)
    real(kind = sgl) :: zero = 0.0
    real(kind = dbl) :: dzero = 0.0_dbl
    logical :: InitialMetaIsNeeded
    logical :: skip_file
    logical :: passed(32)
    type(DateType) :: CurrentTimestamp
    type(DateType) :: FollowingTimestamp
    logical, external :: FileIsRelevantToCurrentPeriod
    logical :: FileEndReached


    !> Initializations
    bDataFound = .false.
    skip_period = .false.
    Raw = error

    !> Timestamp of beginning of current period, as an
    !> initialization to check files contiguity
    call FilenameToTimestamp(FileList(FirstFile)%name, &
        EddyProProj%fname_template, EddyProLog%iso_format, CurrentTimestamp)

    !> Loop on all files relevant to current period
    InitialMetaIsNeeded = MetaIsNeeded
    pN = 0
    if (EddyProProj%biomet_data == 'embedded') nbRecs = 0
    Raw = 0.0_dbl
    CurrentFile = FirstFile
    rawfile_loop: do
        FileEndReached = .false.
        if (CurrentFile > NumFiles) then
            skip_period = .true.
            NextFile = CurrentFile + 1
            return
        end if
        if (CurrentFile > FirstFile) then
            !> Check if current file is relevant to current period.
            if(.not. FileIsRelevantToCurrentPeriod(FileList(CurrentFile)%name, &
                InitialTimestamp, FinalTimestamp)) then
                N = pN
                exit rawfile_loop
            end if

            !> Check file contiguity, if not contiguous, exit cycle
            call FilenameToTimestamp(FileList(CurrentFile)%name, &
                EddyProProj%fname_template, EddyProLog%iso_format, FollowingTimestamp)

            !> If timestamp refers to end of the period, need to set it back to
            !> end of period in this context
            if (EddyProLog%tstamp_end) &
                CurrentTimestamp = CurrentTimestamp + DatafileDateStep

            if (FollowingTimestamp > CurrentTimestamp + DatafileDateStep) then
                N = pN
                NextFile = CurrentFile
                if (N < 1) skip_period = .true.
                return
            end if
            CurrentTimestamp = FollowingTimestamp
            !> If fixed metadata are to be used, does so
            if (EddyProProj%use_extmd_file) then
                LocCol = LocBypassCol
            else
                LocCol = NullCol
            end if
        end if

        !> Some logging
        if (logout) then
            if (CurrentFile > FirstFile)  then
                write(*, *) '          ..' // slash, &
                trim(adjustl(FileList(CurrentFile)%name))
            else
                write(*, *) ' File(s): ..' // slash, &
                trim(adjustl(FileList(CurrentFile)%name))
            end if
        end if

        !> So far the CurrentTimestamp was treated for what it was (for
        !> checking contiguity). Now set anyway at the beginning of file period
        if (EddyProLog%tstamp_end) &
            CurrentTimestamp = CurrentTimestamp - DatafileDateStep

        !> Calculate which portion of the current file shall be imported
        call PortionOfFileInCurrentPeriod(InitialTimestamp, FinalTimestamp, &
            CurrentTimestamp, FirstRecord, LastRecord)
        if (LastRecord < 0) LastRecord = MaxNumFileRecords

        !> Time to allocate the variable holding the amount
        !> of data imported from current file
        allocate(fRaw(LastRecord - FirstRecord + 1, ncol))

        !> Actual file import
        fRaw = error
        select case(EddyProProj%ftype)
            case ('licor_ghg')
                !> for LICOR files, goes to unzipper and file reading
                call ReadLicorGhgArchive(FileList(CurrentFile)%path, &
                    FirstRecord, LastRecord, LocCol, LocBypassCol, MetaIsNeeded, &
                    BiometIsNeeded, EddyProProj%run_mode /= 'md_retrieval', &
                    EddyProProj%run_mode /= 'md_retrieval', &
                    fRaw, size(fRaw, 1), size(fRaw, 2), skip_file, passed, &
                    faulty_col, N, FileEndReached, printout)

                !> File skip control
                if (skip_file .or. (.not.passed(1))) then
                    if (skip_file) call ExceptionHandler(24)
                    if (.not.passed(1)) then
                        call InformOfMetadataProblem(passed, faulty_col)
                        call ExceptionHandler(25)
                    end if
                    N = pN
                    CurrentFile = CurrentFile + 1
                    if (allocated(fRaw)) deallocate(fRaw)
                    cycle rawfile_loop
                end if

            case default
                !> For all other file types
                call ImportNativeData(FileList(CurrentFile)%path, FirstRecord, &
                    LastRecord, LocCol, fRaw, size(fRaw, 1), size(fRaw, 2), &
                    skip_file, N, FileEndReached)

                !> File skip control
                if (skip_file) then
                    call ExceptionHandler(28)
                    N = pN
                    CurrentFile = CurrentFile + 1
                    if (allocated(fRaw)) deallocate(fRaw)
                    cycle rawfile_loop
                end if
        end select

        if (EddyProProj%run_mode /= 'md_retrieval') then
            !> Substitute NaN and Inf with error codes
            where (IsNaN(fRaw(1:N, :)) .or. &
                fRaw(1:N, :) == 1. / zero .or. &
                fRaw(1:N, :) == -1. / zero) fRaw(1:N, :) = error

            !> In case too many data are available (it happens..),
            !> limits to the max reasonable that is such that size
            !> of Raw is not exceeded.
            N = min(N, size(Raw, 1) - pN)

            !> Store raw data
            Raw(pN + 1: pN + N, :) = fRaw(1:N, :)
            pN = pN + N

            !> Biomet data handling
            if (BiometIsNeeded .and. fnbRecs > 0) then

                !> substitute NaN and Inf with error code
                where (IsNaN(fbSet(:, :)) .or. &
                    fbSet(:, :) == 1.0_dbl / dzero .or. &
                    fbSet(:, :) == -1.0_dbl / dzero) &
                    fbSet(:, :) = error

                !> Extend size of bSet to accommodate new data
                if (nbRecs == 0) then
                    if (allocated(bSet)) deallocate(bSet)
                    if (allocated(bTs)) deallocate(bTs)
                    allocate(bSet(fnbRecs, nbVars))
                    allocate(bTs(fnbRecs))
                else
                    allocate(auxbSet(size(bSet, 1)+fnbRecs, size(bSet, 2)))
                    allocate(auxbTs(size(bTs)+fnbRecs))
                    auxbSet(1:size(bSet, 1), :) = bSet(:, :)
                    auxbTs(1:size(bTs)) = bTs(:)
                    deallocate(bSet)
                    deallocate(bTs)
                    allocate(bSet(size(auxbSet, 1), size(auxbSet, 2)))
                    allocate(bTs(size(auxbTs)))
                    bSet = auxbSet
                    bTs = auxbTs
                    deallocate(auxbSet)
                    deallocate(auxbTs)
                end if

                !> Append new biomet data to bSet
                bSet(nbRecs + 1: nbRecs + fnbRecs, :) = fbSet(1:fnbRecs, :)
                bTs(nbRecs + 1: nbRecs + fnbRecs) = fbTs(1:fnbRecs)
                nbRecs = nbRecs + fnbRecs
            end if
        end if

        !> Deallocate fRaw
        if (allocated(fRaw)) deallocate(fRaw)

        !> If LastRecord is smaller than MaxNumFileRecords and the
        !> file was not read to its end, this was surely the last
        !> file to import, so exit cycle
        if (LastRecord < MaxNumFileRecords &
            .and. .not. FileEndReached) exit rawfile_loop

        !> Advance CurrentFile of one unit
        CurrentFile = CurrentFile + 1

        !> Check if current file index is available in Filelist,
        !> otherwise exit cycle
        if(CurrentFile > NumFiles) exit rawfile_loop
    end do rawfile_loop

    !> Define NextFile
    NextFile = CurrentFile
    N = pN
    MetaIsNeeded = InitialMetaIsNeeded
    if (nbRecs > 0) bDataFound = .true.

    !> Define an overall flag to determine whether there are enough
    !> raw data for the current period
    if (N < 1) skip_period = .true.

end subroutine ImportCurrentPeriod
