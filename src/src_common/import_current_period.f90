!***************************************************************************
! import_current_period.f90
! -------------------------
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
subroutine ImportCurrentPeriod(InitialTimestamp, FinalTimestamp, FileList, NumFiles, FirstFile, LocBypassCol, &
    MaxNumFileRecords, MetaIsNeeded, BiometIsNeeded, logout, Raw, nrow, ncol, bRaw, nbrow, nbcol, N, bN, &
    BiometDataExist, skip_period, NextFile, LocCol)

    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol, nbrow, nbcol
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
    integer, intent(out) :: bN
    integer, intent(out) :: NextFile
    real(kind = sgl), intent(out) :: Raw(nrow, ncol)
    real(kind = dbl), intent(out) :: bRaw(nbrow, nbcol)
    logical, intent(out) :: BiometDataExist
    logical, intent(out) :: skip_period
    logical, intent(inout) :: MetaIsNeeded
    !> local variables
    integer :: pN
    integer :: pbN
    integer :: faulty_col
    integer :: CurrentFile
    integer :: FirstRecord
    integer :: LastRecord
    real(kind = sgl), allocatable :: fRaw(:, :)
    real(kind = sgl) :: zero = 0.
    real(kind = dbl) :: dzero = 0d0
    logical :: InitialMetaIsNeeded
    logical :: skip_file
    logical :: passed(32)
    type(DateType) :: CurrentTimestamp
    type(DateType) :: FollowingTimestamp
    logical, external :: FileIsRelevantToCurrentPeriod
    logical :: FileEndReached


    !> Initializations
    BiometDataExist = .false.
    skip_period = .false.
    Raw = error
    bRaw = error

    !> Timestamp of beginning of current period, as an initialization to check files contiguity
    call FilenameToTimestamp(FileList(FirstFile)%name, EddyProProj%fproto, EddyProLog%iso_format, CurrentTimestamp)

    !> Loop on all files relevant to current period
    InitialMetaIsNeeded = MetaIsNeeded
    pN = 0
    pbN = 0
    Raw = 0.d0
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
            if(.not. FileIsRelevantToCurrentPeriod(FileList(CurrentFile)%name, InitialTimestamp, FinalTimestamp)) then
                N = pN
                bN = pbN
                exit rawfile_loop
            end if

            !> Check file contiguity, if not contiguous, exit cycle
            call FilenameToTimestamp(FileList(CurrentFile)%name, EddyProProj%fproto, EddyProLog%iso_format, FollowingTimestamp)
            !> To account for file names
            if (FollowingTimestamp > CurrentTimestamp + DatafileDateStep) then
                N = pN
                bN = pbN
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
                write(*, *) '          ..\', FileList(CurrentFile)%name(1:len_trim(FileList(CurrentFile)%name))
            else
                write(*, *) ' File(s): ..\', FileList(CurrentFile)%name(1:len_trim(FileList(CurrentFile)%name))
            end if
        end if

        !> So far the CurrentTimestamp was treated for what it was (for checking contiguity)
        !> Now set anyway at the beginning of the file period
        if (EddyProLog%tstamp_end) CurrentTimestamp = CurrentTimestamp - DatafileDateStep

        !> Calculate which portion of the current file shall be imported
        call PortionOfFileInCurrentPeriod(InitialTimestamp, FinalTimestamp, CurrentTimestamp, FirstRecord, LastRecord)
        if (LastRecord < 0) LastRecord = MaxNumFileRecords

        !> Time to allocate the variable holding the amount of data imported from current file
        allocate(fRaw(LastRecord - FirstRecord + 1, ncol))

        !> Actual file import
        fRaw = error
        select case(EddyProProj%ftype)
            case ('licor_ghg')
                !> for LICOR files, goes to unzipper and file reading
                call ReadLicorGhgArchive(FileList(CurrentFile)%path, FirstRecord, LastRecord, LocCol, LocBypassCol, MetaIsNeeded, &
                    BiometIsNeeded, EddyProProj%run_mode /= 'md_retrieval', EddyProProj%run_mode /= 'md_retrieval', &
                    fRaw, size(fRaw, 1), size(fRaw, 2), skip_file, passed, faulty_col, N, bN, FileEndReached)

                !> File skip control
                if (skip_file .or. (.not.passed(1))) then
                    if (skip_file) call ExceptionHandler(24)
                    if (.not.passed(1)) then
                        call InformOfMetadataProblem(passed, faulty_col)
                        call ExceptionHandler(25)
                    end if
                    N = pN
                    bN = pbN
                    CurrentFile = CurrentFile + 1
                    if (allocated(fRaw)) deallocate(fRaw)
                    cycle rawfile_loop
                end if
            case default
                !> For all other file types
                call ImportNativeData(FileList(CurrentFile)%path, FirstRecord, LastRecord, &
                    LocCol, fRaw, size(fRaw, 1), size(fRaw, 2), skip_file, N, FileEndReached)

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

            !> In case too many data are available (it happens..), limits to the max reasonable
            !> that is such that size of Raw is not exceeded.
            N = min(N, size(Raw, 1) - pN)

            !> Store raw data
            Raw(pN + 1: pN + N, :) = fRaw(1:N, :)
            pN = pN + N

            !> Biomet date handling
            if (BiometIsNeeded) then
                !> substitute NaN and Inf with error code
                where (IsNaN(BiometSet(1:bN, :)) .or. &
                    BiometSet(1:bN, :) == 1d0 / dzero .or. &
                    BiometSet(1:bN, :) == -1d0 / dzero) BiometSet(1:bN, :) = error
                !> stores all biomet data together
                bRaw(pbN + 1: pbN + bN, :) = BiometSet(1:bN, :)
                ddvec(pbN + 1: pbN + bN)   = dvec(1:bN)
                ttvec(pbN + 1: pbN + bN)   = tvec(1:bN)
                pbN = pbN + bN
            end if
        end if

        !> Deallocate fRaw
        if (allocated(fRaw)) deallocate(fRaw)

        !> If LastRecord is smaller than MaxNumFileRecords and the file was not read to its end,
        !> this was surely the last file to import, so exit cycle
        if (LastRecord < MaxNumFileRecords .and. .not. FileEndReached) exit rawfile_loop

        !> Advance CurrentFile of one unit
        CurrentFile = CurrentFile + 1

        !> Check if current file index is available in Filelist, otherwise exit cycle
        if(CurrentFile > NumFiles) exit rawfile_loop
    end do rawfile_loop

    !> Define NextFile
    NextFile = CurrentFile
    N = pN
    bN = pbN

    MetaIsNeeded = InitialMetaIsNeeded
    if (bN > 0) BiometDataExist = .true.

    !> Define an overall flag to determine whether there are enough raw data for the current period
    if (N < 1) skip_period = .true.
end subroutine ImportCurrentPeriod
