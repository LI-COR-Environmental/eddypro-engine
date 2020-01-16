!***************************************************************************
! eddypro-rp_main.f90
! -------------------
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
! \brief       Program main
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
program EddyproRP
    use m_rp_global_var
    !use netcdf
    !use iso_c_binding
    !use iso_fortran_env
    implicit none

    !> Local variables
    integer :: NumberOfOkPeriods = 0
    integer :: PeriodRecords = 0
    integer :: N2 = 0
    integer :: pcount = 0
    integer :: LatestRawFileIndx = 0
    integer :: NumRawFiles = 0
    integer :: NumberOfPeriods
    integer :: i
    integer :: j
    integer :: SpecRow
    integer :: Nmax
    integer :: Nmin
    integer :: max_nsmpl
    integer :: pfn
    integer :: err_cnt1
    integer :: sec
    integer :: faulty_col
    integer :: STFlg(GHGNumVar)
    integer :: DTFlg(GHGNumVar)
    integer :: month
    integer :: day
    integer :: PeriodActualRecords
    integer :: ton
    integer :: int_doy
    integer :: bLastFile
    integer :: bLastRec
    integer :: TotNumFile
    integer :: NumFileNoRecurse
    integer :: pfStartTimestampIndx
    integer :: pfEndTimestampIndx
    integer :: toStartTimestampIndx
    integer :: toEndTimestampIndx
    integer :: rpStartTimestampIndx
    integer :: rpEndTimestampIndx
    integer :: tlagn(E2NumVar)
    integer :: MaxNumFileRecords
    integer :: NextRawFileIndx
    integer :: InitGas4CalRefCol
    integer :: nCalibEvents
    integer :: NumDynRecords
    integer :: clean
    integer :: dirty
    integer :: latestCleaning
    integer :: NumBiometFiles
    integer :: mkdir_status
    integer :: del_status

    integer, allocatable :: toH2On(:)
    integer, allocatable :: pfNumElem(:)

    real(kind = dbl) :: MissingRecords
    real(kind = dbl) :: float_doy
    real(kind = dbl) :: PP(3, 3)
    real(kind = dbl) :: Mat(3, 3)
    real(kind = dbl) :: Mat2d(2, 2)
    real(kind = dbl) :: pfVec(3)
    real(kind = dbl) :: pfVec2d(2)
    real(kind = dbl) :: PFb2d(2, MaxNumWSect) = 0.d0

    real(kind = dbl), allocatable :: bf(:)
    real(kind = sgl), allocatable :: Raw(:, :)
    real(kind = dbl), allocatable :: E2Set(:, :)
    real(kind = dbl), allocatable :: E2Primes(:, :)
    real(kind = dbl), allocatable :: UserSet(:, :)
    real(kind = dbl), allocatable :: UserPrimes(:, :)
    real(kind = dbl), allocatable :: DiagSet(:, :)
    real(kind = dbl), allocatable :: SpecSet(:, :)
    real(kind = dbl), allocatable :: pfWindBySect(:, :, :)
    real(kind = dbl), allocatable :: pfWind(:, :)

    character(10) :: loggedDate
    character(10) :: date
    character(5) :: time
    character(PathLen) :: suffixOutString
    character(64) :: TmpString1
    character(32) :: char_doy
    character(10) :: tmpDate
    character(5) ::  tmpTime
    character(128) ::  PeriodSkipMessage

    logical :: skip_period
    logical :: passed(32)
    logical :: MetaIsNeeded = .true.
    logical :: EmbBiometDataExist = .false.
    logical :: AddUserStatsHeader = .true.
    logical :: IniFileNotFound
    logical :: initialize
    logical :: initializeBiometOut
    logical :: initializeFluxnetOut
    logical :: InitializeStorage
    logical :: InitOutVarPresence
    logical :: SingMat
    logical :: make_dataset_common
    logical :: make_dataset_rp
    logical :: FilterWhat(E2NumVar)
    logical :: FileEndReached
    logical :: toInit
    logical :: BiometDataFound

    logical, allocatable :: GoPlanarFit(:)

    type (FileListType), allocatable :: RawFileList(:)
    type (FileListType), allocatable :: bFileList(:)
    type (DateType), allocatable :: RawTimeSeries(:)
    type (DateType), allocatable :: MasterTimeSeries(:)
    type (DateType) :: tsStart, tsEnd
    type (DateType) :: LastMetadataTimestamp
    type (DateType) :: tsDatasetStart
    type (DateType) :: tsDatasetEnd
    type (DateType) :: auxStartTimestamp
    type (DateType) :: auxEndTimestamp
    type (DateType) :: SelectedStartTimestamp
    type (DateType) :: SelectedEndTimestamp
    type (StatsType) :: PrevStats
    type (AmbientStateType) :: prevAmbient
    type (QCType) :: StDiff
    type (QCType) :: DtDiff
    type (ColType) :: BypassCol(MaxNumCol)
    type (TestType) :: auxTest
    type(TimeLagOptType), allocatable :: TimelagOpt(:)
    type(TimeLagDatasetType), allocatable :: toSet(:)

    integer, external :: NumOfPeriods
    integer, external :: NumberOfFilesInSubperiod
    real(kind = dbl), external :: LaggedCovarianceNoError
    real(kind = dbl) , external :: Poly6
    integer, external :: CreateDir
    include '../src_common/interfaces.inc'


    !***************************************************************************
    !***************************************************************************
    !****** INITIALIZATION PART COMMON TO ALL SW COMPONENTS ********************
    !***************************************************************************
    !***************************************************************************
    write(*, '(a)') ''
    write(*, '(a)') ' *******************'
    write(*, '(a)') '  Executing EddyPro '
    write(*, '(a)') ' *******************'
    write(*, '(a)') ''

    app = rp_app

    !> Initialize environment
    call InitEnv()

    !> By detault, create FLUXNET output
    EddyProProj%out_fluxnet = .true.

    !> Read setup file
    call ReadIniRP('RawProcess')
    allocate(bf(Meth%spec%nbins + 1))

    !> Add run-mode tag to Timestamp_FilePadding
    call TagRunMode()

    !> EddyPro Express settings
    if (EddyProProj%run_mode == 'express') call ConfigureForExpress()
    if (EddyProProj%run_mode == 'md_retrieval') call ConfigureForMdRetrieval()
    if (EddyProProj%fluxnet_mode) call ConfigureForFluxnet()

    !> Define message for skipped periods
    if (EddyProProj%run_mode /= 'md_retrieval') then
        PeriodSkipMessage = '   Flux averaging period processing time: '
    else
        PeriodSkipMessage = '  Metadata retrieving time: '
    end if

    !> Selects which datasets should be filled with error codes,
    !> based on user selection
    make_dataset_common = EddyProProj%make_dataset
    make_dataset_rp     = EddyProProj%make_dataset

    !> Selects which files to output, considering the selected
    !> spectral correction method
    if (EddyProProj%out_avrg_cosp &
        .or. EddyProProj%out_avrg_spec &
        .or. (EddyProProj%hf_meth /= 'none' &
        .and. EddyProProj%hf_meth /= 'moncrieff_97' &
        .and. EddyProProj%hf_meth /= 'massman_00')) then
        !> in this cases, passage is needed to FCC, so:
        !> don't output files, don't create dataset
        !> don't output metadata
        !> don't calculate spectral correction
        !> don't calculate fluxes 2/3
        !> don't calculate footprint
        EddyProProj%fcc_follows     = .true.
        EddyProProj%out_full        = .false.
        EddyProProj%out_md          = .false.
        make_dataset_common         = .false.
    else
        !> in this cases, does what selected by user
        EddyProProj%fcc_follows  = .false.
    end if

    !> If running in embedded mode, override some settings
    if (EddyProProj%run_env == 'embedded') call ConfigureForEmbedded()
    if (EddyProProj%run_env == 'embedded') RPsetup%out_st = .false.

    !> Create output directory if it does not exist, otherwise is silent
    mkdir_status = CreateDir('"' //trim(adjustl(Dir%main_out)) // '"')

    !> Check on filename template
    call tsValidateTemplate(EddyProProj%fname_template)

    !> Detect number of raw files and allocate RawFileList
    call NumberOfFilesInDir(Dir%main_in, '.'//EddyProProj%fext, .true., &
        EddyProProj%fname_template, TotNumFile, NumFileNoRecurse)

    if (RPsetup%recurse) then
        NumRawFiles = TotNumFile
    else
        NumRawFiles = NumFileNoRecurse
    end if
    allocate(RawFileList(NumRawFiles))

     !> Store names of data files in RawFileList
    call FileListByExt(Dir%main_in, '.'//EddyProProj%fext, .true., .true., &
        EddyProProj%fname_template, EddyProLog%iso_format, .true., &
        RPsetup%recurse, RawFileList, size(RawFileList), .true., indent0)

    if (EddyProProj%use_extmd_file) then
        !> If requested, read external metadata file \n
        !> This is the case with non-GHG files or with GHG files if user \n
        !> explicitly selects an alternative metadata file
        write(*,'(a)', advance = 'no') ' Reading alternative metadata file: "' &
            // AuxFile%metadata(1:len_trim(AuxFile%metadata)) // '"..'
        call ReadMetadataFile(Col, AuxFile%metadata, IniFileNotFound, .true.)
        if (IniFileNotFound) then
            write(*, *)
            call ExceptionHandler(22)
        end if
        !> Retrieve variables to be used (from EddyPro project file) \n
        !> and define user-variables
        call DefineUsedVariables(Col)
        MetaIsNeeded = .false.
        call MetadataFileValidation(Col, passed, faulty_col)
        if (.not. passed(1)) then
            write(*, *)
            call InformOfMetadataProblem(passed, faulty_col)
            call ExceptionHandler(23)
        end if
        write(*,'(a)') ' Done.'
    else
        !> In case of standard GHG processing, without alternative metadata \n
        !> file one GHG file must be opened to read the metadata content for \n
        !> importing information that is necessary before looping on all \n
        !> GHG files.
        allocate(Raw(1, 1))
        BypassCol = NullCol
        if (EddyProProj%ftype == 'licor_ghg') then
            i = 1
            do while (i <= NumRawFiles)
                call ReadLicorGhgArchive(RawFileList(i)%path, -1, -1, Col, &
                    BypassCol, .true., .false., .false., &
                    EddyProProj%run_mode /= 'md_retrieval', &
                    Raw, size(Raw, 1), size(Raw, 2), skip_period, passed, &
                    faulty_col, PeriodRecords, FileEndReached, .false.)
                if (.not. skip_period .and. passed(1)) exit
                i = i + 1
                call InformOfMetadataProblem(passed, faulty_col)
            end do
            if (skip_period .or. (.not. passed(1))) call ExceptionHandler(32)
        end if
        deallocate(Raw)
    end if

    !> MasterSonic-related settings
    call DetectMasterSonic(Col, NumCol)

    !> Override/adjust settings related to the MasterSonic
    call OverrideMasterSonicRelatedSettings()

    !> Now that metadata are read, can set avrg_len in case user didn't
    if (RPsetup%avrg_len <= 0) RPsetup%avrg_len = nint(Metadata%file_length)
    if (EddyProProj%run_mode == 'md_retrieval') &
        RPsetup%avrg_len = nint(Metadata%file_length)

    !> Adjust time constant for planar fit if needed
    if (Meth%det == 'ld') then
        !> If time constant is larger than flux averaging interval,
        !> limit time constant to flux averaging interval and notify
        if (RPsetup%Tconst > RPsetup%avrg_len) then
            call ExceptionHandler(91)
            RPsetup%Tconst = nint(RPsetup%avrg_len * 6d1)
        end if
        !> Default to avrg_len anyway
        if (RPsetup%Tconst <= 0) RPsetup%Tconst = nint(RPsetup%avrg_len * 6d1)
    end if

    !> Some convenient variables
    DatafileDateStep = DateType(0, 0, 0, 0, nint(Metadata%file_length))
    DateStep         = DateType(0, 0, 0, 0, RPsetup%avrg_len)
    MaxNumFileRecords   = nint(Metadata%file_length * 60d0 * Metadata%ac_freq)
    MaxPeriodNumRecords = nint(RPsetup%avrg_len     * 60d0 * Metadata%ac_freq)

    !> Remember bypass columns (or columns detected
    !> from reading a sample GHG file)
    BypassCol = Col

    !> Initialize external biomet data
    if (index(EddyProProj%biomet_data, 'ext_') /= 0) then
        if (EddyProProj%biomet_data == 'ext_dir') then
            write(*,'(a)') ' Reading external biomet file(s) from:'
            write(*,'(a)') '  ' // trim(adjustl(Dir%biomet))
            call NumberOfFilesInDir(Dir%biomet, &
                trim(adjustl(EddyProProj%biomet_tail)), &
                .false., 'none', TotNumFile, NumFileNoRecurse)
            if (EddyProProj%biomet_recurse) then
                NumBiometFiles = TotNumFile
            else
                NumBiometFiles = NumFileNoRecurse
            end if
            write(*, '(a)') ' Done.'
        else
            NumBiometFiles = 1
        end if
        if (.not. allocated(bFileList)) allocate(bFileList(NumBiometFiles))
        call InitExternalBiomet(bFileList, size(bFileList))
    else
        allocate(bFileList(1))
    end if

    !> Open biomet output file
    if (index(EddyProProj%biomet_data, 'ext_') /= 0 .and. nbVars > 0) &
        call InitBiometOut()

    !> Initialize dynamic metadata by reading the file
    !> and figuring out available variables
    if (EddyProProj%use_dynmd_file) call InitDynamicMetadata(NumDynRecords)

    !> Determine potential radiation, based on lat/long info from metadata file
    PotRad = PotentialRadiation(Metadata%lat)

    !> Initialize output files for "user" variables (non-sensitive variables)
    !> if at least one such variable exists
    if (NumUserVar > 0) call InitUserOutFiles()

    !> Retrieve timestamp array in chronological order and
    !> order RawFileList, also in chronological order
    call FilesInChronologicalOrder(RawFileList, size(RawFileList), &
        tsDatasetStart, tsDatasetEnd, '')

    !> Adjust Start/End timestamps to define the boundaries of the
    !> RawTimeSeries. Retrieve the beginning time of first file
    !> and end time of last file.
    if (EddyProLog%tstamp_end) then
        tsDatasetStart = tsDatasetStart - DatafileDateStep
    else
        tsDatasetEnd = tsDatasetEnd + DatafileDateStep
    end if
    call tsRoundToMinute(tsDatasetStart, RPsetup%avrg_len, 'earlier')

    !> Retrieve NumberOfPeriods and allocate RawTimeSeries
    NumberOfPeriods = NumOfPeriods(tsDatasetStart, tsDatasetEnd, DateStep)
    allocate(RawTimeSeries(NumberOfPeriods + 1))

    !> Create timestamp array for full dataset
    call CreateTimeSeries(tsDatasetStart, tsDatasetEnd, DateStep, &
        RawTimeSeries, size(RawTimeSeries), .true.)

    !> Check the dynamic metadata file for calibration data.
    !> If found, builds up time series of absorptance drifts
    if (DriftCorr%method /= 'none') then
        allocate(tsDrifts(NumberOfPeriods + 1))
        allocate(Calib(0:NumDynRecords))  !< elem. 0 is to alloc. start of period
        allocate(tmpCalib(0:NumDynRecords))
        if (EddyProProj%use_dynmd_file) &
            call driftRetrieveCalibrationEvents(nCalibEvents)
    end if

    !> Define exp-binned frequencies extending \n
    !> from f_min = 1/(Flux avrg length) Hz
    !> to f_max = AcFreq/2 (= Nyquist frequency) Hz
    !> It is defined a priori, constant for each data period regardless
    !> of the actual length of the averaging periods
    if (EddyProProj%run_mode /=  'md_retrieval') &
        call BinnedFrequencyVector(bf, Meth%spec%nbins, &
            RPsetup%avrg_len, Metadata%ac_freq)

    !> Allocate array containing all potential data:
    !> rows: all rows potentially needed for current period
    !> columns: all except ignored ones and flag columns
    if (.not. allocated(Raw)) allocate(Raw(MaxPeriodNumRecords, NumAllVar))

    !***************************************************************************
    !***************************************************************************
    !******* TIME LAG OPTIMIZATION IF REQUESTED ********************************
    !***************************************************************************
    !***************************************************************************

    if (trim(adjustl(Meth%tlag)) == 'tlag_opt') then
        if (.not. RPsetup%to_onthefly) then
            call ReadTimelagOptFile(TOSetup%h2o_nclass)
            if (TOSetup%h2o_nclass > 1) &
                TOSetup%h2o_class_size = floor(100d0 / TOSetup%h2o_nclass)
        else
            write(*,'(a)') ' Performing time-lag optimization:'

            if (TOSetup%subperiod) then
                !> Timestamps of start and end of time-lag optimization period
                call DateTimeToDateType(TOSetup%start_date, TOSetup%start_time, auxStartTimestamp)
                call DateTimeToDateType(TOSetup%end_date, TOSetup%end_time, auxEndTimestamp)

                !> In RawTimeSeries, detect indices of first and last files
                !> relevant to time-lag optimization
                call tsExtractSubperiodIndexes(RawTimeSeries, &
                    size(RawTimeSeries), auxStartTimestamp, auxEndTimestamp, &
                    toStartTimestampIndx, toEndTimestampIndx)
                toEndTimestampIndx = toEndTimestampIndx + 1

                if (toStartTimestampIndx == nint(error) &
                    .or. toEndTimestampIndx == nint(error)) &
                    call ExceptionHandler(49)
            else
                toStartTimestampIndx = 1
                toEndTimestampIndx = size(RawTimeSeries)
            end if

            !> Count maximum number of periods for timelag optimization
            write(TmpString1, '(i7)') toEndTimestampIndx - toStartTimestampIndx
            write(*, '(a)') '  Maximum number of flux averaging periods &
                &available for time-lag optimization: ' &
                // trim(adjustl(TmpString1))

            !> Allocate variables that depend upon maximum number of periods
            allocate(TimelagOpt(toEndTimestampIndx - toStartTimestampIndx))

            !> Loop on selected files and calculate relevant statistics
            ton = 0
            month = 0
            day   = 0
            pcount = toStartTimestampIndx - 1
            LatestRawFileIndx = 1
            bLastFile = 1
            bLastRec = 0
            DynamicMetadata = ErrDynamicMetadata
            LastMetadataTimestamp = DateType(0, 0, 0, 0, 0)
            toInit = .true.
            to_periods_loop: do
                pcount = pcount + 1

                !> If embedded metadata are to be used,
                !> reinitialize column information to null
                if (EddyProProj%use_extmd_file) then
                    Col = BypassCol
                else
                    Col = NullCol
                end if

                !> Normal exit instruction: either the last period was
                !> dealt with, or raw files are finished
                if (LatestRawFileIndx > NumRawFiles &
                    .or. pcount >= toEndTimestampIndx) exit to_periods_loop

                !> Define initial/final timestamps
                !> of current period, say [8:00 - 8:30)
                tsStart = RawTimeSeries(pcount)
                tsEnd   = RawTimeSeries(pcount + 1)

                !> Search file containing data starting from the
                !> time closest to tsStart. Searches only from most current
                !> file onward, to avoid wasting time
                call FirstFileOfCurrentPeriod(tsStart, tsEnd, RawFileList, &
                    NumRawFiles, LatestRawFileIndx, NextRawFileIndx, skip_period)

                !> Averaging period advancement
                 if (day /= 0) then
                    if (EddyProProj%caller == 'console') then
                        write(*, '(a)', advance = 'no') '#'
                    else
                        call DisplayProgress('avrg_interval', &
                            '   another small step to the time-lag: ', &
                            tsStart, 'yes')
                    end if
                end if

                !> Daily advancement
                if (day /= tsStart%day &
                    .or. month /= tsStart%month) then
                    month = tsStart%month
                    day   = tsStart%day
                    if (EddyProProj%caller == 'console') then
                        write(*, '(a)')
                        call DisplayProgress('daily','  Importing data for ', &
                            tsStart, 'no')
                    else
                        call DisplayProgress('daily','  Importing data for ', &
                            tsStart, 'yes')
                    end if
                end if

                if (skip_period) cycle to_periods_loop

                !> Import dataset for current period. If using embedded biomet,
                !> also read biomet data. On entrance, NextRawFileIndx contains
                !> the index of the file to start the current period with
                !> On exit, LatestRawFileIndx contains index of latest file used
                call ImportCurrentPeriod(tsStart, tsEnd, &
                    RawFileList, NumRawFiles, NextRawFileIndx, BypassCol, &
                    MaxNumFileRecords, MetaIsNeeded, &
                    EddyProProj%biomet_data == 'embedded', .false., &
                    Raw, size(Raw, 1), size(Raw, 2), PeriodRecords, &
                    EmbBiometDataExist, skip_period, LatestRawFileIndx, Col, &
                    .false.)

                if (skip_period) cycle to_periods_loop

                !> Period skip control with message
                MissingRecords = dfloat(MaxPeriodNumRecords - PeriodRecords) &
                    / dfloat(MaxPeriodNumRecords) * 100d0
                if (PeriodRecords > 0 .and. MissingRecords > RPsetup%max_lack) &
                    cycle to_periods_loop

                !> Filter raw data for user-defined flags
                if (RPsetup%filter_by_raw_flags) &
                    call FilterDatasetForFlags(Col, Raw, &
                        size(Raw, 1), size(Raw, 2))

                !***************************************************************
                !**** RAW FILE IMPORT FINISHES HERE ****************************
                !**** NOW STARTS DATASET DEFINITION ****************************
                !***************************************************************

                !> Allocate arrays for actual data processing
                if (.not. allocated(E2Set)) &
                    allocate(E2Set(PeriodRecords, E2NumVar))
                if (.not. allocated(E2Primes)) &
                    allocate(E2Primes(PeriodRecords, E2NumVar))
                if (.not. allocated(DiagSet)) &
                    allocate(DiagSet(PeriodRecords, MaxNumDiag))

                !> Define EddyPro set of variables for the following processing
                call DefineE2Set(Col, Raw,   size(Raw, 1),     Size(Raw, 2), &
                                    E2Set,   size(E2Set, 1),   Size(E2Set, 2), &
                                    DiagSet, size(DiagSet, 1), Size(DiagSet, 2))

                !> If H2O instrument path type is 'open', doesn't make sense
                !> to use RH classes so set it to 1.
                if (toInit .and. E2Col(h2o)%instr%path_type == 'open') then
                    TOSetup%h2o_nclass = 1
                    toInit = .false.
                end if

                !> Clean up E2Set, eliminating values that are clearly unphysical
                call CleanUpE2Set(E2Set, size(E2Set, 1), size(E2Set, 2))

                !> Define as not present, variables for which
                !> too many values are out-ranged
                call EliminateCorruptedVariables(E2Set, size(E2Set, 1), &
                    size(E2Set, 2), skip_period, .false.)

                !> If either u, v or w have been eliminated, stops processing this period
                if (skip_period) then
                    if(allocated(E2Set)) deallocate(E2Set)
                    if(allocated(E2Primes)) deallocate(E2Primes)
                    if(allocated(DiagSet)) deallocate(DiagSet)
                    cycle to_periods_loop
                end if

                if (.not. any(E2Col(co2:gas4)%present)) then
                    if(allocated(E2Set)) deallocate(E2Set)
                    if(allocated(E2Primes)) deallocate(E2Primes)
                    if(allocated(DiagSet)) deallocate(DiagSet)
                    cycle to_periods_loop
                end if

                !> Update metadata if dynamic metadata are to be used
                if (EddyProProj%use_dynmd_file) &
                    call RetrieveDynamicMetadata(tsEnd, E2Col, size(E2Col))

                !> Retrieve biomet data if they exist
                if (index(EddyProProj%biomet_data, 'ext_') /= 0) then
                    call BiometRetrieveExternalData(bFileList, size(bFileList), &
                        bLastFile, bLastRec, tsStart, &
                        tsEnd, BiometDataFound, .false.)
                elseif (EddyProProj%biomet_data == 'embedded') then
                    call BiometRetrieveEmbeddedData(EmbBiometDataExist, .false.)
                end if

                !> Calculate relative separations between
                !> the analyzers and the anemometer used
                call DefineRelativeSeparations()

                !> Override users choices if needed
                call OverrideSettings()

                !***************************************************************
                !**** DATASET DEFINITION FINISHES HERE *************************
                !**** NOW STARTS RAW DATA REDUCTION ****************************
                !***************************************************************
                !> Interpret diagnostics and filter accordingly
                if (NumDiag > 0) then
                    call InterpretLicorDiagnostics(DiagSet, &
                        size(DiagSet, 1), size(DiagSet, 2))
                    call FilterDatasetForDiagnostics(E2Set, &
                        size(E2Set, 1), size(E2Set, 2), &
                        DiagSet, size(DiagSet, 1), size(DiagSet, 2), &
                        DiagAnemometer, .true.)
                end if
                if(allocated(DiagSet)) deallocate(DiagSet)

                !> Adjust coordinate systems if the case
                call AdjustSonicCoordinates(E2Set, size(E2Set, 1), size(E2Set, 2))

                !> Filter for wind direction if requested
                if (RPSetup%apply_wdf) &
                    call FilterDatasetForWindDirection(E2Set, size(E2Set, 1), size(E2Set, 2))

                !> Calculate basic stats
                call BasicStats(E2Set, size(E2Set, 1), size(E2Set, 2), 1, .false.)
                Stats1 = Stats

                !> Calculate raw screening flags and despike data if requeste
                auxTest = TestType(.true., .false., .false., .true., .false., &
                    .false., .false., .false., .false.)
                call StatisticalScreening(E2Set, &
                    size(E2Set, 1), size(E2Set, 2), auxTest, .false.)

                !> Define as not present, variables for which
                !> too many values are out-ranged
                call EliminateCorruptedVariables(E2Set, &
                    size(E2Set, 1), size(E2Set, 2), skip_period, .false.)

                !> If either u, v or w have been eliminated,
                !> stops processing this period
                if (skip_period) then
                    if(allocated(E2Set)) deallocate(E2Set)
                    if(allocated(E2Primes)) deallocate(E2Primes)
                    cycle to_periods_loop
                end if

                !> Calculate basic stats
                call BasicStats(E2Set, &
                    size(E2Set, 1), size(E2Set, 2), 2, .false.)
                Stats2 = Stats

                !> Apply raw-level cross wind correction
                !> (after Liu et al. 2001), if requested
                if (RPsetup%calib_cw) &
                    call CrossWindCorr(E2Col(u), E2Set, &
                        size(E2Set, 1), size(E2Set, 2), .false.)

                !> Calculate basic stats and output them as requested
                call BasicStats(E2Set, &
                    size(E2Set, 1), size(E2Set, 2), 3, .false.)
                Stats3 = Stats

                !> Angle-of-attack calibration
                call AoaCalibration(E2Set, size(E2Set, 1), size(E2Set, 2))

                !> Gill WindMaster w-boost
                if (RPsetup%calib_wboost) &
                    call ApplyGillWmWBoost(E2Set, size(E2Set, 1), size(E2Set, 2))

                !> Calculate basic stats
                call BasicStats(E2Set, &
                    size(E2Set, 1), size(E2Set, 2), 4, .false.)
                Stats4 = Stats

                !> Apply rotations for tilt correction, if requested
                call TiltCorrection('double_rotation', .false., E2Set, &
                    size(E2Set, 1), size(E2Set, 2), 1, Essentials%yaw, &
                    Essentials%pitch, Essentials%roll, .false.)

                !> Calculate basic stats
                call BasicStats(E2Set, &
                    size(E2Set, 1), size(E2Set, 2), 5, .false.)
                Stats5 = Stats

                !> Convert to mixing ratios (if requested and if the case)
                if (EddyProProj%wpl) &
                    call PointByPointToMixingRatio(E2Set, &
                        size(E2Set, 1), size(E2Set, 2), .false.)

                !> Initialize analytic spectral corrections,
                !> retrieving sensor parameters
                call RetrieveSensorParams()

                !> Adjust min/max time-lags associated to columns, to fit
                !> user settings in the Time lag optimizer dialog
                call AdjustTimelagOptSettings()

                !> Calculate and compensate time-lags
                call TimeLagHandle('maxcov', E2Set, &
                    size(E2Set, 1), size(E2Set, 2), Essentials%actual_timelag, &
                    Essentials%used_timelag, Essentials%def_tlag, .true.)

                !> Calculate basic stats
                call BasicStats(E2Set, &
                    size(E2Set, 1), size(E2Set, 2), 6, .false.)
                Stats6 = Stats

                !> Calculate air and cell parameters
                call AirAndCellParameters()

                !> Apply filter for absolute limits test, if the case
                FilterWhat = .false.
                FilterWhat(co2:gas4) = .true.
                call FilterDatasetForPhysicalThresholds(E2Set, &
                    size(E2Set, 1), size(E2Set, 2), FilterWhat)

                !> Define as not present, variables for which
                !> too many values are out-ranged
                call EliminateCorruptedVariables(E2Set, &
                    size(E2Set, 1), size(E2Set, 2), skip_period, .false.)

                !> If either u, v or w have been eliminated,
                !> stops processing this period
                if (skip_period) then
                    if(allocated(E2Set)) deallocate(E2Set)
                    if(allocated(E2Primes)) deallocate(E2Primes)
                    cycle to_periods_loop
                end if

                call Fluctuations(E2Set, E2Primes, size(E2Set, 1), &
                    size(E2Set, 2), RPsetup%Tconst, Stats6, E2Col)
                if (allocated(E2Set)) deallocate(E2Set)

                !> Calculate basic stats
                call BasicStats(E2Primes, &
                    size(E2Primes, 1), size(E2Primes, 2), 7, .false.)
                if (allocated(E2Primes)) deallocate(E2Primes)
                Stats7 = Stats

                !***************************************************************
                !**** RAW DATA REDUCTION FINISHES HERE. ************************
                !**** TENTATIVE FLUX CALCULATION        ************************
                !***************************************************************

                !> Average mole fractions in [umol mol_a-1] and [mmol mol_a-1]
                call MoleFractionsAndMixingRatios()

                !> Calculate parameters for flux computation
                call FluxParams(.false.)

                !> Calculate fluxes at Level 0
                call Fluxes0_rp(.false.)

                !> Store values if all conditions are met
                ton = ton + 1
                call AddToTimelagOptDataset(TimelagOpt, size(TimelagOpt),ton)

            end do to_periods_loop
            write(*, '(a)')
            write(*, '(a)') ' Done.'

            !*******************************************************************
            !**** RAW DATA REDUCTION FINISHES HERE.    *************************
            !**** NOW STARTS TIME LAG OPT CALCULATIONS *************************
            !*******************************************************************

            !> Adjust time-lag opt dataset to eliminate errors,
            !> so that it's easier to treat them later
            allocate (toSet(ton))
            call FixTimelagOptDataset(TimelagOpt, size(TimelagOpt), &
                toSet, size(toSet), tlagn, size(tlagn))
            if (allocated(TimelagOpt)) deallocate(TimelagOpt)

            allocate(toH2On(TOSetup%h2o_nclass))

            !> Optimize time-lags                                        ******* Improve readability of this subroutine interface
            call OptimizeTimelags(toSet, size(toSet), tlagn, E2NumVar, toH2On, & 
                TOSetup%h2o_nclass, TOSetup%h2o_class_size)

            !> Write time-lag optimization results on output file
            if (.not. (Meth%tlag == 'maxcov')) &
                call WriteOutTimelagOptimization(tlagn, E2NumVar, &
                    toH2On, TOSetup%h2o_nclass, TOSetup%h2o_class_size)

            if (allocated(toH2On)) deallocate(toH2On)
            write(*,'(a)') ' Time-lag optimization session terminated.'
            write(*,'(a)')
        end if
    end if

    !***************************************************************************
    !***************************************************************************
    !********************** PLANAR FIT IF REQUESTED ****************************
    !***************************************************************************
    !***************************************************************************
    if (index(Meth%rot(1:len_trim(Meth%rot)), 'planar_fit') /= 0) then
        if (.not. RPsetup%pf_onthefly) then
            call ReadPlanarFitFile()
            if (.not. allocated(GoPlanarFit)) &
                allocate(GoPlanarFit(PFSetup%num_sec))
            GoPlanarFit = .true.
            secloop2: do sec = 1, PFSetup%num_sec
                do i = 1, 3
                    do j = 1, 3
                        if (PFMat(i, j, sec) == error) then
                            GoPlanarFit(sec) = .false.
                            cycle secloop2
                        end if
                    end do
                end do
            end do secloop2
        else
            write(*,'(a)') ' Performing planar-fit assessment:'

            !> If zero sectors were selected, set to 1 sector by
            !> default and inform
            if (PFSetup%num_sec == 0) then
                call ExceptionHandler(38)
                PFSetup%num_sec = 1
            end if

            !> Allocate variables depending upon number of sectors
            if (.not. allocated(pfNumElem))  &
                allocate(pfNumElem(PFSetup%num_sec))

            if (PFSetup%subperiod) then
                !> Timestamps of start and end of planar fit period
                call DateTimeToDateType(PFSetup%start_date, PFSetup%start_time, &
                    auxStartTimestamp)
                call DateTimeToDateType(PFSetup%end_date, PFSetup%end_time, &
                    auxEndTimestamp)

                !> In RawTimeSeries, detect indexes of first and last files
                !> relevant to planar fit
                call tsExtractSubperiodIndexes(RawTimeSeries, &
                    size(RawTimeSeries), auxStartTimestamp, auxEndTimestamp, &
                    pfStartTimestampIndx, pfEndTimestampIndx)
                pfEndTimestampIndx = pfEndTimestampIndx + 1

                if (pfStartTimestampIndx == nint(error) &
                    .or. pfEndTimestampIndx == nint(error)) &
                    call ExceptionHandler(48)
            else
                pfStartTimestampIndx = 1
                pfEndTimestampIndx = size(RawTimeSeries)
            end if

            !> Count maximum number of periods for planar fit
            write(TmpString1, '(i7)') &
                pfEndTimestampIndx - pfStartTimestampIndx
            write(*, '(a)') '  Maximum number of &
                &flux averaging periods available for planar-fit: ' &
                // trim(adjustl(TmpString1))

            !> Allocate variables that depend upon maximum number of
            !> periods for planar fit
            allocate(pfWind(pfEndTimestampIndx - pfStartTimestampIndx, 3))

            !> Loop on selected files and calculate relevant statistics
            pcount = pfStartTimestampIndx - 1
            pfn = 0
            LatestRawFileIndx = 1
            month = 0
            day   = 0
            DynamicMetadata = ErrDynamicMetadata
            LastMetadataTimestamp = DateType(0, 0, 0, 0, 0)
            pf_periods_loop: do
                pcount = pcount + 1

                !> If embedded metadata are to be used,
                !> reinitialize column information to null
                if (EddyProProj%use_extmd_file) then
                    Col = BypassCol
                else
                    Col = NullCol

                end if

                !> Normal exit instruction: either the last period
                !> was dealt with, or raw files are finished
                if (LatestRawFileIndx > NumRawFiles &
                    .or. pcount >= pfEndTimestampIndx) exit pf_periods_loop

                !> Define initial/final timestamps of
                !> current period, say [8:00 - 8:30)
                tsStart = RawTimeSeries(pcount)
                tsEnd   = RawTimeSeries(pcount + 1)

                !> Search file containing data starting from the time closest to
                !> tsStart. Searches only from most current file
                !> onward, to avoid wasting time
                call FirstFileOfCurrentPeriod(tsStart, tsEnd, &
                    RawFileList, NumRawFiles, LatestRawFileIndx, &
                    NextRawFileIndx, skip_period)

                !> Averaging period advancement
                if (day /= 0) then
                    if (EddyProProj%caller == 'console') then
                        write(*, '(a)', advance = 'no') '#'
                    else
                        call DisplayProgress('avrg_interval', &
                            '   another small step to the planar-fit: ', &
                                tsStart, 'yes')
                    end if
                end if

                !> Daily advancement
                if (day /= tsStart%day &
                    .or. month /= tsStart%month) then
                    month = tsStart%month
                    day   = tsStart%day
                    if (EddyProProj%caller == 'console') then
                        write(*, '(a)')
                        call DisplayProgress('daily', &
                            '  Importing wind data for ', tsStart, 'no')
                    else
                        call DisplayProgress('daily', &
                            '  Importing wind data for ', tsStart, 'yes')
                    end if
                end if

                if (skip_period) cycle pf_periods_loop

                !> Import dataset for current period. If using embedded biomet,
                !> also read biomet data
                !> On entrance, NextRawFileIndx contains the index of the file
                !> to start the current period with.
                !> On exit, LatestRawFileIndx contains the index of
                !> the latest file used
                call ImportCurrentPeriod(tsStart, tsEnd, &
                    RawFileList, NumRawFiles, NextRawFileIndx, BypassCol,  &
                    MaxNumFileRecords, MetaIsNeeded, &
                    .false., .false., &
                    Raw, size(Raw, 1), size(Raw, 2), PeriodRecords, &
                    EmbBiometDataExist, skip_period, LatestRawFileIndx, Col, &
                    .false.)
                if (skip_period) cycle pf_periods_loop

                !> Period skip control with message
                MissingRecords = dfloat(MaxPeriodNumRecords - PeriodRecords) &
                    / dfloat(MaxPeriodNumRecords) * 100d0
                if (PeriodRecords > 0 .and. MissingRecords > RPsetup%max_lack) &
                    cycle pf_periods_loop

                !> Filter raw data for user-defined flags
                if (RPsetup%filter_by_raw_flags) &
                    call FilterDatasetForFlags(Col, &
                        Raw, size(Raw, 1), size(Raw, 2))

                !***************************************************************
                !**** RAW FILE IMPORT FINISHES HERE. ***************************
                !**** NOW STARTS DATASET DEFINITION  ***************************
                !***************************************************************

                !> Allocate arrays for actual data processing
                if (.not. allocated(E2Set)) &
                    allocate(E2Set(PeriodRecords, E2NumVar))
                if (.not. allocated(DiagSet)) &
                    allocate(DiagSet(PeriodRecords, MaxNumDiag))

                !> Define EddyPro set of variables for the following processing
                call DefineE2Set(Col, Raw,   size(Raw, 1),     Size(Raw, 2), &
                                    E2Set,   size(E2Set, 1),   Size(E2Set, 2), &
                                    DiagSet, size(DiagSet, 1), Size(DiagSet, 2))
!                if (allocated(DiagSet))  deallocate(DiagSet)

                !> Clean up E2Set, eliminating values that are clearly un-physical
                call CleanUpE2Set(E2Set, size(E2Set, 1), size(E2Set, 2))

                !> Define as not present, variables for which
                !> too many values are outranged
                call EliminateCorruptedVariables(E2Set, &
                    size(E2Set, 1), size(E2Set, 2), skip_period, .false.)

                !> If either u, v or w have been eliminated,
                !> stops processing this period
                if (skip_period) then
                    if(allocated(E2Set)) deallocate(E2Set)
                    if(allocated(DiagSet)) deallocate(DiagSet)
                    cycle pf_periods_loop
                end if

                !> Update metadata if dynamic metadata are to be used
                if (EddyProProj%use_dynmd_file) &
                    call RetrieveDynamicMetadata(tsEnd, &
                        E2Col, size(E2Col))

                !> Override users choices if needed
                call OverrideSettings()

                !***************************************************************
                !**** DATASET DEFINITION FINISHES HERE. ************************
                !**** NOW STARTS RAW DATA REDUCTION     ************************
                !***************************************************************
                !> Filter only for sonic diagnostics (IRGA is irrelevant in
                !> planar fit)
                if (NumDiag > 0) then
                    call FilterDatasetForDiagnostics(E2Set, size(E2Set, 1), &
                        size(E2Set, 2), DiagSet, &
                        size(DiagSet, 1), size(DiagSet, 2), &
                        DiagAnemometer, .false.)
                end if
                if(allocated(DiagSet)) deallocate(DiagSet)

                !> Adjust coordinate systems if the case
                call AdjustSonicCoordinates(E2Set, &
                    size(E2Set, 1), size(E2Set, 2))

                !> Filter for wind direction if requested
                if (RPSetup%apply_wdf) &
                    call FilterDatasetForWindDirection(E2Set, size(E2Set, 1), size(E2Set, 2))

                !> Calculate basic stats
                call BasicStats(E2Set, size(E2Set, 1), size(E2Set, 2), &
                    1, .false.)
                Stats1 = Stats

                !> Calculate raw screening flags and despike data if requeste
                auxTest = TestType(.true., .false., .false., .true., .false., &
                    .false., .false., .false., .false.)
                call StatisticalScreening(E2Set, &
                    size(E2Set, 1), size(E2Set, 2), auxTest, .false.)

                !> Define as not present, variables for which too
                !> many values are outranged
                call EliminateCorruptedVariables(E2Set, &
                    size(E2Set, 1), size(E2Set, 2), skip_period, .false.)

                !> If either u, v or w have been eliminated,
                !> stops processing this period
                if (skip_period) then
                    if(allocated(E2Set)) deallocate(E2Set)
                    cycle pf_periods_loop
                end if

                !> Calculate basic stats
                call BasicStats(E2Set, size(E2Set, 1), size(E2Set, 2), &
                    2, .false.)
                Stats2 = Stats
                Stats3 = Stats

                !> Angle-of-attack calibration
                call AoaCalibration(E2Set, size(E2Set, 1), size(E2Set, 2))

                !> Gill WindMaster w-boost
                if (RPsetup%calib_wboost) &
                    call ApplyGillWmWBoost(E2Set, size(E2Set, 1), size(E2Set, 2))

                !> Calculate basic stats
                call BasicStats(E2Set, size(E2Set, 1), size(E2Set, 2), &
                    4, .false.)
                Stats4 = Stats
                if (allocated(E2Set)) deallocate(E2Set)

                !> Store statistics needed for planar fit calcuations
                pfn = pfn + 1
                pfWind(pfn, u) = Stats4%Mean(u)
                pfWind(pfn, v) = Stats4%Mean(v)
                pfWind(pfn, w) = Stats4%Mean(w)
            end do pf_periods_loop
            write(*, '(a)')
            write(*, '(a)') ' Done.'

            !*******************************************************************
            !**** RAW DATA REDUCTION FINISHES HERE.  ***************************
            !**** NOW STARTS PLANAR FIT CALCULATIONS ***************************
            !*******************************************************************

            !> Allocate sector-wise wind array
            if (.not. allocated(pfWindBySect)) &
                allocate(pfWindBySect(pfn, 3, PFSetup%num_sec))
            if (.not. allocated(GoPlanarFit)) &
                allocate(GoPlanarFit(PFSetup%num_sec))

            !> Check if wind components are within specified limits
            where (dsqrt(pfWind(1:pfn, u)**2 + pfWind(1:pfn, v)**2) < PFSetup%u_min &
                   .or. dsqrt(pfWind(1:pfn, u)**2 + pfWind(1:pfn, v)**2) > 20d0 &
                   .or. pfWind(1:pfn, w) > PFSetup%w_max)
                pfWind(1:pfn, u) = error
                pfWind(1:pfn, v) = error
                pfWind(1:pfn, w) = error
            end where

            !> Counts and log excluded stats lines
            err_cnt1 = 0
            do i = 1, pfn
                if(pfWind(i, u) == error .or. pfWind(i, v) == error &
                    .or. pfWind(i, w) == error) err_cnt1 = err_cnt1 + 1
            end do
            !if (err_cnt1 /= 0) then
                !Insert here call to ExceptionHandle and notify
                !that at least 1 wind data was excluded
            !end if

            !> Sort wind data according to wind sector
            !> (from pfWind to pfWindBySect)
            call SortWindBySector(pfWind(1:pfn, u:w), pfn, &
                pfNumElem, pfWindBySect)
            deallocate(pfWind)

            !> Some logging
            write(LogInteger, '(i6)') PFSetup%num_sec
            write(*, '(a)') ' Calculating planar fit rotation matrices for ' &
                                  // trim(adjustl(LogInteger)) // ' sector(s).'

            !> Loop over wind sectors
            GoPlanarFit = .true.
            secloop: do sec = 1, PFSetup%num_sec
                write(*, '(a, i2, a)', advance = 'no') '  Sector n.', sec, '..'
                if (PFSetup%wsect_exclude(sec)) then
                    GoPlanarFit(sec) = .false.
                    PFb(:, sec) = error
                    PFMat(:, :, sec) = error
                    call ExceptionHandler(41)
                    cycle secloop
                end if
                if(pfNumElem(sec) < PFSetup%min_per_sec) then
                    GoPlanarFit(sec) = .false.
                    PFb(:, sec) = error
                    PFMat(:, :, sec) = error
                    call ExceptionHandler(33)
                    cycle secloop
                end if

                allocate(pfWind(pfNumElem(sec), 3))
                do i = 1, pfNumElem(sec)
                    pfWind(i, :) = pfWindBySect(i, :, sec)
                end do

                !> Calculate auxiliary variables (see, e.g., Van Dijk et al.2004)
                call PlanarFitAuxParams(pfWind, pfNumElem(sec), Mat, pfVec)
                deallocate(pfWind)

                if (Meth%rot(1:len_trim(Meth%rot)) == 'planar_fit') then
                    !> Invert matrix --> Mat^(-1)
                    call MatrixInversion(Mat, 3, SingMat)

                    !> If singular matrix is found, set results to error
                    if (SingMat) then
                        GoPlanarFit(sec) = .false.
                        PFb(:, sec) = error
                        PFMat(:, :, sec) = error
                        call ExceptionHandler(34)
                        cycle secloop
                    end if

                    !> Calculate plane coefficients: PFb = Mat^(-1) * pfVec
                    PFb(:, sec) = 0d0
                    do i = u, w
                        PFb(:, sec) = PFb(:, sec) + dble(Mat(:, i)) * pfVec(i)
                    end do
                elseif (Meth%rot(1:len_trim(Meth%rot)) &
                        == 'planar_fit_no_bias') then
                    !> Define tensors of 2 elements (out of the 3-elements ones)
                    Mat2d(1,1:2) = Mat(2, 2:3)
                    Mat2d(2,1:2) = Mat(3, 2:3)
                    pfVec2d(1:2) = pfVec(2:3)

                    !> Invert matrix --> Mat^(-1)
                    call MatrixInversion(Mat2d, 2, SingMat)

                    !> If singular matrix is found, set results to error
                    if (SingMat) then
                        GoPlanarFit(sec) = .false.
                        PFb(:, sec) = error
                        PFMat(:, :, sec) = error
                        call ExceptionHandler(34)
                        cycle secloop
                    end if

                    !> Calculate plane coefficients: PFb = Mat^(-1) * pfVec
                    PFb2d(:, sec) = 0d0
                    do i = 1, 2
                        PFb2d(:, sec) = PFb2d(:, sec) &
                            + dble(Mat2d(:, i)) * pfVec2d(i)
                    end do
                    PFb(1, sec) = 0d0
                    PFb(2:3, sec) = PFb2d(1:2, sec)
                end if

                !> Calculate PP (PF rotation matrix, see
                !> Wilczak et al. 2001, BLM)
                call PlanarFitRotationMatrix(sec, PP)

                !> Update sector-wise rotation matrix
                PFMat(:, :, sec) = PP

                write(*, '(a)') ' Done.'
            end do secloop

            !> Fix sectors without calculations, using closest
            !> sector in the angular direction defined by user
            if (index(PFSetup%fix,'clockwise') /= 0) &
                call FixPlanarfitSectors(GoPlanarFit, size(GoPlanarFit))

            !> Write planar fit results on output file
            if (PFSetup%num_sec >= 1) &
                call WriteOutPlanarFit(pfNumElem, PFSetup%num_sec)

            if (allocated (pfWindBySect)) deallocate(pfWindBySect)
            if (allocated (pfNumElem)) deallocate(pfNumElem)
            write(*,'(a)') ' Planar Fit session terminated.'
            write(*,'(a)')
        end if
    else
        if (.not. allocated(GoPlanarFit)) allocate(GoPlanarFit(PFSetup%num_sec))
    end if

    !***************************************************************************
    !***************************************************************************
    !********************** PART COMMON TO NEXT TWO BIG LOOPS ******************
    !***************************************************************************
    !***************************************************************************

    !> Create TimeSeries for actual raw data processing
    if (EddyProProj%subperiod) then
        !> If user selected a sub-period, create time series corresponding to
        !> that period and verify that there is any overlap with RawTimeSeries
        call DateTimeToDateType(EddYProProj%start_date, &
            EddYProProj%start_time, SelectedStartTimestamp)
        call DateTimeToDateType(EddYProProj%end_date, &
            EddYProProj%end_time, SelectedEndTimestamp)

        NumberOfPeriods = NumOfPeriods(SelectedStartTimestamp, &
            SelectedEndTimestamp, DateStep)
        allocate(MasterTimeSeries(NumberOfPeriods + 1))
        call CreateTimeSeries(SelectedStartTimestamp, SelectedEndTimestamp, &
            DateStep, MasterTimeSeries, size(MasterTimeSeries), .false.)

        !> Verify at least partial overlap
        if (MasterTimeSeries(1) > RawTimeSeries(size(RawTimeSeries)) &
            .or. MasterTimeSeries(size(MasterTimeSeries)) < RawTimeSeries(1)) &
            call ExceptionHandler(46)
    else
        allocate(MasterTimeSeries(size(RawTimeSeries)))
        MasterTimeSeries = RawTimeSeries
    end if

    !> ONLY TEMPORARY: you can now actually eliminate rpStartTimestampIndx and
    !> rpEndTimestampIndx
    rpStartTimestampIndx = 1
    rpEndTimestampIndx = size(MasterTimeSeries)

    !***************************************************************************
    !***************************************************************************
    !******************** DEFINITION OF CALIBRATION EVENTS *********************
    !***************************************************************************
    !***************************************************************************
    if (DriftCorr%method /= 'none' .and. nCalibEvents > 0) then
        write(*,'(a)') ' Elaborating IRGA calibration-check history..'

        !> Loop on periods to be processed
        pcount = rpStartTimestampIndx - 1
        LatestRawFileIndx = 1
        latestCleaning = 0
        loggedDate = 'none'
        drift_loop: do
            pcount = pcount + 1

            !> If embedded metadata are to be used,
            !> reinitialize column information to null
            if (EddyProProj%use_extmd_file) then
                Col = BypassCol
            else
                Col = NullCol
            end if

            !> Normal exit instruction: either the last period was
            !> dealt with, or raw files are finished
            if (pcount > rpEndTimestampIndx - 1) exit drift_loop

            !> Normal exit instruction: if all cleaning
            !> events have been processed
            if (latestCleaning >= nCalibEvents) exit drift_loop

            !> Define initial/final timestamps of
            !> current period (say, [8:00 to 8:30))
            tsStart = MasterTimeSeries(pcount)
            tsEnd   = MasterTimeSeries(pcount + 1)

            !> If files are finished, keep going until the end of the selected
            !> period
            if (LatestRawFileIndx > NumRawFiles) cycle drift_loop

            !> Search file containing data starting from the time
            !> closest to tsStart
            !> Searches only from most current file onward, to avoid wasting time
            call FirstFileOfCurrentPeriod(tsStart, &
                tsEnd, RawFileList, NumRawFiles, &
                LatestRawFileIndx, NextRawFileIndx, skip_period)
            if (skip_period) cycle drift_loop

            !> Identify if current period is a cleaning event. If not, skip it
            call tsRelaxedMatch(tsStart, &
                Calib(latestCleaning + 1: nCalibEvents)%ts, &
                nCalibEvents - latestCleaning, &
                datetype(0, 0, 0, 3, 0), 'strictly later', clean)

            !> Identify if current period is closest possible to a cleaning event
            call tsRelaxedMatch(tsStart, &
                Calib(latestCleaning + 1: nCalibEvents)%ts, &
                nCalibEvents - latestCleaning, &
                datetype(0, 0, 0, 3, 0), 'strictly before', dirty)

            !> Cycle if file is not relevant to anything
            if (pcount /= rpStartTimestampIndx &
                .and. clean <= 0 .and. dirty <= 0) then
                LatestRawFileIndx = LatestRawFileIndx + 1
                cycle drift_loop
            end if

            !> Log out if the case
            if (pcount /= rpStartTimestampIndx) then
                call DateTypeToDateTime(tsStart, date, time)
                if (date /= loggedDate) then
                    write(*, '(a)') &
                        '  Calibration-check data found on: ' // date(1:10)
                    loggedDate = date
                end if
            end if

            !> Import dataset for current period. If using embedded biomet,
            !> also read biomet data. On entrance, NextRawFileIndx contains
            !> the index of the file to start the current period with
            !> On exit, LatestRawFileIndx contains the index of the
            !> latest file used
            call ImportCurrentPeriod(tsStart, tsEnd, &
                RawFileList, NumRawFiles, NextRawFileIndx, BypassCol, &
                MaxNumFileRecords, MetaIsNeeded, &
                EddyProProj%biomet_data == 'embedded', .false., Raw, &
                size(Raw, 1), size(Raw, 2), PeriodRecords, EmbBiometDataExist, &
                skip_period, LatestRawFileIndx, Col, .false.)

            !> Period skip control
            if (skip_period) cycle drift_loop

            !> Period skip control with message
            MissingRecords = dfloat(MaxPeriodNumRecords - PeriodRecords) &
                / dfloat(MaxPeriodNumRecords) * 100d0
            if (PeriodRecords > 0 &
                .and. MissingRecords > RPsetup%max_lack) cycle drift_loop

            !> Calculate reference counts
            call ReferenceCounts(dble(Raw), size(Raw, 1), size(Raw, 2))

            !> Special case of first file in the dataset: used to initialize
            !> drift history assuming cleaned instrument at the beginning
            if (pcount == rpStartTimestampIndx) then
                Calib(0)%ts = MasterTimeSeries(rpStartTimestampIndx)
                call DateTypeToDateTime(Calib(0)%ts, Calib(0)%date, Calib(0)%time)
                Calib(0)%ri(co2:h2o) = refCounts(co2:h2o)
                cycle drift_loop
            end if

            !> Case of cleaning event
            !> Assign relevant ri to current Calib dataset
            if (clean > 0) then
                Calib(latestCleaning + clean)%ri(co2:h2o) = refCounts(co2:h2o)
                latestCleaning = latestCleaning + clean
            end if
            !> Case of most dirty file (right before next cleaning event)
            !> Calculate and assign relevant quantities to current Calib dataset
            if (dirty > 0) then
                Calib(latestCleaning)%rf(co2:h2o) = refCounts(co2:h2o)
            end if
        end do drift_loop

        !> Rearrange so that all data needed between
        !> t1 and t2 are stored in Calib(t2)
        tmpCalib = Calib
        do i = 0, nCalibEvents - 1
            !> Calculate number of periods between
            !> each calibration-check event
            Calib(i+1)%numPeriods = &
                NumOfPeriods(Calib(i)%ts, Calib(i+1)%ts, DateStep)

            Calib(i+1)%ri = tmpCalib(i)%ri
            Calib(i+1)%rf = tmpCalib(i)%rf
        end do
        !> For Calib(0) (beginning of dataset), set at clean instrument
        Calib(0)%offset = 0d0
        Calib(0)%ri = error
        Calib(0)%rf = error

!> Only needed and valid for ICOS dataset, where H2O does not
!> start from "clean" but with some offset on day 1, Aug. 3.
!> Artificially set initial ri to the mean value at July 23,
!> when H2O signal was actually "clean", i.e. gives
!> same concentration of LI-7000.
!Calib(1)%ri(h2o) = 34703.78d0

    end if

    !***************************************************************************
    !***************************************************************************
    !***************************** RAW DATA PROCESSING *************************
    !***************************************************************************
    !***************************************************************************

    !> Start loop over all periods contained in the MasterTimeSeries
    NumberOfOkPeriods = 0
    pcount = rpStartTimestampIndx - 1
    LatestRawFileIndx = 1
    bLastFile = 1
    bLastRec = 0
    initialize = .true.
    initializeBiometOut = .true.
    initializeFluxnetOut = .true.
    InitializeStorage = .true.
    InitOutVarPresence = .true.
    DynamicMetadata = ErrDynamicMetadata
    InitGas4CalRefCol = Gas4CalRefCol

    periods_loop: do
        Gas4CalRefCol = InitGas4CalRefCol

        !***********************************************************************
        !**** RAW FILE IMPORT **************************************************
        !***********************************************************************

        if (pcount == rpStartTimestampIndx - 1) then
            !> Some log out
            if (EddyProProj%run_mode /=  'md_retrieval') then
                call hms_delta_print(' Start raw data processing: ', '')
            else
                call hms_delta_print(' Start metadata retrieving: ', '')
            end if
            write(*, '(a)') ' Processing time period:'
            call DateTypeToDateTime(MasterTimeSeries(rpStartTimestampIndx), &
                tmpDate, tmpTime)
            write(*, '(a)') '  Start: ' // tmpDate // ' ' // tmpTime
            call DateTypeToDateTime(MasterTimeSeries(rpEndTimestampIndx - 1) &
                + DateStep , tmpDate, tmpTime)
            write(*, '(a)') '    End: ' // tmpDate // ' ' // tmpTime
            write(TmpString1, '(i7)') &
                rpEndTimestampIndx - rpStartTimestampIndx
            write(*, '(a)') '  Total number of flux averaging periods: ' &
                // trim(adjustl(TmpString1))
            write(*, '(a)')
        end if
        pcount = pcount + 1

        !> If embedded metadata are to be used, reinitialize
        !> column information to null
        if (EddyProProj%use_extmd_file) then
            Col = BypassCol
        else
            Col = NullCol
        end if

        !> Initialize biomet variables to be used in computations
        biomet%val = error

        !> Normal exit instruction: either the last period was dealt with,
        !> or raw files are finished
        if (pcount > rpEndTimestampIndx - 1) exit periods_loop

        !> Define initial/final timestamps of current period (say, 8:00 to 8:29)
        tsStart = MasterTimeSeries(pcount)
        tsEnd   = MasterTimeSeries(pcount + 1)

        !> Associate timestamp of end of the period to current Stats
        call DateTypeToDateTime(tsStart, date, time)
        call DateTypeToDateTime(tsStart, Stats%start_date, Stats%start_time)
        call DateTypeToDateTime(tsEnd, Stats%date, Stats%time)

        !> Some logging
        if (EddyProProj%run_mode /= 'md_retrieval') then
            write(*, '(a)')
            call hms_current_print(' ',': processing new &
                &flux averaging period', .true.)
            write(*, '(a)') ' From: ' &
                // trim(date)   // ' ' // trim(time)
            write(*, '(a)') '   To: ' &
                // trim(Stats%date) // ' ' // trim(Stats%time)
        end if

        !> Define initial part of each output string
        call DateTimeToDOY(Stats%date, Stats%time, int_doy, float_doy)
        write(char_doy, *) float_doy
        call ShrinkString(char_doy)
        suffixOutString =  trim(Stats%date) // ',' // trim(Stats%time) &
                   // ',' // char_doy(1: index(char_doy, '.')+ 4)

        !> Only for external biomet files: retrieve biomet data for current
        !> period and write on output. Even if current period is skipped,
        !> biomet data will be on output.
        if (index(EddyProProj%biomet_data, 'ext_') /= 0) then
            call BiometRetrieveExternalData(bFileList, size(bFileList), &
                bLastFile, bLastRec, tsStart, &
                tsEnd, BiometDataFound, .true.)
            call WriteOutBiomet(suffixOutString, .false.)
        end if

        !> If files are finished, keep going until the end of the selected
        !> period
        if (LatestRawFileIndx > NumRawFiles) then
            if (EddyProProj%run_mode /= 'md_retrieval') then
                call ExceptionHandler(53)
                if (EddYProProj%out_fluxnet) call WriteOutFluxnetOnlyBiomet()
            end if
            call hms_delta_print(PeriodSkipMessage,'')
            cycle periods_loop
        end if

        !> Search file containing data starting from the time
        !> closest to tsStart. Searches only from most current
        !> file onward, to avoid wasting time
        call FirstFileOfCurrentPeriod(tsStart, tsEnd, &
            RawFileList, NumRawFiles, LatestRawFileIndx, &
            NextRawFileIndx, skip_period)

        Essentials%fname = trim(adjustl(RawFileList(NextRawFileIndx)%name))

        suffixOutString =  trim(adjustl(RawFileList(NextRawFileIndx)%name)) &
            // ',' // suffixOutString

        !> Exception handling
        if (skip_period) then
            if (EddyProProj%run_mode /= 'md_retrieval') then
                call ExceptionHandler(53)
                if (EddYProProj%out_fluxnet) call WriteOutFluxnetOnlyBiomet(suffixOutString)
            end if
            call hms_delta_print(PeriodSkipMessage,'')
            cycle periods_loop
        end if

        !> Import dataset for current period. If using embedded biomet, also
        !> read biomet data. On entrance, NextRawFileIndx contains index of
        !> file to start the current period with. On exit,
        !> LatestRawFileIndx contains the index of the latest file used
        call ImportCurrentPeriod(tsStart, tsEnd, RawFileList, &
            NumRawFiles, NextRawFileIndx, BypassCol, MaxNumFileRecords, &
            MetaIsNeeded, EddyProProj%biomet_data == 'embedded', .true., &
            Raw, size(Raw, 1), size(Raw, 2), PeriodRecords, &
            EmbBiometDataExist, skip_period, LatestRawFileIndx, Col, .true.)

        !> If it's running in metadata retriever mode,
        !> create a dummy dataset 1 minute long
        if (EddyProProj%run_mode == 'md_retrieval') then
            PeriodRecords = nint(Metadata%ac_freq * Metadata%file_length * 60d0)
            Raw = 1d0
            NumUserVar = 0
        else
            !> Retrieve biomet data for current period
            if (EddyProProj%biomet_data == 'embedded') then

                !> Retrieve biomet data from the already created bSet
                !> Basically, here only convert units and perform average
                !> over the averaging period
                call BiometRetrieveEmbeddedData(EmbBiometDataExist, .true.)

                !> Open biomet output file in case of embedded biomet files
                if(initializeBiometOut .and. nbVars > 0) then
                    call InitBiometOut()
                    initializeBiometOut  = .false.
                end if
                !> Write biomet output
                call WriteOutBiomet(suffixOutString, .true.)
            end if

            if (.not. allocated(UserCol)) &
                allocate(UserCol(NumUserVar))
            call DefineVars(Col, size(Raw, 2), NumUserVar)

            if (initializeFluxnetOut .and. EddyProProj%out_fluxnet) then
                call InitFluxnetFile_rp()
                initializeFluxnetOut  = .false.
            end if

            !> Period skip control
            if (skip_period) then
                if (EddyProProj%out_fluxnet) call WriteOutFluxnetOnlyBiomet(suffixOutString)
                call hms_delta_print(PeriodSkipMessage,'')
                cycle periods_loop
            end if

            !> Number of valid records imported from raw files
            Essentials%n_in = &
                CountRecordsAndValues(dble(Raw), size(Raw, 1), size(Raw, 2))

            !> Some logging
            write(*, '(a, i6)') '  Number of valid records available for this period: ', Essentials%n_in

            !> Period skip control
            MissingRecords = dfloat(MaxPeriodNumRecords - Essentials%n_in) &
                / dfloat(MaxPeriodNumRecords) * 100d0
            if (Essentials%n_in > 0 .and. MissingRecords > RPsetup%max_lack) then
                if (EddYProProj%out_fluxnet) call WriteOutFluxnetOnlyBiomet(suffixOutString)
                call ExceptionHandler(58)
                call hms_delta_print(PeriodSkipMessage,'')
                cycle periods_loop
            end if

            !> Filter raw data for user-defined flags
            if (RPsetup%filter_by_raw_flags) &
                call FilterDatasetForFlags(Col, Raw, size(Raw, 1), size(Raw, 2))
            Essentials%n_after_custom_flags = &
                CountRecordsAndValues(dble(Raw), size(Raw, 1), size(Raw, 2))

            !> Period skip control
            MissingRecords = dfloat(MaxPeriodNumRecords - Essentials%n_after_custom_flags) &
                / dfloat(MaxPeriodNumRecords) * 100d0
            if (MissingRecords > RPsetup%max_lack) then
                if (EddYProProj%out_fluxnet) call WriteOutFluxnetOnlyBiomet(suffixOutString)
                call ExceptionHandler(58)
                call hms_delta_print(PeriodSkipMessage,'')
                cycle periods_loop
            end if

            !> If drift correction is to be performed with signal strength
            !> proxy, calculate mean refCounts for current period
            if (DriftCorr%method == 'signal_strength') &
                call ReferenceCounts(dble(Raw), size(Raw, 1), size(Raw, 2))
        end if

        !***********************************************************************
        !**** RAW FILE IMPORT FINISHES HERE. NOW STARTS DATASET DEFINITION *****
        !***********************************************************************

        !> Allocate arrays for actual data processing
        if (.not. allocated(E2Set))    &
            allocate(E2Set(PeriodRecords, E2NumVar))
        if (.not. allocated(E2Primes)) &
            allocate(E2Primes(PeriodRecords, E2NumVar))
        if (.not. allocated(DiagSet))  &
            allocate(DiagSet(PeriodRecords, MaxNumDiag))

        !> Define EddyPro set of variables for the following processing
        call DefineE2Set(Col, Raw,   size(Raw, 1),     Size(Raw, 2), &
                            E2Set,   size(E2Set, 1),   Size(E2Set, 2), &
                            DiagSet, size(DiagSet, 1), Size(DiagSet, 2))

        !> Some convenient variables
        if (InitOutVarPresence) then
            OutVarPresent(u:E2NumVar) = E2Col(u:E2NumVar)%present
            InitOutVarPresence = .false.
        end if

        !> Define User set of variables, for main statistics
!        if (NumUserVar > 0) then
            if (.not. allocated(UserSet)) &
                allocate(UserSet(PeriodRecords, NumUserVar))
            if (.not. allocated(UserCol)) &
                allocate(UserCol(NumUserVar))
            if (.not. allocated(UserPrimes)) &
                allocate(UserPrimes(PeriodRecords, NumUserVar))
            call DefineUserSet(Col, Raw, size(Raw, 1), size(Raw, 2), &
                UserSet, size(UserSet, 1), size(UserSet, 2))
!        end if

        RowLags = 0
        if (EddyProProj%run_mode /= 'md_retrieval') then

            !> Update metadata if dynamic metadata are to be used
            if (EddyProProj%use_dynmd_file) &
                call RetrieveDynamicMetadata(tsEnd, E2Col, size(E2Col))

            !> Calculate relative separations between the analyzers
            !> and the anemometer used
            call DefineRelativeSeparations()

            !> Override users choices if needed
            call OverrideSettings()

            !> Determine whether it is day or night-time,
            call AssessDayTime(Stats%date, Stats%time)

            !*******************************************************************
            !**** DATASET DEFINITION FINISHES HERE. ****************************
            !**** STARTS RAW DATA REDUCTION         ****************************
            !*******************************************************************
            !> Interpret diagnostics and filter accordingly
            call InterpretLicorDiagnostics(DiagSet, &
                size(DiagSet, 1), size(DiagSet, 2))
            call FilterDatasetForDiagnostics(E2Set, size(E2Set, 1), &
                size(E2Set, 2), DiagSet, &
                size(DiagSet, 1), size(DiagSet, 2), &
                DiagAnemometer, .true.)
            if(allocated(DiagSet)) deallocate(DiagSet)

            !> Adjust coordinate systems if the case
            call AdjustSonicCoordinates(E2Set, size(E2Set, 1), size(E2Set, 2))

            !> Filter for wind direction if requested
            if (RPSetup%apply_wdf) &
                call FilterDatasetForWindDirection(E2Set, size(E2Set, 1), size(E2Set, 2))

            !> Number of valid records after filtering for wind direction
            Essentials%n_after_wdf = &
                CountRecordsAndValues(E2Set, size(E2Set, 1), size(E2Set, 2))
            PeriodActualRecords = Essentials%n_after_wdf
            
            !> Period skip control
            MissingRecords = dfloat(MaxPeriodNumRecords - Essentials%n_after_wdf) &
                / dfloat(MaxPeriodNumRecords) * 100d0
            if (MissingRecords > RPsetup%max_lack) then
                if (EddYProProj%out_fluxnet) call WriteOutFluxnetOnlyBiomet(suffixOutString)
                if(allocated(E2Set)) deallocate(E2Set)
                if(allocated(E2Primes)) deallocate(E2Primes)
                if(allocated(UserSet)) deallocate(UserSet)
                if(allocated(UserPrimes)) deallocate(UserPrimes)
                call ExceptionHandler(58)
                write(*,*)''
                call hms_delta_print(PeriodSkipMessage,'')
                cycle periods_loop
            end if

            !> Generate cell temperature dataset if the case, using either
            !> (1) native cell temperature, (2) weighted average of ti1 and ti2,
            !> (3) either ti1 or ti2 depending on availability
            call GenerateTcell(E2Set, size(E2Set, 1), size(E2Set, 2))

            !> Filter Tcell to simulate slower response temperature measurement
            !> for conversion to mixing ratio
            !if (RPsetup%tcell_filter_tconst /= 0) &
            !call CRA(E2Set, size(E2Set, 1), size(E2Set, 2), Metadata%ac_freq, &
            !    RPsetup%tcell_filter_tconst, tc)
        end if

        !> Now that variables have been properly assigned, can initialize
        !> main output files. This is done also if run is in
        !> metadata retriever mode
        if(initialize) then
            call InitOutFiles_rp()
            initialize = .false.
        end if

        !> Output first level of stats
        if (EddyProProj%run_mode /= 'md_retrieval') then

            !> ===== 1. RAW DATA AS READ FROM FILES ============================
            !> Output raw dataset first level
            if (RPsetup%out_raw(1)) call OutRawData(Stats%date, Stats%time, &
                E2Set, size(E2Set, 1), size(E2Set, 2), 1)
            !> Calculate basic stats and output them as requested
            call BasicStats(E2Set, size(E2Set, 1), size(E2Set, 2), 1, .true.)
            Stats1 = Stats
            if (RPsetup%out_st(1)) &
                call WriteOutStats(ust1, Stats1, suffixOutString, PeriodRecords)
            if (NumUserVar > 0) then
                call UserBasicStats(UserSet, &
                    size(UserSet, 1), size(UserSet, 2), 1)
                if (RPsetup%out_st(1)) &
                    call WriteOutUserStats(u_user_st1, suffixOutString, &
                        PeriodRecords, AddUserStatsHeader)
                    AddUserStatsHeader = .false.
            end if

            !> Based on mean value, if sonic (or fast) temperature
            !> is out-ranged, search alternative one.
            if (Stats1%Mean(ts) < 220d0 .or. Stats1%Mean(ts) > 340d0) &
                call ReplaceSonicTemperature(E2Set, size(E2Set, 1), &
                    size(E2Set, 2), UserSet, size(UserSet, 1), size(UserSet, 2))

            !> ===== 2. STATISTICAL SCREENING ==================================
            !> Calculate raw screening flags and despike data if requested
            call StatisticalScreening(E2Set, &
                size(E2Set, 1), size(E2Set, 2), Test, .true.)
            if (NumUserVar > 0) call DespikeUserSet(UserSet, &
                size(UserSet, 1), size(UserSet, 2))

            !> Define as not present, variables for which
            !> too many values are out-ranged
            call EliminateCorruptedVariables(E2Set, &
                size(E2Set, 1), size(E2Set, 2), skip_period, .true.)

            !> If either u, v or w have been eliminated,
            !> stops processing this period
                if (skip_period) then
                if (EddYProProj%out_fluxnet) call WriteOutFluxnetOnlyBiomet(suffixOutString)
                if(allocated(E2Set)) deallocate(E2Set)
                if(allocated(E2Primes)) deallocate(E2Primes)
                if(allocated(UserSet)) deallocate(UserSet)
                if(allocated(UserPrimes)) deallocate(UserPrimes)
                call ExceptionHandler(59)
                write(*,*)''
                call hms_delta_print(PeriodSkipMessage,'')
                cycle periods_loop
            end if

            !> If got until here, incrase number of ok periods
            NumberOfOkPeriods = NumberOfOkPeriods + 1
            
            !> Count values available for each variable and value pairs 
            !> available for each main w-covariance
            !>> 
            Essentials%n = ierror
            Essentials%n_wcov = ierror
            !> Wind data
            Essentials%n(w) = &
                CountRecordsAndValues(E2Set, size(E2Set, 1), size(E2Set, 2), w)
            Essentials%n_wcov(u) = &
                CountRecordsAndValues(E2Set, size(E2Set, 1), size(E2Set, 2), w, u)
            !> Gas data
            do j = ts, gas4
                if (E2Col(j)%present) then
                    Essentials%n(j) = &
                        CountRecordsAndValues(E2Set, size(E2Set, 1), size(E2Set, 2), j)
                    Essentials%n_wcov(j) = &
                        CountRecordsAndValues(E2Set, size(E2Set, 1), size(E2Set, 2), w, j)
                end if
            end do

            !> If a 4th gas calibration has to be done (using a 'cal-ref'
            !> column from UserCol) does so. Note that so far the calibration
            !> procedure is fully customized on the needs of a
            !> specific O3 analyzer
            call CalibrateGas4(E2Set, size(E2Set, 1), size(E2Set, 2))

            !> Output raw dataset second level
            if (RPsetup%out_raw(2)) call OutRawData(Stats%date, Stats%time, &
                E2Set, size(E2Set, 1), size(E2Set, 2), 2)
            !> Calculate basic stats and output them as requested
            call BasicStats(E2Set, size(E2Set, 1), size(E2Set, 2), 2, .true.)
            Stats2 = Stats
            if (RPsetup%out_st(2)) &
                call WriteOutStats(ust2, Stats2, suffixOutString, PeriodRecords)
            if (NumUserVar > 0) then
                call UserBasicStats(UserSet, &
                    size(UserSet, 1), size(UserSet, 2), 2)
                if (RPsetup%out_st(2)) &
                    call WriteOutUserStats(u_user_st2, suffixOutString, &
                        PeriodRecords, AddUserStatsHeader)
                    AddUserStatsHeader = .false.
            end if

            !> ===== 3. CROSS-WIND CORRECTION ==================================
            !> Apply raw-level cross wind correction
            !> (after Liu et al. 2001), if requested
            if (RPsetup%calib_cw) then
                call CrossWindCorr(E2Col(u), E2Set, &
                    size(E2Set, 1), size(E2Set, 2), .true.)
            else
                write(*,'(a)') '  Cross-wind correction not requested &
                    &or not applicable'
            end if

            !> Output raw dataset third level
            if (RPsetup%out_raw(3)) call OutRawData(Stats%date, Stats%time, &
                E2Set, size(E2Set, 1), size(E2Set, 2), 3)
            !> Calculate basic stats and output them as requested
            call BasicStats(E2Set, size(E2Set, 1), size(E2Set, 2), 3, .true.)
            Stats3 = Stats
            if (RPsetup%out_st(3)) &
                call WriteOutStats(ust3, Stats3, suffixOutString, PeriodRecords)
            if (NumUserVar > 0) then
                call UserBasicStats(UserSet, &
                    size(UserSet, 1), size(UserSet, 2), 3)
                if (RPsetup%out_st(3)) &
                    call WriteOutUserStats(u_user_st3, suffixOutString, &
                        PeriodRecords, AddUserStatsHeader)
                    AddUserStatsHeader = .false.
            end if

            !> ===== 4. ANGLE OF ATTACK CORRECTION =============================
            !> Angle-of-attack calibration
            call AoaCalibration(E2Set, size(E2Set, 1), size(E2Set, 2))

            !> Gill WindMaster w-boost
            if (RPsetup%calib_wboost) &
                call ApplyGillWmWBoost(E2Set, size(E2Set, 1), size(E2Set, 2))

            !> Output raw dataset forth level
            if (RPsetup%out_raw(4)) call OutRawData(Stats%date, Stats%time, &
                E2Set, size(E2Set, 1), size(E2Set, 2), 4)
            !> Calculate basic stats and output them as requested
            call BasicStats(E2Set, size(E2Set, 1), size(E2Set, 2), 4, .true.)
            Stats4 = Stats
            if (RPsetup%out_st(4)) &
                call WriteOutStats(ust4, Stats4, suffixOutString, PeriodRecords)
            if (NumUserVar > 0) then
                call UserBasicStats(UserSet, &
                    size(UserSet, 1), size(UserSet, 2), 4)
                if (RPsetup%out_st(4)) &
                    call WriteOutUserStats(u_user_st4, suffixOutString, &
                        PeriodRecords, AddUserStatsHeader)
                    AddUserStatsHeader = .false.
            end if

            !> ===== 4.1 CORRECTION OF CALIBRATION DRIFTS ======================
            if (DriftCorr%method /= 'none' .and. nCalibEvents /= 0) &
                call DriftCorrection(E2Set, size(E2Set, 1), size(E2Set, 2), &
                    E2Col, size(E2Col), nCalibEvents, tsStart)

            !> Convert to mixing ratios (if WPL requested, and if the case)
            if (EddyProProj%wpl) &
                call PointByPointToMixingRatio(E2Set, &
                    size(E2Set, 1), size(E2Set, 2), .true.)

            !> ===== 5. TILT CORRECTION ========================================
            !> Apply rotations for tilt correction, if requested
            call TiltCorrection(Meth%rot, GoPlanarFit, E2Set, &
                size(E2Set, 1), size(E2Set, 2), PFSetup%num_sec, &
                Essentials%yaw, Essentials%pitch, Essentials%roll, .true.)

            !> Output raw dataset fifth level
            if (RPsetup%out_raw(5)) call OutRawData(Stats%date, Stats%time, &
                E2Set, size(E2Set, 1), size(E2Set, 2), 5)
            !> Calculate basic stats and output them as requested
            call BasicStats(E2Set, size(E2Set, 1), size(E2Set, 2), 5, .true.)
            Stats5 = Stats
            if (RPsetup%out_st(5)) &
                call WriteOutStats(ust5, Stats5, suffixOutString, PeriodRecords)
            if (NumUserVar > 0) then
                call UserBasicStats(UserSet, &
                    size(UserSet, 1), size(UserSet, 2), 5)
                if (RPsetup%out_st(5)) &
                    call WriteOutUserStats(u_user_st5, suffixOutString, &
                        PeriodRecords, AddUserStatsHeader)
                    AddUserStatsHeader = .false.
            end if

            !> ===== 6. TIMELAG COMPENSATION  ==================================
            !> If available, for files others than GHG, replace flow rate
            !> of LI-7200 provided by user with mean value from raw files
            if (EddyProProj%ftype /= 'licor_ghg' &
                .or. EddyProProj%use_extmd_file) then
                do i = 1, E2NumVar
                    if (NumUserVar > 0) then
                        do j = 1, NumUserVar
                            if (UserCol(j)%var == 'flowrate' &
                                .and. UserCol(j)%instr%model == E2Col(i)%instr%model &
                                .and. UserStats%Mean(j) /= 0d0 &
                                .and. UserStats%Mean(j) /= error) then
                                E2Col(i)%instr%tube_f = UserStats%Mean(j)
                                exit
                            end if
                        end do
                    end if
                end do
            end if

            !> Retrieving instruments parameters
            call RetrieveSensorParams()

            !> Defines plausible nominal timelags and timelag ranges
            call SetTimelags()

            !> Calculate and compensate time-lags
            if (TimeLagOptSelected) Meth%tlag = 'maxcov&default'
            call TimeLagHandle(Meth%tlag(1:len_trim(Meth%tlag)), E2Set, &
                size(E2Set, 1), size(E2Set, 2), Essentials%actual_timelag, &
                Essentials%used_timelag, Essentials%def_tlag, .false.)
            if (TimeLagOptSelected) Meth%tlag = 'tlag_opt'

            !> ===== 6.1 FILTERING MOLAR DENSITY DATA FOR ABSOLUTE LIMITS TEST  ====================
            if (EddyProProj%run_mode /= 'md_retrieval') then
                !> Estimate temperatures, pressures and relevant
                !> air molar volumes
                call AirAndCellParameters()
                if (Test%al .and. RPsetup%filter_al) then
                    !> Apply filter for absolute limits test, if the case
                    FilterWhat = .false.
                    FilterWhat(co2:gas4) = .true.
                    call FilterDatasetForPhysicalThresholds(E2Set, &
                        size(E2Set, 1), size(E2Set, 2), FilterWhat)
                    !> Define as not present, variables for which &
                    !> too many values are out-ranged
                    call EliminateCorruptedVariables(E2Set, &
                        size(E2Set, 1), size(E2Set, 2), skip_period, .true.)
                end if
            end if

            !> Output raw dataset sixth level
            if (RPsetup%out_raw(6)) call OutRawData(Stats%date, Stats%time, &
                E2Set, size(E2Set, 1), size(E2Set, 2), 6)
            !> Calculate basic stats and output them as requested
            call BasicStats(E2Set, size(E2Set, 1), size(E2Set, 2), 6, .true.)
            Stats6 = Stats
            if (RPsetup%out_st(6)) &
                call WriteOutStats(ust6, Stats6, suffixOutString, PeriodRecords)
            if (NumUserVar > 0) then
                call UserBasicStats(UserSet, &
                    size(UserSet, 1), size(UserSet, 2), 6)
                if (RPsetup%out_st(6)) &
                    call WriteOutUserStats(u_user_st6, suffixOutString, &
                        PeriodRecords, AddUserStatsHeader)
                    AddUserStatsHeader = .false.
            end if

            !> Quality check test for stationarity
            call StationarityTest(E2Set, size(E2Set, 1), size(E2Set, 2), StDiff)

            !> Calculate wind speed and maximum wind speed
            call MaxWindSpeed(E2Set, &
                size(E2Set, 1), size(E2Set, 2), Ambient%MWS)
            if (Stats6%mean(u) /= error .and. Stats6%mean(v) /= error &
                .and. Stats6%mean(w) /= error) then
                Ambient%WS = dsqrt(Stats6%mean(u)**2 &
                    + Stats6%mean(v)**2 + Stats6%mean(w)**2)
            else
                Ambient%WS = error
            end if

            !> ===== 6.2 QC tests =============================================
            !> Calculate Kurtosis Index on differenced variables
            call KID(E2Set(:, 1:GHGNumVar), size(E2Set, 1), GHGNumVar)

            !> Calculate Longest Gap Duration
            call LongestGapDuration(E2Set(:, 1:GHGNumVar), size(E2Set, 1), GHGNumVar)

            !> ===== 7. DETRENDING =============================================
            !> Calculate fluctuations based on chosen detrending method
            write(*, '(a)', advance = 'no') '  Detrending..'
            call Fluctuations(E2Set, E2Primes, &
                size(E2Set, 1), size(E2Set, 2), RPsetup%Tconst, Stats, E2Col)
            if (allocated(E2Set)) deallocate(E2Set)
            write(*,'(a)') ' Done.'
            if (NumUserVar > 0) then
                write(*, '(a)', advance = 'no') '  Detrending user set..'
                call UserFluctuations(UserSet, UserPrimes, &
                    size(UserSet, 1), size(UserSet, 2), &
                    RPsetup%Tconst, UserStats, UserCol)
                write(*,'(a)') ' Done.'
                if (allocated(UserSet)) deallocate(UserSet)
            end if

            !> Output raw dataset seventh level
            if (RPsetup%out_raw(7)) &
                call OutRawData(Stats%date, Stats%time, E2Primes, &
                    size(E2Primes, 1), size(E2Primes, 2), 7)
            !> Calculate basic stats and output them as requested
            call BasicStats(E2Primes, &
                size(E2Primes, 1), size(E2Primes, 2), 7, .true.)
            Stats7 = Stats
            if (RPsetup%out_st(7)) &
                call WriteOutStats(ust7, Stats7, suffixOutString, PeriodRecords)
            if (NumUserVar > 0) then
                call UserBasicStats(UserPrimes, &
                    size(UserPrimes, 1), size(UserPrimes, 2), 7)
                if (RPsetup%out_st(7)) &
                    call WriteOutUserStats(u_user_st7, suffixOutString, &
                        PeriodRecords, AddUserStatsHeader)
                    AddUserStatsHeader = .false.
            end if
            if (allocated(UserPrimes)) deallocate(UserPrimes)

            !> ===== 7.1 QC tests =============================================
            !> Fisher's test
            call Fisher(E2Primes(:, 1:GHGNumVar), size(E2Primes, 1), size(E2Primes, 2))

            !> Cross-correlation R^2 test for repeated values 
            call CrossCorrTest(E2Primes(:, 1:GHGNumVar), size(E2Primes, 1), size(E2Primes, 2))

            !> Calculate Mahrt's random error and Nonstationarity ratio anyway.
            call RU_Mahrt_98(E2Primes, size(E2Primes, 1), size(E2Primes, 2))

            !> If requested, estimate random error
            call RandomUncertaintyHandle(E2Primes, size(E2Primes, 1), size(E2Primes, 2))

            !*******************************************************************
            !**** RAW DATA REDUCTION FINISHES HERE *****************************
            !**** CALCULATE AND OUTPUT CO-SPECTRA  *****************************
            !*******************************************************************
            if (RPsetup%do_spectral_analysis) then
                SpecCol = E2Col

                !> Replace gaps with linear interpolation of neighbouring data
                call FixDatasetForSpectra(E2Primes, &
                    size(E2Primes, 1), size(E2Primes, 2), N2)

                !> Set length of dataset by stripping
                !> trailing/leading error codes
                Nmax = maxval(RowLags)
                Nmin = minval(RowLags)
                max_nsmpl = N2 - (Nmax - Nmin)
                if(RPsetup%power_of_two) then
                    !> Calculate power-of-two closest to number of
                    !> available samples
                    call PowerOfTwo(max_nsmpl, SpecRow)
                else
                    !> use all samples
                    SpecRow = max_nsmpl
                end if

                Nmin = - Nmin
                !> Recalculate basic statistics with PeriodRecords=SpecRow
                !> as number of observations
                allocate(SpecSet(SpecRow, E2NumVar))
                SpecSet(1: SpecRow, :) = E2Primes(Nmin + 1: Nmin + SpecRow, :)
                call BasicStats(SpecSet, &
                    size(SpecSet, 1), size(SpecSet, 2), 8, .true.)

                !> Calculate spectra and cospectra and output them all
                call SpectralAnalysis(Stats%date, Stats%time, bf, &
                    SpecSet(:, u:gas4), size(SpecSet, 1), gas4)
                if (allocated(SpecSet)) deallocate(SpecSet)
                if (allocated(E2Primes)) deallocate(E2Primes)

                !> Reset stats to Stats7, after the parenthesis
                !> of spectral analysis
                Stats = Stats7
            else
                Essentials%degH(:) = error
            end if
        end if
        if (allocated(E2Primes)) deallocate(E2Primes)
        if (allocated(UserPrimes)) deallocate(UserPrimes)
        if (allocated(UserSet)) deallocate(UserSet)

        !***********************************************************************
        !**** (CO)SPECTRA CALCULATION FINISHES HERE  ***************************
        !**** NOW STARTS FLUX COMPUTATION/CORRECTION ***************************
        !***********************************************************************
        if (EddyProProj%run_mode /= 'md_retrieval') then

            !> Average mole fractions in [umol mol_a-1] and [mmol mol_a-1]
            call MoleFractionsAndMixingRatios()

            !> Calculate parameters for flux computation
            call FluxParams(.true.)

            if (E2Col(ch4)%Instr%model(1:len_trim(E2Col(ch4)%Instr%model) - 2) &
                == 'li7700') then
                !> Calculate multipliers for LI-7700 spectroscopic correction
                call Multipliers7700(Stats%Pr, Ambient%Ta, Stats%chi(h2o), &
                    Mul7700%A, Mul7700%B, Mul7700%C)
                !> Modify mole fraction and mixing ratio to account for
                !> key(T,P), Eq. 6.13 of LI-7700 manual
                !> Uses multiplies A, because this is equal to key.
                Stats%chi(ch4) = Stats%chi(ch4) * Mul7700%A
                Stats%r(ch4)   = Stats%r(ch4)   * Mul7700%A
            end if

            !> Calculate LI-7500 surface heating correction if requested
            call BurbaTerms()

            !> Calculate fluxes at Level 0
            call Fluxes0_rp(.true.)

            !> As of now, still use CO2 analyzer software version as a proxy for
            !> Logger software version. However, the machinery is in place
            !> for using logger version from [Station], simply remove the
            !> following line.
!            if (E2Col(co2)%instr%sw_ver /= errSwVer) then
            Metadata%logger_swver = E2Col(co2)%instr%sw_ver
!            elseif (E2Col(h2o)%instr%sw_ver /= errSwVer) then
!                Metadata%logger_swver = E2Col(h2o)%instr%sw_ver
!            end if

            if (.not. EddyProProj%fcc_follows) then
                !> Low-pass and high-pass spectral correction factors
                call BandPassSpectralCorrections(E2Col(u)%Instr%height, &
                    Metadata%d, E2Col(u:gas4)%present, Ambient%WS, Ambient%Ta, &
                    Ambient%zL, Metadata%ac_freq, RPsetup%avrg_len, &
                    Metadata%logger_swver, Meth%det, &
                    RPsetup%Tconst, .true., E2Col(u:GHGNumVar)%instr, 1)

                !> Calculate fluxes at Level 1
                call Fluxes1_rp()

                !> Calculate fluxes at Level 2 and Level 3
                call Fluxes23_rp()

                !> Footprint estimation
                foot_model_used = Meth%foot(1:len_trim(Meth%foot))
                call FootprintHandle(Stats%Cov(w, w), Ambient%us, &
                    Ambient%zL, Ambient%WS, Ambient%L, &
                    E2Col(u)%Instr%height, Metadata%d, Metadata%z0)
            else
                Foot = errFootprint
                Flux1 = errFlux
                Flux2 = errFlux
                Flux3 = errFlux
                BPCF = errBPCF
                foot_model_used = 'none'
            end if

            !> Calculate storage terms
            if(InitializeStorage) then
                Stor%H  = error
                Stor%LE = error
                Stor%of = error
                InitializeStorage = .false.
            else
                call Storage(PrevStats, prevAmbient)
            end if
            PrevStats = Stats
            prevAmbient = Ambient
            prevBiomet = biomet

            !> Well developed turbulence conditions test,
            !> after Foken et al. (2004, Handbook of Microm.)
            call DevelopedTurbulenceTest(DtDiff)

            !> flagging the file, after Foken et al. (2004, Handbook of Microm.)
            call QualityFlags(Flux2, StDiff, DtDiff, STFlg, DTFlg, QCFlag, .true.)

            !> Write details on output files if requested
            if(RPsetup%out_qc_details .and. Meth%qcflag /= 'none') &
                call WriteOutQCDetails(suffixOutString, StDiff, DtDiff, STFlg, DTFlg)

            !> Update values of AGC and RSSI as available
            call SetLicorDiagnostics(NumUserVar)
        end if

        !>Write out full output file (main express output)
        if (EddyProProj%out_full) &
            call WriteOutFull(suffixOutString, PeriodRecords, PeriodActualRecords)

        !>Write out full output file (main express output)
        if (EddyProProj%out_md) &
            call WriteOutMetadata(suffixOutString)
            if (EddyProProj%out_fluxnet) &
                call WriteOutFluxnet(StDiff, DtDiff, STFlg, DTFlg)

        if (EddyProProj%run_mode /= 'md_retrieval') then
            call hms_delta_print('  Flux averaging period processing time: ','')
        else
            call hms_delta_print('  Metadata retrieving time: ','')
        end if
        write(*, *)

        if (allocated(UserCol))  deallocate(UserCol)
        if (allocated(E2Set))    deallocate(E2Set)
        if (allocated(E2Primes)) deallocate(E2Primes)
        if (allocated(DiagSet))  deallocate(DiagSet)
        if (allocated(UserSet))  deallocate(UserSet)
    end do periods_loop
    if (allocated(bf)) deallocate(bf)

    !***************************************************************************
    !**** FLUX COMPUTATION FINISHES HERE.                      *****************
    !**** NOW STARTS DATASET CREATION AND OUTPUT FILE HANDLING *****************
    !***************************************************************************
    close(ust1)
    close(ust2)
    close(ust3)
    close(ust4)
    close(ust5)
    close(ust6)
    close(ust7)
    close(u_user_st1)
    close(u_user_st2)
    close(u_user_st3)
    close(u_user_st4)
    close(u_user_st5)
    close(u_user_st6)
    close(u_user_st7)
    close(umd)
    close(uflx)
    close(ufnet_e)
    close(ufnet_b)
    close(uaflx)
    close(uex)
    close(uflxnt)
    close(ubiomet)
    close(uqc)

    !> If no averaging period was performed, return message and cancel tmp files
    if (NumberOfOkPeriods == 0 .and. EddyProProj%run_mode /= 'md_retrieval') then
        !> Delete files in output folder
        del_status = system(trim(comm_del) // ' "' // trim(adjustl(Dir%main_out)) &
            // '*' // Timestamp_FilePadding //'*"'  // comm_err_redirect)

        !> Alerting and closing run
        write(*,'(a)')
        call ExceptionHandler(35)
    end if

    !> Creating datasets from output files
    write(*, '(a)')
    write(*, '(a)') ' Raw data processing terminated. &
        &Creating continuous datasets if necessary..'

    if (make_dataset_common) then
        call CreateDatasetsCommon(MasterTimeSeries, size(MasterTimeSeries), &
            rpStartTimestampIndx, rpEndTimestampIndx, 'RP')
    else
        call RenameTmpFilesCommon()
    end if
    if (make_dataset_rp) then
        call CreateDatasetsRP(MasterTimeSeries, size(MasterTimeSeries), &
            rpStartTimestampIndx, rpEndTimestampIndx)
    else
        call RenameTmpFilesRP()
    end if
    write(*, '(a)') ' Done.'

    !> Edit .eddypro file updating path to ex_file
    call ForceSlash(FLUXNET_Path, .false.)
    call EditIniFile(trim(PrjPath), 'ex_file', &
        trim(FLUXNET_Path(1:index(FLUXNET_Path, '.tmp')-1)))

    if (EddyProProj%run_env /= 'embedded') &
        write(*, '(a)') ' FLUXNET file path: ' &
            // trim(FLUXNET_Path(1:index(FLUXNET_Path, '.tmp')-1))

    !> Copy ".eddypro" file into output folder
    if (.not. EddyProProj%fcc_follows) then
        call CopyFile(trim(adjustl(PrjPath)), &
        trim(adjustl(Dir%main_out)) // 'processing' &
        // Timestamp_FilePadding // '.eddypro')
    end if

    !> Delete tmp folder if running in embedded mode
    if(EddyProProj%run_env == 'desktop') &
        del_status = system(trim(comm_rmdir) // ' "' &
        // trim(adjustl(TmpDir)) // '"')

    if (.not. EddyProProj%fcc_follows) then
        write(*, '(a)') ''
        write(*, '(a)') ' ****************************************************'
        write(*, '(a)') ' Program EddyPro executed gracefully.'
        write(*, '(a)') ' Check results in the selected output directory.     '
        write(*, '(a)') ' ****************************************************'
    end if
    stop ''
end program EddyproRP
