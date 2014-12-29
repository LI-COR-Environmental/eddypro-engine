!***************************************************************************
! eddypro_fcc_main.f90
! --------------------
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
Program EddyproFCC
    use m_fx_global_var
    implicit none

    integer, external :: CreateDir
    integer :: i
    integer :: month
    integer :: nbins
    integer :: open_status
    integer :: nfit(GHGNumVar, 2)
    integer :: day
    integer :: STFlg(GHGNumVar)
    integer :: DTFlg(GHGNumVar)
    integer :: int_doy
    integer :: NumValidExRecords
    integer :: NumExRecords
    integer :: NumberOfPeriods
    integer :: fxStartTimestampIndx
    integer :: fxEndTimestampIndx
    integer :: saStartTimestampIndx
    integer :: saEndTimestampIndx
    integer :: NumFullFiles
    integer :: NumFullFilesNoRecurse
    integer :: NumBinnedFiles
    integer :: NumBinnedFilesNoRecurse
    integer :: fcount
    integer :: nrow_full

    logical :: skip
    logical :: InitializeOuputFiles
    logical :: DoSpectralAssessment
    logical :: skip_spectra
    logical :: skip_cospectra
    logical :: ValidRecord
    logical :: EndOfFileReached
    logical :: exEndReached

    !> Derived type variables
    type(DateType) :: fxStartTimestamp
    type(DateType) :: fxEndTimestamp
    type(DateType) :: saStartTimestamp
    type(DateType) :: saEndTimestamp
    type(DateType) :: binStartTimestamp
    type(DateType) :: binEndTimestamp
    type(DateType) :: SelectedStartTimestamp
    type(DateType) :: SelectedEndTimestamp
    type(Datetype) :: CurrentTimestamp
    type(SpectraSetType) :: BinSpec(MaxNumBins)
    type(SpectraSetType) :: BinCosp(MaxNumBins)
    type(SpectraSetType) :: BinCospForStable(MaxNumBins)
    type(SpectraSetType) :: BinCospForUnstable(MaxNumBins)
    type(InstrumentType) :: AuxInstrument(GHGNumVar)
    type(QCType) :: StDiff
    type(QCType) :: DtDiff
    type(EXType) :: lEx

    !> Allocatable variabled
    type(DateType), allocatable :: MasterTimeSeries(:)
    type(FitSpectraType), allocatable :: FitUnstable(:)
    type(FitSpectraType), allocatable :: FitStable(:)

    !> External functions
    integer, external :: NumOfPeriods
    integer, external :: NumberOfFilesInSubperiod
    real(kind = dbl), external :: func

    !*******************************************************************************
    !*******************************************************************************

    app = fcc_app

    !> Initialize environment
    write(*, '(a)')
    call InitEnv()

    write(*, '(a)') 'Starting flux computation and correction session..'
    write(*, '(a)')

    !> Read ".eddypro" file for both spectral analysis and flux correction
    call ReadIniFX('FluxCorrection')

    !> If running in embedded mode, override some settings
    if (EddyProProj%run_env == 'embedded') call ConfigureForEmbedded('EddyPro-FCC')

    !> Preliminarily read essential files and retrieve a few information
    call InitExVars(fxStartTimestamp, fxEndTimestamp, &
        NumExRecords, NumValidExRecords)

    !> If no good records are found stop execution
    if (NumValidExRecords <= 0) call ExceptionHandler(61)

    !> Retrieve NumberOfPeriods and allocate MasterTimeSeries
    NumberOfPeriods = NumOfPeriods(fxStartTimestamp, fxEndTimestamp, DateStep)
    allocate(MasterTimeSeries(NumberOfPeriods))

    !> Timestamps of selected subperiod
    call DateTimeToDateType(EddYProProj%start_date, EddYProProj%start_time, &
        SelectedStartTimestamp)
    call DateTimeToDateType(EddYProProj%end_date, EddYProProj%end_time, &
        SelectedEndTimestamp)
    SelectedEndTimestamp = SelectedEndTimestamp - DateStep

    !> Create timestamp array for full dataset
    call CreateMasterTimeSeries(fxStartTimestamp, fxEndTimestamp, DateStep, &
        MasterTimeSeries, size(MasterTimeSeries))

    !> Detect indexes of selected sub-period
    call tsExtractSubperiodIndexes(MasterTimeSeries, size(MasterTimeSeries), &
        SelectedStartTimestamp, SelectedEndTimestamp, &
        fxStartTimestampIndx, fxEndTimestampIndx)

    if (fxStartTimestampIndx == nint(error) .or. fxEndTimestampIndx == nint(error)) &
        call ExceptionHandler(50)

    !> For convenience, change timestamps of MasterTimeSeries by increasing of DateStep
    do i = 1, size(MasterTimeSeries)
        MasterTimeSeries(i) = MasterTimeSeries(i) + DateStep
    end do

    !> Initialize import of full cospectra files
    if (FCCsetup%import_full_cospectra) then
        call NumberOfFilesInDir(Dir%full, '.csv', .false., '', &
            NumFullFiles, NumFullFilesNoRecurse)
        allocate(FullFileList(NumFullFiles))

        !> Read names of full cospectra files
        call FileListByExt(Dir%full, '.csv', .true., FullFilePrototype, &
            .true., .true., .true., FullFilelist, size(FullFilelist), .true., '')

        !> Retrieve length of full cospectra for later allocation
        call FullCospectraLength(FullFilelist(1)%path, nrow_full)
    end if

    !****************************************************************
    !****************************************************************
    !*********** SPECTRAL ASSESSMENT SECTION IF REQUESTED ***********
    !****************************************************************
    !****************************************************************

    DoSpectralAssessment = FCCsetup%sa_onthefly .and. FCCsetup%SA%in_situ

    !> If spectral analysis must be performed or cospectral output are requested,
    !> start loop on cospectra files
    if (DoSpectralAssessment .or. EddyProProj%out_avrg_cosp) then
        if (DoSpectralAssessment) write(*, '(a)') &
            ' Starting "spectral assessment" session..'
        if (EddyProProj%out_avrg_cosp) write(*, '(a)') &
            ' Reading (co)spectra for ensemble averaging..'

        !> Convert start/end to timestamps
        call DateTimeToDateType(FCCsetup%SA%start_date, '00:00', saStartTimestamp)
        call DateTimeToDateType(FCCsetup%SA%end_date,   '23:59', saEndTimestamp)

        !> Read names of binned (co)spectra files
        call NumberOfFilesInDir(Dir%binned, '.csv', .false., '', &
            NumBinnedFiles, NumBinnedFilesNoRecurse)
        allocate(BinnedFileList(NumBinnedFiles))

        !> Create list of binned spectra file names
        call FileListByExt(Dir%binned, '.csv', .true., BinnedFilePrototype, &
            .true., .true., .true., BinnedFileList, size(BinnedFileList), &
            .true., ' ')

        !> Set file list in chronological order
        call FilesInChronologicalOrder(BinnedFileList, size(BinnedFileList), &
            binStartTimestamp, binEndTimestamp)

        !> Detect first and last binned files to be used for spectral \n
        !> assessment based on user's dates selection
        call tsExtractSubperiodIndexesFromFilelist(BinnedFileList, &
            size(BinnedFileList), saStartTimestamp, saEndTimestamp, &
            saStartTimestampIndx, saEndTimestampIndx)

        !> Some logging
        write(LogInteger, '(i8)') saEndTimestampIndx - saStartTimestampIndx + 1
        write(*, '(a)') '  Importing and sorting up to ' &
            // trim(adjustl(LogInteger)) // ' binned (co)spectra from files.. '

        !> Create an exponentially spaced frequency array in a range \n
        !> wide enough to accommodate any possible normalized frequency
        dkf(1) = 1d0 / (60d0 * 60d0 * 4d0)     !< 1 / (4 hours in seconds)
        dkf(ndkf + 1) = 200d0 / 2d0            !< Max normalized freq of 200 Hz
        do i = 2, ndkf
            dkf(i) = dkf(1) * dexp(dble(i - 1) &
                * (dlog(dkf(ndkf + 1)) - dlog(dkf(1))) / dble(ndkf))
        end do

        !> Open Ex file to keep it ready for reading
        !> and exit with error in case of problems opening the file
        open(uex, file = AuxFile%ex, status = 'old', iostat = open_status)
        if (open_status /= 0) call ExceptionHandler(60)
        !> Skip header in Ex file
        read(uex, *)

        !> Loop to import binned (co)spectra
        month = 0
        day   = 0
        nfit = 0
        fcount = saStartTimestampIndx - 1
        binned_loop: do
            !> Update file counter
            fcount = fcount + 1

            !> Normal exit instruction
            if (fcount > saEndTimestampIndx) exit binned_loop

            !> Read (co)spectra from file
            call ReadBinnedFile(BinnedFileList(fcount), BinSpec, BinCosp, &
                size(BinSpec), nbins, skip)
            if (skip) cycle binned_loop

            !> Show advancement
            if (day /= BinnedFileList(fcount)%timestamp%Day &
                .or. month /= BinnedFileList(fcount)%timestamp%Month) then
                month = BinnedFileList(fcount)%timestamp%Month
                day   = BinnedFileList(fcount)%timestamp%Day
                call DisplayProgress('daily', &
                    '  Importing binned spectra for ', &
                    BinnedFileList(fcount)%timestamp, 'yes')
            end if

            !> Retrieve ex information for current spectra
            call RetrieveExVarsByTimestamp(uex, &
                BinnedFileList(fcount)%timestamp, lEx, exEndReached, skip)
            if (exEndReached) exit binned_loop
            if (skip) cycle binned_loop

            !> Allocate variables that depend on nbins and initialize them
            if (.not. allocated(MeanBinSpec)) then
                allocate(MeanBinSpec(nbins, MaxGasClasses))
                MeanBinSpec = NullMeanSpec
            end if
            if (.not. allocated(MeanBinCosp)) then
                allocate(MeanBinCosp(nbins, MaxGasClasses))
                MeanBinCosp = NullMeanSpec
                allocate(FitUnstable(nbins * (saEndTimestampIndx - saStartTimestampIndx + 1)))
                FitUnstable = NullFitCosp
                allocate(FitStable  (nbins * (saEndTimestampIndx - saStartTimestampIndx + 1)))
                FitStable   = NullFitCosp
            end if

            !> Eliminate (co)spectra based on user-selected quality criteria
            call CospectraQAQC(BinSpec, BinCosp, size(BinSpec), lEx, &
                BinCospForStable, BinCospForUnstable, skip_spectra, skip_cospectra)

            !> Sort current spectra in relevant classes
            if (.not. skip_spectra) call SpectraSortingAndAveraging(lEx, BinSpec, size(BinSpec), nbins)

            !> Sort current cospectra in time-slot classes
            if (EddyProProj%out_avrg_cosp .and. .not. skip_cospectra) then

                !> Add current cospectra to dataset for regression
                call AddToCospectraFitDataset(lEx, BinCospForStable, BinCospForUnstable, &
                    size(BinCospForStable), nfit, size(nfit, 1), size(nfit, 2), &
                    nbins, FitStable, FitUnstable, size(FitStable))

                !> Sort current cospectra in time slot classes
                call CospectraSortingAndAveraging(BinCosp, size(BinCosp), lEx%time, nbins)
            end if
        end do binned_loop
        close(uex)
        write(*,'(a)') '  Done.'

        !> If cospectra were found for fitting, fit Massman model
        if (EddyProProj%out_avrg_cosp) then
            call FitCospectralModel(nfit, size(nfit, 1), size(nfit, 2), &
                FitStable, FitUnstable, size(FitStable))
        end if

        !> Ensemble average cospectra in stable and unstable cases
        call EnsembleCospectraByStability(nfit, size(nfit, 1), size(nfit, 2), FitStable, FitUnstable, size(FitStable))

        !> Normalize sums for obtaining mean spectra
        call NormalizeMeanSpectraCospectra(nbins)

        !> Detect which average spectra are available for each gas and each class
        call AvailableMeanSpectraCospectra(nbins)

        if (DoSpectralAssessment) then

            !> Determine low-pass TF cut-off frequencies, RH-sorted (H2O) and time-sorted (CO2/CH4/GAS4)
            call FitTFModels(nbins)

            !> Determine analytical relation fc/RH
            call FitRh2Fco()

            !> If necessary, calculate spectral correction factor models as from Ibrom et al. (2007)
            call CorrectionFactorModel(AuxFile%ex, NumExRecords)
        end if

        !> Write everything on output files
        call OutputSpectralAssessmentResults(DoSpectralAssessment, EddyProProj%out_avrg_cosp, nbins)
        write(*,'(a)')

    elseif (.not. FCCsetup%sa_onthefly .and. FCCsetup%SA%in_situ) then

        !> Enter here if an in-situ method was chosen, but results are available
        call ReadSpectralAssessmentFile()

    end if

    !************************************************************************************
    !************************************************************************************
    !*********** MAIN CYCLE ON RESULTS RECORDS RETRIEVED FROM ESSENTIALS FILE ***********
    !************************************************************************************
    !************************************************************************************

    !> Open Ex file to keep it ready for reading
    !> and exit with error in case of problems opening the file
    open(uex, file = AuxFile%ex, status = 'old', iostat = open_status)
    if (open_status /= 0) call ExceptionHandler(60)
    !> Skip first line for header
    read(uex, *)

    month = 0
    day   = 0
    InitializeOuputFiles = .true.
    ex_loop: do i = 1, NumExRecords

        !> Read record from essentials file
        call ReadExRecord('', uex, -1, lEx, ValidRecord, EndOfFileReached)

        !> Initialize presence of key variables for outputting results
        if (InitializeOuputFiles) fcc_var_present(u:GHGNumVar) = lEx%var_present(u:GHGNumVar)

        !> If end of file was reached, exit loop
        if (EndOfFileReached) exit ex_loop

        !> If invalid record was found, cycle loop
        if (.not. ValidRecord) cycle ex_loop

        !> Retrieve timestamp
        call DateTimeToDateType(lEx%date, lEx%time, CurrentTimestamp)

        !> If current timestamp is < start selected timestamp, cycle
        if (CurrentTimestamp < MasterTimeSeries(fxStartTimestampIndx)) cycle ex_loop

        !> If current timestamp is > end selected range, exit
        if (CurrentTimestamp > MasterTimeSeries(fxEndTimestampIndx)) exit ex_loop

        !> Show advancement
        call DateTimetoDOY(lEx%date, lEx%time, int_doy, float_doy)
        if (day /= CurrentTimestamp%day .or. month /= CurrentTimestamp%month) then
            month = CurrentTimestamp%month
            day   = CurrentTimestamp%day
            call DisplayProgress('daily', '  Calculating fluxes for ', CurrentTimestamp, 'yes')
        end if

        !> Band-pass spectral correction factors
        BPCF%of(:) = 1d0

        !> Create aux variables to pass to BandPassSpectralCorrections
        AuxInstrument = NullInstrument
        AuxInstrument(sonic) = lEx%instr(sonic)
        AuxInstrument(co2:gas4) = lEx%instr(ico2:igas4)
        if (.not. allocated(FullFileList)) allocate(FullFileList(1))

        !> Bad pass spectral correction factors
        call BandPassSpectralCorrections(lEx%instr(sonic)%height, lEx%disp_height, &
            lEx%var_present, lEx%WS, lEx%Ta, lEx%zL, lEx%ac_freq, nint(lEx%avrg_length), &
            lEx%det_meth, nint(lEx%det_timec), size(FullFileList), .false., AuxInstrument, &
            FullFileList, nrow_full, lEx, FCCsetup)

        !> Calculate fluxes at Level 1
        call Fluxes1(lEx)

        !> Calculate fluxes at Level 2 and Level 3
        call Fluxes23(lEx)

        !> Calculate footprint estimation
        call FootprintHandle(lEx%var(w), lEx%ustar, lEx%zL, lEx%WS, lEx%L, &
            lEx%instr(sonic)%height, lEx%disp_height, lEx%rough_length)

        !> Calculate quality flags
        StDiff%w_u    = nint(lEx%st_w_u)
        StDiff%w_ts   = nint(lEx%st_w_ts)
        StDiff%w_co2  = nint(lEx%st_w_co2)
        StDiff%w_h2o  = nint(lEx%st_w_h2o)
        StDiff%w_ch4  = nint(lEx%st_w_ch4)
        StDiff%w_gas4 = nint(lEx%st_w_gas4)
        DtDiff%u      = nint(lEx%dt_u)
        DtDiff%w      = nint(lEx%dt_w)
        DtDiff%ts     = nint(lEx%dt_ts)
        call QualityFlags(Flux2, StDiff, DtDiff, STFlg, DTFlg, QCFlag, .false.)

        !> Initialize output files
        if (InitializeOuputFiles) then
            call InitOutFiles(lEx)
            InitializeOuputFiles = .false.
        end if

        !> Write results on output file
        call WriteOutputFiles(lEx)

        !> Write AmeriFlux output if requested
        if (EddyProProj%out_amflux) call WriteAmeriFluxOutput(lEx)
    end do ex_loop
    close(uex)
    close(uflx)
    close(ughgeu)
    close(uaflx)
    close(umd)

    write(*,*)
    write(*,*)
    call sleep(1)

    !> Creating datasets from output files
    if (EddyProProj%make_dataset) then
        call CreateDatasetsCommon(MasterTimeSeries, size(MasterTimeSeries), &
            fxStartTimestampIndx, fxEndTimestampIndx, 'FCC')
    else
        call RenameTmpFilesCommon()
    end if

    !> Delete tmp folder if running in embedded mode
    if(EddyProProj%run_env == 'desktop') call system(trim(comm_rmdir) // ' "' // trim(adjustl(TmpDir)) // '"')

    write(*, '(a)') ''
    write(*, '(a)') ' ****************************************************'
    write(*, '(a)') ' Program EddyPro executed gracefully.'
    write(*, '(a)') ' Check results in the selected output directory.     '
    write(*, '(a)') ' ****************************************************'
    stop ''
end program EddyproFCC
