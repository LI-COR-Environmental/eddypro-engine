!***************************************************************************
! eddypro-fcc_main.f90
! --------------------
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
    integer :: del_status

    character(10) :: sDate, eDate
    character(5) :: sTime, eTime

    logical :: skip
    logical :: InitializeOuputFiles
    logical :: skip_spectra
    logical :: skip_cospectra
    logical :: ValidRecord
    logical :: EndOfFileReached
    logical :: exEndReached

    !> Derived type variables
    type(DateType) :: exStartTimestamp
    type(DateType) :: exEndTimestamp
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
    type(ExType) :: lEx

    !> Allocatable variabled
    type(DateType), allocatable :: exTimeSeries(:)
    type(DateType), allocatable :: MasterTimeSeries(:)
    type(FitSpectraType), allocatable :: FitUnstable(:)
    type(FitSpectraType), allocatable :: FitStable(:)

    !> External functions
    integer, external :: NumOfPeriods
    integer, external :: NumberOfFilesInSubperiod
    real(kind = dbl), external :: func

    include '../src_common/interfaces.inc'

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

    !> Add run-mode tag to Timestamp_FilePadding
    call TagRunMode()

    !> If running in embedded mode, override some settings
    if (EddyProProj%run_env == 'embedded') &
        call ConfigureForEmbedded('EddyPro-FCC')

    !> Preliminarily read essential files and retrieve a few information
    call InitExVars(exStartTimestamp, exEndTimestamp, &
        NumExRecords, NumValidExRecords)

    call ReadExRecord(AuxFile%ex, udf, 1, lEx, ValidRecord, EndOfFileReached)

    !> If no good records are found stop execution
    if (NumValidExRecords <= 0) call ExceptionHandler(61)

    !> Retrieve NumberOfPeriods and allocate exTimeSeries
    NumberOfPeriods = NumOfPeriods(exStartTimestamp, exEndTimestamp, DateStep)
    allocate(exTimeSeries(NumberOfPeriods + 1))
    call CreateTimeSeries(exStartTimestamp, exEndTimestamp, &
        DateStep, exTimeSeries, size(exTimeSeries), .true.)

    !> Define MasterTimeSeries for the period to be considered
    if (EddyProProj%subperiod) then
        call DateTimeToDateType(EddyProProj%start_date, &
            EddyProProj%start_time, SelectedStartTimestamp)
        call DateTimeToDateType(EddyProProj%end_date, &
            EddyProProj%end_time, SelectedEndTimestamp)
        NumberOfPeriods = NumOfPeriods(SelectedStartTimestamp, &
            SelectedEndTimestamp, DateStep)

        allocate(MasterTimeSeries(NumberOfPeriods + 1))
        call CreateTimeSeries(SelectedStartTimestamp, SelectedEndTimestamp, &
            DateStep, MasterTimeSeries, size(MasterTimeSeries), .false.)

        !> Verify at least partial overlap
        if (MasterTimeSeries(1) > exTimeSeries(size(exTimeSeries)) &
            .or. MasterTimeSeries(size(MasterTimeSeries)) < exTimeSeries(1)) &
            call ExceptionHandler(50)
        !> Set start/end indexes in MasterTimeSeries. Note that definition of
        !> fxStartTimestampIndx works, but isn't very nice. Should be improved.
        fxEndTimestampIndx = size(MasterTimeSeries)
        if (EddyProProj%make_dataset) then
            fxStartTimestampIndx = 1
        else
            fxStartTimestampIndx = 2
        end if
    else
        allocate(MasterTimeSeries(size(exTimeSeries)))
        MasterTimeSeries = exTimeSeries
        !> Set start/end indexes in MasterTimeSeries.
        fxStartTimestampIndx = 1
        fxEndTimestampIndx = size(MasterTimeSeries)
    end if

    !> Initialize import of full cospectra files
    if (FCCsetup%import_full_cospectra) then
        call NumberOfFilesInDir(Dir%full, '.csv', .false., '', &
            NumFullFiles, NumFullFilesNoRecurse)
        allocate(FullFileList(NumFullFiles))

        !> Read names of full cospectra files
        call FileListByExt(Dir%full, '.csv', .true., .false., &
            FullFilePrototype, .true., .true., .true., FullFilelist, &
            size(FullFilelist), .true., '')

        !> Retrieve length of full cospectra for later allocation
        call FullCospectraLength(FullFilelist(1)%path, nrow_full)
    end if

    !****************************************************************
    !****************************************************************
    !*********** SPECTRAL ASSESSMENT SECTION IF REQUESTED ***********
    !****************************************************************
    !****************************************************************

    !> If spectral analysis must be performed or co-spectral
    !> outputs are requested, start loop on cospectra files
    if (FCCsetup%pass_thru_spectral_assessment) then
        if (FCCsetup%do_spectral_assessment) write(*, '(a)') &
            ' Starting "spectral assessment" session..'
        if (EddyProProj%out_avrg_cosp .or. EddyProProj%out_avrg_spec) then
            write(*, '(a)') ' Reading (co)spectra from:'
            write(*, '(a)') '  ' // trim(adjustl(Dir%binned))
            write(*, '(a)') ''
        end if

        !> Convert start/end to timestamps
        call DateTimeToDateType(FCCsetup%SA%start_date, &
            FCCsetup%SA%start_time, saStartTimestamp)
        call DateTimeToDateType(FCCsetup%SA%end_date, &
            FCCsetup%SA%end_time, saEndTimestamp)

        !> Spectra hold the filestamp of the end of the period, so increase
        !> start timestamp by DateStep
        saStartTimestamp = saStartTimestamp + DateStep

        !> Read names of binned (co)spectra files
        call NumberOfFilesInDir(Dir%binned, '.csv', .false., '', &
            NumBinnedFiles, NumBinnedFilesNoRecurse)

        if(NumBinnedFiles <= 0) then
            !> Exit loop and set spectral correction to Moncrieff.
            !> Also, cannot create any ensemble spectral output
            EddyProProj%out_avrg_cosp = .false.
            EddyProProj%out_avrg_spec = .false.
            FCCsetup%do_spectral_assessment = .false.
            if (FCCsetup%SA%in_situ) then
                EddyProProj%hf_meth = 'moncrieff_97'
                FCCsetup%SA%in_situ = .false.
            end if
            call ExceptionHandler(89)
            goto 100
        end if
        allocate(BinnedFileList(NumBinnedFiles))

        !> Create list of binned spectra file names
        call FileListByExt(Dir%binned, '.csv', .true., .false., &
            BinnedFilePrototype, .true., .true., .true., &
            BinnedFileList, size(BinnedFileList), .true., ' ')

        !> Set file list in chronological order
        call FilesInChronologicalOrder(BinnedFileList, size(BinnedFileList), &
            binStartTimestamp, binEndTimestamp, ' ')

        !> Detect first and last binned files to be used for spectral \n
        !> assessment based on user's dates selection
        call tsExtractSubperiodIndexesFromFilelist(BinnedFileList, &
            size(BinnedFileList), saStartTimestamp, saEndTimestamp, &
            saStartTimestampIndx, saEndTimestampIndx)

        if(saStartTimestampIndx <= 0 .or. saEndTimestampIndx <= 0) then
            !> Exit loop and set spectral correction to Moncrieff.
            !> Also, cannot create any ensemble spectral output
            EddyProProj%out_avrg_cosp = .false.
            EddyProProj%out_avrg_spec = .false.
            FCCsetup%do_spectral_assessment = .false.
            if (FCCsetup%SA%in_situ) then
                EddyProProj%hf_meth = 'moncrieff_97'
                FCCsetup%SA%in_situ = .false.
            end if
            call ExceptionHandler(90)
            goto 100
        end if

        !> Some logging
        call DateTypeToDateTime(binStartTimestamp - DateStep, sDate, sTime)
        call DateTypeToDateTime(binEndTimestamp, eDate, eTime)
        write(*, '(a)') ''
        write(*, '(a)') '  Period covered by available binned (co)spectra files:'
        write(*, '(a)') '   Start: ' // sDate // ' ' // sTime
        write(*, '(a)') '   End:   ' // eDate // ' ' // eTime

        if (FCCsetup%SA%subperiod) then
            call DateTypeToDateTime(saStartTimestamp, sDate, sTime)
            call DateTypeToDateTime(saEndTimestamp + Datetype(0, 0, 0, 0, 1), eDate, eTime)
            write(*, '(a)') ''
            write(*, '(a)') '  Selected (co)spectra sub-period:'
            write(*, '(a)') '   Start: ' // sDate // ' ' // sTime
            write(*, '(a)') '   End:   ' // eDate // ' ' // eTime
        end if

        write(*, '(a)') ''
        write(LogInteger, '(i8)') saEndTimestampIndx - saStartTimestampIndx + 1
        write(*, '(a)') '  Importing, sorting and ensemble-averaging up to ' &
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

            !> Allocate variables that depend upon nbins and initialize them
            if (.not. allocated(MeanBinSpec)) then
                allocate(MeanBinSpec(nbins, MaxGasClasses))
                allocate(dMeanBinSpec(nbins, MaxGasClasses))
                MeanBinSpec = NullMeanSpec
                dMeanBinSpec = NullMeanSpec
            end if
            if (.not. allocated(MeanBinCosp)) then
                allocate(MeanBinCosp(nbins, MaxGasClasses))
                MeanBinCosp = NullMeanSpec
                allocate(FitUnstable(nbins * &
                    (saEndTimestampIndx - saStartTimestampIndx + 1)))
                FitUnstable = NullFitCosp
                allocate(FitStable  (nbins * &
                    (saEndTimestampIndx - saStartTimestampIndx + 1)))
                FitStable   = NullFitCosp
            end if

            !> Eliminate (co)spectra based on user-selected quality criteria
            call CospectraQAQC(BinSpec, BinCosp, size(BinSpec), lEx, &
                BinCospForStable, BinCospForUnstable, &
                skip_spectra, skip_cospectra)

            !> Sort current spectra in relevant classes
            if (.not. skip_spectra) &
                call SpectraSortingAndAveraging(lEx, BinSpec, &
                    size(BinSpec), nbins)

            !> Sort current cospectra in time-slot classes
            if (EddyProProj%out_avrg_cosp .and. .not. skip_cospectra) then

                !> Add current cospectra to dataset for regression
                call AddToCospectraFitDataset(lEx, BinCospForStable, &
                    BinCospForUnstable, size(BinCospForStable), nfit, &
                    size(nfit, 1), size(nfit, 2), nbins, FitStable, &
                    FitUnstable, size(FitStable))

                !> Sort current cospectra in time slot classes
                call CospectraSortingAndAveraging(BinCosp, size(BinCosp), &
                    lEx%end_time, nbins)
            end if
        end do binned_loop
        close(uex)
        write(*,'(a)') '  Done.'

        if (EddyProProj%out_avrg_cosp) then
            !> If cospectra were found for fitting, fit Massman model
            call FitCospectralModel(nfit, size(nfit, 1), size(nfit, 2), &
                FitStable, FitUnstable, size(FitStable))

            !> Ensemble average cospectra in stable and unstable stratifications
            call EnsembleCospectraByStability(nfit, size(nfit, 1), &
                size(nfit, 2), FitStable, FitUnstable, size(FitStable))

        end if

        !> Normalize sums for obtaining mean spectra
        call NormalizeMeanSpectraCospectra(nbins)

        !> Detect which average spectra are available
        !> for each gas and each class
        call AvailableMeanSpectraCospectra(nbins)

        !> Determine low-pass TF cut-off frequencies, RH-sorted (H2O)
        !> and time-sorted (CO2/CH4/GAS4)
        call FitTFModels(nbins, FCCsetup%do_spectral_assessment)

        !> Spectral attenuation assessment
        if (FCCsetup%do_spectral_assessment) then

            !> Determine analytical relation fc/RH
            call FitRh2Fco()

            !> If necessary, calculate spectral correction factor models
            !> as from Ibrom et al. (2007)
            call CorrectionFactorModel(AuxFile%ex, NumExRecords)
        else
            !> If an in-situ method was chosen, and spectral
            !> assessment file is available, read file
            if (FCCsetup%SA%in_situ) call ReadSpectralAssessmentFile()
        end if

        !> Write number of imported spectra and cospectra on stdout
        if (EddyProProj%out_avrg_spec .or. FCCsetup%do_spectral_assessment) &
            call ReportImportedSpectra(nbins)

        !> Write everything on output files
        call OutputSpectralAssessmentResults(nbins)
        write(*,'(a)')


    else
        !> If an in-situ method was chosen, and spectral
        !> assessment file is available, read file
        if (FCCsetup%SA%in_situ) call ReadSpectralAssessmentFile()
    end if

100 continue

    !***************************************************************************
    !***************************************************************************
    !****** MAIN CYCLE ON RESULTS RECORDS RETRIEVED FROM ESSENTIALS FILE *******
    !***************************************************************************
    !***************************************************************************

    !> Establish present variables
    ! call EstablishPresentVariables()

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
        ! if (InitializeOuputFiles) &
        !     fcc_var_present(u:GHGNumVar) = lEx%var_present(u:GHGNumVar)

        !> If end of file was reached, exit loop
        if (EndOfFileReached) exit ex_loop

        !> If invalid record was found, cycle loop
        if (.not. ValidRecord) cycle ex_loop

        !> Retrieve timestamp
        call DateTimeToDateType(lEx%end_date, lEx%end_time, CurrentTimestamp)

        !> If current timestamp is < start selected timestamp, cycle
        if (CurrentTimestamp < MasterTimeSeries(fxStartTimestampIndx)) &
            cycle ex_loop

        !> If current timestamp is > end selected range, exit
        if (CurrentTimestamp > MasterTimeSeries(fxEndTimestampIndx)) &
            exit ex_loop

        !> Show advancement
        call DateTimetoDOY(lEx%end_date, lEx%end_time, int_doy, float_doy)
        if (day /= CurrentTimestamp%day &
            .or. month /= CurrentTimestamp%month) then
            month = CurrentTimestamp%month
            day   = CurrentTimestamp%day
            call DisplayProgress('daily', '  Calculating fluxes for ', &
                CurrentTimestamp, 'yes')
        end if

        !> Band-pass spectral correction factors
        BPCF%of(:) = 1d0

        !> Create aux variables to pass to BandPassSpectralCorrections
        AuxInstrument = NullInstrument
        AuxInstrument(sonic) = lEx%instr(sonic)
        AuxInstrument(co2:gas4) = lEx%instr(ico2:igas4)
        if (.not. allocated(FullFileList)) allocate(FullFileList(1))

        !> Bad pass spectral correction factors
        call BandPassSpectralCorrections(lEx%instr(sonic)%height, &
            lEx%disp_height, lEx%var_present, lEx%WS, lEx%Ta, lEx%Flux0%zL, &
            lEx%ac_freq, nint(lEx%avrg_length), lEx%logger_swver, &
            lEx%det_meth, nint(lEx%det_timec), .false., AuxInstrument, &
            size(FullFileList), FullFileList, nrow_full, lEx, FCCsetup)

        !> Calculate fluxes at Level 1
        call Fluxes1(lEx)

        !> Calculate fluxes at Level 2 and Level 3
        call Fluxes23(lEx)

        !> Calculate footprint estimation   
        call FootprintHandle(lEx%var(w), lEx%ustar, lEx%zL, lEx%WS, lEx%L, &
            lEx%instr(sonic)%height, lEx%disp_height, lEx%rough_length)

        !> Calculate quality flags
        StDiff%w_u    = nint(lEx%TAU_SS)
        StDiff%w_ts   = nint(lEx%H_SS)
        StDiff%w_co2  = nint(lEx%FC_SS)
        StDiff%w_h2o  = nint(lEx%FH2O_SS)
        StDiff%w_ch4  = nint(lEx%FCH4_SS)
        StDiff%w_gas4 = nint(lEx%FGS4_SS)
        DtDiff%u      = nint(lEx%U_ITC)
        DtDiff%w      = nint(lEx%W_ITC)
        DtDiff%ts     = nint(lEx%TS_ITC)
        call QualityFlags(Flux2, StDiff, DtDiff, STFlg, DTFlg, QCFlag, .false.)

        !> Initialize output files
        if (InitializeOuputFiles) then
            call InitOutFiles(lEx)
            InitializeOuputFiles = .false.
        end if

        !>Write out full output file
        if (EddyProProj%out_full) call WriteOutFull(lEx)
        if (EddyProProj%out_md) call WriteOutMetadata(lEx)
        if (EddyProProj%out_fluxnet) call WriteOutFluxnet(lEx)

    end do ex_loop
    close(uex)
    close(uflx)
    close(ufnet_e)
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
    if(EddyProProj%run_env == 'desktop') &
        del_status = system(trim(comm_rmdir) // ' "' &
        // trim(adjustl(TmpDir)) // '"')

    write(*, '(a)') ''
    write(*, '(a)') ' ****************************************************'
    write(*, '(a)') ' Program EddyPro executed gracefully.'
    write(*, '(a)') ' Check results in the selected output directory.     '
    write(*, '(a)') ' ****************************************************'
    stop ''
end program EddyproFCC
