!***************************************************************************
! output_spectral_assessment_results.f90
! --------------------------------------
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
! \brief       Write results of spectral assessment on output file
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine OutputSpectralAssessmentResults(nbins)
    use m_fx_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nbins
    !> local variables
    integer, external :: CreateDir
    integer :: i
    integer :: pick
    integer :: goodj
    integer :: cls
    integer :: gas
    integer :: month
    integer :: open_status
    integer :: mkdir_status
    real(kind = dbl), external :: func
    real(kind = dbl), external :: kaimal
    character(128) :: Filename
    character(PathLen) :: FilePath
    character(PathLen) :: SpecDir
    character(LongOutstringLen) :: dataline
    character(DatumLen) :: datum
    logical :: proceed
    include '../src_common/interfaces_1.inc'

    !> Create output directory
    mkdir_status = CreateDir(Dir%main_out(1:len_trim(Dir%main_out)))
    SpecDir = Dir%main_out(1:len_trim(Dir%main_out)) // SubDirSpecAn // slash
    mkdir_status = CreateDir('"' // SpecDir(1:len_trim(SpecDir)) // '"')

    !> Select one frequency vector which does not contain only -9999.
    !> If none, it means all spectra are -9999. Basically, spectral
    !> assessment failed.
    goodj = ierror
    ol: do cls = RH10, RH90
        il: do i = 1, nbins - 1
                if (MeanBinSpec(i, cls)%fn(h2o) /= error) then
                    goodj = cls
                    exit ol
                end if
        end do il
    end do ol

    !> SPECTRAL ASSESSMENT
    if (FCCsetup%do_spectral_assessment) then

        if (goodj == ierror) then
            call ExceptionHandler(76)
        else
            write(*,'(a)') ' Writing spectral assessment results on file.. '

            !> Transfer function parameters
            Filename = EddyProProj%id(1:len_trim(EddyProProj%id)) // SA_FilePadding  &
                // Timestamp_FilePadding // TxtExt
            FilePath = SpecDir(1:len_trim(SpecDir)) // Filename(1:len_trim(Filename))
            open(udf, file = FilePath, iostat = open_status)
            if (open_status /= 0) call ExceptionHandler(64)

            write(udf,'(a)') 'Transfer_function_parameters_(TFP)_for_&
                &IIR-shaped_filter_(see_Ibrom_et_al._2007_AFM).'
            write(udf,'(a)') 'fc:_IIR_cut-off_frequency'
            write(udf,'(a)') 'Fn:_normalization_parameter'
            write(udf,'(a)') 'Water_vapour_TFP_are_calculated_for_9_RH_classes.'
            write(udf,'(a)') 'Other_gases_TFP_are_calculated_on_a_monthly_base_&
                &(currently_all_months_together_).'
            write(udf,'(a)') '-----------------------------------------------------&
                &-----------------------------'
            write(udf,'(a)') 'Water vapour TFP              Fn          fc    numerosity'
            write(udf,'(a, 2(f11.5,1x), i13)') 'RH class   5 - 15% = ', &
                RegPar(h2o, RH10)%Fn, RegPar(h2o, RH10)%fc, &
                MeanBinSpec(nbins/2, RH10)%cnt(h2o)
            write(udf,'(a, 2(f11.5,1x), i13)') 'RH class  15 - 25% = ', &
                RegPar(h2o, RH20)%Fn, RegPar(h2o, RH20)%fc, &
                MeanBinSpec(nbins/2, RH20)%cnt(h2o)
            write(udf,'(a, 2(f11.5,1x), i13)') 'RH class  25 - 35% = ', &
                RegPar(h2o, RH30)%Fn, RegPar(h2o, RH30)%fc, &
                MeanBinSpec(nbins/2, RH30)%cnt(h2o)
            write(udf,'(a, 2(f11.5,1x), i13)') 'RH class  35 - 45% = ', &
                RegPar(h2o, RH40)%Fn, RegPar(h2o, RH40)%fc, &
                MeanBinSpec(nbins/2, RH40)%cnt(h2o)
            write(udf,'(a, 2(f11.5,1x), i13)') 'RH class  45 - 55% = ', &
                RegPar(h2o, RH50)%Fn, RegPar(h2o, RH50)%fc, &
                MeanBinSpec(nbins/2, RH50)%cnt(h2o)
            write(udf,'(a, 2(f11.5,1x), i13)') 'RH class  55 - 65% = ', &
                RegPar(h2o, RH60)%Fn, RegPar(h2o, RH60)%fc, &
                MeanBinSpec(nbins/2, RH60)%cnt(h2o)
            write(udf,'(a, 2(f11.5,1x), i13)') 'RH class  65 - 75% = ', &
                RegPar(h2o, RH70)%Fn, RegPar(h2o, RH70)%fc, &
                MeanBinSpec(nbins/2, RH70)%cnt(h2o)
            write(udf,'(a, 2(f11.5,1x), i13)') 'RH class  75 - 85% = ', &
                RegPar(h2o, RH80)%Fn, RegPar(h2o, RH80)%fc, &
                MeanBinSpec(nbins/2, RH80)%cnt(h2o)
            write(udf,'(a, 2(f11.5,1x), i13)') 'RH class  85 - 95% = ', &
                RegPar(h2o, RH90)%Fn, RegPar(h2o, RH90)%fc, &
                MeanBinSpec(nbins/2, RH90)%cnt(h2o)
            write(udf,'(a)') ''

            do gas = co2, gas4
                if (gas == co2)  write(udf,'(a)') 'CO2            TFP            &
                    &Fn          fc'
                if (gas == h2o)  cycle
                if (gas == ch4)  write(udf,'(a)') 'CH4            TFP            &
                    &Fn          fc'
                if (gas == gas4) write(udf,'(a)') g4lab(1:g4l) // '           &
                    &TFP            Fn          fc'

                if (FCCsetup%SA%class(gas, JAN) /= 0) then
                    write(udf,'(a, 2(f11.5,1x))') 'January            = ', &
                        RegPar(gas, FCCsetup%SA%class(gas, JAN))%Fn, &
                        RegPar(gas, FCCsetup%SA%class(gas, JAN))%fc
                else
                    write(udf,'(a, 2(f11.5,1x))') 'January            = ', error, error
                end if
                if (FCCsetup%SA%class(gas, FEB) /= 0) then
                    write(udf,'(a, 2(f11.5,1x))') 'February           = ', &
                        RegPar(gas, FCCsetup%SA%class(gas, FEB))%Fn, &
                        RegPar(gas, FCCsetup%SA%class(gas, FEB))%fc
                else
                    write(udf,'(a, 2(f11.5,1x))') 'February           = ', error, error
                end if
                if (FCCsetup%SA%class(gas, MAR) /= 0) then
                    write(udf,'(a, 2(f11.5,1x))') 'March              = ', &
                        RegPar(gas, FCCsetup%SA%class(gas, MAR))%Fn, &
                        RegPar(gas, FCCsetup%SA%class(gas, MAR))%fc
                else
                    write(udf,'(a, 2(f11.5,1x))') 'March              = ', error, error
                end if
                if (FCCsetup%SA%class(gas, APR) /= 0) then
                    write(udf,'(a, 2(f11.5,1x))') 'April              = ', &
                        RegPar(gas, FCCsetup%SA%class(gas, APR))%Fn, &
                        RegPar(gas, FCCsetup%SA%class(gas, APR))%fc
                else
                    write(udf,'(a, 2(f11.5,1x))') 'April              = ', error, error
                end if
                if (FCCsetup%SA%class(gas, MAY) /= 0) then
                    write(udf,'(a, 2(f11.5,1x))') 'May                = ', &
                        RegPar(gas, FCCsetup%SA%class(gas, MAY))%Fn, &
                        RegPar(gas, FCCsetup%SA%class(gas, MAY))%fc
                else
                    write(udf,'(a, 2(f11.5,1x))') 'May                = ', error, error
                end if
                if (FCCsetup%SA%class(gas, JUN) /= 0) then
                    write(udf,'(a, 2(f11.5,1x))') 'June               = ', &
                        RegPar(gas, FCCsetup%SA%class(gas, JUN))%Fn, &
                        RegPar(gas, FCCsetup%SA%class(gas, JUN))%fc
                else
                    write(udf,'(a, 2(f11.5,1x))') 'June               = ', error, error
                end if
                if (FCCsetup%SA%class(gas, JUL) /= 0) then
                    write(udf,'(a, 2(f11.5,1x))') 'July               = ', &
                        RegPar(gas, FCCsetup%SA%class(gas, JUL))%Fn, &
                        RegPar(gas, FCCsetup%SA%class(gas, JUL))%fc
                else
                    write(udf,'(a, 2(f11.5,1x))') 'July               = ', error, error
                end if
                if (FCCsetup%SA%class(gas, AUG) /= 0) then
                    write(udf,'(a, 2(f11.5,1x))') 'August             = ', &
                        RegPar(gas, FCCsetup%SA%class(gas, AUG))%Fn, &
                        RegPar(gas, FCCsetup%SA%class(gas, AUG))%fc
                else
                    write(udf,'(a, 2(f11.5,1x))') 'August             = ', error, error
                end if
                if (FCCsetup%SA%class(gas, SEP) /= 0) then
                    write(udf,'(a, 2(f11.5,1x))') 'September          = ', &
                        RegPar(gas, FCCsetup%SA%class(gas, SEP))%Fn, &
                        RegPar(gas, FCCsetup%SA%class(gas, SEP))%fc
                else
                    write(udf,'(a, 2(f11.5,1x))') 'September          = ', error, error
                end if
                if (FCCsetup%SA%class(gas, OCT) /= 0) then
                    write(udf,'(a, 2(f11.5,1x))') 'October            = ', &
                        RegPar(gas, FCCsetup%SA%class(gas, OCT))%Fn, &
                        RegPar(gas, FCCsetup%SA%class(gas, OCT))%fc
                else
                    write(udf,'(a, 2(f11.5,1x))') 'October            = ', error, error
                end if
                if (FCCsetup%SA%class(gas, NOV) /= 0) then
                    write(udf,'(a, 2(f11.5,1x))') 'November           = ', &
                        RegPar(gas, FCCsetup%SA%class(gas, NOV))%Fn, &
                        RegPar(gas, FCCsetup%SA%class(gas, NOV))%fc
                else
                    write(udf,'(a, 2(f11.5,1x))') 'November           = ', error, error
                end if
                if (FCCsetup%SA%class(gas, DEC) /= 0) then
                    write(udf,'(a, 2(f11.5,1x))') 'December           = ', &
                        RegPar(gas, FCCsetup%SA%class(gas, DEC))%Fn, &
                        RegPar(gas, FCCsetup%SA%class(gas, DEC))%fc
                else
                    write(udf,'(a, 2(f11.5,1x))') 'December           = ', error, error
                end if
                write(udf, *)
            end do

            !> Exponential fit f_co vs. RH
            write(udf,'(a)') 'RH/fc_exponential_fit_parameters_for_water_vapour&
                &_spectral_corrections'
            write(udf,'(a)') '-----------------------------------'
            write(udf,'(a)') '         exp1         exp2         exp3'
            write(udf,'(3(f13.6))') RegPar(dum, dum)%e1, &
                RegPar(dum, dum)%e2, RegPar(dum, dum)%e3
            write(udf,'(a)') ''
            write(udf,'(a)') ''

            write(udf, '(a)') 'High-pass_correction_factor_model_parameters'
            write(udf, '(a)') 'Model: CF = [c1 * u / (c2 + f_co) + 1] after_&
                &Ibrom_et_al_(2007_AFM)'
            write(udf, '(a)') '---------------------------------------------&
                &----------------------'
            write(udf, '(a)') '                   c1          c2'
            write(udf, '(a, 2(f11.7,1x))') 'unstable = ',UnPar(1), UnPar(2)
            write(udf, '(a, 2(f11.7,1x))') 'stable   = ',StPar(1), StPar(2)
            close(udf)
        end if
    end if

    !> ENSEMBLE AVERAGED SPECTRA
    if (FCCsetup%do_spectral_assessment .or. EddyProProj%out_avrg_spec) then
        !> Average H2O spectra, sorted in RH classes, and predicted spectra
        !> (RHS of eq. 6 in Ibrom et al. 2007, AFM)
        if (goodj == ierror) then
            call ExceptionHandler(77)
        else
            !> Initialize file
            Filename = EddyProProj%id(1:len_trim(EddyProProj%id)) &
                // H2OAvrg_FilePadding // Timestamp_FilePadding // CsvExt
            FilePath = SpecDir(1:len_trim(SpecDir)) // Filename(1:len_trim(Filename))
            open(udf, file = FilePath, iostat = open_status)
            if (open_status /= 0) call ExceptionHandler(64)

            write(udf,'(a)') 'Binned_average_and_predicted_H2O_spectra_sorted_by_RH-class.'
            write(udf,'(a)') ',RH=0.1,,,,RH=0.2,,,,RH=0.3,,,,RH=0.4,,,,RH=0.5,,,,RH=0.6&
                &,,,,RH=0.7,,,,RH=0.8,,,,RH=0.9'

            !> Add number of spectra per class
            dataline = ''
            call AddDatum(dataline, '', separator)
            do cls = RH10, RH90
                call WriteDatumInt(MeanBinSpec(1, cls)%cnt(h2o), datum, &
                    EddyProProj%err_label)
                call AddDatum(dataline, 'n_=_' // datum(1:len_trim(datum)), &
                    separator)
                call AddDatum(dataline, '', separator)
                call AddDatum(dataline, '', separator)
                call AddDatum(dataline, '', separator)
            end do
            write(udf,'(a)') dataline(1:len_trim(dataline) - 1)
            write(udf,'(a)') 'nat_freq,avrg_sp(T),avrg_sp(h2o),denoised_avrg_sp(h2o),pred_sp(h2o)&
                            &,avrg_sp(T),avrg_sp(h2o),denoised_avrg_sp(h2o),pred_sp(h2o)&
                            &,avrg_sp(T),avrg_sp(h2o),denoised_avrg_sp(h2o),pred_sp(h2o)&
                            &,avrg_sp(T),avrg_sp(h2o),denoised_avrg_sp(h2o),pred_sp(h2o)&
                            &,avrg_sp(T),avrg_sp(h2o),denoised_avrg_sp(h2o),pred_sp(h2o)&
                            &,avrg_sp(T),avrg_sp(h2o),denoised_avrg_sp(h2o),pred_sp(h2o)&
                            &,avrg_sp(T),avrg_sp(h2o),denoised_avrg_sp(h2o),pred_sp(h2o)&
                            &,avrg_sp(T),avrg_sp(h2o),denoised_avrg_sp(h2o),pred_sp(h2o)&
                            &,avrg_sp(T),avrg_sp(h2o),denoised_avrg_sp(h2o),pred_sp(h2o)'

            do i = 1, nbins - 1
                call clearstr(dataline)
                if (MeanBinSpec(i, goodj)%fn(h2o) /= error) then
                    call WriteDatumFloat(MeanBinSpec(i, goodj)%fn(h2o), &
                        datum, EddyProProj%err_label)
                    call AddDatum(dataline, datum, separator)
                    do cls = RH10, RH90
                        if (MeanBinSpecAvailable(cls, h2o)) then
                            !> Natural frequency
                            call WriteDatumFloat(MeanBinSpec(i, goodj)%fn(h2o) &
                                * MeanBinSpec(i, cls)%ts(h2o), datum, &
                                EddyProProj%err_label)
                            call AddDatum(dataline, datum, separator)
                            !> Ensemble averaged spectrum
                            call WriteDatumFloat(MeanBinSpec(i, goodj)%fn(h2o) &
                                * MeanBinSpec(i, cls)%of(h2o), datum, &
                                EddyProProj%err_label)
                            call AddDatum(dataline, datum, separator)
                            !> Denoised ensemble averaged spectrum
                            call WriteDatumFloat(MeanBinSpec(i, goodj)%fn(h2o) &
                                * dMeanBinSpec(i, cls)%of(h2o), datum, &
                                EddyProProj%err_label)
                            call AddDatum(dataline, datum, separator)
                            !> Modelled spectrum
                            call WriteDatumFloat(RegPar(h2o, cls)%Fn &
                                * (1d0 / (1d0 + (MeanBinSpec(i, goodj)%fn(h2o) &
                                / RegPar(h2o, cls)%fc)**2 )) &
                                * MeanBinSpec(i, cls)%ts(h2o) &
                                * MeanBinSpec(i, goodj)%fn(h2o), datum, &
                                EddyProProj%err_label)
                            call AddDatum(dataline, datum, separator)

                        else
                            call AddDatum(dataline, &
                                trim(adjustl(EddyProProj%err_label)), separator)
                            call AddDatum(dataline, &
                                trim(adjustl(EddyProProj%err_label)), separator)
                            call AddDatum(dataline, &
                                trim(adjustl(EddyProProj%err_label)), separator)
                            call AddDatum(dataline, &
                                trim(adjustl(EddyProProj%err_label)), separator)
                        end if
                    end do
                    write(udf, '(a)') dataline(1:len_trim(dataline) - 1)
                end if
            end do
            close(udf)

            !> =====================================================================
            !> Average CO2/CH4/GAS4 spectra, sorted by month, and predicted spectra (RHS of
            !> eq. 6 in Ibrom et al. (2007, AFM)

            !> Select one frequency vector which does not contain only -9999.
            !> If none, it means all spectra are -9999, so -9999 is written.
            !> Basically, run failed.
            call GetFnIndex(MeanBinSpec, &
                size(MeanBinSpec, 1), size(MeanBinSpec, 2), goodj, pick)

            !> Write output file if valid goodj and pick were found
            if (goodj > 0 .and. goodj < MaxGasClasses &
                .and. pick > 0 .and. pick < gas4) then
                Filename = EddyProProj%id(1:len_trim(EddyProProj%id)) // PASGAS_Avrg_FilePadding  &
                    // Timestamp_FilePadding // CsvExt
                FilePath = SpecDir(1:len_trim(SpecDir)) // Filename(1:len_trim(Filename))
                open(udf, file = FilePath, iostat = open_status)
                if (open_status /= 0) call ExceptionHandler(64)
                write(udf,'(a)') 'Binned_average_and_predicted_spectra_for_passive_gases'
                !> Add number of spectra per class
                dataline = ''
                call AddDatum(dataline, '', separator)
                do month = JAN, JAN
                    do gas = co2, gas4
                        if (gas /= h2o) then
                            if (FCCsetup%SA%class(gas, month) /= 0) then
                                call WriteDatumInt(MeanBinSpec(1, FCCsetup%SA%class(gas, month))%cnt(gas) &
                                    , datum, EddyProProj%err_label)
                            else
                                call WriteDatumInt(0, datum, EddyProProj%err_label)
                            end if
                            call AddDatum(dataline, 'n_=_' // datum(1:len_trim(datum)), separator)
                            call AddDatum(dataline, '', separator)
                            call AddDatum(dataline, '', separator)
                            call AddDatum(dataline, '', separator)
                        end if
                    end do
                end do
                write(udf,'(a)') dataline(1:len_trim(dataline) - 1)
                write(udf,'(a)') 'nat_freq,&
                    &avrg_sp(T),avrg_sp(co2),denoised_avrg_sp(co2),pred_sp(co2),&
                    &avrg_sp(T),avrg_sp(ch4),denoised_avrg_sp(ch4),pred_sp(ch4),&
                    &avrg_sp(T),avrg_sp(' // g4lab(1:g4l) // '),denoised_avrg_sp(' // g4lab(1:g4l) // '),&
                    &pred_sp(' // g4lab(1:g4l) // ')' !,&

                do i = 1, nbins - 1
                    call clearstr(dataline)
                    if (MeanBinSpec(i, goodj)%fn(pick) /= error) then
                        call WriteDatumFloat(MeanBinSpec(i, goodj)%fn(pick), datum, EddyProProj%err_label)
                        call AddDatum(dataline, datum, separator)
                        do month = JAN, JAN
                            do gas = co2, gas4
                                if (gas == h2o) cycle
                                if (FCCsetup%SA%class(gas, month) /= 0) then
                                    if (MeanBinSpecAvailable(FCCsetup%SA%class(gas, month), gas))then
                                        !> Natural frequency
                                        call WriteDatumFloat(MeanBinSpec(i, goodj)%fn(pick) * MeanBinSpec(i, &
                                            FCCsetup%SA%class(gas, month))%ts(gas), datum, EddyProProj%err_label)
                                        call AddDatum(dataline, datum, separator)
                                        !> Ensemble averaged spectrum
                                        call WriteDatumFloat(MeanBinSpec(i, goodj)%fn(pick) * MeanBinSpec(i, &
                                            FCCsetup%SA%class(gas, month))%of(gas), datum, EddyProProj%err_label)
                                        call AddDatum(dataline, datum, separator)
                                        !> Denoised ensemble averaged spectrum
                                        call WriteDatumFloat(MeanBinSpec(i, goodj)%fn(pick) * dMeanBinSpec(i, &
                                            FCCsetup%SA%class(gas, month))%of(gas), datum, EddyProProj%err_label)
                                        call AddDatum(dataline, datum, separator)
                                        !> Modelled spectrum
                                        call WriteDatumFloat(RegPar(gas, FCCsetup%SA%class(gas, month))%fn &
                                            * (1d0 / (1d0 + (MeanBinSpec(i, goodj)%fn(pick) &
                                            / RegPar(gas, FCCsetup%SA%class(gas, month))%fc)**2 )) &
                                            * MeanBinSpec(i, FCCsetup%SA%class(gas, month))%ts(gas) &
                                            * MeanBinSpec(i, goodj)%fn(pick), datum, EddyProProj%err_label)
                                        call AddDatum(dataline, datum, separator)

                                    else
                                        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
                                        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
                                        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
                                        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
                                    end if
                                else
                                    call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
                                    call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
                                    call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
                                    call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
                                end if
                            end do
                        end do
                        write(udf, '(a)') dataline(1:len_trim(dataline) - 1)
                    end if
                end do
                close(udf)
            else
                call ExceptionHandler(92)
            end if
        end if
    end if

    !> ENSEMBLE AVERAGED COSPECTRA
    if (EddyProProj%out_avrg_cosp) then

        !> =====================================================================
        !> Ensemble cospectra by time of day
        !> Select one frequency vector which does not contain only -9999.
        !> If none, it means all spectra are -9999, so -9999 is written.
        !> Basically, run failed.
        call GetFnIndex(MeanBinCosp, &
            size(MeanBinCosp, 1), size(MeanBinCosp, 2), goodj, pick)


        if (goodj == ierror .or. pick == ierror) then
            call ExceptionHandler(75)
        else
            Filename = trim(adjustl(EddyProProj%id)) &
                // Cosp_FilePadding // Timestamp_FilePadding // CsvExt
            FilePath = trim(adjustl(SpecDir)) // trim(adjustl(Filename))
            open(udf, file = FilePath, iostat = open_status)
            if (open_status /= 0) call ExceptionHandler(64)

            write(udf,'(a)') 'Binned_average_cospectra_sorted_by_time_of_day.'
            write(udf,'(a)') ',00:00-02:59,,,,,03:00-5:59,,,,,06:00-08:59&
                &,,,,,09:00-11:59,,,,,12:00-14:59,,,,,15:00-17:59,,,,,&
                &18:00-20:59,,,,,21:00-23:59'

            !> Add number of cospectra per class
            dataline = ''
            call AddDatum(dataline, '', separator)
            do cls = 1, 8
                do gas = w_ts, w_gas4
                    call WriteDatumInt(MeanBinCosp(1, cls)%cnt(gas), &
                        datum, EddyProProj%err_label)
                    call AddDatum(dataline, 'n_=_' // trim(adjustl(datum)), &
                        separator)
                end do
            end do
            write(udf,'(a)') dataline(1:len_trim(dataline) - 1)
            !> Add header piece
            write(udf,'(a)') &
                'nat_freq,avrg_cosp(w/T),avrg_cosp(w/co2),avrg_cosp(w/h2o),&
                &avrg_cosp(w/ch4),avrg_cosp(w/' // g4lab(1:g4l) //')&
                &,avrg_cosp(w/T),avrg_cosp(w/co2),avrg_cosp(w/h2o),&
                &avrg_cosp(w/ch4),avrg_cosp(w/' // g4lab(1:g4l) //')&
                &,avrg_cosp(w/T),avrg_cosp(w/co2),avrg_cosp(w/h2o),&
                &avrg_cosp(w/ch4),avrg_cosp(w/' // g4lab(1:g4l) //')&
                &,avrg_cosp(w/T),avrg_cosp(w/co2),avrg_cosp(w/h2o),&
                &avrg_cosp(w/ch4),avrg_cosp(w/' // g4lab(1:g4l) //')&
                &,avrg_cosp(w/T),avrg_cosp(w/co2),avrg_cosp(w/h2o),&
                &avrg_cosp(w/ch4),avrg_cosp(w/' // g4lab(1:g4l) //')&
                &,avrg_cosp(w/T),avrg_cosp(w/co2),avrg_cosp(w/h2o),&
                &avrg_cosp(w/ch4),avrg_cosp(w/' // g4lab(1:g4l) //')&
                &,avrg_cosp(w/T),avrg_cosp(w/co2),avrg_cosp(w/h2o),&
                &avrg_cosp(w/ch4),avrg_cosp(w/' // g4lab(1:g4l) //')&
                &,avrg_cosp(w/T),avrg_cosp(w/co2),avrg_cosp(w/h2o),&
                &avrg_cosp(w/ch4),avrg_cosp(w/' // g4lab(1:g4l) //')'

            do i = 1, nbins - 1
                call clearstr(dataline)
                if (MeanBinCosp(i, goodj)%fn(pick) /= error) then
                    call WriteDatumFloat(MeanBinCosp(i, goodj)%fn(pick), &
                        datum, EddyProProj%err_label)
                    call AddDatum(dataline, datum, separator)
                    do cls = 1, 8
                        do gas = w_ts, w_gas4
                            if (MeanBinCospAvailable(cls, gas))then
                                call WriteDatumFloat(MeanBinCosp(i, goodj)%fn(pick) &
                                    * MeanBinCosp(i, cls)%of(gas), datum, &
                                        EddyProProj%err_label)
                                call AddDatum(dataline, datum, separator)
                            else
                                call AddDatum(dataline, &
                                    trim(adjustl(EddyProProj%err_label)), &
                                    separator)
                            end if
                        end do
                    end do
                    write(udf, '(a)') dataline(1:len_trim(dataline) - 1)
                end if
            end do
            close(udf)
        end if

        !> =====================================================================
        !> Stability-sorted cospectra and models,
        !> only if at least one fit succeed
        proceed = .false.
        do gas = w_ts, w_co2
            if ((MassPar(gas, unstable)%a0   /= error .or. &
                MassPar(gas, unstable)%fpeak /= error .or. &
                MassPar(gas, unstable)%mu    /= error) .and. &
                (MassPar(gas, stable)%a0     /= error .or. &
                MassPar(gas, stable)%fpeak   /= error .or. &
                MassPar(gas, stable)%mu      /= error)) then
                proceed = .true.
                exit
            end if
        end do

        !> If not one fit went good, alert on output and exit routine
        if (.not. proceed) then
            call ExceptionHandler(45)
            return
        end if

        Filename = EddyProProj%id(1:len_trim(EddyProProj%id)) // Stability_FilePadding  &
            // Timestamp_FilePadding // CsvExt
        FilePath = SpecDir(1:len_trim(SpecDir)) // Filename(1:len_trim(Filename))

        open(udf, file = FilePath, iostat = open_status)
        if (open_status /= 0) call ExceptionHandler(64)

        write(udf,'(a)') 'Ensemble_cospectra,fitted_Massman_cospectra_and_Kaimal_cospectra.'
        write(udf,'(a)') 'Massman_model:'
        write(udf,'(a)') 'n*Co(ws)/cov(ws)=a0*(fn/fpeak)/(1+(fn/fpeak)^(2mu))^(1.167/mu)'
        write(udf,'(a)') '(note_that_slope_parameter_"m"_is_fixed_to_m=0.75)'
        write(udf,'(a)') 'a0=normalization_factor'
        write(udf,'(a)') 'fpeak=frequency_at_which_cospectrum_attains_highest_value'
        write(udf,'(a)') 'mu=broadness_factor'
        write(udf,'(a)') ''
        write(udf,'(a)') 'For_Kaimal_ideal_cospectra:see_e.g._Moncrieff_et_al._(1997_JoH)_Eqs.12-16.'
        write(udf,'(a)') '(for_stable_stratifications_two_extreme_ideal_cospectra_are_reported_&
            &corresponding_to_z/L=0.01_(slightly_stable)_and_to_z/L=10.0_(very_stable).'
        write(udf,'(a)') '-----------------------------------------------------------------------'
        write(udf,'(a)') 'Massman_model_fit_parameters_for_this_run:'
        write(udf,'(a)') 'w/T (unstable):'
        write(udf,'(a, f10.4)') 'a0,', MassPar(w_ts, unstable)%a0
        write(udf,'(a, f10.4)') 'fpeak,', MassPar(w_ts, unstable)%fpeak
        write(udf,'(a, f10.4)') 'mu,', MassPar(w_ts, unstable)%mu
        write(udf,'(a)') 'w/T (stable):'
        write(udf,'(a, f10.4)') 'a0,', MassPar(w_ts, stable)%a0
        write(udf,'(a, f10.4)') 'fpeak,', MassPar(w_ts, stable)%fpeak
        write(udf,'(a, f10.4)') 'mu,', MassPar(w_ts, stable)%mu
        write(udf,'(a)') 'w/CO2 (unstable):'
        write(udf,'(a, f10.4)') 'a0,', MassPar(w_co2, unstable)%a0
        write(udf,'(a, f10.4)') 'fpeak,', MassPar(w_co2, unstable)%fpeak
        write(udf,'(a, f10.4)') 'mu,', MassPar(w_co2, unstable)%mu
        write(udf,'(a)') 'w/CO2 (stable):'
        write(udf,'(a, f10.4)') 'a0,', MassPar(w_co2, stable)%a0
        write(udf,'(a, f10.4)') 'fpeak,', MassPar(w_co2, stable)%fpeak
        write(udf,'(a, f10.4)') 'mu,', MassPar(w_co2, stable)%mu
        write(udf,'(a)') 'w/H2O (unstable):'
        write(udf,'(a, f10.4)') 'a0,', MassPar(w_h2o, unstable)%a0
        write(udf,'(a, f10.4)') 'fpeak,', MassPar(w_h2o, unstable)%fpeak
        write(udf,'(a, f10.4)') 'mu,', MassPar(w_h2o, unstable)%mu
        write(udf,'(a)') 'w/H2O (stable):'
        write(udf,'(a, f10.4)') 'a0,', MassPar(w_h2o, stable)%a0
        write(udf,'(a, f10.4)') 'fpeak,', MassPar(w_h2o, stable)%fpeak
        write(udf,'(a, f10.4)') 'mu,', MassPar(w_h2o, stable)%mu
        write(udf,'(a)') 'w/CH4 (unstable):'
        write(udf,'(a, f10.4)') 'a0,', MassPar(w_ch4, unstable)%a0
        write(udf,'(a, f10.4)') 'fpeak,', MassPar(w_ch4, unstable)%fpeak
        write(udf,'(a, f10.4)') 'mu,', MassPar(w_ch4, unstable)%mu
        write(udf,'(a)') 'w/CH4 (stable):'
        write(udf,'(a, f10.4)') 'a0,', MassPar(w_ch4, stable)%a0
        write(udf,'(a, f10.4)') 'fpeak,', MassPar(w_ch4, stable)%fpeak
        write(udf,'(a, f10.4)') 'mu,', MassPar(w_ch4, stable)%mu
        write(udf,'(a)') 'w/'// g4lab(1:g4l) // ' (unstable):'
        write(udf,'(a, f10.4)') 'a0,', MassPar(w_gas4, unstable)%a0
        write(udf,'(a, f10.4)') 'fpeak,', MassPar(w_gas4, unstable)%fpeak
        write(udf,'(a, f10.4)') 'mu,', MassPar(w_gas4, unstable)%mu
        write(udf,'(a)') 'w/'// g4lab(1:g4l) // ' (stable):'
        write(udf,'(a, f10.4)') 'a0,', MassPar(w_gas4, stable)%a0
        write(udf,'(a, f10.4)') 'fpeak,', MassPar(w_gas4, stable)%fpeak
        write(udf,'(a, f10.4)') 'mu,', MassPar(w_gas4, stable)%mu
        write(udf,'(a, f10.4)') '-----------------------------------------------------------------------'


        write(udf,'(a)') 'unstable_(-650<L<0),,,,,,,,,,,,,,,,,,,,,,,,,stable(0<L<1000)'
        !> Add header piece
        write(udf,'(a)') 'fn,avrg_cosp(w/T),fit_cosp(w/T),kaimal_cosp,,&
                        &fn,avrg_cosp(w/co2),fit_cosp(w/co2),kaimal_cosp,,&
                        &fn,avrg_cosp(w/h2o),fit_cosp(w/h2o),kaimal_cosp,,&
                        &fn,avrg_cosp(w/ch4),fit_cosp(w/ch4),kaimal_cosp,,&
                        &fn,avrg_cosp(w/' // g4lab(1:g4l) //'),fit_cosp(w/' &
                        // g4lab(1:g4l) //'),kaimal_cosp,,&
                        &fn,avrg_cosp(w/T),fit_cosp(w/T),kaimal_cosp_zL_0.01,kaimal_cosp_zL_10.0,,&
                        &fn,avrg_cosp(w/co2),fit_cosp(w/co2),kaimal_cosp_zL_0.01,kaimal_cosp_zL_10.0,,&
                        &fn,avrg_cosp(w/h2o),fit_cosp(w/h2o),kaimal_cosp_zL_0.01,kaimal_cosp_zL_10.0,,&
                        &fn,avrg_cosp(w/ch4),fit_cosp(w/ch4),kaimal_cosp_zL_0.01,kaimal_cosp_zL_10.0,,&
                        &fn,avrg_cosp(w/' // g4lab(1:g4l) //'),fit_cosp(w/' &
                        // g4lab(1:g4l) //'),kaimal_cosp_zL_0.01,kaimal_cosp_zL_10.0'


        do i = 1, ndkf
            call clearstr(dataline)
            !> Unstable
            do gas = w_ts, w_gas4
                if (MeanStabCospAvailable(unstable, gas))then
                    if (MeanStabilityCosp(i, unstable)%fn(gas) /= error .and. MeanStabilityCosp(i, unstable)%fn(gas) /= 0d0) then
                        call WriteDatumFloat(MeanStabilityCosp(i, unstable)%fn(gas), datum, EddyProProj%err_label)
                        call AddDatum(dataline, datum, separator)
                        !> Ensemble cospectrum
                        if (MeanStabilityCosp(i, unstable)%cnt(gas) > 20) then
                            call WriteDatumFloat(MeanStabilityCosp(i, unstable)%of(gas), datum, EddyProProj%err_label)
                            call AddDatum(dataline, datum, separator)
                        else
                            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
                        end if
                        !> Fitted model cospectrum
                        call WriteDatumFloat(dexp(func(MeanStabilityCosp(i, unstable)%fn(gas), &
                            MassPar(gas, unstable)%a0, MassPar(gas, unstable)%fpeak, &
                            MassPar(gas, unstable)%mu)), datum, EddyProProj%err_label)
                        call AddDatum(dataline, datum, separator)
                        !> Ideal cospectrum
                        call WriteDatumFloat(kaimal(MeanStabilityCosp(i, unstable)%fn(gas), -9999d0, 'unstable'), &
                            datum, EddyProProj%err_label)
                        call AddDatum(dataline, datum, separator)
                        call AddDatum(dataline, '', separator)
                    else
                        call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
                        call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
                        call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
                        call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
                        call AddDatum(dataline, '', separator)
                    end if
                else
                    call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
                    call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
                    call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
                    call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
                    call AddDatum(dataline, '', separator)
                end if
            end do

            !> Stable
            do gas = w_ts, w_gas4
                if (MeanStabCospAvailable(stable, gas))then
                    if (MeanStabilityCosp(i, stable)%fn(gas) /= error .and. MeanStabilityCosp(i, stable)%fn(gas) /= 0d0) then
                            call WriteDatumFloat(MeanStabilityCosp(i, stable)%fn(gas), datum, EddyProProj%err_label)
                            call AddDatum(dataline, datum, separator)
                        !> Ensemble cospectrum
                        if (MeanStabilityCosp(i, stable)%cnt(gas) > 10) then
                            call WriteDatumFloat(MeanStabilityCosp(i, stable)%of(gas), datum, EddyProProj%err_label)
                            call AddDatum(dataline, datum, separator)
                        else
                            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
                        end if
                        !> Fitted model cospectrum
                        call WriteDatumFloat(dexp(func(MeanStabilityCosp(i, stable)%fn(gas), &
                            MassPar(gas, stable)%a0, MassPar(gas, stable)%fpeak, &
                            MassPar(gas, stable)%mu)), datum, EddyProProj%err_label)
                        call AddDatum(dataline, datum, separator)
                        !> Ideal cospectrum
                        call WriteDatumFloat(kaimal(MeanStabilityCosp(i, stable)%fn(gas), 1d-2, 'stable'), &
                            datum, EddyProProj%err_label)
                        call AddDatum(dataline, datum, separator)
                        call WriteDatumFloat(kaimal(MeanStabilityCosp(i, stable)%fn(gas), 1d1 , 'stable'), &
                            datum, EddyProProj%err_label)
                        call AddDatum(dataline, datum, separator)
                        call AddDatum(dataline, '', separator)
                    else
                        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
                        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
                        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
                        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
                        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
                        call AddDatum(dataline, '', separator)
                    end if
                else
                    call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
                    call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
                    call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
                    call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
                    call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
                    call AddDatum(dataline, '', separator)
                end if
            end do
            write(udf, '(a)') dataline(1:len_trim(dataline) - 1)
        end do
        close(udf)
    end if
    write(*,'(a)') ' Done.'

end subroutine OutputSpectralAssessmentResults

!***************************************************************************
!
! \brief
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine GetFnIndex(LocSpec, nrow, ncol, goodj, pick)
    use m_fx_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    integer, intent(in) :: ncol
    type (MeanSpectraType), intent(in) :: LocSpec(nrow, ncol)
    integer, intent(out) :: goodj
    integer, intent(out) :: pick
    !> local variables
    integer :: i
    integer :: cls


    goodj = ierror
    pick = ierror

    ol: do cls = 1, ncol
        il: do i = 1, nrow - 1
                if (LocSpec(i, cls)%fn(ts) > 0d0) then
                    goodj = cls
                    pick = ts
                    exit ol
                end if
        end do il
    end do ol
    if (goodj == ierror) then
        ol1: do cls = 1, ncol
            il1: do i = 1, nrow - 1
                    if (LocSpec(i, cls)%fn(co2) > 0d0) then
                        goodj = cls
                        pick = co2
                        exit ol1
                    end if
            end do il1
        end do ol1
    end if
    if (goodj == ierror) then
        ol2: do cls = 1, ncol
            il2: do i = 1, nrow - 1
                    if (LocSpec(i, cls)%fn(ch4) > 0d0) then
                        goodj = cls
                        pick = ch4
                        exit ol2
                    end if
            end do il2
        end do ol2
    end if
    if (goodj == ierror) then
        ol3: do cls = 1, ncol
            il3: do i = 1, nrow - 1
                    if (LocSpec(i, cls)%fn(gas4) > 0d0) then
                        goodj = cls
                        pick = gas4
                        exit ol3
                    end if
            end do il3
        end do ol3
    end if
end subroutine GetFnIndex
