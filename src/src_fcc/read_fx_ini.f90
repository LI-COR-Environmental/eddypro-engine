!***************************************************************************
! read_fx_ini.f90
! ---------------
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
! \brief       Reads .eddypro file and stores relevant variables
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ReadIniFX(key)
    use m_fx_global_var
    implicit none
    ! in/out variables
    character(*), intent(in) :: key
    !> local variables
    logical :: IniFileNotFound

    write(*,'(a)') ' Reading EddyPro project file: ' &
                     // PrjPath(1:len_trim(PrjPath)) // '..'

    !> parse processing.eddypro file and store [Project] variables,
    !> common to all programs
    call ParseIniFile(PrjPath, 'Project', EPPrjNTags, EPPrjCTags,&
        size(EPPrjNTags), size(EPPrjCTags), &
        SNTagFound, SCTagFound, IniFileNotFound)

    if (IniFileNotFound) call ExceptionHandler(21)
    call WriteProcessingProjectVariables()

    !> parse processing.eddypro file and store all numeric and character tags
    call ParseIniFile(PrjPath, key, SNTags, SCTags, size(SNTags), size(SCTags),&
        SNTagFound, SCTagFound, IniFileNotFound)

    if (IniFileNotFound) call ExceptionHandler(21)
    !> selects only tags needed in this software, and store
    !> them in relevant variables
    call WriteVariablesFX()

    write(*,'(a)')   ' done.'
end subroutine ReadIniFX

!***************************************************************************
!
! \brief       Looks in "SNTags" and "SCTags" and retrieve variables \n
!              used for express processing.
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WriteVariablesFX()
    use m_fx_global_var
    implicit none
    !> local variables
    integer :: i
    integer :: j
    integer :: gas
    integer :: skipped_classes
    integer :: start
    logical :: dirExists


    !> Spectra analysis time period
    if (SCTags(22)%value(1:1) == '1') then
        FCCsetup%SA%start_date = SCTags(1)%value(1:len_trim(SCTags(1)%value))
        FCCsetup%SA%end_date = SCTags(3)%value(1:len_trim(SCTags(3)%value))
    else
        FCCsetup%SA%start_date = '1900-01-01'
        FCCsetup%SA%end_date   = '2100-12-31'
    end if

    !> Essentials files
    AuxFile%ex = trim(adjustl(SCTags(7)%value))
    if (len_trim(AuxFile%ex) == 0) AuxFile%ex = 'none'

    !> Binned (co)spectra folder
    Dir%binned = trim(adjustl(SCTags(8)%value))
    if (len_trim(Dir%binned) == 0) &
        Dir%binned = trim(Dir%main_out) // trim(SubDirBinCospectra)

    !> Full (co)spectra folder
    Dir%full = trim(adjustl(SCTags(9)%value))
    if (len_trim(Dir%full) == 0) &
        Dir%full = trim(Dir%main_out) // trim(SubDirCospectra)

    !> Method of H correction
    select case (SCTags(13)%value(1:1))
        case ('0')
        FCCsetup%H_corr = 'vandijk04'
        case ('1')
        FCCsetup%H_corr = 'kaimal_gaynor91'
    end select

    !> Select low-pass spectra correction method
    FCCsetup%SA%in_situ = .false.
    FCCsetup%SA%ibrom_model = .false.
    FCCsetup%import_full_cospectra = .false.
    select case (EddyProProj%hf_meth)
        case ('horst_97')
            !> Correction after Horst (1997, BLM), in-situ/analytical
            FCCsetup%SA%in_situ = .true.
        case ('ibrom_07')
            !> Correction after Ibrom et al (2007, AFM) fully in-situ
            FCCsetup%SA%in_situ = .true.
            FCCsetup%SA%ibrom_model = .true.
        case ('fratini_12')
            !> Correction after Fratini et al. 2012, fully in-situ
            FCCsetup%SA%add_sonic_lptf = SCTags(17)%value(1:1) == '1'
            FCCsetup%SA%in_situ = .true.
            FCCsetup%import_full_cospectra = .true.
            FCCsetup%SA%ibrom_model = .true.
    end select

    !> Select whether and how to apply correction for sensors
    !> separation from Horst & Lenschow 2009
    FCCsetup%SA%horst_lens09  = 'none'
    select case (SCTags(21)%value(1:1))
        case ('1')
            FCCsetup%SA%horst_lens09  = 'full'
        case ('2')
            FCCsetup%SA%horst_lens09  = 'cross_and_vertical'
    end select

    !> Spectral assessment file path, if declared available
    AuxFile%sa   = 'none'
    if (SCTags(19)%value(1:1) == '0' .and. FCCsetup%SA%in_situ) then
        AuxFile%sa = SCTags(20)%value(1:len_trim(SCTags(20)%value))
        if (len_trim(AuxFile%sa) == 0) AuxFile%sa = 'none'
    end if

    !> Whether to perform spectral assessment on the fly
    !> Either because (1) in-situ spectral corrections were selected with
    !> on-the-fly spectral assessment or spectral assessment file was not found;
    !> (2) or because ensemble averaged spectra or (3) co-spectra have
    !> been requested
    FCCsetup%do_spectral_assessment = &
        FCCsetup%SA%in_situ .and. AuxFile%sa == 'none'

    FCCsetup%pass_thru_spectral_assessment = &
        FCCsetup%do_spectral_assessment &
        .or. EddyProProj%out_avrg_cosp &
        .or. EddyProProj%out_avrg_spec

    !> Check existence of binned cospectra directory if necessary
    if (FCCsetup%pass_thru_spectral_assessment) then
        inquire(file = Dir%binned, exist=dirExists)
        if (.not. dirExists) then
            EddyProProj%out_avrg_cosp = .false.
            EddyProProj%out_avrg_spec = .false.
            FCCsetup%do_spectral_assessment = .false.
            FCCsetup%pass_thru_spectral_assessment = .false.
            if (FCCsetup%SA%in_situ) then
                EddyProProj%hf_meth = 'moncrieff_97'
                FCCsetup%SA%in_situ = .false.
            end if
            call ExceptionHandler(87)
        end if
    end if

    !> Check existence of full cospectra directory if necessary
    if (EddyProProj%hf_meth == 'fratini_12') then
        inquire(file = Dir%full, exist=dirExists)
        if (.not. dirExists) then
            call ExceptionHandler(88)
            EddyProProj%hf_meth = 'moncrieff_97'
            FCCsetup%SA%in_situ = .false.
        end if
    end if

    FCCsetup%SA%lptf = 'none'
    if (EddyProProj%hf_meth == 'custom') then
        !> select low-pass transfer function (LPTF) definition method
        select case (SCTags(15)%value(1:1))
            case ('0')
                !> Transfer function calculated analytically
                !> (Moncrieff et al., 1997)
                FCCsetup%SA%lptf = 'analytic'
            case ('1')
                !> Transfer function calculated in-situ,
                !> after Aubinet et al. 2001 modified to account
                !> for RH-dependency
                FCCsetup%SA%lptf = 'sigma'
            case ('2')
                !> Ttransfer function calculated in-situ,
                !> after Ibrom et al. 2007
                FCCsetup%SA%lptf = 'iir'
            case default
                !> If not specified, set to none
                FCCsetup%SA%lptf = 'none'
        end select
    end if

    !> choose co-spectral models ('none' option not allowed here)
    select case (SCTags(16)%value(1:1))
        case ('0')
            FCCsetup%SA%cosp = 'from_data'
            FCCsetup%SA%cosp_type = 'in_situ'

        case ('1')
            FCCsetup%SA%cosp = 'fitted_from_data'
            FCCsetup%SA%cosp_type = 'in_situ'

        case ('2')
            FCCsetup%SA%cosp = 'moncrieff_97'
            FCCsetup%SA%cosp_type = 'analytic'

        case ('3')
            FCCsetup%SA%cosp = 'horst_97'
            FCCsetup%SA%cosp_type = 'analytic'

        case default
            FCCsetup%SA%cosp = 'none'
            FCCsetup%SA%cosp_type = 'none'
    end select

    if (FCCsetup%SA%lptf == 'sigma' .or. FCCsetup%SA%lptf == 'iir') then
        !> If a data-oriented approach is used for LPTF, an analytic part can be \n
        !> added, correcting for anemometer dynamical response and path-averaging
        FCCsetup%SA%add_sonic_lptf = SCTags(17)%value(1:1) == '1'
    end if

    !> Minimum number of samples per class for spectral averaging
    FCCsetup%SA%min_smpl = idint(dble(SNTags(2)%value))

    !> Minimum and maximum frequencies for transfer functions regression
    i = 3
    do gas = co2, gas4
        FCCsetup%SA%fmin(gas) = dble(SNTags(i)%value)
        FCCsetup%SA%fmax(gas) = dble(SNTags(i+1)%value)
        i = i + 2
    end do

    !> Flux thresholds, used to include/exclude corresponding (co)spectra in
    !> ensemble averages and model fits. Also used to discriminate between
    !> direct method and model in spectral correction after Fratini et al. (2012).
    FCCsetup%SA%min_un_ustar = dble(SNTags(92)%value)
    FCCsetup%SA%min_un_co2   = dble(SNTags(93)%value)
    FCCsetup%SA%min_un_ch4   = dble(SNTags(94)%value)
    FCCsetup%SA%min_un_gas4  = dble(SNTags(95)%value)
    FCCsetup%SA%min_un_LE    = dble(SNTags(96)%value)
    FCCsetup%SA%min_un_H     = dble(SNTags(97)%value)
    FCCsetup%SA%min_st_ustar = dble(SNTags(98)%value)
    FCCsetup%SA%min_st_co2   = dble(SNTags(99)%value)
    FCCsetup%SA%min_st_ch4   = dble(SNTags(100)%value)
    FCCsetup%SA%min_st_gas4  = dble(SNTags(101)%value)
    FCCsetup%SA%min_st_LE    = dble(SNTags(102)%value)
    FCCsetup%SA%min_st_H     = dble(SNTags(103)%value)
    FCCsetup%SA%max_ustar    = dble(SNTags(104)%value)
    FCCsetup%SA%max_co2      = dble(SNTags(105)%value)
    FCCsetup%SA%max_ch4      = dble(SNTags(106)%value)
    FCCsetup%SA%max_gas4     = dble(SNTags(107)%value)
    FCCsetup%SA%max_LE       = dble(SNTags(108)%value)
    FCCsetup%SA%max_H        = dble(SNTags(109)%value)

    !> Whether to use results of Vickers and Mahrt tests to eliminate (co)spectra
    FCCsetup%SA%filter_cosp_by_vm_flags = SCTags(23)%value(1:1) == '1'

    !> Minimum frequency for high-frequency noise detection and elimination
    FCCsetup%SA%hfn_fmin(co2)  = dble(SNTags(16)%value)
    FCCsetup%SA%hfn_fmin(h2o)  = dble(SNTags(17)%value)
    FCCsetup%SA%hfn_fmin(ch4)  = dble(SNTags(18)%value)
    FCCsetup%SA%hfn_fmin(gas4) = dble(SNTags(19)%value)

    !> Assign each month the relevant CO2 class. Max number of groups is 12.
    skipped_classes = 0
    start = 19
    do i = 1, MaxGasClasses
        if (SNTags(start + 2*i - 1)%value > 0d0) then
            do j = nint(SNTags(start + 2*i - 1)%value), nint(SNTags(start + 2*i)%value)
            FCCsetup%SA%class(co2, j) = i
            end do
        else
            skipped_classes = skipped_classes + 1
            cycle
        end if
    end do
    FCCsetup%SA%nclass(co2) = 12 - skipped_classes
    !> Assign each month the relevant CH4 class. Max number of groups is 12.
    skipped_classes = 0
    start = 19 + 24
    do i = 1, MaxGasClasses
        if (SNTags(start + 2*i - 1)%value > 0d0) then
            do j = nint(SNTags(start + 2*i - 1)%value), nint(SNTags(start + 2*i)%value)
            FCCsetup%SA%class(ch4, j) = i
            end do
        else
            skipped_classes = skipped_classes + 1
            cycle
        end if
    end do
    FCCsetup%SA%nclass(ch4) = 12 - skipped_classes
    !> Assign each month the relevant GAS4 class. Max number of groups is 12.
    skipped_classes = 0
    start = 19 + 48
    do i = 1, MaxGasClasses
        if (SNTags(start + 2*i - 1)%value > 0d0) then
            do j = nint(SNTags(start + 2*i - 1)%value), nint(SNTags(start + 2*i)%value)
            FCCsetup%SA%class(gas4, j) = i
            end do
        else
            skipped_classes = skipped_classes + 1
            cycle
        end if
    end do
    FCCsetup%SA%nclass(gas4) = 12 - skipped_classes

    !> adjust Dirs
    call AdjDir(Dir%binned, slash)
    call AdjDir(Dir%full, slash)
    call AdjFilePath(AuxFile%ex, slash)
    call AdjFilePath(AuxFile%sa, slash)
end subroutine WriteVariablesFX
