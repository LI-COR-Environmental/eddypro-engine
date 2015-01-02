!***************************************************************************
! read_ini_rp.f90
! ---------------
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
! \brief       Reads file "processing.eddypro", placed in the predifined program \n
!              folder "prog_folder\ini\" and stores relevant variables
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ReadIniRP(key)
    use m_rp_global_var
    implicit none
    !> in/out variables
    character(*), intent(in) :: key
    !> local variables
    logical :: IniFileNotFound


    write(*,'(a)') ' Reading EddyPro project file: ' &
                     // PrjPath(1:len_trim(PrjPath)) // '..'

    !> parse processing.eddypro file and store [Project] variables, common to all programs
    call ParseIniFile(PrjPath, 'Project', EPPrjNTags, EPPrjCTags,&
        size(EPPrjNTags), size(EPPrjCTags), SNTagFound, SCTagFound, IniFileNotFound)

    if (IniFileNotFound) call ExceptionHandler(21)
    call WriteProcessingProjectVariables()

    !> parse processing.eddypro file and store all numeric and character tags
    call ParseIniFile(PrjPath, key, SNTags, SCTags, size(SNTags), size(SCTags),&
        SNTagFound, SCTagFound, IniFileNotFound)

    if (IniFileNotFound) call ExceptionHandler(21)
    !> selects only tags needed in this software, and store them in relevant variables
    call WriteVariablesRP()

    write(*,'(a)')   ' done.'
end subroutine ReadIniRP

!***************************************************************************
!
! \brief       Looks in "SNTags" and "SCTags" and retrieve variables used for \n
!              express processing.
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WriteVariablesRP()
    use m_rp_global_var
    implicit none
    !> local variables
    integer :: i
    integer :: init_an_flags
    integer :: leap_an_flags
    integer :: init_an_wsect
    integer :: leap_an_wsect
    integer :: init_prof_z
    integer :: hlen
    logical :: proceed

    !> Initializations
    Dir%main_in      = 'none'
    AuxFile%pf       = 'none'

    !> Flags for elimination of individual data points
    leap_an_flags = 3
    init_an_flags = 73 - leap_an_flags
    NumRawFlags = 0
    RPsetup%filter_by_raw_flags = .false.
    do i = 1, MaxNumRawFlags
        if (SNTagFound(init_an_flags + i*leap_an_flags) .and. &
            nint(SNTags(init_an_flags + i*leap_an_flags)%value) > 0) then
            NumRawFlags = NumRawFlags + 1
            RawFlag(NumRawFlags)%col = nint(SNTags(init_an_flags + i*leap_an_flags)%value)
            RawFlag(NumRawFlags)%threshold = dble(SNTags(init_an_flags + i*leap_an_flags + 1)%value)
            RawFlag(NumRawFlags)%upper = .true.
            if (SNTags(init_an_flags + i*leap_an_flags + 2)%value > 0) &
                RawFlag(NumRawFlags)%upper = .false.
        end if
    end do
    if (NumRawFlags > 0) RPsetup%filter_by_raw_flags = .true.

    !> In/out directories
    Dir%main_in = SCTags(1)%value(1:len_trim(SCTags(1)%value))
    if (len_trim(Dir%main_in) == 0) Dir%main_in = 'none'

    !> Everything about raw statistical tests
    Test%sr = SCTags(3)%value(1:1) == '1'
    Test%ar = SCTags(4)%value(1:1) == '1'
    Test%do = SCTags(5)%value(1:1) == '1'
    Test%al = SCTags(6)%value(1:1) == '1'
    Test%sk = SCTags(7)%value(1:1) == '1'
    Test%ds = SCTags(8)%value(1:1) == '1'
    Test%tl = SCTags(9)%value(1:1) == '1'
    Test%aa = SCTags(10)%value(1:1) == '1'
    Test%ns = SCTags(11)%value(1:1) == '1'

    !> method of spike removal
    RPSetup%despike_vm = SCTags(90)%value(1:1) /= '1'

    !> Spike removal test
    sr%num_spk = idint(SNTags(1)%value)
    sr%lim_u = SNTags(2)%value
    sr%hf_lim = SNTags(3)%value
    ar%lim = idint(SNTags(4)%value)
    ar%bins = idint(SNTags(5)%value)
    ar%hf_lim = SNTags(6)%value
    sr%lim_w   = SNTags(54)%value
    sr%lim_co2 = SNTags(55)%value
    sr%lim_h2o = SNTags(56)%value
    sr%lim_ch4 = SNTags(57)%value
    sr%lim_gas4 = SNTags(58)%value

    !> Dropout test
    do%extlim_dw = SNTags(7)%value
    do%hf1_lim = SNTags(8)%value
    do%hf2_lim = SNTags(9)%value

    !> Absolute limits test
    al%u_max = SNTags(10)%value
    al%w_max = SNTags(11)%value
    al%t_min = SNTags(12)%value
    al%t_max = SNTags(13)%value
    al%co2_min = SNTags(14)%value
    al%co2_max = SNTags(15)%value
    al%h2o_min = SNTags(16)%value
    al%h2o_max = SNTags(17)%value
    al%ch4_min = SNTags(59)%value
    al%ch4_max = SNTags(60)%value
    al%gas4_min = SNTags(61)%value
    al%gas4_max = SNTags(62)%value

    !> Skewness and Kurtosis
    sk%hf_skmin = SNTags(18)%value
    sk%hf_skmax = SNTags(19)%value
    sk%sf_skmin = SNTags(20)%value
    sk%sf_skmax = SNTags(21)%value
    sk%hf_kumin = SNTags(22)%value
    sk%hf_kumax = SNTags(23)%value
    sk%sf_kumin = SNTags(24)%value
    sk%sf_kumax = SNTags(25)%value

    !> Discontinuities
    ds%hf_uv = SNTags(26)%value
    ds%hf_w = SNTags(27)%value
    ds%hf_t = SNTags(28)%value
    ds%hf_co2 = SNTags(29)%value
    ds%hf_h2o = SNTags(30)%value
    ds%hf_ch4 = SNTags(63)%value
    ds%hf_gas4 = SNTags(64)%value
    ds%hf_var = SNTags(31)%value
    ds%sf_uv = SNTags(32)%value
    ds%sf_w = SNTags(33)%value
    ds%sf_t = SNTags(34)%value
    ds%sf_co2 = SNTags(35)%value
    ds%sf_h2o = SNTags(36)%value
    ds%sf_ch4 = SNTags(65)%value
    ds%sf_gas4 = SNTags(66)%value
    ds%sf_var = SNTags(37)%value

    !> Timelag
    tl%hf_lim = SNTags(38)%value
    tl%sf_lim = SNTags(39)%value
    tl%def_co2 = SNTags(40)%value
    tl%def_h2o = SNTags(41)%value
    tl%def_ch4 = SNTags(67)%value
    tl%def_n2o = SNTags(68)%value

    !> Angle of attack
    aa%min = SNTags(42)%value
    aa%max = SNTags(43)%value
    aa%lim = SNTags(44)%value

    !> Non-statiorarity of horizontal wind
    ns%hf_lim = SNTags(45)%value

    !> select number of files to process together
    RPsetup%nfiles = nint(SNTags(47)%value)

    !> select angle-of-attack calibration option
    select case (SCTags(12)%value(1:1))
        case ('0')
        RPsetup%calib_aoa = 'none'
        case ('1')
        RPsetup%calib_aoa = 'nakai_12'
        case ('2')
        RPsetup%calib_aoa = 'nakai_06'
    end select

    !> Cross-wind correction
    RPsetup%calib_cw = SCTags(13)%value(1:1) == '1'
    !> select whether to look for raw files in sub-folders (recursive file listing)
    RPsetup%recurse = SCTags(19)%value(1:1) == '1'
    !> select whether to output binned (co)spectra
    RPsetup%out_bin_sp = SCTags(26)%value(1:1) == '1'
    !> select whether to output binned ogives
    RPsetup%out_bin_og = SCTags(51)%value(1:1) == '1'
    !> select whether to convert to mixing ratio
    RPsetup%to_mixing_ratio = SCTags(55)%value(1:1) == '1'

    !> select output file
    RPsetup%out_full_sp(u)   = SCTags(27)%value(1:1) == '1'
    RPsetup%out_full_sp(v)   = SCTags(28)%value(1:1) == '1'
    RPsetup%out_full_sp(w)   = SCTags(29)%value(1:1) == '1'
    RPsetup%out_full_sp(ts)  = SCTags(30)%value(1:1) == '1'
    RPsetup%out_full_sp(co2) = SCTags(31)%value(1:1) == '1'
    RPsetup%out_full_sp(h2o) = SCTags(32)%value(1:1) == '1'
    RPsetup%out_full_sp(ch4) = SCTags(33)%value(1:1) == '1'
    RPsetup%out_full_sp(gas4) = SCTags(34)%value(1:1) == '1'

    RPsetup%out_full_cosp(w_u)   = SCTags(35)%value(1:1) == '1'
    RPsetup%out_full_cosp(w_v)   = SCTags(36)%value(1:1) == '1'
    RPsetup%out_full_cosp(w_ts)  = SCTags(37)%value(1:1) == '1'
    RPsetup%out_full_cosp(w_co2) = SCTags(38)%value(1:1) == '1'
    RPsetup%out_full_cosp(w_h2o) = SCTags(39)%value(1:1) == '1'
    RPsetup%out_full_cosp(w_ch4) = SCTags(40)%value(1:1) == '1'
    RPsetup%out_full_cosp(w_n2o) = SCTags(41)%value(1:1) == '1'

    RPsetup%out_st(1) = SCTags(42)%value(1:1) == '1'
    RPsetup%out_st(2) = SCTags(43)%value(1:1) == '1'
    RPsetup%out_st(3) = SCTags(44)%value(1:1) == '1'
    RPsetup%out_st(4) = SCTags(45)%value(1:1) == '1'
    RPsetup%out_st(5) = SCTags(46)%value(1:1) == '1'
    RPsetup%out_st(6) = SCTags(47)%value(1:1) == '1'
    RPsetup%out_st(7) = SCTags(48)%value(1:1) == '1'

    RPsetup%out_raw(1) = SCTags(68)%value(1:1) == '1'
    RPsetup%out_raw(2) = SCTags(69)%value(1:1) == '1'
    RPsetup%out_raw(3) = SCTags(70)%value(1:1) == '1'
    RPsetup%out_raw(4) = SCTags(71)%value(1:1) == '1'
    RPsetup%out_raw(5) = SCTags(72)%value(1:1) == '1'
    RPsetup%out_raw(6) = SCTags(73)%value(1:1) == '1'
    RPsetup%out_raw(7) = SCTags(74)%value(1:1) == '1'

    RPsetup%out_raw_var(u)   = SCTags(75)%value(1:1) == '1'
    RPsetup%out_raw_var(v)   = SCTags(76)%value(1:1) == '1'
    RPsetup%out_raw_var(w)   = SCTags(77)%value(1:1) == '1'
    RPsetup%out_raw_var(ts)  = SCTags(78)%value(1:1) == '1'
    RPsetup%out_raw_var(co2) = SCTags(79)%value(1:1) == '1'
    RPsetup%out_raw_var(h2o) = SCTags(80)%value(1:1) == '1'
    RPsetup%out_raw_var(ch4) = SCTags(81)%value(1:1) == '1'
    RPsetup%out_raw_var(gas4) = SCTags(82)%value(1:1) == '1'
    RPsetup%out_raw_var(te)  = SCTags(83)%value(1:1) == '1'
    RPsetup%out_raw_var(pe)  = SCTags(84)%value(1:1) == '1'

    !> If no spectral output is selected, identify this situation for skipping
    !> completely the spectral analysis
    RPsetup%do_spectral_analysis = .false.
    if (RPsetup%out_bin_sp .or. RPsetup%out_bin_og &
        .or. RPsetup%out_full_sp(u) .or. RPsetup%out_full_sp(v) .or. RPsetup%out_full_sp(w) .or. RPsetup%out_full_sp(ts) &
        .or. RPsetup%out_full_sp(co2) .or. RPsetup%out_full_sp(h2o) .or. RPsetup%out_full_sp(ch4) .or. RPsetup%out_full_sp(gas4) &
        .or. RPsetup%out_full_cosp(w_u) .or. RPsetup%out_full_cosp(w_v) .or. RPsetup%out_full_cosp(w_ts) &
        .or. RPsetup%out_full_cosp(w_co2) .or. RPsetup%out_full_cosp(w_h2o) .or. RPsetup%out_full_cosp(w_ch4) &
        .or. RPsetup%out_full_cosp(w_n2o)) RPsetup%do_spectral_analysis = .true.

    !> If no variable was selected for output, force out_raw to false regardless of
    !> user setting
    if (.not. (RPsetup%out_raw_var(u) .or. RPsetup%out_raw_var(v) .or. RPsetup%out_raw_var(w) &
    .or. RPsetup%out_raw_var(ts) .or. RPsetup%out_raw_var(co2) .or. RPsetup%out_raw_var(ch4) &
    .or. RPsetup%out_raw_var(gas4) .or. RPsetup%out_raw_var(te) .or. RPsetup%out_raw_var(pe))) then
        RPsetup%out_raw = .false.
    end if

    !> Raw dataset dir
    proceed = .false.
    do i = 1, 7
        if (RPsetup%out_raw(i)) then
            proceed = .true.
            exit
        end if
    end do
    !> Define header of raw dataset files
    raw_out_header = '   '
    hlen = 3
    if (RPsetup%out_raw_var(u))   then
        raw_out_header = raw_out_header(1:hlen) // 'u'
        hlen = hlen + 25
    end if
    if (RPsetup%out_raw_var(v))   then
        raw_out_header = raw_out_header(1:hlen) // 'v'
        hlen = hlen + 25
    end if
    if (RPsetup%out_raw_var(w))   then
        raw_out_header = raw_out_header(1:hlen) // 'w'
        hlen = hlen + 25
    end if
    if (RPsetup%out_raw_var(ts))  then
        raw_out_header = raw_out_header(1:hlen) // 'ts'
        hlen = hlen + 25
    end if
    if (RPsetup%out_raw_var(co2)) then
        raw_out_header = raw_out_header(1:hlen) // 'co2'
        hlen = hlen + 25
    end if
    if (RPsetup%out_raw_var(h2o)) then
        raw_out_header = raw_out_header(1:hlen) // 'h2o'
        hlen = hlen + 25
    end if
    if (RPsetup%out_raw_var(ch4)) then
        raw_out_header = raw_out_header(1:hlen) // 'ch4'
        hlen = hlen + 25
    end if
    if (RPsetup%out_raw_var(gas4)) then
        raw_out_header = raw_out_header(1:hlen) // '4th gas'
        hlen = hlen + 25
    end if
    if (RPsetup%out_raw_var(te))  then
        raw_out_header = raw_out_header(1:hlen) // 'air_t'
        hlen = hlen + 25
    end if
    if (RPsetup%out_raw_var(pe))  then
        raw_out_header = raw_out_header(1:hlen) // 'air_p'
        hlen = hlen + 25
    end if

    !> Output QC details
    RPsetup%out_qc_details = SCTags(85)%value(1:1) == '1'

    !> select detrending method
    select case (SCTags(14)%value(1:1))
        case ('0')
        Meth%det = 'ba'
        case ('1')
        Meth%det = 'ld'
        case ('2')
        Meth%det = 'rm'
        case ('3')
        Meth%det = 'ewa'
        case default
        Meth%det = 'ba'
    end select
    if (Meth%det(1:len_trim(Meth%det)) == 'ld' .or. &
        Meth%det(1:len_trim(Meth%det)) == 'rm' .or. &
        Meth%det(1:len_trim(Meth%det)) == 'ewa') RPsetup%Tconst = nint(SNTags(46)%value)

    !> select rotation method
    select case (SCTags(15)%value(1:1))
        case ('0')
        Meth%rot = 'none'
        case ('1')
        Meth%rot = 'double_rotation'
        case ('2')
        Meth%rot = 'triple_rotation'
        case ('3')
        Meth%rot = 'planar_fit'
        case ('4')
        Meth%rot = 'planar_fit_no_bias'
        case default
        Meth%rot = 'none'
    end select

    !> Planar fit extra settings
    if (index(Meth%rot, 'planar_fit') /= 0) then
        !> Whether to perfom planar fit on the fly or use previous results file
        RPsetup%pf_onthefly = .false.
        if (SCTags(56)%value(1:1) == '1') then
            RPsetup%pf_onthefly = .true.
        else
            AuxFile%pf = SCTags(57)%value(1:len_trim(SCTags(57)%value))
        end if
        !> Whether to subtract b0 from mean w
        RPsetup%pf_subtract_b0 = SCTags(96)%value(1:1) /= '1'
    end if

    !> select time lag handling method
    select case (SCTags(16)%value(1:1))
        case ('0')
        Meth%tlag = 'none'
        case ('1')
        Meth%tlag = 'constant'
        case ('2')
        Meth%tlag = 'maxcov&default'
        case ('3')
        Meth%tlag = 'maxcov'
        case ('4')
        Meth%tlag = 'tlag_opt'
        case default
        Meth%tlag = 'none'
    end select

    !> Time lag optimizer extra settings
    RPsetup%to_onthefly = .false.
    TimeLagOptSelected = .false.
    if (Meth%tlag == 'tlag_opt') then
        TimeLagOptSelected = .true.
        if (SCTags(91)%value(1:1) == '1') then
            RPsetup%to_onthefly = .true.
        else
            AuxFile%to = SCTags(92)%value(1:len_trim(SCTags(92)%value))
        end if
    end if

    !>  tapering window
    select case (SCTags(17)%value(1:1))
        case ('0')
        RPsetup%tap_win = 'squared'
        case ('1')
        RPsetup%tap_win = 'bartlett'
        case ('2')
        RPsetup%tap_win = 'welch'
        case ('3')
        RPsetup%tap_win = 'hamming'
        case ('4')
        RPsetup%tap_win = 'hann'
    end select

    !> number of frequency bins
    Meth%spec%nbins = nint(SNTags(48)%value)

    !> max acceptable lack of data lines in a raw file
    RPsetup%max_lack = SNTags(49)%value

    !> select the averaging length. If zero, files are processed as they are
    RPsetup%avrg_len = nint(SNTags(50)%value)

    !> read wind speed offsets
    RPsetup%offset(u) = 0d0
    RPsetup%offset(v) = 0d0
    RPsetup%offset(w) = 0d0
    RPsetup%offset(u) = dble(SNTags(51)%value)
    RPsetup%offset(v) = dble(SNTags(52)%value)
    RPsetup%offset(w) = dble(SNTags(53)%value)

    !> Planar fit settings
    PFSetup%start_date    = SCTags(49)%value(1:len_trim(SCTags(49)%value))
    PFSetup%end_date      = SCTags(50)%value(1:len_trim(SCTags(50)%value))
    PFSetup%min_per_sec   = nint(SNTags(70)%value)
    PFSetup%w_max         = SNTags(71)%value
    PFSetup%u_min         = SNTags(72)%value
    !> If w_max is found to be < 0.099, it means it has not been set, so it
    !> is forced to the max value, which implies no filtering for w_max.
    if(PFSetup%w_max  <= 0.099d0) PFSetup%w_max = 10d0
    PFSetup%fix = 'clockwise'
    if(SCTags(88)%value(1:1) == '1') PFSetup%fix = 'counterclockwise'
    if(SCTags(88)%value(1:1) == '2') PFSetup%fix = 'double_rotation'

    !> Customization of wind sectors
    if (index(Meth%rot, 'planar_fit') /= 0) then
        leap_an_wsect = 2
        init_an_wsect = 209 - leap_an_wsect
        PFSetup%num_sec = 0
        PFSetup%north_offset = SNTags(208)%value
        do i = 1, MaxNumWSect
            if (SNTagFound(init_an_wsect + i*leap_an_wsect) .and. &
                SNTags(init_an_wsect + i*leap_an_wsect)%value > 0) then
                PFSetup%num_sec = PFSetup%num_sec + 1
                PFSetup%width(PFSetup%num_sec) = SNTags(init_an_wsect + i*leap_an_wsect)%value
                PFSetup%wsect_exclude(PFSetup%num_sec) = nint(SNTags(init_an_wsect + i*leap_an_wsect + 1)%value) == 1
            end if
        end do
        if (PFSetup%num_sec == 0) then
            call ExceptionHandler(40)
            PFSetup%num_sec = 1
        elseif (PFSetup%num_sec == 1) then
            PFSetup%wsect_end(PFSetup%num_sec) = 360
        elseif (PFSetup%num_sec > 1) then
            !> Calculate ending angle of each sector
            do i = 1, PFSetup%num_sec
                PFSetup%wsect_end(i) = nint(sum(PFSetup%width(1:i)))
                if (PFSetup%wsect_end(i) < 0) PFSetup%wsect_end(i) = 360 + PFSetup%wsect_end(i)
            end do
            PFSetup%wsect_end(PFSetup%num_sec) = 360
        end if
    end if

    !> Time lag optimizer settings
    if (RPsetup%to_onthefly) then
        TOSetup%h2o_nclass    = 0
        TOSetup%start_date    = SCTags(93)%value(1:len_trim(SCTags(93)%value))
        TOSetup%end_date      = SCTags(94)%value(1:len_trim(SCTags(94)%value))
        TOSetup%co2_min_flux  = SNTags(194)%value
        TOSetup%ch4_min_flux  = SNTags(195)%value
        TOSetup%gas4_min_flux = SNTags(196)%value
        TOSetup%le_min_flux   = SNTags(197)%value
        TOSetup%pg_range      = SNTags(198)%value
        TOSetup%min_lag(co2)   = SNTags(199)%value
        TOSetup%max_lag(co2)   = SNTags(200)%value
        TOSetup%min_lag(h2o)   = SNTags(201)%value
        TOSetup%max_lag(h2o)   = SNTags(202)%value
        TOSetup%min_lag(ch4)   = SNTags(203)%value
        TOSetup%max_lag(ch4)   = SNTags(204)%value
        TOSetup%min_lag(gas4)  = SNTags(205)%value
        TOSetup%max_lag(gas4)  = SNTags(206)%value
        TOSetup%h2o_nclass    = nint(SNTags(207)%value)
        if (TOSetup%h2o_nclass > 1) then
            TOSetup%h2o_class_size = floor(100d0 / TOSetup%h2o_nclass)
        end if
    end if

    !> Random error estimation settings
    select case (nint(SNTags(281)%value))
        case(1)
            RUsetup%meth = 'finkelstein_sims_01'
        case(2)
            RUsetup%meth = 'mann_lenschow_94'
        case(3)
            RUsetup%meth = 'tbd'
        case default
            RUsetup%meth = 'none'
    end select
    if (RUsetup%meth /= 'none') then
        select case (nint(SNTags(282)%value))
            case(1)
                RUsetup%its_meth = 'cross_0'
            case(2)
                RUsetup%its_meth = 'full_integral'
            case default
                RUsetup%its_meth = 'cross_e'
        end select
        RUsetup%its_sec_fact = nint(SNTags(283)%value)
        RUsetup%tlag_max = nint(SNTags(284)%value)
    end if

    !> Biomet measurements
    select case (SCTags(61)%value(1:len_trim(SCTags(61)%value)))
        case('comma')
            bFileMetadata%separator = ','
        case('semicolon')
            bFileMetadata%separator = ';'
        case('space')
            bFileMetadata%separator = ' '
        case('tab')
            bFileMetadata%separator = char(9)
        case default
            bFileMetadata%separator = SCTags(61)%value(1:1)
    end select
    bFileMetadata%tstamp_ref = SCTags(62)%value(1: len_trim(SCTags(62)%value))
    bFileMetadata%nhead = nint(SNTags(192)%value)

    !> Wheter to filter for spikes and abolute limits
    RPsetup%filter_sr = SCTags(63)%value(1:1) == '1'
    RPsetup%filter_al = SCTags(64)%value(1:1) == '1'

    !> Burba correction params
    RPsetup%bu_corr = 'none'
    if(SCTags(65)%value == '1')  RPsetup%bu_corr = 'yes'
    if(SCTags(65)%value == '-1') RPsetup%bu_corr = 'smart'
    RPsetup%bu_multi = SCTags(66)%value(1:1) == '1'

    !> Whether to use power-of-two samples for FFT
    RPsetup%power_of_two = SCTags(87)%value(1:1) == '1'

    !> Whether to refer wind direction to geographic north
    RPsetup%use_geo_north = SCTags(89)%value(1:1) == '1'

    magnetic_declination = 0d0
    if (RPsetup%use_geo_north) &
        magnetic_declination = nint(SNTags(193)%value)

    !> Biomet measurements numeric params
    bSetup%sel(bTa)   = nint(SNTags(111)%value)
    bSetup%sel(bPa)   = nint(SNTags(112)%value)
    bSetup%sel(bRH)   = nint(SNTags(113)%value)
    bSetup%sel(bPPFD) = nint(SNTags(114)%value)
    bSetup%sel(bLWin) = nint(SNTags(115)%value)
    bSetup%sel(bRg)   = nint(SNTags(116)%value)

    init_prof_z = 120
    bSetup%zT    = error
    bSetup%zCO2  = error
    bSetup%zH2O  = error
    bSetup%zCH4  = error
    bSetup%zGAS4 = error
    do i = 1, MaxProfNodes
        if (nint(SNTags(init_prof_z + i)%value) >= 0) &
            bSetup%zT(i) = nint(SNTags(init_prof_z + i)%value)
        if (nint(SNTags(init_prof_z + MaxProfNodes + i)%value) >= 0) &
            bSetup%zCO2(i) = nint(SNTags(init_prof_z + MaxProfNodes + i)%value)
        if (nint(SNTags(init_prof_z + 2 * MaxProfNodes + i)%value) >= 0) &
            bSetup%zH2O(i)  = nint(SNTags(init_prof_z + 2 * MaxProfNodes + i)%value)
        if (nint(SNTags(init_prof_z + 3 * MaxProfNodes + i)%value) >= 0) &
            bSetup%zCH4(i)  = nint(SNTags(init_prof_z + 3 * MaxProfNodes + i)%value)
        if (nint(SNTags(init_prof_z + 4 * MaxProfNodes + i)%value) >= 0) &
            bSetup%zGAS4(i) = nint(SNTags(init_prof_z + 4 * MaxProfNodes + i)%value)
    end do
    do i = 1, MaxProfNodes - 1
        if (bSetup%zT(i + 1) /= error .and. bSetup%zT(i) /= error) then
            bSetup%dz(1, i) = bSetup%zT(i + 1) - bSetup%zT(i)
        else
            bSetup%dz(1, i) = error
        end if
        bSetup%dz(2, i) = bSetup%zCO2(i + 1)  - bSetup%zCO2(i)
        bSetup%dz(3, i) = bSetup%zH2O(i + 1)  - bSetup%zH2O(i)
        bSetup%dz(4, i) = bSetup%zCH4(i + 1)  - bSetup%zCH4(i)
        bSetup%dz(5, i) = bSetup%zGAS4(i + 1) - bSetup%zGAS4(i)
    end do

    !> Parameters for Burba correction
    !> Multiple linear regressions
    BurbaPar%m(daytime, bot, 1)     =  SNTags(156)%value
    BurbaPar%m(daytime, bot, 2)     =  SNTags(157)%value
    BurbaPar%m(daytime, bot, 3)     =  SNTags(158)%value
    BurbaPar%m(daytime, bot, 4)     =  SNTags(159)%value
    BurbaPar%m(daytime, top, 1)     =  SNTags(160)%value
    BurbaPar%m(daytime, top, 2)     =  SNTags(161)%value
    BurbaPar%m(daytime, top, 3)     =  SNTags(162)%value
    BurbaPar%m(daytime, top, 4)     =  SNTags(163)%value
    BurbaPar%m(daytime, spar, 1)    =  SNTags(164)%value
    BurbaPar%m(daytime, spar, 2)    =  SNTags(165)%value
    BurbaPar%m(daytime, spar, 3)    =  SNTags(166)%value
    BurbaPar%m(daytime, spar, 4)    =  SNTags(167)%value
    BurbaPar%m(nighttime, bot, 1)   =  SNTags(168)%value
    BurbaPar%m(nighttime, bot, 2)   =  SNTags(169)%value
    BurbaPar%m(nighttime, bot, 3)   =  SNTags(170)%value
    BurbaPar%m(nighttime, bot, 4)   =  SNTags(171)%value
    BurbaPar%m(nighttime, top, 1)   =  SNTags(172)%value
    BurbaPar%m(nighttime, top, 2)   =  SNTags(173)%value
    BurbaPar%m(nighttime, top, 3)   =  SNTags(174)%value
    BurbaPar%m(nighttime, top, 4)   =  SNTags(175)%value
    BurbaPar%m(nighttime, spar, 1)  =  SNTags(176)%value
    BurbaPar%m(nighttime, spar, 2)  =  SNTags(177)%value
    BurbaPar%m(nighttime, spar, 3)  =  SNTags(178)%value
    BurbaPar%m(nighttime, spar, 4)  =  SNTags(179)%value
    !> Simple linear regressions
    BurbaPar%l(daytime, bot, 1)     =  SNTags(180)%value
    BurbaPar%l(daytime, bot, 2)     =  SNTags(181)%value
    BurbaPar%l(daytime, top, 1)     =  SNTags(182)%value
    BurbaPar%l(daytime, top, 2)     =  SNTags(183)%value
    BurbaPar%l(daytime, spar, 1)    =  SNTags(184)%value
    BurbaPar%l(daytime, spar, 2)    =  SNTags(185)%value
    BurbaPar%l(nighttime, bot, 1)   =  SNTags(186)%value
    BurbaPar%l(nighttime, bot, 2)   =  SNTags(187)%value
    BurbaPar%l(nighttime, top, 1)   =  SNTags(188)%value
    BurbaPar%l(nighttime, top, 2)   =  SNTags(189)%value
    BurbaPar%l(nighttime, spar, 1)  =  SNTags(190)%value
    BurbaPar%l(nighttime, spar, 2)  =  SNTags(191)%value

    !> Settings related to drift correction
    !> initializations
    DriftCorr%method = 'none'
    select case (nint(SNTags(300)%value))
        case(1)
            DriftCorr%method = 'linear'
        case(2)
            DriftCorr%method = 'signal_strength'
        case default
            DriftCorr%method = 'none'
    end select
    DriftCorr%dir_cal = error
    DriftCorr%inv_cal = error
    DriftCorr%b = error
    DriftCorr%c = error
    !> read values
    DriftCorr%dir_cal(0:6, co2) = SNTags(301:307)%value
    DriftCorr%dir_cal(0:6, h2o) = SNTags(308:314)%value
    DriftCorr%inv_cal(0:6, co2) = SNTags(329:335)%value
    DriftCorr%inv_cal(0:6, h2o) = SNTags(336:342)%value
    DriftCorr%b = SNTags(370)%value
    DriftCorr%c = SNTags(371)%value

    !> If user doesn't want WPL correction, do not even convert to mixing ratio
    if(.not. EddyProProj%wpl) RPsetup%to_mixing_ratio = .false.

    !> adjust paths
    call AdjDir(Dir%main_in, slash)
    call AdjFilePath(AuxFile%pf, slash)
    call AdjFilePath(AuxFile%to, slash)
end subroutine WriteVariablesRP
