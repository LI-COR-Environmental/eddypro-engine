!***************************************************************************
! read_metadata_file.f90
! ----------------------
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
! \brief       Reads metadata file for site and setup information \n
!              expected in EddyPro .metadata format
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ReadMetadataFile(LocCol, MetaFile, IniFileNotFound, printout)
    use m_common_global_var
    implicit none
    !> in/out variables
    character(*), intent(in) :: MetaFile
    type(ColType), intent(inout) :: LocCol(MaxNumCol)
    logical, intent(in) :: printout
    logical, intent(out) :: IniFileNotFound


    !> parse ini file and store all numeric and character tags
    call ParseIniFile(MetaFile, '', ANTags, ACTags, &
        size(ANTags), size(ACTags), ANTagFound, ACTagFound, IniFileNotFound)

    !> selects only tags needed in this software,
    !> and store them in relevant variables
    call WriteEddyProMetadataVariables(LocCol, printout)
end subroutine ReadMetadataFile

!***************************************************************************
!
! \brief       Retrieves relevant variables from the tags found in the \n
!              metadata file expected in EddyPro .metadata format
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WriteEddyProMetadataVariables(LocCol, printout)
    use m_common_global_var
    implicit none
    !> in/out variables
    logical, intent(in) :: printout
    type(ColType), intent(inout) :: LocCol(MaxNumCol)
    !> local variables
    integer :: init_ac_instr
    integer :: init_an_instr
    integer :: leap_ac_instr
    integer :: leap_an_instr
    integer :: init_ac_col
    integer :: init_an_col
    integer :: leap_ac_col
    integer :: leap_an_col
    integer :: NumSkipCol
    integer :: i = 0
    integer :: j = 0
    include 'interfaces.inc'


    if (len_trim(ACTags(1)%value) > 0) then
        Metadata%logger_swver = SwVerFromString(trim(ACTags(1)%value))
    else
        Metadata%logger_swver = errSwVer
    end if

    Metadata%sitename = trim(adjustl(ACTags(9)%value))
    Metadata%canopy_height = 0d0
    Metadata%d = 0d0
    Metadata%z0 = 0d0

    !> Site characteristics
    !> Altitude [m]
    Metadata%alt = dble(ANTags(1)%value)

    !> Altitude cannot be lower than Dead Sea or higher than top of
    !> Mount Everest (which includes reasonable flying altitudes)
    if (Metadata%alt < -428d0 .or. Metadata%alt > 8850d0) then
        if (printout) then
            write(*, '(a)')
            call ExceptionHandler(80)
        end if
        Metadata%alt = 0d0
    end if

    !> Barometric pressure [Pa] (e.g. Campbell & Normann, 1998 -
    !> An introduction to Environmental Biophysics)
    Metadata%bar_press = 1d3 * 101.3d0 *dexp(-Metadata%alt / 8200d0)

    !> Latitude [deg]
    Metadata%lat = dble(ANTags(2)%value)

    !> Latitude cannot be less than -90 or more than 90
    if (Metadata%lat < -90d0 .or. Metadata%lat > 90d0) then
        if (printout) then
            write(*, '(a)')
            call ExceptionHandler(81)
        end if
        Metadata%lat = 0.001d0
    end if

    !> Longitude[deg]
    Metadata%lon = dble(ANTags(3)%value)

    !> Longitude cannot be less than -180 or more than 180
    if (Metadata%lon < -180d0 .or. Metadata%lon > 180d0) then
        if (printout) then
            write(*, '(a)')
            call ExceptionHandler(82)
        end if
        Metadata%lon = 0.001d0
    end if

    !> Canopy height [m]
    Metadata%canopy_height = dble(ANTags(4)%value)

    !> Canopy height cannot be less than 0
    if (Metadata%canopy_height < 0d0) then
        if (printout) then
            write(*, '(a)')
            call ExceptionHandler(83)
        end if
        Metadata%canopy_height = 0d0
    end if

    !> Displacement height [m]
    Metadata%d = dble(ANTags(5)%value)

    !> Displacement height cannot be <= 0 or larger than canopy height
    if (Metadata%d <= 0d0 .or. Metadata%d > Metadata%canopy_height) then
        if (printout) then
            write(*, '(a)')
            call ExceptionHandler(84)
        end if
        Metadata%d = Metadata%canopy_height * 0.67d0
    end if

    !> Roughness length [m]
    Metadata%z0 = dble(ANTags(6)%value)

    !> Roughness length cannot be <= 0 or larger than canopy height
    if (Metadata%z0 <= 0d0 .or. Metadata%z0 > Metadata%canopy_height) then
        if (printout) then
            write(*, '(a)')
            call ExceptionHandler(85)
        end if
        Metadata%z0 = max(Metadata%canopy_height * 0.15d0, 0.001d0)
    end if

    !> Further meta info
    EddyProLog%save_native = ACTags(16)%value(1:1) == '1'
    EddyProLog%timestamp = ACTags(17)%value(1:1) == '1'
    EddyProLog%enable_proc = ACTags(18)%value(1:1) == '1'
    EddyProLog%tstamp_end  = ACTags(20)%value(1:1) == '1'
    select case (ACTags(21)%value(1:1))
        case ('0')
        EddyProLog%native_format = 'tmp_ascii'
        case default
        EddyProLog%native_format = 'none'
    end select

    !> Instruments initializations
    Instr = NullInstrument

    !> Instruments description
    leap_ac_instr = 8
    leap_an_instr = 15
    init_ac_instr = 25 - leap_ac_instr
    init_an_instr = 10 - leap_an_instr
    NumInstruments = 0
    Instr%firm = 'none'
    Instr%model = 'none'
    do i = 1, MaxNumInstruments
        if(ACTagFound(init_ac_instr + i*leap_ac_instr)) then
            NumInstruments = NumInstruments + 1
            Instr(i)%firm = ACTags(init_ac_instr + i*leap_ac_instr)%value &
                (1:len_trim(ACTags(init_ac_instr + i*leap_ac_instr)%value))

            Instr(i)%model = ACTags(init_ac_instr + i*leap_ac_instr + 2)%value &
                (1:len_trim(ACTags(init_ac_instr + i*leap_ac_instr + 2)%value))
            Instr(i)%height = dble(ANTags(init_an_instr + i*leap_an_instr)%value)

            !> For "generic" instruments, retrieve path lengths, time response,
            !> and extinction factors as applicable
            select case (Instr(i)%model(1:len_trim(Instr(i)%model) - 2))
                case ('generic_sonic')
                    Instr(i)%hpath_length = dble(ANTags(init_an_instr + i*leap_an_instr + 8)%value) * 1d-2 !< cm to m
                    Instr(i)%vpath_length = dble(ANTags(init_an_instr + i*leap_an_instr + 9)%value) * 1d-2 !< cm to m
                    Instr(i)%tau = dble(ANTags(init_an_instr + i*leap_an_instr + 10)%value)
                case ('generic_open_path', 'generic_closed_path')
                    Instr(i)%hpath_length = dble(ANTags(init_an_instr + i*leap_an_instr + 8)%value) * 1d-2 !< cm to m
                    Instr(i)%vpath_length = dble(ANTags(init_an_instr + i*leap_an_instr + 9)%value) * 1d-2 !< cm to m
                    Instr(i)%tau = dble(ANTags(init_an_instr + i*leap_an_instr + 10)%value)
                case ('open_path_krypton', 'open_path_lyman', 'closed_path_krypton', 'closed_path_lyman')
                    Instr(i)%hpath_length = dble(ANTags(init_an_instr + i*leap_an_instr + 8)%value) * 1d-2 !< cm to m
                    Instr(i)%vpath_length = dble(ANTags(init_an_instr + i*leap_an_instr + 9)%value) * 1d-2 !< cm to m
                    Instr(i)%tau = dble(ANTags(init_an_instr + i*leap_an_instr + 10)%value)
                    Instr(i)%kw = dble(ANTags(init_an_instr + i*leap_an_instr + 11)%value)
                    Instr(i)%ko = dble(ANTags(init_an_instr + i*leap_an_instr + 12)%value)
            end select

            !> If instrument firm is empty, retrieve it from the instrument model
            if (len(Instr(i)%firm) == 0 .or. Instr(i)%firm == 'none') then
                select case (Instr(i)%model(1:len_trim(Instr(i)%model) - 2))
                    case('li7500', 'li7500a', 'li7500rs', 'li7500ds', 'li7200', &
                        'li7200rs', 'li7700', 'li6262', 'li7000')
                        Instr(i)%firm = 'licor'
                    case('generic_open_path', 'generic_closed_path', 'open_path_krypton', &
                        'open_path_lyman', 'closed_path_krypton', 'closed_path_lyman')
                        Instr(i)%firm = 'other_irga'
                    case('hs_50', 'hs_100', 'r2', 'r3_50', 'r3_100', 'r3a_100', 'wm', 'wmpro')
                        Instr(i)%firm = 'gill'
                    case('usa1_standard', 'usa1_fast', 'usoni3_classa_mp', 'usoni3_cage_mp')
                        Instr(i)%firm = 'metek'
                    case('csat3', 'csat3b')
                        Instr(i)%firm = 'csi'
                    case('81000', '81000v', '81000re', '81000vre')
                        Instr(i)%firm = 'young'
                    case('generic_sonic')
                        Instr(i)%firm = 'other_sonic'
                end select
            end if

            !> If anemoeter model is "generic", set firm to "other" in any case
            if (index(Instr(i)%model, 'generic_sonic') /= 0) Instr(i)%firm = 'other_sonic'

            !> Select whether the instrument is a sonic or a gas analyser
            select case(Instr(i)%firm)
                case('licor', 'other_irga')
                    Instr(i)%category = 'irga'
                case('gill', 'metek', 'young', 'csi', 'other_sonic')
                    Instr(i)%category = 'sonic'
            end select
            !> Instrument software version
            if (len_trim(ACTags(init_ac_instr + i*leap_ac_instr + 1)%value) > 0) &
                Instr(i)%sw_ver_string = &
                    trim(ACTags(init_ac_instr + i*leap_ac_instr + 1)%value)

            !>Replace dash or underscore with dot if the case
            Instr(i)%sw_ver_string = replace(Instr(i)%sw_ver_string, '-', '.', &
                len(Instr(i)%sw_ver_string))
            Instr(i)%sw_ver_string = replace(Instr(i)%sw_ver_string, '_', '.', &
                len(Instr(i)%sw_ver_string))

            !> Create SwVer object from sw version string
            Instr(i)%sw_ver = SwVerFromString(Instr(i)%sw_ver_string)

            !> Retrieve absolute sensor separations (from reference sonic)
            Instr(i)%nsep = dble(ANTags(init_an_instr &
                + i*leap_an_instr + 2)%value) * 1d-2 !< in meters
            Instr(i)%esep = dble(ANTags(init_an_instr &
                + i*leap_an_instr + 3)%value) * 1d-2 !< in meters
            Instr(i)%vsep = dble(ANTags(init_an_instr &
                + i*leap_an_instr + 4)%value) * 1d-2 !< in meters

            !> Retrieve sonic parameters
            if (Instr(i)%category == 'sonic') then
                Instr(i)%wformat = ACTags(init_ac_instr + i*leap_ac_instr + 5)%value &
                    (1:len_trim(ACTags(init_ac_instr + i*leap_ac_instr + 5)%value))
                Instr(i)%wref    = ACTags(init_ac_instr + i*leap_ac_instr + 6)%value &
                    (1:len_trim(ACTags(init_ac_instr + i*leap_ac_instr + 6)%value))
                Instr(i)%head_corr = .false.
                if (ACTags(init_ac_instr + i*leap_ac_instr + 7)%value(1:1) == '1') &
                    Instr(i)%head_corr = .true.
                Instr(i)%north_offset = dble(ANTags(init_an_instr + i*leap_an_instr + 1)%value)
            end if

            !> If format of wind components is not set, it is likely to be u,v,w
            select case (Instr(i)%wformat(1:len_trim(Instr(i)%wformat)))
                case ('uvw', 'polar_w', 'axis')
                    continue
                case default
                    Instr(i)%wformat = 'uvw'
            end select

            !> If sonic is a Gill and wref is unknown, assumes SPAR configuration (the most logic, with
            !> U aligned to North. In case of wrong guess, there is only 30 degree offset in wind speed)
            if (Instr(i)%firm == 'gill') then
                select case (Instr(i)%wref)
                    case ('axis','spar')
                        continue
                    case default
                        Instr(i)%wref = 'spar'
                end select
            end if

            !> Retrieve gas analyser parameters
            if (Instr(i)%category == 'irga') then
                Instr(i)%tube_d = error
                Instr(i)%tube_l = error
                Instr(i)%tube_f = error

                select case (Instr(i)%model(1:len_trim(Instr(i)%model) - 2))
                    case ('li7700', 'li7500', 'li7500a', 'li7500rs', 'li7500ds', &
                        'generic_open_path', 'open_path_krypton', &
                        'open_path_lyman')
                        Instr(i)%path_type = 'open'
                    case default
                        Instr(i)%path_type = 'closed'
                        Instr(i)%tube_d = &
                            dble(ANTags(init_an_instr + i*leap_an_instr + 5)%value) * 1d-3 !< in meters
                        Instr(i)%tube_l = &
                            dble(ANTags(init_an_instr + i*leap_an_instr + 6)%value) * 1d-2  !< in meters
                        Instr(i)%tube_f = &
                            dble(ANTags(init_an_instr + i*leap_an_instr + 7)%value) / 6d4  !< in m+3s-1
                end select
            end if
        end if
    end do

    !> Raw file info
    Metadata%ac_freq = dble(ANTags(7)%value)
    Metadata%file_length = dble(ANTags(8)%value)
    FileInterpreter%header_rows = nint(ANTags(9)%value)
    FileInterpreter%data_label = ACTags(65)%value(1:len_trim(ACTags(65)%value))

    select case (ACTags(23)%value(1:len_trim(ACTags(23)%value)))
        case('comma')
        FileInterpreter%separator = ','
        case('semicolon')
        FileInterpreter%separator = ';'
        case('space')
        FileInterpreter%separator = ' '
        case('tab')
        FileInterpreter%separator = char(9)
        case default
        FileInterpreter%separator = ACTags(23)%value(1:1)
    end select
    FileInterpreter%discard_if_above = .false.
    if (ACTags(24)%value(1:1) == '1') FileInterpreter%discard_if_above = .true.

    FileInterpreter%file_with_text = .false.
    leap_ac_col = 7
    leap_an_col = 8
    init_ac_col = 66 - leap_ac_col
    init_an_col = 85 - leap_an_col
    LocCol%var = 'none'
    LocCol%label = 'none'
    LocCol%measure_type = 'none'
    LocCol%present = .false.
    LocCol%def_tl = 0d0
    LocCol%min_tl = 0d0
    LocCol%max_tl = 0d0
    NumAllVar = 0
    NumCol = 0
    NumSkipCol = 0
    FileWithFlags = .false.
    do i = 1, MaxNumCol
        if(ACTagFound(init_ac_col + i*leap_ac_col)) then
            !> Retrieve variable name
            NumCol = NumCol + 1
            LocCol(i)%var = ACTags(init_ac_col + i*leap_ac_col)%value &
                (1:len_trim(ACTags(init_ac_col + i*leap_ac_col)%value))
            LocCol(i)%label = LocCol(i)%var

            !> Check variable name. If empty, change it to "ignore".
            if (len(LocCol(i)%var) == 0) LocCol%var = 'ignore'

            !> Back compatibility: if a column is declared not_numeric, set flag "file with text" to .true.
            if (index(LocCol(i)%var, 'not_numeric') /= 0) FileInterpreter%file_with_text = .true.

            !> If variable is to be ignored, does so
            if (index(LocCol(i)%var, 'ignore') /= 0 .or. index(LocCol(i)%var, 'not_numeric') /= 0) then
                NumSkipCol = NumSkipCol + 1
                cycle
            end if

            !> Count total number of variables
            NumAllVar = NumAllVar + 1

            !> Read instrument
            LocCol(i)%instr_name = ACTags(init_ac_col + i*leap_ac_col + 3)%value &
                (1:len_trim(ACTags(init_ac_col + i*leap_ac_col + 3)%value))

            !> Associate instrument information to data colums
            do j = 1, NumInstruments
                if (index(LocCol(i)%instr_name, Instr(j)%model(1:len_trim(Instr(j)%model))) /= 0) then
                    LocCol(i)%instr = Instr(j)
                    exit
                end if
            end do

            !> Retrieve remaining column information
            LocCol(i)%measure_type = ACTags(init_ac_col + i*leap_ac_col + 2)%value &
                (1:len_trim(ACTags(init_ac_col + i*leap_ac_col + 2)%value))
            LocCol(i)%unit_in = ACTags(init_ac_col + i*leap_ac_col + 4)%value &
                (1:len_trim(ACTags(init_ac_col + i*leap_ac_col + 4)%value))
            LocCol(i)%conversion_type = ACTags(init_ac_col + i*leap_ac_col + 5)%value &
                (1:len_trim(ACTags(init_ac_col + i*leap_ac_col + 5)%value))
            LocCol(i)%unit_out = ACTags(init_ac_col + i*leap_ac_col + 6)%value &
                (1:len_trim(ACTags(init_ac_col + i*leap_ac_col + 6)%value))
            LocCol(i)%min    = dble(ANTags(init_an_col + i*leap_an_col)%value)
            LocCol(i)%max    = dble(ANTags(init_an_col + i*leap_an_col + 1)%value)
            LocCol(i)%a      = dble(ANTags(init_an_col + i*leap_an_col + 2)%value)
            LocCol(i)%b      = dble(ANTags(init_an_col + i*leap_an_col + 3)%value)
            LocCol(i)%def_tl = dble(ANTags(init_an_col + i*leap_an_col + 4)%value)
            LocCol(i)%min_tl = dble(ANTags(init_an_col + i*leap_an_col + 5)%value)
            LocCol(i)%max_tl = dble(ANTags(init_an_col + i*leap_an_col + 6)%value)

            !> Adjustment for robustness
            if (LocCol(i)%conversion_type /= 'gain_offset' &
                .and. LocCol(i)%conversion_type /= 'zero_fullscale') &
                LocCol(i)%conversion_type = 'none'
        end if
    end do
end subroutine WriteEddyProMetadataVariables
