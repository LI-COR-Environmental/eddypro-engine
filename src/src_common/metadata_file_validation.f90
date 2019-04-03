!***************************************************************************
! metadata_file_validation.f90
! ----------------------------
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
! \brief       Checks information in the metadata file, for a formal and \n
!              partially substantial correctness. Returns a flag such as \n
!              validation "passed" or "not passed".
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine MetadataFileValidation(LocCol, passed, faulty_col)
    use m_common_global_var
    implicit none
    !> in/out variables
    type(ColType), intent(in) :: LocCol(MaxNumCol)
    logical, intent(out) :: passed(32)
    integer, intent(out) :: faulty_col
    !> local variables
    integer :: i
    integer :: j
    logical :: InstrChecked(NumInstruments)
    integer :: present(ts)


    passed = .true.
    InstrChecked = .false.

    !> Check separator for ASCII (virtually impossible to fail, GUI filters)
    if  (EddyProProj%ftype == 'generic_ascii') then
        select case (FileInterpreter%separator)
            case (',',';',' ',char(9))
                continue
            case default
                passed(1) = .false.
                passed(2) = .false.
                return
        end select
    end if

    !> Each Instrument must be validated only if at least one EddyPro
    !> variable is measured with it. Here the instrument check is performed
    !> for columns with property use_it, (which is defined only after reading
    !> the processing project file .eddypro)
    do i = 1, NumCol
        if (LocCol(i)%useit) then
            !> Fast temperature measurements are exception here, because they
            !> can be measured with an un-described instrument.
            if (LocCol(i)%var == 'ts' &
                .and. LocCol(i)%instr%category == 'fast_t_sensor') then
                select case &
                    (LocCol(i)%instr%model(1:len_trim(LocCol(i)%instr%model)-2))
                    case ('li7500', 'li7500a', 'li7500rs', 'li7500ds', 'li7200', &
                        'li7200rs', 'li7700', 'li6262', 'li7000')
                        passed(1) = .false.
                        passed(25) = .false.
                        faulty_col = i
                        return
                    case default
                        cycle
                end select
            end if

            !> Similarly, slow temperature and pressure measurements can
            !> come from any instrument, except from a known sonic.
            if (LocCol(i)%var == 'air_t' .or. LocCol(i)%var == 'air_p') then
                select case &
                    (LocCol(i)%instr%model(1:len_trim(LocCol(i)%instr%model)-2))
                    case ('hs_50', 'hs_100', 'r2', 'r3_50', 'r3_100', &
                        'r3a_100', 'wm', 'wmpro', 'usa1_standard', &
                        'usa1_fast', 'usoni3_classa_mp', 'usoni3_cage_mp', &
                        'csat3', 'csat3b', &
                        '81000', '81000v', '81000re', '81000vre')
                        passed(1) = .false.
                        passed(26) = .false.
                        faulty_col = i
                        return
                    case default
                        cycle
                end select
            end if

            if (LocCol(i)%var == 'ignore' &
                .or. LocCol(i)%var == 'not_numeric') then
                passed(1) = .false.
                passed(24) = .false.
                faulty_col = i
                return
            end if

            !> detects instrument in instrument list and calls actual check
            il: do j = 1, NumInstruments
                if (index(Instr(j)%model, &
                    trim(adjustl(LocCol(i)%Instr%model))) /= 0) then
                    if (.not. InstrChecked(j)) then
                        call InstrumentValidation(Instr(j), LocCol(i), passed)
                        if (.not. passed(1)) then
                            faulty_col = i
                            return
                        end if
                        InstrChecked(j) = .true.
                    end if
                    exit il
                endif
            end do il

            !> If instrument was not found in the list, exit
            if (j == NumInstruments + 1) then
                passed(1) = .false.
                passed(23) = .false.
                faulty_col = i
                return
            end if

            !> Check column
            call ColumnValidation(LocCol(i), passed)
            if (.not. passed(1)) then
                faulty_col = i
                return
            end if
        end if
    end do

    !> Beyond a check of metadata file, checks if the minimum number of
    !> variables needed for flux calculation are present. These are u, v, w and
    !> either sos or ts. Also checks that there are no more than 1
    !> instance of each.
    present = 0
    do i= 1, NumCol
        if (LocCol(i)%var == 'u' .and. LocCol(i)%useit) &
            present(u) = present(u) + 1
        if (LocCol(i)%var == 'v' .and. LocCol(i)%useit) &
            present(v) = present(v) + 1
        if (LocCol(i)%var == 'w' .and. LocCol(i)%useit) &
            present(w) = present(w) + 1
        if ((LocCol(i)%var == 'ts' .or. LocCol(i)%var == 'sos') &
            .and. LocCol(i)%useit) &
            present(ts) = present(ts) + 1
    end do
    do i = u, ts
        if (present(i) /= 1) then
            passed(1) = .false.
            passed(18) = .false.
            return
        end if
    end do
end subroutine MetadataFileValidation

!***************************************************************************
!
! \brief       Checks Instrument description for formal correctness \n
!              Returns the flag "passed"
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine InstrumentValidation(LocInstr, LocCol, passed)
    use m_common_global_var
    implicit none
    !> in/out variables
    type(InstrumentType), intent(in) :: LocInstr
    type(ColType), intent(in) :: LocCol
    logical, intent(out) :: passed(32)

    passed = .true.

    !> Check firm & model
    select case (LocCol%var)
        !> Anemometric variables must come from a sonic anemometer
        case ('u', 'v', 'w', 'ts', 'sos')
            !> check firm
            select case (LocInstr%firm(1:len_trim(LocInstr%firm)))
                case ('gill', 'metek', 'young', 'csi', 'other_sonic')
                    continue
                case default
                    passed(1) = .false.
                    passed(3) = .false.
                    return
            end select
            !> check model
            select case (LocInstr%model(1:len_trim(LocInstr%model)-2))
                case ('hs_50', 'hs_100', 'r2', 'r3_50', 'r3_100', 'r3a_100', 'wm', 'wmpro', &
                      'usa1_standard', 'usa1_fast', 'csat3', 'csat3b', &
                      'usoni3_classa_mp', 'usoni3_cage_mp', &
                      '81000', '81000v', '81000re', '81000vre')
                      continue
                case ('generic_sonic')
                    if (LocInstr%hpath_length * LocInstr%vpath_length * LocInstr%tau == 0) then
                        passed(1) = .false.
                        passed(19) = .false.
                        return
                    end if
                case default
                    passed(1) = .false.
                    passed(4) = .false.
                    return
            end select
        !> Gas concentrations must come from a gas analyzer
        case ('co2', 'h2o', 'ch4', 'n2o')
            !> check firm
            select case (LocInstr%firm(1:len_trim(LocInstr%firm)))
                case ('licor', 'other_irga')
                    continue
                case default
                    passed(1) = .false.
                    passed(5) = .false.
                    return
            end select
            !> check model
            select case (LocInstr%model(1:len_trim(LocInstr%model)-2))
                case ('li7500', 'li7500a', 'li7500rs', 'li7500ds', 'li7200', &
                    'li7200rs', 'li7700', 'li6262', 'li7000')
                    continue
                case ('generic_open_path', 'generic_closed_path')
                    if (LocInstr%hpath_length * LocInstr%vpath_length * LocInstr%tau == 0) then
                        passed(1) = .false.
                        passed(20) = .false.
                        return
                    end if
                case ('open_path_krypton', 'open_path_lyman', 'closed_path_krypton', 'closed_path_lyman')
                    if (LocCol%var /= 'h2o') then
                        passed(1) = .false.
                        passed(21) = .false.
                        return
                    end if
                if (LocInstr%hpath_length * LocInstr%vpath_length * LocInstr%tau * &
                        LocInstr%kw * LocInstr%ko == 0) then
                        passed(1) = .false.
                        passed(22) = .false.
                        return
                    end if
                case default
                    passed(1) = .false.
                    passed(6) = .false.
                    return
            end select
        !> Cell temperature/pressure must come from a closed path gas analyzer
        case ('int_t_1', 'int_t_2', 'cell_t', 'int_p')
            !> check firm
            select case (LocInstr%firm(1:len_trim(LocInstr%firm)))
                case ('licor', 'other_irga')
                    continue
                case default
                    passed(1) = .false.
                    passed(7) = .false.
                    return
            end select
            !> check model
            select case (LocInstr%model(1:len_trim(LocInstr%model)-2))
                case ('li7200', 'li7200rs', 'li6262', 'li7000', 'generic_closed_path')
                    continue
                case default
                    passed(1) = .false.
                    passed(8) = .false.
                    return
            end select
    end select

    !> Check tube properties for closed path IRGAs
    if (LocInstr%path_type == 'closed') then
        if (LocInstr%tube_d <= 0 .or. LocInstr%tube_l <= 0 .or. LocInstr%tube_f <= 0) then
            passed(1) = .false.
            passed(9) = .false.
            return
        end if
    end if
end subroutine InstrumentValidation

!***************************************************************************
!
! \brief       Checks file column description for formal correctness \n
!              Returns the flag "passed"
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ColumnValidation(LocCol, passed)
    use m_common_global_var
    implicit none
    !> in/out variables
    type(ColType), intent(in) :: LocCol
    logical, intent(out) :: passed(32)
    !> local variables
    character(32) :: units

    passed = .true.

    !> Check measure type for gas concentrations
    select case (LocCol%var)
        case ('co2', 'h2o', 'ch4', 'n2o')
            select case (LocCol%measure_type)
                case ('molar_density', 'mole_fraction', 'mixing_ratio')
                    continue
                case default
                    passed(1) = .false.
                    passed(10) = .false.
                    return
            end select
        case default
            continue
    end select

    !> Check output units compatibility, in case conversion is performed or
    !> input units compatibility, in case conversion is not performed
    !> First select which units to check
    if (LocCol%conversion_type /= 'none') then
        units = LocCol%unit_out
    else
        units = LocCol%unit_in
    end if
    select case (LocCol%var)
        case ('co2', 'h2o', 'ch4', 'n2o')
            select case (units)
                case ('ppt', 'ppm', 'ppb', 'mmol_m3', 'umol_m3', 'g_m3', 'mg_m3', 'ug_m3')
                    continue
                case default
                    passed(1) = .false.
                    passed(11) = .false.
                    return
            end select
        case ('u', 'v', 'w', 'sos')
            select case (units)
                case ('m_sec', 'mm_sec', 'cm_sec')
                    continue
                case default
                    passed(1) = .false.
                    passed(12) = .false.
                    return
            end select
        case ('ts', 'int_t_1', 'int_t_2', 'cell_t', 'ait_t')
            select case (units)
                case ('kelvin', 'ckelvin', 'celsius', 'ccelsius')
                    continue
                case default
                    passed(1) = .false.
                    passed(13) = .false.
                    return
            end select
        case ('int_p', 'air_p')
            select case (units)
                case ('pa', 'hpa', 'kpa')
                    continue
                case default
                    passed(1) = .false.
                    passed(14) = .false.
                    return
            end select
    end select

    !> Check in/out ranges (or gain/offset) if conversion type is /= none
    if (LocCol%conversion_type == 'zero_fullscale') then
        if (LocCol%min == 0 .and. LocCol%max == 0) then
            passed(1) = .false.
            passed(15) = .false.
            return
        end if
        if (LocCol%a   == 0 .and. LocCol%b   == 0) then
            passed(1) = .false.
            passed(16) = .false.
            return
        end if
    end if
    if (LocCol%conversion_type == 'gain_offset') then
        if (LocCol%a   == 0) then
            passed(1) = .false.
            passed(17) = .false.
            return
        end if
    end if
end subroutine ColumnValidation
