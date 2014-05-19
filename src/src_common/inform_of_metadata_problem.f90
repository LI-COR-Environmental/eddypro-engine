!***************************************************************************
! inform_of_metadata_problem.f90
! ------------------------------
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
! \brief       Output diagnostic information on why metadata is not valid
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine InformOfMetadataProblem(passed, faulty_col)
    use m_common_global_var
    implicit none
    !> in/out variables
    logical, intent(in) :: passed(32)
    integer, intent(in) :: faulty_col
    !> local variables


    if (.not. passed(2)) then
        write(*, '(a)')  ' Invalid metadata: field separator not supported.'
        call log_msg( ' err=invalid metadata: field separator not supported.')
    end if

    if (.not. passed(3)) then
        write(*, '(a)')  ' Invalid metadata: anemometer firm not recognized for at least one&
            & anemometric variable.'
        write(LogInteger, '(i3)') faulty_col
        call schrinkstring(LogInteger)
        write(*, '(a)')  ' Error detected for column n.: ' // LogInteger(1:len_trim(LogInteger))
        call log_msg( ' err=invalid metadata: anemometer firm not recognized for at least one&
            & anemometric variable..')
    end if

    if (.not. passed(4)) then
        write(*, '(a)')  ' Invalid metadata: anemometer model not recognized for at least one&
            & anemometric variable.'
        write(LogInteger, '(i3)') faulty_col
        call schrinkstring(LogInteger)
        write(*, '(a)')  ' Error detected for column n.: ' // LogInteger(1:len_trim(LogInteger))
        call log_msg( ' err=invalid metadata: anemometer model not recognized for at least one&
            & anemometric variable.')
    end if

    if (.not. passed(5)) then
        write(*, '(a)')  ' Invalid metadata: gas analyzer firm not recognized for at least one&
            & gas concentration or density measurement.'
        write(LogInteger, '(i3)') faulty_col
        call schrinkstring(LogInteger)
        write(*, '(a)')  ' Error detected for column n.: ' // LogInteger(1:len_trim(LogInteger))
        call log_msg( ' err=invalid metadata: gas analyzer firm not recognized for at least one&
            & gas concentration or density measurement.')
    end if

    if (.not. passed(6)) then
        write(*, '(a)')  ' Invalid metadata: gas analyzer model not recognized for at least one&
            & gas concentration or density measurement.'
        write(LogInteger, '(i3)') faulty_col
        call schrinkstring(LogInteger)
        write(*, '(a)')  ' Error detected for column n.: ' // LogInteger(1:len_trim(LogInteger))
        call log_msg( ' err=invalid metadata: gas analyzer model not recognized for at least one&
            & gas concentration or density measurement.')
    end if

    if (.not. passed(7)) then
        write(*, '(a)')  ' Invalid metadata: gas analyzer firm not recognized for at least one&
            & cell temperature or pressure measurement.'
        write(LogInteger, '(i3)') faulty_col
        call schrinkstring(LogInteger)
        write(*, '(a)')  ' Error detected for column n.: ' // LogInteger(1:len_trim(LogInteger))
        call log_msg( ' err=invalid metadata: gas analyzer firm not recognized for at least one&
            & cell temperature or pressure measurement.')
    end if

    if (.not. passed(8)) then
        write(*, '(a)')  ' Invalid metadata: gas analyzer for at least one cell temperature or&
            & pressure is not a closed/enclosed path analyzer.'
        write(LogInteger, '(i3)') faulty_col
        call schrinkstring(LogInteger)
        write(*, '(a)')  ' Error detected for column n.: ' // LogInteger(1:len_trim(LogInteger))
        call log_msg( ' err=invalid metadata: gas analyzer for at least one cell temperature or&
            & pressure is not a closed/enclosed path analyzer.')
    end if

    if (.not. passed(9)) then
        write(*, '(a)')  ' Invalid metadata: attributes of the sampling tube (length, diameter&
            & and flow rate) must be positive numbers.'
        write(LogInteger, '(i3)') faulty_col
        call schrinkstring(LogInteger)
        write(*, '(a)')  ' Error detected for column n.: ' // LogInteger(1:len_trim(LogInteger))
        call log_msg( ' err=invalid metadata: attributes of the sampling tube (length, diameter&
            & and flow rate) must be positive numbers.')
    end if

    if (.not. passed(10)) then
        write(*, '(a)') ' Invalid metadata: unrecognized type of gas measurement. Gas measurement&
            & must be either "molar/mass density", "mole fraction" or "mixing ratio".'
        write(LogInteger, '(i3)') faulty_col
        call schrinkstring(LogInteger)
        write(*, '(a)')  ' Error detected for column n.: ' // LogInteger(1:len_trim(LogInteger))
        call log_msg( ' err=invalid metadata: unrecognized type of gas measurement. Gas measurement&
            & must be either "molar/mass density", "mole fraction" or "mixing ratio".')
    end if

    if (.not. passed(11)) then
        write(*, '(a)') ' Invalid metadata: invalid units for at least one gas measurement.&
            & Valid units for gas measurements are "ppt", "ppm", "ppb", "mmol/m^3", "umol/m^3",&
            & "g/m3", "mg/m3", "' // char(181) // 'g/m^3".'
        write(LogInteger, '(i3)') faulty_col
        call schrinkstring(LogInteger)
        write(*, '(a)')  ' Error detected for column n.: ' // LogInteger(1:len_trim(LogInteger))
        call log_msg( ' err=invalid metadata: invalid units for at least one gas measurement.&
            & Valid units for gas measurements are "ppt", "ppm", "ppb", "mmol/m^3", "umol/m^3".&
            & "g/m3", "mg/m3", "' // char(181) // 'g/m^3".')
    end if

    if (.not. passed(12)) then
        write(*, '(a)') ' Invalid metadata: invalid units for at least one anemometric variable.&
            & Valid units for wind components and speed-of-sound are "m/sec", "mm/sec", "cm/sec".'
        write(LogInteger, '(i3)') faulty_col
        call schrinkstring(LogInteger)
        write(*, '(a)')  ' Error detected for column n.: ' // LogInteger(1:len_trim(LogInteger))
        call log_msg( ' err=invalid metadata: invalid units for at least one anemometric variable.&
            & Valid units for wind components and speed-of-sound are "m/sec", "mm/sec", "cm/sec".')
    end if

    if (.not. passed(13)) then
        write(*, '(a)') ' Invalid metadata: invalid units for at least one temperature reading.&
            & Valid units for temperatures are "K", "cK", "C", "cC".'
        write(LogInteger, '(i3)') faulty_col
        call schrinkstring(LogInteger)
        write(*, '(a)')  ' Error detected for column n.: ' // LogInteger(1:len_trim(LogInteger))
        call log_msg( ' err=invalid metadata: invalid units for at least one temperature reading.&
            & Valid units for temperatures are "K", "cK", "C", "cC".')
    end if

    if (.not. passed(14)) then
        write(*, '(a)') ' Invalid metadata: invalid units for at least one pressure reading.&
            & Valid units for pressure are "Pa", "hPa", "kPa".'
        write(LogInteger, '(i3)') faulty_col
        call schrinkstring(LogInteger)
        write(*, '(a)')  ' Error detected for column n.: ' // LogInteger(1:len_trim(LogInteger))
        call log_msg( ' err=invalid metadata: invalid units for at least one pressure reading.&
            & Valid units for pressure are "Pa", "hPa", "kPa".')
    end if

    if (.not. passed(15)) then
        write(*, '(a)') ' Invalid metadata: if the Zero/Full-Scale conversion is select,&
            & Mimimum and Maximum (input) values cannot be both zero.'
        write(LogInteger, '(i3)') faulty_col
        call schrinkstring(LogInteger)
        write(*, '(a)')  ' Error detected for column n.: ' // LogInteger(1:len_trim(LogInteger))
        call log_msg( ' err=invalid metadata: if the Zero/Full-Scale conversion is select,&
            & Mimimum and Maximum (input) values cannot be both zero.')
    end if

    if (.not. passed(16)) then
        write(*, '(a)') ' Invalid metadata: if the Zero/Full-Scale conversion is select,&
            & A and B (output) values cannot be both zero.'
        write(LogInteger, '(i3)') faulty_col
        call schrinkstring(LogInteger)
        write(*, '(a)')  ' Error detected for column n.: ' // LogInteger(1:len_trim(LogInteger))
        call log_msg( ' err=invalid metadata: if the zero/full-scale conversion is select,&
            & A and B (output) values cannot be both zero.')
    end if

    if (.not. passed(17)) then
        write(*, '(a)') ' Invalid metadata: if the Gain/Offset conversion is selected,&
            & A value (gain) cannot be zero.'
        write(LogInteger, '(i3)') faulty_col
        call schrinkstring(LogInteger)
        write(*, '(a)')  ' Error detected for column n.: ' // LogInteger(1:len_trim(LogInteger))
        call log_msg( ' err=invalid metadata: if the gain/offset conversion is select,&
            & A value (gain) cannot be zero.')
    end if

    if (.not. passed(18)) then
        write(*, '(a)') ' Invalid metadata: at least one among u, v, w and a fast temperature is missing.'
        call log_msg( ' err=invalid metadata: at least one among u, v, w and a fast temperature is missing.')
    end if

    if (.not. passed(19)) then
        write(*, '(a)')  ' Invalid metadata: attributes of a "generic anemometer" (path lengths&
            & and time constant) must be positive numbers.'
        write(LogInteger, '(i3)') faulty_col
        call schrinkstring(LogInteger)
        write(*, '(a)')  ' Error detected for column n.: ' // LogInteger(1:len_trim(LogInteger))
        call log_msg( ' err=invalid metadata: attributes of a "generic anemometer" (path lengths&
            & and time constant) must be positive numbers.')
    end if

    if (.not. passed(20)) then
        write(*, '(a)')  ' Invalid metadata: attributes of a "generic" analyser (path lengths&
            & and time constant) must be positive numbers.'
        write(LogInteger, '(i3)') faulty_col
        call schrinkstring(LogInteger)
        write(*, '(a)')  ' Error detected for column n.: ' // LogInteger(1:len_trim(LogInteger))
        call log_msg( ' err=invalid metadata: attributes of a "generic" analyser (path lengths&
            & and time constant) must be positive numbers.')
    end if

    if (.not. passed(21)) then
        write(*, '(a)')  ' Invalid metadata: krypton or lyman-alpha analysers only&
            & supported for H2O measurements.'
        write(LogInteger, '(i3)') faulty_col
        call schrinkstring(LogInteger)
        write(*, '(a)')  ' Error detected for column n.: ' // LogInteger(1:len_trim(LogInteger))
        call log_msg( ' err=invalid metadata: attributes of a krypton or a lyman-alpha analyser (path lengths&
            & and time constant) must be different from zero.')
    end if

    if (.not. passed(22)) then
        write(*, '(a)')  ' Invalid metadata: attributes of a krypton or a lyman-alpha analyser (path lengths,&
            & time constant and extinction coefficients) must be different from zero.'
        write(LogInteger, '(i3)') faulty_col
        call schrinkstring(LogInteger)
        write(*, '(a)')  ' Error detected for column n.: ' // LogInteger(1:len_trim(LogInteger))
        call log_msg( ' err=invalid metadata: attributes of a krypton or a lyman-alpha analyser (path lengths&
            & and time constant) must be different from zero.')
    end if

    if (.not. passed(23)) then
        write(*, '(a)')  ' Invalid metadata: at least one variable is associated with an inexistent instrument.'
        write(*, '(a)')  ' If you compiled the metadata file with a text editor, check spelling &
            &of instrument models for all variables.'
        write(LogInteger, '(i3)') faulty_col
        call schrinkstring(LogInteger)
        write(*, '(a)')  ' Error detected for column n.: ' // LogInteger(1:len_trim(LogInteger))
        call log_msg( ' err=invalid metadata: attributes of a krypton or a lyman-alpha analyser (path lengths&
            & and time constant) must be different from zero.')
    end if

    if (.not. passed(24)) then
        write(*, '(a)')  ' Invalid metadata: it seems that a data column was selected for flux computation, &
            &which was declared either as "ignore" or "not numeric". &
            &Please check your metadata file and the selected data columns.'
        write(LogInteger, '(i3)') faulty_col
        call schrinkstring(LogInteger)
        write(*, '(a)')  ' Error detected for column n.: ' // LogInteger(1:len_trim(LogInteger))
        call log_msg( ' err=invalid metadata: it seems that a data column was selected for flux computation, &
            &which was declared either as "ignore" or "not numeric". &
            &Please check your metadata file and the selected data columns.')
    end if

    if (.not. passed(25)) then
        write(*, '(a)')  ' Invalid metadata: it seems that a fast ambient temperature measurement was associated &
            &to an LI-COR analyser. No LI-COR analyser provides an ambient temperature measurement suitable for &
            &sensible heat flux computation. Please use sonic temperature instead (i.e. select "None" in the "Fast temperature&
            & reading" Item).'
        write(LogInteger, '(i3)') faulty_col
        call schrinkstring(LogInteger)
        write(*, '(a)')  ' Error detected for column n.: ' // LogInteger(1:len_trim(LogInteger))
        call log_msg( ' err=invalid metadata: it seems that a fast ambient temperature measurement was associated &
            &to an LI-COR analyser. no LI-COR analyser provides an ambient temperature measurement suitable for &
            &sensible heat flux computation. Please use sonic temperature instead (i.e. select "none" in the "fast temperature&
            & reading" Item).')
    end if

    if (.not. passed(26)) then
        write(*, '(a)')  ' Invalid metadata: sonic temperature cannot be selected as a slow ambient temperature&
            &measurement. If you want to use it, just select "none" in the ambient temperature, and EddyPro will automatically&
            &use the sonic temperature, corrected for humidity effects, as ambient temperature.'
        write(LogInteger, '(i3)') faulty_col
        call schrinkstring(LogInteger)
        write(*, '(a)')  ' Error detected for column n.: ' // LogInteger(1:len_trim(LogInteger))
        call log_msg( ' err=invalid metadata: sonic temperature cannot be selected as a slow ambient temperature&
            &measurement. If you want to use it, just select "none" in the ambient temperature, and EddyPro will automatically&
            &use the sonic temperature, corrected for humidity effects, as ambient temperature.')
    end if
end subroutine InformOfMetadataProblem


