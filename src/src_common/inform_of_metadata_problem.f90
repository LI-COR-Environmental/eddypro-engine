!***************************************************************************
! inform_of_metadata_problem.f90
! ------------------------------
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


    write(LogInteger, '(i3)') faulty_col

    if (.not. passed(2)) then
        write(*,*) '  Warning(1001)> Invalid metadata.'
        write(*,*) '  Warning(1001)> Field separator not supported.'
    end if

    if (.not. passed(3)) then
        write(*,*) '  Warning(1001)> Invalid metadata.'
        write(*,*) '  Warning(1001)> Anemometer firm not recognized for at least one anemometric variable.'
        write(*,*) '  Warning(1001)> Problem detected for column n. ' // trim(adjustl(LogInteger))
    end if

    if (.not. passed(4)) then
        write(*,*) '  Warning(1001)> Invalid metadata.'
        write(*,*) '  Warning(1001)> Anemometer model not recognized for at least one anemometric variable.'
        write(*,*) '  Warning(1001)> Problem detected for column n. ' // trim(adjustl(LogInteger))
    end if

    if (.not. passed(5)) then
        write(*,*) '  Warning(1001)> Invalid metadata.'
        write(*,*) '  Warning(1001)> Gas analyser firm not recognized for at least one &
                                     &gas concentration or density measurement.'
        write(*,*) '  Warning(1001)> Problem detected for column n. ' // trim(adjustl(LogInteger))
    end if

    if (.not. passed(6)) then
        write(*,*) '  Warning(1001)> Invalid metadata.'
        write(*,*) '  Warning(1001)> Gas analyser model not recognized for at least one &
                                     &gas concentration or density measurement.'
        write(*,*) '  Warning(1001)> Problem detected for column n. ' // trim(adjustl(LogInteger))
    end if

    if (.not. passed(7)) then
        write(*,*) '  Warning(1001)> Invalid metadata.'
        write(*,*) '  Warning(1001)> Gas analyser firm not recognized for at least one &
                                     &cell temperature or pressure measurement.'
        write(*,*) '  Warning(1001)> Problem detected for column n. ' // trim(adjustl(LogInteger))
    end if

    if (.not. passed(8)) then
        write(*,*) '  Warning(1001)> Invalid metadata.'
        write(*,*) '  Warning(1001)> Gas analyser for at least one cell temperature or &
                                     &pressure is not a closed/enclosed path analyser.'
        write(*,*) '  Warning(1001)> Problem detected for column n. ' // trim(adjustl(LogInteger))
    end if

    if (.not. passed(9)) then
        write(*,*) '  Warning(1001)> Invalid metadata.'
        write(*,*) '  Warning(1001)> Attributes of the sampling tube (length, diameter &
                                     &and flow rate) must be positive numbers.'
        write(*,*) '  Warning(1001)> Problem detected for column n. ' // trim(adjustl(LogInteger))
    end if

    if (.not. passed(10)) then
        write(*,*) '  Warning(1001)> Invalid metadata.'
        write(*,*) '  Warning(1001)> Unrecognised type of gas measurement. Gas measurement&
                                     & must be either "molar/mass density", "mole fraction" or "mixing ratio".'
        write(*,*) '  Warning(1001)> Problem detected for column n. ' // trim(adjustl(LogInteger))
    end if

    if (.not. passed(11)) then
        write(*,*) '  Warning(1001)> Invalid metadata.'
        write(*,*) '  Warning(1001)> Invalid units for at least one gas measurement.'
        write(*,*) '  Warning(1001)> Valid units for gas measurements are "ppt", "ppm", "ppb", "mmol/m^3", &
                                     &"umol/m^3", "g/m3", "mg/m3", "' // char(181) // 'g/m^3".'
        write(*,*) '  Warning(1001)> Problem detected for column n. ' // trim(adjustl(LogInteger))
    end if

    if (.not. passed(12)) then
        write(*,*) '  Warning(1001)> Invalid metadata.'
        write(*,*) '  Warning(1001)> Invalid units for at least one anemometric variable.'
        write(*,*) '  Warning(1001)> Valid units for wind components and speed-of-sound are "m/sec", "mm/sec", &
                                     &"cm/sec".'
        write(*,*) '  Warning(1001)> Problem detected for column n. ' // trim(adjustl(LogInteger))
    end if

    if (.not. passed(13)) then
        write(*,*) '  Warning(1001)> Invalid metadata.'
        write(*,*) '  Warning(1001)> Invalid units for at least one temperature reading.'
        write(*,*) '  Warning(1001)> Valid units for temperatures are "K", "cK", "C", "cC".'
        write(*,*) '  Warning(1001)> Problem detected for column n. ' // trim(adjustl(LogInteger))
    end if

    if (.not. passed(14)) then
        write(*,*) '  Warning(1001)> Invalid metadata.'
        write(*,*) '  Warning(1001)> Invalid units for at least one pressure reading.'
        write(*,*) '  Warning(1001)> Valid units for pressure are "Pa", "hPa", "kPa".'
        write(*,*) '  Warning(1001)> Problem detected for column n. ' // trim(adjustl(LogInteger))
    end if

    if (.not. passed(15)) then
        write(*,*) '  Warning(1001)> Invalid metadata.'
        write(*,*) '  Warning(1001)> If the Zero/Full-Scale conversion is select, &
                                     &Mimimum and Maximum (input) values cannot both be zero.'
        write(*,*) '  Warning(1001)> Problem detected for column n. ' // trim(adjustl(LogInteger))
    end if

    if (.not. passed(16)) then
        write(*,*) '  Warning(1001)> Invalid metadata.'
        write(*,*) '  Warning(1001)> If the Zero/Full-Scale conversion is select, &
                                     &A and B (output) values cannot both be zero.'
        write(*,*) '  Warning(1001)> Problem detected for column n. ' // trim(adjustl(LogInteger))
    end if

    if (.not. passed(17)) then
        write(*,*) '  Warning(1001)> Invalid metadata.'
        write(*,*) '  Warning(1001)> If the Gain/Offset conversion is selected, &
                                     &A value (gain) cannot be zero.'
        write(*,*) '  Warning(1001)> Problem detected for column n. ' // trim(adjustl(LogInteger))
    end if

    if (.not. passed(18)) then
        write(*,*) '  Warning(1001)> Invalid metadata.'
        write(*,*) '  Warning(1001)> At least one among u, v, w and a fast temperature is missing.'
    end if

    if (.not. passed(19)) then
        write(*,*) '  Warning(1001)> Invalid metadata.'
        write(*,*) '  Warning(1001)> Attributes of a "Generic anemometer" (path lengths &
                                     &and time constant) must be positive numbers.'
        write(*,*) '  Warning(1001)> Problem detected for column n. ' // trim(adjustl(LogInteger))
    end if

    if (.not. passed(20)) then
        write(*,*) '  Warning(1001)> Invalid metadata.'
        write(*,*) '  Warning(1001)> Attributes of a "generic analyser" (path lengths &
                                     &and time constant) must be positive numbers.'
        write(*,*) '  Warning(1001)> Problem detected for column n. ' // trim(adjustl(LogInteger))
    end if

    if (.not. passed(21)) then
        write(*,*) '  Warning(1001)> Invalid metadata.'
        write(*,*) '  Warning(1001)> Krypton or lyman-alpha analysers only supported for H2O measurements.'
        write(*,*) '  Warning(1001)> Problem detected for column n. ' // trim(adjustl(LogInteger))
    end if

    if (.not. passed(22)) then
        write(*,*) '  Warning(1001)> Invalid metadata.'
        write(*,*) '  Warning(1001)> Attributes of a krypton or a lyman-alpha analyser (path lengths, &
                                     &time constant and extinction coefficients) must be different from zero.'
        write(*,*) '  Warning(1001)> Problem detected for column n. ' // trim(adjustl(LogInteger))
    end if

    if (.not. passed(23)) then
        write(*,*) '  Warning(1001)> Invalid metadata.'
        write(*,*) '  Warning(1001)> At least one variable is associated with an inexistent instrument.'
        write(*,*) '  Warning(1001)> If you compiled the metadata file with a text editor, check spelling &
                                     &of instrument models for all variables.'
        write(*,*) '  Warning(1001)> Problem detected for column n. ' // trim(adjustl(LogInteger))
    end if

    if (.not. passed(24)) then
        write(*,*) '  Warning(1001)> Invalid metadata.'
        write(*,*) '  Warning(1001)> It seems that a data column was selected for flux computation, &
                                     &which was declared either as "ignore" or "not numeric".'
        write(*,*) '  Warning(1001)> Check the metadata file and the selected data columns.'
        write(*,*) '  Warning(1001)> Problem detected for column n. ' // trim(adjustl(LogInteger))
    end if

    if (.not. passed(25)) then
        write(*,*) '  Warning(1001)> Invalid metadata.'
        write(*,*) '  Warning(1001)> It seems that a fast ambient temperature measurement was associated &
                                     &to an LI-COR analyser.'
        write(*,*) '  Warning(1001)> No LI-COR analyser provides an ambient temperature measurement suitable for &
                                     &sensible heat flux computation.'
        write(*,*) '  Warning(1001)> Use sonic temperature instead (i.e. select "None" in the "Fast temperature &
                                     &reading" item).'
        write(*,*) '  Warning(1001)> Problem detected for column n. ' // trim(adjustl(LogInteger))
    end if

    if (.not. passed(26)) then
        write(*,*) '  Warning(1001)> Invalid metadata.'
        write(*,*) '  Warning(1001)> Sonic temperature cannot be selected as a slow ambient temperature measurement.'
        write(*,*) '  Warning(1001)> If you want to use it, just select "None" in the ambient temperature, and EddyPro &
                                     &will automatically use the sonic temperature, corrected for humidity effects, &
                                     &as ambient temperature.'
        write(*,*) '  Warning(1001)> Problem detected for column n. ' // trim(adjustl(LogInteger))
    end if
end subroutine InformOfMetadataProblem


