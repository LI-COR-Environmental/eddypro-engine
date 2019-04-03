!***************************************************************************
! change_24_into_00.f90
! ---------------------
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
! \brief       Change 24:00 into 00:00 of day after
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine Change24Into00(Datestring, doy_format)
    use m_common_global_var
    implicit none
    !> in/out variables
    logical, intent(in) :: doy_format
    character(*), intent(inout) :: Datestring
    !> local variables
    integer :: int_year
    integer :: int_doy
    integer :: max_days
    integer :: i
    logical :: isleap
    character(5) :: tmp_date


    !> Determine whether current year is leap
    read(Datestring(1:4), '(i4)') int_year
    isleap = leapyear(int_year)
    if (isleap) then
        max_days = 366
    else
        max_days = 365
    end if

    if (doy_format) then
        if (isleap) then
            if (Datestring(10:11) == '24') then
                Datestring(10:11) = '00'
                do i = 1, lndays
                    tmp_date = Datestring(5:6) // '-' // Datestring(7:8)
                    if (DayOfLeapYear(i) == tmp_date) then
                        if (i /= lndays) then
                            Datestring(5:6) = DayOfLeapYear(i + 1)(1:2)
                            Datestring(7:8) = DayOfLeapYear(i + 1)(4:5)
                        else
                            Datestring(5:6) = DayOfLeapYear(1)(1:2)
                            Datestring(7:8) = DayOfLeapYear(1)(4:5)
                            int_year = int_year + 1
                            call int2char(int_year, Datestring(1:4), 0)
                        end if
                        exit
                    end if
                end do
            end if
        else
            if (Datestring(10:11) == '24') then
                Datestring(10:11) = '00'
                do i = 1, ndays
                    tmp_date = Datestring(5:6) // '-' // Datestring(7:8)
                    if (DayOfYear(i) == tmp_date) then
                        if (i /= ndays) then
                            Datestring(5:6) = DayOfYear(i + 1)(1:2)
                            Datestring(7:8) = DayOfYear(i + 1)(4:5)
                        else
                            Datestring(5:6) = DayOfYear(1)(1:2)
                            Datestring(7:8) = DayOfYear(1)(4:5)
                            int_year = int_year + 1
                            call int2char(int_year, Datestring(1:4), 0)
                        end if
                        exit
                    end if
                end do
            end if
        end if
    else
        if (Datestring(10:11) == '24') then
            Datestring(10:11) = '00'
            read(Datestring(5:7), '(i3)') int_doy
            int_doy = int_doy + 1
            if (int_doy == max_days + 1) then
                int_doy = 1
                int_year = int_year + 1
                call int2char(int_year, Datestring(1:4), 0)
            end if
            call int2char(int_doy, Datestring(5:7), 0)
        end if
    end if
end subroutine Change24Into00
