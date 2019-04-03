!***************************************************************************
! sort_wind_by_sector.f90
! -----------------------
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
! \brief       Sort wind data in specified wind sectors
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine SortWindBySector(Wind, nrow, NumElem, WindBySect)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    real(kind = dbl) , intent(in) :: Wind(nrow, 3)
    integer, intent(out) :: NumElem(PFSetup%num_sec)
    real(kind = dbl), intent(out) :: WindBySect(nrow, 3, PFSetup%num_sec)
    !> local variables
    integer :: i = 0
    integer :: sec = 0
    real(kind = dbl) :: WindDir


    write(*, '(a)') ' Sorting wind data by sector..'
    write(LogInteger, '(i2)') PFSetup%num_sec
    write(*, '(a, i1, a)') '  '// trim(adjustl(LogInteger)) &
        // ' wind sector(s) selected.'

    NumElem = 0
    do i = 1, nrow
        if (Wind(i, u) /= error .and. Wind(i, v) /= error &
            .and. Wind(i, w) /= error) then
            !> Calculate wind direction from Wind
            call SingleWindDirection(Wind(i, u:w), &
                E2Col(u)%instr%north_offset + magnetic_declination, WindDir)

            !> Retrieve wind sector
            call WindSector(WindDir, sec)

            !> Assign wind to relevant sector
            if (sec <= 0) cycle
            NumElem(sec) = NumElem(sec) + 1
            WindBySect(NumElem(sec), u:w, sec) = Wind(i, u:w)
        end if
    end do

    write(*, '(a)')   ' Done.'
end subroutine SortWindBySector
