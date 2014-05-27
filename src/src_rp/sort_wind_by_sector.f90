!***************************************************************************
! sort_wind_by_sector.f90
! -----------------------
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
    real (kind = dbl), intent(in) :: Wind(nrow, 3)
    integer, intent(out) :: NumElem(PFSetup%num_sec)
    real(kind = dbl), intent(out) :: WindBySect(nrow, 3, PFSetup%num_sec)
    !> local variables
    integer :: i = 0
    integer :: sec = 0
    real(kind = dbl) :: WindDir


    write(*, '(a)') ' Sorting wind data by sector..'
    write(LogInteger, '(i2)') PFSetup%num_sec
    write(*, '(a, i1, a)') '  '// adjustl(trim(LogInteger)) &
        // ' wind sector(s) selected.'

    NumElem = 0
    do i = 1, nrow
        if (Wind(i, u) /= error .and. Wind(i, v) /= error &
            .and. Wind(i, w) /= error) then
            !> Calculate wind direction from Wind
            call WindDirection(Wind(i, u:w), &
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
