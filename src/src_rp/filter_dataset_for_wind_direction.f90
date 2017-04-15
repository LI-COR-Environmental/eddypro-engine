!***************************************************************************
! filter_data_for_wind_direction.f90
! ----------------------------------
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
! \brief       Filters dataset when wind comes from user-selected directions
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo        
!***************************************************************************
subroutine FilterDatasetForWindDirection(Set, nrow, ncol)
    use m_rp_global_var
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(inout) :: Set(nrow, ncol)
    !> lcoal variables
    integer :: i
    integer :: sec
    real(kind = dbl) :: WD


    Essentials%m_wdf = Essentials%m_wdf + 1
    do i = 1, nrow
        !> Instantaneous wind direction
        call WindDirection(Set(i, u:w), &
            E2Col(u)%instr%north_offset + magnetic_declination, WD)

        !> Set record to error if wind is coming from an excluded sector
        do sec = 1, RPSetup%wdf_num_secs
            if (WD > RPSetup%wdf_start(sec) .and. WD < RPSetup%wdf_end(sec)) then
                Set (i, :) = error
                Essentials%m_wdf = Essentials%m_wdf + 1
            end if 
        end do
    end do
end subroutine FilterDatasetForWindDirection