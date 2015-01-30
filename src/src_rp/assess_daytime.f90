!***************************************************************************
! assess_daytime.f90
! ------------------
! Copyright (C) 2007-2011, Eco2s team, Gerardo Fratini
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
! \brief       Assess whether it is daytime or night-time, based on
!              global radiation, PPFD or potential radiation
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine AssessDaytime(date, time)
    use m_rp_global_var
    implicit none
    !> in/out variables
    character(10), intent(in) :: date
    character(5), intent(in) :: time
    logical, external :: IsDaytime

    Stats%daytime = .false.
    !> Based on Rg
    if (biomet%val(bRg) > 12d0) Stats%daytime = .true.
    !> Based on PPFD
    if (.not. Stats%daytime .and. biomet%val(bPPFD) > 100d0) &
        Stats%daytime = .true.
    !> based on period timestamp and potential radiation
    if (biomet%val(bRg) == error .and. biomet%val(bPPFD) == error) &
        Stats%daytime = IsDaytime(PotRad, date, time)
end subroutine AssessDaytime
