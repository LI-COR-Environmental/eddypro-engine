!***************************************************************************
! wind_sector.f90
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
! \brief       Returns wind sector of given wind direction
! \author      Gerardo Fratini
! \note
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WindSector(wdir, n)
    use m_rp_global_var
    implicit none
    !> in/out variables
    real(kind = dbl), intent(in) :: wdir
    integer, intent(out) :: n
    !> local variables
    real(kind = dbl) :: twdir
    integer :: sec


    !> tilt wind direction according to north offset of planar fit
    twdir = wdir - PFSetup%north_offset
    if (twdir >= 360d0) twdir = twdir - 360d0
    if (twdir < 0d0) twdir = 360d0 + twdir

    !> Now treats everything as if offset is zero
    if(twdir >= 0d0 .and. twdir < dfloat(PFSetup%wsect_end(1))) then
        n = 1
    else
        do sec = 2, PFSetup%num_sec
            if(twdir >= PFSetup%wsect_end(sec - 1) .and. twdir < PFSetup%wsect_end(sec)) then
                n = sec
                exit
            end if
        end do
    end if
end subroutine WindSector
