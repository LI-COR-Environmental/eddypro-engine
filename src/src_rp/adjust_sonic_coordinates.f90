!***************************************************************************
! adjust_sonic_coordinates.f90
! ----------------------------
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
! \file        src/adjust_sonic_coordinates.f90
! \brief       Adjust wind components and offsets in order to reference all
!              wind data to the same system, where:
!              u>0 is wind blowing from South to North
!              v>0 is wind blowing from East to West
!              w>0 is wind going upward
!              wind_dir = 0 is wind from north
!              wind_dir = 90 is wind from east
!              wind_dir = 180 is wind from south
!              wind_dir = 270 is wind from west
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine AdjustSonicCoordinates(Set, nrow, ncol)
    use m_rp_global_var
    !> in/out variables
    integer, intent(in) :: nrow
    integer, intent(in) :: ncol
    real(kind = dbl), intent(inout) :: Set(nrow, ncol)
    !> local variables
    real(kind = dbl) :: windir(nrow)
    real(kind = dbl) :: mag(nrow)
    real(kind = dbl) :: a1(nrow)
    real(kind = dbl) :: a2(nrow)
    real(kind = dbl) :: a3(nrow)
    real(kind = dbl) :: TmpSet(nrow, 2)

    select case (E2Col(u)%instr%wformat)
        case('polar_w')
            !> Polar or axis velocities are converted into u, v
            !> u is wind direction relative to North
            !> v is wind speed (magnitude)
            where(Set(:, u) /= error .and. Set(:, v) /= error)
                windir(:) = Set(:, u)
                mag(:) = Set(:, v)
                Set(:, u) =   mag(:) * dcos(windir(:) * p / 180d0)
                Set(:, v) = - mag(:) * dsin(windir(:) * p / 180d0)
            end where
        case('axis')
            !> Specific to Gill's sonics. See eq. pag. 48 user manual of R3
            !> Conversion is done into U,V,W, apparently converting to AXIS configuration
            !> so wref is forced to AXIS.
            E2Col(u)%instr%wref = 'axis'
            where(Set(:, u) /= error .and. Set(:, v) /= error .and. Set(:, w) /= error)
                a1(:) = Set(:, u)
                a2(:) = Set(:, v)
                a3(:) = Set(:, w)
                Set(:, u) = (2d0 * a1(:) - a2(:) - a3(:)) / 2.1213d0
                Set(:, v) = (a3(:) - a2(:)) / 1.2247d0
                Set(:, w) = (a1 + a2 + a3) / 2.1213d0
            end where
    end select

    !> Compensate for known sonic biases
    where(Set(:, u) /= error .and. Set(:, v) /= error .and. Set(:, w) /= error)
        Set(:, u) = Set(:, u) - RPsetup%offset(u)
        Set(:, v) = Set(:, v) - RPsetup%offset(v)
        Set(:, w) = Set(:, w) - RPsetup%offset(w)
    end where

    !> Directional offset adjustments and conversion to right-handed framework
    select case (E2Col(u)%instr%model(1:len_trim(E2Col(u)%instr%model) - 2))
        case('hs_50', 'hs_100', 'r3_50', 'r3_100', 'r3a_100', 'wm', 'wmpro')
            if (E2Col(u)%instr%wref == 'axis') then
                !> Perform rotation into SPAR configuration, in view of potential
                !> angle of attack correction. Doing the rotation, the north offset
                !> does not need modifying
                where(Set(:, u) /= error .and. Set(:, v) /= error)
                    TmpSet(:, u) = Set(:, u) * dcos(30d0 * p /180d0) - Set(:, v) * dsin(30d0 * p /180d0)
                    TmpSet(:, v) = Set(:, u) * dsin(30d0 * p /180d0) + Set(:, v) * dcos(30d0 * p /180d0)
                    Set(:, u) = TmpSet(:, u)
                    Set(:, v) = TmpSet(:, v)
                end where
            end if
        case('r2')
            !> Perform rotation into SPAR configuration, in view of potential
            !> angle of attack correction. Doing the rotation, the north offset
            !> does not need modifying
            where(Set(:, u) /= error .and. Set(:, v) /= error)
                TmpSet(:, u) = - Set(:, u) * dcos(30d0 * p /180d0) - Set(:, v) * dsin(30d0 * p /180d0)
                TmpSet(:, v) = - Set(:, u) * dsin(30d0 * p /180d0) + Set(:, v) * dcos(30d0 * p /180d0)
                Set(:, u) = TmpSet(:, u)
                Set(:, v) = TmpSet(:, v)
            end where
        case('usa1_standard', 'usa1_fast')
            where(Set(:, v) /= error)
                Set(:, v) = - Set(:, v)
            end where
        case('csat3', 'csat3b')
            E2Col(u)%instr%north_offset = E2Col(u)%instr%north_offset - 180d0
        case('usoni3_cage_mp', 'usoni3_classa_mp')
            E2Col(u)%instr%north_offset = E2Col(u)%instr%north_offset - 90d0
        case('81000', '81000v', '81000re', '81000vre')
            E2Col(u)%instr%north_offset = E2Col(u)%instr%north_offset - 90d0
    end select
end subroutine AdjustSonicCoordinates
