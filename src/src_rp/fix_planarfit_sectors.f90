!***************************************************************************
! fix_planarfit_sectors.f90
! -------------------------
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
! \brief       Replaces error planar fits with closest valid sector rotation matrix
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine FixPlanarfitSectors(GoPlanarFit, N)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    logical, intent(inout) :: GoPlanarFit(N)
    !> Local variables
    integer :: sec
    integer :: sec2
    real(kind = dbl) :: loc_pfmat(3, 3, -N + 1: 2 * N)
    logical :: loc_go(-N + 1: 2 * N)


    !> First, if there is no valid sector, switches to 2D rotations
    do sec = 1, N
        if (GoPlanarFit(sec)) exit
    end do
    if (sec == N + 1) then
        Meth%rot = 'double_rotation'
        call ExceptionHandler(37)
        return
    end if

    !> Define working arrays, only for convenience, for going
    !> clockwise or counter-clockwise easier
    !>
    !> N-1           0   1            N   N+1          2*N
    !>  |--o--x--*---|   |--o--x--*---|    |--o--x--*---|
    !>
    loc_pfmat(:, :, -N + 1: 0)    = PFMat(:, :, 1:N)
    loc_pfmat(:, :, 1: N)         = PFMat(:, :, 1:N)
    loc_pfmat(:, :, N + 1: 2 * N) = PFMat(:, :, 1:N)
    loc_go(-N + 1: 0)    = GoPlanarFit(1:N)
    loc_go(1: N)         = GoPlanarFit(1:N)
    loc_go(N + 1: 2 * N) = GoPlanarFit(1:N)

    do sec = 1, N
        if (.not. GoPlanarFit(sec)) then
            if(PFSetup%fix == 'clockwise') then
                !> Searches clockwise
                do sec2 = sec + 1, 2*N
                    if (loc_go(sec2)) then
                        PFMat(:, :, sec) = loc_pfmat(:, :, sec2)
                        GoPlanarFit(Sec) = .true.
                        exit
                    end if
                end do
            elseif(PFSetup%fix == 'counterclockwise') then
                !> Searches counterclockwise
                do sec2 = sec - 1, - N + 1
                    if (loc_go(sec2)) then
                        PFMat(:, :, sec) = loc_pfmat(:, :, sec2)
                        GoPlanarFit(sec) = .true.
                        exit
                    end if
                end do
            end if
        end if
    end do
end subroutine FixPlanarfitSectors
