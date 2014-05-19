!***************************************************************************
! test_attack_angle.f90
! ---------------------
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
! \brief       Checks for too large wind attack angles \n
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine TestAttackAngle(Set, N)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    real(kind = dbl), intent(inout) :: Set(N, E2NumVar)
    !> local variables
    integer :: i = 0
    integer :: count
    real(kind = dbl) :: aoa = 0.d0
    real(kind = dbl) :: HorVel(N)


    write(*, '(a)', advance = 'no') '   Angle of attack test..'

    !> Calculation of angle of attack for each row
    count = 0
    do i = 1, N
        HorVel(i) = dsqrt((Set(i, U)**2) + (Set(i, V)**2))
        if(HorVel(i) == 0.d0) then
            if (Set(i, W) >= 0.d0) then
                aoa = 90.d0
            end if
            if (Set(i, W) < 0.d0) then
                aoa = - 90.d0
            end if
        else
            aoa = (atan(Set(i, W) / HorVel(i))) * 180.d0 / p
        end if
        if ((aoa < aa%min) .or. (aoa > aa%max)) count = count + 1
    end do
    IntHF%aa = 0
    if (dble(count / N) * 1d2 >= aa%lim) IntHF%aa = 1
    write(*,'(a)') ' done.'
end subroutine TestAttackAngle
