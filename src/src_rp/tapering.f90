!***************************************************************************
! tapering.f90
! ------------
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
! \brief       Defines and applies tapering windows before fft-ing
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine Tapering(Window, xx, N, M, sumw)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    integer, intent(in) :: M
    character(*), intent(in) :: Window
    real(kind = dbl), intent(out) :: sumw
    real(kind = dbl), intent(inout) :: xx(N, M)
    !> local variables
    integer :: i = 0
    integer :: j = 0
    integer :: k
    real(kind = dbl) :: win(N)
    real(kind = dbl) :: sqr_sum
    real(kind = dbl) :: facm, facp
    real(kind = dbl) :: win_han(N)
    real(kind = dbl) :: win_ham(N)
    real(kind = dbl) :: win_bar(N)
    real(kind = dbl) :: win_wel(N)
    real(kind = dbl) :: win_sqr(N)

    write(*, '(a)', advance = 'no') '   Tapering timeseries..'
    !> Selection of window shape
    facm = dble(N)
    facp = 1d0 / dble(N)

    !> Define tapering windows
    do k = 1, N
        win_ham(k) = (0.54d0 + 0.46d0 * dcos(2.d0 * p *(dble(k-1) - facm / 2.d0) * facp))
        win_han(k) = (0.5d0 - 0.5d0 * dcos(2d0 * p * dble(k-1) / (facm - 1d0)))
        win_bar(k) = (1d0 - dabs(((dble(k)-1d0) - facm/2d0) * facp*2d0))
        win_wel(k) = (1d0 - (((dble(k)-1d0) - facm/2d0) * facp*2d0)**2)
        win_sqr(k) = 1d0
    end do

    select case(Window(1:len_trim(Window)))
        case('squared')
            do i = 1, N
                win(i) = win_sqr(i)
            end do
        case('bartlett')
            do i = 1, N
                win(i) = win_bar(i)
            end do
        case('welch')
            do i = 1, N
                win(i) = win_wel(i)
            end do
        case('hamming')
            do i = 1, N
                win(i) = win_ham(i)
            end do
        case('hann')
            do i = 1, N
                win(i) = win_han(i)
            end do
    end select

    !> Apply tapering window to data
    do j = 1, M
        do i = 1, N
            xx(i, j) = xx(i, j) * win(i)
        end do
    end do

    !> Define "window squared and summed" (Numerical Recipes in Fortran, eq. 13.4.11)
    !> that is the normalization factor after FFt-ing and depend on the tapering window
    sqr_sum = 0d0
    do i = 1, N
        sqr_sum = sqr_sum + win(i)**2
    end do
    sumw = sqr_sum * dfloat(N)
    write(*,'(a)') ' Done.'
end subroutine Tapering
