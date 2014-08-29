!***************************************************************************
! binned_frequency_vector.f90
! ---------------------------
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
! \brief       Define binned natural frequency array
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine BinnedFrequencyVector(binned_freq, N, flux_avrg_length, ac_freq)
    use m_rp_global_var
    implicit none
    ! in/out variables
    integer, intent(in) :: N
    integer, intent(in) :: flux_avrg_length
    real(kind = dbl), intent(out) :: binned_freq(N + 1)
    real(kind = dbl), intent(out) :: ac_freq
    ! local variables
    integer :: PeriodSeconds
    integer :: i

    !> Exponentially spaced frequencies extending
    !> from f_min = 1/FileLength Hz --> f_max = AcFreq/2 (= Nyquist frequency)
    if (ac_freq > 0d0) then
        PeriodSeconds = flux_avrg_length * 60
        binned_freq(1) = 1d0 / dfloat(PeriodSeconds)
        binned_freq(N + 1) = ac_freq / 2d0
        do i = 2, N
            binned_freq(i) = binned_freq(1) * dexp(dble(i - 1) * &
                    (dlog(binned_freq(N + 1)) - dlog(binned_freq(1))) / dble(N))
        end do
    else
        call ExceptionHandler(66)
        Meth%spec%lptf = 'none'
    end if
end subroutine BinnedFrequencyVector

