!***************************************************************************
! binned_frequency_vector_fcc.f90
! -------------------------------
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
! \brief       Define natural binned frequency vector
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine BinnedFrequencyVector(lEx, bnf, N)
    use m_fx_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    type(ExType), intent(in) :: lEx
    real(kind = dbl), intent(out) :: bnf(N + 1)
    !> local variables
    integer :: PeriodSeconds
    integer :: i

    !> Exponentially spaced frequencies extending
    !> from f_min = 1/FileLength Hz --> f_max = AcFreq/2 (= Nyquist frequency)
    call log_msg(' inf=creating binned-frequencies array.')
    if (lEx%ac_freq > 0d0) then
        PeriodSeconds = nint(lEx%avrg_length) * 60
        bnf(1) = 1d0 / dfloat(PeriodSeconds)
        bnf(N + 1) = lEx%ac_freq / 2d0
        do i = 2, N
            bnf(i) = bnf(1) * dexp(dble(i - 1) * (dlog(bnf(N + 1)) - dlog(bnf(1))) / dble(N))
        end do
    else
        call log_msg(' err=acquisition frequency seems to be <= 0. &
                     &low-pass spectral correction not applied.')
        call ErrorHandle(2, 0, 5)
        EddyProProj%hf_meth = 'none'
    end if
    call log_msg(' inf=array created correctly.')
end subroutine BinnedFrequencyVector
