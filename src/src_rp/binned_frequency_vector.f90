!***************************************************************************
! binned_frequency_vector.f90
! ---------------------------
! Copyright (C) 2011-2019, LI-COR Biosciences, Inc.  All Rights Reserved.
! Author: Gerardo Fratini
!
! This file is part of EddyPro®.
!
! NON-COMMERCIAL RESEARCH PURPOSES ONLY - EDDYPRO® is licensed for 
! non-commercial academic and government research purposes only, 
! as provided in the EDDYPRO® End User License Agreement. 
! EDDYPRO® may only be used as provided in the End User License Agreement
! and may not be used or accessed for any commercial purposes.
! You may view a copy of the End User License Agreement in the file
! EULA_NON_COMMERCIAL.rtf.
!
! Commercial companies that are LI-COR flux system customers 
! are encouraged to contact LI-COR directly for our commercial 
! EDDYPRO® End User License Agreement.
!
! EDDYPRO® contains Open Source Components (as defined in the 
! End User License Agreement). The licenses and/or notices for the 
! Open Source Components can be found in the file LIBRARIES-ENGINE.txt.
!
! EddyPro® is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
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

