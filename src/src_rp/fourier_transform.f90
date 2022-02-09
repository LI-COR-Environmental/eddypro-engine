!***************************************************************************
! fourier_transform.f90
! ---------------------
! Copyright (C) 2007-2011, Eco2s team, Gerardo Fratini
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
! \brief       Provides a driver to rfftf. Fourier-transform the in/out array "xx".
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine FourierTransform(xx, N, M)

    use m_rp_global_var
    use fftpack, only: dffti, dfftf

    implicit none

    !> in/out variables
    integer, intent(in) :: N
    integer, intent(in) :: M
    real(kind = dbl), intent(inout) :: xx(N, M)

    !> local variables
    integer :: i
    real(dbl) :: xxx(N)
    real(dbl) :: wsave(N*2 + 15)

    write(*, '(a)', advance = 'no') '   FFT-ing..'
    ! call rffti(N, wsave)
    call dffti(N, wsave)
    do i = 1, M
        !> data in 1D vector
        ! xxx(:) = sngl(xx(:, i))
        xxx(:) = xx(:, i)
        !> fast fourier transform
        ! call rfftf(N, xxx, wsave)
        call dfftf(N, xxx, wsave)
        !> replace time data with spectral data
        ! xx(:, i) = dble(xxx(:))
        xx(:, i) = xxx(:)
    end do
    write(*,'(a)') ' Done.'

end subroutine FourierTransform
