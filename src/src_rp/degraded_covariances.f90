!***************************************************************************
! degraded_covariances.f90
! ------------------------
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
! \brief       Calculate degraded temperature timeseries, by filtering in
!              the frequency domain.
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine DegradedCovariances(nf, Cospectrum, N)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    real(kind = dbl), intent(in) :: nf(N/2)
    type (SpectralType), intent(in) :: Cospectrum(N/2 + 1)
    !> local variables
    integer :: i
    integer :: j
    real(kind = dbl) :: cov
    real(kind = dbl) :: iir(N/2)
    real(kind = dbl) :: f_co(NumDegH)
    data f_co(1:NumDegH) &
        /1.626d0, 0.614d0, 0.277d0, 0.133d0, 6.5d-2, 3.2d-2, 1.6d-2, 8d-3, 4d-3/


    !> covariance from unfiltered cospectrum
    cov = 0d0
    do j = 1, N/2
        cov = cov + Cospectrum(j)%of(w_ts)
    end do
    Essentials%degH(NumDegH + 1) = cov / (dble(N) / Metadata%ac_freq)

    !> degraded covariances
    do i = 1, NumDegH
        iir(:) = dsqrt(1d0 / (1d0 + (nf(:)/f_co(i))**2))
        Essentials%degH(i) = 0d0
        do j = 1, N/2
            Essentials%degH(i) = Essentials%degH(i) + Cospectrum(j)%of(w_ts) * iir(j)
        end do
        Essentials%degH(i) = Essentials%degH(i) / (dble(N)/Metadata%ac_freq)
    end do
end subroutine DegradedCovariances
