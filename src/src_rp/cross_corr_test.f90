!***************************************************************************
! cross_corr_test.f90
! -------------------
! Copyright (C) 2018-2020, LI-COR Biosciences, Inc.  All Rights Reserved.
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
! \brief       Compute Cross-Correlation test on time series with and
!              without repeated values, see Vitale et al. (2020), BGD
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine CrossCorrTest(Set, nrow, ncol)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    !> Local variables
    integer :: var
    integer :: i
    integer, parameter :: lagmin = -25
    integer, parameter :: lagmax = 25
    real(kind = dbl) :: mSet(nrow, ncol)
    real(kind = dbl) :: CCF(lagmin: lagmax)
    real(kind = dbl) :: mCCF(lagmin: lagmax)
    real(kind = dbl) :: r
    real(kind = dbl) :: t
    integer(kind = dbl) :: m
    real(kind = dbl) :: cov, sig(1), msig(1)
    real(kind = dbl), external :: LaggedCovarianceNoError

    !> Compute mSet by eliminating repeated values in Set
    write(*, '(a)', advance = 'no') "  Evaluating R2 on CCFs with and without repeated values.."

    !> Set repeated values to error
    mSet = Set
    do var = u, gas4
        if (OutVarPresent(var)) then
            do i = 1, nrow - 1
                if (abs((Set(i + 1, var) - Set(i, var))) < 1d-6) mSet(i, var) = error
            end do
        end if
    end do

    !> Estimate CCFs
    open(123, file='all.txt')
    open(124, file='non_repeated.txt')
    do i = 1, nrow
    write(123, *) Set(i, u:gas4)
    write(124, *) mSet(i, u:gas4)
    end do

    do var = ts, gas4
        if (OutVarPresent(var)) then
            call CrossCorrelation(Set(:, w), Set(:, var), nrow, lagmin, lagmax, CCF)
            call CrossCorrelation(mSet(:, w), mSet(:, var), nrow, lagmin, lagmax, mCCF)

            cov = LaggedCovarianceNoError(CCF, mCCF, size(CCF), 0, error)
            call StDevNoError(CCF, size(CCF), 1, sig, error)
            call StDevNoError(mCCF, size(mCCF), 1, msig, error)
            r = cov / (sig(1) * msig(1))

            ! call unbiased_correlation(CCF, mCCF, lagmax - lagmin + 1, error, 0, r, t, m)
        end if
    end do


    write(*, *) " Done."

end subroutine CrossCorrTest