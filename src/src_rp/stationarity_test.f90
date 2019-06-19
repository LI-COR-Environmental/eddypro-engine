!***************************************************************************
! stationarity_test.f90
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
! \brief       Assess data quality for stationarity (Foken et al. 2004, HoM)
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine StationarityTest(Set, nrow, ncol, StDiff)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, parameter :: ndiv = 6
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    !> local variables
    integer :: i = 0
    integer :: l = 0
    integer :: j = 0
    integer :: subn
    integer :: IntDiff(GHGNumVar, GHGNumVar)
    integer :: IntDiffUstar
    type(QCType) :: StDiff
    real(kind = dbl) :: LocSet(nrow, GHGNumVar)
    real(kind = dbl), allocatable :: SubSet(:, :)
    real(kind = dbl) :: SubCov(GHGNumVar, GHGNumVar)
    real(kind = dbl) :: AvrgCov(GHGNumVar, GHGNumVar)
    real(kind = dbl) :: GlbCov(GHGNumVar, GHGNumVar)
    real(kind = dbl) :: dev
!    real(kind = dbl) :: Mean(GHGNumVar)
!    real(kind = dbl) :: DirYaw
!    real(kind = dbl) :: Yaw(3, 3)
!    real(kind = dbl) :: SinTheta
!    real(kind = dbl) :: CosTheta
    real(kind = dbl) :: GlbUstar
    real(kind = dbl) :: SubUstar


    write(*, '(a)', advance = 'no') '  Performing stationarity test..'

    !> Define LocSet with only variables u to gas4
    LocSet(:, u:GHGNumVar) = Set(:, u:GHGNumVar)

    !> Global covariances from raw data
    call CovarianceMatrixNoError(LocSet, nrow, GHGNumVar, GlbCov, error)

    !> Partial covariances from subsets and their averages
    subn = int(dble(nrow/ndiv))
    allocate(SubSet(subn, GHGNumVar))
    AvrgCov = 0.d0
    do l = 1, ndiv
        do i = 1, subn
            SubSet(i, :) = LocSet(i + (subn*(l - 1)), :)
        end do
        call CovarianceMatrixNoError(SubSet, subn, GHGNumVar, SubCov, error)
        AvrgCov = AvrgCov + SubCov
    end do
    AvrgCov = AvrgCov / dble(ndiv)

    !> Ustar on the whole covariance and on partial covariances
    SubUstar = (AvrgCov(u, w)**2 + AvrgCov(v, w)**2)**0.25d0
    GlbUstar = (GlbCov(u, w)**2  + GlbCov(v, w)**2 )**0.25d0

    !> Differences
    do i = u, GHGNumVar
        do j = u, GHGNumVar
            if (GlbCov(i, j) /= 0d0 .and. GlbCov(i, j) /= error .and. AvrgCov(i, j) /= error) then
                dev = dabs((GlbCov(i, j) - AvrgCov(i, j)) * 1d2 / GlbCov(i, j))
                if (dabs(dev) < 2147483648.d0) then
                    IntDiff(i, j) = int(dev)
                else
                    IntDiff(i, j) = ierror
                end if
            else
                IntDiff(i, j) = ierror
            end if
        end do
    end do

    if (GlbUstar /= 0d0 .and. GlbUstar /= error .and. SubUstar /= error) then
        dev = dabs((GlbUstar - SubUstar) * 1d2 / GlbUstar)
        if (dabs(dev) < 2147483648.d0) then
            IntDiffUstar = int(dev)
        else
            IntDiffUstar = ierror
        end if
    else
        IntDiffUstar = ierror
    end if

    StDiff%u = IntDiff(u, u)
    StDiff%w = IntDiff(w, w)
    StDiff%ts = IntDiff(ts, ts)
    StDiff%co2 = IntDiff(co2, co2)
    StDiff%h2o = IntDiff(h2o, h2o)
    StDiff%ch4 = IntDiff(ch4, ch4)
    StDiff%gas4 = IntDiff(gas4, gas4)

    StDiff%w_u = IntDiffUstar
    StDiff%w_ts = IntDiff(w, ts)
    StDiff%w_co2 = IntDiff(w, co2)
    StDiff%w_h2o = IntDiff(w, h2o)
    StDiff%w_ch4 = IntDiff(w, ch4)
    StDiff%w_gas4 = IntDiff(w, gas4)

    deallocate(SubSet)
    write(*,'(a)') ' Done.'
end subroutine StationarityTest
