!***************************************************************************
! fisher.f90
! ----------
! Copyright (C) 2018-2019, LI-COR Biosciences, Inc.  All Rights Reserved.
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
! \brief       Compute Fisher test on covariances with and without repeated 
!              values
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine fisher(Set, nrow, ncol)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    !> Local variables
    integer :: var
    integer :: i
    integer :: j
    real(kind = dbl) :: Corr(ncol, ncol)
    real(kind = dbl) :: mSet(nrow, ncol)
    real(kind = dbl) :: mCorr(ncol, ncol)


    !> Compute mSet by eliminating repeated values in Set
    write(*, '(a)', advance = 'no') "  Evaluating correlation differences with and without repeated values.."
    mSet = Set
    do var = u, gas4
        do i = 2, nrow
            if ((Set(i, var) - Set(i-1, var)) < 1d-6) mSet(i, var) = error
        end do
    end do
    call CorrelationMatrixNoError(Set, size(Set, 1), size(Set, 2), Corr, error)
    call CorrelationMatrixNoError(mSet, size(mSet, 1), size(mSet, 2), mCorr, error)
    do i = u, gas4
        do j = u, gas4
            if (Corr(i, j) /= error .and. mCorr(i,j) /= error) then
                Essentials%CorrDiff(i, j) = dabs(Corr(i, j) - mCorr(i, j))
            else
                Essentials%CorrDiff(i, j) = error
            end if
        end do
    end do
    write(*, *) " Done."

end subroutine fisher