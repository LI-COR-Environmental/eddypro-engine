!***************************************************************************
! ensemble_cospectra_by_stability.f90
! -----------------------------------
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
! \brief       Ensemble average cospectra in stability classes
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine EnsembleCospectraByStability(nfit, dim1, dim2, FitStable, FitUnstable, fnrow)
    use m_fx_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: dim1
    integer, intent(in) :: dim2
    integer, intent(in) :: fnrow
    integer, intent(in) :: nfit(dim1, dim2)
    type(FitSpectraType), intent(in) :: FitUnstable(fnrow)
    type(FitSpectraType), intent(in) :: FitStable(fnrow)
    !> local variables
    integer :: var
    integer :: i
    integer :: j
    integer :: m

    write(*, '(a)', advance = 'no') ' Sorting cospectra in stability classes.. '
    MeanStabilityCosp = NullMeanSpec
    !> Unstable case
    do var = ts, gas4
        m = nfit(var, unstable)
        do j = 1, ndkf
            do i = 1, m
                if (FitUnstable(i)%fnorm(var) > dkf(j) .and. FitUnstable(i)%fnorm(var) < dkf(j+1)) then
                    MeanStabilityCosp(j, unstable)%cnt(var) = MeanStabilityCosp(j, unstable)%cnt(var) + 1
                    MeanStabilityCosp(j, unstable)%fn(var)  = MeanStabilityCosp(j, unstable)%fn(var)  + FitUnstable(i)%fnorm(var)
                    MeanStabilityCosp(j, unstable)%of(var)  = MeanStabilityCosp(j, unstable)%of(var)  + dlog(FitUnstable(i)%of(var))
                end if
            end do
            if (MeanStabilityCosp(j, unstable)%cnt(var) > 0) then
                MeanStabilityCosp(j, unstable)%of(var) = MeanStabilityCosp(j, unstable)%of(var) &
                                                       / MeanStabilityCosp(j, unstable)%cnt(var)
                MeanStabilityCosp(j, unstable)%of(var) = dexp(MeanStabilityCosp(j, unstable)%of(var))
                MeanStabilityCosp(j, unstable)%fn(var) = MeanStabilityCosp(j, unstable)%fn(var) &
                                                       / MeanStabilityCosp(j, unstable)%cnt(var)
            end if
        end do
    end do

    !> Stable case
    do var = ts, gas4
        m = nfit(var, stable)
        do j = 1, ndkf
            do i = 1, m
                if (FitStable(i)%fnorm(var) > dkf(j) .and. FitStable(i)%fnorm(var) < dkf(j+1)) then
                    MeanStabilityCosp(j, stable)%cnt(var) = MeanStabilityCosp(j, stable)%cnt(var) + 1
                    MeanStabilityCosp(j, stable)%fn(var)  = MeanStabilityCosp(j, stable)%fn(var)  + FitStable(i)%fnorm(var)
                    MeanStabilityCosp(j, stable)%of(var)  = MeanStabilityCosp(j, stable)%of(var)  + FitStable(i)%of(var)
                end if
            end do
            if (MeanStabilityCosp(j, stable)%cnt(var) > 0) then
                MeanStabilityCosp(j, stable)%of(var) = MeanStabilityCosp(j, stable)%of(var) &
                                                       / MeanStabilityCosp(j, stable)%cnt(var)
                MeanStabilityCosp(j, stable)%fn(var) = MeanStabilityCosp(j, stable)%fn(var) &
                                                       / MeanStabilityCosp(j, stable)%cnt(var)
            end if
        end do
    end do
    write(*, '(a)') ' Done.'
end subroutine EnsembleCospectraByStability

