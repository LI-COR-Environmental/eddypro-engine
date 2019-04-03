!***************************************************************************
! add_to_cospectra_fit_dataset.f90
! --------------------------------
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
! \brief       Add current cospectra to dataset for model fitting
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine AddToCospectraFitDataset(lEx, BinCospForStable, BinCospForUnstable, snrow, &
        nfit, nrow, ncol, nbins, FitStable, FitUnstable, fnrow)
    use m_fx_global_var
    implicit none
    !> in/out variables
    type(ExType), intent(in) :: lEx
    integer, intent(in) :: fnrow
    integer, intent(in) :: snrow
    integer, intent(in) :: nrow
    integer, intent(in) :: ncol
    integer, intent(in) :: nbins
    type(SpectraSetType) :: BinCospForStable(snrow)
    type(SpectraSetType) :: BinCospForUnstable(snrow)
    type(FitSpectraType), intent(out) :: FitUnstable(fnrow)
    type(FitSpectraType), intent(out) :: FitStable(fnrow)

    integer, intent(inout) :: nfit(nrow, ncol)
    !> local variables
    integer :: var
    integer :: i
    integer :: sort


    !> Sort by stability regime
    !> Stability conditions as from Kljun et al. (2004)
    sort = 0
    !> Unstable conditions
    if (lEx%L > -650 .and. lEx%L < 0) then
        sort = 1
    !> Stable conditions
    elseif (lEx%L >= 0 .and. lEx%L < 1000) then
        sort = 2
    end if
    if (sort == 0) return


    if (sort == 1) then
        do var = w_ts, w_gas4
            do i = 1, nbins
                if (BinCospForUnstable(i)%fnorm /= error .and. BinCospForUnstable(i)%of(var) > 0d0) then
                    nfit(var, unstable) = nfit(var, unstable) + 1
                    FitUnstable(nfit(var, unstable))%fnorm(var) = BinCospForUnstable(i)%fnorm !< normalized frequency on x
                    FitUnstable(nfit(var, unstable))%of(var) = BinCospForUnstable(i)%fn * BinCospForUnstable(i)%of(var)   !< natural freq * cospectr / covar on y
                endif
            end do
        end do
    elseif(sort == 2) then
        do var = w_ts, w_gas4
            do i = 1, nbins
                if (BinCospForStable(i)%fnorm /= error .and. BinCospForStable(i)%of(var) > 0d0) then
                    nfit(var, stable) = nfit(var, stable) + 1
                    FitStable(nfit(var, stable))%fnorm(var) = BinCospForStable(i)%fnorm !< normalized frequency on x
                    FitStable(nfit(var, stable))%of(var) = BinCospForStable(i)%fn * BinCospForStable(i)%of(var)   !< natural freq * cospectr / covar on y
                endif
            end do
        end do
    end if
end subroutine AddToCospectraFitDataset
