!***************************************************************************
! add_to_cospectra_fit_dataset.f90
! --------------------------------
! Copyright (C) 2011-2015, LI-COR Biosciences
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
