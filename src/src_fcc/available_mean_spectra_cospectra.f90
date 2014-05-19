!***************************************************************************
! available_mean_spectra_cospectra.f90
! ------------------------------------
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
! \brief       Detect which average spectra/cospectra are available for each gas and each class
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine AvailableMeanSpectraCospectra(nbins)
    use m_fx_global_var
    implicit none
    !> In/out variables
    integer, intent(in) :: nbins
    !> local variables
    integer :: gas
    integer :: cls
    integer :: RH
    integer :: bin
    integer :: err_cnt


    MeanBinSpecAvailable = .true.
    MeanBinCospAvailable = .true.
    MeanStabCospAvailable = .true.

    !> Spectra
    do gas = co2, gas4
        if (gas /= h2o) then
            do cls = 1, MaxGasClasses
                err_cnt = 0
                do bin = 1, nbins
                    if (MeanBinSpec(bin, cls)%fnum(gas) == error) err_cnt = err_cnt + 1
                end do
                if (err_cnt == nbins) MeanBinSpecAvailable(cls, gas) = .false.
            end do
        else
            do RH = RH10, RH90
                err_cnt = 0
                do bin = 1, nbins
                    if (MeanBinSpec(bin, RH)%fnum(gas) == error) err_cnt = err_cnt + 1
                end do
                if (err_cnt == nbins) MeanBinSpecAvailable(RH, gas) = .false.
            end do
        end if
    end do

    !> Cospectra
    do gas = ts, gas4
        do cls = 1, 8
            err_cnt = 0
            do bin = 1, nbins
                if (MeanBinCosp(bin, cls)%fnum(gas) == error) err_cnt = err_cnt + 1
            end do
            if (err_cnt == nbins) MeanBinCospAvailable(cls, gas) = .false.
        end do
    end do

    !> Stability cospectra
    do gas = ts, gas4
        do cls = unstable, stable
            err_cnt = 0
            do bin = 1, ndkf
                if (MeanStabilityCosp(bin, cls)%cnt(gas) == 0) err_cnt = err_cnt + 1
            end do
            if (err_cnt == ndkf) MeanStabCospAvailable(cls, gas) = .false.
        end do
    end do

end subroutine AvailableMeanSpectraCospectra
