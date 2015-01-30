!***************************************************************************
! normalize_mean_spectra_cospectra.f90
! ------------------------------------
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
! \brief       Normalize sums to obtain average spectra and cospectra
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine NormalizeMeanSpectraCospectra(nbins)
    use m_fx_global_var
    implicit none
    !> In/out variables
    integer, intent(in) :: nbins
    !> Local variables
    integer :: bin
    integer :: cls


    !> Normalize to complete the averaging process
    do cls = 1, MaxGasClasses
        do bin = 1, nbins
            !> Sorted spectra
            where (MeanBinSpec(bin, cls)%cnt(co2:gas4) >= FCCsetup%SA%min_smpl)
                MeanBinSpec(bin, cls)%fnum(co2:gas4) = &
                    MeanBinSpec(bin, cls)%fnum(co2:gas4) / MeanBinSpec(bin, cls)%cnt(co2:gas4)
                MeanBinSpec(bin, cls)%fn(co2:gas4) = &
                    MeanBinSpec(bin, cls)%fn(co2:gas4)   / dfloat(MeanBinSpec(bin, cls)%cnt(co2:gas4))
                MeanBinSpec(bin, cls)%of(co2:gas4) = &
                    MeanBinSpec(bin, cls)%of(co2:gas4)   / dfloat(MeanBinSpec(bin, cls)%cnt(co2:gas4))
                MeanBinSpec(bin, cls)%ts(co2:gas4) = &
                    MeanBinSpec(bin, cls)%ts(co2:gas4)   / dfloat(MeanBinSpec(bin, cls)%cnt(co2:gas4))
            elsewhere
                MeanBinSpec(bin, cls)%fnum(co2:gas4) = nint(error)
                MeanBinSpec(bin, cls)%fn(co2:gas4) = error
                MeanBinSpec(bin, cls)%of(co2:gas4) = error
                MeanBinSpec(bin, cls)%ts(co2:gas4) = error
            end where
            !> Sorted cospectra
            where (MeanBinCosp(bin, cls)%cnt(ts:gas4) > 0)
                MeanBinCosp(bin, cls)%fnum(ts:gas4) = &
                    MeanBinCosp(bin, cls)%fnum(ts:gas4) / MeanBinCosp(bin, cls)%cnt(ts:gas4)
                MeanBinCosp(bin, cls)%fn(ts:gas4) = &
                    MeanBinCosp(bin, cls)%fn(ts:gas4)   / dfloat(MeanBinCosp(bin, cls)%cnt(ts:gas4))
                MeanBinCosp(bin, cls)%of(ts:gas4) = &
                    MeanBinCosp(bin, cls)%of(ts:gas4)   / dfloat(MeanBinCosp(bin, cls)%cnt(ts:gas4))
            elsewhere
                MeanBinCosp(bin, cls)%fnum(ts:gas4) = nint(error)
                MeanBinCosp(bin, cls)%fn(ts:gas4) = error
                MeanBinCosp(bin, cls)%of(ts:gas4) = error
            end where
        end do
    end do
end subroutine NormalizeMeanSpectraCospectra
