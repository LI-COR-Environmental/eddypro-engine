!***************************************************************************
! normalize_mean_spectra_cospectra.f90
! ------------------------------------
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
