!***************************************************************************
! available_mean_spectra_cospectra.f90
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


    MeanBinSpecAvailable = .true.
    MeanBinCospAvailable = .true.
    MeanStabCospAvailable = .true.

    !> Spectra
    do gas = co2, gas4
        do cls = 1, MaxGasClasses
            if (all(MeanBinSpec(1:nbins, cls)%fnum(gas) == error) &
                .or. all(MeanBinSpec(1:nbins, cls)%fnum(gas) == 0d0)) &
                MeanBinSpecAvailable(cls, gas) = .false.
        end do
    end do

    !> Cospectra
    do gas = ts, gas4
        do cls = 1, 8
            if (all(MeanBinCosp(1:nbins, cls)%fnum(gas) == error) &
                .or. all(MeanBinCosp(1:nbins, cls)%fnum(gas) == 0d0)) then
                MeanBinCospAvailable(cls, gas) = .false.
            else
            end if
        end do
    end do

    !> Stability cospectra
    do gas = ts, gas4
        do cls = unstable, stable
                if (all(MeanStabilityCosp(1:ndkf, cls)%cnt(gas) == error) &
                    .or. all(MeanStabilityCosp(1:ndkf, cls)%cnt(gas) == 0)) &
                    MeanStabCospAvailable(cls, gas) = .false.
        end do
    end do
end subroutine AvailableMeanSpectraCospectra
