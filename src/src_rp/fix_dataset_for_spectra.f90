!***************************************************************************
! fix_dataset_for_spectra.f90
! ---------------------------
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
! \brief       replace gaps with linear interpolation of neighboring data
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine FixDatasetForSpectra(Set, nrow, ncol, nrow2)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    integer, intent(out) :: nrow2
    real(kind = dbl) :: Set(nrow, ncol)
    !> local variables
    integer :: j
    integer :: tnrow


    !> If more than 30% of the data is missing, don't compute spectra
    !> because linear interpolation probably too severly affect spectral shape
    !> This filter is totally arbitrary, only based on anecdotal evidence
    do j = 1, ncol
        if (count(Set(1:nrow, j) == error) > nrow / 3) SpecCol(j)%present = .false.
    end do

    !> nrow2 is the smallest nrow of all columns
    nrow2 = nrow
    do j = 1, GHGNumVar
        if (SpecCol(j)%present) then
            call ReplaceGapWithLinearInterpolation(Set(1:nrow, j), size(Set, 1), tnrow, error)
            if (tnrow < nrow2) nrow2 = tnrow
        end if
    end do
end subroutine FixDatasetForSpectra
