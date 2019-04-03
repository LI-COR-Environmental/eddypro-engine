!*******************************************************************************
! report_imported_spectra.f90
! ---------------------------
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
!*******************************************************************************
!
! \brief       Output on stdout number of spectra and co-spectra imported
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!*******************************************************************************
subroutine ReportImportedSpectra(nbins)
    use m_fx_global_var
    implicit none
    !> In/out variables
    integer, intent(in) :: nbins
    !> Local variables
    integer :: nh2o

    write(*, '(a)') '  Imported gas binned spectra:'
    nh2o = sum(MeanBinSpec(nbins/2, RH10:RH90)%cnt(h2o))
    write(*, '(a, i5)')  '   CO2:   ', MeanBinSpec(nbins/2, 1)%cnt(co2)
    write(*, '(a, i5)')  '   H2O:   ', nh2o
    write(*, '(a, i5)')  '   CH4:   ', MeanBinSpec(nbins/2, 1)%cnt(ch4)
    write(*, '(a, i5)')  '   Gas 4: ', MeanBinSpec(nbins/2, 1)%cnt(gas4)
    write(*, '(a)') '  Done.'
end subroutine ReportImportedSpectra
