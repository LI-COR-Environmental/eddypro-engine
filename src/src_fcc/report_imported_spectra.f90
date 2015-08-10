!*******************************************************************************
! report_imported_spectra.f90
! ---------------------------
! Copyright (C) 2007-2011, Eco2s team, Gerardo Fratini
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
