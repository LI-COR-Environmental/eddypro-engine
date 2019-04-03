!***************************************************************************
! generate_t_cell.f90
! -------------------
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
!***************************************************************************
!
! \brief       Generate LI-7200 cell temperature timeseries based on availability
!              Tcell = Tcell if Tcell is available
!              Tcell = 0.2 * Tin + 0.8 * Tout if Tin and Tout available
!              Tcell = Tin or Tout, if either is missing
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine GenerateTcell(Set, nrow, ncol)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    integer, intent(in) :: ncol
    real(kind = dbl), intent(inout) :: Set(nrow, ncol)
    !> local variables


    if(.not. E2Col(tc)%present) then
        if (E2Col(ti1)%present .and. E2Col(ti2)%present) then
            E2Col(tc) = E2Col(ti1)
            E2Col(tc)%var = 'tc'
            where (Set(1:nrow, ti1) /= error .and. Set(1:nrow, ti2) /= error)
                Set(1:nrow, tc) = Set(1:nrow, ti1) * 0.2d0 + Set(1:nrow, ti2) * 0.8d0
            else where (Set(1:nrow, ti1) /= error)
                Set(1:nrow, tc) = Set(1:nrow, ti1)
            else where (Set(1:nrow, ti2) /= error)
                Set(1:nrow, tc) = Set(1:nrow, ti2)
            else where
                Set(1:nrow, tc) = error
            end where

        elseif(E2Col(ti1)%present) then
            E2Col(tc) = E2Col(ti1)
            E2Col(tc)%var = 'tc'
            Set(1:nrow, tc) = Set(1:nrow, ti1)

        elseif(E2Col(ti2)%present) then
            E2Col(tc) = E2Col(ti2)
            E2Col(tc)%var = 'tc'
            Set(1:nrow, tc) = Set(1:nrow, ti2)
        end if
    end if
    E2Col(ti1)%present = .false.
    E2Col(ti2)%present = .false.
end subroutine GenerateTcell
