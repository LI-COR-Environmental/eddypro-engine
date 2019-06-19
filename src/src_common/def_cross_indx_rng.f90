!***************************************************************************
! def_cross_indx_rng.f90
! ----------------------
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
! \brief
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine DefCrossIndxRng(N, RowLags1, RowLags2, Nmin, Nmax, Nred)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    integer, intent(in) :: RowLags1
    integer, intent(in) :: RowLags2
    integer, intent(out) :: Nmin
    integer, intent(out) :: Nmax
    integer, intent(out) :: Nred

    if (RowLags1 >= 0 .and. RowLags2 >= 0) then
        Nred = N - max(RowLags1, RowLags2)
        Nmin = 1
        Nmax = Nred
    elseif (RowLags1 < 0 .and. RowLags2 >= 0) then
        Nred = N - (RowLags2 - RowLags1)
        Nmin = 1 - RowLags1
        Nmax = N - RowLags2
    elseif (RowLags1 >= 0 .and. RowLags2 < 0) then
        Nred = N - (RowLags1 - RowLags2)
        Nmin = 1 - RowLags2
        Nmax = N - RowLags1
    elseif (RowLags1 < 0 .and. RowLags2 < 0) then
        Nred = N - abs(min(RowLags1, RowLags2))
        Nmin = 1 + abs(min(RowLags1, RowLags2))
        Nmax = N
    end if
end subroutine DefCrossIndxRng
