!***************************************************************************
! ts_round_to_minute.f90
! ----------------------
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
! \brief       Rounds the given timestamp to the provided "precision"
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine tsRoundToMinute(tstamp, approx, where_to)
    use m_common_global_var
    implicit none
    !> In/out variables
    integer, intent(in) :: approx
    character(*), intent(in) :: where_to
    type (DateType), intent(inout) :: tstamp
    !> Local variables
    integer :: base
    integer :: off


    base = approx * (tstamp%minute / approx)
    off  = tstamp%minute - base
    if (off == 0) return

    if (where_to == 'later') then
        tstamp = tstamp + datetype(0, 0, 0, 0, approx - off)
    elseif (where_to == 'earlier') then
        tstamp = tstamp - datetype(0, 0, 0, 0, off)
    elseif(where_to == 'closest') then
        if (off >= approx / 2) then
            tstamp = tstamp + datetype(0, 0, 0, 0, approx - off)
        else
            tstamp = tstamp - datetype(0, 0, 0, 0, off)
        end if
    end if
end subroutine
