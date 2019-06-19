!***************************************************************************
! assess_daytime.f90
! ------------------
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
! \brief       Assess whether it is daytime or night-time, based on
!              global radiation, PPFD or potential radiation
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine AssessDaytime(date, time)
    use m_rp_global_var
    implicit none
    !> in/out variables
    character(10), intent(in) :: date
    character(5), intent(in) :: time
    logical, external :: IsDaytime

    Stats%daytime = .false.
    !> Based on Rg
    if (biomet%val(bRg) > 12d0) Stats%daytime = .true.
    !> Based on PPFD
    if (.not. Stats%daytime .and. biomet%val(bPPFD) > 100d0) &
        Stats%daytime = .true.
    !> based on period timestamp and potential radiation
    if (biomet%val(bRg) == error .and. biomet%val(bPPFD) == error) &
        Stats%daytime = IsDaytime(PotRad, date, time)
end subroutine AssessDaytime
