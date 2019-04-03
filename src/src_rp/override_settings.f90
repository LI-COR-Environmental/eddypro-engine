!***************************************************************************
! override_settings.f90
! ---------------------
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
! \brief       Forces some operations (regardless of user choice) based on instrument
!              models (e.g. CSAT3 no cross-wind correction) and logic
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine OverrideSettings()
    use m_rp_global_var
    implicit none

    !> If biomet measurements are not to be used, they are also not to be output
    if (EddyProProj%biomet_data == 'none') EddyProProj%out_biomet = .false.

    !> if there is no LI-7500 among the instruments, Burba terms should not be calculated
    if (index(E2Col(co2)%Instr%model, 'li7500') == 0 &
        .and. index(E2Col(h2o)%Instr%model,'li7500') == 0) &
        RPsetup%bu_corr = 'none'
end subroutine OverrideSettings
