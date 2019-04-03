!***************************************************************************
! tests_not_performed.f90
! -----------------------
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
! \brief       Sets flags to 99999... for tests not performed
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine TestsNotPerformed()
    use m_rp_global_var
    implicit none

    if(.not.Test%sr) IntHF%sr = 99999999
    if(.not.Test%ar) IntHF%ar = 99999999
    if(.not.Test%do) IntHF%do = 99999999
    if(.not.Test%al) IntHF%al = 99999999
    if(.not.Test%sk) IntHF%sk = 99999999
    if(.not.Test%sk) IntSF%sk = 99999999
    if(.not.Test%ds) IntHF%ds = 99999999
    if(.not.Test%ds) IntSF%ds = 99999999
    if(.not.Test%tl) IntHF%tl = 9999
    if(.not.Test%tl) IntSF%tl = 9999
    if(.not.Test%aa) IntHF%aa = 9
    if(.not.Test%ns) IntHF%ns = 9
end subroutine TestsNotPerformed
