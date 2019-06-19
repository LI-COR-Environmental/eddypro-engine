!***************************************************************************
! int_2_flag.f90
! --------------
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
! \brief       Converts integer flags into characher flags
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine Int2Flags(len)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: len

    call int2char(IntHF%sr, CharHF%sr, len)
    call int2char(IntHF%ar, CharHF%ar, len)
    call int2char(IntHF%do, CharHF%do, len)
    call int2char(IntHF%al, CharHF%al, len)
    call int2char(IntHF%sk, CharHF%sk, len)
    call int2char(IntSF%sk, CharSF%sk, len)
    call int2char(IntHF%ds, CharHF%ds, len)
    call int2char(IntSF%ds, CharSF%ds, len)
    call int2char(IntHF%tl, CharHF%tl, len)
    call int2char(IntSF%tl, CharSF%tl, len)
    call int2char(IntHF%aa, CharHF%aa, len)
    call int2char(IntHF%ns, CharHF%ns, len)
end subroutine Int2Flags
