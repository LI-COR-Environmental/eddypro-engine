!***************************************************************************
! kaimal_models.f90
! -----------------
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
! \brief       Returns Kaimal cospectra as a function of stability \n
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
double precision function kaimal(fnorm, zL, stability)
    use m_common_global_var
    !> in/out variables
    real(kind = dbl), intent(in) :: fnorm
    real(kind = dbl), intent(in) :: zL
    character(*) :: stability
    !> local and returned variables
    real(kind = dbl) :: Ak
    real(kind = dbl) :: Bk

    if (fnorm == 0d0 .or. fnorm == error) then
        kaimal = error
        return
    end if

    if (stability == 'unstable') then
        if (fnorm <= 0.54d0) then
            kaimal = 12.92d0 * fnorm / ((1.d0 + 26.7d0 * fnorm)**1.375d0)
        else
            kaimal = 4.378d0 * fnorm / ((1.d0 +  3.8d0 * fnorm)**2.4d0)
        end if
    elseif (stability == 'stable') then
        Ak = 0.284d0 * ((1.d0 + 6.4d0 * zL)**0.75d0)
        Bk = 2.34d0  * (Ak**(-1.1d0))
        kaimal = fnorm / (Ak + Bk * fnorm**2.1d0)
    else
        kaimal = error
    end if
end function kaimal
