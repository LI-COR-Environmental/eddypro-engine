!***************************************************************************
! kaimal_models.f90
! -----------------
! Copyright (C) 2012-2014, LI-COR Biosciences
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
