!***************************************************************************
! burba_terms.f90
! ---------------
! Copyright (C) 2007-2011, Eco2s team, Gerardo Fratini
! Copyright (C) 2011-2014, LI-COR Biosciences
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
! \brief       Calculate extra sensible heat fluxes (Burba et al. 2008, GCB)
!              for LI-7500(A)
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine BurbaTerms()
    use m_rp_global_var
    implicit none
    !> local variables
    real(kind = dbl), parameter :: r_top = 0.0225d0 ![m]
    real(kind = dbl), parameter :: r_spar = 0.0025d0 ![m]
    real(kind = dbl), parameter :: l_bot = 0.065d0 ![m]
    real(kind = dbl), parameter :: l_top = 0.045d0 ![m]
    real(kind = dbl), parameter :: l_spar = 0.005d0 ![m]
    real(kind = dbl) :: deT_bot
    real(kind = dbl) :: deT_top
    real(kind = dbl) :: deT_spar
    real(kind = dbl) :: delta_bot
    real(kind = dbl) :: delta_top
    real(kind = dbl) :: delta_spar
    real(kind = dbl) :: k_air
    real(kind = dbl) :: Umean


    !> Is user explicitly selects not to perform Burba correction, set terms to zero.
    if (RPsetup%bu_corr == 'none') then
        Burba%h_bot  = 0d0
        Burba%h_top  = 0d0
        Burba%h_spar = 0d0
        return
    end if

    !> Air Conductivity [W m-1 K-1]
    if (Stats%T /= error) then
        k_air = 0.000067d0 * (Stats%T - 273.15d0) + 0.024343d0
    else
        k_air = error
    end if

    !> Mean horizonal wind
    if (Stats%Mean(u) /= error .and. Stats%Mean(v) /= error) then
        Umean = dsqrt(Stats%Mean(u)**2 + Stats%Mean(v)**2)
    else
        Umean = error
    end if

    !> If any parameter for multi linear is not available, set method to simple linear
    if (BiometVar%Rg == error .or. BiometVar%LWin == error) then
        RPsetup%bu_multi = .false.
    end if

    if (RPsetup%bu_multi) then
        !> multiple regression option
        if (Stats%daytime) then
            if (Stats%T /= error .and. BiometVar%Rg /= error .and. Umean /= error) then
                deT_bot  = BurbaPar%m(daytime, bot, 1)  &
                              + BurbaPar%m(daytime, bot, 2)  * (Stats%T - 273.16d0)  &
                              + BurbaPar%m(daytime, bot, 3)  * BiometVar%Rg &
                              + BurbaPar%m(daytime, bot, 4)  * Umean
                deT_top  = BurbaPar%m(daytime, top, 1)  &
                              + BurbaPar%m(daytime, top, 2)  * (Stats%T - 273.16d0) &
                              + BurbaPar%m(daytime, top, 3)  * BiometVar%Rg &
                              + BurbaPar%m(daytime, top, 4)  * Umean
                deT_spar = BurbaPar%m(daytime, spar, 1) &
                              + BurbaPar%m(daytime, spar, 2) * (Stats%T - 273.16d0) &
                              + BurbaPar%m(daytime, spar, 3) * BiometVar%Rg &
                              + BurbaPar%m(daytime, spar, 4) * Umean
            else
                deT_bot  = error
                deT_top  = error
                deT_spar  = error
            end if
        else
            if (Stats%T /= error .and. BiometVar%LWin /= error .and. Umean /= error) then
                deT_bot  = BurbaPar%m(nighttime, bot, 1)  &
                              + BurbaPar%m(nighttime, bot, 2)  * (Stats%T - 273.16d0) &
                              + BurbaPar%m(nighttime, bot, 3)  * BiometVar%LWin &
                              + BurbaPar%m(nighttime, bot, 4)  * Umean
                deT_top  = BurbaPar%m(nighttime, top, 1)  &
                              + BurbaPar%m(nighttime, top, 2)  * (Stats%T - 273.16d0)  &
                              + BurbaPar%m(nighttime, top, 3)  * BiometVar%LWin &
                              + BurbaPar%m(nighttime, top, 4)  * Umean
                deT_spar = BurbaPar%m(nighttime, spar, 1) &
                              + BurbaPar%m(nighttime, spar, 2) * (Stats%T - 273.16d0)  &
                              + BurbaPar%m(nighttime, spar, 3) * BiometVar%LWin &
                              + BurbaPar%m(nighttime, spar, 4) * Umean
            else
                deT_bot  = error
                deT_top  = error
                deT_spar  = error
            end if
        end if
    else
        !> Simple linear regression option
        if(Stats%daytime) then
            if (Stats%T /= error) then
                deT_bot  = BurbaPar%l(daytime, bot, 1) * (Stats%T - 273.16d0)&
                              + BurbaPar%l(daytime, bot, 2) + 273.16d0 - Stats%T
                deT_top  = BurbaPar%l(daytime, top, 1) * (Stats%T - 273.16d0) &
                              + BurbaPar%l(daytime, top, 2) + 273.16d0 - Stats%T
                deT_spar = BurbaPar%l(daytime, spar, 1) * (Stats%T - 273.16d0) &
                              + BurbaPar%l(daytime, spar, 2) + 273.16d0 - Stats%T
            else
                deT_bot  = error
                deT_top  = error
                deT_spar  = error
            end if
        else
            if (Stats%T /= error) then
                deT_bot  = BurbaPar%l(nighttime, bot, 1) * (Stats%T - 273.16d0) &
                              + BurbaPar%l(nighttime, bot, 2) + 273.16d0 - Stats%T
                deT_top  = BurbaPar%l(nighttime, top, 1) * (Stats%T - 273.16d0) &
                              + BurbaPar%l(nighttime, top, 2) + 273.16d0 - Stats%T
                deT_spar = BurbaPar%l(nighttime, spar, 1) * (Stats%T - 273.16d0) &
                              + BurbaPar%l(nighttime, spar, 2) + 273.16d0 - Stats%T
            else
                deT_bot  = error
                deT_top  = error
                deT_spar  = error
            end if
        end if
    end if

    !> U-depedent parameters
    if (Umean /= 0d0 .and. Umean /= error) then
        delta_bot  = 0.004d0  * dsqrt(l_bot / Umean) + 0.004d0
        delta_top  = 0.0028d0 * dsqrt(l_top / Umean) + 0.00025d0 / Umean + 0.0045d0
        delta_spar = 0.0058d0 * dsqrt(l_spar / Umean)
    else
        delta_bot  = error
        delta_top  = error
        delta_spar = error
    endif

    !> Sensible heat fluxes from sensor surfaces (Burba et al. 2008, GCB; Nobel (1983))
    !> Bottom
    if (delta_bot /= 0d0 .and. delta_bot /= error &
        .and.  k_air /= error .and. deT_bot /= error) then
        Burba%h_bot  = k_air * deT_bot / delta_bot
    else
        Burba%h_bot  = 0d0
    end if

    !> Top
    if (delta_top /= 0d0 .and. delta_top /= error &
        .and.  k_air /= error .and. deT_top /= error) then
        Burba%h_top  = k_air * deT_top * (r_top + delta_top) &
                          / (delta_top * r_top)
    else
        Burba%h_top  = 0d0
    end if

    !> Spars
    if (delta_spar /= error .and. k_air /= error &
        .and. deT_spar /= error) then
        Burba%h_spar = 0.15d0 * k_air * deT_spar &
                          / (r_spar * dlog((r_spar + delta_spar) / r_spar))
    else
        Burba%h_spar  = 0d0
    end if
end subroutine BurbaTerms
