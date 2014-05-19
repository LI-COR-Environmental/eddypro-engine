!***************************************************************************
! fluxes0_rp.f90
! --------------
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
! \brief       Calculates uncorrected fluxes (refer to pseudo-code)
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine Fluxes0_rp(printout)
    use m_rp_global_var
    implicit none
    !> in/out variables
    logical, intent(in) :: printout
    !> local variables
    real(kind = dbl) :: Tp

    if (printout) write(*,'(a)', advance = 'no') '  Calculating fluxes Level 0..'

    Flux0 = fluxtype('', '', error, error, error, error, error, error, error, error, error, &
        error, error, error, error, error, error, error, error)

    !> Sensible heat flux, H in [W m-2], Cp in [J Kg-1K-1]
    if (LitePar%RhoCp > 0d0 .and. Stats%Cov(w, ts) /= error) then
        Flux0%H = LitePar%RhoCp * Stats%Cov(w, ts)
    else
        Flux0%H = error
    end if

    !> Random error on sensible heat flux, Cp in [J Kg-1K-1]
    if (RUsetup%meth /= 'none') then
        if (LitePar%RhoCp > 0d0 .and. Essentials%rand_uncer(ts) /= error) &
            Essentials%rand_uncer(ts) = Essentials%rand_uncer(ts) * LitePar%RhoCp
    end if

    !> Internal sensible heat flux, Hint in [W m-2], Cp in [J Kg-1K-1]
    if (LitePar%RhoCp > 0d0) then
        if (Stats%tc_cov_tl_co2 /= error) then
            Flux0%Hi_co2 = LitePar%RhoCp * Stats%tc_cov_tl_co2
        else
            Flux0%Hi_co2 = error
        end if

        if (Stats%tc_cov_tl_h2o /= error) then
            Flux0%Hi_h2o = LitePar%RhoCp * Stats%tc_cov_tl_h2o
        else
            Flux0%Hi_h2o = error
        end if

        if (Stats%tc_cov_tl_ch4 /= error) then
            Flux0%Hi_ch4 = LitePar%RhoCp * Stats%tc_cov_tl_ch4
        else
            Flux0%Hi_ch4 = error
        end if

        if (Stats%tc_cov_tl_gas4 /= error) then
            Flux0%Hi_gas4 = LitePar%RhoCp * Stats%tc_cov_tl_gas4
        else
            Flux0%Hi_gas4 = error
        end if
    else
        Flux0%Hi_co2 = error
        Flux0%Hi_h2o = error
        Flux0%Hi_ch4 = error
        Flux0%Hi_gas4 = error
    end if

    !> Uncorrected co2 flux [umol m-2 s-1]
    if (E2Col(co2)%present) then
        if(E2Col(co2)%measure_type == 'molar_density') then
            if (Stats%Cov(w, co2) /= error) then
                Flux0%co2 = Stats%Cov(w, co2) * 1d3
            else
                Flux0%co2 = error
            end if
            !> Random uncertainty
            if (RUsetup%meth /= 'none') then
                if (Essentials%rand_uncer(co2) /= error) &
                    Essentials%rand_uncer(co2) = Essentials%rand_uncer(co2) * 1d3
            end if
        else if(E2Col(co2)%measure_type == 'mixing_ratio') then
            if (LitePar%Vd > 0d0 .and. Stats%Cov(w, co2) /= error) then
                Flux0%co2 = Stats%Cov(w, co2) / LitePar%Vd
            else
                Flux0%co2 = error
            end if
            !> Random uncertainty
            if (RUsetup%meth /= 'none') then
                if (Essentials%rand_uncer(co2) /= error .and. LitePar%Vd > 0d0) then
                    Essentials%rand_uncer(co2) = Essentials%rand_uncer(co2) / LitePar%Vd
                else
                    Essentials%rand_uncer(co2) = error
                end if
            end if
        else if(E2Col(co2)%measure_type == 'mole_fraction') then
            if (LitePar%Va /= 0d0 .and. LitePar%Va /= error .and. Stats%Cov(w, co2) /= error) then
                Flux0%co2 = Stats%Cov(w, co2) / LitePar%Va
            else
                Flux0%co2 = error
            end if
            !> Random uncertainty
            if (RUsetup%meth /= 'none') then
                if (Essentials%rand_uncer(co2) /= error .and. LitePar%Va > 0d0) then
                    Essentials%rand_uncer(co2) = Essentials%rand_uncer(co2) / LitePar%Va
                else
                    Essentials%rand_uncer(co2) = error
                end if
            end if
        end if
    else
        Flux0%co2 = error
    end if

    !> Uncorrected h2o flux [mmol m-2 s-1]
    if (E2Col(h2o)%present) then
        if(E2Col(h2o)%measure_type == 'molar_density') then
            if (Stats%Cov(w, h2o) /= error) then
                Flux0%h2o = Stats%Cov(w, h2o)
            else
                Flux0%h2o = error
            end if
        else if(E2Col(h2o)%measure_type == 'mole_fraction') then
            if (LitePar%Va > 0d0 .and. Stats%Cov(w, h2o) /= error) then
                Flux0%h2o = Stats%Cov(w, h2o) / LitePar%Va
            else
                Flux0%h2o = error
            end if
            !> Random uncertainty
            if (RUsetup%meth /= 'none') then
                if (Essentials%rand_uncer(h2o) /= error .and. LitePar%Va > 0d0) then
                    Essentials%rand_uncer(h2o) = Essentials%rand_uncer(h2o) / LitePar%Va
                else
                    Essentials%rand_uncer(h2o) = error
                end if
            end if
        else if(E2Col(h2o)%measure_type == 'mixing_ratio') then
            if (LitePar%Vd > 0d0 .and. Stats%Cov(w, h2o) /= error) then
                Flux0%h2o = Stats%Cov(w, h2o) / LitePar%Vd
            else
                Flux0%h2o = error
            end if
            !> Random uncertainty
            if (RUsetup%meth /= 'none') then
                if (Essentials%rand_uncer(h2o) /= error .and. LitePar%Vd > 0d0) then
                    Essentials%rand_uncer(h2o) = Essentials%rand_uncer(h2o) / LitePar%Vd
                else
                    Essentials%rand_uncer(h2o) = error
                end if
            end if
        end if
    else
        Flux0%h2o = error
    end if

    !> Uncorrected ch4 flux [umol m-2 s-1]
    if (E2Col(ch4)%present) then
        if(E2Col(ch4)%measure_type == 'molar_density') then
            if (Stats%Cov(w, ch4) /= error) then
                Flux0%ch4 = Stats%Cov(w, ch4) * 1d3
            else
                Flux0%ch4 = error
            end if
            !> Random uncertainty
            if (RUsetup%meth /= 'none') then
                if (Essentials%rand_uncer(ch4) /= error) &
                    Essentials%rand_uncer(ch4) = Essentials%rand_uncer(ch4) * 1d3
            end if
        else if(E2Col(ch4)%measure_type == 'mixing_ratio') then
            if (LitePar%Vd > 0d0 .and. Stats%Cov(w, ch4) /= error) then
                Flux0%ch4 = Stats%Cov(w, ch4) / LitePar%Vd
            else
                Flux0%ch4 = error
            end if
            !> Random uncertainty
            if (RUsetup%meth /= 'none') then
                if (Essentials%rand_uncer(ch4) /= error .and. LitePar%Vd > 0d0) then
                    Essentials%rand_uncer(ch4) = Essentials%rand_uncer(ch4) / LitePar%Vd
                else
                    Essentials%rand_uncer(ch4) = error
                end if
            end if
        else if(E2Col(ch4)%measure_type == 'mole_fraction') then
            if (LitePar%Va > 0d0 .and. Stats%Cov(w, ch4) /= error) then
                Flux0%ch4 = Stats%Cov(w, ch4) / LitePar%Va
            else
                Flux0%ch4 = error
            end if
            !> Random uncertainty
            if (RUsetup%meth /= 'none') then
                if (Essentials%rand_uncer(ch4) /= error .and. LitePar%Va > 0d0) then
                    Essentials%rand_uncer(ch4) = Essentials%rand_uncer(ch4) / LitePar%Va
                else
                    Essentials%rand_uncer(ch4) = error
                end if
            end if
        end if
    else
        Flux0%ch4 = error
    end if

    !> Uncorrected gas4 flux [umol m-2 s-1]
    if (E2Col(gas4)%present) then
        if(E2Col(gas4)%measure_type == 'molar_density') then
            if (Stats%Cov(w, gas4) /= error) then
                Flux0%gas4 = Stats%Cov(w, gas4) * 1d3
            else
                Flux0%gas4 = error
            end if
            !> Random uncertainty
            if (RUsetup%meth /= 'none') then
                if (Essentials%rand_uncer(gas4) /= error) &
                    Essentials%rand_uncer(gas4) = Essentials%rand_uncer(gas4) * 1d3
            end if
        else if(E2Col(gas4)%measure_type == 'mixing_ratio') then
            if (LitePar%Vd > 0d0 .and. Stats%Cov(w, gas4) /= error) then
                Flux0%gas4 = Stats%Cov(w, gas4) / LitePar%Vd
            else
                Flux0%gas4 = error
            end if
            !> Random uncertainty
            if (RUsetup%meth /= 'none') then
                if (Essentials%rand_uncer(gas4) /= error .and. LitePar%Vd > 0d0) then
                    Essentials%rand_uncer(gas4) = Essentials%rand_uncer(gas4) / LitePar%Vd
                else
                    Essentials%rand_uncer(gas4) = error
                end if
            end if
        else if(E2Col(gas4)%measure_type == 'mole_fraction') then
            if (LitePar%Va > 0d0 .and. Stats%Cov(w, gas4) /= error) then
                Flux0%gas4 = Stats%Cov(w, gas4) / LitePar%Va
            else
                Flux0%gas4 = error
            end if
            !> Random uncertainty
            if (RUsetup%meth /= 'none') then
                if (Essentials%rand_uncer(gas4) /= error .and. LitePar%Va > 0d0) then
                    Essentials%rand_uncer(gas4) = Essentials%rand_uncer(gas4) / LitePar%Va
                else
                    Essentials%rand_uncer(gas4) = error
                end if
            end if
        end if
    else
        Flux0%gas4 = error
    end if

    !> Latent heat flux [W m-2], lambda in [J+1kg-1]
    if (Flux0%h2o /= error .and. LitePar%lambda > 0d0) then
        Flux0%LE = Flux0%h2o * LitePar%lambda * MW(h2o) * 1d-3
    else
        Flux0%LE = error
    end if

    !> Random uncertainty on Latent heat flux, lambda in [J+1kg-1]
    if (RUsetup%meth /= 'none') then
        if (Essentials%rand_uncer(h2o) /= error .and. LitePar%lambda > 0d0) then
            Essentials%rand_uncer_LE = Essentials%rand_uncer(h2o) * LitePar%lambda * MW(h2o) * 1d-3
        else
            Essentials%rand_uncer_LE = error
        end if
    end if

    !> Level 0 evapotranspiration flux [kg m-2 -1]
    if (Flux0%h2o /= error .and. LitePar%lambda > 0d0) then
        Flux0%E = Flux0%h2o * MW(h2o) * 1d-3
    else
        Flux0%E = error
    end if

    !> Level 0 evapotranspiration flux [kg m-2 -1] with H2O covariances at timelags of other scalars
    if (E2Col(h2o)%Instr%path_type == 'closed') then
        if(E2Col(h2o)%measure_type == 'molar_density') then
            if(Stats%h2ocov_tl_co2 /= error) then
                Flux0%E_co2 = Stats%h2ocov_tl_co2 * MW(h2o) * 1d-3
            else
                Flux0%E_co2 = error
            end if
            if(Stats%h2ocov_tl_ch4 /= error) then
                Flux0%E_ch4 = Stats%h2ocov_tl_ch4 * MW(h2o) * 1d-3
            else
                Flux0%E_ch4 = error
            end if
            if(Stats%h2ocov_tl_gas4 /= error) then
                Flux0%E_gas4 = Stats%h2ocov_tl_gas4 * MW(h2o) * 1d-3
            else
                Flux0%E_gas4 = error
            end if
        else if(E2Col(h2o)%measure_type == 'mole_fraction') then
            if (LitePar%Va > 0d0 .and. Stats%h2ocov_tl_co2 /= error) then
                Flux0%E_co2 = Stats%h2ocov_tl_co2  * MW(h2o) * 1d-3 / LitePar%Va
            else
                Flux0%E_co2 = error
            end if
            if (LitePar%Va > 0d0 .and. Stats%h2ocov_tl_ch4 /= error) then
                Flux0%E_ch4 = Stats%h2ocov_tl_ch4  * MW(h2o) * 1d-3 / LitePar%Va
            else
                Flux0%E_ch4 = error
            end if
            if (LitePar%Va > 0d0 .and. Stats%h2ocov_tl_gas4 /= error) then
                Flux0%E_gas4 = Stats%h2ocov_tl_gas4  * MW(h2o) * 1d-3 / LitePar%Va
            else
                Flux0%E_gas4 = error
            end if
        else if (E2Col(h2o)%measure_type == 'mixing_ratio') then
            if (LitePar%Vd > 0d0 .and. Stats%h2ocov_tl_co2 /= error) then
                Flux0%E_co2 = Stats%h2ocov_tl_co2  * MW(h2o) * 1d-3 / LitePar%Vd
            else
                Flux0%E_co2 = error
            end if
            if (LitePar%Vd > 0d0 .and. Stats%h2ocov_tl_ch4 /= error) then
                Flux0%E_ch4 = Stats%h2ocov_tl_ch4  * MW(h2o) * 1d-3 / LitePar%Vd
            else
                Flux0%E_ch4 = error
            end if
            if (LitePar%Vd > 0d0 .and. Stats%h2ocov_tl_gas4 /= error) then
                Flux0%E_gas4 = Stats%h2ocov_tl_gas4  * MW(h2o) * 1d-3 / LitePar%Vd
            else
                Flux0%E_gas4 = error
            end if
        end if
    else
        Flux0%E_co2 = error
        Flux0%E_ch4 = error
        Flux0%E_gas4 = error
    end if

    !> Friction velocity [m s-1]
    if (Stats%Cov(u, w) /= error .and. Stats%Cov(v, w) /= error) then
        LitePar%us = (Stats%Cov(u, w)**2 + Stats%Cov(v, w)**2)**(0.25d0)
    else
        LitePar%us = error
    end if

    !> Momentum flux [kg m-1 s-2], after Van Dijk et al. 2004 Eq. 2.44
    if (RHO%a > 0d0 .and. LitePar%us >= 0d0) then
        Flux0%tau = RHO%a * LitePar%us ** 2d0
    else
        Flux0%tau = error
    end if

    !> Random error on momentum flux
    if (RUsetup%meth /= 'none') then
        if (Essentials%rand_uncer(u) /= error &
            .and. RHO%a > 0d0 .and. LitePar%us >= 0d0) then
            Essentials%rand_uncer(u) = Essentials%rand_uncer(u) * RHO%a
        else
            Essentials%rand_uncer(u) = error
        end if
    end if

    !> Potential temperature
    if (Stats%Pr > 0d0 .and. Stats%T > 0d0) then
        Tp = Stats%T * (1d5 / Stats%Pr)**(.286d0)
    else
        Tp = error
    end if

    !> Monin-Obukhov length (L = - (Tp^ /(k*g))*(ustar**3/(w'Tp')^ in m)
    if (Stats%Cov(w, ts) /= 0d0 .and. Stats%Cov(w, ts) /= error .and. &
        LitePar%us > 0d0 .and. Tp > 0d0) then
        LitePar%L = -Tp * (LitePar%us**3) / (vk * g * Stats%Cov(w, ts))
    else
        LitePar%L = error
    end if
    Essentials%L = LitePar%L

    !> Monin-Obukhov stability parameter (zL = z/L)
    if (LitePar%L /= 0d0 .and. LitePar%L /= error) then
        LitePar%zL = (E2Col(u)%Instr%height - Metadata%d) / LitePar%L
    else
        LitePar%zL = error
    end if
    Essentials%zL = LitePar%zL

    !> SL dynamic temperature(T*), see e.g. Foken and Wichura (1996)
    if (LitePar%us > 0d0 .and. Stats%Cov(w, ts) /= error) then
        LitePar%Ts = - Stats%Cov(w, ts) / LitePar%us
    else
        LitePar%Ts = error
    end if

    !> Bowen ration (Bowen, 1926, Phyis Rev)
    if (Flux0%LE /= 0d0 .and. Flux0%LE /= error) then
        LitePar%Bowen = Flux0%H / Flux0%LE
    else
        LitePar%Bowen = error
    end if
    if (printout) write(*,'(a)')   ' done.'
end subroutine Fluxes0_rp
