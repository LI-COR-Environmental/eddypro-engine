!***************************************************************************
! fluxes23_rp.f90
! ---------------
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
! \brief       Completes flux corrections (Level 2/3). Applies WPL and \n
!              spectral corrections in the appropriate order
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo        Merge with corresponding FCC sub into a common one
!***************************************************************************
subroutine Fluxes23_rp()
    use m_rp_global_var
    implicit none
    !> local variables
    real(kind = dbl) :: Tp
    real(kind = dbl) :: E_nowpl

    write(*,'(a)', advance = 'no') '  Calculating fluxes Level 2 and 3..'

    Flux2 = fluxtype('', '', error, error, error, error, error, error, error, error, error, &
        error, error, error, error, error, error, error, error)
    Flux3 = fluxtype('', '', error, error, error, error, error, error, error, error, error, &
        error, error, error, error, error, error, error, error)

    !> Level 2 end 3 internal sensible heat, do nothing
    Flux2%Hi_co2 = Flux1%Hi_co2
    Flux2%Hi_h2o = Flux1%Hi_h2o
    Flux2%Hi_ch4 = Flux1%Hi_ch4
    Flux2%Hi_gas4 = Flux1%Hi_gas4
    Flux3%Hi_co2 = Flux2%Hi_co2
    Flux3%Hi_h2o = Flux2%Hi_h2o
    Flux3%Hi_ch4 = Flux2%Hi_ch4
    Flux3%Hi_gas4 = Flux2%Hi_gas4

    !> Level 2 evapotranspiration, WPL corrected including Burba if the case
    if (E2Col(h2o)%Instr%path_type == 'open') then
        if (LitePar%RhoCp > 0d0 .and. LitePar%Ta > 0d0 &
            .and. Flux1%E /= error .and. Flux1%H /= error .and. LitePar%sigma /= error) then
            if(index(E2Col(h2o)%Instr%model, 'li7500') /= 0) then
                Flux2%E = (1d0 + mu * LitePar%sigma) * Flux1%E &
                    + (1d0 + mu * LitePar%sigma) * (Flux1%H + Burba%h_top + Burba%h_bot + Burba%h_spar)&
                    * RHO%w / (LitePar%RhoCp * LitePar%Ta)
            else
                Flux2%E = (1d0 + mu * LitePar%sigma) * Flux1%E &
                    + (1d0 + mu * LitePar%sigma) * Flux1%H &
                    * RHO%w / (LitePar%RhoCp * LitePar%Ta)
            end if
        else
            Flux2%E = error
        end if
    else
        select case(E2Col(h2o)%measure_type)
            case ('molar_density', 'mole_fraction')
                if (Flux1%E /= error .and. LitePar%sigma /= error .and. E2Col(h2o)%Va > 0d0 .and. LitePar%Va > 0d0) then
                    if (Flux1%Hi_h2o /= error .and.  Stats%cov(w, pi) /= error) then
                        Flux2%E = (1d0 + mu * LitePar%sigma) * Flux1%E * E2Col(h2o)%Va / LitePar%Va &
                            + (1d0 + mu * LitePar%sigma) * Flux1%Hi_h2o &
                            * RHO%w / (LitePar%RhoCp * LitePar%Tcell) &
                            - (1d0 + mu * LitePar%sigma) * Stats%cov(w, pi) &
                            * RHO%w / (LitePar%Pcell)
                    elseif (Flux1%Hi_h2o /= error) then
                       Flux2%E = (1d0 + mu * LitePar%sigma) * Flux1%E * E2Col(h2o)%Va / LitePar%Va &
                            + (1d0 + mu * LitePar%sigma) * Flux1%Hi_h2o &
                            * RHO%w / (LitePar%RhoCp * LitePar%Tcell)
                    elseif (Stats%cov(w, pi)  /= error) then
                        Flux2%E = (1d0 + mu * LitePar%sigma) * Flux1%E * E2Col(h2o)%Va / LitePar%Va &
                            - (1d0 + mu * LitePar%sigma) * Stats%cov(w, pi) &
                            * RHO%w / (LitePar%Pcell)
                    else
                        Flux2%E = Flux1%E * E2Col(h2o)%Va / LitePar%Va
                    end if
                else
                    Flux2%E = error
                end if
            case ('mixing_ratio')
                Flux2%E = Flux1%E
        end select
    end if

    !> If WPL should not be applied
    if (.not. EddyProProj%wpl) Flux2%E = Flux1%E

    !> But if Flux1 was error, then set also 2 to error
    if (Flux1%E == error) Flux2%E = error

    !> Level 2 h2o and latent heat flux
    if (Flux2%E /= error) then
        Flux2%h2o = Flux2%E * 1d3 / MW(h2o)
        Flux2%LE = Flux2%E * LitePar%lambda
    else
        Flux2%h2o = error
        Flux2%LE  = error
    end if

    !> Level 2 evapotranspiration fluxes with H2O covariances at timelags of other scalars
    !> Do nothing, WPL is deleterious here
    Flux2%E_co2 = Flux1%E_co2
    Flux2%E_ch4 = Flux1%E_ch4
    Flux2%E_gas4 = Flux1%E_gas4

    !> Level 2 Sensible heat
    if (E2Col(ts)%instr%category == 'sonic') then
        !> corrected for humidity, after Van Dyjk et al. (2004) eq. 3.53
        !> revising Schotanus et al. (1983)
        if (Flux1%H /= error) then
            if(Flux0%E /= error .and. Stats%Cov(w, ts) /= error .and. RHO%a > 0d0 .and. &
                LitePar%Q > 0d0 .and. LitePar%RhoCp > 0d0 .and. LitePar%alpha /= error) then
                Flux2%H = Flux1%H &
                    - LitePar%RhoCp * LitePar%alpha * Stats%Mean(ts) * Flux0%E / RHO%a &
                    - LitePar%RhoCp * LitePar%alpha * LitePar%Q * Stats%Cov(w, ts)
                    !> alternative
                    !- LitePar%RhoCp * LitePar%alpha * LitePar%Ta * Flux0%E / RHO%a
            else
                Flux2%H = Flux1%H
            end if
        else
            Flux2%H = error
        end if
    else
        !> Equal Level 1 if Ts is not coming from a sonic, rather from
        !> a fast temperature sensor (such as a thermocouple)
        Flux2%H = Flux1%H
    end if

    !> Map for temperature factor (Van Dijk et al. 2004, eq.3.1)
    !> Currently not applied (to investigate better)
    !if(LitePar%Tmap /= error .and. Flux2%H /= error) then
        !Flux2%H  = Flux2%H * LitePar%Tmap
    !end if

    !> Level 3 sensible heat, spectral corrected
    if(Flux2%H /= error .and. BPCF%of(w_ts) /= error) then
        Flux3%H = Flux2%H * BPCF%of(w_ts)
    else
        Flux3%H = error
    end if

    !> Level 3 for evapotranspiration: for open path, WPL again with corrected H
    !> Starts again from Level 1 of E, Level 2 was only used to calculate H Level 3.
    if(E2Col(h2o)%Instr%path_type == 'open') then
        if (LitePar%RhoCp > 0d0 .and. LitePar%Ta > 0d0 &
            .and. Flux1%E /= error .and. Flux1%H /= error .and. LitePar%sigma /= error) then
            if(index(E2Col(h2o)%Instr%model, 'li7500') /= 0) then
                Flux3%E = (1d0 + mu * LitePar%sigma) * Flux1%E &
                    + (1d0 + mu * LitePar%sigma) * (Flux3%H + Burba%h_top + Burba%h_bot + Burba%h_spar)&
                    * RHO%w / (LitePar%RhoCp * LitePar%Ta)
            else
                Flux3%E = (1d0 + mu * LitePar%sigma) * Flux1%E &
                    + (1d0 + mu * LitePar%sigma) * Flux3%H &
                    * RHO%w / (LitePar%RhoCp * LitePar%Ta)
            end if
        else
            Flux3%E = Flux2%E
        end if
    else
        Flux3%E = Flux2%E
    end if

    !> If WPL should not be applied
    if (.not. EddyProProj%wpl) Flux3%E = Flux2%E

    !> Level 3 latent heat fluxes with H2O covariances at timelags of other scalars
    !> Do nothing
    Flux3%E_co2 = Flux2%E_co2
    Flux3%E_ch4 = Flux2%E_ch4
    Flux3%E_gas4 = Flux2%E_gas4

    !> Level 3 h2o flux and latent heat flux
    if (Flux3%E /= error) then
        Flux3%h2o = Flux3%E * 1d3 / MW(h2o)
        Flux3%LE = Flux3%E * LitePar%lambda
    else
        Flux3%h2o = error
        Flux3%LE  = error
    end if

    !> Calculate E_nowpl for closed and open path systems
    if (E2Col(h2o)%Instr%path_type == 'closed') then
        if (Flux1%E /= error .and. BPCF%of(w_h2o) /= error) then
            E_nowpl = Flux1%E * BPCF%of(w_h2o)
        elseif(Flux1%E /= error) then
            E_nowpl = Flux1%E
        else
            E_nowpl = error
        end if
    else
        E_nowpl = Flux1%E
    end if

    !> Apply spectral correction to h2o and LE Level 3 fluxes for closed path
    if (E2Col(h2o)%Instr%path_type == 'closed' .and. Flux3%E /= error) then
        !> Level 3, spectral correction
        Flux3%h2o = Flux3%h2o * BPCF%of(w_h2o)  !< closed path water flux  was not FR corrected yet
        Flux3%E   = Flux3%E   * BPCF%of(w_h2o)  !< closed path evapotrans. was not FR corrected yet
        Flux3%LE  = Flux3%LE  * BPCF%of(w_h2o)  !< closed path latent heat was not FR corrected yet
    end if

    if (.not. E2Col(h2o)%present) then
        Flux3%h2o = error
        Flux3%E   = error
        Flux3%LE  = error
    end if

    !> Level 2 other gases
    !> co2
    if (E2Col(co2)%Instr%path_type == 'closed') then
        !> Level 2, WPL for closed path, implemented after Ibrom et al. (2007, Tellus, eq. 3a with H contribution from WPL24)
        select case(E2Col(co2)%measure_type)  !< Analitically, it is verified that (E * mu * sigma / rho%w * co2_density) equals (r_c * w'chi_w' / Va)
            case('molar_density')

                !> E, T and P effects
                if (Flux1%co2 /= error .and. LitePar%Va > 0d0 .and. Stats%chi(co2) > 0d0 &
                    .and. Flux3%E_co2 /= error .and. RHO%w > 0d0 &
                    .and. Flux3%Hi_co2 /= error .and. LitePar%RhoCp > 0d0 .and. LitePar%Tcell > 0d0 &
                    .and. Stats%cov(w, pi) /= error .and. LitePar%Pcell /= error) then
                    Flux2%co2 = Flux1%co2 * E2Col(co2)%Va / LitePar%Va &
                        + Flux3%E_co2 * mu * LitePar%sigma / RHO%w &
                        * Stats%chi(co2) / LitePar%Va &
                        + (1d0 + mu * LitePar%sigma) * Flux3%Hi_co2 / (LitePar%RhoCp * LitePar%Tcell) &
                        * Stats%chi(co2) / LitePar%Va &
                        - (1d0 + mu * LitePar%sigma) * Stats%cov(w, pi) / (LitePar%Pcell) &
                        * Stats%chi(co2) / LitePar%Va

                !> Onlt E and T effects
                elseif (Flux1%co2 /= error .and. LitePar%Va > 0d0 .and. Stats%chi(co2) > 0d0 &
                    .and. Flux3%E_co2 /= error .and. RHO%w > 0d0 &
                    .and. Flux3%Hi_co2 /= error .and. LitePar%RhoCp > 0d0 .and. LitePar%Tcell > 0d0) then
                    Flux2%co2 = Flux1%co2 * E2Col(co2)%Va / LitePar%Va &
                        + Flux3%E_co2 * mu * LitePar%sigma / RHO%w &
                        * Stats%chi(co2) / LitePar%Va &
                        + (1d0 + mu * LitePar%sigma) * Flux3%Hi_co2 / (LitePar%RhoCp * LitePar%Tcell) &
                        * Stats%chi(co2) / LitePar%Va

                !> Onlt T and P effects
                elseif (Flux1%co2 /= error .and. LitePar%Va > 0d0 .and. Stats%chi(co2) > 0d0 &
                    .and. Flux3%Hi_co2 /= error .and. LitePar%RhoCp > 0d0 .and. LitePar%Tcell > 0d0 &
                    .and. Stats%cov(w, pi) /= error .and. LitePar%Pcell /= error) then
                    Flux2%co2 = Flux1%co2 * E2Col(co2)%Va / LitePar%Va &
                        + (1d0 + mu * LitePar%sigma) * Flux3%Hi_co2 / (LitePar%RhoCp * LitePar%Tcell) &
                        * Stats%chi(co2) / LitePar%Va &
                        - (1d0 + mu * LitePar%sigma) * Stats%cov(w, pi) / (LitePar%Pcell) &
                        * Stats%chi(co2) / LitePar%Va

                !> Only E and P effects
                elseif (Flux1%co2 /= error .and. LitePar%Va > 0d0 .and. Stats%chi(co2) > 0d0 &
                    .and. Flux3%E_co2 /= error .and. RHO%w > 0d0 &
                    .and. Stats%cov(w, pi) /= error .and. LitePar%Pcell /= error) then
                    Flux2%co2 = Flux1%co2 * E2Col(co2)%Va / LitePar%Va &
                        + Flux3%E_co2 * mu * LitePar%sigma / RHO%w &
                        * Stats%chi(co2) / LitePar%Va &
                        - (1d0 + mu * LitePar%sigma) * Stats%cov(w, pi) / (LitePar%Pcell) &
                        * Stats%chi(co2) / LitePar%Va

                !> Onlt E effects
                elseif (Flux1%co2 /= error .and. LitePar%Va > 0d0 .and. Stats%chi(co2) > 0d0 &
                    .and. Flux3%E_co2 /= error .and. RHO%w > 0d0 ) then
                    Flux2%co2 = Flux1%co2 * E2Col(co2)%Va / LitePar%Va &
                        + Flux3%E_co2 * mu * LitePar%sigma / RHO%w &
                        * Stats%chi(co2) / LitePar%Va

                !> Onlt T effects
                elseif (Flux1%co2 /= error .and. LitePar%Va > 0d0 .and. Stats%chi(co2) > 0d0 &
                    .and. Flux3%Hi_co2 /= error .and. LitePar%RhoCp > 0d0 .and. LitePar%Tcell > 0d0) then
                    Flux2%co2 = Flux1%co2 * E2Col(co2)%Va / LitePar%Va &
                        + (1d0 + mu * LitePar%sigma) * Flux3%Hi_co2 / (LitePar%RhoCp * LitePar%Tcell) &
                        * Stats%chi(co2) / LitePar%Va

                !> Onlt P effects
                elseif (Flux1%co2 /= error .and. LitePar%Va > 0d0 .and. Stats%chi(co2) > 0d0 &
                    .and. Stats%cov(w, pi) /= error .and. LitePar%Pcell /= error) then
                    Flux2%co2 = Flux1%co2 * E2Col(co2)%Va / LitePar%Va &
                        - (1d0 + mu * LitePar%sigma) * Stats%cov(w, pi) / (LitePar%Pcell) &
                        * Stats%chi(co2) / LitePar%Va

                !> No WPL effects
                elseif (Flux1%co2 /= error) then
                    Flux2%co2 = Flux1%co2 * E2Col(co2)%Va / LitePar%Va
                else
                    Flux2%co2 = error
                end if

            case('mole_fraction')   !< Analitically, it is verified that (E * mu * sigma / rho%w * co2_density) equals (r_c * w'chi_w' / Va)

                !> E, T and P effects
                if (Flux1%co2 /= error .and. LitePar%Va > 0d0 .and. Stats%chi(co2) > 0d0 &
                    .and. Flux3%E_co2 /= error .and. RHO%w > 0d0 &
                    .and. Flux3%Hi_co2 /= error .and. LitePar%RhoCp > 0d0 .and. LitePar%Tcell > 0d0 &
                    .and. Stats%cov(w, pi) /= error .and. LitePar%Pcell /= error) then
                    Flux2%co2 = Flux1%co2 &
                        + Flux3%E_co2 * mu * LitePar%sigma / RHO%w &
                        * Stats%chi(co2) / LitePar%Va &
                        + (1d0 + mu * LitePar%sigma) * Flux3%Hi_co2 / (LitePar%RhoCp * LitePar%Tcell) &
                        * Stats%chi(co2) / LitePar%Va &
                        - (1d0 + mu * LitePar%sigma) * Stats%cov(w, pi) / (LitePar%Pcell) &
                        * Stats%chi(co2) / LitePar%Va

                !> Only E and T effects
                elseif (Flux1%co2 /= error .and. LitePar%Va > 0d0 .and. Stats%chi(co2) > 0d0 &
                    .and. Flux3%E_co2 /= error .and. RHO%w > 0d0 &
                    .and. Flux3%Hi_co2 /= error .and. LitePar%RhoCp > 0d0 .and. LitePar%Tcell > 0d0) then
                    Flux2%co2 = Flux1%co2 &
                        + Flux3%E_co2 * mu * LitePar%sigma / RHO%w &
                        * Stats%chi(co2) / LitePar%Va &
                        + (1d0 + mu * LitePar%sigma) * Flux3%Hi_co2 / (LitePar%RhoCp * LitePar%Tcell) &
                        * Stats%chi(co2) / LitePar%Va

                !> Only T and P effects
                elseif (Flux1%co2 /= error .and. LitePar%Va > 0d0 .and. Stats%chi(co2) > 0d0 &
                    .and. Flux3%Hi_co2 /= error .and. LitePar%RhoCp > 0d0 .and. LitePar%Tcell > 0d0 &
                    .and. Stats%cov(w, pi) /= error .and. LitePar%Pcell /= error) then
                    Flux2%co2 = Flux1%co2 &
                        + (1d0 + mu * LitePar%sigma) * Flux3%Hi_co2 / (LitePar%RhoCp * LitePar%Tcell) &
                        * Stats%chi(co2) / LitePar%Va &
                        - (1d0 + mu * LitePar%sigma) * Stats%cov(w, pi) / (LitePar%Pcell) &
                        * Stats%chi(co2) / LitePar%Va

                !> Only E and P effects
                elseif (Flux1%co2 /= error .and. LitePar%Va > 0d0 .and. Stats%chi(co2) > 0d0 &
                    .and. Flux3%E_co2 /= error .and. RHO%w > 0d0 &
                    .and. Stats%cov(w, pi) /= error .and. LitePar%Pcell /= error) then
                    Flux2%co2 = Flux1%co2 &
                        + Flux3%E_co2 * mu * LitePar%sigma / RHO%w &
                        * Stats%chi(co2) / LitePar%Va &
                        - (1d0 + mu * LitePar%sigma) * Stats%cov(w, pi) / (LitePar%Pcell) &
                        * Stats%chi(co2) / LitePar%Va

                !> Only E effects
                elseif (Flux1%co2 /= error .and. LitePar%Va > 0d0 .and. Stats%chi(co2) > 0d0 &
                    .and. Flux3%E_co2 /= error .and. RHO%w > 0d0 ) then
                    Flux2%co2 = Flux1%co2 &
                        + Flux3%E_co2 * mu * LitePar%sigma / RHO%w &
                        * Stats%chi(co2) / LitePar%Va

                !> Onlt T effects
                elseif (Flux1%co2 /= error .and. LitePar%Va > 0d0 .and. Stats%chi(co2) > 0d0 &
                    .and. Flux3%Hi_co2 /= error .and. LitePar%RhoCp > 0d0 .and. LitePar%Tcell > 0d0) then
                    Flux2%co2 = Flux1%co2 &
                        + (1d0 + mu * LitePar%sigma) * Flux3%Hi_co2 / (LitePar%RhoCp * LitePar%Tcell) &
                        * Stats%chi(co2) / LitePar%Va

                !> Onlt P effects
                elseif (Flux1%co2 /= error .and. LitePar%Va > 0d0 .and. Stats%chi(co2) > 0d0 &
                    .and. Stats%cov(w, pi) /= error .and. LitePar%Pcell /= error) then
                    Flux2%co2 = Flux1%co2 &
                        - (1d0 + mu * LitePar%sigma) * Stats%cov(w, pi) / (LitePar%Pcell) &
                        * Stats%chi(co2) / LitePar%Va

                !> No WPL effects
                elseif (Flux1%co2 /= error) then
                    Flux2%co2 = Flux1%co2
                else
                    Flux2%co2 = error
                end if
            case('mixing_ratio')
                Flux2%co2 = Flux1%co2
        end select
    else
        !> Level 2, WPL for open path implemented after e.g. Burba et al. (2008, GCB, eq. 1)
        if(index(E2Col(co2)%Instr%model, 'li7500') /= 0) then
            if (Flux1%co2 /= error .and. Flux3%H /= error .and. Flux3%E /= error &
                .and. RHO%d > 0 .and. LitePar%RhoCp > 0 .and. LitePar%Ta > 0d0 .and. LitePar%sigma /= error) then
                Flux2%co2 = Flux1%co2 + mu * Flux3%E * Stats%d(co2) * 1d3 &
                    / ((1d0 + mu * LitePar%sigma) * RHO%d) &
                    + (Flux3%H + Burba%h_top + Burba%h_bot + Burba%h_spar) * Stats%d(co2) * 1d3 &
                    / (LitePar%RhoCp * LitePar%Ta)
            elseif(Flux1%co2 /= error .and. Flux3%H /= error &
                .and. LitePar%RhoCp >0 .and. LitePar%Ta > 0d0) then
                Flux2%co2 = Flux1%co2 &
                    + (Flux3%H + Burba%h_top + Burba%h_bot + Burba%h_spar) * Stats%d(co2) * 1d3 &
                    / (LitePar%RhoCp * LitePar%Ta)
            elseif(Flux1%co2 /= error .and. Flux3%E /= error &
                .and. RHO%d > 0d0 .and. LitePar%sigma /= error) then
                Flux2%co2 = Flux1%co2 &
                    + mu * Flux3%E * Stats%d(co2) * 1d3 / ((1d0 + mu * LitePar%sigma) * RHO%d)
            elseif(Flux1%co2 /= error) then
                Flux2%co2 = Flux1%co2
            else
                Flux2%co2 = error
            end if
        else
            if (Flux1%co2 /= error .and. Flux3%H /= error .and. Flux3%E /= error &
                .and. RHO%d > 0d0 .and. LitePar%RhoCp > 0d0 .and. LitePar%Ta > 0d0 .and. LitePar%sigma /= error) then
                Flux2%co2 = Flux1%co2 &
                    + mu * Flux3%E * Stats%d(co2) * 1d3 / ((1d0 + mu * LitePar%sigma) * RHO%d) &
                    + Flux3%H * Stats%d(co2) * 1d3 / (LitePar%RhoCp * LitePar%Ta)
            elseif(Flux1%co2 /= error .and. Flux3%H /= error &
                .and. LitePar%RhoCp > 0d0 .and. LitePar%Ta > 0d0) then
                Flux2%co2 = Flux1%co2 &
                    + Flux3%H * Stats%d(co2) * 1d3 / (LitePar%RhoCp * LitePar%Ta)
            elseif(Flux1%co2 /= error .and. Flux3%E /= error &
                .and. RHO%d > 0d0 .and. LitePar%sigma /= error) then
                Flux2%co2 = Flux1%co2 &
                    + mu * Flux3%E * Stats%d(co2) * 1d3 / ((1d0 + mu * LitePar%sigma) * RHO%d)
            elseif(Flux1%co2 /= error) then
                Flux2%co2 = Flux1%co2
            else
                Flux2%co2 = error
            end if
        end if
    end if
    if (.not. E2Col(co2)%present) Flux2%co2 = error

    !> ch4
    if (E2Col(ch4)%Instr%path_type == 'closed') then
        !> Level 2, WPL for closed path, implemented after Ibrom et al. (2007, Tellus, eq. 3a with H contribution from WPL24)
        select case(E2Col(ch4)%measure_type)
            case('molar_density')

                !> E, T and P effects
                if (Flux1%ch4 /= error .and. LitePar%Va > 0d0 .and. Stats%chi(ch4) > 0d0 &
                    .and. Flux3%E_ch4 /= error .and. RHO%w > 0d0 &
                    .and. Flux3%Hi_ch4 /= error .and. LitePar%RhoCp > 0d0 .and. LitePar%Tcell > 0d0 &
                    .and. Stats%cov(w, pi) /= error .and. LitePar%Pcell /= error) then
                    Flux2%ch4 = Flux1%ch4 * E2Col(ch4)%Va / LitePar%Va &
                        + Flux3%E_ch4 * mu * LitePar%sigma / RHO%w &
                        * Stats%chi(ch4) / LitePar%Va &
                        + (1d0 + mu * LitePar%sigma) * Flux3%Hi_ch4 / (LitePar%RhoCp * LitePar%Tcell) &
                        * Stats%chi(ch4) / LitePar%Va &
                        - (1d0 + mu * LitePar%sigma) * Stats%cov(w, pi) / (LitePar%Pcell) &
                        * Stats%chi(ch4) / LitePar%Va

                !> Onlt E and T effects
                elseif (Flux1%ch4 /= error .and. LitePar%Va > 0d0 .and. Stats%chi(ch4) > 0d0 &
                    .and. Flux3%E_ch4 /= error .and. RHO%w > 0d0 &
                    .and. Flux3%Hi_ch4 /= error .and. LitePar%RhoCp > 0d0 .and. LitePar%Tcell > 0d0) then
                    Flux2%ch4 = Flux1%ch4 * E2Col(ch4)%Va / LitePar%Va &
                        + Flux3%E_ch4 * mu * LitePar%sigma / RHO%w &
                        * Stats%chi(ch4) / LitePar%Va &
                        + (1d0 + mu * LitePar%sigma) * Flux3%Hi_ch4 / (LitePar%RhoCp * LitePar%Tcell) &
                        * Stats%chi(ch4) / LitePar%Va

                !> Onlt T and P effects
                elseif (Flux1%ch4 /= error .and. LitePar%Va > 0d0 .and. Stats%chi(ch4) > 0d0 &
                    .and. Flux3%Hi_ch4 /= error .and. LitePar%RhoCp > 0d0 .and. LitePar%Tcell > 0d0 &
                    .and. Stats%cov(w, pi) /= error .and. LitePar%Pcell /= error) then
                    Flux2%ch4 = Flux1%ch4 * E2Col(ch4)%Va / LitePar%Va &
                        + (1d0 + mu * LitePar%sigma) * Flux3%Hi_ch4 / (LitePar%RhoCp * LitePar%Tcell) &
                        * Stats%chi(ch4) / LitePar%Va &
                        - (1d0 + mu * LitePar%sigma) * Stats%cov(w, pi) / (LitePar%Pcell) &
                        * Stats%chi(ch4) / LitePar%Va

                !> Only E and P effects
                elseif (Flux1%ch4 /= error .and. LitePar%Va > 0d0 .and. Stats%chi(ch4) > 0d0 &
                    .and. Flux3%E_ch4 /= error .and. RHO%w > 0d0 &
                    .and. Stats%cov(w, pi) /= error .and. LitePar%Pcell /= error) then
                    Flux2%ch4 = Flux1%ch4 &
                        + Flux3%E_ch4 * mu * LitePar%sigma / RHO%w &
                        * Stats%chi(ch4) / LitePar%Va &
                        - (1d0 + mu * LitePar%sigma) * Stats%cov(w, pi) / (LitePar%Pcell) &
                        * Stats%chi(ch4) / LitePar%Va

                !> Onlt E effects
                elseif (Flux1%ch4 /= error .and. LitePar%Va > 0d0 .and. Stats%chi(ch4) > 0d0 &
                    .and. Flux3%E_ch4 /= error .and. RHO%w > 0d0 ) then
                    Flux2%ch4 = Flux1%ch4 * E2Col(ch4)%Va / LitePar%Va &
                        + Flux3%E_ch4 * mu * LitePar%sigma / RHO%w &
                        * Stats%chi(ch4) / LitePar%Va

                !> Onlt T effects
                elseif (Flux1%ch4 /= error .and. LitePar%Va > 0d0 .and. Stats%chi(ch4) > 0d0 &
                    .and. Flux3%Hi_ch4 /= error .and. LitePar%RhoCp > 0d0 .and. LitePar%Tcell > 0d0) then
                    Flux2%ch4 = Flux1%ch4 * E2Col(ch4)%Va / LitePar%Va &
                        + (1d0 + mu * LitePar%sigma) * Flux3%Hi_ch4 / (LitePar%RhoCp * LitePar%Tcell) &
                        * Stats%chi(ch4) / LitePar%Va

                !> Onlt P effects
                elseif (Flux1%ch4 /= error .and. LitePar%Va > 0d0 .and. Stats%chi(ch4) > 0d0 &
                    .and. Stats%cov(w, pi) /= error .and. LitePar%Pcell /= error) then
                    Flux2%ch4 = Flux1%ch4 * E2Col(ch4)%Va / LitePar%Va &
                        - (1d0 + mu * LitePar%sigma) * Stats%cov(w, pi) / (LitePar%Pcell) &
                        * Stats%chi(ch4) / LitePar%Va

                !> No WPL effects
                elseif (Flux1%ch4 /= error) then
                    Flux2%ch4 = Flux1%ch4 * E2Col(ch4)%Va / LitePar%Va
                else
                    Flux2%ch4 = error
                end if

            case('mole_fraction')

                !> E, T and P effects
                if (Flux1%ch4 /= error .and. LitePar%Va > 0d0 .and. Stats%chi(ch4) > 0d0 &
                    .and. Flux3%E_ch4 /= error .and. RHO%w > 0d0 &
                    .and. Flux3%Hi_ch4 /= error .and. LitePar%RhoCp > 0d0 .and. LitePar%Tcell > 0d0 &
                    .and. Stats%cov(w, pi) /= error .and. LitePar%Pcell /= error) then
                    Flux2%ch4 = Flux1%ch4 &
                        + Flux3%E_ch4 * mu * LitePar%sigma / RHO%w &
                        * Stats%chi(ch4) / LitePar%Va &
                        + (1d0 + mu * LitePar%sigma) * Flux3%Hi_ch4 / (LitePar%RhoCp * LitePar%Tcell) &
                        * Stats%chi(ch4) / LitePar%Va &
                        - (1d0 + mu * LitePar%sigma) * Stats%cov(w, pi) / (LitePar%Pcell) &
                        * Stats%chi(ch4) / LitePar%Va

                !> Onlt E and T effects
                elseif (Flux1%ch4 /= error .and. LitePar%Va > 0d0 .and. Stats%chi(ch4) > 0d0 &
                    .and. Flux3%E_ch4 /= error .and. RHO%w > 0d0 &
                    .and. Flux3%Hi_ch4 /= error .and. LitePar%RhoCp > 0d0 .and. LitePar%Tcell > 0d0) then
                    Flux2%ch4 = Flux1%ch4 &
                        + Flux3%E_ch4 * mu * LitePar%sigma / RHO%w &
                        * Stats%chi(ch4) / LitePar%Va &
                        + (1d0 + mu * LitePar%sigma) * Flux3%Hi_ch4 / (LitePar%RhoCp * LitePar%Tcell) &
                        * Stats%chi(ch4) / LitePar%Va

                !> Onlt T and P effects
                elseif (Flux1%ch4 /= error .and. LitePar%Va > 0d0 .and. Stats%chi(ch4) > 0d0 &
                    .and. Flux3%Hi_ch4 /= error .and. LitePar%RhoCp > 0d0 .and. LitePar%Tcell > 0d0 &
                    .and. Stats%cov(w, pi) /= error .and. LitePar%Pcell /= error) then
                    Flux2%ch4 = Flux1%ch4 &
                        + (1d0 + mu * LitePar%sigma) * Flux3%Hi_ch4 / (LitePar%RhoCp * LitePar%Tcell) &
                        * Stats%chi(ch4) / LitePar%Va &
                        - (1d0 + mu * LitePar%sigma) * Stats%cov(w, pi) / (LitePar%Pcell) &
                        * Stats%chi(ch4) / LitePar%Va

                !> Only E and P effects
                elseif (Flux1%ch4 /= error .and. LitePar%Va > 0d0 .and. Stats%chi(ch4) > 0d0 &
                    .and. Flux3%E_ch4 /= error .and. RHO%w > 0d0 &
                    .and. Stats%cov(w, pi) /= error .and. LitePar%Pcell /= error) then
                    Flux2%ch4 = Flux1%ch4 &
                        + Flux3%E_ch4 * mu * LitePar%sigma / RHO%w &
                        * Stats%chi(ch4) / LitePar%Va &
                        - (1d0 + mu * LitePar%sigma) * Stats%cov(w, pi) / (LitePar%Pcell) &
                        * Stats%chi(ch4) / LitePar%Va

                !> Onlt E effects
                elseif (Flux1%ch4 /= error .and. LitePar%Va > 0d0 .and. Stats%chi(ch4) > 0d0 &
                    .and. Flux3%E_ch4 /= error .and. RHO%w > 0d0 ) then
                    Flux2%ch4 = Flux1%ch4 &
                        + Flux3%E_ch4 * mu * LitePar%sigma / RHO%w &
                        * Stats%chi(ch4) / LitePar%Va

                !> Onlt T effects
                elseif (Flux1%ch4 /= error .and. LitePar%Va > 0d0 .and. Stats%chi(ch4) > 0d0 &
                    .and. Flux3%Hi_ch4 /= error .and. LitePar%RhoCp > 0d0 .and. LitePar%Tcell > 0d0) then
                    Flux2%ch4 = Flux1%ch4 &
                        + (1d0 + mu * LitePar%sigma) * Flux3%Hi_ch4 / (LitePar%RhoCp * LitePar%Tcell) &
                        * Stats%chi(ch4) / LitePar%Va

                !> Onlt P effects
                elseif (Flux1%ch4 /= error .and. LitePar%Va > 0d0 .and. Stats%chi(ch4) > 0d0 &
                    .and. Stats%cov(w, pi) /= error .and. LitePar%Pcell /= error) then
                    Flux2%ch4 = Flux1%ch4 &
                        - (1d0 + mu * LitePar%sigma) * Stats%cov(w, pi) / (LitePar%Pcell) &
                        * Stats%chi(ch4) / LitePar%Va

                !> No WPL effects
                elseif (Flux1%ch4 /= error) then
                    Flux2%ch4 = Flux1%ch4
                else
                    Flux2%ch4 = error
                end if
            case('mixing_ratio')
                Flux2%ch4 = Flux1%ch4
        end select
    else
        !> Level 2, WPL for open path implemented after Webb et al. 1980 (+ 7700 multipliers)
        if (E2Col(ch4)%Instr%model(1:len_trim(E2Col(ch4)%Instr%model) - 2)  == 'li7700') then
            !> Flux formulation including multipliers requires the use of the original WPL formulation
            !> making use of E not corrected for WPL (E_nowpl) and H, both spectrally corrected.
            if (Flux1%ch4 /= error .and. Flux3%H /= error .and. E_nowpl /= error &
                .and. RHO%d > 0d0 .and. LitePar%RhoCp > 0d0 .and. LitePar%Ta > 0d0 .and. LitePar%sigma /= error) then
                Flux2%ch4 = Mul7700%A *(Flux1%ch4 &
                          + Mul7700%B * mu * Stats%d(ch4) * 1d3 * E_nowpl / RHO%d &
                          + Mul7700%C * (1d0 + mu * LitePar%sigma) * Flux3%H * Stats%d(ch4) * 1d3 / (LitePar%RhoCp * LitePar%Ta))
            elseif(Flux1%ch4 /= error .and. Flux3%H /= error &
                .and. LitePar%RhoCp > 0d0 .and. LitePar%Ta > 0d0 .and. LitePar%sigma /= error) then
                Flux2%ch4 = Mul7700%A *(Flux1%ch4 &
                          + Mul7700%C * (1d0 + mu * LitePar%sigma) * Flux3%H * Stats%d(ch4) * 1d3 / (LitePar%RhoCp * LitePar%Ta))
            elseif(Flux1%ch4 /= error .and. E_nowpl /= error .and. RHO%d > 0d0) then
                Flux2%ch4 = Mul7700%A *(Flux1%ch4 &
                          + Mul7700%B * mu * Stats%d(ch4) * 1d3 * E_nowpl / RHO%d)
            elseif (Flux1%ch4 /= error) then
                Flux2%ch4 = Mul7700%A * Flux1%ch4
            else
                Flux2%ch4 = error
            end if
        else
            !> Flux formulation without multipliers follows Burba et al. 2008.
            if (Flux1%ch4 /= error .and. Flux3%H /= error .and. Flux3%E /= error &
                .and. RHO%d > 0d0 .and. LitePar%RhoCp > 0d0 .and. LitePar%Ta > 0d0 .and. LitePar%sigma /= error) then
                Flux2%ch4 = Flux1%ch4 &
                    + mu * Flux3%E * Stats%d(ch4) * 1d3 / ((1d0 + mu * LitePar%sigma) * RHO%d) &
                    + Flux3%H * Stats%d(ch4) * 1d3 / (LitePar%RhoCp * LitePar%Ta)
            elseif(Flux1%ch4 /= error .and. Flux3%H /= error &
                .and. LitePar%RhoCp > 0d0 .and. LitePar%Ta > 0d0) then
                Flux2%ch4 = Flux1%ch4 &
                    + Flux3%H * Stats%d(ch4) * 1d3 / (LitePar%RhoCp * LitePar%Ta)
            elseif(Flux1%ch4 /= error .and. Flux3%E /= error &
                .and. RHO%d > 0d0 .and. LitePar%sigma /= error) then
                Flux2%ch4 = Flux1%ch4 &
                    + mu * Flux3%E * Stats%d(ch4) * 1d3 / ((1d0 + mu * LitePar%sigma) * RHO%d)
            elseif(Flux1%ch4 /= error) then
                Flux2%ch4 = Flux1%ch4
            end if
        end if
    end if
    if (.not. E2Col(ch4)%present) Flux2%ch4 = error

    !> gas4
    if (E2Col(gas4)%Instr%path_type == 'closed') then


        !> Level 2, WPL for closed path, implemented after Ibrom et al. (2007, Tellus, eq. 3a with H contribution from WPL24)
        select case(E2Col(gas4)%measure_type)
            case('molar_density')

                !> E, T and P effects
                if (Flux1%gas4 /= error .and. LitePar%Va > 0d0 .and. Stats%chi(gas4) > 0d0 &
                    .and. Flux3%E_gas4 /= error .and. RHO%w > 0d0 &
                    .and. Flux3%Hi_gas4 /= error .and. LitePar%RhoCp > 0d0 .and. LitePar%Tcell > 0d0 &
                    .and. Stats%cov(w, pi) /= error .and. LitePar%Pcell > 0d0) then
                    Flux2%gas4 = Flux1%gas4 * E2Col(gas4)%Va / LitePar%Va &
                        + Flux3%E_gas4 * mu * LitePar%sigma / RHO%w &
                        * Stats%chi(gas4) / LitePar%Va &
                        + (1d0 + mu * LitePar%sigma) * Flux3%Hi_gas4 / (LitePar%RhoCp * LitePar%Tcell) &
                        * Stats%chi(gas4) / LitePar%Va &
                        - (1d0 + mu * LitePar%sigma) * Stats%cov(w, pi) / (LitePar%Pcell) &
                        * Stats%chi(gas4) / LitePar%Va

                !> Onlt E and T effects
                elseif (Flux1%gas4 /= error .and. LitePar%Va > 0d0 .and. Stats%chi(gas4) > 0d0 &
                    .and. Flux3%E_gas4 /= error .and. RHO%w > 0d0 &
                    .and. Flux3%Hi_gas4 /= error .and. LitePar%RhoCp > 0d0 .and. LitePar%Tcell > 0d0) then
                    Flux2%gas4 = Flux1%gas4 * E2Col(gas4)%Va / LitePar%Va &
                        + Flux3%E_gas4 * mu * LitePar%sigma / RHO%w &
                        * Stats%chi(gas4) / LitePar%Va &
                        + (1d0 + mu * LitePar%sigma) * Flux3%Hi_gas4 / (LitePar%RhoCp * LitePar%Tcell) &
                        * Stats%chi(gas4) / LitePar%Va

                !> Onlt T and P effects
                elseif (Flux1%gas4 /= error .and. LitePar%Va > 0d0 .and. Stats%chi(gas4) > 0d0 &
                    .and. Flux3%Hi_gas4 /= error .and. LitePar%RhoCp > 0d0 .and. LitePar%Tcell > 0d0 &
                    .and. Stats%cov(w, pi) /= error .and. LitePar%Pcell > 0d0) then
                    Flux2%gas4 = Flux1%gas4 * E2Col(gas4)%Va / LitePar%Va &
                        + (1d0 + mu * LitePar%sigma) * Flux3%Hi_gas4 / (LitePar%RhoCp * LitePar%Tcell) &
                        * Stats%chi(gas4) / LitePar%Va &
                        - (1d0 + mu * LitePar%sigma) * Stats%cov(w, pi) / (LitePar%Pcell) &
                        * Stats%chi(gas4) / LitePar%Va

                !> Only E and P effects
                elseif (Flux1%gas4 /= error .and. LitePar%Va > 0d0 .and. Stats%chi(gas4) > 0d0 &
                    .and. Flux3%E_gas4 /= error .and. RHO%w > 0d0 &
                    .and. Stats%cov(w, pi) /= error .and. LitePar%Pcell > 0d0) then
                    Flux2%gas4 = Flux1%gas4 * E2Col(gas4)%Va / LitePar%Va &
                        + Flux3%E_gas4 * mu * LitePar%sigma / RHO%w &
                        * Stats%chi(gas4) / LitePar%Va &
                        - (1d0 + mu * LitePar%sigma) * Stats%cov(w, pi) / (LitePar%Pcell) &
                        * Stats%chi(gas4) / LitePar%Va

                !> Onlt E effects
                elseif (Flux1%gas4 /= error .and. LitePar%Va > 0d0 .and. Stats%chi(gas4) > 0d0 &
                    .and. Flux3%E_gas4 /= error .and. RHO%w > 0d0 ) then
                    Flux2%gas4 = Flux1%gas4 * E2Col(gas4)%Va / LitePar%Va &
                        + Flux3%E_gas4 * mu * LitePar%sigma / RHO%w &
                        * Stats%chi(gas4) / LitePar%Va

                !> Onlt T effects
                elseif (Flux1%gas4 /= error .and. LitePar%Va > 0d0 .and. Stats%chi(gas4) > 0d0 &
                    .and. Flux3%Hi_gas4 /= error .and. LitePar%RhoCp > 0d0 .and. LitePar%Tcell > 0d0) then
                    Flux2%gas4 = Flux1%gas4 * E2Col(gas4)%Va / LitePar%Va &
                        + (1d0 + mu * LitePar%sigma) * Flux3%Hi_gas4 / (LitePar%RhoCp * LitePar%Tcell) &
                        * Stats%chi(gas4) / LitePar%Va

                !> Onlt P effects
                elseif (Flux1%gas4 /= error .and. LitePar%Va > 0d0 .and. Stats%chi(gas4) > 0d0 &
                    .and. Stats%cov(w, pi) /= error .and. LitePar%Pcell > 0d0) then
                    Flux2%gas4 = Flux1%gas4 * E2Col(gas4)%Va / LitePar%Va &
                        - (1d0 + mu * LitePar%sigma) * Stats%cov(w, pi) / (LitePar%Pcell) &
                        * Stats%chi(gas4) / LitePar%Va

                !> No WPL effects
                elseif (Flux1%gas4 /= error) then
                    Flux2%gas4 = Flux1%gas4 * E2Col(gas4)%Va / LitePar%Va
                else
                    Flux2%gas4 = error
                end if

            case('mole_fraction')

                !> E, T and P effects
                if (Flux1%gas4 /= error .and. LitePar%Va > 0d0 .and. Stats%chi(gas4) > 0d0 &
                    .and. Flux3%E_gas4 /= error .and. RHO%w > 0d0 &
                    .and. Flux3%Hi_gas4 /= error .and. LitePar%RhoCp > 0d0 .and. LitePar%Tcell > 0d0 &
                    .and. Stats%cov(w, pi) /= error .and. LitePar%Pcell > 0d0) then
                    Flux2%gas4 = Flux1%gas4 &
                        + Flux3%E_gas4 * mu * LitePar%sigma / RHO%w &
                        * Stats%chi(gas4) / LitePar%Va &
                        + (1d0 + mu * LitePar%sigma) * Flux3%Hi_gas4 / (LitePar%RhoCp * LitePar%Tcell) &
                        * Stats%chi(gas4) / LitePar%Va &
                        - (1d0 + mu * LitePar%sigma) * Stats%cov(w, pi) / (LitePar%Pcell) &
                        * Stats%chi(gas4) / LitePar%Va

                !> Onlt E and T effects
                elseif (Flux1%gas4 /= error .and. LitePar%Va > 0d0 .and. Stats%chi(gas4) > 0d0 &
                    .and. Flux3%E_gas4 /= error .and. RHO%w > 0d0 &
                    .and. Flux3%Hi_gas4 /= error .and. LitePar%RhoCp > 0d0 .and. LitePar%Tcell > 0d0) then
                    Flux2%gas4 = Flux1%gas4 &
                        + Flux3%E_gas4 * mu * LitePar%sigma / RHO%w &
                        * Stats%chi(gas4) / LitePar%Va &
                        + (1d0 + mu * LitePar%sigma) * Flux3%Hi_gas4 / (LitePar%RhoCp * LitePar%Tcell) &
                        * Stats%chi(gas4) / LitePar%Va

                !> Onlt T and P effects
                elseif (Flux1%gas4 /= error .and. LitePar%Va > 0d0 .and. Stats%chi(gas4) > 0d0 &
                    .and. Flux3%Hi_gas4 /= error .and. LitePar%RhoCp > 0d0 .and. LitePar%Tcell > 0d0 &
                    .and. Stats%cov(w, pi) /= error .and. LitePar%Pcell > 0d0) then
                    Flux2%gas4 = Flux1%gas4 &
                        + (1d0 + mu * LitePar%sigma) * Flux3%Hi_gas4 / (LitePar%RhoCp * LitePar%Tcell) &
                        * Stats%chi(gas4) / LitePar%Va &
                        - (1d0 + mu * LitePar%sigma) * Stats%cov(w, pi) / (LitePar%Pcell) &
                        * Stats%chi(gas4) / LitePar%Va

                !> Only E and P effects
                elseif (Flux1%gas4 /= error .and. LitePar%Va > 0d0 .and. Stats%chi(gas4) > 0d0 &
                    .and. Flux3%E_gas4 /= error .and. RHO%w > 0d0 &
                    .and. Stats%cov(w, pi) /= error .and. LitePar%Pcell > 0d0) then
                    Flux2%gas4 = Flux1%gas4 &
                        + Flux3%E_gas4 * mu * LitePar%sigma / RHO%w &
                        * Stats%chi(gas4) / LitePar%Va &
                        - (1d0 + mu * LitePar%sigma) * Stats%cov(w, pi) / (LitePar%Pcell) &
                        * Stats%chi(gas4) / LitePar%Va

                !> Onlt E effects
                elseif (Flux1%gas4 /= error .and. LitePar%Va > 0d0 .and. Stats%chi(gas4) > 0d0 &
                    .and. Flux3%E_gas4 /= error .and. RHO%w > 0d0 ) then
                    Flux2%gas4 = Flux1%gas4 &
                        + Flux3%E_gas4 * mu * LitePar%sigma / RHO%w &
                        * Stats%chi(gas4) / LitePar%Va

                !> Onlt T effects
                elseif (Flux1%gas4 /= error .and. LitePar%Va > 0d0 .and. Stats%chi(gas4) > 0d0 &
                    .and. Flux3%Hi_gas4 /= error .and. LitePar%RhoCp > 0d0 .and. LitePar%Tcell > 0d0) then
                    Flux2%gas4 = Flux1%gas4 &
                        + (1d0 + mu * LitePar%sigma) * Flux3%Hi_gas4 / (LitePar%RhoCp * LitePar%Tcell) &
                        * Stats%chi(gas4) / LitePar%Va

                !> Onlt P effects
                elseif (Flux1%gas4 /= error .and. LitePar%Va > 0d0 .and. Stats%chi(gas4) > 0d0 &
                    .and. Stats%cov(w, pi) /= error .and. LitePar%Pcell > 0d0) then
                    Flux2%gas4 = Flux1%gas4 &
                        - (1d0 + mu * LitePar%sigma) * Stats%cov(w, pi) / (LitePar%Pcell) &
                        * Stats%chi(gas4) / LitePar%Va

                !> No WPL effects
                elseif (Flux1%gas4 /= error) then
                    Flux2%gas4 = Flux1%gas4
                else
                    Flux2%gas4 = error
                end if
            case('mixing_ratio')
                Flux2%gas4 = Flux1%gas4
        end select
    else
        !> Level 2, WPL for open path implemented after e.g. Burba et al. (2008, GCB, eq. 1)
        if (Flux1%gas4 /= error .and. Flux3%H /= error .and. Flux3%E /= error &
            .and. RHO%d * LitePar%RhoCp * LitePar%Ta /= 0d0 .and. LitePar%sigma /= error) then
            Flux2%gas4 = Flux1%gas4 &
                + mu * Flux3%E * Stats%d(gas4) * 1d3 / ((1d0 + mu * LitePar%sigma) * RHO%d) &
                + Flux3%H * Stats%d(gas4) * 1d3 / (LitePar%RhoCp * LitePar%Ta)
        elseif(Flux1%gas4 /= error .and. Flux3%H /= error &
            .and. LitePar%RhoCp > 0d0 .and. LitePar%Ta > 0d0) then
            Flux2%gas4 = Flux1%gas4 &
                + Flux3%H * Stats%d(gas4) * 1d3 / (LitePar%RhoCp * LitePar%Ta)
        elseif(Flux1%gas4 /= error .and. Flux3%E /= error &
            .and. RHO%d > 0d0 .and. LitePar%sigma /= error) then
            Flux2%gas4 = Flux1%gas4 &
                + mu * Flux3%E * Stats%d(gas4) * 1d3 / ((1d0 + mu * LitePar%sigma) * RHO%d)
        elseif(Flux1%gas4 /= error) then
            Flux2%gas4 = Flux1%gas4
        else
            Flux2%gas4 = error
        end if
    end if
    if (.not. E2Col(gas4)%present) Flux2%gas4 = error

    !> If WPL should not be applied
    if (.not. EddyProProj%wpl) then
        Flux2%co2 = Flux1%co2
        Flux2%ch4 = Flux1%ch4
        Flux2%gas4 = Flux1%gas4
    end if

    !> Level 3 other gases. For closed path apply now the spectral correction (e.g. Ibrom et al. 2007)
    !> co2
    if (E2Col(co2)%Instr%path_type == 'closed' .and. Flux2%co2 /= error) then
        !> Level 3, spectral correction
        Flux3%co2 = Flux2%co2 * BPCF%of(w_co2)
    else
        !> Level 3, spectral correction was already applied
        Flux3%co2 = Flux2%co2
    end if
    !> ch4
    if (E2Col(ch4)%Instr%path_type == 'closed' .and. Flux2%ch4 /= error) then
        !> Level 3, spectral correction
        Flux3%ch4 = Flux2%ch4 * BPCF%of(w_ch4)
    else
        !> Level 3, spectral correction was already applied
        Flux3%ch4 = Flux2%ch4
    end if
    !> gas4
    if (E2Col(gas4)%Instr%path_type == 'closed' .and. Flux2%gas4 /= error) then
        !> Level 3, spectral correction
        Flux3%gas4 = Flux2%gas4 * BPCF%of(w_n2o)
    else
        !> Level 3, spectral correction was already applied
        Flux3%gas4 = Flux2%gas4
    end if

    !> Potential temperature
    !> If condition fails, previous value (from Fluxes0) holds for z/LitePar%L
    if (Stats%Pr > 0d0) then
        Tp = LitePar%Ta * (1d5 / Stats%Pr)**(.286d0)
    else
        Tp = error
    end if

    !> Monin-Obukhov length (L = - (Tp^ /(k*g))*(ustar**3/(w'Tp')^ in m)
    if (Flux3%H /= 0d0 .and. Flux3%H /= error .and. &
        LitePar%RhoCp > 0d0 .and. LitePar%us >= 0d0 .and. Tp > error) then
        LitePar%L = -Tp * (LitePar%us**3) / (vk * g * Flux3%H / LitePar%RhoCp)
    else
        LitePar%L = error
    end if

    !> Monin-Obukhov stability parameter (zL = z/L)
    !> If condition fails, previous value (from Fluxes0) holds
    if (LitePar%L /= 0d0 .and. LitePar%L /= error) &
        LitePar%zL = (E2Col(u)%Instr%height - Metadata%d) / LitePar%L

    !> scale temperature(T*)
    !> If condition fails, previous value (from Fluxes0) holds
    if (LitePar%us > 0d0 .and. Flux3%H /= error .and. LitePar%RhoCp > 0d0) &
        LitePar%Ts = Flux3%H / (LitePar%RhoCp * LitePar%us)

    !> Bowen ration (Bowen, 1926, Phyis Rev)
    if (Flux3%LE /= 0d0 .and. Flux3%LE /= error) &
        LitePar%Bowen = Flux3%H / Flux3%LE

    !> Momentum flux
    Flux2%tau = Flux1%tau
    Flux3%tau = Flux1%tau

    !> If fluxes are error, set also time lags to error, just for clarity
    if (Flux2%co2  == error) Essentials%timelag(co2)  = error
    if (Flux2%h2o  == error) Essentials%timelag(h2o)  = error
    if (Flux2%ch4  == error) Essentials%timelag(ch4)  = error
    if (Flux2%gas4 == error) Essentials%timelag(gas4) = error

    write(*,'(a)')   ' done.'
end subroutine Fluxes23_rp
