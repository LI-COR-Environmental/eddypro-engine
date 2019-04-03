!***************************************************************************
! fluxes23.f90
! ------------
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
! \brief       Completes flux correction (Level 2/3). Applies WPL and \n
!              spectral corrections in the appropriate order
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine Fluxes23(lEx)
    use m_fx_global_var
    implicit none
    !> In/out variables
    type(ExType), intent(inout) :: lEx
    !> local variables
    real(kind = dbl) :: Tp
    real(kind = dbl) :: E_nowpl
    real(kind = dbl), parameter :: alpha = 0.51d0


    Flux2 = errFlux
    Flux3 = errFlux

    !> Level 2 end 3 internal sensible heat, do nothing
    Flux2%Hi_co2 = Flux1%Hi_co2
    Flux2%Hi_h2o = Flux1%Hi_h2o
    Flux2%Hi_ch4 = Flux1%Hi_ch4
    Flux2%Hi_gas4 = Flux1%Hi_gas4
    Flux3%Hi_co2 = Flux2%Hi_co2
    Flux3%Hi_h2o = Flux2%Hi_h2o
    Flux3%Hi_ch4 = Flux2%Hi_ch4
    Flux3%Hi_gas4 = Flux2%Hi_gas4

    !> Level 2 evapotranspiration WPL corrected, including Burba if the case
    if (EddyProProj%wpl) then
        if (lEx%instr(ih2o)%path_type == 'open') then
            if (lEx%RhoCp > 0d0 .and. lEx%Ta > 0d0 &
                .and. Flux1%E /= error .and. Flux1%H /= error &
                .and. lEx%sigma /= error) then
                    !> Open-path uses Webb et al. (1980)
                    !> Note that Burba terms are forced to zero
                    !> if analyzer is /= LI-7500
                    Flux2%E = (1d0 + mu * lEx%sigma) * Flux1%E &
                            + (1d0 + mu * lEx%sigma) &
                            * (Flux1%H + lEx%Burba%h_top + lEx%Burba%h_bot + lEx%Burba%h_spar) &
                            * lEx%RHO%w / (lEx%RhoCp * lEx%Ta)
            else
                Flux2%E = error
            end if
        else
            !> Closed-path uses Ibrom et al. (2007) if conversion to mixing
            !> ratio did not already occur (which implies that some variables
            !> were missing)
            select case(lEx%measure_type(h2o))
                case ('molar_density', 'mole_fraction')
                    if (Flux1%E /= error .and. lEx%sigma /= error &
                        .and. lEx%Vcell(h2o) > 0d0 .and. lEx%Va > 0d0) then

                        if (Flux1%Hi_h2o /= error &
                            .and. lEx%cov_w(pi) /= error) then
                            !> Complete formulation, should actually never be
                            !> used cause conversion to mixing ratio should have
                            !> already happened if everything is available
                            Flux2%E = (1d0 + mu * lEx%sigma) * Flux1%E &
                                * lEx%Vcell(h2o) / lEx%Va &
                                + (1d0 + mu * lEx%sigma) * Flux1%Hi_h2o &
                                * lEx%RHO%w / (lEx%RhoCp * lEx%Tcell) &
                                - (1d0 + mu * lEx%sigma) * lEx%cov_w(pi) &
                                * lEx%RHO%w / (lEx%Pcell)

                        elseif (Flux1%Hi_h2o /= error) then
                            !> Correct only for effect of T
                            Flux2%E = (1d0 + mu * lEx%sigma) * Flux1%E &
                                * lEx%Vcell(h2o) / lEx%Va &
                                + (1d0 + mu * lEx%sigma) * Flux1%Hi_h2o &
                                * lEx%RHO%w / (lEx%RhoCp * lEx%Tcell)

                        elseif (lEx%cov_w(pi)  /= error) then
                            !> Correct only for effect of P
                            Flux2%E = (1d0 + mu * lEx%sigma) * Flux1%E &
                                * lEx%Vcell(h2o) / lEx%Va &
                                - (1d0 + mu * lEx%sigma) * lEx%cov_w(pi) &
                                * lEx%RHO%w / (lEx%Pcell)
                        else
                            !> Can't correct for T and P
                            Flux2%E = Flux1%E * lEx%Vcell(h2o) / lEx%Va
                        end if
                    else
                        Flux2%E = error
                    end if
                case ('mixing_ratio')
                    Flux2%E = Flux1%E
            end select
        end if
    else
        !> If WPL should not be applied
        Flux2%E = Flux1%E
    end if
    !> If Flux1 was error, then set also 2 to error
    if (Flux1%E == error) Flux2%E = error

    !> Level 2 h2o and latent heat flux
    if (Flux2%E /= error) then
        Flux2%h2o = Flux2%E * 1d3 / MW(h2o)
        Flux2%ET = Flux2%h2o * h2o_to_ET
        if (lEx%lambda /= error) then
            Flux2%LE = Flux2%E * lEx%lambda
        else
            Flux2%LE = error
        end if
    else
        Flux2%h2o = error
        Flux2%LE  = error
        Flux2%ET  = error
    end if

    !> Level 2 evapotranspiration fluxes with H2O covariances
    !> at time-lags of other scalars. Do nothing, WPL is deleterious here
    Flux2%E_co2 = Flux1%E_co2
    Flux2%E_ch4 = Flux1%E_ch4
    Flux2%E_gas4 = Flux1%E_gas4

    !> Level 2 Sensible heat
    if (lEx%instr(sonic)%category == 'sonic') then
        !> Corrected for humidity, after Van Dyjk et al. (2004) eq. 3.53
        !> revising Schotanus et al. (1983)
        if (Flux1%H /= error) then
            if(lEx%Flux0%E /= error .and. lEx%cov_w(ts) /= error &
                .and. lEx%RHO%a > 0d0 .and. lEx%Q >= 0d0 &
                .and. lEx%RhoCp > 0d0 .and. alpha /= error) then
                Flux2%H = Flux1%H &
                    - lEx%RhoCp * alpha * lEx%Ts * lEx%Flux0%E / lEx%RHO%a &
                    - lEx%RhoCp * alpha * lEx%Q * lEx%cov_w(ts)
                    !> alternative
                    !- lEx%RhoCp * alpha * lEx%Ta * lEx%Flux0%E / lEx%RHO%a
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
    !if(lEx%Tmap /= error .and. Flux2%H /= error) then
        !Flux2%H  = Flux2%H * lEx%Tmap
    !end if

    !> Level 3 sensible heat, spectral corrected
    if(Flux2%H /= error .and. BPCF%of(w_ts) /= error) then
        Flux3%H = Flux2%H * BPCF%of(w_ts)
    else
        Flux3%H = error
    end if

    !> Level 3 for evapotranspiration: for open path, WPL again with corrected H
    !> Starts again from Level 1 of E, Level 2 was only used to calculate H Level 3.
    if(EddyProProj%wpl .and. lEx%instr(ih2o)%path_type == 'open') then
        if (lEx%RhoCp > 0d0 .and. lEx%Ta > 0d0 .and. Flux1%E /= error &
            .and. Flux1%H /= error .and. lEx%sigma /= error) then
            Flux3%E = (1d0 + mu * lEx%sigma) * Flux1%E &
                + (1d0 + mu * lEx%sigma) &
                * (Flux3%H + lEx%Burba%h_top + lEx%Burba%h_bot + lEx%Burba%h_spar)&
                * lEx%RHO%w / (lEx%RhoCp * lEx%Ta)
        else
            Flux3%E = Flux2%E
        end if
    else
        Flux3%E = Flux2%E
    end if

    !> Level 3 latent heat fluxes with H2O covariances at
    !> timelags of other scalars
    !> Do nothing
    Flux3%E_co2 = Flux2%E_co2
    Flux3%E_ch4 = Flux2%E_ch4
    Flux3%E_gas4 = Flux2%E_gas4

    !> Level 3 h2o flux and latent heat flux
    if (Flux3%E /= error) then
        Flux3%h2o = Flux3%E * 1d3 / MW(h2o)
        Flux3%ET = Flux3%h2o * h2o_to_ET
        if (lEx%lambda /= error) then
            Flux3%LE = Flux3%E * lEx%lambda
        else
            Flux3%LE = error
        end if
    else
        Flux3%h2o = error
        Flux3%LE  = error
        Flux3%ET  = error
    end if

    !> Calculate E_nowpl for closed and open path systems
    if (lEx%instr(ih2o)%path_type == 'closed') then
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

    !> Apply spectral correction to h2o and E/LE Level 3 fluxes for closed path
    if (lEx%instr(ih2o)%path_type == 'closed' .and. Flux3%E /= error) then
        !> Level 3, spectral correction
        Flux3%h2o = Flux3%h2o * BPCF%of(w_h2o)
        Flux3%E   = Flux3%E   * BPCF%of(w_h2o)
        Flux3%LE  = Flux3%LE  * BPCF%of(w_h2o)
        Flux3%ET  = Flux3%ET  * BPCF%of(w_h2o)
    end if

    if (.not. lEx%var_present(h2o)) then
        Flux3%h2o = error
        Flux3%E   = error
        Flux3%LE  = error
        Flux3%ET  = error
    end if

    !> Level 2 other gases
    !> co2
    if (lEx%instr(ico2)%path_type == 'closed') then
        !> Level 2, WPL for closed path, implemented after Ibrom et al. (2007),
        !> Tellus, eq. 3a with H contribution from WPL24
        select case(lEx%measure_type(co2))
            !> Analitically, it is verified that:
            !> (E * mu * sigma / rho%w * co2_density)
            !> equals
            !> (r_c * w'chi_w' / Va)
            case('molar_density')

                !> E, T and P effects (should never be actually used)
                if (Flux1%co2 /= error .and. lEx%sigma >= 0d0 .and. lEx%Va > 0d0 .and. lEx%chi(co2) > 0d0 &
                    .and. Flux3%E_co2 /= error .and. lEx%RHO%w > 0d0 &
                    .and. Flux3%Hi_co2 /= error .and. lEx%RhoCp > 0d0 .and. lEx%Tcell > 0d0 &
                    .and. lEx%cov_w(pi) /= error .and. lEx%Pcell /= error) then
                    Flux2%co2 = Flux1%co2 * lEx%Vcell(co2) / lEx%Va &
                        + Flux3%E_co2 * mu * lEx%sigma / lEx%RHO%w &
                        * lEx%chi(co2) / lEx%Va &
                        + (1d0 + mu * lEx%sigma) * Flux3%Hi_co2 / (lEx%RhoCp * lEx%Tcell) &
                        * lEx%chi(co2) / lEx%Va &
                        - (1d0 + mu * lEx%sigma) * lEx%cov_w(pi) / (lEx%Pcell) &
                        * lEx%chi(co2) / lEx%Va

                !> Onlt E and T effects
                elseif (Flux1%co2 /= error .and. lEx%sigma >= 0d0 .and. lEx%Va > 0d0 .and. lEx%chi(co2) > 0d0 &
                    .and. Flux3%E_co2 /= error .and. lEx%RHO%w > 0d0 &
                    .and. Flux3%Hi_co2 /= error .and. lEx%RhoCp > 0d0 .and. lEx%Tcell > 0d0) then
                    Flux2%co2 = Flux1%co2 * lEx%Vcell(co2) / lEx%Va &
                        + Flux3%E_co2 * mu * lEx%sigma / lEx%RHO%w &
                        * lEx%chi(co2) / lEx%Va &
                        + (1d0 + mu * lEx%sigma) * Flux3%Hi_co2 / (lEx%RhoCp * lEx%Tcell) &
                        * lEx%chi(co2) / lEx%Va

                !> Onlt T and P effects
                elseif (Flux1%co2 /= error .and. lEx%sigma >= 0d0 .and. lEx%Va > 0d0 .and. lEx%chi(co2) > 0d0 &
                    .and. Flux3%Hi_co2 /= error .and. lEx%RhoCp > 0d0 .and. lEx%Tcell > 0d0 &
                    .and. lEx%cov_w(pi) /= error .and. lEx%Pcell /= error) then
                    Flux2%co2 = Flux1%co2 * lEx%Vcell(co2) / lEx%Va &
                        + (1d0 + mu * lEx%sigma) * Flux3%Hi_co2 / (lEx%RhoCp * lEx%Tcell) &
                        * lEx%chi(co2) / lEx%Va &
                        - (1d0 + mu * lEx%sigma) * lEx%cov_w(pi) / (lEx%Pcell) &
                        * lEx%chi(co2) / lEx%Va

                !> Only E and P effects
                elseif (Flux1%co2 /= error .and. lEx%sigma >= 0d0 .and. lEx%Va > 0d0 .and. lEx%chi(co2) > 0d0 &
                    .and. Flux3%E_co2 /= error .and. lEx%RHO%w > 0d0 &
                    .and. lEx%cov_w(pi) /= error .and. lEx%Pcell /= error) then
                    Flux2%co2 = Flux1%co2 * lEx%Vcell(co2) / lEx%Va &
                        + Flux3%E_co2 * mu * lEx%sigma / lEx%RHO%w &
                        * lEx%chi(co2) / lEx%Va &
                        - (1d0 + mu * lEx%sigma) * lEx%cov_w(pi) / (lEx%Pcell) &
                        * lEx%chi(co2) / lEx%Va

                !> Onlt E effects
                elseif (Flux1%co2 /= error .and. lEx%sigma >= 0d0 .and. lEx%Va > 0d0 .and. lEx%chi(co2) > 0d0 &
                    .and. Flux3%E_co2 /= error .and. lEx%RHO%w > 0d0 ) then
                    Flux2%co2 = Flux1%co2 * lEx%Vcell(co2) / lEx%Va &
                        + Flux3%E_co2 * mu * lEx%sigma / lEx%RHO%w &
                        * lEx%chi(co2) / lEx%Va

                !> Onlt T effects
                elseif (Flux1%co2 /= error .and. lEx%sigma >= 0d0 .and. lEx%Va > 0d0 .and. lEx%chi(co2) > 0d0 &
                    .and. Flux3%Hi_co2 /= error .and. lEx%RhoCp > 0d0 .and. lEx%Tcell > 0d0) then
                    Flux2%co2 = Flux1%co2 * lEx%Vcell(co2) / lEx%Va &
                        + (1d0 + mu * lEx%sigma) * Flux3%Hi_co2 / (lEx%RhoCp * lEx%Tcell) &
                        * lEx%chi(co2) / lEx%Va

                !> Onlt P effects
                elseif (Flux1%co2 /= error .and. lEx%sigma >= 0d0 .and. lEx%Va > 0d0 .and. lEx%chi(co2) > 0d0 &
                    .and. lEx%cov_w(pi) /= error .and. lEx%Pcell /= error) then
                    Flux2%co2 = Flux1%co2 * lEx%Vcell(co2) / lEx%Va &
                        - (1d0 + mu * lEx%sigma) * lEx%cov_w(pi) / (lEx%Pcell) &
                        * lEx%chi(co2) / lEx%Va

                !> No WPL effects
                elseif (Flux1%co2 /= error) then
                    Flux2%co2 = Flux1%co2 * lEx%Vcell(co2) / lEx%Va
                else
                    Flux2%co2 = error
                end if

            case('mole_fraction')

                !> E, T and P effects (should never be actually used)
                if (Flux1%co2 /= error .and. lEx%sigma >= 0d0 .and. lEx%Va > 0d0 .and. lEx%chi(co2) > 0d0 &
                    .and. Flux3%E_co2 /= error .and. lEx%RHO%w > 0d0 &
                    .and. Flux3%Hi_co2 /= error .and. lEx%RhoCp > 0d0 .and. lEx%Tcell > 0d0 &
                    .and. lEx%cov_w(pi) /= error .and. lEx%Pcell /= error) then
                    Flux2%co2 = Flux1%co2 &
                        + Flux3%E_co2 * mu * lEx%sigma / lEx%RHO%w &
                        * lEx%chi(co2) / lEx%Va &
                        + (1d0 + mu * lEx%sigma) * Flux3%Hi_co2 / (lEx%RhoCp * lEx%Tcell) &
                        * lEx%chi(co2) / lEx%Va &
                        - (1d0 + mu * lEx%sigma) * lEx%cov_w(pi) / (lEx%Pcell) &
                        * lEx%chi(co2) / lEx%Va

                !> Onlt E and T effects
                elseif (Flux1%co2 /= error .and. lEx%sigma >= 0d0 .and. lEx%Va > 0d0 .and. lEx%chi(co2) > 0d0 &
                    .and. Flux3%E_co2 /= error .and. lEx%RHO%w > 0d0 &
                    .and. Flux3%Hi_co2 /= error .and. lEx%RhoCp > 0d0 .and. lEx%Tcell > 0d0) then
                    Flux2%co2 = Flux1%co2 &
                        + Flux3%E_co2 * mu * lEx%sigma / lEx%RHO%w &
                        * lEx%chi(co2) / lEx%Va &
                        + (1d0 + mu * lEx%sigma) * Flux3%Hi_co2 / (lEx%RhoCp * lEx%Tcell) &
                        * lEx%chi(co2) / lEx%Va

                !> Onlt T and P effects
                elseif (Flux1%co2 /= error .and. lEx%sigma >= 0d0 .and. lEx%Va > 0d0 .and. lEx%chi(co2) > 0d0 &
                    .and. Flux3%Hi_co2 /= error .and. lEx%RhoCp > 0d0 .and. lEx%Tcell > 0d0 &
                    .and. lEx%cov_w(pi) /= error .and. lEx%Pcell /= error) then
                    Flux2%co2 = Flux1%co2 &
                        + (1d0 + mu * lEx%sigma) * Flux3%Hi_co2 / (lEx%RhoCp * lEx%Tcell) &
                        * lEx%chi(co2) / lEx%Va &
                        - (1d0 + mu * lEx%sigma) * lEx%cov_w(pi) / (lEx%Pcell) &
                        * lEx%chi(co2) / lEx%Va

                !> Only E and P effects
                elseif (Flux1%co2 /= error .and. lEx%sigma >= 0d0 .and. lEx%Va > 0d0 .and. lEx%chi(co2) > 0d0 &
                    .and. Flux3%E_co2 /= error .and. lEx%RHO%w > 0d0 &
                    .and. lEx%cov_w(pi) /= error .and. lEx%Pcell /= error) then
                    Flux2%co2 = Flux1%co2 &
                        + Flux3%E_co2 * mu * lEx%sigma / lEx%RHO%w &
                        * lEx%chi(co2) / lEx%Va &
                        - (1d0 + mu * lEx%sigma) * lEx%cov_w(pi) / (lEx%Pcell) &
                        * lEx%chi(co2) / lEx%Va

                !> Onlt E effects
                elseif (Flux1%co2 /= error .and. lEx%sigma >= 0d0 .and. lEx%Va > 0d0 .and. lEx%chi(co2) > 0d0 &
                    .and. Flux3%E_co2 /= error .and. lEx%RHO%w > 0d0 ) then
                    Flux2%co2 = Flux1%co2 &
                        + Flux3%E_co2 * mu * lEx%sigma / lEx%RHO%w &
                        * lEx%chi(co2) / lEx%Va

                !> Onlt T effects
                elseif (Flux1%co2 /= error .and. lEx%sigma >= 0d0 .and. lEx%Va > 0d0 .and. lEx%chi(co2) > 0d0 &
                    .and. Flux3%Hi_co2 /= error .and. lEx%RhoCp > 0d0 .and. lEx%Tcell > 0d0) then
                    Flux2%co2 = Flux1%co2 &
                        + (1d0 + mu * lEx%sigma) * Flux3%Hi_co2 / (lEx%RhoCp * lEx%Tcell) &
                        * lEx%chi(co2) / lEx%Va

                !> Onlt P effects
                elseif (Flux1%co2 /= error .and. lEx%sigma >= 0d0 .and. lEx%Va > 0d0 .and. lEx%chi(co2) > 0d0 &
                    .and. lEx%cov_w(pi) /= error .and. lEx%Pcell /= error) then
                    Flux2%co2 = Flux1%co2 &
                        - (1d0 + mu * lEx%sigma) * lEx%cov_w(pi) / (lEx%Pcell) &
                        * lEx%chi(co2) / lEx%Va

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
        if(lEx%instr(ico2)%model(1:len_trim(lEx%instr(ico2)%model) - 2)  == 'li7500') then
            if (Flux1%co2 /= error .and. Flux3%H /= error .and. Flux3%E /= error &
                .and. lEx%RHO%d > 0d0 .and. lEx%RhoCp > 0d0 .and. lEx%Ta > 0d0 .and. lEx%sigma /= error) then
                Flux2%co2 = Flux1%co2 + mu * Flux3%E * lEx%d(co2) * 1d3 &
                    / ((1d0 + mu * lEx%sigma) * lEx%RHO%d) &
                    + (Flux3%H + lEx%Burba%h_top + lEx%Burba%h_bot + lEx%Burba%h_spar) * lEx%d(co2) * 1d3 &
                    / (lEx%RhoCp * lEx%Ta)
            elseif(Flux1%co2 /= error .and. Flux3%H /= error &
                .and. lEx%RhoCp > 0d0 .and. lEx%Ta > 0d0) then
                Flux2%co2 = Flux1%co2 &
                    + (Flux3%H + lEx%Burba%h_top + lEx%Burba%h_bot + lEx%Burba%h_spar) * lEx%d(co2) * 1d3 &
                    / (lEx%RhoCp * lEx%Ta)
            elseif(Flux1%co2 /= error .and. Flux3%E /= error &
                .and. lEx%RHO%d > 0d0 .and. lEx%sigma /= error) then
                Flux2%co2 = Flux1%co2 &
                    + mu * Flux3%E * lEx%d(co2) * 1d3 / ((1d0 + mu * lEx%sigma) * lEx%RHO%d)
            elseif(Flux1%co2 /= error) then
                Flux2%co2 = Flux1%co2
            else
                Flux2%co2 = error
            end if
        else
            if (Flux1%co2 /= error .and. Flux3%H /= error .and. Flux3%E /= error &
                .and. lEx%RHO%d > 0d0 .and. lEx%RhoCp > 0d0 .and. lEx%Ta > 0d0 .and. lEx%sigma /= error) then
                Flux2%co2 = Flux1%co2 &
                    + mu * Flux3%E * lEx%d(co2) * 1d3 / ((1d0 + mu * lEx%sigma) * lEx%RHO%d) &
                    + Flux3%H * lEx%d(co2) * 1d3 / (lEx%RhoCp * lEx%Ta)
            elseif(Flux1%co2 /= error .and. Flux3%H /= error &
                .and. lEx%RhoCp > 0d0 .and. lEx%Ta > 0d0) then
                Flux2%co2 = Flux1%co2 &
                    + Flux3%H * lEx%d(co2) * 1d3 / (lEx%RhoCp * lEx%Ta)
            elseif(Flux1%co2 /= error .and. Flux3%E /= error &
                .and. lEx%RHO%d > 0d0 .and. lEx%sigma /= error) then
                Flux2%co2 = Flux1%co2 &
                    + mu * Flux3%E * lEx%d(co2) * 1d3 / ((1d0 + mu * lEx%sigma) * lEx%RHO%d)
            elseif(Flux1%co2 /= error) then
                Flux2%co2 = Flux1%co2
            else
                Flux2%co2 = error
            end if
        end if
    end if
    if (.not. lEx%var_present(co2)) Flux2%co2 = error

    !> ch4
    if (lEx%instr(ich4)%path_type == 'closed') then
        !> Level 2, WPL for closed path, implemented after Ibrom et al. (2007),
        !> Tellus, eq. 3a with H contribution from WPL24
        select case(lEx%measure_type(ch4))
            case('molar_density')

                !> E, T and P effects (should never be actually used)
                if (Flux1%ch4 /= error .and. lEx%sigma >= 0d0 .and. lEx%Va > 0d0 .and. lEx%chi(ch4) > 0d0 &
                    .and. Flux3%E_ch4 /= error .and. lEx%RHO%w > 0d0 &
                    .and. Flux3%Hi_ch4 /= error .and. lEx%RhoCp > 0d0 .and. lEx%Tcell > 0d0 &
                    .and. lEx%cov_w(pi) /= error .and. lEx%Pcell /= error) then
                    Flux2%ch4 = Flux1%ch4 * lEx%Vcell(ch4) / lEx%Va &
                        + Flux3%E_ch4 * mu * lEx%sigma / lEx%RHO%w &
                        * lEx%chi(ch4) / lEx%Va &
                        + (1d0 + mu * lEx%sigma) * Flux3%Hi_ch4 / (lEx%RhoCp * lEx%Tcell) &
                        * lEx%chi(ch4) / lEx%Va &
                        - (1d0 + mu * lEx%sigma) * lEx%cov_w(pi) / (lEx%Pcell) &
                        * lEx%chi(ch4) / lEx%Va

                !> Onlt E and T effects
                elseif (Flux1%ch4 /= error .and. lEx%sigma >= 0d0 .and. lEx%Va > 0d0 .and. lEx%chi(ch4) > 0d0 &
                    .and. Flux3%E_ch4 /= error .and. lEx%RHO%w > 0d0 &
                    .and. Flux3%Hi_ch4 /= error .and. lEx%RhoCp > 0d0 .and. lEx%Tcell > 0d0) then
                    Flux2%ch4 = Flux1%ch4 * lEx%Vcell(ch4) / lEx%Va &
                        + Flux3%E_ch4 * mu * lEx%sigma / lEx%RHO%w &
                        * lEx%chi(ch4) / lEx%Va &
                        + (1d0 + mu * lEx%sigma) * Flux3%Hi_ch4 / (lEx%RhoCp * lEx%Tcell) &
                        * lEx%chi(ch4) / lEx%Va

                !> Onlt T and P effects
                elseif (Flux1%ch4 /= error .and. lEx%sigma >= 0d0 .and. lEx%Va > 0d0 .and. lEx%chi(ch4) > 0d0 &
                    .and. Flux3%Hi_ch4 /= error .and. lEx%RhoCp > 0d0 .and. lEx%Tcell > 0d0 &
                    .and. lEx%cov_w(pi) /= error .and. lEx%Pcell /= error) then
                    Flux2%ch4 = Flux1%ch4 * lEx%Vcell(ch4) / lEx%Va &
                        + (1d0 + mu * lEx%sigma) * Flux3%Hi_ch4 / (lEx%RhoCp * lEx%Tcell) &
                        * lEx%chi(ch4) / lEx%Va &
                        - (1d0 + mu * lEx%sigma) * lEx%cov_w(pi) / (lEx%Pcell) &
                        * lEx%chi(ch4) / lEx%Va

                !> Only E and P effects
                elseif (Flux1%ch4 /= error .and. lEx%sigma >= 0d0 .and. lEx%Va > 0d0 .and. lEx%chi(ch4) > 0d0 &
                    .and. Flux3%E_ch4 /= error .and. lEx%RHO%w > 0d0 &
                    .and. lEx%cov_w(pi) /= error .and. lEx%Pcell /= error) then
                    Flux2%ch4 = Flux1%ch4 * lEx%Vcell(ch4) / lEx%Va &
                        + Flux3%E_ch4 * mu * lEx%sigma / lEx%RHO%w &
                        * lEx%chi(ch4) / lEx%Va &
                        - (1d0 + mu * lEx%sigma) * lEx%cov_w(pi) / (lEx%Pcell) &
                        * lEx%chi(ch4) / lEx%Va

                !> Onlt E effects
                elseif (Flux1%ch4 /= error .and. lEx%sigma >= 0d0 .and. lEx%Va > 0d0 .and. lEx%chi(ch4) > 0d0 &
                    .and. Flux3%E_ch4 /= error .and. lEx%RHO%w > 0d0 ) then
                    Flux2%ch4 = Flux1%ch4 * lEx%Vcell(ch4) / lEx%Va &
                        + Flux3%E_ch4 * mu * lEx%sigma / lEx%RHO%w &
                        * lEx%chi(ch4) / lEx%Va

                !> Onlt T effects
                elseif (Flux1%ch4 /= error .and. lEx%sigma >= 0d0 .and. lEx%Va > 0d0 .and. lEx%chi(ch4) > 0d0 &
                    .and. Flux3%Hi_ch4 /= error .and. lEx%RhoCp > 0d0 .and. lEx%Tcell > 0d0) then
                    Flux2%ch4 = Flux1%ch4 * lEx%Vcell(ch4) / lEx%Va &
                        + (1d0 + mu * lEx%sigma) * Flux3%Hi_ch4 / (lEx%RhoCp * lEx%Tcell) &
                        * lEx%chi(ch4) / lEx%Va

                !> Onlt P effects
                elseif (Flux1%ch4 /= error .and. lEx%sigma >= 0d0 .and. lEx%Va > 0d0 .and. lEx%chi(ch4) > 0d0 &
                    .and. lEx%cov_w(pi) /= error .and. lEx%Pcell /= error) then
                    Flux2%ch4 = Flux1%ch4 * lEx%Vcell(ch4) / lEx%Va &
                        - (1d0 + mu * lEx%sigma) * lEx%cov_w(pi) / (lEx%Pcell) &
                        * lEx%chi(ch4) / lEx%Va

                !> No WPL effects
                elseif (Flux1%ch4 /= error) then
                    Flux2%ch4 = Flux1%ch4 * lEx%Vcell(ch4) / lEx%Va
                else
                    Flux2%ch4 = error
                end if

            case('mole_fraction')

                !> E, T and P effects
                if (Flux1%ch4 /= error .and. lEx%sigma >= 0d0 .and. lEx%Va > 0d0 .and. lEx%chi(ch4) > 0d0 &
                    .and. Flux3%E_ch4 /= error .and. lEx%RHO%w > 0d0 &
                    .and. Flux3%Hi_ch4 /= error .and. lEx%RhoCp > 0d0 .and. lEx%Tcell > 0d0 &
                    .and. lEx%cov_w(pi) /= error .and. lEx%Pcell /= error) then
                    Flux2%ch4 = Flux1%ch4 &
                        + Flux3%E_ch4 * mu * lEx%sigma / lEx%RHO%w &
                        * lEx%chi(ch4) / lEx%Va &
                        + (1d0 + mu * lEx%sigma) * Flux3%Hi_ch4 / (lEx%RhoCp * lEx%Tcell) &
                        * lEx%chi(ch4) / lEx%Va &
                        - (1d0 + mu * lEx%sigma) * lEx%cov_w(pi) / (lEx%Pcell) &
                        * lEx%chi(ch4) / lEx%Va

                !> Onlt E and T effects
                elseif (Flux1%ch4 /= error .and. lEx%sigma >= 0d0 .and. lEx%Va > 0d0 .and. lEx%chi(ch4) > 0d0 &
                    .and. Flux3%E_ch4 /= error .and. lEx%RHO%w > 0d0 &
                    .and. Flux3%Hi_ch4 /= error .and. lEx%RhoCp > 0d0 .and. lEx%Tcell > 0d0) then
                    Flux2%ch4 = Flux1%ch4 &
                        + Flux3%E_ch4 * mu * lEx%sigma / lEx%RHO%w &
                        * lEx%chi(ch4) / lEx%Va &
                        + (1d0 + mu * lEx%sigma) * Flux3%Hi_ch4 / (lEx%RhoCp * lEx%Tcell) &
                        * lEx%chi(ch4) / lEx%Va

                !> Onlt T and P effects
                elseif (Flux1%ch4 /= error .and. lEx%sigma >= 0d0 .and. lEx%Va > 0d0 .and. lEx%chi(ch4) > 0d0 &
                    .and. Flux3%Hi_ch4 /= error .and. lEx%RhoCp > 0d0 .and. lEx%Tcell > 0d0 &
                    .and. lEx%cov_w(pi) /= error .and. lEx%Pcell /= error) then
                    Flux2%ch4 = Flux1%ch4 &
                        + (1d0 + mu * lEx%sigma) * Flux3%Hi_ch4 / (lEx%RhoCp * lEx%Tcell) &
                        * lEx%chi(ch4) / lEx%Va &
                        - (1d0 + mu * lEx%sigma) * lEx%cov_w(pi) / (lEx%Pcell) &
                        * lEx%chi(ch4) / lEx%Va

                !> Only E and P effects
                elseif (Flux1%ch4 /= error .and. lEx%sigma >= 0d0 .and. lEx%Va > 0d0 .and. lEx%chi(ch4) > 0d0 &
                    .and. Flux3%E_ch4 /= error .and. lEx%RHO%w > 0d0 &
                    .and. lEx%cov_w(pi) /= error .and. lEx%Pcell /= error) then
                    Flux2%ch4 = Flux1%ch4 &
                        + Flux3%E_ch4 * mu * lEx%sigma / lEx%RHO%w &
                        * lEx%chi(ch4) / lEx%Va &
                        - (1d0 + mu * lEx%sigma) * lEx%cov_w(pi) / (lEx%Pcell) &
                        * lEx%chi(ch4) / lEx%Va

                !> Onlt E effects
                elseif (Flux1%ch4 /= error .and. lEx%sigma >= 0d0 .and. lEx%Va > 0d0 .and. lEx%chi(ch4) > 0d0 &
                    .and. Flux3%E_ch4 /= error .and. lEx%RHO%w > 0d0 ) then
                    Flux2%ch4 = Flux1%ch4 &
                        + Flux3%E_ch4 * mu * lEx%sigma / lEx%RHO%w &
                        * lEx%chi(ch4) / lEx%Va

                !> Onlt T effects
                elseif (Flux1%ch4 /= error .and. lEx%sigma >= 0d0 .and. lEx%Va > 0d0 .and. lEx%chi(ch4) > 0d0 &
                    .and. Flux3%Hi_ch4 /= error .and. lEx%RhoCp > 0d0 .and. lEx%Tcell > 0d0) then
                    Flux2%ch4 = Flux1%ch4 &
                        + (1d0 + mu * lEx%sigma) * Flux3%Hi_ch4 / (lEx%RhoCp * lEx%Tcell) &
                        * lEx%chi(ch4) / lEx%Va

                !> Onlt P effects
                elseif (Flux1%ch4 /= error .and. lEx%sigma >= 0d0 .and. lEx%Va > 0d0 .and. lEx%chi(ch4) > 0d0 &
                    .and. lEx%cov_w(pi) /= error .and. lEx%Pcell /= error) then
                    Flux2%ch4 = Flux1%ch4 &
                        - (1d0 + mu * lEx%sigma) * lEx%cov_w(pi) / (lEx%Pcell) &
                        * lEx%chi(ch4) / lEx%Va

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
        if (lEx%instr(ich4)%model(1:len_trim(lEx%instr(ich4)%model) - 2)  == 'li7700') then
            !> Flux formulation including multipliers requires the use of the original WPL formulation
            !> making use of E not corrected for WPL (E_nowpl) and H, both spectrally corrected.
            if (Flux1%ch4 /= error .and. Flux3%H /= error .and. E_nowpl /= error &
                .and. lEx%RHO%d > 0d0 .and. lEx%RhoCp > 0d0 .and. lEx%Ta > 0d0 .and. lEx%sigma /= error) then
                Flux2%ch4 = lEx%Mul7700%A *(Flux1%ch4 &
                          + lEx%Mul7700%B * mu * lEx%d(ch4) * 1d3 * E_nowpl / lEx%RHO%d &
                          + lEx%Mul7700%C * (1d0 + mu * lEx%sigma) * Flux3%H * lEx%d(ch4) * 1d3 / (lEx%RhoCp * lEx%Ta))
            elseif(Flux1%ch4 /= error .and. Flux3%H /= error &
                .and. lEx%RhoCp > 0d0 .and. lEx%Ta > 0d0 .and. lEx%sigma /= error) then
                Flux2%ch4 = lEx%Mul7700%A *(Flux1%ch4 &
                          + lEx%Mul7700%C * (1d0 + mu * lEx%sigma) * Flux3%H * lEx%d(ch4) * 1d3 / (lEx%RhoCp * lEx%Ta))
            elseif(Flux1%ch4 /= error .and. E_nowpl /= error .and. lEx%RHO%d > 0d0) then
                Flux2%ch4 = lEx%Mul7700%A *(Flux1%ch4 &
                          + lEx%Mul7700%B * mu * lEx%d(ch4) * 1d3 * E_nowpl / lEx%RHO%d)
            elseif (Flux1%ch4 /= error) then
                Flux2%ch4 = lEx%Mul7700%A * Flux1%ch4
            else
                Flux2%ch4 = error
            end if
        else
            !> Flux formulation without multipliers follows Burba et al. 2008.
            if (Flux1%ch4 /= error .and. Flux3%H /= error .and. Flux3%E /= error &
                .and. lEx%RHO%d > 0d0 .and. lEx%RhoCp > 0d0 .and. lEx%Ta > 0d0 .and. lEx%sigma /= error) then
                Flux2%ch4 = Flux1%ch4 &
                    + mu * Flux3%E * lEx%d(ch4) * 1d3 / ((1d0 + mu * lEx%sigma) * lEx%RHO%d) &
                    + Flux3%H * lEx%d(ch4) * 1d3 / (lEx%RhoCp * lEx%Ta)
            elseif(Flux1%ch4 /= error .and. Flux3%H /= error &
                .and. lEx%RhoCp > 0d0 .and. lEx%Ta > 0d0) then
                Flux2%ch4 = Flux1%ch4 &
                    + Flux3%H * lEx%d(ch4) * 1d3 / (lEx%RhoCp * lEx%Ta)
            elseif(Flux1%ch4 /= error .and. Flux3%E /= error &
                .and. lEx%RHO%d > 0d0 .and. lEx%sigma /= error) then
                Flux2%ch4 = Flux1%ch4 &
                    + mu * Flux3%E * lEx%d(ch4) * 1d3 / ((1d0 + mu * lEx%sigma) * lEx%RHO%d)
            elseif(Flux1%ch4 /= error) then
                Flux2%ch4 = Flux1%ch4
            end if
        end if
    end if
    if (.not. lEx%var_present(ch4)) Flux2%ch4 = error

    !> gas4
    if (lEx%instr(igas4)%path_type == 'closed') then
        !> Level 2, WPL for closed path, implemented after Ibrom et al. (2007),
        !> Tellus, eq. 3a with H contribution from WPL24
        select case(lEx%measure_type(gas4))
            case('molar_density')

                !> E, T and P effects (should never be actually used)
                if (Flux1%gas4 /= error .and. lEx%sigma >= 0d0 .and. lEx%Va > 0d0 .and. lEx%chi(n2o) > 0d0 &
                    .and. Flux3%E_gas4 /= error .and. lEx%RHO%w > 0d0 &
                    .and. Flux3%Hi_gas4 /= error .and. lEx%RhoCp > 0d0 .and. lEx%Tcell > 0d0 &
                    .and. lEx%cov_w(pi) /= error .and. lEx%Pcell /= error) then
                    Flux2%gas4 = Flux1%gas4 * lEx%Vcell(n2o) / lEx%Va &
                        + Flux3%E_gas4 * mu * lEx%sigma / lEx%RHO%w &
                        * lEx%chi(n2o) / lEx%Va &
                        + (1d0 + mu * lEx%sigma) * Flux3%Hi_gas4 / (lEx%RhoCp * lEx%Tcell) &
                        * lEx%chi(n2o) / lEx%Va &
                        - (1d0 + mu * lEx%sigma) * lEx%cov_w(pi) / (lEx%Pcell) &
                        * lEx%chi(n2o) / lEx%Va

                !> Onlt E and T effects
                elseif (Flux1%gas4 /= error .and. lEx%sigma >= 0d0 .and. lEx%Va > 0d0 .and. lEx%chi(n2o) > 0d0 &
                    .and. Flux3%E_gas4 /= error .and. lEx%RHO%w > 0d0 &
                    .and. Flux3%Hi_gas4 /= error .and. lEx%RhoCp > 0d0 .and. lEx%Tcell > 0d0) then
                    Flux2%gas4 = Flux1%gas4 * lEx%Vcell(n2o) / lEx%Va &
                        + Flux3%E_gas4 * mu * lEx%sigma / lEx%RHO%w &
                        * lEx%chi(n2o) / lEx%Va &
                        + (1d0 + mu * lEx%sigma) * Flux3%Hi_gas4 / (lEx%RhoCp * lEx%Tcell) &
                        * lEx%chi(n2o) / lEx%Va

                !> Onlt T and P effects
                elseif (Flux1%gas4 /= error .and. lEx%sigma >= 0d0 .and. lEx%Va > 0d0 .and. lEx%chi(n2o) > 0d0 &
                    .and. Flux3%Hi_gas4 /= error .and. lEx%RhoCp > 0d0 .and. lEx%Tcell > 0d0 &
                    .and. lEx%cov_w(pi) /= error .and. lEx%Pcell /= error) then
                    Flux2%gas4 = Flux1%gas4 * lEx%Vcell(n2o) / lEx%Va &
                        + (1d0 + mu * lEx%sigma) * Flux3%Hi_gas4 / (lEx%RhoCp * lEx%Tcell) &
                        * lEx%chi(n2o) / lEx%Va &
                        - (1d0 + mu * lEx%sigma) * lEx%cov_w(pi) / (lEx%Pcell) &
                        * lEx%chi(n2o) / lEx%Va

                !> Only E and P effects
                elseif (Flux1%gas4 /= error .and. lEx%sigma >= 0d0 .and. lEx%Va > 0d0 .and. lEx%chi(n2o) > 0d0 &
                    .and. Flux3%E_gas4 /= error .and. lEx%RHO%w > 0d0 &
                    .and. lEx%cov_w(pi) /= error .and. lEx%Pcell /= error) then
                    Flux2%gas4 = Flux1%gas4 * lEx%Vcell(n2o) / lEx%Va &
                        + Flux3%E_gas4 * mu * lEx%sigma / lEx%RHO%w &
                        * lEx%chi(n2o) / lEx%Va &
                        - (1d0 + mu * lEx%sigma) * lEx%cov_w(pi) / (lEx%Pcell) &
                        * lEx%chi(n2o) / lEx%Va

                !> Onlt E effects
                elseif (Flux1%gas4 /= error .and. lEx%sigma >= 0d0 .and. lEx%Va > 0d0 .and. lEx%chi(n2o) > 0d0 &
                    .and. Flux3%E_gas4 /= error .and. lEx%RHO%w > 0d0 ) then
                    Flux2%gas4 = Flux1%gas4 * lEx%Vcell(n2o) / lEx%Va &
                        + Flux3%E_gas4 * mu * lEx%sigma / lEx%RHO%w &
                        * lEx%chi(n2o) / lEx%Va

                !> Onlt T effects
                elseif (Flux1%gas4 /= error .and. lEx%sigma >= 0d0 .and. lEx%Va > 0d0 .and. lEx%chi(n2o) > 0d0 &
                    .and. Flux3%Hi_gas4 /= error .and. lEx%RhoCp > 0d0 .and. lEx%Tcell > 0d0) then
                    Flux2%gas4 = Flux1%gas4 * lEx%Vcell(n2o) / lEx%Va &
                        + (1d0 + mu * lEx%sigma) * Flux3%Hi_gas4 / (lEx%RhoCp * lEx%Tcell) &
                        * lEx%chi(n2o) / lEx%Va

                !> Onlt P effects
                elseif (Flux1%gas4 /= error .and. lEx%sigma >= 0d0 .and. lEx%Va > 0d0 .and. lEx%chi(n2o) > 0d0 &
                    .and. lEx%cov_w(pi) /= error .and. lEx%Pcell /= error) then
                    Flux2%gas4 = Flux1%gas4 * lEx%Vcell(n2o) / lEx%Va &
                        - (1d0 + mu * lEx%sigma) * lEx%cov_w(pi) / (lEx%Pcell) &
                        * lEx%chi(n2o) / lEx%Va

                !> No WPL effects
                elseif (Flux1%gas4 /= error) then
                    Flux2%gas4 = Flux1%gas4 * lEx%Vcell(n2o) / lEx%Va
                else
                    Flux2%gas4 = error
                end if

            case('mole_fraction')

                !> E, T and P effects (should never be actually used)
                if (Flux1%gas4 /= error .and. lEx%sigma >= 0d0 .and. lEx%Va > 0d0 .and. lEx%chi(n2o) > 0d0 &
                    .and. Flux3%E_gas4 /= error .and. lEx%RHO%w > 0d0 &
                    .and. Flux3%Hi_gas4 /= error .and. lEx%RhoCp > 0d0 .and. lEx%Tcell > 0d0 &
                    .and. lEx%cov_w(pi) /= error .and. lEx%Pcell /= error) then
                    Flux2%gas4 = Flux1%gas4 &
                        + Flux3%E_gas4 * mu * lEx%sigma / lEx%RHO%w &
                        * lEx%chi(n2o) / lEx%Va &
                        + (1d0 + mu * lEx%sigma) * Flux3%Hi_gas4 / (lEx%RhoCp * lEx%Tcell) &
                        * lEx%chi(n2o) / lEx%Va &
                        - (1d0 + mu * lEx%sigma) * lEx%cov_w(pi) / (lEx%Pcell) &
                        * lEx%chi(n2o) / lEx%Va

                !> Onlt E and T effects
                elseif (Flux1%gas4 /= error .and. lEx%sigma >= 0d0 .and. lEx%Va > 0d0 .and. lEx%chi(n2o) > 0d0 &
                    .and. Flux3%E_gas4 /= error .and. lEx%RHO%w > 0d0 &
                    .and. Flux3%Hi_gas4 /= error .and. lEx%RhoCp > 0d0 .and. lEx%Tcell > 0d0) then
                    Flux2%gas4 = Flux1%gas4 &
                        + Flux3%E_gas4 * mu * lEx%sigma / lEx%RHO%w &
                        * lEx%chi(n2o) / lEx%Va &
                        + (1d0 + mu * lEx%sigma) * Flux3%Hi_gas4 / (lEx%RhoCp * lEx%Tcell) &
                        * lEx%chi(n2o) / lEx%Va

                !> Onlt T and P effects
                elseif (Flux1%gas4 /= error .and. lEx%sigma >= 0d0 .and. lEx%Va > 0d0 .and. lEx%chi(n2o) > 0d0 &
                    .and. Flux3%Hi_gas4 /= error .and. lEx%RhoCp > 0d0 .and. lEx%Tcell > 0d0 &
                    .and. lEx%cov_w(pi) /= error .and. lEx%Pcell /= error) then
                    Flux2%gas4 = Flux1%gas4 &
                        + (1d0 + mu * lEx%sigma) * Flux3%Hi_gas4 / (lEx%RhoCp * lEx%Tcell) &
                        * lEx%chi(n2o) / lEx%Va &
                        - (1d0 + mu * lEx%sigma) * lEx%cov_w(pi) / (lEx%Pcell) &
                        * lEx%chi(n2o) / lEx%Va

                !> Only E and P effects
                elseif (Flux1%gas4 /= error .and. lEx%sigma >= 0d0 .and. lEx%Va > 0d0 .and. lEx%chi(n2o) > 0d0 &
                    .and. Flux3%E_gas4 /= error .and. lEx%RHO%w > 0d0 &
                    .and. lEx%cov_w(pi) /= error .and. lEx%Pcell /= error) then
                    Flux2%gas4 = Flux1%gas4 &
                        + Flux3%E_gas4 * mu * lEx%sigma / lEx%RHO%w &
                        * lEx%chi(n2o) / lEx%Va &
                        - (1d0 + mu * lEx%sigma) * lEx%cov_w(pi) / (lEx%Pcell) &
                        * lEx%chi(n2o) / lEx%Va

                !> Onlt E effects
                elseif (Flux1%gas4 /= error .and. lEx%sigma >= 0d0 .and. lEx%Va > 0d0 .and. lEx%chi(n2o) > 0d0 &
                    .and. Flux3%E_gas4 /= error .and. lEx%RHO%w > 0d0 ) then
                    Flux2%gas4 = Flux1%gas4 &
                        + Flux3%E_gas4 * mu * lEx%sigma / lEx%RHO%w &
                        * lEx%chi(n2o) / lEx%Va

                !> Onlt T effects
                elseif (Flux1%gas4 /= error .and. lEx%sigma >= 0d0 .and. lEx%Va > 0d0 .and. lEx%chi(n2o) > 0d0 &
                    .and. Flux3%Hi_gas4 /= error .and. lEx%RhoCp > 0d0 .and. lEx%Tcell > 0d0) then
                    Flux2%gas4 = Flux1%gas4 &
                        + (1d0 + mu * lEx%sigma) * Flux3%Hi_gas4 / (lEx%RhoCp * lEx%Tcell) &
                        * lEx%chi(n2o) / lEx%Va

                !> Onlt P effects
                elseif (Flux1%gas4 /= error .and. lEx%sigma >= 0d0 .and. lEx%Va > 0d0 .and. lEx%chi(n2o) > 0d0 &
                    .and. lEx%cov_w(pi) /= error .and. lEx%Pcell /= error) then
                    Flux2%gas4 = Flux1%gas4 &
                        - (1d0 + mu * lEx%sigma) * lEx%cov_w(pi) / (lEx%Pcell) &
                        * lEx%chi(n2o) / lEx%Va

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
            .and. lEx%RHO%d > 0d0 .and. lEx%RhoCp > 0d0 .and. lEx%Ta > 0d0 .and. lEx%sigma /= error) then
            Flux2%gas4 = Flux1%gas4 &
                + mu * Flux3%E * lEx%d(gas4) * 1d3 / ((1d0 + mu * lEx%sigma) * lEx%RHO%d) &
                + Flux3%H * lEx%d(gas4) * 1d3 / (lEx%RhoCp * lEx%Ta)
        elseif(Flux1%gas4 /= error .and. Flux3%H /= error &
            .and. lEx%RhoCp > 0d0 .and. lEx%Ta > 0d0) then
            Flux2%gas4 = Flux1%gas4 &
                + Flux3%H * lEx%d(gas4) * 1d3 / (lEx%RhoCp * lEx%Ta)
        elseif(Flux1%gas4 /= error .and. Flux3%E /= error &
            .and. lEx%RHO%d > 0d0 .and. lEx%sigma /= error) then
            Flux2%gas4 = Flux1%gas4 &
                + mu * Flux3%E * lEx%d(gas4) * 1d3 / ((1d0 + mu * lEx%sigma) * lEx%RHO%d)
        elseif(Flux1%gas4 /= error) then
            Flux2%gas4 = Flux1%gas4
        else
            Flux2%gas4 = error
        end if
    end if
    if (.not. lEx%var_present(gas4)) Flux2%gas4 = error

    !> If WPL should not be applied
    if (.not. EddyProProj%wpl) then
        Flux2%co2 = Flux1%co2
        Flux2%ch4 = Flux1%ch4
        Flux2%gas4 = Flux1%gas4
    end if

    !> Level 3 other gases. For closed path apply now the spectral correction (e.g. Ibrom et al. 2007)
    !> co2
    if (lEx%instr(ico2)%path_type == 'closed' .and. Flux2%co2 /= error) then
        !> Level 3, spectral correction
        Flux3%co2 = Flux2%co2 * BPCF%of(w_co2)
    else
        !> Level 3, spectral correction was already applied
        Flux3%co2 = Flux2%co2
    end if
    !> ch4
    if (lEx%instr(ich4)%path_type == 'closed' .and. Flux2%ch4 /= error) then
        !> Level 3, spectral correction
        Flux3%ch4 = Flux2%ch4 * BPCF%of(w_ch4)
    else
        !> Level 3, spectral correction was already applied
        Flux3%ch4 = Flux2%ch4
    end if
    !> n2o
    if (lEx%instr(igas4)%path_type == 'closed' .and. Flux2%gas4 /= error) then
        !> Level 3, spectral correction
        Flux3%gas4 = Flux2%gas4 * BPCF%of(w_gas4)
    else
        !> Level 3, spectral correction was already applied
        Flux3%gas4 = Flux2%gas4
    end if

    !> Potential temperature
    !> If condition fails, previous value (from Fluxes0) holds for z/Ambient%L
    if (lEx%Pa > 0d0) then
        Tp = lEx%Ta * (1d5 / lEx%Pa)**(.286d0)
    else
        Tp = error
    end if

    !> Momentum flux and friction velocity
    Flux2%tau = Flux1%tau
    Flux3%tau = Flux1%tau
    Flux2%ustar = Flux1%ustar
    Flux3%ustar = Flux1%ustar

    !> Monin-Obukhov length (L = - (Tp^ /(k*g))*(ustar**3/(w'Tp')^ in m)
    if (Flux3%H /= 0d0 .and. Flux3%H /= error .and. &
        lEx%RhoCp > 0d0 .and. Flux3%ustar >= 0d0 .and. Tp > 0d0) then
        lEx%L = -Tp * (Flux3%ustar**3) / (vk * g * Flux3%H / lEx%RhoCp)
    else
        lEx%L = error
    end if

    !> Monin-Obukhov stability parameter (zL = z/L)
    !> If condition fails, previous value (from Fluxes0) holds
    if (lEx%L /= 0d0 .and. lEx%L /= error) &
        lEx%zL = (lEx%instr(sonic)%height - lEx%disp_height) / lEx%L

    !> scale temperature(T*)
    !> If condition fails, previous value (from Fluxes0) holds
    if (Flux3%ustar > 0d0 .and. Flux3%H /= error .and. lEx%RhoCp > 0d0) &
        lEx%Tstar = Flux3%H / (lEx%RhoCp * Flux3%ustar)

    !> Bowen ration (Bowen, 1926, Phyis Rev)
    if (Flux3%LE /= 0d0 .and. Flux3%LE /= error .and. Flux3%H /= error) then
        lEx%Bowen = Flux3%H / Flux3%LE
    else
        lEx%Bowen = error
    end if
end subroutine Fluxes23
