!***************************************************************************
! flux_params.f90
! ---------------
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
! \brief       Calculate micromet and auxilary params useful for \n
!              flux computation and correction, and for user analysis
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine FluxParams(printout)
    use m_rp_global_var
    implicit none
    !> in/out variables
    logical, intent(in) :: printout
    !> local variables
    ! real(kind = dbl) :: Ma
    real(kind = dbl) :: Cpd
    real(kind = dbl) :: Cpv

    if (printout) write(*,'(a)', advance = 'no') &
        '  Calculating auxiliary variables..'

    Ambient%alpha = 0.51d0

    !> Water vapour partial pressure at saturation [Pa]
    !> (this formula gives same results as that in Buck (1981),
    !> cited in Campbell and Norman (1998) - Environmental Biophysics
    if (Stats%T > 0d0) then
        Ambient%es = (dexp(77.345d0 + 0.0057d0 * Stats%T &
                      - 7235.d0 / Stats%T)) / Stats%T**(8.2d0)
    else
        Ambient%es = error
    end if

    if (biomet%val(bRH) > 0d0 .and. biomet%val(bRH) < RHmax) then
        !> If meteo RH is available, uses it for all slow parameters,
        !> including redefining chi, r and d of H2O
        Stats%RH = biomet%val(bRH)
        !> Water vapour partial pressure [Pa]
        if (Ambient%es /= error) then
            Ambient%e = Stats%RH * 1d-2 * Ambient%es
        else
            Ambient%e = error
        end if
        !> vapor pressure deficit [Pa]
        if (Ambient%e /= error) then
            Ambient%VPD = Ambient%es - Ambient%e
        else
            Ambient%VPD = error
        end if
        !> Water vapour mass density [kg_w m-3]
        if (Ambient%e /= error .and. Stats%T /= error .and. Stats%T /= 0d0) then
            RHO%w = Ambient%e / (Rw * Stats%T)
        else
            RHO%w = error
        end if
        !> Water vapour concentrations and densities
        !> Water vapour mole fractions
        if (RHO%w /= error .and. Ambient%Va /= error) then
            Stats%chi(h2o) = RHO%w * Ambient%Va / MW(h2o) * 1d3
            !> Water vapour mixing ratio
            Stats%r(h2o)   = Stats%chi(h2o) / (1.d0 - Stats%chi(h2o) * 1d-3)
            !> Water vapour molar density
            if (E2Col(h2o)%instr%path_type == 'closed') then
                if (E2Col(h2o)%Va > 0d0) then
                    Stats%d(h2o) = Stats%chi(h2o) / E2Col(h2o)%Va
                else
                    Stats%d(h2o) = error
                    Stats%r(h2o) = error
                    Stats%chi(h2o) = error
                end if
            else
                Stats%d(h2o) = Stats%chi(h2o) / Ambient%Va
            end if
        else
            Stats%chi(h2o) = error
            Stats%r(h2o) = error
            Stats%d(h2o) = error
        end if
    else
        !> If meteo RH is not available or out of range, uses H2O from raw data
        !> Molecular weight of wet air:
        !> Ma = chi(h2o) * MW(h2o) + chi(dry_air) * Md
        !> if chi(dry_air) = 1 - chi(h2o) (assumes chi(h2o) in mmol mol_a-1)
        ! if (Stats%chi(h2o) > 0d0) then
        !     Ma = (Stats%chi(h2o) * 1d-3) * MW(h2o) &
        !        + (1d0 - Stats%chi(h2o) * 1d-3) * Md
        ! else
        !     Ma = error
        ! end if

        !> Water vapour mass density [kg_w m-3]
        !> from mole fraction [mmol_w / mol_a]
        !> (good also when native is molar density)
        if (Stats%chi(h2o) > 0d0 .and. Ambient%Va > 0d0) then
            RHO%w = (Stats%chi(h2o) / Ambient%Va) * MW(h2o) * 1d-3
        else
            RHO%w = error
        end if

        if (Stats%T > 0d0 .and. RHO%w >= 0d0) then
            !> Water vapour partial pressure [Pa]
            Ambient%e  =  RHO%w * Rw * Stats%T
            if (Ambient%es > 0d0) then
                !> Relative huimidity [%]
                Stats%RH = Ambient%e * 1d2 / Ambient%es
                !> vapor pressure deficit [hPa]
                Ambient%VPD = Ambient%es - Ambient%e
            else
                Stats%RH    = 0d0
                Ambient%VPD = 0d0
            end if
        else
            Ambient%e = error
            Ambient%es = error
            Stats%RH = error
            Ambient%VPD = error
        end if
        if (Stats%RH < 0d0 .or. Stats%RH > RHmax) then
            Stats%RH = error
            Ambient%VPD = error
        end if
        if (Stats%RH > 100d0 .and. Stats%RH < RHmax) then
            Stats%RH = 100d0 !< RH slightly higher than 100% is set to 100%
            Ambient%VPD = 0d0 !< RH slightly higher than 100%, VPD is set to 0
        end if
    end if

    !> Dew-point temperature [K], after Campbell and Norman (1998)
    !> Environmental Biophysics. Here e is in Pa, thus it must be divided
    !> by 10^3 to get kPa as in the formula.
    if (Ambient%e > 0d0) then
        Ambient%Td = (240.97d0 * dlog(Ambient%e * 1d-3 /0.611d0) &
            / (17.502d0 - dlog(Ambient%e * 1d-3 / 0.611d0))) + 273.15d0
    else
        Ambient%Td = error
    end if

    !> Dry air partial pressure [Pa], as:
    !> ambient P minus water vapor partial pressure
    if (Stats%Pr > 0d0) then
        if (Ambient%e > 0d0) then
            Ambient%p_d =  Stats%Pr - Ambient%e
        else
            Ambient%p_d = Stats%Pr
        end if
    else
        Ambient%p_d = error
    end if

    !> Molar volume of dry air [m+3 mol_d-1] after Ibrom et al. (2007, Tellus B)
    if (Ambient%p_d > 0d0) then
        Ambient%Vd = (Stats%Pr * Ambient%Va) / Ambient%p_d
    else
        Ambient%Vd = error
    end if

    !> Density of dry air [kg_d m-3]
    if (Stats%T > 0d0) then
        RHO%d = Ambient%p_d / (Rd * Stats%T)
    else
        RHO%d = error
    end if

    !> Dry air heat capacity at costant pressure [J+1kg-1K-1],
    !> as a function of temperature
    Cpd = 1005d0 + (Stats%T - 273.15d0 + 23.12d0)**2 / 3364d0

    !> Density of wet air [kg_a m-3]
    if (RHO%d > 0d0) then
        if (RHO%w >= 0d0) then
            RHO%a = RHO%d + RHO%w
            !> alternative: analytically derived from RHO%a = Pa * Ma / (Ru * T)
            !> Gives identical result.
            !RHO%a = (Stats%Pr - (1d0-MW(h2o)/Md) * Ambient%e) / (Rd*Stats%T)
        else
            RHO%a = RHO%d
        end if
    else
        RHO%a = error
    end if

    !> Specific humidity [kg_w kg_a-1]
    if (RHO%a > 0d0 .and. RHO%w >= 0d0) then
        Ambient%Q = RHO%w / RHO%a
    else
        Ambient%Q = error
    end if

    !> Air temperature = sonic temperature (corrected for side-wind) \n
    !> corrected for humidity (T in K), or = biomet T
    !> Condition is posed on either biomet T or raw air T (Mean(te))
    if (Stats%Mean(te) > 0d0 .or. biomet%val(bTa) > 0d0) then
        Ambient%Ta = Stats%T
        if (Stats%Mean(ts) > 0d0) then
            !> temperature mapping factor (Van Dijk et al. 2004, eq.3.1)
            Ambient%Tmap = Ambient%Ta / Stats%Mean(ts)
        else
            Ambient%Tmap = error
        end if
    elseif (E2Col(ts)%instr%category == 'fast_t_sensor') then
        !> If Ts was actually from a fast temperature sensor,
        !> do not apply Q correction
        Ambient%Ta = Stats%Mean(ts)
        Ambient%Tmap = 1d0
    else
        if (Ambient%Q > 0d0 .and. Ambient%alpha /= error &
            .and. Stats%Mean(ts) > 0d0) then
            Ambient%Ta = Stats%Mean(ts) / (1.d0 + Ambient%alpha * Ambient%Q)
            Ambient%Tmap = Ambient%Ta / Stats%Mean(ts)
        else
            Ambient%Ta = Stats%Mean(ts)
            Ambient%Tmap = 1d0
        end if

        !> Iterate the calculation of main quantities,
        !> after having better estimated air T
        if (Ambient%Ta > 0d0) then
            Ambient%es = (dexp(77.345d0 + 0.0057d0 * Ambient%Ta &
                          - 7235.d0 / Ambient%Ta)) / Ambient%Ta**(8.2d0)

            if (RHO%w >= 0d0) then
                Ambient%e  =  RHO%w * Rw * Ambient%Ta
                Stats%RH = Ambient%e * 1d2 / Ambient%es
                Ambient%VPD = Ambient%es - Ambient%e
                if (Stats%RH < 0d0 .or. Stats%RH > RHmax) then
                    Stats%RH = error
                    Ambient%VPD = error
                end if
                if (Stats%RH > 100d0 .and. Stats%RH < RHmax) then
                    Stats%RH = 100d0 !< RH slightly higher than 100% is set to 100%
                    Ambient%VPD = 0d0 !< RH slightly higher than 100% VPD is set to 0
                end if
                Ambient%Td = (240.97d0 * dlog(Ambient%e * 1d-3 /0.611d0) &
                    / (17.502d0 - dlog(Ambient%e * 1d-3 / 0.611d0))) + 273.15d0
                Ambient%p_d =  Stats%Pr - Ambient%e
                if (Ambient%p_d > 0d0) then
                    Ambient%Vd = (Stats%Pr * Ambient%Va) / Ambient%p_d
                else
                    Ambient%Vd = error
                end if
                RHO%d = Ambient%p_d / (Rd * Ambient%Ta)
                RHO%a = RHO%d + RHO%w
                Cpd = 1005d0 + (Ambient%Ta - 273.15d0 + 23.12d0)**2 / 3364d0
                Ambient%Q = RHO%w / RHO%a
                if (E2Col(ts)%instr%category == 'fast_t_sensor') then
                    !> If Ts was actually from a fast temperature sensor,
                    !> do not apply Q correction
                    Ambient%Ta = Stats%Mean(ts)
                    Ambient%Tmap = 1d0
                else
                    Ambient%Ta = Stats%Mean(ts) &
                        / (1.d0 + Ambient%alpha * Ambient%Q)
                    Ambient%Tmap = Ambient%Ta / Stats%Mean(ts)
                end if
            else
                Ambient%e = error
                Stats%RH = error
                Ambient%Q = error
                Ambient%Td = error
                Ambient%p_d = Stats%Pr
                Ambient%Vd = (Stats%Pr * Ambient%Va) / Ambient%p_d
                RHO%d = Ambient%p_d / (Rd * Ambient%Ta)
                RHO%a = RHO%d
                Cpd = 1005d0 + (Ambient%Ta - 273.15d0 + 23.12d0)**2 &
                    / 3364d0
                Ambient%Ta  = Stats%Mean(ts)
                Ambient%Tmap = 1
            end if
        end if
    end if

    !> Cell Temperature, if applicable
    if (Stats%Mean(tc) > 0d0) then
        Ambient%Tcell = Stats%Mean(tc)
    else
        Ambient%Tcell = Ambient%Ta
    end if

    !> Water vapour heat capacity at costant pressure [J+1kg-1K-1],
    !> as a function of temperature and RH
    Cpv = 1859d0 + 0.13d0 * Stats%RH &
        + (0.193d0 + 5.6d-3 * Stats%RH) * (Ambient%Ta - 273.15d0) &
        + (1d-3 + 5d-5 * Stats%RH) * (Ambient%Ta - 273.15d0)**2
    !> RhoAir by Cp (this is wet air Cp), in [J+1K-1m-3]
    if (RHO%d > 0d0 .and. RHO%w >= 0d0) then
            Ambient%RhoCp = Cpv * RHO%w + Cpd * Rho%d
            !> Alternative formulation (gives identical result)
            !Ambient%RhoCp = RHO%a * (Cpd * (1d0 - Ambient%Q) + Cpv * Ambient%Q)
        elseif (RHO%d > 0d0) then
            !> If RHO%d exists but RHO%w not, RhoCp is
            !> calculated as referred to dry air
            Ambient%RhoCp = Cpd * Rho%d
        else
            Ambient%RhoCp = error
    end if

    !> Specific heat of evaporation [J kg-1 K-1]
    if (Ambient%Ta > 0d0) then
        !> Gives same result of:
        !> lambda = − 0.0000614342*T^3 + 0.00158927*T^2 − 2.36418*T
        !> + 2500.79 in a large range (-30 to 50 °C)
        Ambient%lambda = (3147.5d0 - 2.37d0 * Ambient%Ta) * 1d3
    else
        Ambient%lambda = error
    end if

    !> water to dry air density ratio [adim.]
    if (RHO%d >0d0 .and. RHO%w > 0d0) then
        Ambient%sigma = RHO%w / RHO%d
    else
        Ambient%sigma = error
    end if

    if (printout) write(*,'(a)') ' Done.'
end subroutine FluxParams

