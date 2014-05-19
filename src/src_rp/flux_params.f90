!***************************************************************************
! flux_params.f90
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
    real(kind = dbl) :: Ma
    real(kind = dbl) :: Cpd
    real(kind = dbl) :: Cpv

    if (printout) write(*,'(a)', advance = 'no') '  Calculating auxiliary variables..'

    LitePar%alpha = 0.51d0

    !> Water vapour partial pressure at saturation [Pa] (this formula gives same results as
    !> that in Buck (1981), cited in Campbell and Norman (1998) - Environmental Biophysics
    if (Stats%T > 0d0) then
        LitePar%es = (dexp(77.345d0 + 0.0057d0 * Stats%T &
                      - 7235.d0 / Stats%T)) / Stats%T**(8.2d0)
    else
        LitePar%es = error
    end if

    if (Stats%mRH > 0d0 .and. Stats%mRH < RHmax) then
        !> If meteo RH is available, uses it for all slow parameters, including redefining
        !> chi, r and d of H2O
        Stats%RH = Stats%mRH
        !> Water vapour partial pressure [Pa]
        if (LitePar%es /= error) then
            LitePar%e = Stats%RH * 1d-2 * LitePar%es
        else
            LitePar%e = error
        end if
        !> vapor pressure deficit [hPa]
        if (LitePar%e /= error) then
            LitePar%VPD = LitePar%es - LitePar%e
        else
            LitePar%VPD = error
        end if
        !> Water vapour mass density [kg_w m-3]
        if (LitePar%e /= error .and. Stats%T /= error .and. Stats%T /= 0d0) then
            RHO%w = LitePar%e / (Rw * Stats%T)
        else
            RHO%w = error
        end if
        !> Water vapour concentrations and densities
        !> Water vapour mole fractions
        if (RHO%w /= error .and. LitePar%Va /= error) then
            Stats%chi(h2o) = RHO%w * LitePar%Va / MW(h2o) * 1d3
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
                Stats%d(h2o) = Stats%chi(h2o) / LitePar%Va
            end if
        else
            Stats%chi(h2o) = error
            Stats%r(h2o) = error
            Stats%d(h2o) = error
        end if
    else
        !> If meteo RH is not available or out of range, uses H2O from raw data
        !> Molecular weight of wet air --> Ma = chi(h2o) * MW(h2o) + chi(dry_air) * Md
        !> if chi(dry_air) = 1 - chi(h2o) (assumes chi(h2o) in mmol mol_a-1)
        if (Stats%chi(h2o) > 0d0) then
            Ma = (Stats%chi(h2o) * 1d-3) * MW(h2o) &
               + (1d0 - Stats%chi(h2o) * 1d-3) * Md
        else
            Ma = error
        end if

        !> Water vapour mass density [kg_w m-3]
        !> from mole fraction [mmol_w / mol_a] (good also when native is molar density)
        if (Stats%chi(h2o) > 0d0 .and. LitePar%Va > 0d0) then
            RHO%w = (Stats%chi(h2o) / LitePar%Va) * MW(h2o) * 1d-3
        else
            RHO%w = error
        end if

        if (Stats%T > 0d0 .and. RHO%w >= 0d0) then
            !> Water vapour partial pressure [Pa]
            LitePar%e  =  RHO%w * Rw * Stats%T
            if (LitePar%es > 0d0) then
                !> Relative huimidity [%]
                Stats%RH = LitePar%e * 1d2 / LitePar%es
                !> vapor pressure deficit [hPa]
                LitePar%VPD = LitePar%es - LitePar%e
            else
                Stats%RH    = 0d0
                LitePar%VPD = 0d0
            end if
        else
            LitePar%e = error
            LitePar%es = error
            Stats%RH = error
            LitePar%VPD = error
        end if
        if (Stats%RH < 0d0 .or. Stats%RH > RHmax) then
            Stats%RH = error
            LitePar%VPD = error
        end if
        if (Stats%RH > 100d0 .and. Stats%RH < RHmax) then
            Stats%RH = 100d0 !< RH slightly higher than 100% is set to 100%
            LitePar%VPD = 0d0 !< RH slightly higher than 100%, VPD is set to 0
        end if
    end if

    !> Dew-point temperature [K], after Campbell and Norman (1998) - Environmental Biophysics
    !> here e is in Pa, thus it must be divided by 10^3 to get kPa as in the formula.
    if (LitePar%e > 0d0) then
        LitePar%Td = (240.97d0 * dlog(LitePar%e * 1d-3 /0.611d0) / (17.502d0 - dlog(LitePar%e * 1d-3 / 0.611d0))) + 273.16d0
    else
        LitePar%Td = error
    end if

    !> Dry air partial pressure [Pa], as ambient P - water vapour partial pressure
    if (Stats%Pr > 0d0) then
        if (LitePar%e > 0d0) then
            LitePar%p_d =  Stats%Pr - LitePar%e
        else
            LitePar%p_d = Stats%Pr
        end if
    else
        LitePar%p_d = error
    end if

    !> Molar volume of dry air [m+3 mol_d-1] after Ibrom et al. (2007, Tellus B)
    if (LitePar%p_d > 0d0) then
        LitePar%Vd = (Stats%Pr * LitePar%Va) / LitePar%p_d
    else
        LitePar%Vd = error
    end if

    !> Density of dry air [kg_d m-3]
    if (Stats%T > 0d0) then
        RHO%d = LitePar%p_d / (Rd * Stats%T)
    else
        RHO%d = error
    end if

    !> Dry air heat capacity at costant pressure [J+1kg-1K-1], as a function of temperature
    Cpd = 1005d0 + (Stats%T - 273.16d0 + 23.12d0)**2 / 3364d0

    !> Density of wet air [kg_a m-3]
    if (RHO%d > 0d0) then
        if (RHO%w >= 0d0) then
            RHO%a = RHO%d + RHO%w
            !> alternative: analytically derived from RHO%a = Pa * Ma / (Ru * T), gives identical result
            !RHO%a = (Stats%Pr - (1d0 - MW(h2o) / Md) * LitePar%e) / (Rd * Stats%T)
        else
            RHO%a = RHO%d
        end if
    else
        RHO%a = error
    end if

    !> Specific humidity [kg_w kg_a-1]
    if (RHO%a > 0d0 .and. RHO%w >= 0d0) then
        LitePar%Q = RHO%w / RHO%a
    else
        LitePar%Q = error
    end if

    !> Air temperature = sonic temperature (corrected for side-wind) \n
    !> corrected for humidity (T in K), or = meteo T
    !> Condition is posed on either external meteo T, (mT) or raw air T (Mean(te))
    if (Stats%Mean(te) > 0d0 .or. Stats%mT > 0d0) then
        LitePar%Ta = Stats%T
        if (Stats%Mean(ts) > 0d0) then
            !> temperature mapping factor (Van Dijk et al. 2004, eq.3.1)
            LitePar%Tmap = LitePar%Ta / Stats%Mean(ts)
        else
            LitePar%Tmap = error
        end if
    elseif (E2Col(ts)%instr%category == 'fast_t_sensor') then
            !> If Ts was actually from a fast temperature sensor, do not apply Q correction
            LitePar%Ta = Stats%Mean(ts)
            LitePar%Tmap = 1d0
    else
        if (LitePar%Q > 0d0 .and. LitePar%alpha /= error .and. Stats%Mean(ts) > 0d0) then
            LitePar%Ta = Stats%Mean(ts) / (1.d0 + LitePar%alpha * LitePar%Q)
            LitePar%Tmap = LitePar%Ta / Stats%Mean(ts)
        else
            LitePar%Ta = Stats%Mean(ts)
            LitePar%Tmap = 1d0
        end if

        !> Iterate the calculation of main quantities, after having better estimated air T
        if (LitePar%Ta > 0d0) then
            LitePar%es = (dexp(77.345d0 + 0.0057d0 * LitePar%Ta &
                          - 7235.d0 / LitePar%Ta)) / LitePar%Ta**(8.2d0)

            if (RHO%w >= 0d0) then
                LitePar%e  =  RHO%w * Rw * LitePar%Ta
                Stats%RH = LitePar%e * 1d2 / LitePar%es
                LitePar%VPD = LitePar%es - LitePar%e
                if (Stats%RH < 0d0 .or. Stats%RH > RHmax) then
                    Stats%RH = error
                    LitePar%VPD = error
                end if
                if (Stats%RH > 100d0 .and. Stats%RH < RHmax) then
                    Stats%RH = 100d0 !< RH slightly higher than 100% is set to 100%
                    LitePar%VPD = 0d0 !< RH slightly higher than 100%, VPD is set to 0
                end if
                LitePar%Td = (240.97d0 * dlog(LitePar%e * 1d-3 /0.611d0) / (17.502d0 - dlog(LitePar%e * 1d-3 / 0.611d0))) + 273.16d0
                LitePar%p_d =  Stats%Pr - LitePar%e
                if (LitePar%p_d > 0d0) then
                    LitePar%Vd = (Stats%Pr * LitePar%Va) / LitePar%p_d
                else
                    LitePar%Vd = error
                end if
                RHO%d = LitePar%p_d / (Rd * LitePar%Ta)
                RHO%a = RHO%d + RHO%w
                Cpd = 1005d0 + (LitePar%Ta - 273.16d0 + 23.12d0)**2 / 3364d0
                LitePar%Q = RHO%w / RHO%a
                if (E2Col(ts)%instr%category == 'fast_t_sensor') then
                    !> If Ts was actually from a fast temperature sensor, do not apply Q correction
                    LitePar%Ta = Stats%Mean(ts)
                    LitePar%Tmap = 1d0
                else
                    LitePar%Ta = Stats%Mean(ts) / (1.d0 + LitePar%alpha * LitePar%Q)
                    LitePar%Tmap = LitePar%Ta / Stats%Mean(ts)
                end if
            else
                LitePar%e   = error
                Stats%RH    = error
                LitePar%Q   = error
                LitePar%Td  = error
                LitePar%p_d = Stats%Pr
                LitePar%Vd  = (Stats%Pr * LitePar%Va) / LitePar%p_d
                RHO%d       = LitePar%p_d / (Rd * LitePar%Ta)
                RHO%a       = RHO%d
                Cpd         = 1005d0 + (LitePar%Ta - 273.16d0 + 23.12d0)**2 / 3364d0
                LitePar%Ta  = Stats%Mean(ts)
                LitePar%Tmap = 1
            end if
        end if
    end if

    !> Cell Temperature, if applicable
    if (Stats%Mean(tc) > 0d0) then
        LitePar%Tcell = Stats%Mean(tc)
    else
        LitePar%Tcell = LitePar%Ta
    end if

    !> Water vapour heat capacity at costant pressure [J+1kg-1K-1], as a function of temperature and RH
    Cpv = 1859d0 + 0.13d0 * Stats%RH + (0.193d0 + 5.6d-3 * Stats%RH) * (LitePar%Ta - 273.16d0) + &
          (1d-3 + 5d-5 * Stats%RH) * (LitePar%Ta - 273.16d0)**2
    !> RhoAir by Cp (this is wet air Cp), in [J+1K-1m-3]
    if (RHO%d > 0d0 .and. RHO%w >= 0d0) then
            LitePar%RhoCp = Cpv * RHO%w + Cpd * Rho%d
            !> Alternative formulation (gives identical result)
            !LitePar%RhoCp = RHO%a * (Cpd * (1d0 - LitePar%Q) + Cpv * LitePar%Q)
        elseif (RHO%d > 0d0) then
            !> If RHO%d exists but RHO%w not, RhoCp is calculated as referred to dry air
            LitePar%RhoCp = Cpd * Rho%d
        else
            LitePar%RhoCp = error
    end if

    !> Specific heat of evaporation [J kg-1 K-1]
    if (LitePar%Ta > 0d0) then
        !> give same result of lambda = − 0.0000614342*T^3 + 0.00158927*T^2 − 2.36418*T + 2500.79 in a large range (-30 to 50 °C)
        LitePar%lambda = (3147.5d0 - 2.37d0 * LitePar%Ta) * 1d3
    else
        LitePar%lambda = error
    end if

    !> water to dry air density ratio [adim.]
    if (RHO%d >0d0 .and. RHO%w > 0d0) then
        LitePar%sigma = RHO%w / RHO%d
    else
        LitePar%sigma = error
    end if

    if (printout) write(*,'(a)') ' done.'
end subroutine FluxParams

