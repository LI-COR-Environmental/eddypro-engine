!***************************************************************************
! set_timelags.f90
! ----------------
! Copyright (C) 2011-2015, LI-COR Biosciences
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
! \brief       Sets nominal timelag and timelag ranges, based on instruments \n
!              and EC setup (tube, flow rate, sensor separation), but only if
!              values in the metadata file are inplausible.
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine SetTimelags()
    use m_rp_global_var
    implicit none
    !> local variables
    integer :: gas
    integer :: cls
    real(kind = dbl) :: lRH
    real(kind = dbl) :: mult(GHGNumVar)
    real(kind = dbl) :: tube_time(GHGNumVar)
    real(kind = dbl) :: tube_volume(GHGNumVar)
    real(kind = dbl) :: cell_time(GHGNumVar)
    real(kind = dbl) :: cell_volume(GHGNumVar)
    real(kind = dbl) :: safety

    !> Multiplier
    mult(:) = 2d0 !< For passive gases
    mult(h2o) = 10d0 !< For active gases
    safety = 0.3d0  !< Safety margin for min/max setting, should nominal tlag be very close to zero

    !> set time-lags to optimized values if selected so by user
    if (meth%tlag == 'tlag_opt') then
        do gas = co2, gas4
            if (E2Col(gas)%present) then
                if (gas /= h2o) then
                    !> Passive gases
                    E2Col(gas)%def_tl = toPasGas(gas)%def
                    E2Col(gas)%min_tl = toPasGas(gas)%min
                    E2Col(gas)%max_tl = toPasGas(gas)%max
                else
                    !> For water vapor, if requested adjust time-lag to current RH
                    !> either taken from meteo or estimated locally from raw data
                    if (TOSetup%h2o_nclass > 1) then
                        if (biomet%val(bRH) > 0d0 .and. biomet%val(bRH) < RHmax) then
                            !> If meteo RH is available, uses that one
                            lRH = biomet%val(bRH)
                        else
                            !> If meteo RH is not available, calculate one
                            call LocalRhEstimate(lRH)
                        end if
                        do cls = 1, TOSetup%h2o_nclass
                            if (lRH >= (cls - 1) * TOSetup%h2o_class_size .and. lRH <= cls * TOSetup%h2o_class_size) then
                                E2Col(gas)%def_tl = toH2O(cls)%def
                                E2Col(gas)%min_tl = toH2O(cls)%min
                                E2Col(gas)%max_tl = toH2O(cls)%max
                            end if
                        end do
                    end if
                end if
            end if
        end do
    else
        do gas = co2, gas4
            if (E2Col(gas)%instr%path_type == 'closed') then
                if (E2Col(gas)%def_tl == 0d0) then
                    tube_volume(gas) = (p * (E2Col(gas)%instr%tube_d / 2d0)**2 * E2Col(gas)%instr%tube_l)
                    tube_time(gas) = tube_volume(gas) / E2Col(gas)%instr%tube_f
                    cell_volume(gas) = (p * (E2Col(gas)%instr%hpath_length / 2d0)**2 * E2Col(gas)%instr%vpath_length)
                    cell_time(gas) = cell_volume(gas) / E2Col(gas)%instr%tube_f
                    E2Col(gas)%def_tl = tube_time(gas) + cell_time(gas)
                end if
                if (E2Col(gas)%min_tl == 0d0) E2Col(gas)%min_tl = &
                                              max(0d0, E2Col(gas)%def_tl - 2d0)
                if (E2Col(gas)%max_tl == 0d0) E2Col(gas)%max_tl = &
                                              E2Col(gas)%def_tl + mult(gas) * E2Col(gas)%def_tl + safety

            elseif (E2Col(gas)%instr%path_type == 'open') then
                if (E2Col(gas)%min_tl == 0d0) &
                    E2Col(gas)%min_tl = - dsqrt(E2Col(gas)%instr%hsep**2 + E2Col(gas)%instr%vsep**2) * 2d0 - safety
                if (E2Col(gas)%max_tl == 0d0) &
                    E2Col(gas)%max_tl = dsqrt(E2Col(gas)%instr%hsep**2 + E2Col(gas)%instr%vsep**2) * 2d0 + safety
            end if
        end do
    end if
end subroutine SetTimelags

!***************************************************************************
!
! \brief       Provides an RH estimate for the purpose of time-lag optimization \n
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine LocalRhEstimate(lRH)
    use m_rp_global_var
    implicit none
    !> in/out variables
    real(kind = dbl), intent(out) :: lRH
    !> local variables
    real(kind = dbl) :: lChi
    real(kind = dbl) :: locT
    real(kind = dbl) :: locPr
    real(kind = dbl) :: locVa
    real(kind = dbl) :: lRHOw
    real(kind = dbl) :: locES
    real(kind = dbl) :: Ma

    !> Air temperature and pressure estimates
    !> Last true condition determines which temperature is used
    locT = Stats%Mean(ts)
    if(Stats%Mean(te)  > 220d0 .and. Stats%Mean(te) < 340d0) locT = Stats%Mean(te)
    if(biomet%val(bTa) > 220d0 .and. biomet%val(bTa) < 340d0) locT = biomet%val(bTa)

    !> Last true condition determines which pressure is used
    locPr = Metadata%bar_press
    if(Stats%Mean(pe)  > 40000 .and. Stats%Mean(pe)  < 110000) locPr = Stats%Mean(pe)
    if(biomet%val(bPa) > 40000 .and. biomet%val(bPa) < 110000) locPr = biomet%val(bPa)

    !> Ambient air molar volume [m+3 mol-1] and air mass density [kg m-3]
    if (locPr > 0d0 .and. locT /= error) then
        locVa = Ru * locT / locPr
    else
        locVa = error
    end if

    !> First calculate stuff for H2O
    lChi = error
    select case (E2Col(h2o)%measure_type)
        case ('mixing_ratio')
            !> If water vapour is mixing ratio, convert to mole fraction
            lChi = Stats%Mean(h2o) / (1d0 + Stats%Mean(h2o) / 1d3)
        case('mole_fraction')
            lChi = Stats%Mean(h2o)
        case('molar_density')
            !> If water vapour is molar density [mmol_w m-3] calculate mole fraction
            !> [mmol_w mol_a-1] by multiplication by air mole volume [m+3 mol_a-1]
            if (locVa > 0d0) then
                lChi = Stats%Mean(h2o) * locVa
            else
                lChi = error
            end if
        case default
            lChi = error
    end select

    !> If meteo RH is not available or out of range, uses H2O from raw data
    !> Molecular weight of wet air --> Ma = chi(h2o) * MW(h2o) + chi(dry_air) * Md
    !> if chi(dry_air) = 1 - chi(h2o) (assumes chi(h2o) in mmol mol_a-1)
    if (lChi > 0d0) then
        Ma = (lChi * 1d-3) * MW(h2o) + (1d0 - lChi * 1d-3) * Md
    else
        Ma = error
    end if

    !> Water vapour mass density [kg_w m-3]
    !> from mole fraction [mmol_w / mol_a] (good also when native is molar density)
    if (lChi > 0d0 .and. locVa > 0d0) then
        lRHOw = (lChi / locVa) * MW(h2o) * 1d-3
    else
        lRHOw = error
    end if

    !> Water vapour partial pressure at saturation [Pa] (this formula gives same results as
    !> that in Buck (1981), cited in Campbell and Norman (1998) - Environmental Biophysics
    if (locT > 0d0) then
        locES = (dexp(77.345d0 + 0.0057d0 * locT - 7235.d0 / locT)) / locT**(8.2d0)
    else
        locES = error
    end if

    !> Relative humidity
    lRH = error
    if (locT > 0d0 .and. lRHOw >= 0d0 .and. locES > 0) then
        !> Relative huimidity [%]
        lRH = (lRHOw * Rw * locT) * 1d2 / locES
    end if
end subroutine LocalRhEstimate

