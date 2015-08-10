!***************************************************************************
! molefractions_and_mixingratios.f90
! ----------------------------------
! Copyright (C) 2007-2011, Eco2s team, Gerardo Fratini
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
! \brief       Calculate average mole fractions, mixing ratios and
!              molar density
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo        change name (include densities)
!***************************************************************************
subroutine MoleFractionsAndMixingRatios()
    use m_rp_global_var
    implicit none
    !> local variables
    integer :: gas
    real(kind=dbl) :: LocVa(GHGNumVar)


    !> Initialization
    do gas = co2, gas4
        Stats%d(gas) = error
        Stats%chi(gas) = error
        Stats%r(gas) = error
        !> Set correct air molar volume, depending
        !> on instrument path type (closed or open)
        if (E2Col(gas)%present &
            .and. E2Col(gas)%instr%path_type == 'closed') then
            LocVa(gas) = E2Col(gas)%Va
        else
            LocVa(gas) = Ambient%Va
        end if
    end do

    !> First calculate stuff for H2O
    select case (E2Col(h2o)%measure_type)
        case ('mixing_ratio')
            !> If water vapour is mixing ratio, convert to mole fraction
            Stats%r(h2o)   = Stats%Mean(h2o)
            Stats%chi(h2o) = Stats%Mean(h2o) / (1d0 + Stats%Mean(h2o) / 1d3)
            if (LocVa(h2o) > 0d0) then
                Stats%d(h2o) = Stats%chi(h2o) / LocVa(h2o)
            else
                Stats%d(h2o) = error
            end if
        case('mole_fraction')
            !> If water vapour is already mole fraction, takes mean value
            Stats%chi(h2o) = Stats%Mean(h2o)
            Stats%r(h2o)   = Stats%Mean(h2o) / (1.d0 - Stats%Mean(h2o) / 1d3)
            if (LocVa(h2o) > 0d0) then
                Stats%d(h2o) = Stats%chi(h2o) / LocVa(h2o)
            else
                Stats%d(h2o) = error
            end if
        case('molar_density')
            !> If water vapour is molar density [mmol_w m-3]
            !> calculate mole fraction [mmol_w mol_a-1] by multiplication by
            !> air mole volume [m+3 mol_a-1]
            Stats%d(h2o) = Stats%Mean(h2o)
            if (LocVa(h2o) > 0d0) then
                Stats%chi(h2o) = Stats%Mean(h2o) * LocVa(h2o)
                Stats%r(h2o)   = Stats%chi(h2o) / (1.d0 - Stats%chi(h2o) / 1d3)
            else
                Stats%chi(h2o) = error
                Stats%r(h2o) = error
            end if
        case default
            Stats%d(h2o) = error
            Stats%r(h2o) = error
            Stats%chi(h2o) = error
    end select

    !> Calculate average mole fractions and mixing ratios where applicable
    do gas = co2, gas4
        if (gas /= h2o) then
            select case (E2Col(gas)%measure_type)
                case('mixing_ratio')
                    Stats%r(gas)   = Stats%Mean(gas)
                    if (Stats%r(h2o) /= error) then
                        Stats%chi(gas) = Stats%Mean(gas) &
                            / (1.d0 + Stats%r(h2o) * 1d-3)
                    else
                        Stats%chi(gas) = Stats%r(gas)
                    end if
                    if (LocVa(gas) > 0d0) then
                        Stats%d(gas) = Stats%chi(gas) / LocVa(gas) * 1d-3
                    else
                        Stats%d(gas) = error
                    end if
                case('mole_fraction')
                    Stats%chi(gas) = Stats%Mean(gas)
                    if (Stats%chi(h2o) /= error) then
                        Stats%r(gas) = Stats%chi(gas) &
                            / (1.d0 - Stats%chi(h2o) * 1d-3)
                    else
                        Stats%r(gas) = Stats%chi(gas)
                    end if
                    if (LocVa(gas) > 0d0) then
                        Stats%d(gas) = Stats%chi(gas) / LocVa(gas) * 1d-3
                    else
                        Stats%d(gas) = error
                    end if
                case('molar_density')
                    Stats%d(gas) = Stats%Mean(gas)
                    if (LocVa(gas) > 0d0) then
                        Stats%chi(gas) = Stats%Mean(gas) * LocVa(gas) * 1d3
                        if (Stats%chi(h2o) /= error) then
                            Stats%r(gas) = Stats%chi(gas) &
                                / (1.d0 - Stats%chi(h2o) * 1d-3)
                        else
                            Stats%r(gas) = Stats%chi(gas)
                        end if
                    else
                        Stats%chi(h2o) = error
                        Stats%r(h2o) = error
                    end if
            end select
        end if
    end do
end subroutine MoleFractionsAndMixingRatios
