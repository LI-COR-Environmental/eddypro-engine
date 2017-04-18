!***************************************************************************
! point_by_point_to_mixing_ratio.f90
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
! \brief       Convert mole fractions (moles per mole wet air) or \n
!              molar density to mixing ratios (moles per mole dry air)
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine PointByPointToMixingRatio(Set, nrow, ncol, printout)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    logical, intent(in) :: printout
    real(kind = dbl), intent(inout) :: Set(nrow, ncol)
    !> Local variables
    real(kind = dbl) :: H2Omf(nrow)
    real(kind = dbl) :: Va(nrow)
    logical :: cellVaAvailable


    !> Point-by-point, accurate conversion to mixing ratio cannot be
    !> performed if H2O is not measured by a closed-path analyzer
    if (.not. E2Col(h2o)%present &
        .or. E2Col(h2o)%instr%path_type /= 'closed') return

    !> Calculates time series of air molar volume [kg+1m-3] in the cell of the
    !> of the instrument for which T and P are available.
    if (E2Col(tc)%present .and. E2Col(pi)%present &
        .and. E2Col(tc)%instr%model == E2Col(h2o)%instr%model &
        .and. E2Col(pi)%instr%model == E2Col(h2o)%instr%model) then
        where (Set(:, pi) > 0d0 .and. Set(:, tc) > 0d0)
            Va(:) = Ru * Set(:, tc) / Set(:, pi)
        elsewhere
            Va(:) = error
        end where
    else
        Va(:) = error
    end if

    if (all(Va(:) == error)) then
        cellVaAvailable = .false.
    else
        cellVaAvailable = .true.
    end if
    
    !> Locally transform h2o into mole fraction [mmol/mol]
    select case (E2Col(h2o)%measure_type)
        case ('mixing_ratio')
            where(Set(:, h2o) /=  error)
                H2Omf(:) = Set(:, h2o) / (1d0 + Set(:, h2o) / 1d3)
            elsewhere
                H2Omf(:) = error
            end where
        case ('mole_fraction')
            H2Omf(:) = Set(:, h2o)
        case ('molar_density')
            where (Va(:) /= error .and. Set(:, h2o) /= error)
                H2Omf(:) = Set(:, h2o) * Va(:)
            elsewhere
                H2Omf(:) = error
            end where
    end select

    !> If there is any scalar expressed as mole_fraction or molar density, 
    !> measured by the same analyzer of H2O, convert it into mixing ratio using
    !> water vapor mole fraction as calculated above
    !> (H2O mole fraction is expressed as  mmol_w / mol_a)
    !> and using cell air molar volume
    if (printout) write(*, '(a)', advance = 'no') &
        '  WPL step: converting into mixing ratios wherever possible..'

    !> Special case of H2O
    if (E2Col(h2o)%measure_type == 'mole_fraction') then
        E2Col(h2o)%measure_type = 'mixing_ratio'
        where(H2Omf(:) /= error .and. Set(:, h2o) /= error)
            Set(:, h2o) = Set(:, h2o) / (1.d0 - H2Omf(:) * 1d-3)
        elsewhere
            Set(:, h2o) = error
        endwhere
    elseif (E2Col(h2o)%measure_type == 'molar_density' .and. cellVaAvailable) then
        E2Col(h2o)%measure_type = 'mixing_ratio'
        where(Va(:) /= error .and. H2Omf(:) /= error)
            Set(:, h2o) = Set(:, h2o) * Va(:) / (1.d0 - H2Omf(:) * 1d-3)
        elsewhere
            Set(:, h2o) = error
        endwhere
    end if

    !> CO2
    if (E2Col(co2)%instr%model == E2Col(h2o)%instr%model) then
        if (E2Col(co2)%measure_type == 'mole_fraction')then
            E2Col(co2)%measure_type = 'mixing_ratio'
            where(H2Omf(:) /= error .and. Set(:, co2) /= error)
                Set(:, co2) = Set(:, co2) / (1.d0 - H2Omf(:) * 1d-3)
            elsewhere
                Set(:, co2) = error
            endwhere
        elseif (E2Col(co2)%measure_type == 'molar_density' .and. cellVaAvailable) then
            E2Col(co2)%measure_type = 'mixing_ratio'
            where(Va(:) /= error .and. H2Omf(:) /= error &
                .and. Set(:, co2) /= error)
                Set(:, co2) = Set(:, co2) * Va(:) * 1d3 &
                    / (1.d0 - H2Omf(:) * 1d-3)
            elsewhere
                Set(:, co2) = error
            endwhere
        end if
    end if

    !> CH4
    if (E2Col(ch4)%instr%model == E2Col(h2o)%instr%model) then
        if (E2Col(ch4)%measure_type == 'mole_fraction') then
            E2Col(ch4)%measure_type = 'mixing_ratio'
            where(H2Omf(:) /= error .and. Set(:, ch4) /= error)
                Set(:, ch4) = Set(:, ch4) / (1.d0 - H2Omf(:) * 1d-3)
            elsewhere
                Set(:, ch4) = error
            endwhere
        elseif (E2Col(ch4)%measure_type == 'molar_density' .and. cellVaAvailable) then
            E2Col(ch4)%measure_type = 'mixing_ratio'
            where(Va(:) /= error .and. H2Omf(:) /= error &
                .and. Set(:, ch4) /= error)
                Set(:, ch4) = Set(:, ch4) * Va(:) * 1d3 &
                    / (1.d0 - H2Omf(:) * 1d-3)
            elsewhere
                Set(:, ch4) = error
            endwhere
        end if
    end if

    if (E2Col(gas4)%instr%model == E2Col(h2o)%instr%model) then
        if (E2Col(gas4)%measure_type == 'mole_fraction') then
            E2Col(gas4)%measure_type = 'mixing_ratio'
            where(H2Omf(:) /= error .and. Set(:, gas4) /= error)
                Set(:, gas4) = Set(:, gas4) / (1.d0 - H2Omf(:) * 1d-3)
            elsewhere
                Set(:, gas4) = error
            endwhere
        elseif (E2Col(gas4)%measure_type == 'molar_density' .and. cellVaAvailable) then
            E2Col(gas4)%measure_type = 'mixing_ratio'
            where(Va(:) /= error .and. H2Omf(:) /= error &
            .and. Set(:, gas4) /= error)
                Set(:, gas4) = Set(:, gas4) * Va(:) * 1d3 &
                    / (1.d0 - H2Omf(:) * 1d-3)
            elsewhere
                Set(:, gas4) = error
            endwhere
        end if
    end if
    if (printout) write(*,'(a)') ' Done.'
end subroutine PointByPointToMixingRatio
