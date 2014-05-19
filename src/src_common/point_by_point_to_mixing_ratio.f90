!***************************************************************************
! point_by_point_to_mixing_ratio.f90
! ----------------------------------
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
! \brief       Convert mole fractions (moles per mole wet air) or molar density \n
!              to mixing ratios (moles per mole dry air)
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
    real(kind = dbl) :: TmpH2O(nrow)
    real(kind = dbl) :: Va(nrow)
    logical :: make_d2r


    !> Define a local variable for deciding if molar density can be transformed to mixing ratio
    !> and if so calculates time series of air molar volume [kg+1m-3]
    make_d2r = .false.
    if (E2Col(tc)%present .and. E2Col(pi)%present .and. E2Col(h2o)%instr%path_type == 'closed') then
        make_d2r = .true.
        where (Set(1:nrow, pi) /= 0d0)
            Va(1:nrow) = Ru * Set(1:nrow, tc) / Set(1:nrow, pi)
        elsewhere
            Va(1:nrow) = error
        end where
    end if

    !> Locally transform h2o into mole fraction [mmol/mol]
    select case (E2Col(h2o)%measure_type)
        case ('mixing_ratio')
            TmpH2O(:) = Set(:, h2o) / (1d0 + Set(:, h2o) / 1d3)
        case ('mole_fraction')
            TmpH2O(:) = Set(:, h2o)
        case ('molar_density')
            if (make_d2r) then
                where (Va(1:nrow) /= error)
                    TmpH2O(1:nrow) = Set(1:nrow, h2o) * Va(1:nrow)
                elsewhere
                    TmpH2O(1:nrow) = error
                end where
            else
                return
            end if
    end select

    !> If there is any scalar expressed as mole_fraction, coming from the same analyser as H2O
    !> convert it into mixing ratio using water vapor mole fraction as calculated above
    !> H2O mole fraction is expressed as  mmol_w / mol_a
    if (printout) write(*, '(a)', advance = 'no') '  Converting into mixing ratios..'
    if (E2Col(co2)%instr%model == E2Col(h2o)%instr%model) then
        if (E2Col(co2)%measure_type == 'mole_fraction')then
            E2Col(co2)%measure_type = 'mixing_ratio'
            where(TmpH2O(1:nrow) /= error .and. Set(1:nrow, co2) /= error)
                Set(:, co2) = Set(:, co2) / (1.d0 - TmpH2O(:) * 1d-3)
            elsewhere
                Set(1:nrow, co2) = error
            endwhere
        elseif (E2Col(co2)%measure_type == 'molar_density') then
            if (make_d2r) then
                E2Col(co2)%measure_type = 'mixing_ratio'
                where(Va(1:nrow) /= error .and. TmpH2O(1:nrow) /= error .and. Set(1:nrow, co2) /= error)
                    Set(1:nrow, co2) = Set(1:nrow, co2) * Va(1:nrow) * 1d3 / (1.d0 - TmpH2O(1:nrow) * 1d-3)
                elsewhere
                    Set(1:nrow, co2) = error
                endwhere
            else
                continue
            end if
        end if
    end if

    if (E2Col(h2o)%measure_type == 'mole_fraction') then
        E2Col(h2o)%measure_type = 'mixing_ratio'
        where(TmpH2O(1:nrow) /= error .and. Set(1:nrow, h2o) /= error)
        Set(:, h2o) = Set(:, h2o) / (1.d0 - TmpH2O(:) * 1d-3)
        elsewhere
            Set(1:nrow, h2o) = error
        endwhere
    elseif (E2Col(h2o)%measure_type == 'molar_density') then
        if (make_d2r) then
            E2Col(h2o)%measure_type = 'mixing_ratio'
            where(Va(1:nrow) /= error .and. TmpH2O(1:nrow) /= error)
                Set(1:nrow, h2o) = Set(1:nrow, h2o) * Va(1:nrow) / (1.d0 - TmpH2O(1:nrow) * 1d-3)
            elsewhere
                Set(1:nrow, h2o) = error
            endwhere
        else
            continue
        end if
    end if

    if (E2Col(ch4)%instr%model == E2Col(h2o)%instr%model) then
        if (E2Col(ch4)%measure_type == 'mole_fraction') then
            E2Col(ch4)%measure_type = 'mixing_ratio'
            where(TmpH2O(1:nrow) /= error .and. Set(1:nrow, ch4) /= error)
                Set(:, ch4) = Set(:, ch4) / (1.d0 - TmpH2O(:) * 1d-3)
            elsewhere
                Set(1:nrow, ch4) = error
            endwhere
        elseif (E2Col(ch4)%measure_type == 'molar_density') then
            if (make_d2r) then
                E2Col(ch4)%measure_type = 'mixing_ratio'
                where(Va(1:nrow) /= error .and. TmpH2O(1:nrow) /= error)
                    Set(1:nrow, ch4) = Set(1:nrow, ch4) * Va(1:nrow) * 1d3 / (1.d0 - TmpH2O(1:nrow) * 1d-3)
                elsewhere
                    Set(1:nrow, ch4) = error
                endwhere
            else
                continue
            end if
        end if
    end if

    if (E2Col(gas4)%instr%model == E2Col(h2o)%instr%model) then
        if (E2Col(gas4)%measure_type == 'mole_fraction') then
            E2Col(gas4)%measure_type = 'mixing_ratio'
            where(TmpH2O(1:nrow) /= error .and. Set(1:nrow, gas4) /= error)
                Set(:, gas4) = Set(:, gas4) / (1.d0 - TmpH2O(:) * 1d-3)
            elsewhere
                Set(1:nrow, gas4) = error
            endwhere
        elseif (E2Col(gas4)%measure_type == 'molar_density') then
            if (make_d2r) then
                E2Col(gas4)%measure_type = 'mixing_ratio'
                where(Va(1:nrow) /= error .and. TmpH2O(1:nrow) /= error)
                    Set(1:nrow, gas4) = Set(1:nrow, gas4) * Va(1:nrow) * 1d3 / (1.d0 - TmpH2O(1:nrow) * 1d-3)
                elsewhere
                    Set(1:nrow, gas4) = error
                endwhere
            else
                continue
            end if
        end if
    end if
    if (printout) write(*,'(a)') ' done.'
end subroutine PointByPointToMixingRatio
