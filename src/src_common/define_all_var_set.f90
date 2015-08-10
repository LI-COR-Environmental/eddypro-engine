!***************************************************************************
! define_all_var_set.f90
! ----------------------
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
! \brief       Creates a dataset with with all variables in phisical and \n
!              standardized units
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine DefineAllVarSet(LocCol, fRaw, nrow, ncol, N)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    integer, intent(in) :: ncol
    integer, intent(in) :: N
    type(ColType), intent(inout) :: LocCol(MaxNumCol)
    real(kind = sgl), intent(inout) :: fRaw(nrow, ncol)
    !> local variables
    integer :: j
    real(kind = sgl) :: DumVec(N)

    !> Converts units, according to information in the metadata file
    !> Physical, non-standard units are detected and converted into standard
    !> Volt (or any other "non-physical") are converted as well, using conversion parameters
    do j = 1, NumAllVar
        select case (LocCol(j)%var)
            !> Wind components, taken to [m s-1]
            case('u', 'v', 'w', 'sos')
                if (LocCol(j)%conversion_type == 'none') then
                    select case(LocCol(j)%unit_in(1:len_trim(LocCol(j)%unit_in)))
                        case ('m_sec')
                            if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) /= 'sos') cycle
                        case ('mm_sec')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 1e-3
                            end where
                            if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) /= 'sos') cycle
                        case ('cm_sec')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 1e-2
                            end where
                            if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) /= 'sos') cycle
                        case default
                            cycle
                    end select
                else
                    DumVec(1:N) = fRaw(1:N, j)
                    call LinearConversion(LocCol, DumVec(1:N), N, j)
                    fRaw(1:N, j) = DumVec(1:N)
                    select case(LocCol(j)%unit_out(1:len_trim(LocCol(j)%unit_out)))
                        case ('m_sec')
                            if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) /= 'sos') cycle
                        case ('mm_sec')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 1e-3
                            end where
                            if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) /= 'sos') cycle
                        case ('cm_sec')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 1e-2
                            end where
                            if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) /= 'sos') cycle
                        case default
                            if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) /= 'sos') cycle
                        end select
                end if
                !> if speed-of-sound, it is taken from [m s-1] to sonic temperature [K]
                if(LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'sos') then
                    where(fRaw(1:N, j) /= error)
                        fRaw(1:N, j) = (fRaw(1:N, j))**2 / 403.
                    end where
                    cycle
                end if

            !> Temperatures (K)
            case('ts', 'cell_t', 'int_t_1', 'int_t_2', 'air_t')
                if (LocCol(j)%conversion_type == 'none') then
                    select case(LocCol(j)%unit_in(1:len_trim(LocCol(j)%unit_in)))
                        case ('kelvin')
                            cycle
                        case ('ckelvin')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 1e-2
                            end where
                            cycle
                        case ('celsius')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) + 273.15
                            end where
                        case ('ccelsius')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 1e-2 + 273.15
                            end where
                            cycle
                        case default
                            cycle
                    end select
                else
                    DumVec(1:N) = fRaw(1:N, j)
                    call LinearConversion(LocCol, DumVec(1:N), N, j)
                    fRaw(1:N, j) = DumVec(1:N)
                    select case(LocCol(j)%unit_out(1:len_trim(LocCol(j)%unit_out)))
                        case ('kelvin')
                            cycle
                        case ('ckelvin')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 1e-2
                            end where
                            cycle
                        case ('celsius')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) + 273.15
                            end where
                        case ('ccelsius')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 1e-2 + 273.15
                            end where
                            cycle
                        case default
                            cycle
                    end select
                end if

            !> Pressures (Pa)
            case('int_p', 'air_p')
                if (LocCol(j)%conversion_type == 'none') then
                    select case(LocCol(j)%unit_in(1:len_trim(LocCol(j)%unit_in)))
                        case ('pa')
                            cycle
                        case ('hpa')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 1e2
                            end where
                            cycle
                        case ('kpa')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 1e3
                            end where
                            cycle
                    end select
                else
                    DumVec(1:N) = fRaw(1:N, j)
                    call LinearConversion(LocCol, DumVec(1:N), N, j)
                    fRaw(1:N, j) = DumVec(1:N)
                    select case(LocCol(j)%unit_out(1:len_trim(LocCol(j)%unit_out)))
                        case ('pa')
                            cycle
                        case ('hpa')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 1e2
                            end where
                            cycle
                        case ('kpa')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 1e3
                            end where
                            cycle
                        case default
                            cycle
                    end select
                end if

            !> Concentrations of trace gases
            case('co2', 'ch4', 'n2o')
                if (LocCol(j)%conversion_type == 'none') then
                    select case(LocCol(j)%unit_in(1:len_trim(LocCol(j)%unit_in)))
                        case ('mmol_m3', 'ppm')
                            cycle
                        case ('ppt')
                            fRaw(1:N, j) = fRaw(1:N, j) * 1e3
                        case ('umol_m3', 'ppb')
                            fRaw(1:N, j) = fRaw(1:N, j) * 1e-3
                        case ('g_m3')
                            select case (LocCol(j)%var)
                                case('co2')
                                    where(fRaw(1:N, j) /= error)
                                        fRaw(1:N, j) = fRaw(1:N, j) / MW(co2)
                                    end where
                                case('ch4')
                                    where(fRaw(1:N, j) /= error)
                                        fRaw(1:N, j) = fRaw(1:N, j) / MW(ch4)
                                    end where
                                case('n2o')
                                    where(fRaw(1:N, j) /= error)
                                        fRaw(1:N, j) = fRaw(1:N, j) / MW(gas4)
                                    end where
                            end select
                        case ('mg_m3')
                            select case (LocCol(j)%var)
                                case('co2')
                                    where(fRaw(1:N, j) /= error)
                                        fRaw(1:N, j) = fRaw(1:N, j) / MW(co2) * 1e-3
                                    end where
                                case('ch4')
                                    where(fRaw(1:N, j) /= error)
                                        fRaw(1:N, j) = fRaw(1:N, j) / MW(ch4) * 1e-3
                                    end where
                                case('n2o')
                                    where(fRaw(1:N, j) /= error)
                                        fRaw(1:N, j) = fRaw(1:N, j) / MW(gas4) * 1e-3
                                    end where
                            end select
                        case ('ug_m3')
                            select case (LocCol(j)%var)
                                case('co2')
                                    where(fRaw(1:N, j) /= error)
                                        fRaw(1:N, j) = fRaw(1:N, j) / MW(co2) * 1e-6
                                    end where
                                case('ch4')
                                    where(fRaw(1:N, j) /= error)
                                        fRaw(1:N, j) = fRaw(1:N, j) / MW(ch4) * 1e-6
                                    end where
                                case('n2o')
                                    where(fRaw(1:N, j) /= error)
                                        fRaw(1:N, j) = fRaw(1:N, j) / MW(gas4) * 1e-6
                                    end where
                            end select
                        case default
                            cycle
                    end select
                else
                    DumVec(1:N) = fRaw(1:N, j)
                    call LinearConversion(LocCol, DumVec(1:N), N, j)
                    fRaw(1:N, j) = DumVec(1:N)
                    select case(LocCol(j)%unit_out(1:len_trim(LocCol(j)%unit_out)))
                        case ('ppt')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 1e3
                            end where
                        case ('umol_m3', 'ppb')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 1e-3
                            end where
                        case ('g_m3')
                            select case (LocCol(j)%var)
                                case('co2')
                                    where(fRaw(1:N, j) /= error)
                                        fRaw(1:N, j) = fRaw(1:N, j) / MW(co2)
                                    end where
                                case('ch4')
                                    where(fRaw(1:N, j) /= error)
                                        fRaw(1:N, j) = fRaw(1:N, j) / MW(ch4)
                                    end where
                                case('n2o')
                                    where(fRaw(1:N, j) /= error)
                                        fRaw(1:N, j) = fRaw(1:N, j) / MW(gas4)
                                    end where
                            end select
                        case ('mg_m3')
                            select case (LocCol(j)%var)
                                case('co2')
                                    where(fRaw(1:N, j) /= error)
                                        fRaw(1:N, j) = fRaw(1:N, j) / MW(co2) * 1e-3
                                    end where
                                case('ch4')
                                    where(fRaw(1:N, j) /= error)
                                        fRaw(1:N, j) = fRaw(1:N, j) / MW(ch4) * 1e-3
                                    end where
                                case('n2o')
                                    where(fRaw(1:N, j) /= error)
                                        fRaw(1:N, j) = fRaw(1:N, j) / MW(gas4) * 1e-3
                                    end where
                            end select
                        case ('ug_m3')
                            select case (LocCol(j)%var)
                                case('co2')
                                    where(fRaw(1:N, j) /= error)
                                        fRaw(1:N, j) = fRaw(1:N, j) / MW(co2) * 1e-6
                                    end where
                                case('ch4')
                                    where(fRaw(1:N, j) /= error)
                                        fRaw(1:N, j) = fRaw(1:N, j) / MW(ch4) * 1e-6
                                    end where
                                case('n2o')
                                    where(fRaw(1:N, j) /= error)
                                        fRaw(1:N, j) = fRaw(1:N, j) / MW(gas4) * 1e-6
                                    end where
                            end select
                        case default
                            cycle
                    end select
                end if

            !> Concentrations of water vapour
            case('h2o')
                if (LocCol(j)%conversion_type == 'none') then
                    select case(LocCol(j)%unit_in(1:len_trim(LocCol(j)%unit_in)))
                        case ('mmol_m3', 'ppt')
                            cycle
                        case ('umol_m3', 'ppm')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 1e-3
                            end where
                        case ('ppb')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 1e-6
                            end where
                        case ('g_m3')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) / MW(h2o)
                            end where
                        case ('mg_m3')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) / MW(h2o) * 1e-3
                            end where
                        case ('ug_m3')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) / MW(h2o) * 1e-6
                            end where
                        case default
                            cycle
                    end select
                else
                    DumVec(1:N) = fRaw(1:N, j)
                    call LinearConversion(LocCol, DumVec(1:N), N, j)
                    fRaw(1:N, j) = DumVec(1:N)
                    select case(LocCol(j)%unit_out(1:len_trim(LocCol(j)%unit_out)))
                        case ('umol_m3', 'ppm')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 1e-3
                            end where
                        case ('ppb')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 1e-6
                            end where
                        case ('g_m3')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) / MW(h2o)
                            end where
                        case ('mg_m3')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) / MW(h2o) * 1e-3
                            end where
                        case ('ug_m3')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) / MW(h2o) * 1e-6
                            end where
                        case default
                            cycle
                    end select
                end if

            !> Flow rates (m+3s-1)
            case('flowrate')
                if (LocCol(j)%conversion_type == 'none') then
                    select case(LocCol(j)%unit_in(1:len_trim(LocCol(j)%unit_in)))
                        case ('m3_s')
                            cycle
                        case ('cm3_s')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 1e-6
                            end where
                            cycle
                        case ('lit_m')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) *  1.66666667e-5
                            end where
                            cycle
                        case ('ft3_s')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 0.028316846592e0
                            end where
                            cycle
                        case ('in3_s')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 1.6387064e-5
                            end where
                            cycle
                        case default
                            cycle
                    end select
                else
                    DumVec(1:N) = fRaw(1:N, j)
                    call LinearConversion(LocCol, DumVec(1:N), N, j)
                    fRaw(1:N, j) = DumVec(1:N)
                    select case(LocCol(j)%unit_out(1:len_trim(LocCol(j)%unit_out)))
                        case ('m3_s')
                            cycle
                        case ('cm3_s')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 1e-6
                            end where
                            cycle
                        case ('lit_m')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) *  1.66666667e-5
                            end where
                            cycle
                        case ('ft3_s')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 0.028316846592
                            end where
                            cycle
                        case ('in3_s')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 1.6387064e-5
                            end where
                            cycle
                    end select
                end if

            !> All other variables (custom variables)
            case default
                DumVec(1:N) = fRaw(1:N, j)
                call LinearConversion(LocCol, DumVec(1:N), N, j)
                fRaw(1:N, j) = DumVec(1:N)
        end select
    end do
end subroutine DefineAllVarSet

!***************************************************************************
!
! \brief       Performs a linear rescaling as either gain/offset or min/max.
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine LinearConversion(LocCol, Vec, nrow, j)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    integer, intent(in) :: j
    type(ColType), intent(in) :: LocCol(MaxNumCol)
    real(kind = sgl), intent(inout) :: Vec(nrow)

    select case (LocCol(j)%conversion_type(1:len_trim(LocCol(j)%conversion_type)))
        case ('zero_fullscale')
            where(Vec(:) /= error)
                Vec(:) = sngl(((LocCol(j)%b - LocCol(j)%a) / &
                    (LocCol(j)%max - LocCol(j)%min)) * (dble(Vec(:)) - LocCol(j)%min) + LocCol(j)%a)
            end where
        case ('gain_offset')
            where(Vec(:) /= error)
                Vec(:) = sngl(LocCol(j)%a * dble(Vec(:)) + LocCol(j)%b)
            end where
        case default
            continue
    end select
end subroutine LinearConversion
