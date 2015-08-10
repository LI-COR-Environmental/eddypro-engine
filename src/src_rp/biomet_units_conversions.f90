!***************************************************************************
! biomet_units_conversions.f90
! ----------------------------
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
! \brief       Convert input units into standard units
! \author      Gerardo Fratini
! \note
!              Radiations (Rg, Rn, Rd, Rr, LWin, LWout, Ruva, Ruvb) \n
!              are not expected to need unit conversion
!              Photons flux densities (PPFD, PPFDd, PPFDr, PPFDbc, APAR) \n
!              are not expected to need unit conversion
!              Alb, PRI, SWC, SHF are not expected to need unit conversion
! \sa
! \bug
! \deprecated
! \test
!***************************************************************************
subroutine BiometStandardEddyProUnits()
    use m_rp_global_var
    implicit none
    !> local variables
    integer :: i

    !> Temperatures
    do i = 1, nbVars
        select case(trim(bVars(i)%nature))
            case('TEMPERATURE')
                select case(bVars(i)%unit_in)
                    case('C','°C')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) + 273.16d0
                        end where
                    case('F','°F')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = (bSet(:, i) - 32d0) * 5d0 / 9d0 &
                                + 273.16d0
                        end where
                    case('CK')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 1d-2
                        end where
                    case('CC','C°C')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 1d-2 + 273.16d0
                        end where
                    case('CF','C°F')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = (bSet(:, i) * 1d-2 - 32d0) * 5d0 / 9d0 &
                                + 273.16d0
                        end where
                    case default
                end select

            case('RELATIVE_HUMIDITY')
                select case(bVars(i)%unit_in)
                    case('NUMBER','#','DIMENSIONLESS')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 1d2
                        end where
                    case default
                        continue
                end select

            case('PRESSURE')
                select case(bVars(i)%unit_in)
                    case('HPA')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 1d2
                        end where
                    case('KPA')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 1d3
                        end where
                    case('MMHG', 'TORR')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 133.32d0
                        end where
                    case('PSI')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 6894.6d0
                        end where
                    case('BAR')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 1d5
                        end where
                    case('ATM')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 0.980665d5
                        end where
                    case default
                        continue
                end select

!            !> Precipitation is converted to [m]
!            case('PRECIPITATION')
!                select case(bVars(i)%unit_in)
!                    case('NM')
!                        where (bSet(:, i) /= error)
!                            bSet(:, i) = bSet(:, i) * 1d-6
!                        end where
!                    case('UM')
!                        where (bSet(:, i) /= error)
!                            bSet(:, i) = bSet(:, i) * 1d-3
!                        end where
!                    case('CM')
!                        where (bSet(:, i) /= error)
!                            bSet(:, i) = bSet(:, i) * 10
!                        end where
!                    case('M')
!                        where (bSet(:, i) /= error)
!                            bSet(:, i) = bSet(:, i) * 1d3
!                        end where
!                    case default
!                        continue
!                end select

            !> Lengths
            !> converted to [m]
            case('LENGTH', 'PRECIPITATION')
                select case(bVars(i)%unit_in)
                    case('NM')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 1d-9
                        end where
                    case('UM')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 1d-6
                        end where
                    case('MM')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 1d-3
                        end where
                    case('CM')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 1d-2
                        end where
                    case('KM')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 1d3
                        end where
                    case default
                        continue
                end select

            case('CONCENTRATION')
                select case(trim(adjustl(bVars(i)%label)))
                !> CO2 is converted to PPM
                case ('CO2')
                    select case(bVars(i)%unit_in)
                        case('PPT', 'MMOL/MOL', 'MMOL+1MOL-1', 'MMOLMOL-1')
                            where (bSet(:, i) /= error)
                                bSet(:, i) = bSet(:, i) * 1d3
                            end where
                        case('PPB', 'NMOL/MOL', 'NMOL+1MOL-1', 'NMOLMOL-1')
                            where (bSet(:, i) /= error)
                                bSet(:, i) = bSet(:, i) * 1d-3
                            end where
                        case default
                            continue
                    end select
                !> H2O is converted to PPT
                case ('H2O')
                    select case(bVars(i)%unit_in)
                        case('PPM', 'UMOL/MOL', 'UMOL+1MOL-1', 'UMOLMOL-1')
                            where (bSet(:, i) /= error)
                                bSet(:, i) = bSet(:, i) * 1d-3
                            end where
                        case('PPB', 'NMOL/MOL', 'NMOL+1MOL-1', 'NMOLMOL-1')
                            where (bSet(:, i) /= error)
                                bSet(:, i) = bSet(:, i) * 1d-6
                            end where
                        case default
                            continue
                    end select
                !> Any other gas is converted to PPB
                case default
                    select case(bVars(i)%unit_in)
                        case('PPT', 'MMOL/MOL', 'MMOL+1MOL-1', 'MMOLMOL-1')
                            where (bSet(:, i) /= error)
                                bSet(:, i) = bSet(:, i) * 1d6
                            end where
                        case('PPM', 'UMOL/MOL', 'UMOL+1MOL-1', 'UMOLMOL-1')
                            where (bSet(:, i) /= error)
                                bSet(:, i) = bSet(:, i) * 1d3
                            end where
                        case default
                            continue
                    end select
                end select

            case('SPEED')
                select case(bVars(i)%unit_in)
                    case('CM+1S-1','CM/S','CMS^-1','CMS-1')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 1d-2
                        end where
                    case('MM+1S-1','MM/S','MMS^-1','MMS-1')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 1d-3
                        end where
                    case default
                        continue
                end select

            case('ANGULAR_DIRECTION')
            case('FLOW')
        end select
    end do
end subroutine BiometStandardEddyProUnits

!***************************************************************************
!
! \brief       Start from EddyPro-standardized units to create
!              dataset with FLUXNET-standardized units
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
!***************************************************************************
subroutine BiometStandardFluxnetUnits()
    use m_rp_global_var
    implicit none
    !> local variables
    integer :: i


    !> Most variables will have same units..
    bAggrFluxnet = bAggr

    !> Change units as needed
    do i = 1, nbVars
        !> All temperatures converted to [degC]
        if (trim(bVars(i)%nature) == 'TEMPERATURE' .and. bAggr(i) /= error) &
            bAggrFluxnet(i) = bAggr(i) - 273.16d0
        !> Air pressure is converted to [kPa]
        if (bVars(i)%fluxnet_base_name == 'PA' .and. bAggr(i) /= error) &
            bAggrFluxnet(i) = bAggr(i) * 1d-3
        !> VPD is converted to [hPa]
        if (bVars(i)%fluxnet_base_name == 'VPD' .and. bAggr(i) /= error) &
            bAggrFluxnet(i) = bAggr(i) * 1d-2
        !> All precipitations are converted to [mm]
        if (trim(bVars(i)%nature) == 'PRECIPITATION' .and. bAggr(i) /= error) &
            bAggrFluxnet(i) = bAggr(i) * 1d3
        !> Snow depth is converted to [cm]
        if (bVars(i)%fluxnet_base_name == 'SNOW_D' .and. bAggr(i) /= error) &
            bAggrFluxnet(i) = bAggr(i) * 1d2
        !> Water table depth is converted to [cm]
        if (bVars(i)%fluxnet_base_name == 'WATER_TABLE_DEPTH' .and. bAggr(i) /= error) &
            bAggrFluxnet(i) = bAggr(i) * 1d2
        !> SWC is converted to [%]
        if (bVars(i)%fluxnet_base_name == 'SWC' .and. bAggr(i) /= error) &
            bAggrFluxnet(i) = bAggr(i) * 1d2
        !> RUNOFF is converted to [mm]
        if (bVars(i)%fluxnet_base_name == 'RUNOFF' .and. bAggr(i) /= error) &
            bAggrFluxnet(i) = bAggr(i) * 1d3
        !> THROUGHFALL is converted to [mm]
        if (bVars(i)%fluxnet_base_name == 'THROUGHFALL' .and. bAggr(i) /= error) &
            bAggrFluxnet(i) = bAggr(i) * 1d3
        !> DBH is converted to [cm]
        if (bVars(i)%fluxnet_base_name == 'DBH' .and. bAggr(i) /= error) &
            bAggrFluxnet(i) = bAggr(i) * 1d2
    end do
end subroutine BiometStandardFluxnetUnits
