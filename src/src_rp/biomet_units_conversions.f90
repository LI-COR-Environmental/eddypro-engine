!***************************************************************************
! biomet_units_conversions.f90
! ----------------------------
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


    do i = 1, nbVars
        select case(trim(bVars(i)%nature))
            case('TEMPERATURE')
                select case(bVars(i)%unit_in)
                    case('C','�C')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) + 273.15d0
                        end where
                    case('F','�F')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = (bSet(:, i) - 32d0) * 5d0 / 9d0 &
                                + 273.15d0
                        end where
                    case('CK')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 1d-2
                        end where
                    case('CC','C�C')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 1d-2 + 273.15d0
                        end where
                    case('CF','C�F')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = (bSet(:, i) * 1d-2 - 32d0) * 5d0 / 9d0 &
                                + 273.15d0
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
    bAggrFluxnet = bAggrEddyPro

    !> Change units as needed
    do i = 1, nbVars
        !> All temperatures converted to [degC]
        if (trim(bVars(i)%nature) == 'TEMPERATURE' .and. bAggrEddyPro(i) /= error) &
            bAggrFluxnet(i) = bAggrEddyPro(i) - 273.15d0
        !> Air pressure is converted to [kPa]
        if (bVars(i)%fluxnet_base_name == 'PA' .and. bAggrEddyPro(i) /= error) &
            bAggrFluxnet(i) = bAggrEddyPro(i) * 1d-3
        !> VPD is converted to [hPa]
        if (bVars(i)%fluxnet_base_name == 'VPD' .and. bAggrEddyPro(i) /= error) &
            bAggrFluxnet(i) = bAggrEddyPro(i) * 1d-2
        !> All precipitations are converted to [mm]
        if (trim(bVars(i)%nature) == 'PRECIPITATION' .and. bAggrEddyPro(i) /= error) &
            bAggrFluxnet(i) = bAggrEddyPro(i) * 1d3
        !> Snow depth is converted to [cm]
        if (bVars(i)%fluxnet_base_name == 'SNOW_D' .and. bAggrEddyPro(i) /= error) &
            bAggrFluxnet(i) = bAggrEddyPro(i) * 1d2
        !> Water table depth is converted to [cm]
        if (bVars(i)%fluxnet_base_name == 'WATER_TABLE_DEPTH' .and. bAggrEddyPro(i) /= error) &
            bAggrFluxnet(i) = bAggrEddyPro(i) * 1d2
        !> SWC is converted to [%]
        if (bVars(i)%fluxnet_base_name == 'SWC' .and. bAggrEddyPro(i) /= error) &
            bAggrFluxnet(i) = bAggrEddyPro(i) * 1d2
        !> RUNOFF is converted to [mm]
        if (bVars(i)%fluxnet_base_name == 'RUNOFF' .and. bAggrEddyPro(i) /= error) &
            bAggrFluxnet(i) = bAggrEddyPro(i) * 1d3
        !> THROUGHFALL is converted to [mm]
        if (bVars(i)%fluxnet_base_name == 'THROUGHFALL' .and. bAggrEddyPro(i) /= error) &
            bAggrFluxnet(i) = bAggrEddyPro(i) * 1d3
        !> DBH is converted to [cm]
        if (bVars(i)%fluxnet_base_name == 'DBH' .and. bAggrEddyPro(i) /= error) &
            bAggrFluxnet(i) = bAggrEddyPro(i) * 1d2
    end do
end subroutine BiometStandardFluxnetUnits
