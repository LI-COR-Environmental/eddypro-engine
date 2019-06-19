!***************************************************************************
! filter_dataset_for_physical_thresholds.f90
! ------------------------------------------
! Copyright (C) 2007-2011, Eco2s team, Gerardo Fratini
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
! \brief       Apply "absolute limits" filter to variables selected in SelectWhat
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine FilterDatasetForPhysicalThresholds(Set, N, M, FilterWhat)
    use m_rp_global_var
    !> in/out variables
    integer, intent(in) :: N
    integer, intent(in) :: M
    logical, intent(in) :: FilterWhat(M)
    real(kind = dbl), intent(inout) :: Set(N, M)


    !> Eliminate for values out of bounds in the whole dataset
    !> co2, expressed as [mmol m-3] if molar_density, [umol mol-1] otherwise
    if (E2Col(co2)%present .and. FilterWhat(co2)) then
        if(E2Col(co2)%measure_type == 'molar_density') then
            do i = 1, N
                if( (Set(i, co2) * Ambient%Va * 1d3 < al%co2_min) .or. &
                    (Set(i, co2) * Ambient%Va * 1d3 > al%co2_max)) Set(i, co2) = error
            end do
        else
            do i = 1, N
                if( (Set(i, co2) < al%co2_min) .or. &
                (Set(i, co2) > al%co2_max)) Set(i, co2) = error
            end do
        end if
    end if

    !> h2o, expressed as [mmol m-3] if molar_density, [mmol mol-1] otherwise
    if (E2Col(h2o)%present .and. FilterWhat(h2o)) then
        if(E2Col(h2o)%measure_type == 'molar_density') then
            do i = 1, N
                if( (Set(i, h2o) * Ambient%Va < al%h2o_min) .or. &
                (Set(i, h2o) * Ambient%Va > al%h2o_max)) Set(i, h2o) = error
            end do
        else
            do i = 1, N
                if( (Set(i, h2o) < al%h2o_min) .or. &
                (Set(i, h2o) > al%h2o_max)) Set(i, h2o) = error
            end do
        end if
    end if

    !> ch4, expressed as [mmol m-3] if molar_density, [umol mol-1] otherwise
    if (E2Col(ch4)%present .and. FilterWhat(ch4)) then
        if(E2Col(ch4)%measure_type == 'molar_density') then
            do i = 1, N
                if( (Set(i, ch4) * Ambient%Va * 1d3 < al%ch4_min) .or. &
                (Set(i, ch4) * Ambient%Va * 1d3 > al%ch4_max)) Set(i, ch4) = error
            end do
        else
            do i = 1, N
                if( (Set(i, ch4) < al%ch4_min) .or. &
                (Set(i, ch4) > al%ch4_max)) Set(i, ch4) = error
            end do
        end if
    end if

    !> gas4, expressed as [mmol m-3] if molar_density, [umol mol-1] otherwise
    if (E2Col(gas4)%present .and. FilterWhat(gas4)) then
        if(E2Col(gas4)%measure_type == 'molar_density') then
            do i = 1, N
                if( (Set(i, gas4) * Ambient%Va * 1d3 < al%gas4_min) .or. &
                (Set(i, gas4) * Ambient%Va * 1d3 > al%gas4_max)) Set(i, gas4) = error
            end do
        else
            do i = 1, N
                if( (Set(i, gas4) < al%gas4_min) .or. &
                (Set(i, gas4) > al%gas4_max)) Set(i, gas4) = error
            end do
        end if
    end if
end subroutine FilterDatasetForPhysicalThresholds
