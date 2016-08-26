!***************************************************************************
! writeout_timelag_optimization.f90
! ---------------------------------
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
! \brief       Write time-lag optimization results on output file \n
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WriteOutTimelagOptimization(actn, M, h2o_n, ncls, cls_size)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: M
    integer, intent(in) :: ncls
    integer, intent(in) :: actn(M)
    real(kind = dbl), intent(in) :: cls_size
    integer, intent(out) :: h2o_n(ncls)
    integer, external :: CreateDir
    !> local variables
    integer :: cls
    integer :: open_status = 1
    character(4) :: min
    character(4) :: max
    character(9) :: txt


    !> Create output file
    TimelagOpt_Path = Dir%main_out(1:len_trim(Dir%main_out)) &
              // EddyProProj%id(1:len_trim(EddyProProj%id)) &
              // TimelagOpt_FilePadding // Timestamp_FilePadding // TxtExt
    open(uto, file = TimelagOpt_Path, iostat = open_status, encoding = 'utf-8')

    !> Write on output file time-lag optimization results
    write(uto, '(a)') 'Time-lag_optimisation_results'
    write(uto, '(a, f7.2)') 'Plausibility_range_[timefolds_standard_deviation]:',TOSetup%pg_range
    write(uto, '(a, a)') 'Beginning_of_timelag_optimization_period: ', TOSetup%start_date
    write(uto, '(a, a)') 'End_of_timelag_optimization_period: ', TOSetup%end_date
    write(uto, '(a)')

    if (E2Col(co2)%present) then
        write(uto, '(a, i5)') 'Number_of_timelags_used_for_co2:', actn(co2)
        write(uto, '(a, f6.2)') 'Median_co2_timelag_[s]:', toPasGas(co2)%def
        write(uto, '(a, f6.2)') 'Mimimum_co2_timelag_[s]:', toPasGas(co2)%min
        write(uto, '(a, f6.2)') 'Maximum_co2_timelag_[s]:', toPasGas(co2)%max
        write(uto, '(a)')
    end if

    if (E2Col(h2o)%present .and. ncls <= 1) then
        write(uto, '(a, i5)') 'Number_of_timelags_used_for_h2o:', actn(h2o)
        write(uto, '(a, f6.2)') 'Median_h2o_timelag_[s]:', toPasGas(h2o)%def
        write(uto, '(a, f6.2)') 'Mimimum_h2o_timelag_[s]:', toPasGas(h2o)%min
        write(uto, '(a, f6.2)') 'Maximum_h2o_timelag_[s]:', toPasGas(h2o)%max
        write(uto, '(a)')
    end if

    if (E2Col(ch4)%present) then
        write(uto, '(a, i4)') 'Number_of_timelags_used_for_ch4:', actn(ch4)
        write(uto, '(a, f6.2)') 'Median_ch4_timelag_[s]:', toPasGas(ch4)%def
        write(uto, '(a, f6.2)') 'Mimimum_ch4_timelag_[s]:', toPasGas(ch4)%min
        write(uto, '(a, f6.2)') 'Maximum_ch4_timelag_[s]:', toPasGas(ch4)%max
        write(uto, '(a)')
    end if

    if (E2Col(gas4)%present) then
        write(uto, '(a, i4)') 'Number_of_timelags_used_for_4th_gas:', actn(gas4)
        write(uto, '(a, f6.2)') 'Median_4th_gas_timelag_[s]:' , toPasGas(gas4)%def
        write(uto, '(a, f6.2)') 'Mimimum_4th_gas_timelag_[s]:', toPasGas(gas4)%min
        write(uto, '(a, f6.2)') 'Maximum_4th_gas_timelag_[s]:', toPasGas(gas4)%max
        write(uto, '(a)')
    end if

    if (E2Col(h2o)%present .and. ncls > 1) then
        write(uto, '(a, i4)') 'H2O_timelag_determinations_as_a_function_of_relative_humidity'
        write(uto, '(a, i4)') 'Classes with numerosity < 30 are inferred (see software documentation)'
        write(uto,'(a)')             'class     RH-range       med_h2o       min_h2o       max_h2o     class_num'
        do cls = 1, ncls
            write(min, '(i4)') nint((cls - 1) * cls_size)
            call ShrinkString(min)
            write(max, '(i4)') nint(cls * cls_size)
            call ShrinkString(max)
            txt = min(1:len_trim(min)) // ' - ' // max(1:len_trim(max)) // '%'
            write(uto,'(i5, 5x, a9, 3(f13.2,1x), i13)') cls,  txt, toH2O(cls)%def, toH2O(cls)%min, toH2O(cls)%max, h2o_n(cls)
        end do
    end if
    close(uto)
    write(*,'(a)') '  Results written on file: ' &
        // TimelagOpt_Path(1:len_trim(TimelagOpt_Path))
end subroutine WriteOutTimelagOptimization

