!***************************************************************************
! read_timelag_opt_file.f90
! -------------------------
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
! \brief       Read time-lag optimization file and import relevant parameters
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ReadTimelagOptFile(ncls)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer :: ncls
    !> local variables
    integer :: open_status
    integer :: read_status
    character(500) :: strg


    !> Open planar fit file and read rotation matrices
    write(*,'(a)') ' Reading time-lag optimization file: ' // AuxFile%to(1:len_trim(AuxFile%to))
    open(udf, file = AuxFile%to, status = 'old', iostat = open_status)

    if (open_status == 0) then
        write(*, '(a)') '  Time lag optimization file found, retrieving content..'
        do
            !> co2
            read(udf, '(a)', iostat = read_status) strg
            if (read_status /= 0) exit
            if(index(strg, 'Median_co2_timelag_[s]') /= 0) then
                read(strg(index(strg, ':')+1:len_trim(strg)), '(f6.2)') toPasGas(co2)%def
                cycle
            end if
            if(index(strg, 'Mimimum_co2_timelag_[s]') /= 0) then
                read(strg(index(strg, ':')+1:len_trim(strg)), '(f6.2)') toPasGas(co2)%min
                cycle
            end if
            if(index(strg, 'Maximum_co2_timelag_[s]') /= 0) then
                read(strg(index(strg, ':')+1:len_trim(strg)), '(f6.2)') toPasGas(co2)%max
                cycle
            end if
            !> h2o
            if(index(strg, 'Median_h2o_timelag_[s]') /= 0) then
                ncls = 0
                read(strg(index(strg, ':')+1:len_trim(strg)), '(f6.2)') toPasGas(h2o)%def
                cycle
            end if
            if(index(strg, 'Mimimum_h2o_timelag_[s]') /= 0) then
                read(strg(index(strg, ':')+1:len_trim(strg)), '(f6.2)') toPasGas(h2o)%min
                cycle
            end if
            if(index(strg, 'Maximum_h2o_timelag_[s]') /= 0) then
                read(strg(index(strg, ':')+1:len_trim(strg)), '(f6.2)') toPasGas(h2o)%max
                cycle
            end if
            !> ch4
            if(index(strg, 'Median_ch4_timelag_[s]') /= 0) then
                read(strg(index(strg, ':')+1:len_trim(strg)), '(f6.2)') toPasGas(ch4)%def
                cycle
            end if
            if(index(strg, 'Mimimum_ch4_timelag_[s]') /= 0) then
                read(strg(index(strg, ':')+1:len_trim(strg)), '(f6.2)') toPasGas(ch4)%min
                cycle
            end if
            if(index(strg, 'Maximum_ch4_timelag_[s]') /= 0) then
                read(strg(index(strg, ':')+1:len_trim(strg)), '(f6.2)') toPasGas(ch4)%max
                cycle
            end if
            !> 4th gas
            if(index(strg, 'Median_4th_gas_timelag_[s]') /= 0) then
                read(strg(index(strg, ':')+1:len_trim(strg)), '(f6.2)') toPasGas(gas4)%def
                cycle
            end if
            if(index(strg, 'Mimimum_4th_gas_timelag_[s]') /= 0) then
                read(strg(index(strg, ':')+1:len_trim(strg)), '(f6.2)') toPasGas(gas4)%min
                cycle
            end if
            if(index(strg, 'Maximum_4th_gas_timelag_[s]') /= 0) then
                read(strg(index(strg, ':')+1:len_trim(strg)), '(f6.2)') toPasGas(gas4)%max
                cycle
            end if

            !> h2o as a function of RH
            ncls = 0
            if (index(strg, 'H2O_timelag_determinations_as_a_function') /= 0) then
                !> Skip one line
                read(udf, *)
                read(udf, *)
                !> Read as many classes as available
                do
                    read(udf, '(a)', iostat = read_status) strg
                    if (read_status /= 0 .or. index(strg, '%') == 0) exit
                    ncls = ncls + 1
                    read(strg(20:len_trim(strg)), '(3(f14.2))') &
                        toH2O(ncls)%def,  toH2O(ncls)%min,  toH2O(ncls)%max
                end do
                exit
            end if
        end do
    else
       !> If the specified file is not found or is empty, switches to covariance maximization without default
        Meth%tlag = 'maxcov'
        call ExceptionHandler(39)
    end if
    write(*,'(a)')   '  Done.'
end subroutine ReadTimelagOptFile
