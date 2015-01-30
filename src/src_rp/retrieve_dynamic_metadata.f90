!***************************************************************************
! retrieve_dynamic_metadata.f90
! -----------------------------
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
! \brief       Reads external dynamic (time-varying) metadata file. If dynamic
!              metadata are found, they are used instead of those in
!              alternative or embedded metadata files
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine RetrieveDynamicMetadata(FinalTimestamp, LocCol, ncol)
    use m_rp_global_var
    implicit none
    !> In/out variables
    integer, intent(in) :: ncol
    type (DateType), intent(in) :: FinalTimestamp
    type (ColType), intent(out) :: LocCol(ncol)
    !> local variables
    integer :: cnt
    integer :: open_status
    integer :: read_status
    integer :: var_num
    integer :: sepa
    character(10) :: date
    character(5) :: time
    character(32) :: mdStringVars(256)
    character(32) :: mdCurrentStringVars(256)
    character(LongInstringLen) :: dataline
    type (DateType) :: mdCurrentTimestamp


    write(*,'(a)', advance = 'no') '  Retrieving dynamic metadata..'

    !> Open dynamic metadata file
    open(udf, file = AuxFile%DynMD, status = 'old', iostat = open_status)
    if (open_status /= 0) then
        write(*,*)
        call ExceptionHandler(68)
        return
    end if

    !> Skip header
    read(udf, '(a)') dataline

    !> Start cycle on remaining records to find an appropriate one
    mdCurrentStringVars = 'none'
    cnt = 0
    record_loop: do
        dataline = ''
        var_num = 0
        !> Read date and time if present and exit on error
        read(udf, '(a)', iostat = read_status) dataline
        if (read_status /= 0) then
            if (cnt == 0) exit record_loop
        else
            if (len(trim(dataline)) == 0) cycle record_loop

            !> Store data in a temporary array as text
            mdStringVars = 'none'
            do
                sepa = index(dataline, separator)
                if (sepa == 0) sepa = len_trim(dataline) + 1
                if (len_trim(dataline) == 0) exit
                var_num = var_num + 1
                mdStringVars(var_num) = dataline(1:sepa - 1)
                dataline = dataline(sepa + 1: len_trim(dataline))
            end do

            !> Retrieve date and time and check suitability to current time period
            if (DynamicMetadataOrder(dynmd_date) /= nint(error)) &
                read(mdStringVars(DynamicMetadataOrder(dynmd_date)), *) date
            if (DynamicMetadataOrder(dynmd_time) /= nint(error)) &
                read(mdStringVars(DynamicMetadataOrder(dynmd_time)), *) time

            !> Retrieve timestamp from date/time
            call DateTimeToDateType(date, time, mdCurrentTimestamp)

            !> Check suitability of Metadata for current period
            !> Normal case
            if (mdCurrentTimestamp < FinalTimestamp) then
                cnt = cnt + 1
                mdCurrentStringVars = mdStringVars
                cycle record_loop
            end if
            if (cnt == 0) exit record_loop
        end if

        !> If it gets here, means that it's time to
        !> retrieve Dynamic Metadata from previous dataline
        call ReadMetadataFromTextVars(mdCurrentStringVars, size(mdStringVars))

        !> Normal loop exit instruction (if CurrentTimestamp
        !> is beyond FinalTimestamp)
        if (mdCurrentTimestamp >= FinalTimestamp) exit record_loop

        !> If it arrived here, means that Dynamic Metadata
        !> were found and retrieved. Can exit cycle
        exit record_loop
    end do record_loop
    close(udf)

    !> Correct dynamic metadata in case of plain mistakes
    call FixDynamicMetadata()

    !> Now update Metadata and E2Col with current DynamicMetadata,
    !> whether or not they have been updated at the current round
    call ExtractUsableMetadataFromDynamic(LocCol, size(LocCol))

    write(*, '(a)') ' Done.'
end subroutine RetrieveDynamicMetadata

!*******************************************************************************
!
! \brief       Determine whether metadata are relevant to current period
!              Metadata are relevant if:
!              - if they are "in the past", if they are
!                closer to InitialTimestamp than the last used
!              - if they are "in the future", if they overlap with current
!                period at any extent
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!*******************************************************************************
!logical function MetadataAreRelevantToCurrentPeriod(InitialTimestamp, FinalTimestamp, &
!    mdCurrentTimestamp, LastMetadataTimestamp)
!    use m_rp_global_var
!    implicit none
!    > In/out variables
!    type (DateType), intent(in) :: InitialTimestamp
!    type (DateType), intent(in) :: FinalTimestamp
!    type (DateType), intent(in) :: mdCurrentTimestamp
!    type (DateType), intent(in) :: LastMetadataTimestamp
!
!
!    MetadataAreRelevantToCurrentPeriod = .false.
!    if (mdCurrentTimestamp < InitialTimestamp) then
!        !> If CurrentMetada is earlier than InitialTimestamp,
!        !> check if its lag is smaller than
!        !> that of LastMetadata used
!        if (abs(Timelag(InitialTimestamp, mdCurrentTimestamp)) &
!            <= abs(TimeLag(InitialTimestamp, LastMetadataTimestamp))) &
!            MetadataAreRelevantToCurrentPeriod = .true.
!    else
!        !> If CurrentMetada is later or equal to InitialTimestamp,
!        !> check if it overlaps with current period
!        if (mdCurrentTimestamp < FinalTimestamp) &
!           MetadataAreRelevantToCurrentPeriod = .true.
!    end if
!
!end function MetadataAreRelevantToCurrentPeriod

!***************************************************************************
!
! \brief       Extract values of variables stored as text
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ReadMetadataFromTextVars(mdStringVars, nrow)
    use m_rp_global_var
    implicit none
    !> In/out variables
    integer, intent(in) :: nrow
    character(*), intent(in) :: mdStringVars(nrow)

    !> Site location
    if (DynamicMetadataOrder(altitude) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(altitude)), *) &
        DynamicMetadata%alt
    if (DynamicMetadataOrder(latitude) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(latitude)), *) &
        DynamicMetadata%lat
    if (DynamicMetadataOrder(longitude) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(longitude)), *) &
        DynamicMetadata%lon

    !> File length and acquisition frequency
    if (DynamicMetadataOrder(file_length) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(file_length)), *) &
        DynamicMetadata%file_length
    if (DynamicMetadataOrder(acquisition_frequency) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(acquisition_frequency)), *) &
        DynamicMetadata%ac_freq

    !> Canopy height, displacement height and roughness length
    if (DynamicMetadataOrder(canopy_height) /= nint(error)) then
        read(mdStringVars(DynamicMetadataOrder(canopy_height)), *) &
        DynamicMetadata%canopy_height
        DynamicMetadata%d = 0.66d0 * DynamicMetadata%canopy_height
        DynamicMetadata%z0 = 0.2d0 * DynamicMetadata%canopy_height
    end if
    if (DynamicMetadataOrder(displacement_height) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(displacement_height)), *) &
        DynamicMetadata%d
    if (DynamicMetadataOrder(roughness_length) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(roughness_length)), *) &
        DynamicMetadata%z0

    !> Master sonic info
    if (DynamicMetadataOrder(master_sonic_manufacturer) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(master_sonic_manufacturer)), *) &
        DynamicMetadata%instr(u)%firm
    if (DynamicMetadataOrder(master_sonic_model) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(master_sonic_model)), *) &
        DynamicMetadata%instr(u)%model
    if (DynamicMetadataOrder(master_sonic_height) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(master_sonic_height)), *) &
        DynamicMetadata%instr(u)%height
    if (DynamicMetadataOrder(master_sonic_wformat) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(master_sonic_wformat)), *) &
        DynamicMetadata%instr(u)%wformat
    if (DynamicMetadataOrder(master_sonic_wref) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(master_sonic_wref)), *) &
        DynamicMetadata%instr(u)%wref
    if (DynamicMetadataOrder(master_sonic_north_offset) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(master_sonic_north_offset)), *) &
        DynamicMetadata%instr(u)%north_offset
    if (DynamicMetadataOrder(master_sonic_hpath_length) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(master_sonic_hpath_length)), *) &
        DynamicMetadata%instr(u)%hpath_length
    if (DynamicMetadataOrder(master_sonic_vpath_length) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(master_sonic_vpath_length)), *) &
        DynamicMetadata%instr(u)%vpath_length
    if (DynamicMetadataOrder(master_sonic_tau) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(master_sonic_tau)), *) &
        DynamicMetadata%instr(u)%tau

    !> co2 irga
    if (DynamicMetadataOrder(co2_irga_manufacturer) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(co2_irga_manufacturer)), *) &
        DynamicMetadata%instr(co2)%firm
    if (DynamicMetadataOrder(co2_irga_model) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(co2_irga_model)), *) &
        DynamicMetadata%instr(co2)%model
    if (DynamicMetadataOrder(co2_measure_type) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(co2_measure_type)), *) &
        DynamicMetadata%measure_type(co2)
    if (DynamicMetadataOrder(co2_irga_northward_separation) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(co2_irga_northward_separation)), *) &
        DynamicMetadata%instr(co2)%nsep
    if (DynamicMetadataOrder(co2_irga_eastward_separation) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(co2_irga_eastward_separation)), *) &
        DynamicMetadata%instr(co2)%esep
    if (DynamicMetadataOrder(co2_irga_vertical_separation) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(co2_irga_vertical_separation)), *) &
        DynamicMetadata%instr(co2)%vsep
    if (DynamicMetadataOrder(co2_irga_tube_length) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(co2_irga_tube_length)), *) &
        DynamicMetadata%instr(co2)%tube_l
    if (DynamicMetadataOrder(co2_irga_tube_diameter) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(co2_irga_tube_diameter)), *) &
        DynamicMetadata%instr(co2)%tube_d
    if (DynamicMetadataOrder(co2_irga_tube_flowrate) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(co2_irga_tube_flowrate)), *) &
        DynamicMetadata%instr(co2)%tube_f
    if (DynamicMetadataOrder(co2_irga_kw) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(co2_irga_kw)), *) &
        DynamicMetadata%instr(co2)%kw
    if (DynamicMetadataOrder(co2_irga_ko) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(co2_irga_ko)), *) &
        DynamicMetadata%instr(co2)%ko
    if (DynamicMetadataOrder(co2_irga_hpath_length) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(co2_irga_hpath_length)), *) &
        DynamicMetadata%instr(co2)%hpath_length
    if (DynamicMetadataOrder(co2_irga_vpath_length) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(co2_irga_vpath_length)), *) &
        DynamicMetadata%instr(co2)%vpath_length
    if (DynamicMetadataOrder(co2_irga_tau) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(co2_irga_tau)), *) &
        DynamicMetadata%instr(co2)%tau

    !> h2o irga
    if (DynamicMetadataOrder(h2o_irga_manufacturer) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(h2o_irga_manufacturer)), *) &
        DynamicMetadata%instr(h2o)%firm
    if (DynamicMetadataOrder(h2o_irga_model) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(h2o_irga_model)), *) &
        DynamicMetadata%instr(h2o)%model
    if (DynamicMetadataOrder(h2o_measure_type) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(h2o_measure_type)), *) &
        DynamicMetadata%measure_type(h2o)
    if (DynamicMetadataOrder(h2o_irga_northward_separation) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(h2o_irga_northward_separation)), *) &
        DynamicMetadata%instr(h2o)%nsep
    if (DynamicMetadataOrder(h2o_irga_eastward_separation) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(h2o_irga_eastward_separation)), *) &
        DynamicMetadata%instr(h2o)%esep
    if (DynamicMetadataOrder(h2o_irga_vertical_separation) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(h2o_irga_vertical_separation)), *) &
        DynamicMetadata%instr(h2o)%vsep
    if (DynamicMetadataOrder(h2o_irga_tube_length) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(h2o_irga_tube_length)), *) &
        DynamicMetadata%instr(h2o)%tube_l
    if (DynamicMetadataOrder(h2o_irga_tube_diameter) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(h2o_irga_tube_diameter)), *) &
        DynamicMetadata%instr(h2o)%tube_d
    if (DynamicMetadataOrder(h2o_irga_tube_flowrate) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(h2o_irga_tube_flowrate)), *) &
        DynamicMetadata%instr(h2o)%tube_f
    if (DynamicMetadataOrder(h2o_irga_kw) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(h2o_irga_kw)), *) &
        DynamicMetadata%instr(h2o)%kw
    if (DynamicMetadataOrder(h2o_irga_ko) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(h2o_irga_ko)), *) &
        DynamicMetadata%instr(h2o)%ko
    if (DynamicMetadataOrder(h2o_irga_hpath_length) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(h2o_irga_hpath_length)), *) &
        DynamicMetadata%instr(h2o)%hpath_length
    if (DynamicMetadataOrder(h2o_irga_vpath_length) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(h2o_irga_vpath_length)), *) &
        DynamicMetadata%instr(h2o)%vpath_length
    if (DynamicMetadataOrder(h2o_irga_tau) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(h2o_irga_tau)), *) &
        DynamicMetadata%instr(h2o)%tau

    !> ch4 irga
    if (DynamicMetadataOrder(ch4_irga_manufacturer) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(ch4_irga_manufacturer)), *) &
        DynamicMetadata%instr(ch4)%firm
    if (DynamicMetadataOrder(ch4_irga_model) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(ch4_irga_model)), *) &
        DynamicMetadata%instr(ch4)%model
    if (DynamicMetadataOrder(ch4_measure_type) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(ch4_measure_type)), *) &
        DynamicMetadata%measure_type(ch4)
    if (DynamicMetadataOrder(ch4_irga_northward_separation) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(ch4_irga_northward_separation)), *) &
        DynamicMetadata%instr(ch4)%nsep
    if (DynamicMetadataOrder(ch4_irga_eastward_separation) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(ch4_irga_eastward_separation)), *) &
        DynamicMetadata%instr(ch4)%esep
    if (DynamicMetadataOrder(ch4_irga_vertical_separation) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(ch4_irga_vertical_separation)), *) &
        DynamicMetadata%instr(ch4)%vsep
    if (DynamicMetadataOrder(ch4_irga_tube_length) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(ch4_irga_tube_length)), *) &
        DynamicMetadata%instr(ch4)%tube_l
    if (DynamicMetadataOrder(ch4_irga_tube_diameter) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(ch4_irga_tube_diameter)), *) &
        DynamicMetadata%instr(ch4)%tube_d
    if (DynamicMetadataOrder(ch4_irga_tube_flowrate) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(ch4_irga_tube_flowrate)), *) &
        DynamicMetadata%instr(ch4)%tube_f
    if (DynamicMetadataOrder(ch4_irga_kw) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(ch4_irga_kw)), *) &
        DynamicMetadata%instr(ch4)%kw
    if (DynamicMetadataOrder(ch4_irga_ko) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(ch4_irga_ko)), *) &
        DynamicMetadata%instr(ch4)%ko
    if (DynamicMetadataOrder(ch4_irga_hpath_length) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(ch4_irga_hpath_length)), *) &
        DynamicMetadata%instr(ch4)%hpath_length
    if (DynamicMetadataOrder(ch4_irga_vpath_length) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(ch4_irga_vpath_length)), *) &
        DynamicMetadata%instr(ch4)%vpath_length
    if (DynamicMetadataOrder(ch4_irga_tau) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(ch4_irga_tau)), *) &
        DynamicMetadata%instr(ch4)%tau

    !> 4th gas irga
    if (DynamicMetadataOrder(gas4_irga_manufacturer) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(gas4_irga_manufacturer)), *) &
        DynamicMetadata%instr(gas4)%firm
    if (DynamicMetadataOrder(gas4_irga_model) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(gas4_irga_model)), *) &
        DynamicMetadata%instr(gas4)%model
    if (DynamicMetadataOrder(gas4_measure_type) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(gas4_measure_type)), *) &
        DynamicMetadata%instr(gas4)%nsep
    if (DynamicMetadataOrder(gas4_irga_northward_separation) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(gas4_irga_northward_separation)), *) &
        DynamicMetadata%instr(gas4)%nsep
    if (DynamicMetadataOrder(gas4_irga_eastward_separation) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(gas4_irga_eastward_separation)), *) &
        DynamicMetadata%instr(gas4)%esep
    if (DynamicMetadataOrder(gas4_irga_vertical_separation) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(gas4_irga_vertical_separation)), *) &
        DynamicMetadata%instr(gas4)%vsep
    if (DynamicMetadataOrder(gas4_irga_tube_length) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(gas4_irga_tube_length)), *) &
        DynamicMetadata%instr(gas4)%tube_l
    if (DynamicMetadataOrder(gas4_irga_tube_diameter) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(gas4_irga_tube_diameter)), *) &
        DynamicMetadata%instr(gas4)%tube_d
    if (DynamicMetadataOrder(gas4_irga_tube_flowrate) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(gas4_irga_tube_flowrate)), *) &
        DynamicMetadata%instr(gas4)%tube_f
    if (DynamicMetadataOrder(gas4_irga_kw) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(gas4_irga_kw)), *) &
        DynamicMetadata%instr(gas4)%kw
    if (DynamicMetadataOrder(gas4_irga_ko) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(gas4_irga_ko)), *) &
        DynamicMetadata%instr(gas4)%ko
    if (DynamicMetadataOrder(gas4_irga_hpath_length) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(gas4_irga_hpath_length)), *) &
        DynamicMetadata%instr(gas4)%hpath_length
    if (DynamicMetadataOrder(gas4_irga_vpath_length) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(gas4_irga_vpath_length)), *) &
        DynamicMetadata%instr(gas4)%vpath_length
    if (DynamicMetadataOrder(gas4_irga_tau) /= nint(error)) &
        read(mdStringVars(DynamicMetadataOrder(gas4_irga_tau)), *) &
        DynamicMetadata%instr(gas4)%tau

end subroutine ReadMetadataFromTextVars

!***************************************************************************
!
! \brief       Fills necessary missing information on dynamic metadata
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine FixDynamicMetadata()
    use m_rp_global_var
    implicit none
    !> local variables
    integer :: j

    do j = 1, E2NumVar
        !> Adjust instrument firms
        select case (DynamicMetadata%instr(j)%model(1:len_trim(DynamicMetadata%instr(j)%model) - 2))
            case ('li6262','li7000','li7200','li7500','li7500a','li7700')
                DynamicMetadata%instr(j)%firm = 'licor'
            case ('generic_open_path', 'generic_closed_path', &
                    'open_path_krypton', 'open_path_lyman', &
                    'closed_path_krypton', 'closed_path_lyman')
                DynamicMetadata%instr(j)%firm = 'other_irga'
            case('hs_50', 'hs_100', 'r2', 'r3_50', 'r3_100', &
                'r3a_100', 'wm', 'wmpro')
                DynamicMetadata%instr(j)%firm = 'gill'
            case('usa1_standard', 'usa1_fast')
                DynamicMetadata%instr(j)%firm = 'metek'
            case('csat3')
                DynamicMetadata%instr(j)%firm = 'csi'
            case('81000')
                DynamicMetadata%instr(j)%firm = 'young'
            case('generic_sonic')
                DynamicMetadata%instr(j)%firm = 'other_sonic'
        end select
        !> Determine whether the instrument is a sonic or a gas analyser
        select case(DynamicMetadata%instr(j)%firm)
            case('licor', 'other_irga')
                DynamicMetadata%instr(j)%category = 'irga'
            case('gill', 'metek', 'young', 'csi', 'other_sonic')
                DynamicMetadata%instr(j)%category = 'sonic'
        end select
        !> If format of wind components is not set, it is likely to be u,v,w
        select case (DynamicMetadata%instr(j)%wformat)
            case ('uvw', 'polar_w', 'axis')
                continue
            case ('none')
                continue
            case default
                DynamicMetadata%instr(j)%wformat = 'uvw'
        end select
        !> If sonic is a Gill and wref is unknown, assumes SPAR configuration
        !> (the most logic, with U aligned to North. In case of wrong guess,
        !> there is only 30 degree offset in wind speed)
        if (DynamicMetadata%instr(j)%firm == 'gill') then
            select case (DynamicMetadata%instr(j)%wref)
                case ('axis','spar')
                    continue
                case ('none')
                    continue
                case default
                    DynamicMetadata%instr(j)%wref = 'spar'
            end select
        end if

        if (DynamicMetadata%instr(j)%category == 'irga') then
            select case (DynamicMetadata%instr(j)%model(1:len_trim(DynamicMetadata%instr(j)%model) - 2))
                case ('li7700', 'li7500', 'li7500a', 'generic_open_path', &
                    'open_path_krypton', 'open_path_lyman')
                    DynamicMetadata%instr(j)%path_type = 'open'
                case ('none')
                    continue
                case default
                    DynamicMetadata%instr(j)%path_type = 'closed'
            end select
        end if

        !> Retrieve gas analyser parameters
        if (DynamicMetadata%instr(j)%tube_d /= error) &
            DynamicMetadata%instr(j)%tube_d = &
            DynamicMetadata%instr(j)%tube_d * 1d-3 !< in meters
        if (DynamicMetadata%instr(j)%tube_l /= error) &
            DynamicMetadata%instr(j)%tube_l = &
            DynamicMetadata%instr(j)%tube_l * 1d-2  !< in meters
        if (DynamicMetadata%instr(j)%tube_f /= error) &
            DynamicMetadata%instr(j)%tube_f = &
            DynamicMetadata%instr(j)%tube_f / 6d4  !< in m+3s-1

        if (DynamicMetadata%instr(j)%nsep /= error) &
            DynamicMetadata%instr(j)%nsep = &
            DynamicMetadata%instr(j)%nsep * 1d-2 !< in meters
        if (DynamicMetadata%instr(j)%esep /= error) &
            DynamicMetadata%instr(j)%esep = &
            DynamicMetadata%instr(j)%esep * 1d-2 !< in meters
        if (DynamicMetadata%instr(j)%vsep /= error) &
            DynamicMetadata%instr(j)%vsep = &
            DynamicMetadata%instr(j)%vsep * 1d-2 !< in meters
    end do

end subroutine FixDynamicMetadata

!***************************************************************************
!
! \brief       Pick dynamic metadata according to timestamp info
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ExtractUsableMetadataFromDynamic(LocCol, ncol)
    use m_rp_global_var
    implicit none
    !> In/out variables
    integer, intent(in) :: ncol
    type(ColType), intent(inout) :: LocCol(ncol)
    !> Local variables
    integer :: gas
    integer :: gas2
    logical :: instr_updated(GHGNumVar)


    !> General
    if (DynamicMetadata%lat /= error) Metadata%lat = DynamicMetadata%lat
    if (DynamicMetadata%lon /= error) Metadata%lon = DynamicMetadata%lon
    if (DynamicMetadata%alt /= error) Metadata%alt = DynamicMetadata%alt
    if (DynamicMetadata%canopy_height /= error) Metadata%canopy_height = DynamicMetadata%canopy_height
    if (DynamicMetadata%d /= error) Metadata%d = DynamicMetadata%d
    if (DynamicMetadata%z0 /= error) Metadata%z0 = DynamicMetadata%z0
    if (DynamicMetadata%ac_freq /= error) Metadata%ac_freq = DynamicMetadata%ac_freq
    if (DynamicMetadata%file_length /= error) Metadata%file_length = DynamicMetadata%file_length
    !> Sonic
    if (DynamicMetadata%instr(u)%model /= 'none') LocCol(u:ts)%instr%model = DynamicMetadata%instr(u)%model
    if (DynamicMetadata%instr(u)%wformat /= 'none') LocCol(u:ts)%instr%wformat = DynamicMetadata%instr(u)%wformat
    if (DynamicMetadata%instr(u)%wref /= 'none') LocCol(u:ts)%instr%wref = DynamicMetadata%instr(u)%wref
    if (DynamicMetadata%instr(u)%height /= error) LocCol(u:ts)%instr%height = DynamicMetadata%instr(u)%height
    if (DynamicMetadata%instr(u)%north_offset /= error) LocCol(u:ts)%instr%north_offset = DynamicMetadata%instr(u)%north_offset
    !> Gases
    instr_updated = .false.
    do gas = co2, gas4
        if (DynamicMetadata%instr(gas)%path_type /= 'none') then
            instr_updated(gas) = .true.
            LocCol(gas)%instr%path_type = DynamicMetadata%instr(gas)%path_type
        end if
        if (DynamicMetadata%instr(gas)%model /= 'none') then
            instr_updated(gas) = .true.
            LocCol(gas)%instr%model = DynamicMetadata%instr(gas)%model
        end if
        if (DynamicMetadata%instr(gas)%nsep /= error) then
            instr_updated(gas) = .true.
            LocCol(gas)%instr%nsep = DynamicMetadata%instr(gas)%nsep
        end if
        if (DynamicMetadata%instr(gas)%esep /= error) then
            instr_updated(gas) = .true.
            LocCol(gas)%instr%esep = DynamicMetadata%instr(gas)%esep
        end if
        if (DynamicMetadata%instr(gas)%vsep /= error) then
            instr_updated(gas) = .true.
            LocCol(gas)%instr%vsep = DynamicMetadata%instr(gas)%vsep
        end if
        if (DynamicMetadata%instr(gas)%tube_l /= error) then
            instr_updated(gas) = .true.
            LocCol(gas)%instr%tube_l = DynamicMetadata%instr(gas)%tube_l
        end if
        if (DynamicMetadata%instr(gas)%tube_d /= error) then
            instr_updated(gas) = .true.
            LocCol(gas)%instr%tube_d = DynamicMetadata%instr(gas)%tube_d
        end if
        if (DynamicMetadata%instr(gas)%tube_f /= error) then
            instr_updated(gas) = .true.
            LocCol(gas)%instr%tube_f = DynamicMetadata%instr(gas)%tube_f
        end if
        if (DynamicMetadata%instr(gas)%kw /= error) then
            instr_updated(gas) = .true.
            LocCol(gas)%instr%kw = DynamicMetadata%instr(gas)%kw
        end if
        if (DynamicMetadata%instr(gas)%ko /= error) then
            instr_updated(gas) = .true.
            LocCol(gas)%instr%ko = DynamicMetadata%instr(gas)%ko
        end if
        if (DynamicMetadata%instr(gas)%hpath_length /= error) then
            instr_updated(gas) = .true.
            LocCol(gas)%instr%hpath_length = DynamicMetadata%instr(gas)%hpath_length
        end if
        if (DynamicMetadata%instr(gas)%vpath_length /= error) then
            instr_updated(gas) = .true.
            LocCol(gas)%instr%vpath_length = DynamicMetadata%instr(gas)%vpath_length
        end if
        if (DynamicMetadata%instr(gas)%tau /= error) then
            instr_updated(gas) = .true.
            LocCol(gas)%instr%tau = DynamicMetadata%instr(gas)%tau
        end if
    end do

    !> Consideration of gases from same analyser
    do gas = co2, gas4
        do gas2 = co2, gas4
            if ((LocCol(gas2)%instr%model == LocCol(gas)%instr%model) .and. instr_updated(gas)) &
                LocCol(gas2)%instr = LocCol(gas)%instr
        end do
    end do

end subroutine ExtractUsableMetadataFromDynamic
