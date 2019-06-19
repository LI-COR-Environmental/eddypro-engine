!***************************************************************************
! write_out_fluxnet_only_biomet.f90
! ---------------------------------
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
! \brief       Write line to FLUXNET output with only biomet data, if available
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WriteOutFluxnetOnlyBiomet()
    use m_rp_global_var
    implicit none
    !> local variables
    integer :: i
    integer :: indx
    integer :: int_doy
    real(kind = dbl) :: float_doy
    character(32) :: char_doy
    character(LongOutstringLen) :: dataline
    character(14) :: tsIso
    real(kind = dbl), allocatable :: bAggrOut(:)
    include '../src_common/interfaces.inc'

    call clearstr(dataline)

    !> Start/end imestamps
    tsIso = Stats%start_date(1:4) // Stats%start_date(6:7) // Stats%start_date(9:10) &
                // Stats%start_time(1:2) // Stats%start_time(4:5)
    call AddDatum(dataline, trim(adjustl(tsIso)), separator)
    tsIso = Stats%date(1:4) // Stats%date(6:7) // Stats%date(9:10) &
                // Stats%time(1:2) // Stats%time(4:5)
    call AddDatum(dataline, trim(adjustl(tsIso)), separator)

    !> DOYs
    !>  Start
    call DateTimeToDOY(Stats%start_date, Stats%start_time, int_doy, float_doy)
    write(char_doy, *) float_doy
    call AddDatum(dataline, trim(adjustl(char_doy(1: index(char_doy, '.')+ 4))), separator)
    !>  End
    call DateTimeToDOY(Stats%date, Stats%time, int_doy, float_doy)
    write(char_doy, *) float_doy
    call AddDatum(dataline, trim(adjustl(char_doy(1: index(char_doy, '.')+ 4))), separator)

    !> Not enough data
    call AddDatum(dataline, 'not_enough_data', separator)

    !> Potential Radiations
    indx = DateTimeToHalfHourNumber(Stats%date, Stats%time)
    call AddFloatDatumToDataline(PotRad(indx), dataline, EddyProProj%err_label)

    !> Daytime
    if (Stats%daytime) then
        call AddDatum(dataline, '0', separator)
    else
        call AddDatum(dataline, '1', separator)
    endif

    !> Write error codes in place of fixed columns
    do i = 1, 456
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    end do
    !> Write error codes in place of custom variables
    do i = 1, NumUserVar + 1
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    end do

    !> write all aggregated biomet values in FLUXNET units
    call AddIntDatumToDataline(nbVars, dataline, EddyProProj%err_label)
    if (nbVars > 0) then
        if (.not. allocated(bAggrOut)) allocate(bAggrOut(size(bAggr)))
        if (EddyProProj%fluxnet_standardize_biomet) then
            bAggrOut = bAggrFluxnet
        else
            bAggrOut = bAggr
        end if

        do i = 1, nbVars
            call AddFloatDatumToDataline(bAggrOut(i), dataline, EddyProProj%err_label)
        end do
    end if
    write(uflxnt, '(a)') dataline(1:len_trim(dataline) - 1)

end subroutine WriteOutFluxnetOnlyBiomet
