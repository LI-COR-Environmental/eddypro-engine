!***************************************************************************
! write_out_icos_only_biomet.f90
! ------------------------------
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
! \brief       Write line to ICOS output with only biomet data, if available
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WriteOutIcosOnlyBiomet()
    use m_rp_global_var
    implicit none
    !> local variables
    integer :: i
    character(LongOutstringLen) :: dataline
    character(14) :: tsIso
    real(kind = dbl), allocatable :: bAggrOut(:)
    include '../src_common/interfaces.inc'

    !> write Essentials output file (csv) for communication
    !> with Fluxes
    call clearstr(dataline)

    !> Start/end imestamps
    tsIso = Stats%start_date(1:4) // Stats%start_date(6:7) // Stats%start_date(9:10) &
                // Stats%start_time(1:2) // Stats%start_time(4:5)
    call AddDatum(dataline, trim(adjustl(tsIso)), separator)
    tsIso = Stats%date(1:4) // Stats%date(6:7) // Stats%date(9:10) &
                // Stats%time(1:2) // Stats%time(4:5)
    call AddDatum(dataline, trim(adjustl(tsIso)), separator)

    !> Write error codes in place of fixed columns
    do i = 1, 471
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
        if (EddyProProj%icos_standardize_biomet) then
            bAggrOut = bAggrFluxnet
        else
            bAggrOut = bAggr
        end if

        do i = 1, nbVars
            call AddFloatDatumToDataline(bAggrOut(i), dataline, EddyProProj%err_label)
        end do
    end if
    write(uicos, '(a)') dataline(1:len_trim(dataline) - 1)

end subroutine WriteOutIcosOnlyBiomet
