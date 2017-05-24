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
subroutine WriteOutIcosOnlyBiomet(init_string)
    use m_rp_global_var
    implicit none
    !> in/out variables
    character(*), intent(in) :: init_string
    !> local variables
    integer :: i
    character(LongOutstringLen) :: dataline
    character(DatumLen) :: datum
    character(64) :: tmp_init_string
    character(14) :: iso_basic
    include '../src_common/interfaces.inc'

    !> write Essentials output file (csv) for communication
    !> with Fluxes
    call clearstr(dataline)

    !> Timestamp
    tmp_init_string = &
        init_string(index(init_string, ',') +1: &
                    index(init_string, ',', .true.) - 1)
    iso_basic = tmp_init_string(1:4) // tmp_init_string(6:7) &
        // tmp_init_string(9:10) // tmp_init_string(12:13)  &
        // tmp_init_string(15:16) // '00'
    call AddDatum(dataline, trim(adjustl(iso_basic)), separator)

    !> Write error codes in place of fixed columns
    do i = 1, 481
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    end do
    !> Write error codes in place of custom variables
    do i = 1, NumUserVar + 1
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    end do

    !> write all aggregated biomet values in FLUXNET units
    write(datum, *) nbVars
    call AddDatum(dataline, datum, separator)
    do i = 1, nbVars
        call WriteDatumFloat(bAggrFluxnet(i), datum, '-9999.')
        call AddDatum(dataline, datum, separator)
    end do
    write(uicos, '(a)') dataline(1:len_trim(dataline) - 1)

end subroutine WriteOutIcosOnlyBiomet
