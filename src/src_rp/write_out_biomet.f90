!***************************************************************************
! write_out_biomet.f90
! --------------------
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
! \brief       Write biomet data on (temporary) output files
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WriteOutBiomet(init_string, embedded)
    use m_rp_global_var
    implicit none
    !> in/out variables
    character(*), intent(in) :: init_string
    logical, intent(in) :: embedded
    !> local variables
    integer :: i
    character(LongOutstringLen) :: dataline
    character(DatumLen) :: datum
    character(len=len(init_string)) :: prefix
    character(64) :: tmp_init_string
    character(14) :: iso_basic
    include '../src_common/interfaces.inc'

    !>==========================================================================
    !> EddyPro's BIOMET output
    if (embedded) then
        prefix = init_string(index(init_string, ',') + 1: &
                             len_trim(init_string))
    else
        prefix = trim(init_string)
    end if

    if (EddyProProj%out_biomet .and. nbVars > 0) then
        call clearstr(dataline)
        call AddDatum(dataline, trim(adjustl(prefix)), separator)

        do i = 1, nbVars
            call WriteDatumFloat(bAggrEddyPro(i), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        end do
        write(ubiomet, '(a)') dataline(1:len_trim(dataline) - 1)
    end if

    !>==========================================================================
    !> FLUXNET BIOMET output
    if (EddyProProj%out_fluxnet_biomet) then
        !> Edit init_string to fit FLUXNET format
        !> Strip DOY
        if (embedded) then
            tmp_init_string = init_string(index(init_string, ',') + 1: &
                index(init_string, ',', .true.) - 1)
        else
            tmp_init_string = &
                init_string(1: index(init_string, ',', .true.) - 1)
        end if

        !> derive ISO basic format timestamp
        iso_basic = tmp_init_string(1:4) // tmp_init_string(6:7) &
            // tmp_init_string(9:10) // tmp_init_string(12:13)  &
            // tmp_init_string(15:16) // '00'

        call clearstr(dataline)
        call AddDatum(dataline, trim(adjustl(iso_basic)), separator)

        !> All aggregated biomet values in FLUXNET units
        do i = 1, nbVars
            call WriteDatumFloat(bAggrFluxnet(i), datum, '-9999.')
            call AddDatum(dataline, datum, separator)
        end do
        write(ufnet_b, '(a)') dataline(1:len_trim(dataline) - 1)
    end if
end subroutine WriteOutBiomet
