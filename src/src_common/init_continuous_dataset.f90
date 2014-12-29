!***************************************************************************
! init_continuous_dataset.f90
! ---------------------------
! Copyright (C) 2007-2011, Eco2s team, Gerardo Fratini
! Copyright (C) 2011-2014, LI-COR Biosciences
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
! \brief       Initialize continuous dataset by creating an error string \n
!              to replace missing averaging periods in a continuous dataset
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine InitContinuousDataset(PathIn, ErrString, hnrow)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: hnrow
    character(*), intent(in) :: PathIn
    character(*), intent(out) :: ErrString
    !> local variabes
    integer :: sepa
    integer :: io_status
    integer :: read_status
    integer :: i
    character(MaxOutstringLen) :: Datastring


    !> Open file
    open(udf, file = PathIn(1:len_trim(PathIn)), status = 'old', iostat = io_status)
    if (io_status /= 0 ) then
        write(*,'(a)')
        write(*,'(a)') '  A problem occurred while opening file: ', PathIn(1:len_trim(PathIn))
        write(*,'(a)') '   File not imported.'
        return
    end if

    !> Skip header
    do i = 1, hnrow
        read(udf,*)
    end do

    !> Now import first valid Datastring
    do
        read(udf, '(a)', iostat = read_status) Datastring
        if (read_status > 0) cycle
        exit
    end do

    ErrString = ','
    !> skips everything before the end of the date/time/DOY
    sepa = index(Datastring, ',')
    Datastring = Datastring(sepa + 1: len_trim(Datastring))
    sepa = index(Datastring, ':') + 3
    Datastring = Datastring(sepa + 1: len_trim(Datastring))
    sepa = index(Datastring, ',')
    Datastring = Datastring(sepa + 1: len_trim(Datastring))

    !> Now Datastring starts from first sensible datum
    !> Substitute datum with errors
    do
        sepa = index(Datastring, ',')
        if (sepa /= 0) then
            ErrString = ErrString(1:len_trim(ErrString)) // EddyProProj%err_label(1:len_trim(EddyProProj%err_label)) // ','
            Datastring = Datastring(sepa + 1: len_trim(Datastring))
        else
            exit
        end if
    end do
    !> adds the last one
    ErrString = ErrString(1:len_trim(ErrString)) // EddyProProj%err_label(1:len_trim(EddyProProj%err_label))
    close(udf)
end subroutine InitContinuousDataset

