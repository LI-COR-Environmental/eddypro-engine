!***************************************************************************
! edit_ini_file.f90
! -----------------
!Copyright (C) 2014, LI-COR Biosciences
!
!This file is part of EddyPro (TM).
!
!EddyPro (TM) is free software: you can redistribute it and/or modify
!it under the terms of the GNU General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.
!
!EddyPro (TM) is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License
!along with EddyPro (TM).  If not, see <http://www.gnu.org/licenses/>.
!***************************************************************************
!
! \brief       Modify INI file, changing the value for the specified tag
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine EditIniFile(fname, tag, newval)
    use m_common_global_var
    !> In/out variables
    character(*), intent(in) :: fname
    character(*), intent(in) :: tag
    character(*), intent(in) :: newval
    !> Local variables
    integer :: ierr
    integer :: sepa
    integer :: io_error
    character(256) :: tfname
    character(512) :: row
    character(64) :: currtag


    !> Open file to be edited and temp file
    open(10, file=trim(fname), iostat=io_error)
    tfname = trim(fname) // '.tmp'
    open(11, file=trim(tfname), status='new', iostat=io_error)

    !> Copy whole file into temp file including modification
    do
        read(10, '(a)', iostat=ierr) row
        if (ierr < 0) exit
        sepa = index(row, '=')
        if (sepa /= 0) then
            currtag = row(1:sepa-1)
            if (currtag == trim(tag)) row = trim(tag) // '='// trim(newval)
        end if
        write(11, '(a)') trim(adjustl(row))
    end do
    close(10, status='DELETE')
    close(11)

    !> Input file has been deleted, now change name of tmp file into old input file name
    call rename(trim(adjustl(tfname)), trim(adjustl(fname)), status = io_error)
end subroutine EditIniFile

