!***************************************************************************
! copy_file.f90
! -------------
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
! \brief       Copy a text file by reading/writing content line by line
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine CopyFile(ifname, ofname)
    use m_common_global_var
    !> In/out variables
    character(*), intent(in) :: ifname
    character(*), intent(in) :: ofname
    !> Local variables
    integer :: io_error
    character(LongInstringLen) :: dataline

    !> Open existing file
    open(10, file = trim(ifname), status = 'old', iostat = io_error)
    if (io_error /= 0) return

    !> Open new file
    open(11, file = trim(ofname), status = 'new', iostat = io_error)
    if (io_error /= 0) return

    !> Copy file line by line
    io_error = 0
    do
        read(10, '(a)', iostat = io_error) dataline
        if(io_error /= 0) exit
        write(11, '(a)') trim(dataline)
    end do
end subroutine CopyFile
