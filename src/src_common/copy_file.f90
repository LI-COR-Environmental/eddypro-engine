!***************************************************************************
! copy_file.f90
! -------------
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
!EddyPro® is free software: you can redistribute it and/or modify
!it under the terms of the GNU General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.
!EddyPro® is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!GNU General Public License for more details.
!You should have received a copy of the GNU General Public License
!along with EddyPro®.  If not, see <http://www.gnu.org/licenses/>.
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
