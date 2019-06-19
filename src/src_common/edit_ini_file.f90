!***************************************************************************
! edit_ini_file.f90
! -----------------
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
    character(PathLen) :: tfname
    character(ShortInstringLen) :: dataline
    character(64) :: currtag


    !> Open file to be edited and temp file
    open(10, file=trim(fname), iostat=io_error)
    tfname = trim(fname) // '.tmp'
    open(11, file=trim(tfname), status='new', iostat=io_error)

    !> Copy whole file into temp file including modification
    do
        read(10, '(a)', iostat=ierr) dataline
        if (ierr < 0) exit
        sepa = index(dataline, '=')
        if (sepa /= 0) then
            currtag = dataline(1:sepa-1)
            if (currtag == trim(tag)) dataline = trim(tag) // '='// trim(newval)
        end if
        write(11, '(a)') trim(adjustl(dataline))
    end do
    close(10, status='DELETE')
    close(11)

    !> Input file has been deleted, now change name of tmp file into old input file name
    call rename(trim(adjustl(tfname)), trim(adjustl(fname)), status = io_error)
end subroutine EditIniFile

