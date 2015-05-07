!***************************************************************************
! set_os_environment.f90
! ----------------------
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
! \brief       Set constants depending on Operating system
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine SetOSEnvironment()
    use m_common_global_var
    implicit none


    select case(OS(1:len_trim(OS)))
        case('win')
            slash = '\'
            escape = '/'
            comm_err_redirect = ' 2> nul'
            comm_out_redirect = '  > nul'
            comm_del          = 'del '
            comm_rmdir        = 'rmdir /s /q'
            comm_7zip         = '7z.exe '
            comm_7zip_x_opt   = 'x -y '
            comm_copy         = 'copy '
            comm_move         = 'move '
            comm_force_opt    = '/Y '
            comm_dir          = 'dir  /O:D /B '
        case('linux')
            slash = '/'
            escape = '\'
            comm_err_redirect = ' 2> /dev/null'
            comm_out_redirect = ' > /dev/null'
            comm_del          = 'rm '
            comm_rmdir        = 'rm -r -f'
            comm_7zip         = './7za '
            comm_7zip_x_opt   = 'x -y '
            comm_copy         = 'cp '
            comm_move         = 'mv '
            comm_force_opt    = '-y '
            comm_dir          = 'ls '
!            comm_dir          = 'find -iname '
        case('mac')
            slash = '/'
            escape = '\'
            comm_err_redirect = ' 2> /dev/null'
            comm_out_redirect = ' > /dev/null'
            comm_del          = 'rm '
            comm_rmdir        = 'rm -r -f'
            comm_7zip         = './7za '
            comm_7zip_x_opt   = 'x -y '
            comm_copy         = 'cp '
            comm_move         = 'mv '
            comm_force_opt    = '-y '
            comm_dir          = 'ls '
!            comm_dir          = 'find -iname '
    end select
end subroutine SetOSEnvironment
