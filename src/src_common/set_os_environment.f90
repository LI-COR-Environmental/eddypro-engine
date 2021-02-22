!***************************************************************************
! set_os_environment.f90
! ----------------------
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
            comm_7zip_out     = '-o'
            comm_copy         = 'copy '
            comm_move         = 'move '
            comm_force_opt    = '/Y '
            comm_dir          = 'dir  /O:D /B '
            comm_mkdir        = 'mkdir'
        case('linux')
            slash = '/'
            escape = '\'
            comm_err_redirect = ' 2> /dev/null'
            comm_out_redirect = ' > /dev/null'
            comm_del          = 'rm '
            comm_rmdir        = 'rm -r -f'
            comm_7zip         = '7za '
            comm_7zip_x_opt   = 'x -y '
            comm_7zip_out     = '-o'
            comm_copy         = 'cp '
            comm_move         = 'mv '
            comm_force_opt    = '-y '
            comm_dir          = 'ls '
!            comm_dir          = 'find -iname '
            comm_mkdir        = 'mkdir -p'
        case('mac')
            slash = '/'
            escape = '\'
            comm_err_redirect = ' 2> /dev/null'
            comm_out_redirect = ' > /dev/null'
            comm_del          = 'rm '
            comm_rmdir        = 'rm -r -f'
            ! comm_7zip         = './7za '
            ! comm_7zip_x_opt   = 'x -y '
            ! comm_7zip_out     = '-o '
            comm_7zip         = 'unzip '
            comm_7zip_x_opt   = '-o '
            comm_7zip_out     = '-d '
            comm_copy         = 'cp '
            comm_move         = 'mv '
            comm_force_opt    = '-y '
            comm_dir          = 'ls '
!            comm_dir          = 'find -iname '
            comm_mkdir        = 'mkdir -p'
    end select
end subroutine SetOSEnvironment
