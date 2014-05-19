!***************************************************************************
! init_env.f90
! ------------
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
! \brief       Initialize environment variables
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine InitEnv()
    use m_common_global_var
    implicit none
    include 'version_and_date.inc'
    !> local variables
    integer, parameter :: n_cmdpar = 3
    integer :: i
    integer :: make_dir
    integer :: aux
    character(32) :: timestring
    character(8) :: switch
    character(256) :: value
    character(32) :: tmpDirPadding
    character(3), parameter :: OS_default = 'win'
    integer, external :: CreateDir


    !> Store current timestamp information
    call hms_current_string(timestring)
    if(timestring(12:12) == ' ') timestring(12:12) = '0'
    Timestamp_FilePadding = '_' // timestring(1:10) // 'T' &
        // timestring(12:13) // timestring(15:16) // timestring(18:19)

    tmpDirPadding = '_' // timestring(1:10) // 'T' &
        // timestring(12:13) // timestring(15:16) // timestring(18:19) // timestring(21:23)

    !> Get command lines switches and values
    OS = ''
    homedir = ''
    EddyProProj%run_env = ''
    do i = 1, n_cmdpar
        call get_command_argument(2*(i-1)+1, switch)
        call get_command_argument(2*(i-1)+2, value)
        if (switch(1:len_trim(switch)) == '-s') then
            !> Switch for "system", the host operating system
            OS = trim(value)
            if (OS(1:1) == '-') OS = ''
        elseif (switch(1:len_trim(switch)) == '-e') then
            !> Switch for "environment", the 'home' working directory
            homedir = trim(value)
            if (homedir(1:1) == '-') homedir = ''
        elseif (switch(1:len_trim(switch)) == '-m') then
                !> Switch for "mode", whether "embedded" or "desktop" mode
                EddyProProj%run_env = trim(value)
                if (EddyProProj%run_env(1:1) == '-') EddyProProj%run_env = ''
        elseif (switch(1:len_trim(switch)) == '-v') then
            call InformOfSoftwareVersion(sw_ver, build_date)
        elseif (switch(1:len_trim(switch)) == '-h') then
            call CommandLineHelp(sw_ver, build_date)
        end if
    end do

    !> Set OS-dependent parameters
    if (len_trim(OS) == 0) OS = OS_default
    call SetOSEnvironment()
    if (len_trim(homedir) == 0) homedir = '..'
    if (len_trim(EddyProProj%run_env) == 0) EddyProProj%run_env = 'desktop'

    !> Define default unit number (udf), run specific
    call hms_current_hms(aux, aux, aux, udf)
    if (udf < 200) udf = udf + 200
    udf2 = udf + 1

    !> Define path of key eddypro files/dirs
    call AdjDir(homedir, slash)
    IniDir = trim(homedir) // 'ini' // slash
    LogDir = trim(homedir) // 'log' // slash
    PrjPath = trim(IniDir) // trim(PrjFile)

    !> Define TmpDir differently if it's in desktop or embedded mode
    if (EddyProProj%run_env == 'desktop') then
        TmpDir = trim(homedir) // 'tmp' // slash // 'tmp' // trim(adjustl(tmpDirPadding)) // slash
    else
        TmpDir = trim(homedir) // 'tmp' // slash
    end if

    select case (app)
        case ('EddyPro-RP')
        LogPath = trim(LogDir) // trim(LogFileRP)
        case ('EddyPro-FCC')
        LogPath = trim(LogDir) // trim(LogFileFX)
    end select

    !> Create log dir in case it doesn't exist (for use from command line)
    make_dir = CreateDir('"' // trim(LogDir) // '"')
    make_dir = CreateDir('"' // trim(TmpDir) // '"')

    !> Initialize log file with append disabled.
    call log_startup(trim(LogPath), .false.)
    call log_configure("timestamp", .false.)
    call log_configure("writeonstdout", .false.)
    call log_configure("stoponerror", .false.)
    call log_delimiter(LOG_LEVEL_VOLUME)
    call log_msg(' inf=executing ' // trim(app))
    call log_delimiter(LOG_LEVEL_VOLUME)
    LogString = ' inf=hosting operating system: ' // trim(OS)
    call log_msg(LogString)
    call clearstr(LogString)
end subroutine InitEnv

!***************************************************************************
! \file        src/init_env.f90
! \brief       Provide software version on demand (switch "-v")
! \version     4.2.0
! \date        2013-07-02
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine InformOfSoftwareVersion(sw_ver, build_date)
    use m_common_global_var
    implicit none
    !> In/out variables
    character(*), intent(in) :: sw_ver
    character(*), intent(in) :: build_date


    write (*, '(a)') ' ' // trim(adjustl(app)) // ', version ' // trim(adjustl(sw_ver)) // &
        &', build ' // trim(adjustl(build_date)) // '.'
    stop
end subroutine InformOfSoftwareVersion

!***************************************************************************
! \file        src/init_env.f90
! \brief       Provide software version on demand (switch "-v")
! \version     4.2.0
! \date        2013-07-02
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine CommandLineHelp(sw_ver, build_date)
    use m_common_global_var
    implicit none
    !> In/out variables
    character(*), intent(in) :: sw_ver
    character(*), intent(in) :: build_date


    write(*, '(a)') ' Help for ' // trim(adjustl(app))
    write(*, '(a)') ' --------------------'
    write (*, '(a)') ' ' // trim(adjustl(app)) // ', version ' // trim(adjustl(sw_ver)) // &
        &', build ' // trim(adjustl(build_date)) // '.'
    write(*,*)
    write(*, '(a)') ' USAGE: eddypro_rp [OPTION [ARG]]'
    write(*,*)
    write(*, '(a)') ' OPTIONS:'
    write(*, '(a)') '   -s [win | linux | mac]    Operating system; if not provided assumes "win"'
    write(*, '(a)') '   -m [embbeded | desktop]   Running mode; if not provided assumes "desktop"'
    write(*, '(a)') '   -e [DIRECTORY]            Working directory, to be provided in embedded mode; if not provided assumes \.'
    write(*, '(a)') '   -h                        Display this help and exit'
    write(*, '(a)') '   -v                        Output version information and exit'
    stop
end subroutine CommandLineHelp
