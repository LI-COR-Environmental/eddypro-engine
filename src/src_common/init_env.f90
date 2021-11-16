!***************************************************************************
! init_env.f90
! ------------
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
    integer :: i
    integer :: make_dir
    integer :: io_status
    integer :: aux
    character(32) :: timestring
    character(256) :: switch
    character(PathLen) :: projPath
    character(256) :: arg
    character(32) :: tmpDirPadding
    real         :: realrand
    character(5) :: strrand
    character(PathLen) :: tempdir
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
    EddyProProj%caller = ''
    tempdir = ''
    projPath = ''
    i = 1
    arg_loop: do
        !> Read switch
        call get_command_argument(i, value=switch, status=io_status)
        if (io_status > 0 .or. len_trim(switch) == 0) exit arg_loop
        i = i + 1
        call get_command_argument(i, value=arg, status=io_status)
        i = i + 1

        if (switch(1:1) == '-') then
            select case(trim(adjustl(switch)))

                !> Switch for "system", the host operating system
                case('-s', '--system')
                    if (io_status > 0 .or. len_trim(switch) == 0) exit arg_loop
                    OS = trim(arg)
                    if (OS(1:1) == '-') OS = ''

                !> Switch for "environment", the 'home' working directory
                case('-e', '--environment')
                    if (io_status > 0 .or. len_trim(switch) == 0) exit arg_loop
                    homedir = trim(arg)
                    if (homedir(1:1) == '-') homedir = ''

                !> Switch for "mode", whether "embedded" or "desktop" mode
                case('-m', '--mode')
                    if (io_status > 0 .or. len_trim(switch) == 0) exit arg_loop
                    EddyProProj%run_env = trim(arg)
                    if (EddyProProj%run_env(1:1) == '-') EddyProProj%run_env = ''

                !> Switch for "caller", whether "gui" or "console"
                case('-c', '--caller')
                    if (io_status > 0 .or. len_trim(switch) == 0) exit arg_loop
                    EddyProProj%caller = trim(arg)
                    if (EddyProProj%caller(1:1) == '-') EddyProProj%caller = ''

                !> Switch for "temporary directory", the temporary working directory
                case('-t', '--tmpdir')
                    if (io_status > 0 .or. len_trim(switch) == 0) exit arg_loop
                    tempdir = trim(arg)
                    if (tempdir(1:1) == '-') tempdir = ''

                !> Software version
                case('-v', '--version')
                    call InformOfSoftwareVersion(sw_ver, build_date)

                !> Minimal command line help
                case('-h', '--help')
                    call CommandLineHelp(sw_ver, build_date)
            end select
        else
            projPath = trim(switch)
            if (index(projPath, '.eddypro') == 0) projPath = ''
        end if
    end do arg_loop

    !> Set OS-dependent parameters
    if (len_trim(OS) == 0) OS = OS_default
    call SetOSEnvironment()

    !> Default values if args are not passed
    if (len_trim(homedir) == 0) homedir = '..'
    if (len_trim(EddyProProj%run_env) == 0) EddyProProj%run_env = 'desktop'
    if (len_trim(EddyProProj%caller) == 0)  EddyProProj%caller  = 'console'
    if (len_trim(tempdir) == 0) then
        tempdir = homedir
        if (EddyProProj%run_env == 'desktop') &
            tempdir = trim(tempdir) // 'tmp' // slash
    endif

    !> Define default unit number (udf), run specific
    call hms_current_hms(aux, aux, aux, udf)
    if (udf < 200) udf = udf + 200
    udf2 = udf + 1

    !> Define path of key eddypro files/dirs
    call AdjDir(homedir, slash)
    IniDir = trim(homedir) // 'ini' // slash
    if (projPath == '') then
        PrjPath = trim(IniDir) // trim(PrjFile)
    else
        PrjPath = projPath
    end if

    !> Define TmpDir differently if it's in desktop or embedded mode
    call random_number(realrand)
    write(strrand, '(i0.5)') int(realrand*10000)
    if (EddyProProj%run_env == 'desktop') then
        TmpDir = trim(tempdir) // 'tmp' // trim(adjustl(tmpDirPadding)) &
            // '_' // strrand // slash
    else
        TmpDir = trim(tempdir) // 'tmp' // '_' // strrand // slash
    end if

    !> Create TmpDir in case it doesn't exist (for use from command line)
    make_dir = CreateDir('"' // trim(TmpDir) // '"')
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
!
! \brief       Provide software version on demand (switch "-v")
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
    !> Local variables
    character(11) :: prog


    if (app == 'EddyPro-RP') then
        prog = 'eddypro_rp'
    else
        prog = 'eddypro_fcc'
    end if

    write(*, '(a)') ' Help for ' // trim(adjustl(app))
    write(*, '(a)') ' --------------------'
    write (*, '(a)') ' ' // trim(adjustl(app)) // ', version ' // trim(adjustl(sw_ver)) // &
        &', build ' // trim(adjustl(build_date)) // '.'
    write(*,*)
    write(*, '(a)') ' USAGE: ' // trim(prog) // ' [OPTION [ARG]] [PROJ_FILE]'
    write(*,*)
    write(*, '(a)') ' OPTIONS:'
    write(*, '(a)') '   [-s | --system [win | linux | mac]]  Operating system; if not provided assumes "win"'
    write(*, '(a)') '   [-m | --mode [embedded | desktop]]   Running mode; if not provided assumes "desktop"'
    write(*, '(a)') '   [-c | --caller [gui | console]]      Caller; if not provided assumes "console"'
    write(*, '(a)') '   [-e | --environment [DIRECTORY]]     Working directory, to be provided in embedded mode;&
                                                             & if not provided assumes \.'
    write(*, '(a)') '   [-h | --help]                        Display this help and exit'
    write(*, '(a)') '   [-v | --version]                     Output version information and exit'
    write(*, '(a)')
    write(*, '(a)') ' PROJ_FILE                              Path of project (*.eddypro) file;&
                                                             & if not provided, assumes ..\ini\processing.eddypro'
    stop
end subroutine CommandLineHelp
