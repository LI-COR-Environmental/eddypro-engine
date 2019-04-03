!***************************************************************************
! init_user_outfiles.f90
! ----------------------
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
! \brief       Initialize output files containing results for insensitive \n
!              variables
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine InitUserOutFiles()
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, external :: CreateDir
    !> local variables
    integer :: mkdir_status = 1      ! initializing to false
    integer :: open_status = 1      ! initializing to false
    integer :: dot
    integer :: i
    character(PathLen) :: Test_Path
    logical :: proceed


    !> create sub-directory
    proceed = .false.
    do i = 1, 7
        if (RPsetup%out_st(i)) then
            proceed = .true.
            exit
        end if
    end do
    if (proceed) then
        UserStatsDir = Dir%main_out(1:len_trim(Dir%main_out)) // SubDirUserStats // slash
        mkdir_status = CreateDir('"' // UserStatsDir(1:len_trim(UserStatsDir)) // '"')
    end if

    !> Statistics files Level 1
    if (RPsetup%out_st(1)) then
        Test_Path = UserStatsDir(1:len_trim(UserStatsDir)) &
                  // EddyProProj%id(1:len_trim(EddyProProj%id)) &
                  // UserStats1_FilePadding // Timestamp_FilePadding // CsvExt
        dot = index(Test_Path, CsvExt, .true.) - 1
        UserSt1_Path = Test_Path(1:dot) // CsvTmpExt
        open(u_user_st1, file = UserSt1_Path, iostat = open_status, encoding = 'utf-8')
        write(u_user_st1, '(a)') 'first_statistics:_on_raw_data'
    end if

    !> Statistics files Level 2
    if (RPsetup%out_st(2)) then
        Test_Path = UserStatsDir(1:len_trim(UserStatsDir)) &
                  // EddyProProj%id(1:len_trim(EddyProProj%id)) &
                  // UserStats2_FilePadding // Timestamp_FilePadding // CsvExt
        dot = index(Test_Path, CsvExt, .true.) - 1
        UserSt2_Path = Test_Path(1:dot) // CsvTmpExt
        open(u_user_st2, file = UserSt2_Path, iostat = open_status, encoding = 'utf-8')
        write(u_user_st2, '(a)') 'second_statistics:_on_raw_data_after_despiking'
    end if

    !> Statistics files Level 3
    if (RPsetup%out_st(3)) then
        Test_Path = UserStatsDir(1:len_trim(UserStatsDir)) &
                  // EddyProProj%id(1:len_trim(EddyProProj%id)) &
                  // UserStats3_FilePadding // Timestamp_FilePadding // CsvExt
        dot = index(Test_Path, CsvExt, .true.) - 1
        UserSt3_Path = Test_Path(1:dot) // CsvTmpExt
        open(u_user_st3, file = UserSt3_Path, iostat = open_status, encoding = 'utf-8')
        write(u_user_st3, '(a)') 'third_statistics:_on_raw_data_after_despiking_and_cross-wind_correction'
    end if

    !> Statistics files Level 4
    if (RPsetup%out_st(4)) then
        Test_Path = UserStatsDir(1:len_trim(UserStatsDir)) &
                  // EddyProProj%id(1:len_trim(EddyProProj%id)) &
                  // UserStats4_FilePadding // Timestamp_FilePadding // CsvExt
        dot = index(Test_Path, CsvExt, .true.) - 1
        UserSt4_Path = Test_Path(1:dot) // CsvTmpExt
        open(u_user_st4, file = UserSt4_Path, iostat = open_status, encoding = 'utf-8')
        write(u_user_st4, '(a)') 'forth_statistics:_on_raw_data_after_despiking_cross_wind_correction&
            &_and_angle-of-attack_correction'
    end if

    !> Statistics files Level 5
    if (RPsetup%out_st(5)) then
        Test_Path = UserStatsDir(1:len_trim(UserStatsDir)) &
                  // EddyProProj%id(1:len_trim(EddyProProj%id)) &
                  // UserStats5_FilePadding // Timestamp_FilePadding // CsvExt
        dot = index(Test_Path, CsvExt, .true.) - 1
        UserSt5_Path = Test_Path(1:dot) // CsvTmpExt
        open(u_user_st5, file = UserSt5_Path, iostat = open_status, encoding = 'utf-8')
        write(u_user_st5, '(a)') 'fifth_statistics:_on_raw_data_after_despiking_cross_wind_correction&
            &_angle-of-attack_correction_and_tilt_correction'
    end if

    !> Statistics files Level 6
    if (RPsetup%out_st(6)) then
        Test_Path = UserStatsDir(1:len_trim(UserStatsDir)) &
                  // EddyProProj%id(1:len_trim(EddyProProj%id)) &
                  // UserStats6_FilePadding // Timestamp_FilePadding // CsvExt
        dot = index(Test_Path, CsvExt, .true.) - 1
        UserSt6_Path = Test_Path(1:dot) // CsvTmpExt
        open(u_user_st6, file = UserSt6_Path, iostat = open_status, encoding = 'utf-8')
        write(u_user_st6, '(a)') 'sixth_statistics:_on_raw_data_after_despiking_cross_wind_correction&
            &_angle-of-attack_correction_tilt_correction_and_time-lag_compensation'
    end if

    !> Statistics files Level 7
    if (RPsetup%out_st(7)) then
        Test_Path = UserStatsDir(1:len_trim(UserStatsDir)) &
                  // EddyProProj%id(1:len_trim(EddyProProj%id)) &
                  // UserStats7_FilePadding // Timestamp_FilePadding // CsvExt
        dot = index(Test_Path, CsvExt, .true.) - 1
        UserSt7_Path = Test_Path(1:dot) // CsvTmpExt
        open(u_user_st7, file = UserSt7_Path, iostat = open_status, encoding = 'utf-8')
        write(u_user_st7, '(a)') 'seventh_statistics:_on_raw_data_after_despiking_cross_wind_correction&
            &_angle-of-attack_correction_tilt_correction_time-lag_compensation_and_detrending'
    end if
end subroutine InitUserOutFiles
