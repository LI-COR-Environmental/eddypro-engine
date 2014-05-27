!***************************************************************************
! init_user_outfiles.f90
! ----------------------
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
    character(256) :: Test_Path
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
            &_angle-of-attack_correction_and_double_rotation'
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
            &_angle-of-attack_correction_double_rotation_and_time-lag_compensation'
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
            &_angle-of-attack_correction_double_rotation_time-lag_compensation_and_detrending'
    end if
end subroutine InitUserOutFiles
