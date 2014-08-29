!***************************************************************************
! create_datasets_rp.f90
! ----------------------
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
! \brief       Create continuous datasets from gapped ones, only for files
!              that are created solely by RP
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine CreateDatasetsRP(MasterTimeSeries, nrow, rpStartIndx, rpEndIndx)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    integer, intent(in) :: rpStartIndx
    integer, intent(in) :: rpEndIndx
    type (DateType), intent(in) :: MasterTimeSeries(nrow)
    !> local variables
    integer :: del_status
    integer :: tmp_indx
    integer :: move_status = 1
    character(256) :: OutFile


    !> L1 statistics
    if (RPsetup%out_st(1)) then
        write(*,'(a)', advance = 'no') '  Creating Level 1 Statistics dataset..'
        call MakeDataset(St1_Path(1:len_trim(St1_Path)), &
            MasterTimeSeries, size(MasterTimeSeries), rpStartIndx, rpEndIndx, .true., 2)
        write(*,'(a)') ' Done.'
    end if

    !> L2 statistics
    if (RPsetup%out_st(2)) then
        write(*,'(a)', advance = 'no') '  Creating Level 2 Statistics dataset..'
        call MakeDataset(St2_Path(1:len_trim(St2_Path)), &
            MasterTimeSeries, size(MasterTimeSeries), rpStartIndx, rpEndIndx, .true., 2)
        write(*,'(a)') ' Done.'
    end if

    !> L3 statistics
    if (RPsetup%out_st(3)) then
        write(*,'(a)', advance = 'no') '  Creating Level 3 Statistics dataset..'
        call MakeDataset(St3_Path(1:len_trim(St3_Path)), &
            MasterTimeSeries, size(MasterTimeSeries), rpStartIndx, rpEndIndx, .true., 2)
        write(*,'(a)') ' Done.'
    end if

    !> L4 statistics
    if (RPsetup%out_st(4)) then
        write(*,'(a)', advance = 'no') '  Creating Level 4 Statistics dataset..'
        call MakeDataset(St4_Path(1:len_trim(St4_Path)), &
            MasterTimeSeries, size(MasterTimeSeries), rpStartIndx, rpEndIndx, .true., 2)
        write(*,'(a)') ' Done.'
    end if

    !> L5 statistics
    if (RPsetup%out_st(5)) then
        write(*,'(a)', advance = 'no') '  Creating Level 5 Statistics dataset..'
        call MakeDataset(St5_Path(1:len_trim(St5_Path)), &
            MasterTimeSeries, size(MasterTimeSeries), rpStartIndx, rpEndIndx, .true., 2)
        write(*,'(a)') ' Done.'
    end if

    !> L6 statistics
    if (RPsetup%out_st(6)) then
        write(*,'(a)', advance = 'no') '  Creating Level 6 Statistics dataset..'
        call MakeDataset(St6_Path(1:len_trim(St6_Path)), &
            MasterTimeSeries, size(MasterTimeSeries), rpStartIndx, rpEndIndx, .true., 2)
        write(*,'(a)') ' Done.'
    end if

    !> L7 statistics
    if (RPsetup%out_st(7)) then
        write(*,'(a)', advance = 'no') '  Creating Level 7 Statistics dataset..'
        call MakeDataset(St7_Path(1:len_trim(St7_Path)), &
            MasterTimeSeries, size(MasterTimeSeries), rpStartIndx, rpEndIndx, .true., 2)
        write(*,'(a)') ' Done.'
    end if

    if (NumUserVar > 0) then
        !> L1 to L7 user statistics
        if (RPsetup%out_st(1)) then
            write(*,'(a)', advance = 'no') '  Creating Level 1 Statistics dataset for user variables..'
            call MakeDataset(UserSt1_Path(1:len_trim(UserSt1_Path)), &
                MasterTimeSeries, size(MasterTimeSeries), rpStartIndx, rpEndIndx, .true., 2)
            write(*,'(a)') ' Done.'
        end if

        if (RPsetup%out_st(2)) then
            write(*,'(a)', advance = 'no') '  Creating Level 2 Statistics dataset for user variables..'
            call MakeDataset(UserSt2_Path(1:len_trim(UserSt2_Path)), &
                MasterTimeSeries, size(MasterTimeSeries), rpStartIndx, rpEndIndx, .true., 2)
            write(*,'(a)') ' Done.'
        end if

        if (RPsetup%out_st(3)) then
            write(*,'(a)', advance = 'no') '  Creating Level 3 Statistics dataset for user variables..'
            call MakeDataset(UserSt3_Path(1:len_trim(UserSt3_Path)), &
                MasterTimeSeries, size(MasterTimeSeries), rpStartIndx, rpEndIndx, .true., 2)
            write(*,'(a)') ' Done.'
        end if

        if (RPsetup%out_st(4)) then
            write(*,'(a)', advance = 'no') '  Creating Level 4 Statistics dataset for user variables..'
            call MakeDataset(UserSt4_Path(1:len_trim(UserSt4_Path)), &
                MasterTimeSeries, size(MasterTimeSeries), rpStartIndx, rpEndIndx, .true., 2)
            write(*,'(a)') ' Done.'
        end if

        if (RPsetup%out_st(5)) then
            write(*,'(a)', advance = 'no') '  Creating Level 5 Statistics dataset for user variables..'
            call MakeDataset(UserSt5_Path(1:len_trim(UserSt5_Path)), &
                MasterTimeSeries, size(MasterTimeSeries), rpStartIndx, rpEndIndx, .true., 2)
            write(*,'(a)') ' Done.'
        end if

        if (RPsetup%out_st(6)) then
            write(*,'(a)', advance = 'no') '  Creating Level 6 Statistics dataset for user variables..'
            call MakeDataset(UserSt6_Path(1:len_trim(UserSt6_Path)), &
                MasterTimeSeries, size(MasterTimeSeries), rpStartIndx, rpEndIndx, .true., 2)
            write(*,'(a)') ' Done.'
        end if

        if (RPsetup%out_st(7)) then
            write(*,'(a)', advance = 'no') '  Creating Level 7 Statistics dataset for user variables..'
            call MakeDataset(UserSt7_Path(1:len_trim(UserSt7_Path)), &
                MasterTimeSeries, size(MasterTimeSeries), rpStartIndx, rpEndIndx, .true., 2)
            write(*,'(a)') ' Done.'
        end if
    end if

    !> Essentials file is not filled (useless waste of time)
    if (EddyProProj%out_essentials) then
        tmp_indx = index(Essentials_Path, TmpExt)
        OutFile = Essentials_Path(1: tmp_indx - 1)
        move_status = system(comm_move // '"' // Essentials_Path(1:len_trim(Essentials_Path)) // '" "' &
                    // OutFile(1:len_trim(OutFile)) // '"' // comm_out_redirect // comm_err_redirect)
    end if


    !> QC file
    if(RPsetup%out_qc_details .and. Meth%qcflag /= 'none') then
        write(*,'(a)', advance = 'no') '  Creating QC details dataset..'
        call MakeDataset(QCdetails_Path(1:len_trim(QCdetails_Path)), &
            MasterTimeSeries, size(MasterTimeSeries), rpStartIndx, rpEndIndx, .true., 3)
        write(*,'(a)') ' Done.'
    end if

    !> Biomet measurements file
    if (EddyProProj%out_biomet .and. NumBiometVar > 0) then
        write(*,'(a)', advance = 'no') '  Creating Biomet dataset..'
        call MakeDataset(Slow_Path(1:len_trim(Slow_Path)), &
            MasterTimeSeries, size(MasterTimeSeries), rpStartIndx, rpEndIndx, .false., 2)
        write(*,'(a)') ' Done.'
    end if

    !> Remove temporary output file
    if (len_trim(QCdetails_Path) /= 0 .and. RPsetup%out_qc_details .and. Meth%qcflag /= 'none') &
        del_status = system(comm_del // '"' // QCdetails_Path(1:len_trim(QCdetails_Path)) // '"')

    if (len_trim(Slow_Path) /= 0 .and. EddyProProj%out_biomet .and. NumBiometVar > 0) &
        del_status = system(comm_del // '"' // Slow_Path(1:len_trim(Slow_Path)) // '"')

    if (len_trim(St1_Path) /= 0 .and. RPsetup%out_st(1)) &
        del_status = system(comm_del // '"' // St1_Path(1:len_trim(St1_Path)) // '"')
    if (len_trim(St2_Path) /= 0 .and. RPsetup%out_st(2)) &
        del_status = system(comm_del // '"' // St2_Path(1:len_trim(St2_Path)) // '"')
    if (len_trim(St3_Path) /= 0 .and. RPsetup%out_st(3)) &
        del_status = system(comm_del // '"' // St3_Path(1:len_trim(St3_Path)) // '"')
    if (len_trim(St4_Path) /= 0 .and. RPsetup%out_st(4)) &
        del_status = system(comm_del // '"' // St4_Path(1:len_trim(St4_Path)) // '"')
    if (len_trim(St5_Path) /= 0 .and. RPsetup%out_st(5)) &
        del_status = system(comm_del // '"' // St5_Path(1:len_trim(St5_Path)) // '"')
    if (len_trim(St6_Path) /= 0 .and. RPsetup%out_st(6)) &
        del_status = system(comm_del // '"' // St6_Path(1:len_trim(St6_Path)) // '"')
    if (len_trim(St7_Path) /= 0 .and. RPsetup%out_st(7)) &
        del_status = system(comm_del // '"' // St7_Path(1:len_trim(St7_Path)) // '"')

    if (NumUserVar > 0) then
        if (len_trim(UserSt1_Path) /= 0 .and. RPsetup%out_st(1)) &
            del_status = system(comm_del // '"' // UserSt1_Path(1:len_trim(UserSt1_Path)) // '"')
        if (len_trim(UserSt2_Path) /= 0 .and. RPsetup%out_st(2)) &
            del_status = system(comm_del // '"' // UserSt2_Path(1:len_trim(UserSt2_Path)) // '"')
        if (len_trim(UserSt3_Path) /= 0 .and. RPsetup%out_st(3)) &
            del_status = system(comm_del // '"' // UserSt3_Path(1:len_trim(UserSt3_Path)) // '"')
        if (len_trim(UserSt4_Path) /= 0 .and. RPsetup%out_st(4)) &
            del_status = system(comm_del // '"' // UserSt4_Path(1:len_trim(UserSt4_Path)) // '"')
        if (len_trim(UserSt5_Path) /= 0 .and. RPsetup%out_st(5)) &
            del_status = system(comm_del // '"' // UserSt5_Path(1:len_trim(UserSt5_Path)) // '"')
        if (len_trim(UserSt6_Path) /= 0 .and. RPsetup%out_st(6)) &
            del_status = system(comm_del // '"' // UserSt6_Path(1:len_trim(UserSt6_Path)) // '"')
        if (len_trim(UserSt7_Path) /= 0 .and. RPsetup%out_st(7)) &
            del_status = system(comm_del // '"' // UserSt7_Path(1:len_trim(UserSt7_Path)) // '"')
    end if
end subroutine CreateDatasetsRP
