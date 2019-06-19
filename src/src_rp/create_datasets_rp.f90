!***************************************************************************
! create_datasets_rp.f90
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
subroutine CreateDatasetsRP(TimeSeries, nrow, StartIndx, EndIndx)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    integer, intent(in) :: StartIndx
    integer, intent(in) :: EndIndx
    type (DateType), intent(in) :: TimeSeries(nrow)
    !> local variables
    integer :: del_status


    !> L1 statistics
    if (RPsetup%out_st(1)) then
        write(*,'(a)', advance = 'no') '  Creating Level 1 Statistics dataset..'
        call MakeDataset(St1_Path(1:len_trim(St1_Path)), &
            TimeSeries, size(TimeSeries), &
            StartIndx, EndIndx, .true., 2)
        write(*,'(a)') ' Done.'
    end if

    !> L2 statistics
    if (RPsetup%out_st(2)) then
        write(*,'(a)', advance = 'no') '  Creating Level 2 Statistics dataset..'
        call MakeDataset(St2_Path(1:len_trim(St2_Path)), &
            TimeSeries, size(TimeSeries), &
            StartIndx, EndIndx, .true., 2)
        write(*,'(a)') ' Done.'
    end if

    !> L3 statistics
    if (RPsetup%out_st(3)) then
        write(*,'(a)', advance = 'no') '  Creating Level 3 Statistics dataset..'
        call MakeDataset(St3_Path(1:len_trim(St3_Path)), &
            TimeSeries, size(TimeSeries), &
            StartIndx, EndIndx, .true., 2)
        write(*,'(a)') ' Done.'
    end if

    !> L4 statistics
    if (RPsetup%out_st(4)) then
        write(*,'(a)', advance = 'no') '  Creating Level 4 Statistics dataset..'
        call MakeDataset(St4_Path(1:len_trim(St4_Path)), &
            TimeSeries, size(TimeSeries), &
            StartIndx, EndIndx, .true., 2)
        write(*,'(a)') ' Done.'
    end if

    !> L5 statistics
    if (RPsetup%out_st(5)) then
        write(*,'(a)', advance = 'no') '  Creating Level 5 Statistics dataset..'
        call MakeDataset(St5_Path(1:len_trim(St5_Path)), &
            TimeSeries, size(TimeSeries), &
            StartIndx, EndIndx, .true., 2)
        write(*,'(a)') ' Done.'
    end if

    !> L6 statistics
    if (RPsetup%out_st(6)) then
        write(*,'(a)', advance = 'no') '  Creating Level 6 Statistics dataset..'
        call MakeDataset(St6_Path(1:len_trim(St6_Path)), &
            TimeSeries, size(TimeSeries), &
            StartIndx, EndIndx, .true., 2)
        write(*,'(a)') ' Done.'
    end if

    !> L7 statistics
    if (RPsetup%out_st(7)) then
        write(*,'(a)', advance = 'no') '  Creating Level 7 Statistics dataset..'
        call MakeDataset(St7_Path(1:len_trim(St7_Path)), &
            TimeSeries, size(TimeSeries), &
            StartIndx, EndIndx, .true., 2)
        write(*,'(a)') ' Done.'
    end if

    if (NumUserVar > 0) then
        !> L1 to L7 user statistics
        if (RPsetup%out_st(1)) then
            write(*,'(a)', advance = 'no') &
                '  Creating Level 1 Statistics dataset for user variables..'
            call MakeDataset(UserSt1_Path(1:len_trim(UserSt1_Path)), &
                TimeSeries, size(TimeSeries), &
                StartIndx, EndIndx, .true., 2)
            write(*,'(a)') ' Done.'
        end if

        if (RPsetup%out_st(2)) then
            write(*,'(a)', advance = 'no') &
                '  Creating Level 2 Statistics dataset for user variables..'
            call MakeDataset(UserSt2_Path(1:len_trim(UserSt2_Path)), &
                TimeSeries, size(TimeSeries), &
                StartIndx, EndIndx, .true., 2)
            write(*,'(a)') ' Done.'
        end if

        if (RPsetup%out_st(3)) then
            write(*,'(a)', advance = 'no') &
                '  Creating Level 3 Statistics dataset for user variables..'
            call MakeDataset(UserSt3_Path(1:len_trim(UserSt3_Path)), &
                TimeSeries, size(TimeSeries), &
                StartIndx, EndIndx, .true., 2)
            write(*,'(a)') ' Done.'
        end if

        if (RPsetup%out_st(4)) then
            write(*,'(a)', advance = 'no') &
                '  Creating Level 4 Statistics dataset for user variables..'
            call MakeDataset(UserSt4_Path(1:len_trim(UserSt4_Path)), &
                TimeSeries, size(TimeSeries), &
                StartIndx, EndIndx, .true., 2)
            write(*,'(a)') ' Done.'
        end if

        if (RPsetup%out_st(5)) then
            write(*,'(a)', advance = 'no') &
                '  Creating Level 5 Statistics dataset for user variables..'
            call MakeDataset(UserSt5_Path(1:len_trim(UserSt5_Path)), &
                TimeSeries, size(TimeSeries), &
                StartIndx, EndIndx, .true., 2)
            write(*,'(a)') ' Done.'
        end if

        if (RPsetup%out_st(6)) then
            write(*,'(a)', advance = 'no') &
                '  Creating Level 6 Statistics dataset for user variables..'
            call MakeDataset(UserSt6_Path(1:len_trim(UserSt6_Path)), &
                TimeSeries, size(TimeSeries), &
                StartIndx, EndIndx, .true., 2)
            write(*,'(a)') ' Done.'
        end if

        if (RPsetup%out_st(7)) then
            write(*,'(a)', advance = 'no') &
                '  Creating Level 7 Statistics dataset for user variables..'
            call MakeDataset(UserSt7_Path(1:len_trim(UserSt7_Path)), &
                TimeSeries, size(TimeSeries), &
                StartIndx, EndIndx, .true., 2)
            write(*,'(a)') ' Done.'
        end if
    end if

    !> QC file
    if(RPsetup%out_qc_details .and. Meth%qcflag /= 'none') then
        write(*,'(a)', advance = 'no') '  Creating QC details dataset..'
        call MakeDataset(QCdetails_Path(1:len_trim(QCdetails_Path)), &
            TimeSeries, size(TimeSeries), &
            StartIndx, EndIndx, .true., 3)
        write(*,'(a)') ' Done.'
    end if

    !> Biomet measurements file
    if (EddyProProj%out_biomet .and. nbVars > 0) then
        write(*,'(a)', advance = 'no') '  Creating Biomet dataset..'
        call MakeDataset(Biomet_Path(1:len_trim(Biomet_Path)), &
            TimeSeries, size(TimeSeries), &
            StartIndx, EndIndx, .false., 2)
        write(*,'(a)') ' Done.'
    end if

    !> Remove temporary output file
    if (len_trim(QCdetails_Path) /= 0 &
        .and. RPsetup%out_qc_details .and. Meth%qcflag /= 'none') &
        del_status = system(comm_del // '"' &
        // QCdetails_Path(1:len_trim(QCdetails_Path)) // '"')

    if (len_trim(Biomet_Path) /= 0 &
        .and. EddyProProj%out_biomet .and. nbVars > 0) &
        del_status = system(comm_del // '"' &
        // Biomet_Path(1:len_trim(Biomet_Path)) // '"')

    if (len_trim(St1_Path) /= 0 .and. RPsetup%out_st(1)) &
        del_status = system(comm_del // '"' &
        // St1_Path(1:len_trim(St1_Path)) // '"')
    if (len_trim(St2_Path) /= 0 .and. RPsetup%out_st(2)) &
        del_status = system(comm_del // '"' &
        // St2_Path(1:len_trim(St2_Path)) // '"')
    if (len_trim(St3_Path) /= 0 .and. RPsetup%out_st(3)) &
        del_status = system(comm_del // '"' &
        // St3_Path(1:len_trim(St3_Path)) // '"')
    if (len_trim(St4_Path) /= 0 .and. RPsetup%out_st(4)) &
        del_status = system(comm_del // '"' &
        // St4_Path(1:len_trim(St4_Path)) // '"')
    if (len_trim(St5_Path) /= 0 .and. RPsetup%out_st(5)) &
        del_status = system(comm_del // '"' &
        // St5_Path(1:len_trim(St5_Path)) // '"')
    if (len_trim(St6_Path) /= 0 .and. RPsetup%out_st(6)) &
        del_status = system(comm_del // '"' &
        // St6_Path(1:len_trim(St6_Path)) // '"')
    if (len_trim(St7_Path) /= 0 .and. RPsetup%out_st(7)) &
        del_status = system(comm_del // '"' &
        // St7_Path(1:len_trim(St7_Path)) // '"')

    if (NumUserVar > 0) then
        if (len_trim(UserSt1_Path) /= 0 .and. RPsetup%out_st(1)) &
            del_status = system(comm_del // '"' &
            // UserSt1_Path(1:len_trim(UserSt1_Path)) // '"')
        if (len_trim(UserSt2_Path) /= 0 .and. RPsetup%out_st(2)) &
            del_status = system(comm_del // '"' &
            // UserSt2_Path(1:len_trim(UserSt2_Path)) // '"')
        if (len_trim(UserSt3_Path) /= 0 .and. RPsetup%out_st(3)) &
            del_status = system(comm_del // '"' &
            // UserSt3_Path(1:len_trim(UserSt3_Path)) // '"')
        if (len_trim(UserSt4_Path) /= 0 .and. RPsetup%out_st(4)) &
            del_status = system(comm_del // '"' &
            // UserSt4_Path(1:len_trim(UserSt4_Path)) // '"')
        if (len_trim(UserSt5_Path) /= 0 .and. RPsetup%out_st(5)) &
            del_status = system(comm_del // '"' &
            // UserSt5_Path(1:len_trim(UserSt5_Path)) // '"')
        if (len_trim(UserSt6_Path) /= 0 .and. RPsetup%out_st(6)) &
            del_status = system(comm_del // '"' &
            // UserSt6_Path(1:len_trim(UserSt6_Path)) // '"')
        if (len_trim(UserSt7_Path) /= 0 .and. RPsetup%out_st(7)) &
            del_status = system(comm_del // '"' &
            // UserSt7_Path(1:len_trim(UserSt7_Path)) // '"')
    end if
end subroutine CreateDatasetsRP
