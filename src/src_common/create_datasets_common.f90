!***************************************************************************
! create_datasets_common.f90
! --------------------------
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
! \brief       Create continuous datasets from gapped ones, for common output files
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine CreateDatasetsCommon(TimeSeries, nrow, StartIndx, EndIndx)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    integer, intent(in) :: StartIndx
    integer, intent(in) :: EndIndx
    type (DateType), intent(in) :: TimeSeries(nrow)
    !> local variables
    integer :: del_status
    integer :: tmp_indx
    integer :: move_status = 1
    character(PathLen) :: OutPath


    !> Full out file
    if (EddyProProj%out_full) then
        write(*,'(a)', advance = 'no') '  Closing Full Output file..'
        call MakeDataset(FullOut_Path(1:len_trim(FullOut_Path)), &
            TimeSeries, size(TimeSeries), &
            StartIndx, EndIndx, .true., 3)
        write(*,'(a)') ' Done.'
    end if

    !> FLUXNET file - NEVER filled. Only renamed.
    if (EddyProProj%out_fluxnet) then
        write(*,'(a)', advance = 'no') &
            '  Closing FLUXNET output file..'
        tmp_indx = index(FLUXNET_Path, TmpExt)
        OutPath = FLUXNET_Path(1: tmp_indx - 1)
        move_status = system(comm_move // '"' &
            // FLUXNET_Path(1:len_trim(FLUXNET_Path)) // '" "' &
            // OutPath(1:len_trim(OutPath)) // '"' &
            // comm_out_redirect // comm_err_redirect)
            write(*,'(a)') ' Done.'
    end if

    !> Metadata file
    if (EddyProProj%out_md) then
        write(*,'(a)', advance = 'no') &
            '  Closing Metadata file..'
        call MakeDataset(Metadata_Path(1:len_trim(Metadata_Path)), &
            TimeSeries, size(TimeSeries), &
            StartIndx, EndIndx, .true., 1)
        write(*,'(a)') ' Done.'
    end if

    !> Remove temporary output file
    if (len_trim(FullOut_Path) /= 0 .and. EddyProProj%out_full) &
        del_status = system(comm_del // '"' &
        // FullOut_Path(1:len_trim(FullOut_Path)) // '"')

    if (len_trim(Metadata_Path) /= 0 .and. EddyProProj%out_md) &
        del_status = system(comm_del // '"' &
        // Metadata_Path(1:len_trim(Metadata_Path)) // '"')
end subroutine CreateDatasetsCommon
