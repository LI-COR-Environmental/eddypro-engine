!***************************************************************************
! rename_tmp_files_common.f90
! ---------------------------
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
! \brief       Remove tmp extension from temporary file names
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine RenameTmpFilesCommon()
    use m_common_global_var
    implicit none
    !> local variables
    integer :: tmp_indx
    integer :: move_status = 1
    character(PathLen) :: OutPath

    write(*,'(a)', advance = 'no') ' Closing COMMON output files..'

    !> Full out file
    if (EddyProProj%out_full) then
        tmp_indx = index(FullOut_Path, TmpExt)
        OutPath = FullOut_Path(1: tmp_indx - 1)
        move_status = system(comm_move // '"' &
        // FullOut_Path(1:len_trim(FullOut_Path)) // '" "' &
            // OutPath(1:len_trim(OutPath)) // '"' &
            // comm_out_redirect // comm_err_redirect)
    end if

    !> Metadata
    if (EddyProProj%out_md) then
        tmp_indx = index(Metadata_Path, TmpExt)
        OutPath = Metadata_Path(1: tmp_indx - 1)
        move_status = system(comm_move // '"' &
            // Metadata_Path(1:len_trim(Metadata_Path)) // '" "' &
            // OutPath(1:len_trim(OutPath)) // '"' &
            // comm_out_redirect // comm_err_redirect)
    end if

    !> FLUXNET file
    if (EddyProProj%out_fluxnet) then
        tmp_indx = index(FLUXNET_Path, TmpExt)
        OutPath = FLUXNET_Path(1: tmp_indx - 1)
        move_status = system(comm_move // '"' &
            // FLUXNET_Path(1:len_trim(FLUXNET_Path)) // '" "' &
            // OutPath(1:len_trim(OutPath)) // '"' &
            // comm_out_redirect // comm_err_redirect)
    end if

    write(*,'(a)') ' Done.'
end subroutine RenameTmpFilesCommon
