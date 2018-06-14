!***************************************************************************
! create_datasets_common.f90
! --------------------------
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
        write(*,'(a)', advance = 'no') '  Creating Full Output dataset..'
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
        write(*,'(a)', advance = 'no') '  Creating Metadata dataset..'
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
