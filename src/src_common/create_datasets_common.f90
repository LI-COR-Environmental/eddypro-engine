!***************************************************************************
! create_datasets_common.f90
! --------------------------
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
! \brief       Create continuous datasets from gapped ones, for common output files
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine CreateDatasetsCommon(MasterTimeSeries, nrow, rpStartIndx, rpEndIndx)
    use m_common_global_var
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


    !> Full out file
    if (EddyProProj%out_full) then
        write(*,'(a)', advance = 'no') '  Creating Full Output dataset..'
        call MakeDataset(FullOut_Path(1:len_trim(FullOut_Path)), &
            MasterTimeSeries, size(MasterTimeSeries), rpStartIndx, rpEndIndx, .true., 3)
        write(*,'(a)') ' Done.'
    end if

    !> GHG-EUROPE file - it is NEVER filled. Only renamed.
    if (EddyProProj%out_ghg_eu) then
        tmp_indx = index(GHGEUROPE_Path, TmpExt)
        OutFile = GHGEUROPE_Path(1: tmp_indx - 1)
        move_status = system(comm_move // '"' // GHGEUROPE_Path(1:len_trim(GHGEUROPE_Path)) // '" "' &
            // OutFile(1:len_trim(OutFile)) // '"' // comm_out_redirect // comm_err_redirect)
!        write(*,'(a)', advance = 'no') '  Creating GHG-EUROPE-style dataset..'
!        call MakeDataset(GHGEUROPE_Path(1:len_trim(GHGEUROPE_Path)), &
!            MasterTimeSeries, size(MasterTimeSeries), rpStartIndx, rpEndIndx, .true., 3)
!        write(*,'(a)') ' Done.'
    end if

    !> AmeriFlux file
    if (EddyProProj%out_amflux) then
        !> Ameriflux_Path
        tmp_indx = index(Ameriflux_Path, TmpExt)
        OutFile = Ameriflux_Path(1: tmp_indx - 1)
        move_status = system(comm_move // '"' // Ameriflux_Path(1:len_trim(Ameriflux_Path)) // '" "' &
                    // OutFile(1:len_trim(OutFile)) // '"' // comm_out_redirect // comm_err_redirect)
    end if

    !> Metadata file
    if (EddyProProj%out_md) then
        write(*,'(a)', advance = 'no') '  Creating Metadata dataset..'
        call MakeDataset(Metadata_Path(1:len_trim(Metadata_Path)), &
            MasterTimeSeries, size(MasterTimeSeries), rpStartIndx, rpEndIndx, .true., 1)
        write(*,'(a)') ' Done.'
    end if

    !> Remove temporary output file
    if (len_trim(FullOut_Path) /= 0 .and. EddyProProj%out_full) &
        del_status = system(comm_del // '"' // FullOut_Path(1:len_trim(FullOut_Path)) // '"')

!    if (len_trim(GHGEUROPE_Path) /= 0 .and. EddyProProj%out_ghg_eu) &
!        del_status = system(comm_del // '"' // GHGEUROPE_Path(1:len_trim(GHGEUROPE_Path)) // '"')

    if (len_trim(Metadata_Path) /= 0 .and. EddyProProj%out_md) &
        del_status = system(comm_del // '"' // Metadata_Path(1:len_trim(Metadata_Path)) // '"')
end subroutine CreateDatasetsCommon
