!***************************************************************************
! rename_tmp_files_common.f90
! ---------------------------
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

    !> FLUXNET (fluxes) file
    if (EddyProProj%out_fluxnet_eddy) then
        tmp_indx = index(FLUXNET_EDDY_Path, TmpExt)
        OutPath = FLUXNET_EDDY_Path(1: tmp_indx - 1)
        move_status = system(comm_move // '"' &
            // FLUXNET_EDDY_Path(1:len_trim(FLUXNET_EDDY_Path)) // '" "' &
            // OutPath(1:len_trim(OutPath)) // '"' &
            // comm_out_redirect // comm_err_redirect)
    end if

    !> ICOS file
    if (EddyProProj%out_icos) then
        tmp_indx = index(ICOS_Path, TmpExt)
        OutPath = ICOS_Path(1: tmp_indx - 1)
        move_status = system(comm_move // '"' &
            // ICOS_Path(1:len_trim(ICOS_Path)) // '" "' &
            // OutPath(1:len_trim(OutPath)) // '"' &
            // comm_out_redirect // comm_err_redirect)
    end if

    !> AmeriFlux file
    if (EddyProProj%out_amflux) then
        tmp_indx = index(AmeriFlux_Path, TmpExt)
        OutPath = AmeriFlux_Path(1: tmp_indx - 1)
        move_status = system(comm_move // '"' &
            // AmeriFlux_Path(1:len_trim(AmeriFlux_Path)) // '" "' &
            // OutPath(1:len_trim(OutPath)) // '"' &
            // comm_out_redirect // comm_err_redirect)
    end if
    write(*,'(a)') ' Done.'
end subroutine RenameTmpFilesCommon
