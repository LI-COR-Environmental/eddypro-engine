!***************************************************************************
! configure_for_embedded.f90
! --------------------------
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
! \brief       Overwrite settings of the .eddypro file, for use \n
!              in embedded mode
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ConfigureForEmbedded()
    use m_common_global_var
    implicit none
    !> Local variables
    character(CommLen) :: comm
    character(ShortInstringLen) :: dataline
    integer :: dir_status
    integer :: io_status

    select case (app)
        !> EddyPro-RP
        case ('EddyPro-RP')
            Dir%main_in  = trim(homedir) // 'raw_files' // slash

            !> Retrieve planar fit file name if needed
            if (index(Meth%rot, 'planar_fit') /= 0) then

                !> Retrieve planar fit file name from /ini folder
                comm = 'find "' // trim(homedir) // 'ini' // slash &
                    // '" -iname *_planar_fit_*' // ' > ' // '"' &
                    // trim(adjustl(TmpDir)) // 'pf_flist.tmp" ' &
                    // comm_err_redirect
                dir_status = system(comm)
                open(udf, file = trim(adjustl(TmpDir)) &
                    // 'pf_flist.tmp', iostat = io_status)
                AuxFile%pf = 'none'
                if (io_status == 0) then
                    read(udf, '(a128)', iostat = io_status) dataline
                    if(io_status == 0) then
                        AuxFile%pf = trim(adjustl(dataline))
                        call StripFileName(AuxFile%pf)
                    end if
                end if
                close(udf)
            end if

            !> Retrieve time-lag optimization file name if needed
            if (index(Meth%tlag, 'tlag_opt') /= 0) then

                !> Retrieve planar fit file name from /ini folder
                comm = 'find "' // trim(homedir) // 'ini' // slash &
                    // '" -iname *_timelag_opt_*'// ' > ' // '"' &
                    // trim(adjustl(TmpDir)) // 'to_flist.tmp" ' &
                    // comm_err_redirect
                dir_status = system(comm)
                open(udf, file = trim(adjustl(TmpDir)) &
                    // 'to_flist.tmp', iostat = io_status)
                AuxFile%to = 'none'
                if (io_status == 0) then
                    read(udf, '(a128)', iostat = io_status) dataline
                    if(io_status == 0) then
                        AuxFile%to = trim(adjustl(dataline))
                        call StripFileName(AuxFile%to)
                    end if
                end if
                close(udf)
            end if

            !> Delete all temporary files
            call system(comm_del // '"' // trim(adjustl(TmpDir)) // '"*.tmp ' &
                // comm_err_redirect)

            !> Selection of output files
            if (EddyProProj%fcc_follows) then
                EddyProProj%out_essentials = .true.
            else
                EddyProProj%out_essentials = .false.
            end if
            EddyProProj%out_fluxnet  = .false.
            EddyProProj%out_md      = .false.
            if (EddyProProj%biomet_data /= 'none') then
                EddyProProj%out_biomet = .true.
            else
                EddyProProj%out_biomet = .false.
            end if


        !> EddyPro-FCC
        case ('EddyPro-FCC')
            !> Retrieve essentials file name from /output folder
            comm = 'find "' // trim(homedir) // 'output' // slash // &
                '" -iname *_essentials_*' // ' > ' // trim(adjustl(TmpDir)) &
                // 'ex_flist.tmp ' // comm_err_redirect
            dir_status = system(comm)
            open(udf, file = trim(adjustl(TmpDir)) &
                // 'ex_flist.tmp', iostat = io_status)
            AuxFile%ex = 'none'
            if (io_status == 0) then
                read(udf, '(a128)', iostat = io_status) dataline
                if(io_status == 0) then
                    AuxFile%ex = trim(adjustl(dataline))
                    call StripFileName(AuxFile%ex)
                end if
            end if
            close(udf)

            !> Retrieve spectral assessment file name if needed
            if (EddyProProj%hf_meth =='fratini_12' .or. &
                EddyProProj%hf_meth =='horst_97' .or. &
                EddyProProj%hf_meth =='ibrom_07') then

                !> Retrieve spectral assessment file name from /ini folder
                comm = 'find "' // trim(homedir) // 'ini' // slash &
                    // '" -iname *_spectral_assessment_*' // ' > ' &
                    // trim(adjustl(TmpDir)) // 'sa_flist.tmp ' &
                    // comm_err_redirect
                dir_status = system(comm)
                open(udf, file = trim(adjustl(TmpDir)) &
                    // 'sa_flist.tmp', iostat = io_status)
                AuxFile%sa = 'none'
                if (io_status == 0) then
                    read(udf, '(a128)', iostat = io_status) dataline
                    if(io_status == 0) then
                        AuxFile%sa = trim(adjustl(dataline))
                        call StripFileName(AuxFile%sa)
                    end if
                end if
                close(udf)
            end if

            !> Delet all temporary files
            call system(comm_del // '"' // trim(adjustl(TmpDir)) &
                // '"*.tmp ' // comm_err_redirect)

            !> Selection of output files
            EddyProProj%out_fluxnet  = .false.
    end select

    !> Common settings
    Dir%main_out = trim(homedir) // 'output' // slash

end subroutine ConfigureForEmbedded
