!***************************************************************************
! configure_for_embedded.f90
! --------------------------
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
    character(128) :: sa_fname
    integer :: dir_status
    integer :: io_status
    integer :: ix

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

            ! EddyProProj%out_fluxnet  = .false.
            EddyProProj%out_md      = .false.
            if (EddyProProj%biomet_data /= 'none') then
                EddyProProj%out_biomet = .true.
            else
                EddyProProj%out_biomet = .false.
            end if


        !> EddyPro-FCC
        case ('EddyPro-FCC')
            !> Retrieve FLUXNET file name from /output folder
            comm = 'find "' // trim(homedir) // 'output' // slash // &
                '" -iname *_fluxnet_*' // ' > ' // trim(adjustl(TmpDir)) &
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

                ! Retrieve file name from project file
                ix = index(AuxFile%sa, slash, back=.true.)
                sa_fname = AuxFile%sa(ix+1: len_trim(AuxFile%sa))
                ! File path is $HOME/ini/sa_fname
                AuxFile%sa = trim(homedir) // 'ini' // slash // trim(sa_fname)
                
                !> Retrieve spectral assessment file name from /ini folder
                ! comm = 'find "' // trim(homedir) // 'ini' // slash &
                !     // '" -iname *spectral_assessment*' // ' > ' &
                !     // trim(adjustl(TmpDir)) // 'sa_flist.tmp ' &
                !     // comm_err_redirect
                ! dir_status = system(comm)
                ! open(udf, file = trim(adjustl(TmpDir)) &
                !     // 'sa_flist.tmp', iostat = io_status)
                ! AuxFile%sa = 'none'
                ! if (io_status == 0) then
                !     read(udf, '(a128)', iostat = io_status) dataline
                !     if(io_status == 0) then
                !         AuxFile%sa = trim(adjustl(dataline))
                !         call StripFileName(AuxFile%sa)
                !     end if
                ! end if
                ! close(udf)

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
