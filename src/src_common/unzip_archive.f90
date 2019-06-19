!***************************************************************************
! unzip_archive.f90
! -----------------
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
! \brief       Unzip archive, containing at max a status, a metadata and a data file
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine UnZipArchive(ZipFile, MetaExt, DataExt, MetaFile, DataFile, &
    BiometFile, BiometMetaFile, skip_file)
    use m_common_global_var
    implicit none
    !> in/out variables
    character(*), intent(in) :: ZipFile
    character(*), intent(in) :: MetaExt
    character(*), intent(in) :: DataExt
    character(*), intent(out) :: MetaFile
    character(*), intent(out) :: DataFile
    character(*), intent(out) :: BiometFile
    character(*), intent(out) :: BiometMetaFile
    logical, intent(out) :: skip_file
    !> local variables
    integer :: i
    integer :: io_status
    integer :: del_status
    integer :: unzip_status
    integer :: dir_status
    character(CommLen) :: comm
    character(128) :: dataline
    character(128) :: TmpString


    call clearstr(MetaFile)
    call clearstr(DataFile)
    call clearstr(BiometFile)
    call clearstr(BiometMetaFile)

    !> Delete residual files in tmp folder
    comm = trim(comm_del) // ' "' // trim(adjustl(TmpDir)) &
        // '"*.*' // trim(adjustl(DataExt)) &
        // ' ' // comm_err_redirect
    del_status = system(comm)
    comm = trim(comm_del) // ' "' // trim(adjustl(TmpDir)) &
        // '"*.tmp' // ' ' // comm_err_redirect
    del_status = system(comm)

    !> Extract files from archive
    comm = trim(comm_7zip) // ' ' // trim(comm_7zip_x_opt) &
        // ' "' // ZipFile(1:len_trim(ZipFile)) // '" -o"' &
        // trim(adjustl(TmpDir)) // '"' &
        // comm_out_redirect // comm_err_redirect

    unzip_status = system(comm)
    if (unzip_status /= 0) then
        call ExceptionHandler(14)
        skip_file = .true.
        return
    end if
    call clearstr(comm)

    !> Metadata files
    comm = trim(adjustl(comm_dir)) // ' "' &
        // trim(adjustl(TmpDir)) // '"*.' &
        // trim(adjustl(MetaExt))  // &
        ' > "' // trim(adjustl(TmpDir)) // 'meta_flist.tmp" ' &
        // comm_err_redirect
    dir_status = system(comm)

    open(udf, file = trim(adjustl(TmpDir)) // 'meta_flist.tmp', &
        iostat = io_status)
    MetaFile = 'none'
    BiometMetaFile = 'none'
    if (io_status == 0) then
        do i = 1, 2
            read(udf, '(a256)', iostat = io_status) dataline
            if(io_status == 0) then
                if (index(dataline, '-biomet.metadata') /= 0) then
                    BiometMetaFile = dataline(1:len_trim(dataline))
                    call StripFileName(BiometMetaFile)
                else
                    MetaFile = dataline(1:len_trim(dataline))
                    call StripFileName(MetaFile)
                end if
            end if
        end do
    end if
    close(udf, status = 'delete')
    TmpString = MetaFile
    call basename(TmpString, MetaFile, slash)
    TmpString = BiometMetaFile
    call basename(TmpString, BiometMetaFile, slash)

    !> Raw data file
    comm = trim(adjustl(comm_dir)) // ' "' &
        // trim(adjustl(TmpDir)) // '"*.' &
        // trim(adjustl(DataExt))  // &
        ' > "' // trim(adjustl(TmpDir)) // 'data_flist.tmp" ' &
        // comm_err_redirect
    dir_status = system(comm)

    open(udf, file = trim(adjustl(TmpDir)) // 'data_flist.tmp', &
        iostat = io_status)
    DataFile = 'none'
    BiometFile = 'none'
    if (io_status == 0) then
        do i = 1, 2
            read(udf, '(a128)', iostat = io_status) dataline
            if(io_status == 0) then
                if (index(dataline, '-biomet.data') /= 0) then
                    BiometFile = dataline(1:len_trim(dataline))
                    call StripFileName(BiometFile)
                else
                    DataFile = dataline(1:len_trim(dataline))
                    call StripFileName(DataFile)
                end if
            end if
        end do
    end if
    close(udf, status = 'delete')
    TmpString = DataFile
    call basename(TmpString, DataFile, slash)
    TmpString = BiometFile
    call basename(TmpString, BiometFile, slash)
end subroutine UnZipArchive
