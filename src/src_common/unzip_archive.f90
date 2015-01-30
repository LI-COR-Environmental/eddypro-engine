!***************************************************************************
! unzip_archive.f90
! -----------------
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
        // trim(adjustl(TmpDir)) // '"'&
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
            read(udf, '(a128)', iostat = io_status) dataline
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

    !> Fast metadata and data file names definition, based on zip file name
    !> Less robust
    !ZipFileNameBody = ZipFile(index(ZipFile, slash, .true.) + 1: len_trim(ZipFile) - 3)
    !MetaFile = ZipFileNameBody(1:len_trim(ZipFileNameBody)) // MetaExt
    !DataFile = ZipFileNameBody(1:len_trim(ZipFileNameBody))  // DataExt
    !BiometFile = ZipFileNameBody(1:len_trim(ZipFileNameBody))  // BiometExt
    !BiometMetaFile = ZipFileNameBody(1:len_trim(ZipFileNameBody))  // BiometMetaExt
end subroutine UnZipArchive
