!***************************************************************************
! read_ext_biomet_files.f90
! -------------------------
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
! \brief       Interface to reading biomet file, to handle multiple files and \n
!              the file queue. \n
! \author      Gerardo Fratini
! \sa
! \bug
! \deprecated
! \test
!***************************************************************************
subroutine ReadExtBiometFiles(BiometFileList, NumBiometFiles, &
        last_nfl, last_nrec, CurrentTimestamp, BiometDataExist, nbiomet, printout)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: NumBiometFiles
    type(FileListType), intent(in) :: BiometFileList(NumBiometFiles)
    type(DateType), intent(in) :: CurrentTimestamp
    integer, intent(out) :: nbiomet
    logical, intent(out) :: BiometDataExist
    integer, intent(inout) :: last_nfl
    integer, intent(inout) :: last_nrec
    logical, intent(in) :: printout
    !> local variables
    integer :: i
    integer :: j
    integer :: ii
    integer :: jj
    integer :: iii
    integer :: nfl
    integer :: nrec
    integer :: open_status
    integer :: read_status
    integer :: var_num
    integer :: sepa
    integer :: ncstm
    character(32) :: text_vars(1000)
    character(1024) :: datastring
    character(64) :: tstamp_string
    type (DateType) :: tol
    type (DateType) :: win
    type (DateType) :: biomet_date
    logical :: BiometPeriodHooked
    logical :: skip_init


    !> Initializations
    skip_init = .true.

    !> Define time periods for which biomet data are needed
    !> MOVE TO A GLOBAL NAMESPACE (add to bFileMetadata and probably
    !> set in main)
    win = datetype(0, 0, 0, 0, RPsetup%avrg_len)
    tol = datetype(0, 0, 0, 0, bFileMetadata%time_step / 2)

    !> Initialization
    i = 0
    nbiomet = 0
    BiometPeriodHooked = .false.
    file_loop: do nfl = last_nfl, NumBiometFiles
        nrec = 0

        !> Open biomet measurement file(s) and read data,
        !> selecting those in plausible ranges
        if (printout) then
            if (nfl == last_nfl) write(*, '(a)') &
                '  Searching biomet data in file: '
            write(*, '(a)') '   ' &
                // trim(adjustl(BiometFileList(nfl)%path))
        end if

        call ReadBiometFile(BiometFileList(nfl)%path, last_nrec, skip_init, &

            skip_file)




                    !> Retrieve timestamp from tstamp_string
                    call BiometTimestamp(trim(adjustl(bFileMetadata%tsPattern)), &
                        biomet_date)


                    !else
                        if (BiometPeriodHooked) then
                            last_nrec = nrec - 1
                            close(udf)
                            exit file_loop
                        else
                            i = i - 1
                            if (biomet_date > CurrentTimestamp + tol) then
                                last_nrec = nrec - 1
                                close(udf)
                                exit file_loop
                            else
                                cycle rec_loop
                            end if
                        end if
                    end if
                end do ol
            end do rec_loop
        else
            call ExceptionHandler(2)
        end if
    end do file_loop
    nbiomet = i - 1

    if (nbiomet < 1) then
        BiometDataExist = .false.
        if (printout) write(*,'(a)') '  No valid biomet records found for &
            &this averaging period. Continuing without biomet data.'
    else
        if (printout) then
            write(LogInteger, '(i6)') nbiomet
            write(*,'(a)') '  ' // trim(adjustl(LogInteger)) &
                // ' biomet record(s) imported correctly for this averaging period.'
        end if
    end if

    !> Adjust units as needed
    call BiometStandardUnits(nbiomet)

    !> Adjust timesteps of Biomet data if needed
    call AdjustBiometTimestamps(nbiomet)
end subroutine ReadExtBiometFiles


!***************************************************************************
!
! \brief       Retrieve time step intrinsic in Biomet file (in minutes)
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine AdjustBiometTimestamps(N)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    !> local variables
    integer :: i
    type(datetype) :: BioometTimestamp

    select case (BiometSetup%tstamp_ref)
        case ('begin')
            do i = 1, N
                call DateTimeToDateType(Biomet(i)%date, Biomet(i)%time, BioometTimestamp)
                BioometTimestamp = BioometTimestamp + datetype(0, 0, 0, 0, BiometSetup%tstep)
                call DateTypeToDateTime(BioometTimestamp, Biomet(i)%date, Biomet(i)%time)
            end do

        case ('middle')
            do i = 1, N
                call DateTimeToDateType(Biomet(i)%date, Biomet(i)%time, BioometTimestamp)
                BioometTimestamp = BioometTimestamp + datetype(0, 0, 0, 0, BiometSetup%tstep/2)
                call DateTypeToDateTime(BioometTimestamp, Biomet(i)%date, Biomet(i)%time)
            end do
        case ('end')
        return
    end select
end subroutine AdjustBiometTimestamps
