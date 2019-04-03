!***************************************************************************
! init_biomet_out.f90
! -------------------
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
! \brief       Initializes biomet output file
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine InitBiometOut()
    use m_rp_global_var
    use iso_fortran_env
    implicit none

    !> Local variables
    integer :: open_status
    integer :: i
    integer :: dot
    character(PathLen) :: Test_Path
    character(LongOutstringLen) :: header1
    character(LongOutstringLen) :: header2
    character(LongOutstringLen) :: head1_utf8
    character(LongOutstringLen) :: head2_utf8


    !>==========================================================================
    !> EddyPro's biomet output
    if (EddyProProj%out_biomet .and. nbVars > 0) then
        Test_Path = Dir%main_out(1:len_trim(Dir%main_out)) &
                  // EddyProProj%id(1:len_trim(EddyProProj%id)) &
                  // Biomet_FilePadding // Timestamp_FilePadding // CsvExt
        dot = index(Test_Path, CsvExt, .true.) - 1
        Biomet_Path = Test_Path(1:dot) // CsvTmpExt
        open(ubiomet, file = Biomet_Path, iostat = open_status, encoding = 'utf-8')

        !> Initialize string to void
        call Clearstr(header1)
        call Clearstr(header2)
        call Clearstr(head1_utf8)
        call Clearstr(head2_utf8)

        call AddDatum(header1,'date,time,DOY', separator)
        call AddDatum(header2,'[yyyy-mm-dd],[HH:MM],[ddd.ddd]', separator)

        do i = 1, nbVars
            call AddDatum(header1, trim(bVars(i)%label), separator)
            call AddDatum(header2,trim(bVars(i)%pretty_unit_out), separator)
        end do

        call latin1_to_utf8(header1, head1_utf8)
        call latin1_to_utf8(header2, head2_utf8)

        !> Write on output file
        write(ubiomet, '(a)') head1_utf8(1:len_trim(head1_utf8) - 1)
        write(ubiomet, '(a)') head2_utf8(1:len_trim(head2_utf8) - 1)
    end if

end subroutine InitBiometOut
