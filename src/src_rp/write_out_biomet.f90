!***************************************************************************
! write_out_biomet.f90
! --------------------
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
! \brief       Write biomet data on (temporary) output files
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WriteOutBiomet(init_string, embedded)
    use m_rp_global_var
    implicit none
    !> in/out variables
    character(*), intent(in) :: init_string
    logical, intent(in) :: embedded
    !> local variables
    integer :: i
    character(LongOutstringLen) :: dataline
    character(DatumLen) :: datum
    character(len=len(init_string)) :: prefix
    include '../src_common/interfaces.inc'

    !>==========================================================================
    !> EddyPro's BIOMET output
    if (embedded) then
        prefix = init_string(index(init_string, ',') + 1: &
                             len_trim(init_string))
    else
        prefix = trim(init_string)
    end if

    if (EddyProProj%out_biomet .and. nbVars > 0) then
        call clearstr(dataline)
        call AddDatum(dataline, trim(adjustl(prefix)), separator)

        do i = 1, nbVars
            call WriteDatumFloat(bAggrEddyPro(i), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        end do
        write(ubiomet, '(a)') dataline(1:len_trim(dataline) - 1)
    end if

end subroutine WriteOutBiomet
