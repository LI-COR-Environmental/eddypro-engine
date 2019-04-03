!***************************************************************************
! wind_sector.f90
! ---------------
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
! \brief       Returns wind sector of given wind direction
! \author      Gerardo Fratini
! \note
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WindSector(wdir, n)
    use m_rp_global_var
    implicit none
    !> in/out variables
    real(kind = dbl), intent(in) :: wdir
    integer, intent(out) :: n
    !> local variables
    real(kind = dbl) :: twdir
    integer :: sec


    !> tilt wind direction according to north offset of planar fit
    twdir = wdir - PFSetup%north_offset
    if (twdir >= 360d0) twdir = twdir - 360d0
    if (twdir < 0d0) twdir = 360d0 + twdir

    !> Now treats everything as if offset is zero
    if(twdir >= 0d0 .and. twdir < dfloat(PFSetup%wsect_end(1))) then
        n = 1
    else
        do sec = 2, PFSetup%num_sec
            if(twdir >= PFSetup%wsect_end(sec - 1) .and. twdir < PFSetup%wsect_end(sec)) then
                n = sec
                exit
            end if
        end do
    end if
end subroutine WindSector
