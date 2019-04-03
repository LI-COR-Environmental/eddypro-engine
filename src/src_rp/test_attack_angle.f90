!***************************************************************************
! test_attack_angle.f90
! ---------------------
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
! \brief       Checks for too large wind attack angles \n
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine TestAttackAngle(Set, N)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    real(kind = dbl), intent(inout) :: Set(N, E2NumVar)
    !> local variables
    integer :: i = 0
    integer :: cnt
    integer :: ocnt
    real(kind = dbl) :: aoa = 0.d0
    real(kind = dbl) :: hwind


    write(*, '(a)', advance = 'no') '   Angle of attack test..'

    !> Calculation of angle of attack for each row
    cnt = 0
    ocnt = 0
    do i = 1, N
        if (Set(i, u) == error &
            .or. Set(i, v) == error .or. Set(i, w) == error) cycle
        cnt = cnt + 1
        hwind = dsqrt(Set(i, u)**2 + Set(i, v)**2)
        if(hwind == 0d0 .and. Set(i, w) /= 0d0) then
            ocnt = ocnt + 1
        elseif (hwind /= 0d0) then
            aoa = (atan(Set(i, w) / hwind)) * 180.d0 / p
            if (aoa < aa%min .or. aoa > aa%max) ocnt = ocnt + 1
        end if
    end do
    Essentials%aa_s = dble(ocnt) / dble(cnt) * 1d2
    IntHF%aa = 0
    if (dble(ocnt) / dble(cnt) * 1d2 >= aa%lim) IntHF%aa = 1
    write(*,'(a)') ' Done.'
end subroutine TestAttackAngle
