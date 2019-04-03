!***************************************************************************
! statistical_screening.f90
! -------------------------
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
! \brief       Driver to statistical screening according to
!              Vicker and Mahrt (1997, JAOT)
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine StatisticalScreening(Set, nrow, ncol, Tests, printout)
    use m_rp_global_var
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    logical, intent(in) :: printout
    type(TestType), intent(in) :: Tests
    real(kind = dbl), intent(inout) :: Set(nrow, ncol)

    if (printout) write(*, '(a)') '  Raw level statistical screening..'

    !> Absolute limits (al)
    if (Tests%al) call TestAbsoluteLimits(Set, nrow, printout)

    !> Spike count/removal (sr)
    if (Tests%sr) then
        if (RPSetup%despike_vickers97) then
            call TestSpikeDetectionVickers97(Set, nrow, printout)
        else
            call TestSpikeDetectionMauder13(Set, nrow, printout)
        end if
    end if

    !> Amplitude resolution and dropouts (ar, do)
    if (Tests%ar .or. Tests%do) call TestAmpResDropOut(Set, nrow)

    !> Skewness and kurtosis (sk)
    if (Tests%sk) call TestHigherMoments(Set, nrow)

    !> Discontinuities (ds)
    if (Tests%ds) call TestDiscontinuities(Set, nrow)

    !> Time lag (tl)
    if (Tests%tl) call TestTimeLag(Set, nrow)

    !> Angle of attack (aa)
    if (Tests%aa) call TestAttackAngle(Set, nrow)

    !> Non-steady horizontal wind speed (ns)
    if (Tests%ns) call TestNonSteadyWind(Set, nrow)

    !> Set flags to 9 for tests not performed
    call TestsNotPerformed()

    call Int2Flags(9)
    if (printout) write(*,'(a)') '  Done.'
end subroutine StatisticalScreening
