!***************************************************************************
! statistical_screening.f90
! -------------------------
! Copyright (C) 2007-2011, Eco2s team, Gerardo Fratini
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
! \brief       Driver to statistical screening according to Vicker and Mahrt (1997, JAOT)
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine StatisticalScreening(Set, nrow, ncol, LocTest, printout)
    use m_rp_global_var
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    logical, intent(in) :: printout
    type(TestType), intent(in) :: LocTest
    real(kind = dbl), intent(inout) :: Set(nrow, ncol)

    if (printout) write(*, '(a)') '  Raw level statistical screening..'

    !> Spike count/removal (sr)
    if (LocTest%sr) then
        if (RPSetup%despike_vm) then
            call TestSpikeDetection(Set, nrow, printout)
        else
            call TestSpikeDetectionComingSoon(Set, nrow, printout)
        end if
    end if

    !> Amplitude resolution and dropouts (ar, do)
    if (LocTest%ar .or. LocTest%do) call TestAmpResDropOut(Set, nrow)

    !> Absolute limits (al)
    if (LocTest%al) call TestAbsoluteLimits(Set, nrow, printout)

    !> Skewness and kurtosis (sk)
    if (LocTest%sk) call TestHigherMoments(Set, nrow)

    !> Discontinuities (ds)
    if (LocTest%ds) call TestDiscontinuities(Set, nrow)

    !> Time lag (tl)
    if (LocTest%tl) call TestTimeLag(Set, nrow)

    !> Angle of attack (aa)
    if (LocTest%aa) call TestAttackAngle(Set, nrow)

    !> Non-steady horizontal wind speed (ns)
    if (LocTest%ns) call TestNonSteadyWind(Set, nrow)

    !> Set flags to 9 for tests not performed
    call TestsNotPerformed()

    call Int2Flags(9)
    if (printout) write(*,'(a)') '  Done.'
end subroutine StatisticalScreening
