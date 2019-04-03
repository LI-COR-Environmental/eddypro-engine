!***************************************************************************
! user_fluctuations.f90
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
! \brief       Calculate fluctuations around a trend, defined by the \n
!              chosen detrending method
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine UserFluctuations(Set, Primes, nrow, ncol, Tconst, LocStats, LocCol)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    integer, intent(in) :: Tconst
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    type(UserStatsType), intent(in) :: LocStats
    type(ColType), intent(in) :: LocCol(ncol)
    real(kind = dbl), intent(out) :: Primes(nrow, ncol)
    !> local variables
    integer :: var


    select case(Meth%det(1:len_trim(Meth%det)))
        case('ld')
            call LinDetrend(Set, Primes, Tconst, LocCol, nrow, ncol)
        case('rm')
            call RunningMean(Set, Primes, Tconst, LocCol, nrow, ncol)
        case('ew')
            call ExpWeightAvrg(Set, Primes, Tconst, LocCol, nrow, ncol)
        case('ba')
            do var = 1, ncol
                if (LocCol(var)%present) then
                    where(Set(1:nrow, var) /= error)
                        Primes(1:nrow, var) = Set(1:nrow, var) - LocStats%Mean(var)
                    elsewhere
                        Primes(1:nrow, var) = error
                    end where
                end if
            end do
    end select
end subroutine UserFluctuations
