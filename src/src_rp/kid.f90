!***************************************************************************
! kid.f90
! -------
! Copyright (C) 2018-2019, LI-COR Biosciences, Inc.  All Rights Reserved.
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
! \brief       Compute Kurtosis index on differenced variables
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine KID(Set, nrow, ncol)

    use m_rp_global_var
    use mo_timelag_handle, only: VariableStochasticDetrending

    implicit none

    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    !> Local variables
    integer :: var
    real(kind = dbl) :: Primes(nrow, ncol)

    do var = u, ts
        call VariableStochasticDetrending(Set(:, var), Primes(:, var), nrow)
        call KurtosisNoError(Primes(:, var), nrow, 1, Essentials%KID(var), error)
        Essentials%ZCD(var) = count(abs(Primes(:, var)) < 1d-6)
    end do

    do var = co2, gas4
        call VariableStochasticDetrending(Set(:, var), Primes(:, var), nrow)
        call KurtosisNoError(Primes(:, var), nrow, 1, Essentials%KID(var), error)
        if (E2Col(var)%present) then
            Essentials%ZCD(var) = count(abs(Primes(:, var)) < 1d-6)
        else
            Essentials%ZCD(var) = ierror
        end if
    end do

end subroutine KID
