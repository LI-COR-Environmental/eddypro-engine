!***************************************************************************
! cleanup_e2_set.f90
! ------------------
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
! \brief       Set implausible values to error code in EddyPro main dataset \n
!              variables
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine CleanUpE2Set(Set, nrow, ncol)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: ncol
    integer, intent(in) :: nrow
    real(kind = dbl), intent(inout) :: Set(nrow, ncol)


    !> General
    where(Set(:, :) < -300d0)
        Set(:, :) = error
    end where
    !> Specific for wind components
    where(Set(:, u:w) < -50d0 .or. Set(:, u:w) > 50d0)
        Set(:, u:w) = error
    end where
    !> Specific for temperatures
    where(Set(:, ts) < 200d0 .or. Set(:, ts) > 350d0)
        Set(:, ts) = error
    end where
    where(Set(:, ti1) < 200d0 .or. Set(:, ti1) > 350d0)
        Set(:, ti1) = error
    end where
    where(Set(:, ti2) < 200d0 .or. Set(:, ti2) > 350d0)
        Set(:, ti2) = error
    end where
    where(Set(:, te) < 200d0 .or. Set(:, te) > 350d0)
        Set(:, te) = error
    end where
    !> Specific for concentrations
    !> TBD
end subroutine CleanUpE2Set
