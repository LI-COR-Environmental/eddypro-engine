!***************************************************************************
! fix_timelag_opt_dataset.f90
! ---------------------------
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
! \brief       Eliminate error codes for easier following processing
!              Needs to create a new RH column for each gas
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine FixTimelagOptDataset(TimelagOpt, nrow, toSet, ton, actn, tlncol)
    use m_rp_global_var
    implicit none
    !> In/out variables
    integer, intent(in) :: nrow
    integer, intent(in) :: ton
    integer, intent(in) :: tlncol
    type(TimeLagOptType), intent(in):: TimelagOpt(nrow)
    type(TimeLagDatasetType), intent(out):: toSet(ton)
    integer, intent(out) :: actn(tlncol)
    !> Local variables
    integer :: i
    integer :: gas

    toSet = TimelagDatasetType(0d0, 0d0)
    actn = 0
    do i = 1, ton
        do gas = co2, gas4
            if(gas == h2o .and. TimelagOpt(i)%RH == error) cycle
            if(TimelagOpt(i)%tlag(gas) /= error) then
                actn(gas) = actn(gas) + 1
                toSet(actn(gas))%tlag(gas) = TimelagOpt(i)%tlag(gas)
                if (gas == h2o) toSet(actn(gas))%RH = TimelagOpt(i)%RH
            end if
        end do
    end do
end subroutine FixTimelagOptDataset

