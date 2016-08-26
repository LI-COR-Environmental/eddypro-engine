!***************************************************************************
! user_timelag_handle.f90
! -----------------------
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
! \brief       Calculates time-lags (in terms of number of records) \n
!              for all non-sensitive variables
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine UserTimeLagHandle(TlagMeth, UserSet, unrow, uncol, E2W, nrow)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: unrow, uncol
    integer, intent(in) :: nrow
    real(kind = dbl), intent(in) :: E2W(nrow)
    character(*), intent(in) :: TlagMeth
    real(kind = dbl), intent(inout) :: UserSet(unrow, uncol)
    !> local variables
    integer :: i = 0
    integer :: j = 0
    integer :: min_rl(uncol)
    integer :: max_rl(uncol)
    integer :: def_rl(uncol)
    real(kind = dbl) :: FirstCol(unrow)
    real(kind = dbl) :: SecondCol(unrow)
    real(kind = dbl) :: TLag(uncol)
    real(kind = dbl) :: TmpUserSet(unrow, uncol)
    logical :: DefTlagUsed(uncol)

    write(*, '(a)', advance = 'no') '  Compensating user variables time-lags..'

    !> for User scalars, initialise auxiliary vars to zero
    def_rl(:) = 0
    min_rl(:) = 0
    max_rl(:) = 0
    !> Define "row-lags" for scalars, using timelags retrieved from metadata file
    where (UserCol(1:uncol)%var /= 'none')
        def_rl(1:uncol) = nint(UserCol(1:uncol)%def_tl * Metadata%ac_freq)
        min_rl(1:uncol) = nint(UserCol(1:uncol)%min_tl * Metadata%ac_freq)
        max_rl(1:uncol) = nint(UserCol(1:uncol)%max_tl * Metadata%ac_freq)
    endwhere

    !> Calculate actual time-lags according to the chosen method
    select case(Meth%tlag(1:len_trim(Meth%tlag)))
        case ('constant')
            !> Constant timelags are set equal to default values (user selected)
            UserRowLags(1:uncol) = def_rl(1:uncol)
            TLag(1:uncol) = UserCol(1:uncol)%def_tl
        case ('maxcov', 'maxcov&default')
            !> Covariance maximization method, with or without default
            do j = 1, uncol
                if (UserCol(j)%var /= 'none') then
                    FirstCol(:)  = E2W(:)
                    SecondCol(:) = UserSet(:, j)
                    call CovMax(TlagMeth, def_rl(j), min_rl(j), max_rl(j), &
                        FirstCol, SecondCol, size(FirstCol), TLag(j), UserRowLags(j), DefTlagUsed(j))
                else
                    UserRowLags(j) = 0
                    TLag(j) = 0d0
                end if
            end do
        case ('none')
            !> not compensating for timelags
            UserRowLags(1:uncol) = 0
            TLag(1:uncol) = 0d0
    end select

    !> Align data according to relevant time-lags, filling remaining with error code.
    do j = 1, uncol
        if (UserRowLags(j) >= 0) then
            !> For positive lags
            do i = 1, unrow - UserRowLags(j)
                TmpUserSet(i, j) = UserSet(i + UserRowLags(j), j)
            end do
            do i = unrow - UserRowlags(j) + 1, unrow
                TmpUserSet(i, j) = error
            end do
        else
            !> For negative lags
            do i = 1, abs(UserRowLags(j))
                TmpUserSet(i, j) = error
            end do
            do i = abs(UserRowLags(j)) + 1, unrow
                TmpUserSet(i, j) = UserSet(i + UserRowLags(j), j)
            end do
        end if
    end do
    UserSet = TmpUserSet
    write(*,'(a)') ' Done.'
end subroutine UserTimeLagHandle
