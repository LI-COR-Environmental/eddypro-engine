!***************************************************************************
! sort.f90
! --------
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
!>           
!>  \brief     Sorts X into ascending order - Quicksort
!>             Quicksort chooses a "pivot" in the set, and explores the
!>             array from both ends, looking for a value > pivot with the
!>             increasing index, for a value <= pivot with the decreasing
!>             index, and swapping them when it has found one of each.
!>             The array is then subdivided in 2 ([3]) subsets:
!>             { values <= pivot} {pivot} {values > pivot}
!>             One then call recursively the program to sort each subset.
!>             When the size of the subarray is small enough, one uses an
!>             insertion sort that is faster for very small sets.
!>             Michel Olagnon - Apr. 2000
!>             http://www.fortran-2000.com/rank/refsor.f90
!
! \author      Patched by Gerardo Fratini from original code from Alan Miller,
!              available in the public domain at:
!              http://www.fortran-2000.com/rank/refsor.f90
!              See also here:
!              https://jblevins.org/mirror/amiller/
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine sort (X, N, XX)
    use m_common_global_var
    implicit none

    integer, intent (in) :: N
    real (kind=dbl), intent (inout) :: X(N)
    real (kind=dbl), intent (inout) :: XX(N)

    XX = X
    call D_subsor(XX, N, 1, Size(X))
    call D_inssor(XX, N)
    return
end subroutine sort

!***************************************************************************
!>  Sorts X from IDEB1 to IFIN1
Recursive subroutine D_subsor(X, N, IDEB1, IFIN1)
    use m_common_global_var
    implicit none

    integer, intent (in) :: N
    real(kind=dbl), intent(inout) :: X(N)
    integer, intent (in) :: IDEB1, IFIN1

    integer, parameter :: NINS = 16 ! Max for insertion sort
    integer :: ICRS, IDEB, IDCR, IFIN, IMIL
    real(kind=dbl) :: XPIV, XWRK
    IDEB = IDEB1
    IFIN = IFIN1

    !>  If we don't have enough values to make it worth while, we leave
    !>  them unsorted, and the final insertion sort will take care of them
    if ((IFIN - IDEB) > NINS) then
        IMIL = (IDEB+IFIN) / 2

        !>  One chooses a pivot, median of 1st, last, and middle values
        if (X(IMIL) < X(IDEB)) then
            XWRK = X(IDEB)
            X(IDEB) = X(IMIL)
            X(IMIL) = XWRK
        end if
        if (X(IMIL) > X(IFIN)) then
            XWRK = X(IFIN)
            X(IFIN) = X(IMIL)
            X(IMIL) = XWRK
            if (X(IMIL) < X(IDEB)) then
                XWRK = X(IDEB)
                X(IDEB) = X(IMIL)
                X(IMIL) = XWRK
            end if
        end if
        XPIV = X(IMIL)

        !  One exchanges values to put those > pivot in the end and
        !  those <= pivot at the beginning
        ICRS = IDEB
        IDCR = IFIN
        ECH2: do
            do
                ICRS = ICRS + 1
                if (ICRS >= IDCR) then

                    !  the first  >  pivot is IDCR
                    !  the last   <= pivot is ICRS-1
                    !  Note: If one arrives here on the first iteration, then
                    !        the pivot is the maximum of the set, the last value is equal
                    !        to it, and one can reduce by one the size of the set to process,
                    !        as if X(IFIN) > XPIV
                    Exit ECH2
                end if
                if (X(ICRS) > XPIV) Exit
            end do
            do
                if (X(IDCR) <= XPIV) Exit
                IDCR = IDCR - 1
                if (ICRS >= IDCR) then
                    !  The last value < pivot is always ICRS-1
                    Exit ECH2
                end if
            end do
            XWRK = X(IDCR)
            X(IDCR) = X(ICRS)
            X(ICRS) = XWRK
        end do ECH2
        !>  One now sorts each of the two sub-intervals
        call D_subsor (X, N, IDEB1, ICRS-1)
        call D_subsor (X, N, IDCR, IFIN1)
    end if
    return
end subroutine D_subsor
 
 
subroutine D_inssor (X, N)
    use m_common_global_var
    implicit none

    integer, intent (in) :: N
    real(kind=dbl), intent (inout) :: X(N)
    integer :: ICRS, IDCR
    real(kind=dbl) :: XWRK

    do ICRS = 2, Size (X)
        XWRK = X(ICRS)
        if (XWRK >= X(ICRS-1)) Cycle
        X(ICRS) = X(ICRS-1)
        do IDCR = ICRS - 2, 1, - 1
            if (XWRK >= X(IDCR)) Exit
            X(IDCR+1) = X(IDCR)
        end do
        X(IDCR+1) = XWRK
    end do
    return
end subroutine D_inssor
