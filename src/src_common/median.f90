!***************************************************************************
! median.f90
! ----------
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
! \brief       Find the median of vec(1), ... , vec(n), using as much of the
!              quicksort algorithm as is needed to isolate it.
! \author      Patched by Gerardo Fratini from original code from Alan Miller,
!              available in the public domain at:
!              http://jblevins.org/mirror/amiller/median.f90
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine median(ivec, n, med)
    use m_common_global_var
    implicit none
    !> In/out variables
    integer, intent(in) :: n
    real(kind = dbl), intent(in), dimension(:) :: ivec(n)
    real(kind = dbl), intent(out) :: med
    !> local variables
    real(kind = dbl), dimension(:) :: vec(n)
    real(kind = dbl) :: temp, xhi, xlo, xmax, xmin
    logical :: odd
    integer :: hi, lo, nby2, nby2p1, mid, i, j, k


    vec = ivec
    nby2 = n / 2
    nby2p1 = nby2 + 1
    odd = .true.

    !> HI & LO are position limits encompassing the median.
    if (n == 2 * nby2) odd = .false.
    lo = 1
    hi = n
    if (n < 3) then
        if (n < 1) then
        med = 0.0
        return
    end if
    med = vec(1)
    if (n == 1) return
        med = 0.5*(med + vec(2))
        return
    end if

    !> Find median of 1st, middle & last values.
10  mid = (lo + hi)/2
    med = vec(mid)
    xlo = vec(lo)
    xhi = vec(hi)
    if (xhi < xlo) then ! Swap xhi & xlo
      temp = xhi
      xhi = xlo
      xlo = temp
    end if
    if (med > xhi) then
        med = xhi
    elseif(med < xlo) then
        med = xlo
    end if

    !> The basic quicksort algorithm to move all values <= the sort key (med)
    !> to the left-hand end, and all higher values to the other end.
    i = lo
    j = hi
50  do
        if (vec(i) >= med) exit
        i = i + 1
    end do
    do
        if (vec(j) <= med) exit
        j = j - 1
    end do
    if (i < j) then
        temp = vec(i)
        vec(i) = vec(j)
        vec(j) = temp
        i = i + 1
        j = j - 1
        !> Decide which half the median is in.
        if (i <= j) go to 50
    end if

    if (.NOT. odd) then
        if (j == nby2 .AND. i == nby2p1) go to 130
        if (j < nby2) lo = i
        if (i > nby2p1) hi = j
        if (i /= j) go to 100
        if (i == nby2) lo = nby2
        if (j == nby2p1) hi = nby2p1
        else
        if (j < nby2p1) lo = i
        if (i > nby2p1) hi = j
        if (i /= j) go to 100

        ! Test whether median has been isolated.
        if (i == nby2p1) return
    end if
100 if (lo < hi - 1) go to 10

    if (.NOT. odd) then
        med = 0.5*(vec(nby2) + vec(nby2p1))
        return
    end if
    temp = vec(lo)
    if (temp > vec(hi)) then
        vec(lo) = vec(hi)
        vec(hi) = temp
    end if
    med = vec(nby2p1)
    return

    !> Special case, N even, J = N/2 & I = J + 1, so the median is
    !> between the two halves of the series.   Find max. of the first
    !> half & min. of the second half, then average.
130 xmax = vec(1)
    do k = lo, j
        xmax = MAX(xmax, vec(k))
    end do
    xmin = vec(n)
    do k = i, hi
        xmin = MIN(xmin, vec(k))
    end do
    med = 0.5*(xmin + xmax)
    return
end subroutine median

!***************************************************************************
! median.f90
! ----------
! Copyright (C) 2007-2011, Eco2s team, Gerardo Fratini
! Copyright (C) 2011-2019, LI-COR Biosciences, Inc.  All Rights Reserved.
! Author: Gerardo Fratini
!
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
! \ brief      Calculate Quartiles using SAS Method 5
!              This method is the default method of SAS and is based on the 
!              empirical distribution function. 
!              Based on discussion in this paper 
!              http://www.haiweb.org/medicineprices/manual/quartiles_iTSS.pdf
! \author      Patched by Gerardo Fratini from original code
!              available in the public domain at:
!              http://fortranwiki.org/fortran/show/Quartiles

! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
double precision function quantile_sas5(x, N, qin)
    use m_common_global_var
    implicit none
    integer, intent(in) :: N
    !> In/out variables
    ! real(kind = dbl), intent(in) :: x(N)
    real(kind = dbl), intent(inout) :: x(N)
    real(kind = dbl), intent(in) :: qin
    !> Local variables
    real(kind = dbl) :: xx(N)
    real(kind = dbl), parameter :: tol = 1d-8
    real(kind = dbl) :: a,b,c
    real(kind = dbl) :: diff
    integer :: ib


    a = N * qin
    b = mod(a, 1d0)
    c = a - b

    !> Sort array

    !> This works on Mac but not Win (gfortran issue?)
    ! call sort(x, size(x), xx)


    !> This works on Mac and Win
    call HPSORT(size(x), x)
    xx = x

    !> Find quantile qin
    ib = int(c)
    diff = b - 0d0
    if (diff <= tol) then
        quantile_sas5 = (xx(ib+1) + xx(ib)) / 2d0
    else
        quantile_sas5 = xx(ib+1)
    end if
end function quantile_sas5