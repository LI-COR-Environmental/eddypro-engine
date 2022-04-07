!***************************************************************************
! qn_estimator.f90
! ----------------
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
!--------------------------------------------------------------------
!  Note from the original implementation, which can be found at 
!  http://ftp.uni-bayreuth.de/math/statlib/general/snqn
!
!  This file contains a Fortran function for the robust estimator
!  of scale denoted as Qn, proposed in Rousseeuw and Croux (1993).
!  The estimator has a high breakdown point and a smooth and bounded
!  influence function. The algorithm given here is very fast (running
!  in O(nlogn) time) and needs only O(n) storage space.
!
!  Rousseeuw, P.J. and Croux, C. (1993), "Alternatives to the
!     Median Absolute Deviation," Journal of the American
!     Statistical Association, Vol. 88, 1273-1283.
!
!  This software may be used and copied freely, provided
!  reference is made to the abovementioned paper.
!
!***************************************************************************
!
! \brief       Compute the scale estimator Qn 
!              after Rousseeuw, P.J. and Croux, C. (1993)
! \author      Adapted by Gerardo Fratini from original code
! \note        Original code:
!              http://ftp.uni-bayreuth.de/math/statlib/general/snqn
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
!
!   Efficient algorithm for the scale estimator:
!
!       Qn = dn * 2.2219 * {|x_i-x_j|; i<j}_(k)
!
!   Parameters of the function Qn :
!       x  : real array containing the observations
!       n  : number of observations (n >=2)
!
!   The function Qn uses the procedures:
!      whimed(a,iw,n): finds the weighted high median of an array
!                      a of length n, using the array iw (also of
!                      length n) with positive integer weights.
!      qn_sort(x,n,y) : sorts an array x of length n, and stores the
!                    result in an array y (of size at least n)
!      pull(a,n,k) : finds the k-th order statistic of an
!                    array a of length n
!
double precision function Qn(x, n)
    use m_numeric_kinds
    implicit none
    !> in/out variables
    integer, intent(in) :: n
    real(kind = dbl), intent(in) :: x(n)
    !> local variables
    real(kind = dbl) :: y(500), work(500)
    real(kind = dbl) :: dn, trial
    integer :: left(500), right(500), weight(500), Q(500), P(500)
    integer :: i, j, jj, h, k, knew, jhelp, nL, nR, sumQ, sumP
    logical :: found

    double precision, external :: whimed
    double precision, external :: pull

    h = n/2 + 1
    k = h * (h-1) / 2
    call qn_sort(x,n,y)
    do i = 1, n
        left(i) = n - i + 2
        right(i) = n
    end do

    jhelp = n*(n+1)/2
    knew = k + jhelp
    nL = jhelp
    nR = n * n
    found = .false.

200 continue

    if ( (nR-nL > n) .and. (.not.found) ) then
        j=1
        do i = 2, n
            if (left(i) <= right(i)) then
                weight(j) = right(i) - left(i) + 1
                jhelp = left(i) + weight(j) / 2
                work(j) = y(i) - y(n + 1 - jhelp)
                j = j + 1
            endif
        end do

        trial = whimed(work, weight, j - 1)
        j=0

        do i = n, 1, -1
45          continue
            if ( (j < n) .and. ((y(i) - y(n-j)) < trial) ) then
                j = j + 1
                goto 45
            endif
            P(i) = j
        end do
        j = n + 1

        do i = 1, n
55          continue
            if ( (y(i) - y(n - j + 2) ) > trial ) then
                j = j - 1
                goto 55
            endif
            Q(i) = j
        end do
        sumP = 0
        sumQ = 0

        do i = 1, n
            sumP = sumP + P(i)
            sumQ = sumQ + Q(i) -1
        end do
              
        if (knew <= sumP) then
            do i = 1, n
                right(i) = P(i)
            end do
            nR = sumP
        else
            if (knew > sumQ) then
                do i = 1, n
                    left(i)=Q(i)
                end do
                nL = sumQ
            else
                Qn = trial
                found = .true.
            endif
        endif
        goto 200
    endif


    if (.not.found) then
        j = 1
        do i = 2, n
            if ( left(i) <= right(i) ) then
                do jj = left(i), right(i)
                    work(j) = y(i) - y(n - jj + 1)
                    j = j + 1
                end do
            endif
        end do
        Qn = pull(work, j-1, knew-nL)
    endif

    dn = -9999.
    if (n <= 9) then
        if (n == 2) dn=0.399
        if (n == 3) dn=0.994
        if (n == 4) dn=0.512
        if (n == 5) dn=0.844
        if (n == 6) dn=0.611
        if (n == 7) dn=0.857
        if (n == 8) dn=0.669
        if (n == 9) dn=0.872
    else
        if (mod(n, 2) == 1) dn = n / (n + 1.4)
        if (mod(n, 2) == 0) dn = n / (n + 3.8)
    endif

    if (dn /= -9999.) then
        Qn = dn * 2.2219 * Qn
    else
        Qn = -9999.
    end if

end function Qn

!***************************************************************************
!
! \brief       Whimed
! \author      Adapted by Gerardo Fratini from original code
! \note        Original code:
!              http://ftp.uni-bayreuth.de/math/statlib/general/snqn
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
!
!  Algorithm to compute the weighted high median in O(n) time.
!
!  The whimed is defined as the smallest a(j) such that the sum
!  of the weights of all a(i) <= a(j) is strictly greater than
!  half of the total weight.
!
!  Parameters of this function:
!        a: real array containing the observations
!        n: number of observations
!       iw: array of integer weights of the observations.
!
!  This function uses the function pull.
!
!  The size of acand, iwcand must be at least n.
!
double precision function whimed(a, iw, n)
    use m_numeric_kinds
    implicit none
    !> in/out variables
    integer, intent(in) :: n
    real(kind = dbl), intent(inout) :: a(n)
    integer, intent(inout) :: iw(n)
    !> local variables
    integer :: iwcand(500)
    real(kind = dbl) :: acand(500)
    real(kind = dbl) :: trial, wtotal, wleft, wright, wrest, wmid
    integer :: i, kcand, nn

    double precision, external :: pull

    nn = n
    wtotal = 0
    do i = 1, nn
        wtotal = wtotal + iw(i)
    end do
    wrest=0

100 continue
    trial = pull(a, nn, nn / 2 + 1)
    wleft = 0
    wmid = 0
    wright = 0
    do i = 1, nn
        if (a(i) < trial) then
            wleft = wleft + iw(i)
        else
            if (a(i) > trial) then
                wright = wright + iw(i)
            else
                wmid = wmid + iw(i)
            endif
        endif
    end do
          
    if ( (2 * wrest + 2 * wleft) > wtotal ) then
        kcand=0
        do i = 1, nn
            if (a(i) < trial) then
                kcand = kcand + 1
                acand(kcand) = a(i)
                iwcand(kcand) = iw(i)
            endif
        end do
        nn = kcand
    else
        if ( (2 * wrest + 2 * wleft + 2 * wmid) > wtotal ) then
            whimed = trial
            return
        else
            kcand=0
            do i = 1, nn
                if(a(i).gt.trial) then
                    kcand = kcand + 1
                    acand(kcand) = a(i)
                    iwcand(kcand) = iw(i)
                endif
            end do
            nn=kcand
            wrest=wrest+wleft+wmid
        endif
    endif
    do i = 1, nn
          a(i) = acand(i)
          iw(i) = iwcand(i)
    end do
    go to 100
end function whimed

!***************************************************************************
!
! \brief       Sorts an array a of length n<=1000, and stores
!              the result in the array b.
! \author      Adapted by Gerardo Fratini from original code
! \note        Original code:
!              http://ftp.uni-bayreuth.de/math/statlib/general/snqn
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine qn_sort(a, n, b)
    use m_numeric_kinds
    implicit none
    !> in/out variables
    integer, intent(in) :: n
    real(kind = dbl), intent(in) :: a(n)
    real(kind = dbl), intent(inout) :: b(n)
    !> local variables
    integer :: jlv(1000), jrv(1000)
    real(kind = dbl) :: amm, xx
    integer :: j, jss, jtwe, jr, jndl, jnc

    do j = 1, n
          b(j)=a(j)
    end do
    jss = 1
    jlv(1) = 1
    jrv(1) = n
10  jndl = jlv(jss)
    jr = jrv(jss)
    jss = jss - 1
20  jnc = jndl
    j = jr
    jtwe = (jndl + jr) / 2
    xx = b(jtwe)
30  if (b(jnc) >= xx) goto 40
    jnc = jnc + 1
    goto 30
40  if ( xx >= b(j)) goto 50
    j = j - 1
    goto 40
50  if (jnc > j) goto 60
    amm = b(jnc)
    b(jnc) = b(j)
    b(j) = amm
    jnc = jnc + 1
    j = j - 1
60  if (jnc <= j) goto 30
    if ( (j - jndl) < (jr - jnc) ) goto 80
    if (jndl >= j) goto 70
    jss = jss + 1
    jlv(jss) = jndl
    jrv(jss) = j
70  jndl = jnc
    goto 100
80  if (jnc >= jr) goto 90
    jss = jss + 1
    jlv(jss) = jnc
    jrv(jss) = jr
90  jr = j
100 if ( jndl < jr) goto 20
    if ( jss /= 0) goto 10
end subroutine qn_sort

!***************************************************************************
!
! \brief       Finds the kth order statistic of an array a of length n<=1000
! \author      Adapted by Gerardo Fratini from original code
! \note        Original code:
!              http://ftp.uni-bayreuth.de/math/statlib/general/snqn
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
double precision function pull(a, n, k)
    use m_numeric_kinds
    implicit none
    !> in/out variables
    integer, intent(in) :: n
    integer, intent(in) :: k
    real(kind = dbl), intent(in) :: a(n)
    !> local variables
    real(kind = dbl) :: b(1000)
    real(kind = dbl) :: ax, buffer
    integer :: lr, l, jnc, j


    do j = 1, n
        b(j) = a(j)
    end do
    l = 1
    lr = n
20  if (l >= lr) goto 90
    ax = b(k)
    jnc = l
    j = lr
30  if (jnc > j) goto 80
40  if (b(jnc) >= ax) goto 50
    jnc = jnc + 1
    goto 40
50  if (b(j) <= ax)goto 60
    j = j - 1
    goto 50
60  if (jnc > j) goto 70
    buffer = b(jnc)
    b(jnc) = b(j)
    b(j) = buffer
    jnc = jnc+1
    j = j - 1
70  goto 30
80  if (j < k) l = jnc
    if (k < jnc) lr = j
    goto 20
90  pull = b(k)
end function pull

