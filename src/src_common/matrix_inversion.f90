!***************************************************************************
! matrix_inversion.f90
! --------------------
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
! \brief       Linear equation solution by Gauss-Jordan elimination; matrix inversion  \n
!              original definition: gauss - j(Mat, n, np, b, m, mp) \n
!              b(1:n, 1:m) is an input matrix containing the m right-hand side vectors, \n
!              on output, b(1:n, 1:m) is replaced by the corresponding set of solution vectors.
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine MatrixInversion(Mat, n, SingMat)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: n
    logical, intent(out) :: SingMat
    real(kind = dbl), intent(inout) :: Mat(n, n)
    integer, parameter :: m = 1
    integer, parameter :: mp = 1
    integer, parameter :: np = 3
    integer, parameter :: NMAX = 50
    !> local variables
    integer :: icol = 0
    integer :: irow = 0
    integer :: i
    integer :: j
    integer :: k
    integer :: l
    integer :: ll
    integer :: indxc(NMAX)
    integer :: indxr(NMAX)
    integer :: ipiv(NMAX)
    real(kind = dbl) :: b(np, mp)
    real(kind = dbl) :: big
    real(kind = dbl) :: dumm
    real(kind = dbl) :: pivinv


    b = 1.d0
    SingMat = .false.
    do j = 1, n
        ipiv(j) = 0
    end do
    do i = 1, n
        big = 0.d0
        !> Looking for the pivot element
        do j = 1, n
            if(ipiv(j) /= 1)then
                do k = 1, n
                    if (ipiv(k) == 0) then
                    if (abs(Mat(j, k)) >= big)then
                    big = abs(Mat(j, k))
                    irow = j
                    icol = k
                    end if
                    end if
                end do
            end if
        end do
        ipiv(icol) = ipiv(icol) + 1
        !> Interchanging rows, if needed, to put the pivot element on the diagonal
        if (irow /= icol) then
            do l = 1, n
                dumm = Mat(irow, l)
                Mat(irow, l) = Mat(icol, l)
                Mat(icol, l) = dumm
            end do
            do l = 1, m
                dumm = b(irow, l)
                b(irow, l) = b(icol, l)
                b(icol, l) = dumm
            end do
        end if
        !> Dividing the pivot row by the pivot element
        indxr(i) = irow
        indxc(i) = icol
        if (Mat(icol, icol) == 0.d0) then
            SingMat = .true.
        end if
        pivinv = 1.d0 / Mat(icol, icol)
        Mat(icol, icol) = 1.d0
        do l = 1, n
            Mat(icol, l) = Mat(icol, l)*pivinv
        end do
        do l = 1, m
            b(icol, l) = b(icol, l)*pivinv
        end do
        !> Reducing the rows (except for the pivot one)
        do ll = 1, n
            if(ll /= icol)then
                dumm = Mat(ll, icol)
                Mat(ll, icol) = 0.d0
                do l = 1, n
                    Mat(ll, l) = Mat(ll, l) - Mat(icol, l)*dumm
                end do
                do l = 1, m
                    b(ll, l) = b(ll, l) - b(icol, l)*dumm
                end do
            end if
        end do
    end do
    !> Unscrambling the solution in view of the column interchanges
    do l = n, 1, - 1
        if(indxr(l) /= indxc(l))then
            do k = 1, n
                dumm = Mat(k, indxr(l))
                Mat(k, indxr(l)) = Mat(k, indxc(l))
                Mat(k, indxc(l)) = dumm
            end do
        end if
    end do
end subroutine MatrixInversion
