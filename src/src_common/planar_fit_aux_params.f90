!***************************************************************************
! planarfit_aux_params.f90
! ------------------------
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
! \brief       Calculate auxiliary parameters used in planar fit calculations
!              as from Wilczak et al. (2001, BLM). MAtrix and vector definitions
!              from van Dijk et al. 2004, eq. 3.40
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine PlanarFitAuxParams(Wind, NumRow, Mat, Vec)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: NumRow
    real(kind = dbl), intent(in) :: Wind(NumRow, 3)
    real(kind = dbl), intent(out) :: Mat(3, 3)
    real(kind = dbl), intent(out) :: Vec(3)
    !> local variables
    real(kind = dbl) :: Mean(3)
    real(kind = dbl) :: SqMean(2)
    real(kind = dbl) :: UV
    real(kind = dbl) :: UW
    real(kind = dbl) :: VW
    integer :: i = 0


    Mean   = 0.d0
    SqMean = 0.d0
    UV       = 0.d0
    UW       = 0.d0
    VW     = 0.d0
    do i = 1, NumRow
        Mean(u)   = Mean(u)   + Wind(i, u)
        Mean(v)   = Mean(v)   + Wind(i, v)
        Mean(w)   = Mean(w)   + Wind(i, w)
        SqMean(u) = SqMean(u) + Wind(i, u)**2
        SqMean(v) = SqMean(v) + Wind(i, v)**2
        UV        = UV        + Wind(i, u)*Wind(i, v)
        UW        = UW        + Wind(i, u)*Wind(i, w)
        VW        = VW        + Wind(i, v)*Wind(i, w)
    end do

    Mean   = Mean   / dble(NumRow)
    SqMean = SqMean / dble(NumRow)
    UV     = UV     / dble(NumRow)
    UW     = UW     / dble(NumRow)
    VW     = VW     / dble(NumRow)

    !> Mat assuming b0 (=b(1) here) differento from 0,
    !> according to Wiczak et al. 2001.
    Mat(u, u) = 1d0
    Mat(u, v) = Mean(u)
    Mat(u, w) = Mean(v)
    Mat(v, u) = Mat(u, v)
    Mat(v, v) = SqMean(u)
    Mat(v, w) = UV
    Mat(w, u) = Mat(u, w)
    Mat(w, v) = Mat(v, w)
    Mat(w, w) = SqMean(v)

    Vec(u) = Mean(w)
    Vec(v) = UW
    Vec(w) = VW
end subroutine PlanarFitAuxParams
