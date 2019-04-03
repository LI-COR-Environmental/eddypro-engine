!***************************************************************************
! planar_fit_aux_params.f90
! -------------------------
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
