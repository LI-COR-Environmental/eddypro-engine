!***************************************************************************
! m_numeric_kinds.f90
! -------------------
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
! \brief       Define numeric kinds
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
module m_numeric_kinds

    implicit none

    !> Declare parameters
    integer, parameter :: short = kind(2)
    !integer, parameter :: int   = kind(4)
    integer, parameter :: long  = kind(8)
    integer, parameter :: sgl   = kind(0.0)
    integer, parameter :: dbl   = kind(0.0d0)
    integer, parameter :: utf8 = selected_char_kind('ISO_10646')

    integer, parameter :: i1 = selected_int_kind(2)
    integer, parameter :: i2 = selected_int_kind(4)
    integer, parameter :: i4 = selected_int_kind(9)
    integer, parameter :: i8 = selected_int_kind(18)
    integer, parameter :: sp = selected_real_kind(6,37)
    integer, parameter :: dp = selected_real_kind(15,307)

end module m_numeric_kinds
