!***************************************************************************
! configure_for_express.f90
! -------------------------
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
! \brief       Set all entries to predefined values, valid for Express
!              This bypasses all user settings in terms of processing choices
! \author      Gerardo Fratini
! \note        Angle of attack correction is set later in main,
!              because it requires master_sonic
!              which is not know at this stage, when running in embedded mode
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ConfigureForFluxnet
    use m_rp_global_var
    implicit none

    !> Error code must be -9999
    EddyProProj%err_label = '-9999'
    

end subroutine ConfigureForFluxnet
