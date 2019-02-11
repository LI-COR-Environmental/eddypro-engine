!***************************************************************************
! wind_direction.f90
! ------------------
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
! \brief       Calculates wind direction from single wind components pair and offset
! \brief       It's essentially a vector direction + offset
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine SingleWindDirection(Wind, offset, WindDir)
    use m_common_global_var
    implicit none
    !> in/out variables
    real(kind = dbl), intent(in) :: Wind(3)
    real(kind = dbl), intent(in) :: offset
    real(kind = dbl), intent(out) :: WindDir


    if (Wind(U) /= error .and. Wind(V) /= error)  then 
        !> Calculate raw wind direction from wind vector
        WindDir = 180 - (datan2(Wind(V), Wind(U)) * 180d0 / p)

        !> accounts for user-supplied anemometer mis-alignment
        WindDir = WindDir + offset

        !> wrap within 0 - 360
        if (WindDir >= 360d0) WindDir = WindDir - 360d0
        if (WindDir < 0d0)   WindDir = 360d0 + WindDir
    else
        WindDir = error
    end if
end subroutine SingleWindDirection

!***************************************************************************
!
! \brief       Calculates mean wind direction and compensates offset
! \author      Gerardo Fratini
! \note        
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine AverageWindDirection(Set, nrow, ncol, offset, WindDir, err_float)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(in) :: err_float
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    real(kind = dbl), intent(in) :: offset
    real(kind = dbl), intent(out) :: WindDir
    !> Local variables
    real(kind = dbl):: wd(nrow)
    integer :: i


    !> Compute raw-level wind-direction
    do i = 1, nrow
        call SingleWindDirection(Set(i, u:w), offset, wd(i))
    end do
    
    !> Compute mean wind direction
    call AngularAverageNoError(wd, nrow, 1, WindDir, err_float)

end subroutine AverageWindDirection


!***************************************************************************
!
! \brief       Calculates standard deviation of wind direction
! \author      Gerardo Fratini
! \note        
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WindDirectionStDev(Set, nrow, ncol, WindDirStDev, err_float)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(in) :: err_float
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    real(kind = dbl), intent(out) :: WindDirStDev
    !> Local variables
    real(kind = dbl):: wd(nrow)
    integer :: i


    !> Compute raw-level wind-direction
    do i = 1, nrow
        call SingleWindDirection(Set(i, u:w), 0d0, wd(i))
    end do
    
    !> Compute mean wind direction
    call AngularStDevApproxNoError(wd, nrow, 1, WindDirStDev, err_float)

end subroutine WindDirectionStDev
