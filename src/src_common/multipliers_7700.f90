!***************************************************************************
! multipliers_7700.f90
! --------------------
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
! \brief       Calculate multipliers for spectroscopic corrections of LI-7700
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine Multipliers7700(AirPress, AirTemp, H2OMoleFrac, A, B, C)
    use m_common_global_var
    implicit none
    !> in/out variables
    real(kind = dbl), intent(in) :: AirPress
    real(kind = dbl), intent(in) :: AirTemp
    real(kind = dbl), intent(in) :: H2OMoleFrac
    real(kind = dbl) :: A, B, C
    !> local variables
    real(kind = dbl) :: Peq
    real(kind = dbl) :: key
    real(kind = dbl) :: coll_t
    real(kind = dbl) :: coll_p

    !> Control on inputs. In case of problems set multipliers to 1 and exit subroutine
    if(H2OMoleFrac == error .or. AirTemp == error .or. AirPress == error) then
        A = 1d0
        B = 1d0
        C = 1d0
        return
    end if

    !> Calculate equivalent pressure (in kPa)
    Peq = AirPress * (1d0 + 0.455d0 * H2OMoleFrac * 1d-3) * 1d-3

    !> Retrieve parameters from LUT, by linear interpolation
    call LinearLookUp(KeyTable,  size(KeyTable, 1),  size(KeyTable, 2),  Peq, AirTemp - 273.15d0, .true., key)
    call LinearLookUp(KeyTTable, size(KeyTTable, 1), size(KeyTTable, 2), Peq, AirTemp - 273.15d0, .true., coll_t)
    call LinearLookUp(KeyPTable, size(KeyPTable, 1), size(KeyPTable, 2), Peq, AirTemp - 273.15d0, .true., coll_p)

    !> Calculate multipliers
    A = key
    B = 1d0 + (1d0 - 1.46d0 * H2OMoleFrac * 1d-3) * coll_p
    C = 1d0 + (1d0 - H2OMoleFrac * 1d-3) * coll_t !+ H2OMoleFrac  * 1d-3 * (B - 1d0)
end subroutine Multipliers7700

!***************************************************************************
!
! \brief       Calculate bi-linear interpolation, for use with LUT
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine LinearLookUp(lTable, nx, ny, inx, iny, xfirst, outz)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nx, ny
    real(kind = dbl) :: lTable(nx, ny)
    real(kind = dbl) :: inx, iny, outz
    logical :: xfirst
    !> local variables
    real(kind = dbl) :: min_x, max_x, min_y, max_y
    real(kind = dbl) :: tlz, trz, tcz, blz, brz, bcz, lcz, rcz

    !> Retrieve bounding values (up/down, both in x and y)
    call LookUp2D(lTable, nx, ny, inx, iny, min_x, max_x, min_y, max_y, tlz, trz, blz, brz)

    if (min_x == error) then
        outz = error
        return
    end if

    !> Calculate linear interpolation
    if (xfirst) then
        !> first interpolate in x
        tcz = (trz - tlz) / (max_y - min_y) * (iny - min_y) + tlz
        bcz = (brz - blz) / (max_y - min_y) * (iny - min_y) + blz
        !> now interpolate in y
        outz = (bcz - tcz) / (max_x - min_x) * (inx - min_x) + tcz
    else
        !> first interpolate in y
        lcz = (blz - tlz) / (max_x - min_x) * (inx - min_x) + tlz
        rcz = (brz - trz) / (max_x - min_x) * (inx - min_x) + trz
        !> now interpolate in x
        outz = (rcz - lcz) / (max_y - min_y) * (iny - min_y) + lcz
    end if

end subroutine LinearLookUp

!***************************************************************************
!
! \brief       Retrieve bounding values from look-up table
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine LookUp2D(lTable, nx, ny, inx, iny, min_x, max_x, min_y, max_y, tlz, trz, blz, brz)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nx, ny
    real(kind = dbl), intent(in) :: inx, iny !< input values
    real(kind = dbl), intent(out) :: min_x, max_x, min_y, max_y
    real(kind = dbl), intent(out) :: tlz, trz, blz, brz
    real(kind = dbl) :: lTable(nx, ny)
    !> local variables
    integer :: i
    integer :: j
    integer :: indxi
    integer :: indxj

    indxi = nint(error)
    indxj = nint(error)

    do i = 2, nx - 1
        if (inx > lTable(i, 1) .and. inx <= lTable(i+1, 1)) then
            min_x = lTable(i, 1)
            max_x = lTable(i+1, 1)
            indxi = i
            exit
        end if
    end do
    do j = 1, ny - 1
        if (iny > lTable(1, j) .and. iny <= lTable(1, j+1)) then
            min_y = lTable(1, j)
            max_y = lTable(1, j+1)
            indxj = j
            exit
        end if
    end do
    if (indxi /= error .and. indxj /= error) then
        tlz = lTable(indxi, indxj)
        trz = lTable(indxi, indxj + 1)
        blz = lTable(indxi + 1, indxj)
        brz = lTable(indxi + 1, indxj + 1)
    else
        min_x = error
    end if
end subroutine LookUp2D
