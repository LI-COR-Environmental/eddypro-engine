!***************************************************************************
! flush_output.f90
! ----------------
!Copyright (C) 2011, LI-COR Biosciences
!
!This file is part of EddyPro (TM).
!
!EddyPro (TM) is free software: you can redistribute it and/or modify
!it under the terms of the GNU General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.
!
!EddyPro (TM) is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License
!along with EddyPro (TM).  If not, see <http://www.gnu.org/licenses/>.
!***************************************************************************

!***************************************************************************
! \file        src/flush_output.f90
! \brief       Write stored output results on output files
! \version     3.0.0
! \date        2011-10-21
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine FlushOutput(N)
    use m_fx_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    !> local variables
    integer :: j

    do j = 1, N
        if (EddyProProj%out_full) &
            write(uflx, '(a)')   chunk_rich(j)(1:len_trim(chunk_rich(j)) - 1)
        if (EddyProProj%out_ghg_eu) &
            write(ughgeu, '(a)') chunk_eu(j)(1:len_trim(chunk_eu(j)) - 1)
        if (EddyProProj%out_amflux) &
            write(uaflx,*) chunk_amflux(j)(1:len_trim(chunk_amflux(j)) - 1)
        if (EddyProProj%out_md) &
            write(umd,*) chunk_md(j)(1:len_trim(chunk_md(j)) - 1)
    end do
end subroutine FlushOutput
