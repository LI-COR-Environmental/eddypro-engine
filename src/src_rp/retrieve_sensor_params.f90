!***************************************************************************
! retrieve_sensor_params.f90
! --------------------------
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
! \brief       Retrieves sensor parameters, basing on firm and model
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo        Specify gas analyser's acquisition frequencies
!***************************************************************************
subroutine RetrieveSensorParams()
    use m_rp_global_var
    implicit none
    !> local variables
    integer :: gas

    !> anemometer parameters
    select case(E2Col(u)%Instr%model(1:len_trim(E2Col(u)%Instr%model) - 2))
        case('hs_50')
            E2Col(u:ts)%Instr%vpath_length = 0.125d0
            E2Col(u:ts)%Instr%hpath_length = 0.110d0
            E2Col(u:ts)%Instr%tau = 1d0 / 50d0
        case('hs_100')
            E2Col(u:ts)%Instr%vpath_length = 0.125d0
            E2Col(u:ts)%Instr%hpath_length = 0.110d0
            E2Col(u:ts)%Instr%tau = 1d0 / 100d0
        case('r2')
            E2Col(u:ts)%Instr%vpath_length = 0.125d0
            E2Col(u:ts)%Instr%hpath_length = 0.110d0
            E2Col(u:ts)%Instr%tau = 1d0 / 50.d0
        case('r3_50')
            E2Col(u:ts)%Instr%vpath_length = 0.125d0
            E2Col(u:ts)%Instr%hpath_length = 0.110d0
            E2Col(u:ts)%Instr%tau = 1d0 / 50d0
        case('r3_100')
            E2Col(u:ts)%Instr%vpath_length = 0.125d0
            E2Col(u:ts)%Instr%hpath_length = 0.110d0
            E2Col(u:ts)%Instr%tau = 1d0 / 100d0
        case('r3a_100')
            E2Col(u:ts)%Instr%vpath_length = 0.125d0
            E2Col(u:ts)%Instr%hpath_length = 0.110d0
            E2Col(u:ts)%Instr%tau = 1d0 / 100d0
        case('wm')
            E2Col(u:ts)%Instr%vpath_length = 0.125d0
            E2Col(u:ts)%Instr%hpath_length = 0.110d0
            E2Col(u:ts)%Instr%tau = 1d0 / 32d0
        case('wmpro')
            E2Col(u:ts)%Instr%vpath_length = 0.125d0
            E2Col(u:ts)%Instr%hpath_length = 0.110d0
            E2Col(u:ts)%Instr%tau = 1d0 / 40d0
        case('usa1_standard')
            E2Col(u:ts)%Instr%vpath_length = 0.175d0    !< to be adjusted
            E2Col(u:ts)%Instr%hpath_length = 0.175d0    !< to be adjusted
            E2Col(u:ts)%Instr%tau = 1d0 / 25d0
        case('usa1_fast')
            E2Col(u:ts)%Instr%vpath_length = 0.175d0    !< to be adjusted
            E2Col(u:ts)%Instr%hpath_length = 0.175d0    !< to be adjusted
            E2Col(u:ts)%Instr%tau = 1d0 / 50d0
        case('csat3', 'csat3b')
            E2Col(u:ts)%Instr%vpath_length = 0.115d0
            E2Col(u:ts)%Instr%hpath_length = 0.058d0
            E2Col(u:ts)%Instr%tau = 1d0 / 60d0
        case('81000', '81000v', '81000re', '81000vre')
            E2Col(u:ts)%Instr%vpath_length = 0.150d0
            E2Col(u:ts)%Instr%tau = 1d0 / 40d0
    end select

    !> gas analyzer parameters
    do gas = co2, gas4
        if(E2Col(gas)%present) then
            select case (E2Col(gas)%Instr%model(1:len_trim(E2Col(gas)%Instr%model) - 2))
                case('li7500')
                    E2Col(gas)%Instr%vpath_length = 0.127d0
                    E2Col(gas)%Instr%hpath_length = 0.0095d0
                    E2Col(gas)%Instr%tau = 1d-1
                case('li7500a','li7500rs')
                    E2Col(gas)%Instr%vpath_length = 0.127d0
                    E2Col(gas)%Instr%hpath_length = 0.0095d0
                    E2Col(gas)%Instr%tau = 1d-1
                case('li7200', 'li7200rs')
                    E2Col(gas)%Instr%vpath_length = 0.127d0
                    E2Col(gas)%Instr%hpath_length = 0.0127d0
                    E2Col(gas)%Instr%tau = 1d-1
                case('li7700')
                    E2Col(gas)%Instr%vpath_length = 0.50d0
                    E2Col(gas)%Instr%hpath_length = 0.0635d0
                    E2Col(gas)%Instr%tau = 1d-1
                case('li6262')
                    E2Col(gas)%Instr%vpath_length = 0.152d0
                    E2Col(gas)%Instr%hpath_length = 0.013d0
                    E2Col(gas)%Instr%tau = 1d-1
                case('li7000')
                    E2Col(gas)%Instr%vpath_length = 0.152d0
                    E2Col(gas)%Instr%hpath_length = 0.0095d0
                    E2Col(gas)%Instr%tau = 1d-1
            end select
        end if
    end do
end subroutine RetrieveSensorParams

