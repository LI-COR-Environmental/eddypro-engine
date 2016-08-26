!***************************************************************************
! override_settings.f90
! ---------------------
! Copyright (C) 2007-2011, Eco2s team, Gerardo Fratini
! Copyright (C) 2011-2016, LI-COR Biosciences
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
! \brief       Forces some operations (regardless of user choice) based on instrument
!              models (e.g. CSAT3 no cross-wind correction) and logic
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine OverrideSettings()
    use m_rp_global_var
    implicit none

    !> If biomet measurements are not to be used, they are also not to be output
    if (EddyProProj%biomet_data == 'none') EddyProProj%out_biomet = .false.
    if (EddyProProj%biomet_data == 'none') EddyProProj%out_fluxnet_biomet = .false.

    !> if there is no LI-7500 among the instruments, Burba terms should not be calculated
    if (index(E2Col(co2)%Instr%model, 'li7500') == 0 &
        .and. index(E2Col(h2o)%Instr%model,'li7500') == 0) &
        RPsetup%bu_corr = 'none'
end subroutine OverrideSettings
