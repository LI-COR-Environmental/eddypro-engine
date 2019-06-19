!***************************************************************************
! developed_turbulence_test.f90
! -----------------------------
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
! \brief       Assess data quality for well developed turbulent conditions
!              after Foken et al. 2005 (Handbood of Micromet.)
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine DevelopedTurbulenceTest(DtDiff)
    use m_rp_global_var
    implicit none
    ! in/out variables
    type(QCType), intent(out)   :: DtDiff
    ! local variables
    real(kind = dbl) :: sigW_mod
    real(kind = dbl) :: sigU_mod
    real(kind = dbl) :: sigT_mod
    real(kind = dbl) :: sigW_meas
    real(kind = dbl) :: sigU_meas
    real(kind = dbl) :: sigT_meas
    real(kind = dbl) :: Fcor
    real(kind = dbl) :: zplus
    real(kind = dbl) :: omg

    ! Stability ranges, after Goeckede et al. (2004, AFM)
    ! For velocity statistics:
    !               -3           -0.2              0         0.4
    !_________________________________________________________________________
    ! -oo            |             |               |          |

    ! For scalar statistics:
    !                -1        -0.0625      0       0.02
    !_________________________________________________________________________
    ! -oo            |             |        |        |

    !> Coriolis factor
    omg = 2d0 * p / (24d0 * 60d0 * 60d0)
    Fcor = dabs(2d0 * omg * dsin(Metadata%lat * p / 180d0))
    !> z+, Thomas & Foken (2002)
    zplus = 1d0

    !> Modeled characteristics after Gockede et al. 2004, AFM, Table 1
    if (Ambient%zL < -0.2d0) then
    !if (Ambient%zL > -3d0 .and. Ambient%zL < -0.2d0) then  !< bounds the test to z/L > -3 as from Goeckede et al. 2004
        sigW_mod = 1.3d0 * (1d0 - 2d0 * Ambient%zL)**(1d0/3d0)
        sigU_mod = 4.15d0 * dabs(Ambient%zL)**(1d0/8d0)
    !elseif (Ambient%zL > -0.2d0 .and. Ambient%zL < 0.4d0) then  !< bounds the test to z/L < 0.4 as from Goeckede et al. 2004
    elseif (Ambient%zL > -0.2d0) then
        sigW_mod = 0.21d0 * dlog(Fcor * zplus / Ambient%us) + 3.1d0
        sigU_mod = 0.44d0 * dlog(Fcor * zplus / Ambient%us) + 6.3d0
    else
        sigW_mod = error
        sigU_mod = error
    end if

!    !> ITC for w and u after Mauder and Foken (TK3 manual)
!    if (Ambient%zL < -0.032d0) then
!        sigW_mod = 2.00d0 * dabs(Ambient%zL)**(1d0/8d0)
!        sigU_mod = 4.15d0 * dabs(Ambient%zL)**(1d0/8d0)
!    else
!        sigW_mod = 2.7d0
!        sigU_mod = 1.3d0
!    end if
!
!    !> Neutral range definition by Thomas and Foken (2002) has priority over Foken and Wichura (1996)
!    if (Ambient%zL > -0.2d0 .and. Ambient%zL < 0.4d0) then
!        sigW_mod = 0.21d0 * dlog(Fcor * zplus / Ambient%us) + 3.1d0
!        sigU_mod = 0.44d0 * dlog(Fcor * zplus / Ambient%us) + 6.3d0
!    else
!        sigW_mod = error
!        sigU_mod = error
!    end if

    !> ITC for Ts
    if(Ambient%zL < -1) then
        sigT_mod = dabs(Ambient%zL)**(-1d0/3d0)
    elseif(Ambient%zL >= -1.d0 .and. Ambient%zL < -0.0625d0) then
        sigT_mod = dabs(Ambient%zL)**(-1d0/4d0)
    elseif(Ambient%zL >= -0.0625d0 .and. Ambient%zL < 0.02d0) then
        sigT_mod = 0.5d0 * (dabs(Ambient%zL)**(-0.5d0))
    elseif(Ambient%zL > 0.02d0) then
        sigT_mod = 1.4d0 * dabs(Ambient%zL)**(-1d0/4d0)
    else
        sigT_mod = error
    end if

    !> Measured normalized standard deviations
    if(Ambient%us /= 0d0 .and. Ambient%us /= error .and. Stats%Cov(u, u) /= error) then
        sigU_meas = dsqrt(Stats%Cov(u, u)) / Ambient%us
    else
        sigU_meas = error
    end if
    if(Ambient%us /= 0d0 .and. Ambient%us /= error .and. Stats%Cov(w, w) /= error) then
        sigW_meas = dsqrt(Stats%Cov(w, w)) / Ambient%us
    else
        sigW_meas = error
    end if
    if(Ambient%Ts /= 0d0 .and. Ambient%Ts /= error .and. Stats%Cov(ts, ts) /= error) then
        sigT_meas = dabs(dsqrt(Stats%Cov(ts, ts)) / Ambient%Ts)
    else
        sigT_meas = error
    end if

    !> Relative differences
    if (sigU_mod /= 0d0 .and. sigU_mod /= error .and. sigU_meas /= error) then
        DtDiff%u = idint(dabs((sigU_mod - sigU_meas) / sigU_mod) * 1d2)
    else
        DtDiff%u = idint(error)
    end if

    if (sigW_mod /= 0d0 .and. sigW_mod /= error .and. sigW_meas /= error) then
        DtDiff%w = idint(dabs((sigW_mod - sigW_meas) / sigW_mod) * 1d2)
    else
        DtDiff%w = idint(error)
    end if

    if (sigT_mod /= 0d0 .and. sigT_mod /= error .and. sigT_meas /= error) then
        DtDiff%ts = idint(dabs((sigT_mod - sigT_meas) / sigT_mod) * 1d2)
    else
        DtDiff%ts = idint(error)
    end if
end subroutine DevelopedTurbulenceTest
