!***************************************************************************
! qc_flags_subs.f90
! -----------------
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
! \brief       Quality flagging according to the chosen method
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine QualityFlags(lFlux2, StDiff, DtDiff, STFlg, DTFlg, lQCFlag, printout)
    use m_common_global_var
    implicit none
    !> in/out variables
    type(FluxType), intent(in) :: lFlux2
    type(QCType), intent(in) :: StDiff
    type(QCType), intent(in) :: DtDiff
    logical, intent(in) :: printout
    type(QCType), intent(out)   :: lQCFlag
    integer, intent(out) :: STFlg(GHGNumVar)
    integer, intent(out) :: DTFlg(GHGNumVar)


    if (printout) write(*,'(a)', advance = 'no') '  Calculating quality flags..'

    !> Stationarity flags
    call PartialFlagLF(StDiff%w_co2, STFlg(w_co2))
    call PartialFlagLF(StDiff%w_h2o, STFlg(w_h2o))
    call PartialFlagLF(StDiff%w_ch4, STFlg(w_ch4))
    call PartialFlagLF(StDiff%w_gas4, STFlg(w_gas4))
    call PartialFlagLF(StDiff%w_ts,  STFlg(w_ts))
    call PartialFlagLF(StDiff%w_u,   STFlg(w_u))
    !> Developed turbulence flags
    call PartialFlagLF(DtDiff%u, DTFlg(u))
    call PartialFlagLF(DtDiff%w, DTFlg(w))
    call PartialFlagLF(DtDiff%ts, DTFlg(ts))
    DTFlg(u)  = max(DTFlg(u),  DTFlg(w))

    if (.not. EddyProProj%fcc_follows) then
        select case(Meth%qcflag(1:len_trim(Meth%qcflag)))
            case ('none')
                lQCFlag%tau = nint(error)
                lQCFlag%H = nint(error)
                lQCFlag%co2 = nint(error)
                lQCFlag%h2o = nint(error)
                lQCFlag%ch4 = nint(error)
                lQCFlag%gas4 = nint(error)
            case ('mauder_foken_04')
                !> Combined flags according to Mauder and Foken (2004)
                call GTK2Flag(STFlg(w_u),   DTFlg(u), lQCFlag%tau)
                call GTK2Flag(STFlg(w_ts),  DTFlg(w), lQCFlag%H)
                call GTK2Flag(STFlg(w_co2), DTFlg(w), lQCFlag%co2)
                call GTK2Flag(STFlg(w_h2o), DTFlg(w), lQCFlag%h2o)
                call GTK2Flag(STFlg(w_ch4), DTFlg(w), lQCFlag%ch4)
                call GTK2Flag(STFlg(w_gas4), DTFlg(w), lQCFlag%gas4)
            case ('foken_03')
                !> Combined flags according to Foken (2003), retrieved from Foken et al. (2004, HoM)
                call FokenFlag(STFlg(w_u),   DTFlg(u), lQCFlag%tau)
                call FokenFlag(STFlg(w_ts),  DTFlg(w), lQCFlag%H)
                call FokenFlag(STFlg(w_co2), DTFlg(w), lQCFlag%co2)
                call FokenFlag(STFlg(w_h2o), DTFlg(w), lQCFlag%h2o)
                call FokenFlag(STFlg(w_ch4), DTFlg(w), lQCFlag%ch4)
                call FokenFlag(STFlg(w_gas4), DTFlg(w), lQCFlag%gas4)
            case ('goeckede_06')
                !> Combined flags according to Goeckede et al. (2006)
                call GoeckedeFlag(STFlg(w_u),   DTFlg(u), lQCFlag%tau)
                call GoeckedeFlag(STFlg(w_ts),  DTFlg(w), lQCFlag%H)
                call GoeckedeFlag(STFlg(w_co2), DTFlg(w), lQCFlag%co2)
                call GoeckedeFlag(STFlg(w_h2o), DTFlg(w), lQCFlag%h2o)
                call GoeckedeFlag(STFlg(w_ch4), DTFlg(w), lQCFlag%ch4)
                call GoeckedeFlag(STFlg(w_gas4), DTFlg(w), lQCFlag%gas4)
        end select
        !> If fluxes are set to error, set to error also the quality flags
        if (lFlux2%H    == error) lQCFlag%H    = nint(error)
        if (lFlux2%h2o  == error) lQCFlag%h2o  = nint(error)
        if (lFlux2%co2  == error) lQCFlag%co2  = nint(error)
        if (lFlux2%ch4  == error) lQCFlag%ch4  = nint(error)
        if (lFlux2%gas4 == error) lQCFlag%gas4 = nint(error)
    else
        lQCFlag%tau = nint(error)
        lQCFlag%H = nint(error)
        lQCFlag%co2 = nint(error)
        lQCFlag%h2o = nint(error)
        lQCFlag%ch4 = nint(error)
        lQCFlag%gas4 = nint(error)
    end if

    if (printout) write(*, '(a)') ' Done.'
end subroutine QualityFlags

!***************************************************************************
!
! \brief       Partial flags, after Foken et al. (2004, Handbook of Microm.)
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine PartialFlagLF(val, flag)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: val
    integer, intent(out) :: flag

    select case (val)
        case (0:15)
            flag = 1
        case (16:30)
            flag = 2
        case (31:50)
            flag = 3
        case (51:75)
            flag = 4
        case (76:100)
            flag = 5
        case (101:250)
            flag = 6
        case (251:500)
            flag = 7
        case (501:1000)
            flag = 8
        case (1001:)
            flag = 9
        case default
            flag = ierror
    end select
end subroutine PartialFlagLF

!***************************************************************************
!
! \brief
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
!subroutine GhgEuropeFlagLF(STFlg, DTFlg, OAFlag)
!    use m_common_global_var
!    implicit none
!    !> in/out variables
!    integer, intent(in) :: STFlg
!    integer, intent(in) :: DTFlg
!    integer, intent(out) :: OAFlag
!
!
!    if (STFlg == idint(error) .or. DTFlg == idint(error)) then
!        OAFlag = idint(error)
!        return
!    end if
!    OAFlag = nint(error)
!
!    !> Flags combination: it might need a feedback from Foken's group.
!    if( STFlg >= 6 .and. DTFlg >= 6)                                        OAFlag = 2
!    if((STFlg >= 1 .and. STFlg <= 5) .and. (DTFlg >= 1 .and. DTFlg <= 5))   OAFlag = 1
!    if((STFlg >= 1 .and. STFlg <= 5) .and. (DTFlg >= 1 .and. DTFlg <= 2))   OAFlag = 0
!    if((STFlg >= 1 .and. STFlg <= 2) .and. (DTFlg >= 1 .and. DTFlg <= 5))   OAFlag = 0
!end subroutine GhgEuropeFlagLF

!***************************************************************************
!
! \brief       Final flags, after Mauder and Foken (2004), TK2 documentation
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine GTK2Flag(STFlg, DTFlg, OAFlag)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(inout) :: STFlg
    integer, intent(inout) :: DTFlg
    integer, intent(out) :: OAFlag


    !> stat test < 30  ==> stat flag <= 2
    !> itc test  < 30  ==> itc flag  <= 2
    !> stat test < 100 ==> stat flag <= 5
    !> itc test  < 100 ==> itc  flag <= 5
    if (STFlg == ierror .or. DTFlg == ierror) then
        OAFlag = 2
        return
    end if

    OAFlag = 2
    if((STFlg >= 1 .and. STFlg <= 2) .and. (DTFlg >= 1 .and. DTFlg <= 2)) then
        OAFlag = 0
    elseif ((STFlg >= 1 .and. STFlg <= 5) .and. (DTFlg >= 1 .and. DTFlg <= 5)) then
        OAFlag = 1
    else
        OAFlag = 2
    end if
end subroutine GTK2Flag

!***************************************************************************
!
! \brief       Partial flags, after Gockede et al. (2004, AFM)
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine GoeckedeFlag(STFlg, DTFlg, OAFlag)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(inout) :: STFlg
    integer, intent(inout) :: DTFlg
    integer, intent(out) :: OAFlag

    if (STFlg == idint(error) .or. DTFlg == idint(error)) then
        OAFlag = 5
        return
    end if

    !> From Table 2 Goecked et al. (2004, AFM)
    !> Note that the range (STFlg >= 5 .and. STFlg <= 9) .and. (DTFlg >= 5 .and. DTFlg <= 6) is not
    !> provided in the paper, here it is included in the flag 5 (last).
    !> Other ranges not supported are left with error codes
    OAFlag = 5
    if((STFlg >= 1 .and. STFlg <= 2) .and. (DTFlg >= 1 .and. DTFlg <= 2))   OAFlag = 1
    if((STFlg >= 1 .and. STFlg <= 2) .and. (DTFlg >= 3 .and. DTFlg <= 4))   OAFlag = 2
    if((STFlg >= 3 .and. STFlg <= 4) .and. (DTFlg >= 3 .and. DTFlg <= 4))   OAFlag = 3
    if((STFlg >= 3 .and. STFlg <= 4) .and. (DTFlg >= 5 .and. DTFlg <= 6))   OAFlag = 4
    if((STFlg >= 5 .and. STFlg <= 9) .and. (DTFlg >= 5 .and. DTFlg <= 9))   OAFlag = 5
end subroutine GoeckedeFlag

!***************************************************************************
!
! \brief       Partial flags, after Foken et al. (2003) as retrieved from
!              Foken et al. 2004, Handbook of Micrometeorology, Table 9.4
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine FokenFlag(STFlg, DTFlg, OAFlag)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(inout) :: STFlg
    integer, intent(inout) :: DTFlg
    integer, intent(out) :: OAFlag

    if (STFlg == idint(error) .or. DTFlg == idint(error)) then
        OAFlag = 9
        return
    end if

    !> From Table 9.4 Foken et al. 2004, HoM
    !> Note that values 7 and 8 of the overall flag was somewhat interpreted from the text
    !> of the Table, which is ambiguous (ranges of flags 7 and 8 overlap with previous ones)
    OAFlag = 9
    if(STFlg == 1                    .and. (DTFlg >= 1 .and. DTFlg <= 2))   OAFlag = 1
    if(STFlg == 2                    .and. (DTFlg >= 1 .and. DTFlg <= 2))   OAFlag = 2
    if((STFlg >= 1 .and. STFlg <= 2) .and. (DTFlg >= 3 .and. DTFlg <= 4))   OAFlag = 3
    if((STFlg >= 3 .and. STFlg <= 4) .and. (DTFlg >= 1 .and. DTFlg <= 2))   OAFlag = 4
    if((STFlg >= 1 .and. STFlg <= 4) .and. (DTFlg >= 3 .and. DTFlg <= 5))   OAFlag = 5
    if(STFlg == 5                    .and. (DTFlg >= 1 .and. DTFlg <= 5))   OAFlag = 6

    if(STFlg == 6                    .and. (DTFlg >= 1 .and. DTFlg <= 6))   OAFlag = 7
    if((STFlg >= 1 .and. STFlg <= 6) .and. DTFlg == 6)                      OAFlag = 7

    if((STFlg >= 7 .and. STFlg <= 8) .and. (DTFlg >= 1 .and. DTFlg <= 8))   OAFlag = 8
    if((STFlg >= 1 .and. STFlg <= 8) .and. (DTFlg >= 7 .and. DTFlg <= 8))   OAFlag = 8

    if(STFlg == 9                    .or.  DTFlg == 9)                      OAFlag = 9
end subroutine FokenFlag
