!***************************************************************************
! fix_planarfit_sectors.f90
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
! \brief       Replaces error planar fits with closest valid sector rotation matrix
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine FixPlanarfitSectors(GoPlanarFit, N)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    logical, intent(inout) :: GoPlanarFit(N)
    !> Local variables
    integer :: sec
    integer :: sec2
    real(kind = dbl) :: loc_pfmat(3, 3, -N + 1: 2 * N)
    logical :: loc_go(-N + 1: 2 * N)
    real(kind = dbl)  :: loc_pfb(3, -N + 1: 2 * N)


    !> First, if there is no valid sector, switches to 2D rotations
    do sec = 1, N
        if (GoPlanarFit(sec)) exit
    end do
    if (sec == N + 1) then
        Meth%rot = 'double_rotation'
        call ExceptionHandler(37)
        return
    end if

    !> Define working arrays, only for convenience, for going
    !> clockwise or counter-clockwise easier
    !>
    !> N-1           0   1            N   N+1          2*N
    !>  |--o--x--*---|   |--o--x--*---|    |--o--x--*---|
    !>
    loc_pfmat(:, :, -N + 1: 0)    = PFMat(:, :, 1:N)
    loc_pfmat(:, :, 1: N)         = PFMat(:, :, 1:N)
    loc_pfmat(:, :, N + 1: 2 * N) = PFMat(:, :, 1:N)
    loc_go(-N + 1: 0)    = GoPlanarFit(1:N)
    loc_go(1: N)         = GoPlanarFit(1:N)
    loc_go(N + 1: 2 * N) = GoPlanarFit(1:N)
    loc_pfb(:, -N + 1: 0)    = PFb(:, 1:N)
    loc_pfb(:, 1: N)         = PFb(:, 1:N)
    loc_pfb(:, N + 1: 2 * N) = PFb(:, 1:N)

    do sec = 1, N
        if (.not. GoPlanarFit(sec)) then
            if(PFSetup%fix == 'clockwise') then
                !> Searches clockwise
                do sec2 = sec + 1, 2*N
                    if (loc_go(sec2)) then
                        PFMat(:, :, sec) = loc_pfmat(:, :, sec2)
                        PFb(:, sec) = loc_pfb(:, sec2)
                        GoPlanarFit(Sec) = .true.
                        exit
                    end if
                end do
            elseif(PFSetup%fix == 'counterclockwise') then
                !> Searches counterclockwise
                do sec2 = sec - 1, - N + 1
                    if (loc_go(sec2)) then
                        PFMat(:, :, sec) = loc_pfmat(:, :, sec2)
                        PFb(:, sec) = loc_pfb(:, sec2)
                        GoPlanarFit(sec) = .true.
                        exit
                    end if
                end do
            end if
        end if
    end do
end subroutine FixPlanarfitSectors
