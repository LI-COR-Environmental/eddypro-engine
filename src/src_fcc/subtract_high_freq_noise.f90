!***************************************************************************
! subtract_high_freq_noise.f90
! ----------------------------
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
! \brief       Regresses high-frequency noise and subtracts it from \n
!              spectrum
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine SubtractHighFreqNoise(lSpec, nrow, ncol, nlong, N, nclass, nbins)
    use m_fx_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    integer, intent(in) :: nbins
    integer, intent(in) :: nclass
    integer, intent(in) :: nrow, ncol
    integer, intent(in) :: nlong(N, nclass)
    type(LongSpectraType), intent(inout) :: lSpec(nrow, ncol)
    !> local variables
    integer :: gas
    integer :: cls
    integer :: imin
    integer :: ipt
    integer :: i
    real(kind = dbl) :: sumf, sumff, avgf
    real(kind = dbl) :: sumy, sumfy, avgy
    real(kind = dbl) :: m, q
    real(kind = dbl) :: noise(nrow)

    !> Linear regression of high frequency spectral range
    do gas = co2, gas4
        if (FCCsetup%SA%hfn_fmin(gas) > 0d0) then
            do cls = 1, nclass
                if (MeanBinSpecAvailable(cls, gas)) then
                    !> first, detect index of minimum frequency
                    imin = 0
                    do i = 1, nlong(gas, cls)
                        if (lSpec(i, cls)%fn(gas) > FCCsetup%SA%hfn_fmin(gas)) then
                            imin = i
                            exit
                        end if
                    end do
                    if (imin == 0) cycle
                    ipt = nlong(gas, cls) - imin + 1

                    !> Calculate linear regression of high freq spectrum
                    sumf  = 0d0
                    sumy  = 0d0
                    sumff  = 0d0
                    sumfy  = 0d0
                    do i = imin, nlong(gas, cls)
                        sumf  = sumf  + lSpec(i, cls)%fn(gas)
                        sumy  = sumy  + lSpec(i, cls)%of(gas)
                        sumff = sumff + lSpec(i, cls)%fn(gas)**2
                        sumfy = sumfy + lSpec(i, cls)%fn(gas) * lSpec(i, cls)%of(gas)
                    end do
                    avgf = sumf / dfloat(ipt)
                    avgy = sumy / dfloat(ipt)
                    m = (sumfy - sumf * avgy) / (sumff - sumf * avgf)
                    q = avgy - m * avgf
                    !> Define noise
                    noise(1:nlong(gas, cls)) = m * lSpec(1:nlong(gas, cls), cls)%fn(gas) + q

                    !> Subtract noise from spectra (both long and short)
                    do i = 1, nlong(gas, cls)
                        !if (lSpec(i, cls)%fn(gas) >= FCCsetup%SA%hfn_fmin(gas) * 0.3d0)
                            lSpec(i, cls)%of(gas) = lSpec(i, cls)%of(gas) - noise(i)
                    end do
                    do i = 1, nbins
                        !if (MeanBinSpec(i, cls)%fn(gas) >= FCCsetup%SA%hfn_fmin(gas) * 0.3d0)
                        MeanBinSpec(i, cls)%of(gas) =  MeanBinSpec(i, cls)%of(gas) &
                            -  (m * MeanBinSpec(i, cls)%fn(gas) + q)
                    end do
                end if
            end do
        end if
    end do
end subroutine SubtractHighFreqNoise
