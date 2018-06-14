!***************************************************************************
! m_index_parameters.f90
! ----------------------
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
! \brief       Declare parameters for array indexes.
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
module m_index_parameters
    use m_numeric_kinds
    implicit none
    save


    !> other useful labels
    integer, parameter :: dum = 1

    integer, parameter :: diag72 = 1
    integer, parameter :: diag75 = 2
    integer, parameter :: diag77 = 3
    integer, parameter :: diagAnem = 4
    integer, parameter :: diagStaA = 5
    integer, parameter :: diagStaD = 6

    integer, parameter :: Fn = 1
    integer, parameter :: fc = 2
    integer, parameter :: f_half = 3
    integer, parameter :: exp1 = 4
    integer, parameter :: exp2 = 5
    integer, parameter :: exp3 = 6
    integer, parameter :: pow1 = 7
    integer, parameter :: pow2 = 8
    integer, parameter :: offsetpow1 = 9
    integer, parameter :: offsetpow2 = 10
    integer, parameter :: offsetpow3 = 11

    integer, parameter :: wc = 1
    integer, parameter :: wh = 2
    integer, parameter :: wu = 3

    integer, parameter :: RH10 = 1
    integer, parameter :: RH20 = 2
    integer, parameter :: RH30 = 3
    integer, parameter :: RH40 = 4
    integer, parameter :: RH50 = 5
    integer, parameter :: RH60 = 6
    integer, parameter :: RH70 = 7
    integer, parameter :: RH80 = 8
    integer, parameter :: RH90 = 9

    integer, parameter :: JAN = 1
    integer, parameter :: FEB = 2
    integer, parameter :: MAR = 3
    integer, parameter :: APR = 4
    integer, parameter :: MAY = 5
    integer, parameter :: JUN = 6
    integer, parameter :: JUL = 7
    integer, parameter :: AUG = 8
    integer, parameter :: SEP = 9
    integer, parameter :: OCT = 10
    integer, parameter :: NOV = 11
    integer, parameter :: DEC = 12

    integer, parameter :: daytime = 1
    integer, parameter :: nighttime = 2
    integer, parameter :: bot = 1
    integer, parameter :: top = 2
    integer, parameter :: spar = 3

    integer, parameter :: dynmd_date                        = 1
    integer, parameter :: dynmd_time                        = 2
    integer, parameter :: latitude                          = 3
    integer, parameter :: longitude                         = 4
    integer, parameter :: altitude                          = 5
    integer, parameter :: file_length                       = 6
    integer, parameter :: acquisition_frequency             = 7
    integer, parameter :: canopy_height                     = 8
    integer, parameter :: displacement_height               = 9
    integer, parameter :: roughness_length                  = 10
    integer, parameter :: master_sonic_manufacturer         = 11
    integer, parameter :: master_sonic_model                = 12
    integer, parameter :: master_sonic_height               = 13
    integer, parameter :: master_sonic_wformat              = 14
    integer, parameter :: master_sonic_wref                 = 15
    integer, parameter :: master_sonic_north_offset         = 16
    integer, parameter :: master_sonic_hpath_length         = 17
    integer, parameter :: master_sonic_vpath_length         = 18
    integer, parameter :: master_sonic_tau                  = 19
    integer, parameter :: co2_irga_manufacturer             = 20
    integer, parameter :: co2_irga_model                    = 21
    integer, parameter :: co2_measure_type                  = 22
    integer, parameter :: co2_irga_northward_separation     = 23
    integer, parameter :: co2_irga_eastward_separation      = 24
    integer, parameter :: co2_irga_vertical_separation      = 25
    integer, parameter :: co2_irga_tube_length              = 26
    integer, parameter :: co2_irga_tube_diameter            = 27
    integer, parameter :: co2_irga_tube_flowrate            = 28
    integer, parameter :: co2_irga_kw                       = 29
    integer, parameter :: co2_irga_ko                       = 30
    integer, parameter :: co2_irga_hpath_length             = 31
    integer, parameter :: co2_irga_vpath_length             = 32
    integer, parameter :: co2_irga_tau                      = 33
    integer, parameter :: h2o_irga_manufacturer             = 34
    integer, parameter :: h2o_irga_model                    = 35
    integer, parameter :: h2o_measure_type                  = 36
    integer, parameter :: h2o_irga_northward_separation     = 37
    integer, parameter :: h2o_irga_eastward_separation      = 38
    integer, parameter :: h2o_irga_vertical_separation      = 39
    integer, parameter :: h2o_irga_tube_length              = 40
    integer, parameter :: h2o_irga_tube_diameter            = 41
    integer, parameter :: h2o_irga_tube_flowrate            = 42
    integer, parameter :: h2o_irga_kw                       = 43
    integer, parameter :: h2o_irga_ko                       = 44
    integer, parameter :: h2o_irga_hpath_length             = 45
    integer, parameter :: h2o_irga_vpath_length             = 46
    integer, parameter :: h2o_irga_tau                      = 47
    integer, parameter :: ch4_irga_manufacturer             = 48
    integer, parameter :: ch4_irga_model                    = 49
    integer, parameter :: ch4_measure_type                  = 50
    integer, parameter :: ch4_irga_northward_separation     = 51
    integer, parameter :: ch4_irga_eastward_separation      = 52
    integer, parameter :: ch4_irga_vertical_separation      = 53
    integer, parameter :: ch4_irga_tube_length              = 54
    integer, parameter :: ch4_irga_tube_diameter            = 55
    integer, parameter :: ch4_irga_tube_flowrate            = 56
    integer, parameter :: ch4_irga_kw                       = 57
    integer, parameter :: ch4_irga_ko                       = 58
    integer, parameter :: ch4_irga_hpath_length             = 59
    integer, parameter :: ch4_irga_vpath_length             = 60
    integer, parameter :: ch4_irga_tau                      = 61
    integer, parameter :: gas4_irga_manufacturer            = 62
    integer, parameter :: gas4_irga_model                   = 63
    integer, parameter :: gas4_measure_type                 = 64
    integer, parameter :: gas4_irga_northward_separation    = 65
    integer, parameter :: gas4_irga_eastward_separation     = 66
    integer, parameter :: gas4_irga_vertical_separation     = 67
    integer, parameter :: gas4_irga_tube_length             = 68
    integer, parameter :: gas4_irga_tube_diameter           = 69
    integer, parameter :: gas4_irga_tube_flowrate           = 70
    integer, parameter :: gas4_irga_kw                      = 71
    integer, parameter :: gas4_irga_ko                      = 72
    integer, parameter :: gas4_irga_hpath_length            = 73
    integer, parameter :: gas4_irga_vpath_length            = 74
    integer, parameter :: gas4_irga_tau                     = 75


    !> Biomet measurements
    integer, parameter :: bTa   = 1
    integer, parameter :: bPa   = 2
    integer, parameter :: bRH   = 3
    integer, parameter :: bPPFD = 4
    integer, parameter :: bLWin = 5
    integer, parameter :: bRg   = 6

    !> logical file units
    integer ::  udf
    integer ::  udf2
    integer, parameter ::  uqc  = 106
    integer, parameter :: uspk  = 107
    integer, parameter ::  umd  = 108
    integer, parameter :: uflx  = 109
    integer, parameter :: udtc  = 110
    integer, parameter :: uvdtc = 111
    integer, parameter :: ubpcf = 112
    integer, parameter ::  ust1 = 113
    integer, parameter ::  ust2 = 114
    integer, parameter ::  ust3 = 115
    integer, parameter ::  ust4 = 116
    integer, parameter ::  ust5 = 117
    integer, parameter ::  ust6 = 118
    integer, parameter ::  ust7 = 119
    integer, parameter :: u_user_st1 = 120
    integer, parameter :: u_user_st2 = 121
    integer, parameter :: u_user_st3 = 122
    integer, parameter :: u_user_st4 = 123
    integer, parameter :: u_user_st5 = 124
    integer, parameter :: u_user_st6 = 125
    integer, parameter :: u_user_st7 = 126
    integer, parameter ::     ufnet_e = 127
    integer, parameter ::     ufnet_b = 128
    integer, parameter ::     uaflx  = 129
    integer, parameter ::     uex    = 130
    integer, parameter ::     ubiomet  = 131
    integer, parameter ::     uflxnt  = 132

    integer, parameter :: uout  = 142
    integer, parameter :: upf   = 143
    integer, parameter :: upf1  = 144
    integer, parameter :: upf2  = 145
    integer, parameter :: uxls  = 146
    integer, parameter :: utxt  = 147
    integer, parameter :: umed  = 148
    integer, parameter :: u1    = 149
    integer, parameter :: u2    = 150
    integer, parameter :: u3    = 151
    integer, parameter :: u4    = 152
    integer, parameter :: u5    = 153
    integer, parameter :: ucosp = 154
    integer, parameter :: utf   = 155
    integer, parameter :: unat  = 156
    integer, parameter :: udata = 157
    integer, parameter :: uaux  = 158
    integer, parameter :: uoff  = 159
    integer, parameter :: urh   = 160
    integer, parameter :: uto   = 161
end module m_index_parameters
