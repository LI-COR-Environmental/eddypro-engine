!***************************************************************************
! exception_handler.f90
! ---------------------
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
! \brief       Manages error and Warning instances, possibly aborting execution
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ExceptionHandler(error_code)
    implicit none
    !> in/out variables
    integer :: error_code


    select case (error_code)
        case(0)
            write(*,*) ' Fatal error(0)> Occurred while retrieving user home path.'
            write(*,*) ' Fatal error(0)> EddyPro is not able to locate project files and must terminate.'
            write(*,*) ' Fatal error(0)> Program execution aborted.'
            stop 1
        case(1)
            write(*,*) ' Fatal error(1)> Temporary file "flist.tmp" not created.'
            write(*,*) ' Fatal error(1)> Most likely, no files matching the raw file name format were found.'
            write(*,*) ' Fatal error(1)> Program execution aborted.'
            stop 1
        case(2)
            write(*,*) ' Error(2)> Occurred while scanning file for biomet data.'
            write(*,*) ' Error(2)> EddyPro could not retrieve biomet data from this file.'
            write(*,*) ' Error(2)> File will be ignored.'
        case(3)
            write(*,*) ' Error(3)> Occurred while reading metadata file in GHG archive. File not found or empty.'
            write(*,*) ' Error(3)> GHG archive appears to be corrupted or invalid and will thus be skipped.'
        case(4)
            write(*,*) ' Error(4)> Occurred while reading raw data file in GHG archive. File not found or empty.'
            write(*,*) ' Error(4)> GHG archive appears to be corrupted or invalid and will thus be skipped.'
        case(5)
            write(*,*) ' Error(5)> Occurred while reading biomet data file or corresponding metadata file in GHG archive.'
            write(*,*) ' Error(5)> File(s) not found or empty. Biomet data not used for this time period.'
        case(6)
            write(*,*) ' Error(6)> Occurred while opening current raw data file. File is empty.'
            write(*,*) ' Error(6)> Raw data file skipped.'
        case(7)
            write(*,*)
            write(*,*) ' Error(7)> Occurred while opening INI-format file. Looking for a solution..'
        case(8)
            write(*,*) ' Fatal error(8)> No files matching the specified extension were found in the selected folder.'
            write(*,*) ' Fatal error(8)> Program execution aborted.'
            stop 1
        case(14)
            write(*,*) ' Error(14)> Occurred while unzipping GHG archive. File skipped.'
        case(20)
            write(*,*) '  Fatal error(20)> Incorrect or unsupported "Raw file name format".'
            write(*,*) '  Fatal error(20)> Check it and try again, or rename raw files appropriately (see software documentation).'
            write(*,*) '  Fatal error(20)> Program execution aborted.'
            stop 1
        case(21)
            write(*,*) '  Fatal error(21)> Configuration file (project file with extension *.eddypro) not found.'
            write(*,*) '  Fatal error(21)> Program execution aborted.'
            stop 1
        case(22)
            write(*,*) '  Fatal error(22)> "Alternative metadata file" not found.'
            write(*,*) '  Fatal error(22)> Program execution aborted.'
            stop 1
        case(23)
            write(*,*) '  Fatal error(23)> Occurred while validating "Alternative metadata file".'
            write(*,*) '  Fatal error(23)> Check entries in the "Metadata file editor" and try again.'
            write(*,*) '  Fatal error(23)> Program execution aborted.'
            stop 1
        case(24)
            write(*,*) '  Error(24)> Occurred while opening GHG file.'
            write(*,*) '  Error(24)> File skipped.'
        case(25)
            write(*,*) '  Error(25)> Occurred while validating embedded metadata file.'
            write(*,*) '  Error(25)> GHG file skipped.'
        case(28)
            write(*,*) '  Error(28)> Occurred while opening raw data file.'
            write(*,*) '  Error(28)> Raw file skipped.'
        case(29)
            write(*,*) ' Alert(29)> Occurred while reading rotation matrices from auxiliary planar-fit file.'
            write(*,*) ' Alert(29)> Axis rotation method switched to "Double rotation".'
        case(30)
            write(*,*) ' Alert(30)> Occurred while opening auxiliary file for planar-fit. File not found or empty.'
            write(*,*) ' Alert(30)> Axis rotation method switched to "Double rotation".'
        case(31)
            write(*,*) ' Fatal error(31)> Unrecognized TOB1 data format. Select either IEEE4 or FP2 data formats and try again.'
            write(*,*) ' Fatal error(31)> Program execution aborted.'
            stop 1
        case(32)
            write(*,*) ' Fatal error(32)> No valid GHG file was found in the selected "Raw data directory".'
            write(*,*) ' Fatal error(32)> The problem could also be due to corrupted archives or invalid metadata files.'
            write(*,*) ' Fatal error(32)> Program execution aborted.'
            stop 1
        case(33)
            write(*,*) ' Error(33)> Number of wind records for this sector is less than requested.'
            write(*,*) ' Error(33)> Planar-fit rotation matrix not calculated for this sector.'
        case(34)
            write(*,*) ' Error(34)> Occurred while calculating planar-fit rotations &
                                     &for this sector: singular matrix found.'
            write(*,*) ' Error(34)> Planar-fit rotation matrix not calculated for this sector.'
        case(35)
            write(*,*) ' Fatal error(35)> Oops! Something went wrong. EddyPro was not able to process any raw file.'
            write(*,*) ' Fatal error(35)> Output files not created.'
            write(*,*) ' Fatal error(35)> Program execution aborted.'
            stop 1
        case(36)
            write(*,*) ' Fatal error(36)> No "Output directory" was selected. Select an "Output directory" before running EddyPro.'
            write(*,*) ' Fatal error(36)> Program execution aborted.'
            stop 1
        case(37)
            write(*,*) ' Alert(37)> No valid planar-fit rotation matrix found for any wind sector.'
            write(*,*) ' Alert(37)> Axis rotation method switched to "Double rotation".'
        case(38)
            write(*,*) ' Alert(38)> No sectors selected for planar fit.'
            write(*,*) ' Alert(38)> Forcing to 1 sector of 360 degrees.'
        case(39)
            write(*,*) ' Alert(39)> Error while opening auxiliary file for time-lag optimization. File not found or empty.'
            write(*,*) ' Alert(39)> Time-lag detection method switched to "Covariance maximization".'
        case(40)
            write(*,*) ' Alert(40)> Occurred while retrieving wind sectors configuration for planar-fit.'
            write(*,*) ' Alert(40)> Forcing to 1 sector of 360 degrees.'
        case(41)
            write(*,*) ' Warning(41)> Wind-sector excluded by user. No planar-fit matrix calculated for this sector.'
        case(42)
            write(*,*) ' Error(42)> Method for random uncertainty estimation not recognized.'
            write(*,*) ' Error(42)> Random uncertainty not calculated.'
        case(43)
            write(*,*) ' Alert(43)> Time-lag optimization failed.'
            write(*,*) ' Alert(43)> Switching to method "Covariance maximization" for time-lag detection.'
        case(44)
            write(*,*) ' Error(44)> Occurred while reading or interpreting biomet file.'
            write(*,*) ' Error(44)> Biomet data not used for this period.'
        case(45)
            write(*,*) ' Error(45)> Not enough valid co-spectra were found for fitting models, or fitting procedure failed.'
            write(*,*) ' Error(45)> Stability-sorted ensemble averaged cospectra outputs not created.'
        case(46)
            write(*,*) '  Fatal error(46)> The dataset does not contain any raw file &
                                           &corresponding to the selected sub-period.'
            write(*,*) '  Fatal error(46)> Select another sub-period or a different "Raw data directory".'
            write(*,*) '  Fatal error(46)> Try also un-checking the option "Select sub-period" or'
            write(*,*) '  Fatal error(46)> checking the option "Search in sub-folders", in the "Basic settings page".'
            write(*,*) '  Fatal error(46)> Program execution aborted.'
            stop 1
        case(48)
            write(*,*) '  Fatal error(48)> The dataset does not contain any raw file &
                                           &corresponding to the sub-period selected for planar-fit.'
            write(*,*) '  Fatal error(48)> Select another sub-period or another "Raw data directory".'
            write(*,*) '  Fatal error(48)> Try also un-checking the option "Select sub-period" in the planar-fit dialogue.'
            write(*,*) '  Fatal error(48)> Program execution aborted.'
            stop 1
        case(49)
            write(*,*) '  Fatal error(49)> The dataset does not contain any raw file &
                                           &corresponding to the sub-period selected for time-lag optimization.'
            write(*,*) '  Fatal error(49)> Select another sub-period or another "Raw data directory".'
            write(*,*) '  Fatal error(49)> Try also un-checking the option "Select sub-period" &
                                           &in the time-lag optimization dialogue.'
            write(*,*) '  Fatal error(49)> Program execution aborted.'
            stop 1
        case(50)
            write(*,*) '  Fatal error(50)> "Essentials" files does not contain any results &
                                           &corresponding to the selected sub-period.'
            write(*,*) '  Fatal error(50)> Select another sub-period.'
            write(*,*) '  Fatal error(50)> Try also un-checking the option "Select sub-period" in "Basic settings" page.'
            write(*,*) '  Fatal error(50)> Program execution aborted.'
            stop 1
        case(51)
            write(*,*) '  Error(51)> Calculation of tube attenuation parameter (lambda) failed for at least one gas.'
            write(*,*) '  Error(51)> Processing continues but tube attenuation is not included in the spectral correction.'
        case(52)
            write(*,*) '  Fatal error(52)> An unexpected, unrecognized internal problem occurred.'
            write(*,*) '  Fatal error(52)> Program execution aborted.'
            stop 1
        case(53)
            write(*,*) '  Warning(53)> No raw data file relevant to current averaging period was found.'
            write(*,*) '  Warning(53)> Skipping to next averaging period.'
        case(54)
            write(*,*) '  Error(54)> Occurred while opening or reading SLT (binary) file.'
            write(*,*) '  Error(54)> Raw file skipped.'
        case(55)
            write(*,*) '  Error(55)> Occurred while opening or reading generic binary file.'
            write(*,*) '  Error(55)> Raw file skipped.'
        case(56)
            write(*,*) '  Error(56)> Occurred while opening or reading TOB1 (binary) file.'
            write(*,*) '  Error(56)> Raw file skipped.'
        case(57)
            write(*,*) '  Error(57)> Occurred while opening or reading ASCII file.'
            write(*,*) '  Error(57)> Raw file skipped.'
        case(58)
            write(*,*) ' Warning(58)> Available samples not enough for an averaging period.'
            write(*,*) ' Warning(58)> Skipping to next averaging period.'
        case(59)
            write(*,*) '  Error(59)> At least one wind component appears to be corrupted (too many implausible values).'
            write(*,*) '  Error(59)> This may also be the result of data exclusion by the "Absolute limits" test or by a'
            write(*,*) '  Error(59)> custom-designed "Flag" in the "Basic Settings" page.'
            write(*,*) '  Error(59)> If the problem occurs for many or all raw files, check those settings.'
            write(*,*) '  Error(59)> Skipping to next averaging period.'
        case(60)
            write(*,*) ' Fatal error(60)> Occurred while opening or reading "essentials" file.'
            write(*,*) ' Fatal error(60)> Execution will be aborted, but you may be able to avoid re-processing raw data'
            write(*,*) ' Fatal error(60)> by fixing the problem with the "essentials" file and using the option'
            write(*,*) ' Fatal error(60)> "Previous data directory" in the "Basic Settings" page.'
            write(*,*) ' Fatal error(60)> Please consult software documentation.'
            write(*,*) ' Fatal error(60)> Program execution aborted.'
            stop 1
        case(61)
            write(*,*) ' Fatal error(61)> No valid data records found in the "essentials" file.'
            write(*,*) ' Fatal error(61)> Program execution aborted.'
            stop 1
        case(62)
            write(*,*) ' Error(62)> Occurred while reading binned (co)spectra file.'
            write(*,*) ' Error(62)> File skipped.'
        case(63)
            write(*,*) ' Error(63)> Occurred while reading full (co)spectra file.'
            write(*,*) ' Error(63)> File skipped.'
        case(64)
            write(*,*) ' Error(64)> Occurred while creating output file.'
            write(*,*) ' Error(64)> Some spectral assessment results will not be written on output file.'
        case(65)
            write(*,*) ' Alert(65)> Occurred while reading auxiliary "spectral assessment" file.'
            write(*,*) ' Alert(65)> High-frequency spectral correction method switched to Moncrieff et al. (1997).'
        case(66)
            write(*,*) ' Alert(66)> Acquisition frequency appears to be set to a value <= zero.'
            write(*,*) ' Alert(66)> High-frequency spectral corrections cannot be calculated.'
            write(*,*) ' Alert(66)> Proceeding without spectral corrections.'
        case(67)
            write(*,*) ' Error(67)> Occurred while opening output file.'
            write(*,*) ' Error(67)> Output dataset not created.'
        case(68)
            write(*,*) ' Error(68)> Occurred while reading dynamic metadata file. File not found or empty.'
            write(*,*) ' Error(68)> Dynamic metadata not used in this run.'
        case(69)
            write(*,*) ' Error(69)> There is a problem with results of the spectral assessment.'
            write(*,*) ' Error(69)> High-frequency spectral correction method switched to Moncrieff et al. (1997).'
        case(70)
            write(*,*) ' Error(70)> Inconsistent number of variables in biomet files.'
            write(*,*) ' Error(70)> EddyPro cannot resolve the conflict and will thus proceed without using biomet data.'
        case(71)
            write(*,*) ' Error(71)> No valid biomet record imported.'
            write(*,*) ' Error(71)> EddyPro will proceed without using biomet data.'
        case(72)
            write(*,*) '  Warning(72)> No valid biomet record found for this period.'
        case(73)
            write(*,*) '  Error(73)> The label of at least one biomet variable misses'
            write(*,*) '  Error(73)> or has incomplete positional qualifier ("_x_y_z" suffix).'
            write(*,*) '  Error(73)> EddyPro will proceed without using biomet data.'
        case(74)
            write(*,*) '  Error(74)> No valid binned (co)spectra files were found'
            write(*,*) '  Error(74)> EddyPro cannot perform spectral asssessment, nor '
            write(*,*) '  Error(74)> create ensemble averaged (co)spectra. If the case, spectral '
            write(*,*) '  Error(74)> correction method will be switched to Moncrieff et al. (2007).'
        case(75)
            write(*,*) ' Error(75)> Not enough valid co-spectra were found for making ensemble averages.'
            write(*,*) ' Error(75)> Time-sorted ensemble averaged cospectra outputs not created.'
        case(76)
            write(*,*) ' Error(76)> EddyPro could not calculate ensemble spectra.'
            write(*,*) ' Error(76)> Spectral assessment failed. Spectral assessment file not created.'
        case(77)
            write(*,*) ' Error(77)> EddyPro could not calculate ensemble spectra.'
            write(*,*) ' Error(77)> Ensemble averaged spectral output not created.'
        case(78)
            write(*,*) ' Fatal error(78)> No files matching the expected template were found in the selected folder.'
            write(*,*) ' Fatal error(78)> Program execution aborted.'
            stop 1
        case(79)
            write(*,*) ' Error(79)> Inconsistent variable labels or units in biomet files.'
            write(*,*) ' Error(79)> EddyPro cannot resolve the conflict and will thus proceed without using biomet data.'
        case(80)
            write(*,*) ' Warning(80)> Implausible altitude value detected. Altitude defaulted to zero.'
        case(81)
            write(*,*) ' Warning(81)> Implausible latitude value detected. Latitude defaulted to zero.'
        case(82)
            write(*,*) ' Warning(82)> Implausible longitude value detected. Longitude defaulted to zero.'
        case(83)
            write(*,*) ' Warning(83)> Implausible canopy height value detected. Canopy height defaulted to zero.'
        case(84)
            write(*,*) ' Warning(84)> Implausible displacement height value detected.'
            write(*,*) ' Warning(84)> Displacement height defaulted to 0.67 times the canopy height.'
        case(85)
            write(*,*) ' Warning(85)> Implausible roughness length value detected.'
            write(*,*) ' Warning(85)> Roughness length  defaulted to 0.15 times the canopy height,'
            write(*,*) ' Warning(85)> or to 1.0mm if canopy height is zero.'
        case(86)
            write(*,*) ' Fatal error(86)> Could not retrieve files from directory. Either directory does not exist'
            write(*,*) ' Fatal error(86)> or it does not contain files matching the selected requirements.'
            write(*,*) ' Fatal error(86)> Program execution aborted.'
            stop 1
        case(87)
            write(*,*) '  Error(87)> Entered or inferred "Binned co-spectra files directory" does not exist.'
            write(*,*) '  Error(87)> EddyPro cannot perform spectral assessment, calculate ensemble averaged spectra'
            write(*,*) '  Error(87)> or calculate ensemble averaged co-spectra.'
            write(*,*) '  Error(87)> Continuing by switching to Moncrieff et al. (1997) spectral corrections'
            write(*,*) '  Error(87)> if an in-situ method was selected, and ignoring selection of ensemble'
            write(*,*) '  Error(87)> averaged spectra or co-spectra outputs'
        case(88)
            write(*,*) '  Error(88)> Entered or inferred "Full co-spectra files directory" does not exist.'
            write(*,*) '  Error(88)> EddyPro cannot use spectral correction method of Fratini et al. (2012)'
            write(*,*) '  Error(88)> Continuing by switching to Moncrieff et al. (1997)'
        case(89)
            write(*,*) '  Error(89)> Entered or inferred "Binned co-spectra files directory" does not contain any valid files.'
            write(*,*) '  Error(89)> EddyPro cannot perform spectral assessment, calculate ensemble averaged spectra'
            write(*,*) '  Error(89)> or calculate ensemble averaged co-spectra.'
            write(*,*) '  Error(89)> Continuing by switching to Moncrieff et al. (1997) spectral corrections'
            write(*,*) '  Error(89)> if an in-situ method was selected, and ignoring selection of ensemble'
            write(*,*) '  Error(89)> averaged spectra or co-spectra outputs'
        case(90)
            write(*,*) '  Error(90)> Entered or inferred "Binned co-spectra files directory" does not contain'
            write(*,*) '  Error(90)> any files corresponding to the selected start/end period.'
            write(*,*) '  Error(90)> EddyPro cannot perform spectral assessment, calculate ensemble averaged spectra'
            write(*,*) '  Error(90)> or calculate ensemble averaged co-spectra.'
            write(*,*) '  Error(90)> EddyPro will continue, switching to Moncrieff et al. (1997) spectral corrections'
            write(*,*) '  Error(90)> if an in-situ method was selected, and ignoring selection of ensemble'
            write(*,*) '  Error(90)> averaged spectra or co-spectra outputs.'
        case(91)
            write(*,*) '  Alert(91)> Time constant of linear detrending cannot be larger than flux averaging interval.'
            write(*,*) '  Alert(91)> Automatically set time constant equal to flux averaging interval.'
        case(92)
            write(*,*) ' Error(92)> Occurred while writing passive gases ensemble spectra on output file.'
            write(*,*) ' Error(92)> File not created.'
        case(93)
            write(*,*) ' Error(93)> Embedded biomet selected, but not processing GHG files. This is not possible.'
            write(*,*) ' Error(93)> EddyPro will proceed without using biomet data.'
        case(94)
            write(*,*) ' Warning(94)> The selected or inferred angle-of-attack correction method is not applicable'
            write(*,*) ' Warning(94)> to data collected with selected sonic anemometer.'
            write(*,*) ' Warning(94)> Continuing without applying any angle-of-attack correction.'
        case(95)
            write(*,*) ' Warning(95)> The selected "w-boost" correction is not applicable'
            write(*,*) ' Warning(95)> to data collected with selected sonic anemometer.'
            write(*,*) ' Warning(95)> Continuing without applying "w-boost" correction.'
    end select
end subroutine ExceptionHandler
