include '../src_common/interfaces_1.inc'

interface
    subroutine BandPassSpectralCorrections(measuring_height, displ_height,&
        loc_var_present, wind_speed, t_air, zL, ac_frequency, avrg_length, logger_sw_ver,&
        detrending_method, detrending_time_constant, printout, LocInstr, &
        nfull, LocFileList, nrow_full, lEx, LocSetup)
        use m_common_global_var
        implicit none
        !> In/out variables
        integer, intent(in) :: avrg_length
        integer, intent(in) :: detrending_time_constant
        integer, intent(in) :: nfull
        real(kind = dbl), intent(in) :: measuring_height
        real(kind = dbl), intent(in) :: displ_height
        logical, intent(in) :: loc_var_present(GHGNumVar)
        real(kind = dbl), intent(in) :: wind_speed
        real(kind = dbl), intent(in) :: t_air
        real(kind = dbl), intent(in) :: zL
        real(kind = dbl), intent(in) :: ac_frequency
        character(2), intent(in) :: detrending_method
        logical, intent(in) :: printout
        type(InstrumentType), intent(in) :: LocInstr(GHGNumVar)
        type(SwVerType) :: logger_sw_ver
        !> Optional variables
        integer, optional, intent(in):: nrow_full
        type(FileListType), optional, intent(in) :: LocFileList(nfull)
        type(ExType), optional, intent(in) :: lEx
        type(FCCsetupType), optional, intent(in) :: LocSetup
    end subroutine BandPassSpectralCorrections
end interface