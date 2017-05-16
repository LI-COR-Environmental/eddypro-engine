interface
    subroutine fcn(m, npar, x, fvec, fjac, iflag)
        implicit none
        integer, parameter :: dbl   = kind(0.0d0)
        integer, intent(in)            :: m, npar
        real(kind = dbl), intent(in)    :: x(:)
        real(kind = dbl), intent(inout) :: fvec(:)
        real(kind = dbl), intent(out)   :: fjac(:,:)
        integer, intent(inout)         :: iflag
    end subroutine fcn
end interface

interface
    function SwVerFromString(string)
        use m_common_global_var
        implicit none
        !> in/out variables
        character(*), intent(in) :: string
        type(SwVerType) :: SwVerFromString
        !> Local variables
        character(4) :: chunk
    end function
end interface

interface
    function CompareSwVer(v1, v2)
        use m_common_global_var
        implicit none
        !> In/out variables
        type(SwVerType), intent(in) :: v1
        type(SwVerType), intent(in) :: v2
        logical :: CompareSwVer
    end function
end interface

interface
    function EqualSwVer(v1, v2)
        use m_common_global_var
        implicit none
        !> In/out variables
        type(SwVerType), intent(in) :: v1
        type(SwVerType), intent(in) :: v2
        logical :: EqualSwVer
    end function
end interface

interface
    function DefaultSonicOutputRate(model)
        use m_common_global_var
        implicit none
        !> In/out variables
        character(*), intent(in) :: model
        integer :: DefaultSonicOutputRate
    end function
end interface

interface
    function sinc(x, N)
        use m_numeric_kinds
        implicit none
        !> In/out variables
        integer, intent(in) :: N
        real(kind = dbl) , intent(in) :: x(N)
        real(kind = dbl)  :: sinc(N)
    end function
end interface

interface
    function asymmetric_linear_transformation(x, N, pgain, poffset, ngain, noffset)
        use m_numeric_kinds
        implicit none
        !> In/out variables
        integer, intent(in) :: N
        real(kind = dbl) , intent(inout) :: x(N)
        real(kind = dbl) , intent(in) :: pgain
        real(kind = dbl) , intent(in) :: poffset
        real(kind = dbl) , intent(in) :: ngain
        real(kind = dbl) , intent(in) :: noffset
        real(kind = dbl)  :: asymmetric_linear_transformation(N)
    end function
end interface

interface
    function CountRecordsAndValues(Set, nrow, ncol, var1, var2)
        use m_common_global_var
        implicit none
        !> In/out variables
        integer, intent(in) :: nrow, ncol
        real(kind = dbl), intent(in) :: Set(nrow, ncol)
        integer, optional, intent(in) :: var1, var2
        integer :: CountRecordsAndValues
    end function
end interface