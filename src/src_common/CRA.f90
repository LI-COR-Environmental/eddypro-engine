subroutine CRA(Set, nrow, ncol, fs, tconst, fcol)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    integer, intent(in) :: ncol
    integer, intent(in) :: tconst
    integer, intent(in) :: fcol
    real(kind = dbl), intent(in) :: fs
    real(kind = dbl), intent(inout) :: Set(nrow, ncol)
    !> local variables
    integer :: i
    integer :: np
    integer :: nnp
    real(kind = dbl) tmp(nrow)


    np = nint(fs * tconst)

    !> Centered average for first part of array
    do i = 1, np/2
        nnp = (i-1)*2
        call AverageNoError(Set(i-nnp/2:i+nnp/2, fcol), nnp+1, 1, tmp(i), error)
    end do
    !> Centered average for inner points
    do i = np/2+1, nrow-np/2
        call AverageNoError(Set(i-np/2:i+np/2, fcol), np+1, 1, tmp(i), error)
    end do
    !> Centered average for last part of array
    do i = nrow-np/2+1, nrow
        nnp = (nrow-i)*2
        call AverageNoError(Set(i-nnp/2:i+nnp/2, fcol), nnp+1, 1, tmp(i), error)
    end do

    Set(:, fcol) = tmp
end subroutine CRA
