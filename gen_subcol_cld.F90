! ****************************************************************************
! Generate stochastic cloudy/clear subcolumns using algorithm described in
! Raisanen et al. (2004).
!
! Inputs:
!     cf(npts, nlev)        cloud fraction of layers
!
! Outputs:
!     cb(npts, ncol, nlev)  cloudy/clear flag (1 -> cloud, 0 -> clear)
!
! NOTE:   Inputs MUST be from TOA to surface!
! AUTHOR: Benjamin R. Hillman
! DATE:   12 March 2016
! ****************************************************************************
subroutine gen_subcol_cld(npts, ncol, nlev, cf, alpha, cb, seed)
    !f2py integer, intent(in) :: npts, ncol, nlev
    !f2py real, intent(in) :: cf(npts, nlev), alpha(npts,nlev)
    !f2py real, intent(out) :: cb(npts, ncol, nlev)

    implicit none

    ! inputs
    integer, intent(in) :: npts, ncol, nlev
    real, intent(in) :: cf(npts,nlev)  ! cloud fraction
    real, intent(in) :: alpha(npts,nlev)  ! overlap parameter
    integer, intent(in) :: seed  ! seed for random number generator

    ! outputs
    real, intent(out) :: cb(npts,ncol,nlev)  ! subcolumn cloudy/clear flag
      
    ! local variables
    real, dimension(npts,ncol,nlev) :: x, rn1, rn2, rn3
    integer :: i, j, k

    ! initialize psuedo-random number generator
    if (seed /= 0) then
        print *, 'Initializing random seed from gen_subcol_cld'
        call init_random_seed(seed)
    end if

    ! generate random numbers
    call random_number(rn1)
    call random_number(rn2)
    call random_number(rn3)

    ! calculate x(i, k)
    x(:, :, 1) = rn1(:, :, 1)
    do k = 2, nlev
        where (rn2(:, :, k) <= spread(alpha(:, k), 2, ncol))
            x(:, :, k) = x(:, :, k - 1)
        elsewhere
            x(:, :, k) = rn3(:, :, k)
        endwhere
    end do

    ! cloudy binary flag
    where (x > (1.0 - spread(cf, 2, ncol)))
        cb = 1.0
    elsewhere
        cb = 0.0
    endwhere
end subroutine
