! Inputs
!     npts
!     ncol
!     nlev
!     padj      Type of precipitation adjustment to perform;
!               1 -> Remove from least cloudy, add to most cloudy columns
!               2 -> Remove / add precip from randomly selected columns
!               3 -> Remove / add precip systematically
!     pf
!     cb
!         
subroutine adjust_precip(npts, ncol, nlev, pf, cb, pb, padj, seed)
    !f2py integer, intent(in) :: npts, ncol, nlev
    !f2py real, intent(in) :: pf(npts, nlev), cb(npts, ncol, nlev)
    !f2py real, intent(inout) :: pb(npts, ncol, nlev)
    !f2py integer, intent(in) :: padj=1
    !f2py integer, intent(in) :: seed=-1

    use sort
    implicit none

    ! inputs
    integer, intent(in) :: npts, ncol, nlev
    real, intent(in) :: pf(npts, nlev)
    real, intent(in) :: cb(npts, ncol, nlev)
    real, intent(inout) :: pb(npts, ncol, nlev)
    integer, intent(in) :: padj  ! type of adjustment
    integer, intent(in) :: seed

    ! local namespace
    integer :: i, j, k  ! loop indices
    real :: cn(npts, ncol), x
    integer :: ci(npts, ncol)

    ! seed random number generator
    if (seed /= 0) then
        call init_random_seed(seed)
    end if

    ! get number of cloudy levels in each subcolumn
    cn = sum(cb, 3)

    ! loop over points and levels
    if (padj == 1) then
        do i = 1, npts
            ! sort cloudy counts; ci returns indices to original array, so if 
            ! we want the least cloudy column, take ci(1), if we want the most 
            ! cloudy column, take ci(ncol), etc.
            call quick_sort(cn(i, :), ci(i, :))

            ! loop over levels
            do k = 1, nlev
                ! remove precip if we have too much
                j = 1
                do while(sum(pb(i, :, k)) / ncol > pf(i, k))
                    ! remove precip from LEAST cloudy columns first!
                    pb(i, ci(i, j), k) = 0.0
                    j = j + 1
                end do

                ! add precip if we have too little
                j = ncol
                do while(sum(pb(i, :, k)) / ncol < pf(i, k))
                    ! add precip from MOST cloudy columns first!
                    pb(i, ci(i, j), k) = 1.0
                    j = j - 1
                end do
            end do
        end do
    else if (padj == 2) then
        do i = 1, npts
            ! loop over levels
            do k = 1, nlev
                ! remove precip if we have too much
                do while(sum(pb(i, :, k)) / ncol > pf(i, k))
                    call random_number(x)
                    j = int(x * ncol) + 1
                    pb(i, j, k) = 0.0
                end do

                ! add precip if we have too little
                do while(sum(pb(i, :, k)) / ncol < pf(i, k))
                    call random_number(x)
                    j = int(x * ncol) + 1
                    pb(i, j, k) = 1.0
                end do
            end do
        end do 
    else if (padj == 3) then
        do i = 1, npts
            ! remove precip systematically
            do k = 1, nlev
                ! remove precip if we have too much
                j = 1
                do while(sum(pb(i, :, k)) / ncol > pf(i, k))
                    pb(i, j, k) = 0.0
                    j = j + 1
                end do

                ! add precip if we have too little
                j = ncol
                do while(sum(pb(i, :, k)) / ncol < pf(i, k))
                    pb(i, j, k) = 1.0
                    j = j - 1
                end do
            end do
        end do
    else
        print *, 'ERROR: padj flag value not recognized!'
        stop
    end if 
end subroutine
