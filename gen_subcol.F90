! ****************************************************************************
! Generate stochastic subcolumns of cloud and precipitation condensate, given
! gridbox-mean values as would be specified by a global climate model. This
! subroutine calls a number of other subroutines, and is the "master" entry
! point for my generalized subcolumn generator.
!
! Inputs:
!   npts        Number of points (geolocation)
!   ncol        Number of subcolumns to generate
!   nlev        Number of vertical levels (height)
!   cf          Gridbox-average cloud fraction
!   pf          Gridbox-average precipitation fraction
!   qc_mean     Gridbox-average of IN-CLOUD cloud liquid condensate
!   qi_mean     Gridbox-average of IN-CLOUD cloud ice condensate
!   qpc_mean    Gridbox-average of IN-PRECIP precipitating liquid condensate
!   qpi_mean    Gridbox-average of IN-PRECIP precipitating ice condensate
!   qc_var      Gridbox-variance of IN-CLOUD subcolumn cloud liquid condensate
!   qi_var      Gridbox-variance of IN-CLOUD subcolumn cloud ice condensate
!   qpc_var     Gridbox-variance of IN-PRECIP precipitating liquid condensate
!   qpi_var     Gridbox-variance of IN-PRECIP precipitating ice condensate
!   alphac      Overlap parameter for cloud occurrence
!   alphap      Overlap parameter for precipitation occurrence
!               (not currently used)
!   qc_rho      Rank correlation between cloudy parts of adjacent layers
!   qi_rho      Rank correlation between cloudy parts of adjacent layers
!   qpc_rho     Rank correlation between precipitating parts of adjacent layers
!   qpi_rho     Rank correlation between precipitating parts of adjacent layers
!   padj        Flag to adjust generated precipitation by precipitation 
!               fracion
!   scale_means Flag to scale the gridbox-averaged in-cloud and in-precip
!               condensate amounts by the new (generated) cloud and precip
!               fractions in order to reproduce the desired gridbox-average
!               condensate (as opposed to reproducing the in-cloud or in-
!               precip averages)
!   seed        Seed for random number generator. If negative, will seed with
!               system time from operating system, otherwise seeds with value
!               given. Useful for repeating experiments with same psuedo-
!               random numbers.
!
! Outputs:
!   qc_out      Subcolumn cloud liquid condensate
!   qi_out      Subcolumn cloud ice condensate
!   qpc_out     Subcolumn precipitating liquid condensate
!   qpi_out     Subcolumn precipitating ice condensate
! ****************************************************************************
subroutine gen_subcol(npts, ncol, nlev, cf, pf, &
                      qc_mean, qi_mean, qpc_mean, qpi_mean, &
                      qc_var, qi_var, qpc_var, qpi_var, &
                      alphac, alphap, &
                      qc_rho, qi_rho, qpc_rho, qpi_rho, &
                      qc_out, qi_out, qpc_out, qpi_out, &
                      padj, scale_means, seed, verbose)
    !f2py integer, intent(in) :: npts, ncol, nlev
    !f2py real, intent(in) :: cf(npts, nlev), pf(npts, nlev)
    !f2py real, intent(in) :: qc_mean(npts, nlev), qi_mean(npts, nlev)
    !f2py real, intent(in) :: qpc_mean(npts, nlev), qpi_mean(npts, nlev)
    !f2py real, intent(in) :: qc_var(npts, nlev), qi_var(npts, nlev)
    !f2py real, intent(in) :: qpc_var(npts, nlev), qpi_var(npts, nlev)
    !f2py real, intent(in) :: alphac(npts, nlev), alphap(npts, nlev)
    !f2py real, intent(in) :: qc_rho(npts, nlev), qi_rho(npts, nlev)
    !f2py real, intent(in) :: qpc_rho(npts, nlev), qpi_rho(npts, nlev)
    !f2py integer, optional, intent(in) :: padj=0, seed=-1, scale_means=1,
    !f2py integer, optional, intent(in) :: verbose=0
    !f2py real, intent(out) :: qc_out(npts, ncol, nlev)
    !f2py real, intent(out) :: qi_out(npts,ncol,nlev)
    !f2py real, intent(out) :: qpc_out(npts, ncol, nlev)
    !f2py real, intent(out) :: qpi_out(npts,ncol,nlev)

    implicit none

    ! inputs
    integer, intent(in) :: npts, ncol, nlev
    real, intent(in), dimension(npts, nlev) :: cf
    real, intent(in), dimension(npts, nlev) :: pf
    real, intent(in), dimension(npts, nlev) :: qc_mean
    real, intent(in), dimension(npts, nlev) :: qi_mean
    real, intent(in), dimension(npts, nlev) :: qpc_mean
    real, intent(in), dimension(npts, nlev) :: qpi_mean
    real, intent(in), dimension(npts, nlev) :: qc_var
    real, intent(in), dimension(npts, nlev) :: qi_var
    real, intent(in), dimension(npts, nlev) :: qpc_var
    real, intent(in), dimension(npts, nlev) :: qpi_var
    real, intent(in), dimension(npts, nlev) :: alphac   
    real, intent(in), dimension(npts, nlev) :: alphap   
    real, intent(in), dimension(npts, nlev) :: qc_rho, qi_rho  
    real, intent(in), dimension(npts, nlev) :: qpc_rho, qpi_rho  
    integer, intent(in) :: padj, seed, scale_means, verbose

    ! outputs
    real, intent(out), dimension(npts, ncol, nlev) :: qc_out
    real, intent(out), dimension(npts, ncol, nlev) :: qi_out
    real, intent(out), dimension(npts, ncol, nlev) :: qpc_out
    real, intent(out), dimension(npts, ncol, nlev) :: qpi_out

    ! local namespace
    real :: fracout(npts, ncol, nlev), cf_new(npts, nlev)
    real :: precout(npts, ncol, nlev), pf_new(npts, nlev)
    real :: qc_scaled(npts, nlev), qi_scaled(npts, nlev)
    real :: qpc_scaled(npts, nlev), qpi_scaled(npts, nlev)

    ! seed random number generator
    if (seed /= 0) then
        if (verbose > 0) print *, 'Initializing random number from gen_subcol &
                                 &with', seed
        call init_random_seed(seed)
    end if

    ! generate subcolumn cloud occurrence
    if (all(alphac < 0)) then
        if (verbose > 0) print *, 'Generating subcolumn cloud occurrence &
                                  &using SCOPS'
        call scops(npts, nlev, ncol, cf, 0.0*cf, 3, fracout, 0, 0)
    else
        if (verbose > 0) print *, 'Generating subcolumn cloud occurrence &
                                  &using R04'
        call gen_subcol_cld(npts, ncol, nlev, cf, alphac, fracout, 0)
    end if

    ! generate subcolumn precip occurrence
    if (all(alphap < 0)) then
        if (verbose > 0) print *, 'Generating subcolumn precip occurrence &
                                  &using PREC_SCOPS'
        call prec_scops(npts, nlev, ncol, pf, 0.0*pf, fracout, precout)
    else
        if (verbose > 0) print *, 'Generating subcolumn precip occurrence &
                                  &using R04'
        if (verbose > 0) print *, 'WARNING: this has not been tested &
                                  &extensively, and may not produce &
                                  &desirable overlap with cloud!'
        call gen_subcol_cld(npts, ncol, nlev, pf, alphap, precout, 0)
    end if

    ! constrain subcolumn precip occurrence
    if (padj > 0) then
        if (verbose > 0) print *, 'Adjusting precipitation occurrence'
        call adjust_precip(npts, ncol, nlev, pf, fracout, precout, padj, 0)
    end if

    ! scale the gridbox-mean condensate amounts by new occurrence fractions?
    if (scale_means > 0) then
        ! calculate new cloud and precipitation fractions
        cf_new = sum(fracout, dim=2) / ncol
        pf_new = sum(precout, dim=2) / ncol

        ! scale means
        qc_scaled = qc_mean * cf / cf_new
        qi_scaled = qi_mean * cf / cf_new
        qpc_scaled = qpc_mean * pf / pf_new
        qpi_scaled = qpi_mean * pf / pf_new
    else
        qc_scaled = qc_mean
        qi_scaled = qi_mean
        qpc_scaled = qpc_mean
        qpi_scaled = qpi_mean
    end if

    ! generate subcolumn cloud condensate
    if (all(qc_var <= 0) .and. all(qi_var <= 0)) then
        if (verbose > 0) print *, 'Generating homogeneous cloud condensate'
        ! homogeneous condensate
        where (fracout > 0)
            qc_out = spread(qc_scaled, 2, ncol)
            qi_out = spread(qi_scaled, 2, ncol)
        elsewhere
            qc_out = 0.0
            qi_out = 0.0
        endwhere
    else
        ! heterogeneous condensate following R04
        if (verbose > 0) print *, 'Generating heterogeneous cloud condensate'
        call gen_subcol_var( &
            npts, ncol, nlev, fracout, &
            qc_scaled, qc_var, qc_rho, qc_out, 0 &
        )
        call gen_subcol_var( &
            npts, ncol, nlev, fracout, &
            qi_scaled, qi_var, qi_rho, qi_out, 0 &
        )
    end if

    ! generate subcolumn precip condensate
    if (all(qpi_var <= 0) .and. all(qpi_var <= 0)) then
        ! homogeneous condensate
        if (verbose > 0) print *, 'Generating homogeneous precipitation &
                                  &condensate'
        where (precout > 0)
            qpc_out = spread(qpc_scaled, 2, ncol)
            qpi_out = spread(qpi_scaled, 2, ncol)
        elsewhere
            qpc_out = 0.0
            qpi_out = 0.0
        endwhere
    else
        ! heterogeneous condensate following R04
        if (verbose > 0) print *, 'Generating heterogeneous precipitation &
                                  &condensate using R04'
        call gen_subcol_var( &
            npts, ncol, nlev, precout, &
            qpc_scaled, qpc_var, qpc_rho, qpc_out, 0 & 
        )
        call gen_subcol_var( &
            npts, ncol, nlev, precout, &
            qpi_scaled, qpi_var, qpi_rho, qpi_out, 0 & 
        )
    end if
end subroutine
