module QUADPACK
    use integrationDataTypes
    implicit none
    
    
    
    contains
    
    
    
    subroutine aaaa

!*****************************************************************************80
!
!! AAAA is a dummy subroutine with QUADPACK documentation in its comments.
!
! 1. introduction
!
!    quadpack is a fortran subroutine package for the numerical
!    computation of definite one-dimensional integrals. it originated
!    from a joint project of r. piessens and e. de doncker (appl.
!    math. and progr. div.- k.u.leuven, belgium), c. ueberhuber (inst.
!    fuer math.- techn.u.wien, austria), and d. kahaner (nation. bur.
!    of standards- washington d.c., u.s.a.).
!
! 2. survey
!
!    - qags : is an integrator based on globally adaptive interval
!             subdivision in connection with extrapolation (de doncker,
!             1978) by the epsilon algorithm (wynn, 1956).
!
!    - qagp : serves the same purposes as qags, but also allows
!             for eventual user-supplied information, i.e. the
!             abscissae of internal singularities, discontinuities
!             and other difficulties of the integrand function.
!             the algorithm is a modification of that in qags.
!
!    - qagi : handles integration over infinite intervals. the
!             infinite range is mapped onto a finite interval and
!             then the same strategy as in qags is applied.
!
!    - qawo : is a routine for the integration of cos(omega*x)*f(x)
!             or sin(omega*x)*f(x) over a finite interval (a,b).
!             omega is is specified by the user
!             the rule evaluation component is based on the
!             modified clenshaw-curtis technique.
!             an adaptive subdivision scheme is used connected with
!             an extrapolation procedure, which is a modification
!             of that in qags and provides the possibility to deal
!             even with singularities in f.
!
!    - qawf : calculates the fourier cosine or fourier sine
!             transform of f(x), for user-supplied interval (a,
!             infinity), omega, and f. the procedure of qawo is
!             used on successive finite intervals, and convergence
!             acceleration by means of the epsilon algorithm (wynn,
!             1956) is applied to the series of the integral
!             contributions.
!
!    - qaws : integrates w(x)*f(x) over (a,b) with a < b finite,
!             and   w(x) = ((x-a)**alfa)*((b-x)**beta)*v(x)
!             where v(x) = 1 or log(x-a) or log(b-x)
!                            or log(x-a)*log(b-x)
!             and   -1 < alfa, -1 < beta.
!             the user specifies a, b, alfa, beta and the type of
!             the function v.
!             a globally adaptive subdivision strategy is applied,
!             with modified clenshaw-curtis integration on the
!             subintervals which contain a or b.
!
!    - qawc : computes the cauchy principal value of f(x)/(x-c)
!             over a finite interval (a,b) and for
!             user-determined c.
!             the strategy is globally adaptive, and modified
!             clenshaw-curtis integration is used on the subranges
!             which contain the point x = c.
!
!  each of the routines above also has a "more detailed" version
!    with a name ending in e, as qage.  these provide more
!    information and control than the easier versions.
!
!
!   the preceeding routines are all automatic.  that is, the user
!      inputs his problem and an error tolerance.  the routine
!      attempts to perform the integration to within the requested
!      absolute or relative error.
!   there are, in addition, a number of non-automatic integrators.
!      these are most useful when the problem is such that the
!      user knows that a fixed rule will provide the accuracy
!      required.  typically they return an error estimate but make
!      no attempt to satisfy any particular input error request.
!
!      qk15
!      qk21_x
!      qk31
!      qk41
!      qk51
!      qk61
!           estimate the integral on [a,b] using 15, 21,..., 61
!           point rule and return an error estimate.
!      qk15i 15 point rule for (semi)infinite interval.
!      qk15w 15 point rule for special singular weight functions.
!      qc25c 25 point rule for cauchy principal values
!      qc25o 25 point rule for sin/cos integrand.
!      qmomo integrates k-th degree chebychev polynomial times
!            function with various explicit singularities.
!
! 3. guidelines for the use of quadpack
!
!    here it is not our purpose to investigate the question when
!    automatic quadrature should be used. we shall rather attempt
!    to help the user who already made the decision to use quadpack,
!    with selecting an appropriate routine or a combination of
!    several routines for handling his problem.
!
!    for both quadrature over finite and over infinite intervals,
!    one of the first questions to be answered by the user is
!    related to the amount of computer time he wants to spend,
!    versus his -own- time which would be needed, for example, for
!    manual subdivision of the interval or other analytic
!    manipulations.
!
!    (1) the user may not care about computer time, or not be
!        willing to do any analysis of the problem. especially when
!        only one or a few integrals must be calculated, this attitude
!        can be perfectly reasonable. in this case it is clear that
!        either the most sophisticated of the routines for finite
!        intervals, qags, must be used, or its analogue for infinite
!        intervals, qagi. these routines are able to cope with
!        rather difficult, even with improper integrals.
!        this way of proceeding may be expensive. but the integrator
!        is supposed to give you an answer in return, with additional
!        information in the case of a failure, through its error
!        estimate and flag. yet it must be stressed that the programs
!        cannot be totally reliable.
!
!    (2) the user may want to examine the integrand function.
!        if bad local difficulties occur, such as a discontinuity, a
!        singularity, derivative singularity or high peak at one or
!        more points within the interval, the first advice is to
!        split up the interval at these points. the integrand must
!        then be examinated over each of the subintervals separately,
!        so that a suitable integrator can be selected for each of
!        them. if this yields problems involving relative accuracies
!        to be imposed on -finite- subintervals, one can make use of
!        qagp, which must be provided with the positions of the local
!        difficulties. however, if strong singularities are present
!        and a high accuracy is requested, application of qags on the
!        subintervals may yield a better result.
!
!        for quadrature over finite intervals we thus dispose of qags
!        and
!        - qng for well-behaved integrands,
!        - qag for functions with an oscillating behavior of a non
!          specific type,
!        - qawo for functions, eventually singular, containing a
!          factor cos(omega*x) or sin(omega*x) where omega is known,
!        - qaws for integrands with algebraico-logarithmic end point
!          singularities of known type,
!        - qawc for cauchy principal values.
!
!        remark
!
!        on return, the work arrays in the argument lists of the
!        adaptive integrators contain information about the interval
!        subdivision process and hence about the integrand behavior:
!        the end points of the subintervals, the local integral
!        contributions and error estimates, and eventually other
!        characteristics. for this reason, and because of its simple
!        globally adaptive nature, the routine qag in particular is
!        well-suited for integrand examination. difficult spots can
!        be located by investigating the error estimates on the
!        subintervals.
!
!        for infinite intervals we provide only one general-purpose
!        routine, qagi. it is based on the qags algorithm applied
!        after a transformation of the original interval into (0,1).
!        yet it may eventuate that another type of transformation is
!        more appropriate, or one might prefer to break up the
!        original interval and use qagi only on the infinite part
!        and so on. these kinds of actions suggest a combined use of
!        different quadpack integrators. note that, when the only
!        difficulty is an integrand singularity at the finite
!        integration limit, it will in general not be necessary to
!        break up the interval, as qagi deals with several types of
!        singularity at the boundary point of the integration range.
!        it also handles slowly convergent improper integrals, on
!        the condition that the integrand does not oscillate over
!        the entire infinite interval. if it does we would advise
!        to sum succeeding positive and negative contributions to
!        the integral -e.g. integrate between the zeros- with one
!        or more of the finite-range integrators, and apply
!        convergence acceleration eventually by means of quadpack
!        subroutine qelg which implements the epsilon algorithm.
!        such quadrature problems include the fourier transform as
!        a special case. yet for the latter we have an automatic
!        integrator available, qawf.
!
  return
end subroutine
subroutine qag ( f_ptr, dat, a, b, epsabs, epsrel, key, result, abserr, neval, ier )

!*****************************************************************************80
!
!! QAG approximates an integral over a finite interval.
!
!  Discussion:
!
!    The routine calculates an approximation RESULT to a definite integral   
!      I = integral of F over (A,B),
!    hopefully satisfying
!      || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
!
!    QAG is a simple globally adaptive integrator using the strategy of 
!    Aind (Piessens, 1973).  It is possible to choose between 6 pairs of
!    Gauss-Kronrod quadrature formulae for the rule evaluation component. 
!    The pairs of high degree of precision are suitable for handling
!    integration difficulties due to a strongly oscillating integrand.
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real ( kind = 8 ) F, the name of the function routine, of the form
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!    which evaluates the integrand function.
!
!    Input, real ( kind = 8 ) A, B, the limits of integration.
!
!    Input, real ( kind = 8 ) EPSABS, EPSREL, the absolute and relative accuracy requested.
!
!    Input, integer ( kind = 8 ) KEY, chooses the order of the local integration rule:
!    1,  7 Gauss points, 15 Gauss-Kronrod points,
!    2, 10 Gauss points, 21 Gauss-Kronrod points,
!    3, 15 Gauss points, 31 Gauss-Kronrod points,
!    4, 20 Gauss points, 41 Gauss-Kronrod points,
!    5, 25 Gauss points, 51 Gauss-Kronrod points,
!    6, 30 Gauss points, 61 Gauss-Kronrod points.
!
!    Output, real ( kind = 8 ) RESULT, the estimated value of the integral.
!
!    Output, real ( kind = 8 ) ABSERR, an estimate of || I - RESULT ||.
!
!    Output, integer ( kind = 8 ) NEVAL, the number of times the integral was evaluated.
!
!    Output, integer ( kind = 8 ) IER, return code.
!    0, normal and reliable termination of the routine.  It is assumed that the 
!      requested accuracy has been achieved.
!    1, maximum number of subdivisions allowed has been achieved.  One can 
!      allow more subdivisions by increasing the value of LIMIT in QAG. 
!      However, if this yields no improvement it is advised to analyze the
!      integrand to determine the integration difficulties.  If the position
!      of a local difficulty can be determined, such as a singularity or
!      discontinuity within the interval) one will probably gain from 
!      splitting up the interval at this point and calling the integrator 
!      on the subranges.  If possible, an appropriate special-purpose 
!      integrator should be used which is designed for handling the type 
!      of difficulty involved.
!    2, the occurrence of roundoff error is detected, which prevents the
!      requested tolerance from being achieved.
!    3, extremely bad integrand behavior occurs at some points of the
!      integration interval.
!    6, the input is invalid, because EPSABS < 0 and EPSREL < 0.
!
!  Local parameters:
!
!    LIMIT is the maximum number of subintervals allowed in
!    the subdivision process of QAGE.
!
  implicit none

  integer ( kind = 8 ), parameter :: limit = 500

  real ( kind = 8 ) a
  real ( kind = 8 ) abserr
  real ( kind = 8 ) alist(limit)
  real ( kind = 8 ) b
  real ( kind = 8 ) blist(limit)
  real ( kind = 8 ) elist(limit)
  real ( kind = 8 ) epsabs
  real ( kind = 8 ) epsrel
  !real ( kind = 8 ), external :: f
  integer ( kind = 8 ) ier
  integer ( kind = 8 ) iord(limit)
  integer ( kind = 8 ) key
  integer ( kind = 8 ) last
  integer ( kind = 8 ) neval
  real ( kind = 8 ) result
  real ( kind = 8 ) rlist(limit)
  class(dataCollectionBase), target :: dat
  procedure (f_int_dat), intent(in), pointer :: f_ptr

  call qage ( f_ptr, dat, a, b, epsabs, epsrel, key, limit, result, abserr, neval, &
    ier, alist, blist, rlist, elist, iord, last )

  return
end subroutine
subroutine qage ( f_ptr, dat, a, b, epsabs, epsrel, key, limit, result, abserr, neval, &
  ier, alist, blist, rlist, elist, iord, last )

!*****************************************************************************80
!
!! QAGE estimates a definite integral.
!
!  Discussion:
!
!    The routine calculates an approximation RESULT to a definite integral   
!      I = integral of F over (A,B),
!    hopefully satisfying
!      || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real ( kind = 8 ) F, the name of the function routine, of the form
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!    which evaluates the integrand function.
!
!    Input, real ( kind = 8 ) A, B, the limits of integration.
!
!    Input, real ( kind = 8 ) EPSABS, EPSREL, the absolute and relative accuracy requested.
!
!    Input, integer ( kind = 8 ) KEY, chooses the order of the local integration rule:
!    1,  7 Gauss points, 15 Gauss-Kronrod points,
!    2, 10 Gauss points, 21 Gauss-Kronrod points,
!    3, 15 Gauss points, 31 Gauss-Kronrod points,
!    4, 20 Gauss points, 41 Gauss-Kronrod points,
!    5, 25 Gauss points, 51 Gauss-Kronrod points,
!    6, 30 Gauss points, 61 Gauss-Kronrod points.
!
!    Input, integer ( kind = 8 ) LIMIT, the maximum number of subintervals that
!    can be used.
!
!    Output, real ( kind = 8 ) RESULT, the estimated value of the integral.
!
!    Output, real ( kind = 8 ) ABSERR, an estimate of || I - RESULT ||.
!
!    Output, integer ( kind = 8 ) NEVAL, the number of times the integral was evaluated.
!
!    Output, integer ( kind = 8 ) IER, return code.
!    0, normal and reliable termination of the routine.  It is assumed that the 
!      requested accuracy has been achieved.
!    1, maximum number of subdivisions allowed has been achieved.  One can 
!      allow more subdivisions by increasing the value of LIMIT in QAG. 
!      However, if this yields no improvement it is advised to analyze the
!      integrand to determine the integration difficulties.  If the position
!      of a local difficulty can be determined, such as a singularity or
!      discontinuity within the interval) one will probably gain from 
!      splitting up the interval at this point and calling the integrator 
!      on the subranges.  If possible, an appropriate special-purpose 
!      integrator should be used which is designed for handling the type 
!      of difficulty involved.
!    2, the occurrence of roundoff error is detected, which prevents the
!      requested tolerance from being achieved.
!    3, extremely bad integrand behavior occurs at some points of the
!      integration interval.
!    6, the input is invalid, because EPSABS < 0 and EPSREL < 0.
!
!    Workspace, real ( kind = 8 ) ALIST(LIMIT), BLIST(LIMIT), contains in entries 1 
!    through LAST the left and right ends of the partition subintervals.
!
!    Workspace, real ( kind = 8 ) RLIST(LIMIT), contains in entries 1 through LAST
!    the integral approximations on the subintervals.
!
!    Workspace, real ( kind = 8 ) ELIST(LIMIT), contains in entries 1 through LAST
!    the absolute error estimates on the subintervals.
!
!    Output, integer ( kind = 8 ) IORD(LIMIT), the first K elements of which are pointers 
!    to the error estimates over the subintervals, such that
!    elist(iord(1)), ..., elist(iord(k)) form a decreasing sequence, with
!    k = last if last <= (limit/2+2), and k = limit+1-last otherwise.
!
!    Output, integer ( kind = 8 ) LAST, the number of subintervals actually produced 
!    in the subdivision process.
!
!  Local parameters:
!
!    alist     - list of left end points of all subintervals
!                       considered up to now
!    blist     - list of right end points of all subintervals
!                       considered up to now
!    elist(i)  - error estimate applying to rlist(i)
!    maxerr    - pointer to the interval with largest error estimate
!    errmax    - elist(maxerr)
!    area      - sum of the integrals over the subintervals
!    errsum    - sum of the errors over the subintervals
!    errbnd    - requested accuracy max(epsabs,epsrel*abs(result))
!    *****1    - variable for the left subinterval
!    *****2    - variable for the right subinterval
!    last      - index for subdivision
!
  implicit none

  integer ( kind = 8 ) limit

  real ( kind = 8 ) a
  real ( kind = 8 ) abserr
  real ( kind = 8 ) alist(limit)
  real ( kind = 8 ) area
  real ( kind = 8 ) area1
  real ( kind = 8 ) area12
  real ( kind = 8 ) area2
  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  real ( kind = 8 ) b
  real ( kind = 8 ) blist(limit)
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) c
  real ( kind = 8 ) defabs
  real ( kind = 8 ) defab1
  real ( kind = 8 ) defab2
  real ( kind = 8 ) elist(limit)
  real ( kind = 8 ) epsabs
  real ( kind = 8 ) epsrel
  real ( kind = 8 ) errbnd
  real ( kind = 8 ) errmax
  real ( kind = 8 ) error1
  real ( kind = 8 ) error2
  real ( kind = 8 ) erro12
  real ( kind = 8 ) errsum
  !real ( kind = 8 ), external :: f
  procedure (f_int_dat), intent(in), pointer :: f_ptr
  integer ( kind = 8 ) ier
  integer ( kind = 8 ) iord(limit)
  integer ( kind = 8 ) iroff1
  integer ( kind = 8 ) iroff2
  integer ( kind = 8 ) key
  integer ( kind = 8 ) keyf
  integer ( kind = 8 ) last
  integer ( kind = 8 ) maxerr
  integer ( kind = 8 ) neval
  integer ( kind = 8 ) nrmax
  real ( kind = 8 ) resabs
  real ( kind = 8 ) result
  real ( kind = 8 ) rlist(limit)
  class(dataCollectionBase), target :: dat
!
!  Test on validity of parameters.
!
  ier = 0
  neval = 0
  last = 0
  result = 0.0E+00
  abserr = 0.0E+00
  alist(1) = a
  blist(1) = b
  rlist(1) = 0.0E+00
  elist(1) = 0.0E+00
  iord(1) = 0

  if ( epsabs < 0.0E+00 .and. epsrel < 0.0E+00 ) then
    ier = 6
    return
  end if
!
!  First approximation to the integral.
!
  keyf = key
  keyf = max ( keyf, 1 )
  keyf = min ( keyf, 6 )

  c = keyf
  neval = 0

  if ( keyf == 1 ) then
    call qk15 ( f_ptr, dat, a, b, result, abserr, defabs, resabs )
  else if ( keyf == 2 ) then
    call qk21_x ( f_ptr, dat, a, b, result, abserr, defabs, resabs )
  else if ( keyf == 3 ) then
    call qk31 ( f_ptr, dat, a, b, result, abserr, defabs, resabs )
  else if ( keyf == 4 ) then
    call qk41 ( f_ptr, dat, a, b, result, abserr, defabs, resabs )
  else if ( keyf == 5 ) then
    call qk51 ( f_ptr, dat, a, b, result, abserr, defabs, resabs )
  else if ( keyf == 6 ) then
    call qk61 ( f_ptr, dat, a, b, result, abserr, defabs, resabs )
  end if

  last = 1
  rlist(1) = result
  elist(1) = abserr
  iord(1) = 1
!
!  Test on accuracy.
!
  errbnd = max ( epsabs, epsrel * abs ( result ) )

  if ( abserr <= 5.0E+01 * epsilon ( defabs ) * defabs .and. &
    errbnd < abserr ) then
    ier = 2
  end if

  if ( limit == 1 ) then
    ier = 1
  end if

  if ( ier /= 0 .or. &
    ( abserr <= errbnd .and. abserr /= resabs ) .or. &
    abserr == 0.0E+00 ) then

    if ( keyf /= 1 ) then
      neval = (10*keyf+1) * (2*neval+1)
    else
      neval = 30 * neval + 15
    end if

    return

  end if
!
!  Initialization.
!
  errmax = abserr
  maxerr = 1
  area = result
  errsum = abserr
  nrmax = 1
  iroff1 = 0
  iroff2 = 0

  do last = 2, limit
!
!  Bisect the subinterval with the largest error estimate.
!
    a1 = alist(maxerr)
    b1 = 0.5E+00 * ( alist(maxerr) + blist(maxerr) )
    a2 = b1
    b2 = blist(maxerr)

    if ( keyf == 1 ) then
      call qk15 ( f_ptr, dat, a1, b1, area1, error1, resabs, defab1 )
    else if ( keyf == 2 ) then
      call qk21_x ( f_ptr, dat, a1, b1, area1, error1, resabs, defab1 )
    else if ( keyf == 3 ) then
      call qk31 ( f_ptr, dat, a1, b1, area1, error1, resabs, defab1 )
    else if ( keyf == 4 ) then
      call qk41 ( f_ptr, dat, a1, b1, area1, error1, resabs, defab1)
    else if ( keyf == 5 ) then
      call qk51 ( f_ptr, dat, a1, b1, area1, error1, resabs, defab1 )
    else if ( keyf == 6 ) then
      call qk61 ( f_ptr, dat, a1, b1, area1, error1, resabs, defab1 )
    end if

    if ( keyf == 1 ) then
      call qk15 ( f_ptr, dat, a2, b2, area2, error2, resabs, defab2 )
    else if ( keyf == 2 ) then
      call qk21_x ( f_ptr, dat, a2, b2, area2, error2, resabs, defab2 )
    else if ( keyf == 3 ) then
      call qk31 ( f_ptr, dat, a2, b2, area2, error2, resabs, defab2 )
    else if ( keyf == 4 ) then
      call qk41 ( f_ptr, dat, a2, b2, area2, error2, resabs, defab2 )
    else if ( keyf == 5 ) then
      call qk51 ( f_ptr, dat, a2, b2, area2, error2, resabs, defab2 )
    else if ( keyf == 6 ) then
      call qk61 ( f_ptr, dat, a2, b2, area2, error2, resabs, defab2 )
    end if
!
!  Improve previous approximations to integral and error and
!  test for accuracy.
!
    neval = neval + 1
    area12 = area1 + area2
    erro12 = error1 + error2
    errsum = errsum + erro12 - errmax
    area = area + area12 - rlist(maxerr)

    if ( defab1 /= error1 .and. defab2 /= error2 ) then

      if ( abs ( rlist(maxerr) - area12 ) <= 1.0E-05 * abs ( area12 ) &
        .and. 9.9E-01 * errmax <= erro12 ) then
        iroff1 = iroff1 + 1
      end if

      if ( 10 < last .and. errmax < erro12 ) then
        iroff2 = iroff2 + 1
      end if

    end if

    rlist(maxerr) = area1
    rlist(last) = area2
    errbnd = max ( epsabs, epsrel * abs ( area ) )
!
!  Test for roundoff error and eventually set error flag.
!
    if ( errbnd < errsum ) then

      if ( 6 <= iroff1 .or. 20 <= iroff2 ) then
        ier = 2
      end if
!
!  Set error flag in the case that the number of subintervals
!  equals limit.
!
      if ( last == limit ) then
        ier = 1
      end if
!
!  Set error flag in the case of bad integrand behavior
!  at a point of the integration range.
!
      if ( max ( abs ( a1 ), abs ( b2 ) ) <= ( 1.0E+00 + c * 1.0E+03 * &
        epsilon ( a1 ) ) * ( abs ( a2 ) + 1.0E+04 * tiny ( a2 ) ) ) then
        ier = 3
      end if

    end if
!
!  Append the newly-created intervals to the list.
!
    if ( error2 <= error1 ) then
      alist(last) = a2
      blist(maxerr) = b1
      blist(last) = b2
      elist(maxerr) = error1
      elist(last) = error2
    else
      alist(maxerr) = a2
      alist(last) = a1
      blist(last) = b1
      rlist(maxerr) = area2
      rlist(last) = area1
      elist(maxerr) = error2
      elist(last) = error1
    end if
!
!  Call QSORT to maintain the descending ordering
!  in the list of error estimates and select the subinterval
!  with the largest error estimate (to be bisected next).
!
    call qsort ( limit, last, maxerr, errmax, elist, iord, nrmax )
 
    if ( ier /= 0 .or. errsum <= errbnd ) then
      exit
    end if

  end do
!
!  Compute final result.
!
  result = sum ( rlist(1:last) )

  abserr = errsum

  if ( keyf /= 1 ) then
    neval = ( 10 * keyf + 1 ) * ( 2 * neval + 1 )
  else
    neval = 30 * neval + 15
  end if

  return
end subroutine
subroutine qagi ( f_ptr, dat, bound, inf, epsabs, epsrel, result, abserr, neval, ier )

!*****************************************************************************80
!
!! QAGI estimates an integral over a semi-infinite or infinite interval.
!
!  Discussion:
!
!    The routine calculates an approximation RESULT to a definite integral   
!      I = integral of F over (A, +Infinity), 
!    or 
!      I = integral of F over (-Infinity,A)
!    or 
!      I = integral of F over (-Infinity,+Infinity),
!    hopefully satisfying
!      || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real ( kind = 8 ) F, the name of the function routine, of the form
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!    which evaluates the integrand function.
!
!    Input, real ( kind = 8 ) BOUND, the value of the finite endpoint of the integration
!    range, if any, that is, if INF is 1 or -1.
!
!    Input, integer ( kind = 8 ) INF, indicates the type of integration range.
!    1:  (  BOUND,    +Infinity),
!    -1: ( -Infinity,  BOUND),
!    2:  ( -Infinity, +Infinity).
!
!    Input, real ( kind = 8 ) EPSABS, EPSREL, the absolute and relative accuracy requested.
!
!    Output, real ( kind = 8 ) RESULT, the estimated value of the integral.
!
!    Output, real ( kind = 8 ) ABSERR, an estimate of || I - RESULT ||.
!
!    Output, integer ( kind = 8 ) NEVAL, the number of times the integral was evaluated.
!
!    Output, integer ( kind = 8 ) IER, error indicator.
!    0, normal and reliable termination of the routine.  It is assumed that 
!      the requested accuracy has been achieved.
!    > 0,  abnormal termination of the routine.  The estimates for result
!      and error are less reliable.  It is assumed that the requested
!      accuracy has not been achieved.
!    1, maximum number of subdivisions allowed has been achieved.  One can 
!      allow more subdivisions by increasing the data value of LIMIT in QAGI
!      (and taking the according dimension adjustments into account).
!      However, if this yields no improvement it is advised to analyze the
!      integrand in order to determine the integration difficulties.  If the
!      position of a local difficulty can be determined (e.g. singularity,
!      discontinuity within the interval) one will probably gain from
!      splitting up the interval at this point and calling the integrator 
!      on the subranges.  If possible, an appropriate special-purpose 
!      integrator should be used, which is designed for handling the type
!      of difficulty involved.
!    2, the occurrence of roundoff error is detected, which prevents the
!      requested tolerance from being achieved.  The error may be
!      under-estimated.
!    3, extremely bad integrand behavior occurs at some points of the
!      integration interval.
!    4, the algorithm does not converge.  Roundoff error is detected in the
!      extrapolation table.  It is assumed that the requested tolerance
!      cannot be achieved, and that the returned result is the best which 
!      can be obtained.
!    5, the integral is probably divergent, or slowly convergent.  It must 
!      be noted that divergence can occur with any other value of IER.
!    6, the input is invalid, because INF /= 1 and INF /= -1 and INF /= 2, or
!      epsabs < 0 and epsrel < 0.  result, abserr, neval are set to zero.
!
!  Local parameters:
!
!            the dimension of rlist2 is determined by the value of
!            limexp in QEXTR.
!
!           alist     - list of left end points of all subintervals
!                       considered up to now
!           blist     - list of right end points of all subintervals
!                       considered up to now
!           rlist(i)  - approximation to the integral over
!                       (alist(i),blist(i))
!           rlist2    - array of dimension at least (limexp+2),
!                       containing the part of the epsilon table
!                       which is still needed for further computations
!           elist(i)  - error estimate applying to rlist(i)
!           maxerr    - pointer to the interval with largest error
!                       estimate
!           errmax    - elist(maxerr)
!           erlast    - error on the interval currently subdivided
!                       (before that subdivision has taken place)
!           area      - sum of the integrals over the subintervals
!           errsum    - sum of the errors over the subintervals
!           errbnd    - requested accuracy max(epsabs,epsrel*
!                       abs(result))
!           *****1    - variable for the left subinterval
!           *****2    - variable for the right subinterval
!           last      - index for subdivision
!           nres      - number of calls to the extrapolation routine
!           numrl2    - number of elements currently in rlist2. if an
!                       appropriate approximation to the compounded
!                       integral has been obtained, it is put in
!                       rlist2(numrl2) after numrl2 has been increased
!                       by one.
!           small     - length of the smallest interval considered up
!                       to now, multiplied by 1.5
!           erlarg    - sum of the errors over the intervals larger
!                       than the smallest interval considered up to now
!           extrap    - logical variable denoting that the routine
!                       is attempting to perform extrapolation. i.e.
!                       before subdividing the smallest interval we
!                       try to decrease the value of erlarg.
!           noext     - logical variable denoting that extrapolation
!                       is no longer allowed (true-value)
!
  implicit none

  integer ( kind = 8 ), parameter :: limit = 500

  real ( kind = 8 ) abseps
  real ( kind = 8 ) abserr
  real ( kind = 8 ) alist(limit)
  real ( kind = 8 ) area
  real ( kind = 8 ) area1
  real ( kind = 8 ) area12
  real ( kind = 8 ) area2
  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  real ( kind = 8 ) blist(limit)
  real ( kind = 8 ) boun
  real ( kind = 8 ) bound
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) correc
  real ( kind = 8 ) defabs
  real ( kind = 8 ) defab1
  real ( kind = 8 ) defab2
  real ( kind = 8 ) dres
  real ( kind = 8 ) elist(limit)
  real ( kind = 8 ) epsabs
  real ( kind = 8 ) epsrel
  real ( kind = 8 ) erlarg
  real ( kind = 8 ) erlast
  real ( kind = 8 ) errbnd
  real ( kind = 8 ) errmax
  real ( kind = 8 ) error1
  real ( kind = 8 ) error2
  real ( kind = 8 ) erro12
  real ( kind = 8 ) errsum
  real ( kind = 8 ) ertest
  logical extrap
  !real ( kind = 8 ), external :: f
  procedure (f_int_dat), intent(in), pointer :: f_ptr 
  integer ( kind = 8 ) id
  integer ( kind = 8 ) ier
  integer ( kind = 8 ) ierro
  integer ( kind = 8 ) inf
  integer ( kind = 8 ) iord(limit)
  integer ( kind = 8 ) iroff1
  integer ( kind = 8 ) iroff2
  integer ( kind = 8 ) iroff3
  integer ( kind = 8 ) jupbnd
  integer ( kind = 8 ) k
  integer ( kind = 8 ) ksgn
  integer ( kind = 8 ) ktmin
  integer ( kind = 8 ) last
  integer ( kind = 8 ) maxerr
  integer ( kind = 8 ) neval
  logical noext
  integer ( kind = 8 ) nres
  integer ( kind = 8 ) nrmax
  integer ( kind = 8 ) numrl2
  real ( kind = 8 ) resabs
  real ( kind = 8 ) reseps
  real ( kind = 8 ) result
  real ( kind = 8 ) res3la(3)
  real ( kind = 8 ) rlist(limit)
  real ( kind = 8 ) rlist2(52)
  real ( kind = 8 ) small
  class( dataCollectionBase), target :: dat
!
!  Test on validity of parameters.
!
  ier = 0
  neval = 0
  last = 0
  result = 0.0E+00
  abserr = 0.0E+00
  alist(1) = 0.0E+00
  blist(1) = 1.0E+00
  rlist(1) = 0.0E+00
  elist(1) = 0.0E+00
  iord(1) = 0

  if ( epsabs < 0.0E+00 .and. epsrel < 0.0E+00 ) then
    ier = 6
    return
  end if
!
!  First approximation to the integral.
!
!  Determine the interval to be mapped onto (0,1).
!  If INF = 2 the integral is computed as i = i1+i2, where
!  i1 = integral of f over (-infinity,0),
!  i2 = integral of f over (0,+infinity).
!
  if ( inf == 2 ) then
    boun = 0.0E+00
  else
    boun = bound
  end if

  call qk15i ( f_ptr, dat, boun, inf, 0.0E+00, 1.0E+00, result, abserr, defabs, resabs )
!
!  Test on accuracy.
!
  last = 1
  rlist(1) = result
  elist(1) = abserr
  iord(1) = 1
  dres = abs ( result )
  errbnd = max ( epsabs, epsrel * dres )

  if ( abserr <= 100.0E+00 * epsilon ( defabs ) * defabs .and. &
    errbnd < abserr ) then
    ier = 2
  end if

  if ( limit == 1 ) then
    ier = 1
  end if

  if ( ier /= 0 .or. (abserr <= errbnd .and. abserr /= resabs ) .or. &
    abserr == 0.0E+00 ) go to 130
!
!  Initialization.
!
  rlist2(1) = result
  errmax = abserr
  maxerr = 1
  area = result
  errsum = abserr
  abserr = huge ( abserr )
  nrmax = 1
  nres = 0
  ktmin = 0
  numrl2 = 2
  extrap = .false.
  noext = .false.
  ierro = 0
  iroff1 = 0
  iroff2 = 0
  iroff3 = 0

  if ( ( 1.0E+00 - 5.0E+01 * epsilon ( defabs ) ) * defabs <= dres ) then
    ksgn = 1
  else
    ksgn = -1
  end if

  do last = 2, limit
!
!  Bisect the subinterval with nrmax-th largest error estimate.
!
    a1 = alist(maxerr)
    b1 = 5.0E-01 * ( alist(maxerr) + blist(maxerr) )
    a2 = b1
    b2 = blist(maxerr)
    erlast = errmax
    call qk15i ( f_ptr, dat, boun, inf, a1, b1, area1, error1, resabs, defab1 )
    call qk15i ( f_ptr, dat, boun, inf, a2, b2, area2, error2, resabs, defab2 )
!
!  Improve previous approximations to integral and error
!  and test for accuracy.
!
    area12 = area1 + area2
    erro12 = error1 + error2
    errsum = errsum + erro12 - errmax
    area = area + area12 - rlist(maxerr)

    if ( defab1 /= error1 .and. defab2 /= error2 ) then

      if ( abs ( rlist(maxerr) - area12 ) <= 1.0E-05 * abs ( area12 ) &
        .and. 9.9E-01 * errmax <= erro12 ) then

        if ( extrap ) then
          iroff2 = iroff2 + 1
        end if

        if ( .not. extrap ) then
          iroff1 = iroff1 + 1
        end if

      end if

      if ( 10 < last .and. errmax < erro12 ) then
        iroff3 = iroff3 + 1
      end if

    end if

    rlist(maxerr) = area1
    rlist(last) = area2
    errbnd = max ( epsabs, epsrel * abs ( area ) )
!
!  Test for roundoff error and eventually set error flag.
!
    if ( 10 <= iroff1 + iroff2 .or. 20 <= iroff3 ) then
      ier = 2
    end if

    if ( 5 <= iroff2 ) then
      ierro = 3
    end if
!
!  Set error flag in the case that the number of subintervals equals LIMIT.
!
    if ( last == limit ) then
      ier = 1
    end if
!
!  Set error flag in the case of bad integrand behavior
!  at some points of the integration range.
!
    if ( max ( abs(a1), abs(b2) ) <= (1.0E+00 + 1.0E+03 * epsilon ( a1 ) ) * &
    ( abs(a2) + 1.0E+03 * tiny ( a2 ) )) then
      ier = 4
    end if
!
!  Append the newly-created intervals to the list.
!
    if ( error2 <= error1 ) then
      alist(last) = a2
      blist(maxerr) = b1
      blist(last) = b2
      elist(maxerr) = error1
      elist(last) = error2
    else
      alist(maxerr) = a2
      alist(last) = a1
      blist(last) = b1
      rlist(maxerr) = area2
      rlist(last) = area1
      elist(maxerr) = error2
      elist(last) = error1
    end if
!
!  Call QSORT to maintain the descending ordering
!  in the list of error estimates and select the subinterval
!  with NRMAX-th largest error estimate (to be bisected next).
!
    call qsort ( limit, last, maxerr, errmax, elist, iord, nrmax )

    if ( errsum <= errbnd ) go to 115

    if ( ier /= 0 ) then
      exit
    end if

    if ( last == 2 ) then
      small = 3.75E-01
      erlarg = errsum
      ertest = errbnd
      rlist2(2) = area
      cycle
    end if

    if ( noext ) then
      cycle
    end if

    erlarg = erlarg - erlast

    if ( small < abs ( b1 - a1 ) ) then
      erlarg = erlarg + erro12
    end if
!
!  Test whether the interval to be bisected next is the
!  smallest interval.
!
    if ( .not. extrap ) then

      if ( small < abs ( blist(maxerr) - alist(maxerr) ) ) then
        cycle
      end if

      extrap = .true.
      nrmax = 2

    end if

    if ( ierro == 3 .or. erlarg <= ertest ) then
      go to 60
    end if
!
!  The smallest interval has the largest error.
!  before bisecting decrease the sum of the errors over the
!  larger intervals (erlarg) and perform extrapolation.
!
    id = nrmax
    jupbnd = last

    if ( (2+limit/2) < last ) then
      jupbnd = limit + 3 - last
    end if

    do k = id, jupbnd
      maxerr = iord(nrmax)
      errmax = elist(maxerr)
      if ( small < abs ( blist(maxerr) - alist(maxerr) ) ) then
        go to 90
      end if
      nrmax = nrmax + 1
    end do
!
!  Extrapolate.
!
60  continue

    numrl2 = numrl2 + 1
    rlist2(numrl2) = area
    call qextr ( numrl2, rlist2, reseps, abseps, res3la, nres ) 
    ktmin = ktmin+1

    if ( 5 < ktmin .and. abserr < 1.0E-03 * errsum ) then
      ier = 5
    end if

    if ( abseps < abserr ) then

      ktmin = 0
      abserr = abseps
      result = reseps
      correc = erlarg
      ertest = max ( epsabs, epsrel * abs(reseps) )

      if ( abserr <= ertest ) then
        exit
      end if

    end if
!
!  Prepare bisection of the smallest interval.
!
    if ( numrl2 == 1 ) then
      noext = .true.
    end if

    if ( ier == 5 ) then
      exit
    end if

    maxerr = iord(1)
    errmax = elist(maxerr)
    nrmax = 1
    extrap = .false.
    small = small * 5.0E-01
    erlarg = errsum

90  continue

  end do
!
!  Set final result and error estimate.
!
  if ( abserr == huge ( abserr ) ) then
    go to 115
  end if

  if ( ( ier + ierro ) == 0 ) then
    go to 110
  end if

  if ( ierro == 3 ) then
    abserr = abserr + correc
  end if

  if ( ier == 0 ) then
    ier = 3
  end if

  if ( result /= 0.0E+00 .and. area /= 0.0E+00) then
    go to 105
  end if

  if ( errsum < abserr ) then
    go to 115
  end if

  if ( area == 0.0E+00 ) then
    go to 130
  end if

  go to 110

105 continue

  if ( errsum / abs ( area ) < abserr / abs ( result )  ) then
    go to 115
  end if
!
!  Test on divergence
!
110 continue

  if ( ksgn == (-1) .and. &
  max ( abs(result), abs(area) ) <=  defabs * 1.0E-02) go to 130

  if ( 1.0E-02 > (result/area) .or. &
    (result/area) > 1.0E+02 .or. &
    errsum > abs(area)) then
    ier = 6
  end if

  go to 130
!
!  Compute global integral sum.
!
  115 continue

  result = sum ( rlist(1:last) )

  abserr = errsum
  130 continue

  neval = 30 * last - 15
  if ( inf == 2 ) then
    neval = 2 * neval
  end if

  if ( 2 < ier ) then
    ier = ier - 1
  end if

  return
end subroutine
subroutine qagp ( f_ptr, dat, a, b, npts2, points, epsabs, epsrel, result, abserr, &
  neval, ier )

!*****************************************************************************80
!
!! QAGP computes a definite integral.
!
!  Discussion:
!
!    The routine calculates an approximation RESULT to a definite integral   
!      I = integral of F over (A,B),
!    hopefully satisfying
!      || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
!
!    Interior break points of the integration interval,
!    where local difficulties of the integrand may occur, such as
!    singularities or discontinuities, are provided by the user.
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real ( kind = 8 ) F, the name of the function routine, 
!    of the form
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!    which evaluates the integrand function.
!
!    Input, real ( kind = 8 ) A, B, the limits of integration.
!
!    Input, integer ( kind = 8 ) NPTS2, the number of user-supplied break points within 
!    the integration range, plus 2.  NPTS2 must be at least 2.
!
!    Input/output, real ( kind = 8 ) POINTS(NPTS2), contains the user provided interior
!    breakpoints in entries 1 through NPTS2-2.  If these points are not
!    in ascending order on input, they will be sorted.
!
!    Input, real ( kind = 8 ) EPSABS, EPSREL, the absolute and relative accuracy requested.
!
!    Output, real ( kind = 8 ) RESULT, the estimated value of the integral.
!
!    Output, real ( kind = 8 ) ABSERR, an estimate of || I - RESULT ||.
!
!    Output, integer ( kind = 8 ) NEVAL, the number of times the integral was evaluated.
!
!    Output, integer ( kind = 8 ) IER, return flag.
!                     ier = 0 normal and reliable termination of the
!                             routine. it is assumed that the requested
!                             accuracy has been achieved.
!                     ier > 0 abnormal termination of the routine.
!                             the estimates for integral and error are
!                             less reliable. it is assumed that the
!                             requested accuracy has not been achieved.
!                     ier = 1 maximum number of subdivisions allowed
!                             has been achieved. one can allow more
!                             subdivisions by increasing the data value
!                             of limit in qagp(and taking the according
!                             dimension adjustments into account).
!                             however, if this yields no improvement
!                             it is advised to analyze the integrand
!                             in order to determine the integration
!                             difficulties. if the position of a local
!                             difficulty can be determined (i.e.
!                             singularity, discontinuity within the
!                             interval), it should be supplied to the
!                             routine as an element of the vector
!                             points. if necessary, an appropriate
!                             special-purpose integrator must be used,
!                             which is designed for handling the type
!                             of difficulty involved.
!                         = 2 the occurrence of roundoff error is
!                             detected, which prevents the requested
!                             tolerance from being achieved.
!                             the error may be under-estimated.
!                         = 3 extremely bad integrand behavior occurs
!                             at some points of the integration
!                             interval.
!                         = 4 the algorithm does not converge. roundoff
!                             error is detected in the extrapolation
!                             table. it is presumed that the requested
!                             tolerance cannot be achieved, and that
!                             the returned result is the best which
!                             can be obtained.
!                         = 5 the integral is probably divergent, or
!                             slowly convergent. it must be noted that
!                             divergence can occur with any other value
!                             of ier > 0.
!                         = 6 the input is invalid because
!                             npts2 < 2 or
!                             break points are specified outside
!                             the integration range or
!                             epsabs < 0 and epsrel < 0,
!                             or limit < npts2.
!                             result, abserr, neval are set to zero.
!
!  Local parameters:
!
!            the dimension of rlist2 is determined by the value of
!            limexp in QEXTR (rlist2 should be of dimension
!            (limexp+2) at least).
!
!           alist     - list of left end points of all subintervals
!                       considered up to now
!           blist     - list of right end points of all subintervals
!                       considered up to now
!           rlist(i)  - approximation to the integral over
!                       (alist(i),blist(i))
!           rlist2    - array of dimension at least limexp+2
!                       containing the part of the epsilon table which
!                       is still needed for further computations
!           elist(i)  - error estimate applying to rlist(i)
!           maxerr    - pointer to the interval with largest error
!                       estimate
!           errmax    - elist(maxerr)
!           erlast    - error on the interval currently subdivided
!                       (before that subdivision has taken place)
!           area      - sum of the integrals over the subintervals
!           errsum    - sum of the errors over the subintervals
!           errbnd    - requested accuracy max(epsabs,epsrel*
!                       abs(result))
!           *****1    - variable for the left subinterval
!           *****2    - variable for the right subinterval
!           last      - index for subdivision
!           nres      - number of calls to the extrapolation routine
!           numrl2    - number of elements in rlist2. if an appropriate
!                       approximation to the compounded integral has
!                       obtained, it is put in rlist2(numrl2) after
!                       numrl2 has been increased by one.
!           erlarg    - sum of the errors over the intervals larger
!                       than the smallest interval considered up to now
!           extrap    - logical variable denoting that the routine
!                       is attempting to perform extrapolation. i.e.
!                       before subdividing the smallest interval we
!                       try to decrease the value of erlarg.
!           noext     - logical variable denoting that extrapolation is
!                       no longer allowed (true-value)
!
  implicit none

  integer ( kind = 8 ), parameter :: limit = 500

  real ( kind = 8 ) a
  real ( kind = 8 ) abseps
  real ( kind = 8 ) abserr
  real ( kind = 8 ) alist(limit)
  real ( kind = 8 ) area
  real ( kind = 8 ) area1
  real ( kind = 8 ) area12
  real ( kind = 8 ) area2
  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  real ( kind = 8 ) b
  real ( kind = 8 ) blist(limit)
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) correc
  real ( kind = 8 ) defabs
  real ( kind = 8 ) defab1
  real ( kind = 8 ) defab2
  real ( kind = 8 ) dres
  real ( kind = 8 ) elist(limit)
  real ( kind = 8 ) epsabs
  real ( kind = 8 ) epsrel
  real ( kind = 8 ) erlarg
  real ( kind = 8 ) erlast
  real ( kind = 8 ) errbnd
  real ( kind = 8 ) errmax
  real ( kind = 8 ) error1
  real ( kind = 8 ) erro12
  real ( kind = 8 ) error2
  real ( kind = 8 ) errsum
  real ( kind = 8 ) ertest
  logical extrap
  !real ( kind = 8 ), external :: f
  procedure (f_int_dat), intent(in), pointer :: f_ptr 
  integer ( kind = 8 ) i
  integer ( kind = 8 ) id
  integer ( kind = 8 ) ier
  integer ( kind = 8 ) ierro
  integer ( kind = 8 ) ind1
  integer ( kind = 8 ) ind2
  integer ( kind = 8 ) iord(limit)
  integer ( kind = 8 ) iroff1
  integer ( kind = 8 ) iroff2
  integer ( kind = 8 ) iroff3
  integer ( kind = 8 ) j
  integer ( kind = 8 ) jlow
  integer ( kind = 8 ) jupbnd
  integer ( kind = 8 ) k
  integer ( kind = 8 ) ksgn
  integer ( kind = 8 ) ktmin
  integer ( kind = 8 ) last
  integer ( kind = 8 ) levcur
  integer ( kind = 8 ) level(limit)
  integer ( kind = 8 ) levmax
  integer ( kind = 8 ) maxerr
  integer ( kind = 8 ) ndin(40)
  integer ( kind = 8 ) neval
  integer ( kind = 8 ) nint
  logical noext
  integer ( kind = 8 ) npts
  integer ( kind = 8 ) npts2
  integer ( kind = 8 ) nres
  integer ( kind = 8 ) nrmax
  integer ( kind = 8 ) numrl2
  real ( kind = 8 ) points(40)
  real ( kind = 8 ) pts(40)
  real ( kind = 8 ) resa
  real ( kind = 8 ) resabs
  real ( kind = 8 ) reseps
  real ( kind = 8 ) result
  real ( kind = 8 ) res3la(3)
  real ( kind = 8 ) rlist(limit)
  real ( kind = 8 ) rlist2(52)
  real ( kind = 8 ) sign
  real ( kind = 8 ) temp
  class(dataCollectionBase), target :: dat
!
!  Test on validity of parameters.
!
  ier = 0
  neval = 0
  last = 0
  result = 0.0E+00
  abserr = 0.0E+00
  alist(1) = a
  blist(1) = b
  rlist(1) = 0.0E+00
  elist(1) = 0.0E+00
  iord(1) = 0
  level(1) = 0
  npts = npts2 - 2

  if ( npts2 < 2 ) then
    ier = 6
    return
  else if ( limit <= npts .or. ( epsabs < 0.0E+00 .and. &
    epsrel < 0.0E+00) ) then
    ier = 6
    return
  end if
!
!  If any break points are provided, sort them into an
!  ascending sequence.
!
  if ( b < a ) then
    sign = -1.0E+00
  else
    sign = +1.0E+00
  end if

  pts(1) = min ( a, b )

  do i = 1, npts
    pts(i+1) = points(i)
  end do

  pts(npts+2) = max ( a, b )
  nint = npts+1
  a1 = pts(1)

  if ( npts /= 0 ) then

    do i = 1, nint
      do j = i+1, nint+1
        if ( pts(j) < pts(i) ) then
          temp   = pts(i)
          pts(i) = pts(j)
          pts(j) = temp
        end if
      end do
    end do

    if ( pts(1) /= min ( a, b ) .or. pts(nint+1) /= max ( a, b ) ) then
      ier = 6
      return
    end if

  end if
!
!  Compute first integral and error approximations.
!
  resabs = 0.0E+00

  do i = 1, nint

    b1 = pts(i+1)
    call qk21_x ( f_ptr, dat, a1, b1, area1, error1, defabs, resa )
    abserr = abserr + error1
    result = result + area1
    ndin(i) = 0

    if ( error1 == resa .and. error1 /= 0.0E+00 ) then
      ndin(i) = 1
    end if

    resabs = resabs + defabs
    level(i) = 0
    elist(i) = error1
    alist(i) = a1
    blist(i) = b1
    rlist(i) = area1
    iord(i) = i
    a1 = b1

  end do

  errsum = 0.0E+00

  do i = 1, nint
    if ( ndin(i) == 1 ) then
      elist(i) = abserr
    end if
    errsum = errsum + elist(i)
  end do
!
!  Test on accuracy.
!
  last = nint
  neval = 21 * nint
  dres = abs ( result )
  errbnd = max ( epsabs, epsrel * dres )

  if ( abserr <= 1.0E+02 * epsilon ( resabs ) * resabs .and. &
    abserr > errbnd ) then
    ier = 2
  end if

  if ( nint /= 1 ) then

    do i = 1, npts

      jlow = i+1
      ind1 = iord(i)

      do j = jlow, nint
        ind2 = iord(j)
        if ( elist(ind1) <= elist(ind2) ) then
          ind1 = ind2
          k = j
        end if
      end do

      if ( ind1 /= iord(i) ) then
        iord(k) = iord(i)
        iord(i) = ind1
      end if

    end do

    if ( limit < npts2 ) then
      ier = 1
    end if

  end if

  if ( ier /= 0 .or. abserr <= errbnd ) then
    return
  end if
!
!  Initialization
!
  rlist2(1) = result
  maxerr = iord(1)
  errmax = elist(maxerr)
  area = result
  nrmax = 1
  nres = 0
  numrl2 = 1
  ktmin = 0
  extrap = .false.
  noext = .false.
  erlarg = errsum
  ertest = errbnd
  levmax = 1
  iroff1 = 0
  iroff2 = 0
  iroff3 = 0
  ierro = 0
  abserr = huge ( abserr )

  if ( dres >= ( 1.0E+00 - 0.5E+00 * epsilon ( resabs ) ) * resabs ) then
    ksgn = 1
  else
    ksgn = -1
  end if

  do last = npts2, limit
!
!  Bisect the subinterval with the NRMAX-th largest error estimate.
!
    levcur = level(maxerr) + 1
    a1 = alist(maxerr)
    b1 = 0.5E+00 * ( alist(maxerr) + blist(maxerr) )
    a2 = b1
    b2 = blist(maxerr)
    erlast = errmax
    call qk21_x ( f_ptr, dat, a1, b1, area1, error1, resa, defab1 )
    call qk21_x ( f_ptr, dat, a2, b2, area2, error2, resa, defab2 )
!
!  Improve previous approximations to integral and error
!  and test for accuracy.
!
    neval = neval + 42
    area12 = area1 + area2
    erro12 = error1 + error2
    errsum = errsum + erro12 -errmax
    area = area + area12 - rlist(maxerr)

    if ( defab1 /= error1 .and. defab2 /= error2 ) then

      if ( abs ( rlist ( maxerr ) - area12 ) <= 1.0E-05 * abs(area12) .and. &
        erro12 >= 9.9E-01 * errmax ) then

        if ( extrap ) then
          iroff2 = iroff2+1
        else
          iroff1 = iroff1+1
        end if

      end if

      if ( last > 10 .and. erro12 > errmax ) then
        iroff3 = iroff3 + 1
      end if

    end if

    level(maxerr) = levcur
    level(last) = levcur
    rlist(maxerr) = area1
    rlist(last) = area2
    errbnd = max ( epsabs, epsrel * abs ( area ) )
!
!  Test for roundoff error and eventually set error flag.
!
    if ( 10 <= iroff1 + iroff2 .or. 20 <= iroff3 ) then
      ier = 2
    end if

    if ( 5 <= iroff2 ) then
      ierro = 3
    end if
!
!  Set error flag in the case that the number of subintervals
!  equals limit.
!
    if ( last == limit ) then
      ier = 1
    end if
!
!  Set error flag in the case of bad integrand behavior
!  at a point of the integration range
!
    if ( max ( abs(a1), abs(b2)) <= ( 1.0E+00 + 1.0E+03 * epsilon ( a1 ) )* &
    ( abs ( a2 ) + 1.0E+03 * tiny ( a2 ) ) ) then
      ier = 4
    end if
!
!  Append the newly-created intervals to the list.
!
    if ( error2 <= error1 ) then
      alist(last) = a2
      blist(maxerr) = b1
      blist(last) = b2
      elist(maxerr) = error1
      elist(last) = error2
    else
      alist(maxerr) = a2
      alist(last) = a1
      blist(last) = b1
      rlist(maxerr) = area2
      rlist(last) = area1
      elist(maxerr) = error2
      elist(last) = error1
    end if
!
!  Call QSORT to maintain the descending ordering
!  in the list of error estimates and select the subinterval
!  with nrmax-th largest error estimate (to be bisected next).
!
    call qsort ( limit, last, maxerr, errmax, elist, iord, nrmax )

    if ( errsum <= errbnd ) then
      go to 190
    end if

    if ( ier /= 0 ) then
      exit
    end if

    if ( noext ) then
      cycle
    end if

    erlarg = erlarg - erlast

    if ( levcur+1 <= levmax ) then
      erlarg = erlarg + erro12
    end if
!
!  Test whether the interval to be bisected next is the
!  smallest interval.
!
    if ( .not. extrap ) then

      if ( level(maxerr)+1 <= levmax ) then
        cycle
      end if

      extrap = .true.
      nrmax = 2

    end if
!
!  The smallest interval has the largest error.
!  Before bisecting decrease the sum of the errors over the
!  larger intervals (erlarg) and perform extrapolation.
!
    if ( ierro /= 3 .and. erlarg > ertest ) then

      id = nrmax
      jupbnd = last
      if ( last > (2+limit/2) ) then
        jupbnd = limit+3-last
      end if

      do k = id, jupbnd
        maxerr = iord(nrmax)
        errmax = elist(maxerr)
        if ( level(maxerr)+1 <= levmax ) then
          go to 160
        end if
        nrmax = nrmax + 1
      end do

    end if
!
!  Perform extrapolation.
!
    numrl2 = numrl2 + 1
    rlist2(numrl2) = area

    if ( numrl2 <= 2 ) then
      go to 155
    end if

    call qextr ( numrl2, rlist2, reseps, abseps, res3la, nres )
    ktmin = ktmin+1

    if ( 5 < ktmin .and. abserr < 1.0E-03 * errsum ) then
      ier = 5
    end if

    if ( abseps < abserr ) then

      ktmin = 0
      abserr = abseps
      result = reseps
      correc = erlarg
      ertest = max ( epsabs, epsrel * abs(reseps) )

      if ( abserr < ertest ) then
        exit
      end if

    end if
!
!  Prepare bisection of the smallest interval.
!
    if ( numrl2 == 1 ) then
      noext = .true.
    end if

    if ( 5 <= ier ) then
      exit
    end if

155 continue

    maxerr = iord(1)
    errmax = elist(maxerr)
    nrmax = 1
    extrap = .false.
    levmax = levmax + 1
    erlarg = errsum

160 continue

  end do
!
!  Set the final result.
!
  if ( abserr == huge ( abserr ) ) then
    go to 190
  end if

  if ( ( ier + ierro ) == 0 ) then
    go to 180
  end if

  if ( ierro == 3 ) then
    abserr = abserr + correc
  end if

  if ( ier == 0 ) then
    ier = 3
  end if

  if ( result /= 0.0E+00 .and. area /= 0.0E+00 ) then
    go to 175
  end if

  if ( errsum < abserr ) then
    go to 190
  end if

  if ( area == 0.0E+00 ) then
    go to 210
  end if

  go to 180

175 continue

  if ( abserr / abs(result) > errsum / abs(area) ) then
    go to 190
  end if
!
!  Test on divergence.
!
  180 continue

  if ( ksgn == (-1) .and. max ( abs(result),abs(area)) <=  &
    resabs*1.0E-02 ) go to 210

  if ( 1.0E-02 > (result/area) .or. (result/area) > 1.0E+02 .or. &
    errsum > abs(area) ) then
    ier = 6
  end if

  go to 210
!
!  Compute global integral sum.
!
190 continue

  result = sum ( rlist(1:last) )

  abserr = errsum

210 continue

  if ( 2 < ier ) then
    ier = ier - 1
  end if

  result = result * sign

  return
end subroutine
subroutine qags ( f_ptr, dat, a, b, epsabs, epsrel, result, abserr, neval, ier )

!*****************************************************************************80
!
!! QAGS estimates the integral of a function.
!
!  Discussion:
!
!    The routine calculates an approximation RESULT to a definite integral   
!      I = integral of F over (A,B),
!    hopefully satisfying
!      || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real ( kind = 8 ) F, the name of the function routine, of the form
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!    which evaluates the integrand function.
!
!    Input, real ( kind = 8 ) A, B, the limits of integration.
!
!    Input, real ( kind = 8 ) EPSABS, EPSREL, the absolute and relative accuracy requested.
!
!    Output, real ( kind = 8 ) RESULT, the estimated value of the integral.
!
!    Output, real ( kind = 8 ) ABSERR, an estimate of || I - RESULT ||.
!
!    Output, integer ( kind = 8 ) NEVAL, the number of times the integral was evaluated.
!
!    Output, integer ( kind = 8 ) IER, error flag.
!                     ier = 0 normal and reliable termination of the
!                             routine. it is assumed that the requested
!                             accuracy has been achieved.
!                     ier > 0 abnormal termination of the routine
!                             the estimates for integral and error are
!                             less reliable. it is assumed that the
!                             requested accuracy has not been achieved.
!                         = 1 maximum number of subdivisions allowed
!                             has been achieved. one can allow more sub-
!                             divisions by increasing the data value of
!                             limit in qags (and taking the according
!                             dimension adjustments into account).
!                             however, if this yields no improvement
!                             it is advised to analyze the integrand
!                             in order to determine the integration
!                             difficulties. if the position of a
!                             local difficulty can be determined (e.g.
!                             singularity, discontinuity within the
!                             interval) one will probably gain from
!                             splitting up the interval at this point
!                             and calling the integrator on the sub-
!                             ranges. if possible, an appropriate
!                             special-purpose integrator should be used,
!                             which is designed for handling the type
!                             of difficulty involved.
!                         = 2 the occurrence of roundoff error is detec-
!                             ted, which prevents the requested
!                             tolerance from being achieved.
!                             the error may be under-estimated.
!                         = 3 extremely bad integrand behavior occurs
!                             at some  points of the integration
!                             interval.
!                         = 4 the algorithm does not converge. roundoff
!                             error is detected in the extrapolation
!                             table. it is presumed that the requested
!                             tolerance cannot be achieved, and that the
!                             returned result is the best which can be
!                             obtained.
!                         = 5 the integral is probably divergent, or
!                             slowly convergent. it must be noted that
!                             divergence can occur with any other value
!                             of ier.
!                         = 6 the input is invalid, because
!                             epsabs < 0 and epsrel < 0,
!                             result, abserr and neval are set to zero.
!
!  Local Parameters:
!
!           alist     - list of left end points of all subintervals
!                       considered up to now
!           blist     - list of right end points of all subintervals
!                       considered up to now
!           rlist(i)  - approximation to the integral over
!                       (alist(i),blist(i))
!           rlist2    - array of dimension at least limexp+2 containing
!                       the part of the epsilon table which is still
!                       needed for further computations
!           elist(i)  - error estimate applying to rlist(i)
!           maxerr    - pointer to the interval with largest error
!                       estimate
!           errmax    - elist(maxerr)
!           erlast    - error on the interval currently subdivided
!                       (before that subdivision has taken place)
!           area      - sum of the integrals over the subintervals
!           errsum    - sum of the errors over the subintervals
!           errbnd    - requested accuracy max(epsabs,epsrel*
!                       abs(result))
!           *****1    - variable for the left interval
!           *****2    - variable for the right interval
!           last      - index for subdivision
!           nres      - number of calls to the extrapolation routine
!           numrl2    - number of elements currently in rlist2. if an
!                       appropriate approximation to the compounded
!                       integral has been obtained it is put in
!                       rlist2(numrl2) after numrl2 has been increased
!                       by one.
!           small     - length of the smallest interval considered
!                       up to now, multiplied by 1.5
!           erlarg    - sum of the errors over the intervals larger
!                       than the smallest interval considered up to now
!           extrap    - logical variable denoting that the routine is
!                       attempting to perform extrapolation i.e. before
!                       subdividing the smallest interval we try to
!                       decrease the value of erlarg.
!           noext     - logical variable denoting that extrapolation
!                       is no longer allowed (true value)
!
  implicit none

  integer ( kind = 8 ), parameter :: limit = 500

  real ( kind = 8 ) a
  real ( kind = 8 ) abseps
  real ( kind = 8 ) abserr
  real ( kind = 8 ) alist(limit)
  real ( kind = 8 ) area
  real ( kind = 8 ) area1
  real ( kind = 8 ) area12
  real ( kind = 8 ) area2
  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  real ( kind = 8 ) b
  real ( kind = 8 ) blist(limit)
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) correc
  real ( kind = 8 ) defabs
  real ( kind = 8 ) defab1
  real ( kind = 8 ) defab2
  real ( kind = 8 ) dres
  real ( kind = 8 ) elist(limit)
  real ( kind = 8 ) epsabs
  real ( kind = 8 ) epsrel
  real ( kind = 8 ) erlarg
  real ( kind = 8 ) erlast
  real ( kind = 8 ) errbnd
  real ( kind = 8 ) errmax
  real ( kind = 8 ) error1
  real ( kind = 8 ) error2
  real ( kind = 8 ) erro12
  real ( kind = 8 ) errsum
  real ( kind = 8 ) ertest
  logical extrap
  !real ( kind = 8 ), external :: f
  procedure (f_int_dat), intent(in), pointer :: f_ptr 
  integer ( kind = 8 ) id
  integer ( kind = 8 ) ier
  integer ( kind = 8 ) ierro
  integer ( kind = 8 ) iord(limit)
  integer ( kind = 8 ) iroff1
  integer ( kind = 8 ) iroff2
  integer ( kind = 8 ) iroff3
  integer ( kind = 8 ) jupbnd
  integer ( kind = 8 ) k
  integer ( kind = 8 ) ksgn
  integer ( kind = 8 ) ktmin
  integer ( kind = 8 ) last
  logical noext
  integer ( kind = 8 ) maxerr
  integer ( kind = 8 ) neval
  integer ( kind = 8 ) nres
  integer ( kind = 8 ) nrmax
  integer ( kind = 8 ) numrl2
  real ( kind = 8 ) resabs
  real ( kind = 8 ) reseps
  real ( kind = 8 ) result
  real ( kind = 8 ) res3la(3)
  real ( kind = 8 ) rlist(limit)
  real ( kind = 8 ) rlist2(52)
  real ( kind = 8 ) small
  class(dataCollectionBase), target :: dat
!
!  The dimension of rlist2 is determined by the value of
!  limexp in QEXTR (rlist2 should be of dimension
!  (limexp+2) at least).
!
!  Test on validity of parameters.
!
  ier = 0
  neval = 0
  last = 0
  result = 0.0E+00
  abserr = 0.0E+00
  alist(1) = a
  blist(1) = b
  rlist(1) = 0.0E+00
  elist(1) = 0.0E+00

  if ( epsabs < 0.0E+00 .and. epsrel < 0.0E+00 ) then
    ier = 6
    return
  end if
!
!  First approximation to the integral.
!
  ierro = 0
  call qk21_x ( f_ptr, dat, a, b, result, abserr, defabs, resabs )
!
!  Test on accuracy.
!
  dres = abs ( result )
  errbnd = max ( epsabs, epsrel * dres )
  last = 1
  rlist(1) = result
  elist(1) = abserr
  iord(1) = 1

  if ( abserr <= 1.0E+02 * epsilon ( defabs ) * defabs .and. &
    abserr > errbnd ) then
    ier = 2
  end if

  if ( limit == 1 ) then
    ier = 1
  end if

  if ( ier /= 0 .or. (abserr <= errbnd .and. abserr /= resabs ) .or. &
    abserr == 0.0E+00 ) go to 140
!
!  Initialization.
!
  rlist2(1) = result
  errmax = abserr
  maxerr = 1
  area = result
  errsum = abserr
  abserr = huge ( abserr )
  nrmax = 1
  nres = 0
  numrl2 = 2
  ktmin = 0
  extrap = .false.
  noext = .false.
  iroff1 = 0
  iroff2 = 0
  iroff3 = 0

  if ( dres >= (1.0E+00 - 5.0E+01* epsilon ( defabs ) ) * defabs ) then
    ksgn = 1
  else
    ksgn = -1
  end if

  do last = 2, limit
!
!  Bisect the subinterval with the nrmax-th largest error estimate.
!
    a1 = alist(maxerr)
    b1 = 5.0E-01 * ( alist(maxerr) + blist(maxerr) )
    a2 = b1
    b2 = blist(maxerr)
    erlast = errmax
    call qk21_x ( f_ptr, dat, a1, b1, area1, error1, resabs, defab1 )
    call qk21_x ( f_ptr, dat, a2, b2, area2, error2, resabs, defab2 )
!
!  Improve previous approximations to integral and error
!  and test for accuracy.
!
    area12 = area1+area2
    erro12 = error1+error2
    errsum = errsum+erro12-errmax
    area = area+area12-rlist(maxerr)

    if ( defab1 == error1 .or. defab2 == error2 ) go to 15

    if ( abs ( rlist(maxerr) - area12) > 1.0E-05 * abs(area12) &
      .or. erro12 < 9.9E-01 * errmax ) go to 10

    if ( extrap ) then
      iroff2 = iroff2+1
    else
      iroff1 = iroff1+1
    end if

10  continue

    if ( last > 10 .and. erro12 > errmax ) then
      iroff3 = iroff3+1
    end if

15  continue

    rlist(maxerr) = area1
    rlist(last) = area2
    errbnd = max ( epsabs, epsrel*abs(area) )
!
!  Test for roundoff error and eventually set error flag.
!
    if ( iroff1+iroff2 >= 10 .or. iroff3 >= 20 ) then
      ier = 2
    end if

    if ( iroff2 >= 5 ) then
      ierro = 3
    end if
!
!  Set error flag in the case that the number of subintervals
!  equals limit.
!
    if ( last == limit ) then
      ier = 1
    end if
!
!  Set error flag in the case of bad integrand behavior
!  at a point of the integration range.
!
    if ( max ( abs(a1),abs(b2)) <= (1.0E+00+1.0E+03* epsilon ( a1 ) )* &
      (abs(a2)+1.0E+03* tiny ( a2 ) ) ) then
      ier = 4
    end if
!
!  Append the newly-created intervals to the list.
!
    if ( error2 <= error1 ) then
      alist(last) = a2
      blist(maxerr) = b1
      blist(last) = b2
      elist(maxerr) = error1
      elist(last) = error2
    else
      alist(maxerr) = a2
      alist(last) = a1
      blist(last) = b1
      rlist(maxerr) = area2
      rlist(last) = area1
      elist(maxerr) = error2
      elist(last) = error1
    end if
!
!  Call QSORT to maintain the descending ordering
!  in the list of error estimates and select the subinterval
!  with nrmax-th largest error estimate (to be bisected next).
!
    call qsort ( limit, last, maxerr, errmax, elist, iord, nrmax )

    if ( errsum <= errbnd ) go to 115

    if ( ier /= 0 ) then
      exit
    end if

    if ( last == 2 ) go to 80
    if ( noext ) go to 90

    erlarg = erlarg-erlast

    if ( abs(b1-a1) > small ) then
      erlarg = erlarg+erro12
    end if
!
!  Test whether the interval to be bisected next is the
!  smallest interval.
!
    if ( .not. extrap ) then
      if ( abs(blist(maxerr)-alist(maxerr)) > small ) go to 90
      extrap = .true.
      nrmax = 2
    end if

!40  continue
!
!  The smallest interval has the largest error.
!  Before bisecting decrease the sum of the errors over the
!  larger intervals (erlarg) and perform extrapolation.
!
    if ( ierro /= 3 .and. erlarg > ertest ) then

      id = nrmax
      jupbnd = last

      if ( last > (2+limit/2) ) then
        jupbnd = limit+3-last
      end if

      do k = id, jupbnd
        maxerr = iord(nrmax)
        errmax = elist(maxerr)
        if ( abs(blist(maxerr)-alist(maxerr)) > small ) then
          go to 90
        end if
        nrmax = nrmax+1
      end do

    end if
!
!  Perform extrapolation.
!
!60  continue

    numrl2 = numrl2+1
    rlist2(numrl2) = area
    call qextr ( numrl2, rlist2, reseps, abseps, res3la, nres )
    ktmin = ktmin+1

    if ( ktmin > 5 .and. abserr < 1.0E-03 * errsum ) then
      ier = 5
    end if

    if ( abseps < abserr ) then

      ktmin = 0
      abserr = abseps
      result = reseps
      correc = erlarg
      ertest = max ( epsabs,epsrel*abs(reseps))

      if ( abserr <= ertest ) then
        exit
      end if

    end if
!
!  Prepare bisection of the smallest interval.
!
    if ( numrl2 == 1 ) then
      noext = .true.
    end if

    if ( ier == 5 ) then
      exit
    end if

    maxerr = iord(1)
    errmax = elist(maxerr)
    nrmax = 1
    extrap = .false.
    small = small * 5.0E-01
    erlarg = errsum
    go to 90

80  continue

    small = abs ( b - a ) * 3.75E-01
    erlarg = errsum
    ertest = errbnd
    rlist2(2) = area

90  continue

  end do
!
!  Set final result and error estimate.
!
  if ( abserr == huge ( abserr ) ) then
    go to 115
  end if

  if ( ier + ierro == 0 ) then
    go to 110
  end if

  if ( ierro == 3 ) then
    abserr = abserr + correc
  end if

  if ( ier == 0 ) then
    ier = 3
  end if

  if ( result /= 0.0E+00 .and. area /= 0.0E+00 ) then
    go to 105
  end if

  if ( abserr > errsum ) go to 115
  if ( area == 0.0E+00 ) go to 130
  go to 110

105 continue

  if ( abserr/abs(result) > errsum/abs(area) ) go to 115
!
!  Test on divergence.
!
110 continue

  if ( ksgn == (-1) .and. max ( abs(result),abs(area)) <=  &
   defabs*1.0E-02 ) go to 130

  if ( 1.0E-02 > (result/area) .or. (result/area) > 1.0E+02 &
   .or. errsum > abs(area) ) then
    ier = 6
  end if

  go to 130
!
!  Compute global integral sum.
!
115 continue

  result = sum ( rlist(1:last) )

  abserr = errsum

130 continue
 
  if ( 2 < ier ) then
    ier = ier - 1
  end if

140 continue

  neval = 42*last-21

  return
end subroutine qags

!:: needed for doing double or triple integrals
subroutine qags_x ( f_ptr, dat, a, b, epsabs, epsrel, result, abserr, neval, ier )

!*****************************************************************************80
!
!! QAGS estimates the integral of a function.
!
!  Discussion:
!
!    The routine calculates an approximation RESULT to a definite integral   
!      I = integral of F over (A,B),
!    hopefully satisfying
!      || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real ( kind = 8 ) F, the name of the function routine, of the form
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!    which evaluates the integrand function.
!
!    Input, real ( kind = 8 ) A, B, the limits of integration.
!
!    Input, real ( kind = 8 ) EPSABS, EPSREL, the absolute and relative accuracy requested.
!
!    Output, real ( kind = 8 ) RESULT, the estimated value of the integral.
!
!    Output, real ( kind = 8 ) ABSERR, an estimate of || I - RESULT ||.
!
!    Output, integer ( kind = 8 ) NEVAL, the number of times the integral was evaluated.
!
!    Output, integer ( kind = 8 ) IER, error flag.
!                     ier = 0 normal and reliable termination of the
!                             routine. it is assumed that the requested
!                             accuracy has been achieved.
!                     ier > 0 abnormal termination of the routine
!                             the estimates for integral and error are
!                             less reliable. it is assumed that the
!                             requested accuracy has not been achieved.
!                         = 1 maximum number of subdivisions allowed
!                             has been achieved. one can allow more sub-
!                             divisions by increasing the data value of
!                             limit in qags (and taking the according
!                             dimension adjustments into account).
!                             however, if this yields no improvement
!                             it is advised to analyze the integrand
!                             in order to determine the integration
!                             difficulties. if the position of a
!                             local difficulty can be determined (e.g.
!                             singularity, discontinuity within the
!                             interval) one will probably gain from
!                             splitting up the interval at this point
!                             and calling the integrator on the sub-
!                             ranges. if possible, an appropriate
!                             special-purpose integrator should be used,
!                             which is designed for handling the type
!                             of difficulty involved.
!                         = 2 the occurrence of roundoff error is detec-
!                             ted, which prevents the requested
!                             tolerance from being achieved.
!                             the error may be under-estimated.
!                         = 3 extremely bad integrand behavior occurs
!                             at some  points of the integration
!                             interval.
!                         = 4 the algorithm does not converge. roundoff
!                             error is detected in the extrapolation
!                             table. it is presumed that the requested
!                             tolerance cannot be achieved, and that the
!                             returned result is the best which can be
!                             obtained.
!                         = 5 the integral is probably divergent, or
!                             slowly convergent. it must be noted that
!                             divergence can occur with any other value
!                             of ier.
!                         = 6 the input is invalid, because
!                             epsabs < 0 and epsrel < 0,
!                             result, abserr and neval are set to zero.
!
!  Local Parameters:
!
!           alist     - list of left end points of all subintervals
!                       considered up to now
!           blist     - list of right end points of all subintervals
!                       considered up to now
!           rlist(i)  - approximation to the integral over
!                       (alist(i),blist(i))
!           rlist2    - array of dimension at least limexp+2 containing
!                       the part of the epsilon table which is still
!                       needed for further computations
!           elist(i)  - error estimate applying to rlist(i)
!           maxerr    - pointer to the interval with largest error
!                       estimate
!           errmax    - elist(maxerr)
!           erlast    - error on the interval currently subdivided
!                       (before that subdivision has taken place)
!           area      - sum of the integrals over the subintervals
!           errsum    - sum of the errors over the subintervals
!           errbnd    - requested accuracy max(epsabs,epsrel*
!                       abs(result))
!           *****1    - variable for the left interval
!           *****2    - variable for the right interval
!           last      - index for subdivision
!           nres      - number of calls to the extrapolation routine
!           numrl2    - number of elements currently in rlist2. if an
!                       appropriate approximation to the compounded
!                       integral has been obtained it is put in
!                       rlist2(numrl2) after numrl2 has been increased
!                       by one.
!           small     - length of the smallest interval considered
!                       up to now, multiplied by 1.5
!           erlarg    - sum of the errors over the intervals larger
!                       than the smallest interval considered up to now
!           extrap    - logical variable denoting that the routine is
!                       attempting to perform extrapolation i.e. before
!                       subdividing the smallest interval we try to
!                       decrease the value of erlarg.
!           noext     - logical variable denoting that extrapolation
!                       is no longer allowed (true value)
!
  implicit none

  integer ( kind = 8 ), parameter :: limit = 500

  real ( kind = 8 ) a
  real ( kind = 8 ) abseps
  real ( kind = 8 ) abserr
  real ( kind = 8 ) alist(limit)
  real ( kind = 8 ) area
  real ( kind = 8 ) area1
  real ( kind = 8 ) area12
  real ( kind = 8 ) area2
  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  real ( kind = 8 ) b
  real ( kind = 8 ) blist(limit)
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) correc
  real ( kind = 8 ) defabs
  real ( kind = 8 ) defab1
  real ( kind = 8 ) defab2
  real ( kind = 8 ) dres
  real ( kind = 8 ) elist(limit)
  real ( kind = 8 ) epsabs
  real ( kind = 8 ) epsrel
  real ( kind = 8 ) erlarg
  real ( kind = 8 ) erlast
  real ( kind = 8 ) errbnd
  real ( kind = 8 ) errmax
  real ( kind = 8 ) error1
  real ( kind = 8 ) error2
  real ( kind = 8 ) erro12
  real ( kind = 8 ) errsum
  real ( kind = 8 ) ertest
  logical extrap
  !real ( kind = 8 ), external :: f
  procedure (f_int_dat), intent(in), pointer :: f_ptr 
  integer ( kind = 8 ) id
  integer ( kind = 8 ) ier
  integer ( kind = 8 ) ierro
  integer ( kind = 8 ) iord(limit)
  integer ( kind = 8 ) iroff1
  integer ( kind = 8 ) iroff2
  integer ( kind = 8 ) iroff3
  integer ( kind = 8 ) jupbnd
  integer ( kind = 8 ) k
  integer ( kind = 8 ) ksgn
  integer ( kind = 8 ) ktmin
  integer ( kind = 8 ) last
  logical noext
  integer ( kind = 8 ) maxerr
  integer ( kind = 8 ) neval
  integer ( kind = 8 ) nres
  integer ( kind = 8 ) nrmax
  integer ( kind = 8 ) numrl2
  real ( kind = 8 ) resabs
  real ( kind = 8 ) reseps
  real ( kind = 8 ) result
  real ( kind = 8 ) res3la(3)
  real ( kind = 8 ) rlist(limit)
  real ( kind = 8 ) rlist2(52)
  real ( kind = 8 ) small
  
  class(dataCollectionBase), target :: dat
  
  !f_ptr => f
  
!
!  The dimension of rlist2 is determined by the value of
!  limexp in QEXTR (rlist2 should be of dimension
!  (limexp+2) at least).
!
!  Test on validity of parameters.
!
  ier = 0
  neval = 0
  last = 0
  result = 0.0E+00
  abserr = 0.0E+00
  alist(1) = a
  blist(1) = b
  rlist(1) = 0.0E+00
  elist(1) = 0.0E+00

  if ( epsabs < 0.0E+00 .and. epsrel < 0.0E+00 ) then
    ier = 6
    return
  end if
!
!  First approximation to the integral.
!
  ierro = 0
  call qk21_x ( f_ptr, dat, a, b, result, abserr, defabs, resabs )
!
!  Test on accuracy.
!
  dres = abs ( result )
  errbnd = max ( epsabs, epsrel * dres )
  last = 1
  rlist(1) = result
  elist(1) = abserr
  iord(1) = 1

  if ( abserr <= 1.0E+02 * epsilon ( defabs ) * defabs .and. &
    abserr > errbnd ) then
    ier = 2
  end if

  if ( limit == 1 ) then
    ier = 1
  end if

  if ( ier /= 0 .or. (abserr <= errbnd .and. abserr /= resabs ) .or. &
    abserr == 0.0E+00 ) go to 140
!
!  Initialization.
!
  rlist2(1) = result
  errmax = abserr
  maxerr = 1
  area = result
  errsum = abserr
  abserr = huge ( abserr )
  nrmax = 1
  nres = 0
  numrl2 = 2
  ktmin = 0
  extrap = .false.
  noext = .false.
  iroff1 = 0
  iroff2 = 0
  iroff3 = 0

  if ( dres >= (1.0E+00 - 5.0E+01* epsilon ( defabs ) ) * defabs ) then
    ksgn = 1
  else
    ksgn = -1
  end if

  do last = 2, limit
!
!  Bisect the subinterval with the nrmax-th largest error estimate.
!
    a1 = alist(maxerr)
    b1 = 5.0E-01 * ( alist(maxerr) + blist(maxerr) )
    a2 = b1
    b2 = blist(maxerr)
    erlast = errmax
    call qk21_x ( f_ptr, dat, a1, b1, area1, error1, resabs, defab1 )
    call qk21_x ( f_ptr, dat, a2, b2, area2, error2, resabs, defab2 )
!
!  Improve previous approximations to integral and error
!  and test for accuracy.
!
    area12 = area1+area2
    erro12 = error1+error2
    errsum = errsum+erro12-errmax
    area = area+area12-rlist(maxerr)

    if ( defab1 == error1 .or. defab2 == error2 ) go to 15

    if ( abs ( rlist(maxerr) - area12) > 1.0E-05 * abs(area12) &
      .or. erro12 < 9.9E-01 * errmax ) go to 10

    if ( extrap ) then
      iroff2 = iroff2+1
    else
      iroff1 = iroff1+1
    end if

10  continue

    if ( last > 10 .and. erro12 > errmax ) then
      iroff3 = iroff3+1
    end if

15  continue

    rlist(maxerr) = area1
    rlist(last) = area2
    errbnd = max ( epsabs, epsrel*abs(area) )
!
!  Test for roundoff error and eventually set error flag.
!
    if ( iroff1+iroff2 >= 10 .or. iroff3 >= 20 ) then
      ier = 2
    end if

    if ( iroff2 >= 5 ) then
      ierro = 3
    end if
!
!  Set error flag in the case that the number of subintervals
!  equals limit.
!
    if ( last == limit ) then
      ier = 1
    end if
!
!  Set error flag in the case of bad integrand behavior
!  at a point of the integration range.
!
    if ( max ( abs(a1),abs(b2)) <= (1.0E+00+1.0E+03* epsilon ( a1 ) )* &
      (abs(a2)+1.0E+03* tiny ( a2 ) ) ) then
      ier = 4
    end if
!
!  Append the newly-created intervals to the list.
!
    if ( error2 <= error1 ) then
      alist(last) = a2
      blist(maxerr) = b1
      blist(last) = b2
      elist(maxerr) = error1
      elist(last) = error2
    else
      alist(maxerr) = a2
      alist(last) = a1
      blist(last) = b1
      rlist(maxerr) = area2
      rlist(last) = area1
      elist(maxerr) = error2
      elist(last) = error1
    end if
!
!  Call QSORT to maintain the descending ordering
!  in the list of error estimates and select the subinterval
!  with nrmax-th largest error estimate (to be bisected next).
!
    call qsort ( limit, last, maxerr, errmax, elist, iord, nrmax )

    if ( errsum <= errbnd ) go to 115

    if ( ier /= 0 ) then
      exit
    end if

    if ( last == 2 ) go to 80
    if ( noext ) go to 90

    erlarg = erlarg-erlast

    if ( abs(b1-a1) > small ) then
      erlarg = erlarg+erro12
    end if
!
!  Test whether the interval to be bisected next is the
!  smallest interval.
!
    if ( .not. extrap ) then
      if ( abs(blist(maxerr)-alist(maxerr)) > small ) go to 90
      extrap = .true.
      nrmax = 2
    end if

!40  continue
!
!  The smallest interval has the largest error.
!  Before bisecting decrease the sum of the errors over the
!  larger intervals (erlarg) and perform extrapolation.
!
    if ( ierro /= 3 .and. erlarg > ertest ) then

      id = nrmax
      jupbnd = last

      if ( last > (2+limit/2) ) then
        jupbnd = limit+3-last
      end if

      do k = id, jupbnd
        maxerr = iord(nrmax)
        errmax = elist(maxerr)
        if ( abs(blist(maxerr)-alist(maxerr)) > small ) then
          go to 90
        end if
        nrmax = nrmax+1
      end do

    end if
!
!  Perform extrapolation.
!
!60  continue

    numrl2 = numrl2+1
    rlist2(numrl2) = area
    call qextr ( numrl2, rlist2, reseps, abseps, res3la, nres )
    ktmin = ktmin+1

    if ( ktmin > 5 .and. abserr < 1.0E-03 * errsum ) then
      ier = 5
    end if

    if ( abseps < abserr ) then

      ktmin = 0
      abserr = abseps
      result = reseps
      correc = erlarg
      ertest = max ( epsabs,epsrel*abs(reseps))

      if ( abserr <= ertest ) then
        exit
      end if

    end if
!
!  Prepare bisection of the smallest interval.
!
    if ( numrl2 == 1 ) then
      noext = .true.
    end if

    if ( ier == 5 ) then
      exit
    end if

    maxerr = iord(1)
    errmax = elist(maxerr)
    nrmax = 1
    extrap = .false.
    small = small * 5.0E-01
    erlarg = errsum
    go to 90

80  continue

    small = abs ( b - a ) * 3.75E-01
    erlarg = errsum
    ertest = errbnd
    rlist2(2) = area

90  continue

  end do
!
!  Set final result and error estimate.
!
  if ( abserr == huge ( abserr ) ) then
    go to 115
  end if

  if ( ier + ierro == 0 ) then
    go to 110
  end if

  if ( ierro == 3 ) then
    abserr = abserr + correc
  end if

  if ( ier == 0 ) then
    ier = 3
  end if

  if ( result /= 0.0E+00 .and. area /= 0.0E+00 ) then
    go to 105
  end if

  if ( abserr > errsum ) go to 115
  if ( area == 0.0E+00 ) go to 130
  go to 110

105 continue

  if ( abserr/abs(result) > errsum/abs(area) ) go to 115
!
!  Test on divergence.
!
110 continue

  if ( ksgn == (-1) .and. max ( abs(result),abs(area)) <=  &
   defabs*1.0E-02 ) go to 130

  if ( 1.0E-02 > (result/area) .or. (result/area) > 1.0E+02 &
   .or. errsum > abs(area) ) then
    ier = 6
  end if

  go to 130
!
!  Compute global integral sum.
!
115 continue

  result = sum ( rlist(1:last) )

  abserr = errsum

130 continue
 
  if ( 2 < ier ) then
    ier = ier - 1
  end if

140 continue

  neval = 42*last-21

  return
end subroutine qags_x

subroutine qags_x_vec ( f_ptr, dat, a, b, epsabs, epsrel, result, abserr, neval, ier )

!*****************************************************************************80
!
!! QAGS estimates the integral of a function.
!
!  Discussion:
!
!    The routine calculates an approximation RESULT to a definite integral   
!      I = integral of F over (A,B),
!    hopefully satisfying
!      || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real ( kind = 8 ) F, the name of the function routine, of the form
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!    which evaluates the integrand function.
!
!    Input, real ( kind = 8 ) A, B, the limits of integration.
!
!    Input, real ( kind = 8 ) EPSABS, EPSREL, the absolute and relative accuracy requested.
!
!    Output, real ( kind = 8 ) RESULT, the estimated value of the integral.
!
!    Output, real ( kind = 8 ) ABSERR, an estimate of || I - RESULT ||.
!
!    Output, integer ( kind = 8 ) NEVAL, the number of times the integral was evaluated.
!
!    Output, integer ( kind = 8 ) IER, error flag.
!                     ier = 0 normal and reliable termination of the
!                             routine. it is assumed that the requested
!                             accuracy has been achieved.
!                     ier > 0 abnormal termination of the routine
!                             the estimates for integral and error are
!                             less reliable. it is assumed that the
!                             requested accuracy has not been achieved.
!                         = 1 maximum number of subdivisions allowed
!                             has been achieved. one can allow more sub-
!                             divisions by increasing the data value of
!                             limit in qags (and taking the according
!                             dimension adjustments into account).
!                             however, if this yields no improvement
!                             it is advised to analyze the integrand
!                             in order to determine the integration
!                             difficulties. if the position of a
!                             local difficulty can be determined (e.g.
!                             singularity, discontinuity within the
!                             interval) one will probably gain from
!                             splitting up the interval at this point
!                             and calling the integrator on the sub-
!                             ranges. if possible, an appropriate
!                             special-purpose integrator should be used,
!                             which is designed for handling the type
!                             of difficulty involved.
!                         = 2 the occurrence of roundoff error is detec-
!                             ted, which prevents the requested
!                             tolerance from being achieved.
!                             the error may be under-estimated.
!                         = 3 extremely bad integrand behavior occurs
!                             at some  points of the integration
!                             interval.
!                         = 4 the algorithm does not converge. roundoff
!                             error is detected in the extrapolation
!                             table. it is presumed that the requested
!                             tolerance cannot be achieved, and that the
!                             returned result is the best which can be
!                             obtained.
!                         = 5 the integral is probably divergent, or
!                             slowly convergent. it must be noted that
!                             divergence can occur with any other value
!                             of ier.
!                         = 6 the input is invalid, because
!                             epsabs < 0 and epsrel < 0,
!                             result, abserr and neval are set to zero.
!
!  Local Parameters:
!
!           alist     - list of left end points of all subintervals
!                       considered up to now
!           blist     - list of right end points of all subintervals
!                       considered up to now
!           rlist(i)  - approximation to the integral over
!                       (alist(i),blist(i))
!           rlist2    - array of dimension at least limexp+2 containing
!                       the part of the epsilon table which is still
!                       needed for further computations
!           elist(i)  - error estimate applying to rlist(i)
!           maxerr    - pointer to the interval with largest error
!                       estimate
!           errmax    - elist(maxerr)
!           erlast    - error on the interval currently subdivided
!                       (before that subdivision has taken place)
!           area      - sum of the integrals over the subintervals
!           errsum    - sum of the errors over the subintervals
!           errbnd    - requested accuracy max(epsabs,epsrel*
!                       abs(result))
!           *****1    - variable for the left interval
!           *****2    - variable for the right interval
!           last      - index for subdivision
!           nres      - number of calls to the extrapolation routine
!           numrl2    - number of elements currently in rlist2. if an
!                       appropriate approximation to the compounded
!                       integral has been obtained it is put in
!                       rlist2(numrl2) after numrl2 has been increased
!                       by one.
!           small     - length of the smallest interval considered
!                       up to now, multiplied by 1.5
!           erlarg    - sum of the errors over the intervals larger
!                       than the smallest interval considered up to now
!           extrap    - logical variable denoting that the routine is
!                       attempting to perform extrapolation i.e. before
!                       subdividing the smallest interval we try to
!                       decrease the value of erlarg.
!           noext     - logical variable denoting that extrapolation
!                       is no longer allowed (true value)
!
  implicit none

  integer ( kind = 8 ), parameter :: limit = 500

  real ( kind = 8 ) a
  real ( kind = 8 ) abseps
  real ( kind = 8 ) abserr
  real ( kind = 8 ) alist(limit)
  real ( kind = 8 ) area
  real ( kind = 8 ) area1
  real ( kind = 8 ) area12
  real ( kind = 8 ) area2
  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  real ( kind = 8 ) b
  real ( kind = 8 ) blist(limit)
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) correc
  real ( kind = 8 ) defabs
  real ( kind = 8 ) defab1
  real ( kind = 8 ) defab2
  real ( kind = 8 ) dres
  real ( kind = 8 ) elist(limit)
  real ( kind = 8 ) epsabs
  real ( kind = 8 ) epsrel
  real ( kind = 8 ) erlarg
  real ( kind = 8 ) erlast
  real ( kind = 8 ) errbnd
  real ( kind = 8 ) errmax
  real ( kind = 8 ) error1
  real ( kind = 8 ) error2
  real ( kind = 8 ) erro12
  real ( kind = 8 ) errsum
  real ( kind = 8 ) ertest
  logical extrap
  !real ( kind = 8 ), external :: f
  procedure (f_int_dat_vec), intent(in), pointer :: f_ptr 
  integer ( kind = 8 ) id
  integer ( kind = 8 ) ier
  integer ( kind = 8 ) ierro
  integer ( kind = 8 ) iord(limit)
  integer ( kind = 8 ) iroff1
  integer ( kind = 8 ) iroff2
  integer ( kind = 8 ) iroff3
  integer ( kind = 8 ) jupbnd
  integer ( kind = 8 ) k
  integer ( kind = 8 ) ksgn
  integer ( kind = 8 ) ktmin
  integer ( kind = 8 ) last
  logical noext
  integer ( kind = 8 ) maxerr
  integer ( kind = 8 ) neval
  integer ( kind = 8 ) nres
  integer ( kind = 8 ) nrmax
  integer ( kind = 8 ) numrl2
  real ( kind = 8 ) resabs
  real ( kind = 8 ) reseps
  real ( kind = 8 ) result
  real ( kind = 8 ) res3la(3)
  real ( kind = 8 ) rlist(limit)
  real ( kind = 8 ) rlist2(52)
  real ( kind = 8 ) small
  
  class(dataCollectionBase), target :: dat
  
  !f_ptr => f
  
!
!  The dimension of rlist2 is determined by the value of
!  limexp in QEXTR (rlist2 should be of dimension
!  (limexp+2) at least).
!
!  Test on validity of parameters.
!
  ier = 0
  neval = 0
  last = 0
  result = 0.0E+00
  abserr = 0.0E+00
  alist(1) = a
  blist(1) = b
  rlist(1) = 0.0E+00
  elist(1) = 0.0E+00

  if ( epsabs < 0.0E+00 .and. epsrel < 0.0E+00 ) then
    ier = 6
    return
  end if
!
!  First approximation to the integral.
!
  ierro = 0
  call qk21_x_vec ( f_ptr, dat, a, b, result, abserr, defabs, resabs )
!
!  Test on accuracy.
!
  dres = abs ( result )
  errbnd = max ( epsabs, epsrel * dres )
  last = 1
  rlist(1) = result
  elist(1) = abserr
  iord(1) = 1

  if ( abserr <= 1.0E+02 * epsilon ( defabs ) * defabs .and. &
    abserr > errbnd ) then
    ier = 2
  end if

  if ( limit == 1 ) then
    ier = 1
  end if

  if ( ier /= 0 .or. (abserr <= errbnd .and. abserr /= resabs ) .or. &
    abserr == 0.0E+00 ) go to 140
!
!  Initialization.
!
  rlist2(1) = result
  errmax = abserr
  maxerr = 1
  area = result
  errsum = abserr
  abserr = huge ( abserr )
  nrmax = 1
  nres = 0
  numrl2 = 2
  ktmin = 0
  extrap = .false.
  noext = .false.
  iroff1 = 0
  iroff2 = 0
  iroff3 = 0

  if ( dres >= (1.0E+00 - 5.0E+01* epsilon ( defabs ) ) * defabs ) then
    ksgn = 1
  else
    ksgn = -1
  end if

  do last = 2, limit
!
!  Bisect the subinterval with the nrmax-th largest error estimate.
!
    a1 = alist(maxerr)
    b1 = 5.0E-01 * ( alist(maxerr) + blist(maxerr) )
    a2 = b1
    b2 = blist(maxerr)
    erlast = errmax
    call qk21_x_vec ( f_ptr, dat, a1, b1, area1, error1, resabs, defab1 )
    call qk21_x_vec ( f_ptr, dat, a2, b2, area2, error2, resabs, defab2 )
!
!  Improve previous approximations to integral and error
!  and test for accuracy.
!
    area12 = area1+area2
    erro12 = error1+error2
    errsum = errsum+erro12-errmax
    area = area+area12-rlist(maxerr)

    if ( defab1 == error1 .or. defab2 == error2 ) go to 15

    if ( abs ( rlist(maxerr) - area12) > 1.0E-05 * abs(area12) &
      .or. erro12 < 9.9E-01 * errmax ) go to 10

    if ( extrap ) then
      iroff2 = iroff2+1
    else
      iroff1 = iroff1+1
    end if

10  continue

    if ( last > 10 .and. erro12 > errmax ) then
      iroff3 = iroff3+1
    end if

15  continue

    rlist(maxerr) = area1
    rlist(last) = area2
    errbnd = max ( epsabs, epsrel*abs(area) )
!
!  Test for roundoff error and eventually set error flag.
!
    if ( iroff1+iroff2 >= 10 .or. iroff3 >= 20 ) then
      ier = 2
    end if

    if ( iroff2 >= 5 ) then
      ierro = 3
    end if
!
!  Set error flag in the case that the number of subintervals
!  equals limit.
!
    if ( last == limit ) then
      ier = 1
    end if
!
!  Set error flag in the case of bad integrand behavior
!  at a point of the integration range.
!
    if ( max ( abs(a1),abs(b2)) <= (1.0E+00+1.0E+03* epsilon ( a1 ) )* &
      (abs(a2)+1.0E+03* tiny ( a2 ) ) ) then
      ier = 4
    end if
!
!  Append the newly-created intervals to the list.
!
    if ( error2 <= error1 ) then
      alist(last) = a2
      blist(maxerr) = b1
      blist(last) = b2
      elist(maxerr) = error1
      elist(last) = error2
    else
      alist(maxerr) = a2
      alist(last) = a1
      blist(last) = b1
      rlist(maxerr) = area2
      rlist(last) = area1
      elist(maxerr) = error2
      elist(last) = error1
    end if
!
!  Call QSORT to maintain the descending ordering
!  in the list of error estimates and select the subinterval
!  with nrmax-th largest error estimate (to be bisected next).
!
    call qsort ( limit, last, maxerr, errmax, elist, iord, nrmax )

    if ( errsum <= errbnd ) go to 115

    if ( ier /= 0 ) then
      exit
    end if

    if ( last == 2 ) go to 80
    if ( noext ) go to 90

    erlarg = erlarg-erlast

    if ( abs(b1-a1) > small ) then
      erlarg = erlarg+erro12
    end if
!
!  Test whether the interval to be bisected next is the
!  smallest interval.
!
    if ( .not. extrap ) then
      if ( abs(blist(maxerr)-alist(maxerr)) > small ) go to 90
      extrap = .true.
      nrmax = 2
    end if

!40  continue
!
!  The smallest interval has the largest error.
!  Before bisecting decrease the sum of the errors over the
!  larger intervals (erlarg) and perform extrapolation.
!
    if ( ierro /= 3 .and. erlarg > ertest ) then

      id = nrmax
      jupbnd = last

      if ( last > (2+limit/2) ) then
        jupbnd = limit+3-last
      end if

      do k = id, jupbnd
        maxerr = iord(nrmax)
        errmax = elist(maxerr)
        if ( abs(blist(maxerr)-alist(maxerr)) > small ) then
          go to 90
        end if
        nrmax = nrmax+1
      end do

    end if
!
!  Perform extrapolation.
!
!60  continue

    numrl2 = numrl2+1
    rlist2(numrl2) = area
    call qextr ( numrl2, rlist2, reseps, abseps, res3la, nres )
    ktmin = ktmin+1

    if ( ktmin > 5 .and. abserr < 1.0E-03 * errsum ) then
      ier = 5
    end if

    if ( abseps < abserr ) then

      ktmin = 0
      abserr = abseps
      result = reseps
      correc = erlarg
      ertest = max ( epsabs,epsrel*abs(reseps))

      if ( abserr <= ertest ) then
        exit
      end if

    end if
!
!  Prepare bisection of the smallest interval.
!
    if ( numrl2 == 1 ) then
      noext = .true.
    end if

    if ( ier == 5 ) then
      exit
    end if

    maxerr = iord(1)
    errmax = elist(maxerr)
    nrmax = 1
    extrap = .false.
    small = small * 5.0E-01
    erlarg = errsum
    go to 90

80  continue

    small = abs ( b - a ) * 3.75E-01
    erlarg = errsum
    ertest = errbnd
    rlist2(2) = area

90  continue

  end do
!
!  Set final result and error estimate.
!
  if ( abserr == huge ( abserr ) ) then
    go to 115
  end if

  if ( ier + ierro == 0 ) then
    go to 110
  end if

  if ( ierro == 3 ) then
    abserr = abserr + correc
  end if

  if ( ier == 0 ) then
    ier = 3
  end if

  if ( result /= 0.0E+00 .and. area /= 0.0E+00 ) then
    go to 105
  end if

  if ( abserr > errsum ) go to 115
  if ( area == 0.0E+00 ) go to 130
  go to 110

105 continue

  if ( abserr/abs(result) > errsum/abs(area) ) go to 115
!
!  Test on divergence.
!
110 continue

  if ( ksgn == (-1) .and. max ( abs(result),abs(area)) <=  &
   defabs*1.0E-02 ) go to 130

  if ( 1.0E-02 > (result/area) .or. (result/area) > 1.0E+02 &
   .or. errsum > abs(area) ) then
    ier = 6
  end if

  go to 130
!
!  Compute global integral sum.
!
115 continue

  result = sum ( rlist(1:last) )

  abserr = errsum

130 continue
 
  if ( 2 < ier ) then
    ier = ier - 1
  end if

140 continue

  neval = 42*last-21

  return
end subroutine qags_x_vec

!::needed for doing double or triple integrals
subroutine qags_y ( f_ptr, dat, a, b, epsabs, epsrel, result, abserr, neval, ier )

!*****************************************************************************80
!
!! QAGS estimates the integral of a function.
!
!  Discussion:
!
!    The routine calculates an approximation RESULT to a definite integral   
!      I = integral of F over (A,B),
!    hopefully satisfying
!      || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real ( kind = 8 ) F, the name of the function routine, of the form
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!    which evaluates the integrand function.
!
!    Input, real ( kind = 8 ) A, B, the limits of integration.
!
!    Input, real ( kind = 8 ) EPSABS, EPSREL, the absolute and relative accuracy requested.
!
!    Output, real ( kind = 8 ) RESULT, the estimated value of the integral.
!
!    Output, real ( kind = 8 ) ABSERR, an estimate of || I - RESULT ||.
!
!    Output, integer ( kind = 8 ) NEVAL, the number of times the integral was evaluated.
!
!    Output, integer ( kind = 8 ) IER, error flag.
!                     ier = 0 normal and reliable termination of the
!                             routine. it is assumed that the requested
!                             accuracy has been achieved.
!                     ier > 0 abnormal termination of the routine
!                             the estimates for integral and error are
!                             less reliable. it is assumed that the
!                             requested accuracy has not been achieved.
!                         = 1 maximum number of subdivisions allowed
!                             has been achieved. one can allow more sub-
!                             divisions by increasing the data value of
!                             limit in qags (and taking the according
!                             dimension adjustments into account).
!                             however, if this yields no improvement
!                             it is advised to analyze the integrand
!                             in order to determine the integration
!                             difficulties. if the position of a
!                             local difficulty can be determined (e.g.
!                             singularity, discontinuity within the
!                             interval) one will probably gain from
!                             splitting up the interval at this point
!                             and calling the integrator on the sub-
!                             ranges. if possible, an appropriate
!                             special-purpose integrator should be used,
!                             which is designed for handling the type
!                             of difficulty involved.
!                         = 2 the occurrence of roundoff error is detec-
!                             ted, which prevents the requested
!                             tolerance from being achieved.
!                             the error may be under-estimated.
!                         = 3 extremely bad integrand behavior occurs
!                             at some  points of the integration
!                             interval.
!                         = 4 the algorithm does not converge. roundoff
!                             error is detected in the extrapolation
!                             table. it is presumed that the requested
!                             tolerance cannot be achieved, and that the
!                             returned result is the best which can be
!                             obtained.
!                         = 5 the integral is probably divergent, or
!                             slowly convergent. it must be noted that
!                             divergence can occur with any other value
!                             of ier.
!                         = 6 the input is invalid, because
!                             epsabs < 0 and epsrel < 0,
!                             result, abserr and neval are set to zero.
!
!  Local Parameters:
!
!           alist     - list of left end points of all subintervals
!                       considered up to now
!           blist     - list of right end points of all subintervals
!                       considered up to now
!           rlist(i)  - approximation to the integral over
!                       (alist(i),blist(i))
!           rlist2    - array of dimension at least limexp+2 containing
!                       the part of the epsilon table which is still
!                       needed for further computations
!           elist(i)  - error estimate applying to rlist(i)
!           maxerr    - pointer to the interval with largest error
!                       estimate
!           errmax    - elist(maxerr)
!           erlast    - error on the interval currently subdivided
!                       (before that subdivision has taken place)
!           area      - sum of the integrals over the subintervals
!           errsum    - sum of the errors over the subintervals
!           errbnd    - requested accuracy max(epsabs,epsrel*
!                       abs(result))
!           *****1    - variable for the left interval
!           *****2    - variable for the right interval
!           last      - index for subdivision
!           nres      - number of calls to the extrapolation routine
!           numrl2    - number of elements currently in rlist2. if an
!                       appropriate approximation to the compounded
!                       integral has been obtained it is put in
!                       rlist2(numrl2) after numrl2 has been increased
!                       by one.
!           small     - length of the smallest interval considered
!                       up to now, multiplied by 1.5
!           erlarg    - sum of the errors over the intervals larger
!                       than the smallest interval considered up to now
!           extrap    - logical variable denoting that the routine is
!                       attempting to perform extrapolation i.e. before
!                       subdividing the smallest interval we try to
!                       decrease the value of erlarg.
!           noext     - logical variable denoting that extrapolation
!                       is no longer allowed (true value)
!
  implicit none

  integer ( kind = 8 ), parameter :: limit = 500

  real ( kind = 8 ) a
  real ( kind = 8 ) abseps
  real ( kind = 8 ) abserr
  real ( kind = 8 ) alist(limit)
  real ( kind = 8 ) area
  real ( kind = 8 ) area1
  real ( kind = 8 ) area12
  real ( kind = 8 ) area2
  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  real ( kind = 8 ) b
  real ( kind = 8 ) blist(limit)
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) correc
  real ( kind = 8 ) defabs
  real ( kind = 8 ) defab1
  real ( kind = 8 ) defab2
  real ( kind = 8 ) dres
  real ( kind = 8 ) elist(limit)
  real ( kind = 8 ) epsabs
  real ( kind = 8 ) epsrel
  real ( kind = 8 ) erlarg
  real ( kind = 8 ) erlast
  real ( kind = 8 ) errbnd
  real ( kind = 8 ) errmax
  real ( kind = 8 ) error1
  real ( kind = 8 ) error2
  real ( kind = 8 ) erro12
  real ( kind = 8 ) errsum
  real ( kind = 8 ) ertest
  logical extrap
  !real ( kind = 8 ), external :: f
  procedure (f_int_dat),intent(in), pointer :: f_ptr 
  integer ( kind = 8 ) id
  integer ( kind = 8 ) ier
  integer ( kind = 8 ) ierro
  integer ( kind = 8 ) iord(limit)
  integer ( kind = 8 ) iroff1
  integer ( kind = 8 ) iroff2
  integer ( kind = 8 ) iroff3
  integer ( kind = 8 ) jupbnd
  integer ( kind = 8 ) k
  integer ( kind = 8 ) ksgn
  integer ( kind = 8 ) ktmin
  integer ( kind = 8 ) last
  logical noext
  integer ( kind = 8 ) maxerr
  integer ( kind = 8 ) neval
  integer ( kind = 8 ) nres
  integer ( kind = 8 ) nrmax
  integer ( kind = 8 ) numrl2
  real ( kind = 8 ) resabs
  real ( kind = 8 ) reseps
  real ( kind = 8 ) result
  real ( kind = 8 ) res3la(3)
  real ( kind = 8 ) rlist(limit)
  real ( kind = 8 ) rlist2(52)
  real ( kind = 8 ) small
  class(dataCollectionBase), target :: dat
!
!  The dimension of rlist2 is determined by the value of
!  limexp in QEXTR (rlist2 should be of dimension
!  (limexp+2) at least).
!
!  Test on validity of parameters.
!
 !f_ptr => f

  ier = 0
  neval = 0
  last = 0
  result = 0.0E+00
  abserr = 0.0E+00
  alist(1) = a
  blist(1) = b
  rlist(1) = 0.0E+00
  elist(1) = 0.0E+00

  if ( epsabs < 0.0E+00 .and. epsrel < 0.0E+00 ) then
    ier = 6
    return
  end if
!
!  First approximation to the integral.
!
  ierro = 0
  call qk21_y ( f_ptr, dat, a, b, result, abserr, defabs, resabs )
!
!  Test on accuracy.
!
  dres = abs ( result )
  errbnd = max ( epsabs, epsrel * dres )
  last = 1
  rlist(1) = result
  elist(1) = abserr
  iord(1) = 1

  if ( abserr <= 1.0E+02 * epsilon ( defabs ) * defabs .and. &
    abserr > errbnd ) then
    ier = 2
  end if

  if ( limit == 1 ) then
    ier = 1
  end if

  if ( ier /= 0 .or. (abserr <= errbnd .and. abserr /= resabs ) .or. &
    abserr == 0.0E+00 ) go to 140
!
!  Initialization.
!
  rlist2(1) = result
  errmax = abserr
  maxerr = 1
  area = result
  errsum = abserr
  abserr = huge ( abserr )
  nrmax = 1
  nres = 0
  numrl2 = 2
  ktmin = 0
  extrap = .false.
  noext = .false.
  iroff1 = 0
  iroff2 = 0
  iroff3 = 0

  if ( dres >= (1.0E+00 - 5.0E+01* epsilon ( defabs ) ) * defabs ) then
    ksgn = 1
  else
    ksgn = -1
  end if

  do last = 2, limit
!
!  Bisect the subinterval with the nrmax-th largest error estimate.
!
    a1 = alist(maxerr)
    b1 = 5.0E-01 * ( alist(maxerr) + blist(maxerr) )
    a2 = b1
    b2 = blist(maxerr)
    erlast = errmax
    call qk21_y ( f_ptr, dat, a1, b1, area1, error1, resabs, defab1 )
    call qk21_y ( f_ptr, dat, a2, b2, area2, error2, resabs, defab2 )
!
!  Improve previous approximations to integral and error
!  and test for accuracy.
!
    area12 = area1+area2
    erro12 = error1+error2
    errsum = errsum+erro12-errmax
    area = area+area12-rlist(maxerr)

    if ( defab1 == error1 .or. defab2 == error2 ) go to 15

    if ( abs ( rlist(maxerr) - area12) > 1.0E-05 * abs(area12) &
      .or. erro12 < 9.9E-01 * errmax ) go to 10

    if ( extrap ) then
      iroff2 = iroff2+1
    else
      iroff1 = iroff1+1
    end if

10  continue

    if ( last > 10 .and. erro12 > errmax ) then
      iroff3 = iroff3+1
    end if

15  continue

    rlist(maxerr) = area1
    rlist(last) = area2
    errbnd = max ( epsabs, epsrel*abs(area) )
!
!  Test for roundoff error and eventually set error flag.
!
    if ( iroff1+iroff2 >= 10 .or. iroff3 >= 20 ) then
      ier = 2
    end if

    if ( iroff2 >= 5 ) then
      ierro = 3
    end if
!
!  Set error flag in the case that the number of subintervals
!  equals limit.
!
    if ( last == limit ) then
      ier = 1
    end if
!
!  Set error flag in the case of bad integrand behavior
!  at a point of the integration range.
!
    if ( max ( abs(a1),abs(b2)) <= (1.0E+00+1.0E+03* epsilon ( a1 ) )* &
      (abs(a2)+1.0E+03* tiny ( a2 ) ) ) then
      ier = 4
    end if
!
!  Append the newly-created intervals to the list.
!
    if ( error2 <= error1 ) then
      alist(last) = a2
      blist(maxerr) = b1
      blist(last) = b2
      elist(maxerr) = error1
      elist(last) = error2
    else
      alist(maxerr) = a2
      alist(last) = a1
      blist(last) = b1
      rlist(maxerr) = area2
      rlist(last) = area1
      elist(maxerr) = error2
      elist(last) = error1
    end if
!
!  Call QSORT to maintain the descending ordering
!  in the list of error estimates and select the subinterval
!  with nrmax-th largest error estimate (to be bisected next).
!
    call qsort ( limit, last, maxerr, errmax, elist, iord, nrmax )

    if ( errsum <= errbnd ) go to 115

    if ( ier /= 0 ) then
      exit
    end if

    if ( last == 2 ) go to 80
    if ( noext ) go to 90

    erlarg = erlarg-erlast

    if ( abs(b1-a1) > small ) then
      erlarg = erlarg+erro12
    end if
!
!  Test whether the interval to be bisected next is the
!  smallest interval.
!
    if ( .not. extrap ) then
      if ( abs(blist(maxerr)-alist(maxerr)) > small ) go to 90
      extrap = .true.
      nrmax = 2
    end if

!40  continue
!
!  The smallest interval has the largest error.
!  Before bisecting decrease the sum of the errors over the
!  larger intervals (erlarg) and perform extrapolation.
!
    if ( ierro /= 3 .and. erlarg > ertest ) then

      id = nrmax
      jupbnd = last

      if ( last > (2+limit/2) ) then
        jupbnd = limit+3-last
      end if

      do k = id, jupbnd
        maxerr = iord(nrmax)
        errmax = elist(maxerr)
        if ( abs(blist(maxerr)-alist(maxerr)) > small ) then
          go to 90
        end if
        nrmax = nrmax+1
      end do

    end if
!
!  Perform extrapolation.
!
!60  continue

    numrl2 = numrl2+1
    rlist2(numrl2) = area
    call qextr ( numrl2, rlist2, reseps, abseps, res3la, nres )
    ktmin = ktmin+1

    if ( ktmin > 5 .and. abserr < 1.0E-03 * errsum ) then
      ier = 5
    end if

    if ( abseps < abserr ) then

      ktmin = 0
      abserr = abseps
      result = reseps
      correc = erlarg
      ertest = max ( epsabs,epsrel*abs(reseps))

      if ( abserr <= ertest ) then
        exit
      end if

    end if
!
!  Prepare bisection of the smallest interval.
!
    if ( numrl2 == 1 ) then
      noext = .true.
    end if

    if ( ier == 5 ) then
      exit
    end if

    maxerr = iord(1)
    errmax = elist(maxerr)
    nrmax = 1
    extrap = .false.
    small = small * 5.0E-01
    erlarg = errsum
    go to 90

80  continue

    small = abs ( b - a ) * 3.75E-01
    erlarg = errsum
    ertest = errbnd
    rlist2(2) = area

90  continue

  end do
!
!  Set final result and error estimate.
!
  if ( abserr == huge ( abserr ) ) then
    go to 115
  end if

  if ( ier + ierro == 0 ) then
    go to 110
  end if

  if ( ierro == 3 ) then
    abserr = abserr + correc
  end if

  if ( ier == 0 ) then
    ier = 3
  end if

  if ( result /= 0.0E+00 .and. area /= 0.0E+00 ) then
    go to 105
  end if

  if ( abserr > errsum ) go to 115
  if ( area == 0.0E+00 ) go to 130
  go to 110

105 continue

  if ( abserr/abs(result) > errsum/abs(area) ) go to 115
!
!  Test on divergence.
!
110 continue

  if ( ksgn == (-1) .and. max ( abs(result),abs(area)) <=  &
   defabs*1.0E-02 ) go to 130

  if ( 1.0E-02 > (result/area) .or. (result/area) > 1.0E+02 &
   .or. errsum > abs(area) ) then
    ier = 6
  end if

  go to 130
!
!  Compute global integral sum.
!
115 continue

  result = sum ( rlist(1:last) )

  abserr = errsum

130 continue
 
  if ( 2 < ier ) then
    ier = ier - 1
  end if

140 continue

  neval = 42*last-21

  return
end subroutine qags_y

subroutine qags_y_vec ( f_ptr, dat, a, b, epsabs, epsrel, result, abserr, neval, ier )

!*****************************************************************************80
!
!! QAGS estimates the integral of a function.
!
!  Discussion:
!
!    The routine calculates an approximation RESULT to a definite integral   
!      I = integral of F over (A,B),
!    hopefully satisfying
!      || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real ( kind = 8 ) F, the name of the function routine, of the form
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!    which evaluates the integrand function.
!
!    Input, real ( kind = 8 ) A, B, the limits of integration.
!
!    Input, real ( kind = 8 ) EPSABS, EPSREL, the absolute and relative accuracy requested.
!
!    Output, real ( kind = 8 ) RESULT, the estimated value of the integral.
!
!    Output, real ( kind = 8 ) ABSERR, an estimate of || I - RESULT ||.
!
!    Output, integer ( kind = 8 ) NEVAL, the number of times the integral was evaluated.
!
!    Output, integer ( kind = 8 ) IER, error flag.
!                     ier = 0 normal and reliable termination of the
!                             routine. it is assumed that the requested
!                             accuracy has been achieved.
!                     ier > 0 abnormal termination of the routine
!                             the estimates for integral and error are
!                             less reliable. it is assumed that the
!                             requested accuracy has not been achieved.
!                         = 1 maximum number of subdivisions allowed
!                             has been achieved. one can allow more sub-
!                             divisions by increasing the data value of
!                             limit in qags (and taking the according
!                             dimension adjustments into account).
!                             however, if this yields no improvement
!                             it is advised to analyze the integrand
!                             in order to determine the integration
!                             difficulties. if the position of a
!                             local difficulty can be determined (e.g.
!                             singularity, discontinuity within the
!                             interval) one will probably gain from
!                             splitting up the interval at this point
!                             and calling the integrator on the sub-
!                             ranges. if possible, an appropriate
!                             special-purpose integrator should be used,
!                             which is designed for handling the type
!                             of difficulty involved.
!                         = 2 the occurrence of roundoff error is detec-
!                             ted, which prevents the requested
!                             tolerance from being achieved.
!                             the error may be under-estimated.
!                         = 3 extremely bad integrand behavior occurs
!                             at some  points of the integration
!                             interval.
!                         = 4 the algorithm does not converge. roundoff
!                             error is detected in the extrapolation
!                             table. it is presumed that the requested
!                             tolerance cannot be achieved, and that the
!                             returned result is the best which can be
!                             obtained.
!                         = 5 the integral is probably divergent, or
!                             slowly convergent. it must be noted that
!                             divergence can occur with any other value
!                             of ier.
!                         = 6 the input is invalid, because
!                             epsabs < 0 and epsrel < 0,
!                             result, abserr and neval are set to zero.
!
!  Local Parameters:
!
!           alist     - list of left end points of all subintervals
!                       considered up to now
!           blist     - list of right end points of all subintervals
!                       considered up to now
!           rlist(i)  - approximation to the integral over
!                       (alist(i),blist(i))
!           rlist2    - array of dimension at least limexp+2 containing
!                       the part of the epsilon table which is still
!                       needed for further computations
!           elist(i)  - error estimate applying to rlist(i)
!           maxerr    - pointer to the interval with largest error
!                       estimate
!           errmax    - elist(maxerr)
!           erlast    - error on the interval currently subdivided
!                       (before that subdivision has taken place)
!           area      - sum of the integrals over the subintervals
!           errsum    - sum of the errors over the subintervals
!           errbnd    - requested accuracy max(epsabs,epsrel*
!                       abs(result))
!           *****1    - variable for the left interval
!           *****2    - variable for the right interval
!           last      - index for subdivision
!           nres      - number of calls to the extrapolation routine
!           numrl2    - number of elements currently in rlist2. if an
!                       appropriate approximation to the compounded
!                       integral has been obtained it is put in
!                       rlist2(numrl2) after numrl2 has been increased
!                       by one.
!           small     - length of the smallest interval considered
!                       up to now, multiplied by 1.5
!           erlarg    - sum of the errors over the intervals larger
!                       than the smallest interval considered up to now
!           extrap    - logical variable denoting that the routine is
!                       attempting to perform extrapolation i.e. before
!                       subdividing the smallest interval we try to
!                       decrease the value of erlarg.
!           noext     - logical variable denoting that extrapolation
!                       is no longer allowed (true value)
!
  implicit none

  integer ( kind = 8 ), parameter :: limit = 500

  real ( kind = 8 ) a
  real ( kind = 8 ) abseps
  real ( kind = 8 ) abserr
  real ( kind = 8 ) alist(limit)
  real ( kind = 8 ) area
  real ( kind = 8 ) area1
  real ( kind = 8 ) area12
  real ( kind = 8 ) area2
  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  real ( kind = 8 ) b
  real ( kind = 8 ) blist(limit)
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) correc
  real ( kind = 8 ) defabs
  real ( kind = 8 ) defab1
  real ( kind = 8 ) defab2
  real ( kind = 8 ) dres
  real ( kind = 8 ) elist(limit)
  real ( kind = 8 ) epsabs
  real ( kind = 8 ) epsrel
  real ( kind = 8 ) erlarg
  real ( kind = 8 ) erlast
  real ( kind = 8 ) errbnd
  real ( kind = 8 ) errmax
  real ( kind = 8 ) error1
  real ( kind = 8 ) error2
  real ( kind = 8 ) erro12
  real ( kind = 8 ) errsum
  real ( kind = 8 ) ertest
  logical extrap
  !real ( kind = 8 ), external :: f
  procedure (f_int_dat_vec), intent(in), pointer :: f_ptr 
  integer ( kind = 8 ) id
  integer ( kind = 8 ) ier
  integer ( kind = 8 ) ierro
  integer ( kind = 8 ) iord(limit)
  integer ( kind = 8 ) iroff1
  integer ( kind = 8 ) iroff2
  integer ( kind = 8 ) iroff3
  integer ( kind = 8 ) jupbnd
  integer ( kind = 8 ) k
  integer ( kind = 8 ) ksgn
  integer ( kind = 8 ) ktmin
  integer ( kind = 8 ) last
  logical noext
  integer ( kind = 8 ) maxerr
  integer ( kind = 8 ) neval
  integer ( kind = 8 ) nres
  integer ( kind = 8 ) nrmax
  integer ( kind = 8 ) numrl2
  real ( kind = 8 ) resabs
  real ( kind = 8 ) reseps
  real ( kind = 8 ) result
  real ( kind = 8 ) res3la(3)
  real ( kind = 8 ) rlist(limit)
  real ( kind = 8 ) rlist2(52)
  real ( kind = 8 ) small
  
  class(dataCollectionBase), target :: dat
  
  !f_ptr => f
  
!
!  The dimension of rlist2 is determined by the value of
!  limexp in QEXTR (rlist2 should be of dimension
!  (limexp+2) at least).
!
!  Test on validity of parameters.
!
  ier = 0
  neval = 0
  last = 0
  result = 0.0E+00
  abserr = 0.0E+00
  alist(1) = a
  blist(1) = b
  rlist(1) = 0.0E+00
  elist(1) = 0.0E+00

  if ( epsabs < 0.0E+00 .and. epsrel < 0.0E+00 ) then
    ier = 6
    return
  end if
!
!  First approximation to the integral.
!
  ierro = 0
  call qk21_y_vec ( f_ptr, dat, a, b, result, abserr, defabs, resabs )
!
!  Test on accuracy.
!
  dres = abs ( result )
  errbnd = max ( epsabs, epsrel * dres )
  last = 1
  rlist(1) = result
  elist(1) = abserr
  iord(1) = 1

  if ( abserr <= 1.0E+02 * epsilon ( defabs ) * defabs .and. &
    abserr > errbnd ) then
    ier = 2
  end if

  if ( limit == 1 ) then
    ier = 1
  end if

  if ( ier /= 0 .or. (abserr <= errbnd .and. abserr /= resabs ) .or. &
    abserr == 0.0E+00 ) go to 140
!
!  Initialization.
!
  rlist2(1) = result
  errmax = abserr
  maxerr = 1
  area = result
  errsum = abserr
  abserr = huge ( abserr )
  nrmax = 1
  nres = 0
  numrl2 = 2
  ktmin = 0
  extrap = .false.
  noext = .false.
  iroff1 = 0
  iroff2 = 0
  iroff3 = 0

  if ( dres >= (1.0E+00 - 5.0E+01* epsilon ( defabs ) ) * defabs ) then
    ksgn = 1
  else
    ksgn = -1
  end if

  do last = 2, limit
!
!  Bisect the subinterval with the nrmax-th largest error estimate.
!
    a1 = alist(maxerr)
    b1 = 5.0E-01 * ( alist(maxerr) + blist(maxerr) )
    a2 = b1
    b2 = blist(maxerr)
    erlast = errmax
    call qk21_y_vec ( f_ptr, dat, a1, b1, area1, error1, resabs, defab1 )
    call qk21_y_vec ( f_ptr, dat, a2, b2, area2, error2, resabs, defab2 )
!
!  Improve previous approximations to integral and error
!  and test for accuracy.
!
    area12 = area1+area2
    erro12 = error1+error2
    errsum = errsum+erro12-errmax
    area = area+area12-rlist(maxerr)

    if ( defab1 == error1 .or. defab2 == error2 ) go to 15

    if ( abs ( rlist(maxerr) - area12) > 1.0E-05 * abs(area12) &
      .or. erro12 < 9.9E-01 * errmax ) go to 10

    if ( extrap ) then
      iroff2 = iroff2+1
    else
      iroff1 = iroff1+1
    end if

10  continue

    if ( last > 10 .and. erro12 > errmax ) then
      iroff3 = iroff3+1
    end if

15  continue

    rlist(maxerr) = area1
    rlist(last) = area2
    errbnd = max ( epsabs, epsrel*abs(area) )
!
!  Test for roundoff error and eventually set error flag.
!
    if ( iroff1+iroff2 >= 10 .or. iroff3 >= 20 ) then
      ier = 2
    end if

    if ( iroff2 >= 5 ) then
      ierro = 3
    end if
!
!  Set error flag in the case that the number of subintervals
!  equals limit.
!
    if ( last == limit ) then
      ier = 1
    end if
!
!  Set error flag in the case of bad integrand behavior
!  at a point of the integration range.
!
    if ( max ( abs(a1),abs(b2)) <= (1.0E+00+1.0E+03* epsilon ( a1 ) )* &
      (abs(a2)+1.0E+03* tiny ( a2 ) ) ) then
      ier = 4
    end if
!
!  Append the newly-created intervals to the list.
!
    if ( error2 <= error1 ) then
      alist(last) = a2
      blist(maxerr) = b1
      blist(last) = b2
      elist(maxerr) = error1
      elist(last) = error2
    else
      alist(maxerr) = a2
      alist(last) = a1
      blist(last) = b1
      rlist(maxerr) = area2
      rlist(last) = area1
      elist(maxerr) = error2
      elist(last) = error1
    end if
!
!  Call QSORT to maintain the descending ordering
!  in the list of error estimates and select the subinterval
!  with nrmax-th largest error estimate (to be bisected next).
!
    call qsort ( limit, last, maxerr, errmax, elist, iord, nrmax )

    if ( errsum <= errbnd ) go to 115

    if ( ier /= 0 ) then
      exit
    end if

    if ( last == 2 ) go to 80
    if ( noext ) go to 90

    erlarg = erlarg-erlast

    if ( abs(b1-a1) > small ) then
      erlarg = erlarg+erro12
    end if
!
!  Test whether the interval to be bisected next is the
!  smallest interval.
!
    if ( .not. extrap ) then
      if ( abs(blist(maxerr)-alist(maxerr)) > small ) go to 90
      extrap = .true.
      nrmax = 2
    end if

!40  continue
!
!  The smallest interval has the largest error.
!  Before bisecting decrease the sum of the errors over the
!  larger intervals (erlarg) and perform extrapolation.
!
    if ( ierro /= 3 .and. erlarg > ertest ) then

      id = nrmax
      jupbnd = last

      if ( last > (2+limit/2) ) then
        jupbnd = limit+3-last
      end if

      do k = id, jupbnd
        maxerr = iord(nrmax)
        errmax = elist(maxerr)
        if ( abs(blist(maxerr)-alist(maxerr)) > small ) then
          go to 90
        end if
        nrmax = nrmax+1
      end do

    end if
!
!  Perform extrapolation.
!
!60  continue

    numrl2 = numrl2+1
    rlist2(numrl2) = area
    call qextr ( numrl2, rlist2, reseps, abseps, res3la, nres )
    ktmin = ktmin+1

    if ( ktmin > 5 .and. abserr < 1.0E-03 * errsum ) then
      ier = 5
    end if

    if ( abseps < abserr ) then

      ktmin = 0
      abserr = abseps
      result = reseps
      correc = erlarg
      ertest = max ( epsabs,epsrel*abs(reseps))

      if ( abserr <= ertest ) then
        exit
      end if

    end if
!
!  Prepare bisection of the smallest interval.
!
    if ( numrl2 == 1 ) then
      noext = .true.
    end if

    if ( ier == 5 ) then
      exit
    end if

    maxerr = iord(1)
    errmax = elist(maxerr)
    nrmax = 1
    extrap = .false.
    small = small * 5.0E-01
    erlarg = errsum
    go to 90

80  continue

    small = abs ( b - a ) * 3.75E-01
    erlarg = errsum
    ertest = errbnd
    rlist2(2) = area

90  continue

  end do
!
!  Set final result and error estimate.
!
  if ( abserr == huge ( abserr ) ) then
    go to 115
  end if

  if ( ier + ierro == 0 ) then
    go to 110
  end if

  if ( ierro == 3 ) then
    abserr = abserr + correc
  end if

  if ( ier == 0 ) then
    ier = 3
  end if

  if ( result /= 0.0E+00 .and. area /= 0.0E+00 ) then
    go to 105
  end if

  if ( abserr > errsum ) go to 115
  if ( area == 0.0E+00 ) go to 130
  go to 110

105 continue

  if ( abserr/abs(result) > errsum/abs(area) ) go to 115
!
!  Test on divergence.
!
110 continue

  if ( ksgn == (-1) .and. max ( abs(result),abs(area)) <=  &
   defabs*1.0E-02 ) go to 130

  if ( 1.0E-02 > (result/area) .or. (result/area) > 1.0E+02 &
   .or. errsum > abs(area) ) then
    ier = 6
  end if

  go to 130
!
!  Compute global integral sum.
!
115 continue

  result = sum ( rlist(1:last) )

  abserr = errsum

130 continue
 
  if ( 2 < ier ) then
    ier = ier - 1
  end if

140 continue

  neval = 42*last-21

  return
end subroutine qags_y_vec

subroutine qawc ( f_ptr, dat, a, b, c, epsabs, epsrel, result, abserr, neval, ier )

!*****************************************************************************80
!
!! QAWC computes a Cauchy principal value.
!
!  Discussion:
!
!    The routine calculates an approximation RESULT to a Cauchy principal
!    value 
!      I = integral of F*W over (A,B),
!    with
!      W(X) = 1 / (X-C),
!    with C distinct from A and B, hopefully satisfying
!      || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real ( kind = 8 ) F, the name of the function routine, of the form
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!    which evaluates the integrand function.
!
!    Input, real ( kind = 8 ) A, B, the limits of integration.
!
!    Input, real ( kind = 8 ) C, a parameter in the weight function, which must
!    not be equal to A or B.
!
!    Input, real ( kind = 8 ) EPSABS, EPSREL, the absolute and relative accuracy requested.
!
!    Output, real ( kind = 8 ) RESULT, the estimated value of the integral.
!
!    Output, real ( kind = 8 ) ABSERR, an estimate of || I - RESULT ||.
!
!    Output, integer ( kind = 8 ) NEVAL, the number of times the integral was evaluated.
!
!            ier    - integer ( kind = 8 )
!                     ier = 0 normal and reliable termination of the
!                             routine. it is assumed that the requested
!                             accuracy has been achieved.
!                     ier > 0 abnormal termination of the routine
!                             the estimates for integral and error are
!                             less reliable. it is assumed that the
!                             requested accuracy has not been achieved.
!                     ier = 1 maximum number of subdivisions allowed
!                             has been achieved. one can allow more sub-
!                             divisions by increasing the data value of
!                             limit in qawc (and taking the according
!                             dimension adjustments into account).
!                             however, if this yields no improvement it
!                             is advised to analyze the integrand in
!                             order to determine the integration
!                             difficulties. if the position of a local
!                             difficulty can be determined (e.g.
!                             singularity, discontinuity within the
!                             interval one will probably gain from
!                             splitting up the interval at this point
!                             and calling appropriate integrators on the
!                             subranges.
!                         = 2 the occurrence of roundoff error is detec-
!                             ted, which prevents the requested
!                             tolerance from being achieved.
!                         = 3 extremely bad integrand behavior occurs
!                             at some points of the integration
!                             interval.
!                         = 6 the input is invalid, because
!                             c = a or c = b or
!                             epsabs < 0 and epsrel < 0,
!                             result, abserr, neval are set to zero.
!
!  Local parameters:
!
!    LIMIT is the maximum number of subintervals allowed in the
!    subdivision process of qawce. take care that limit >= 1.
!
implicit none
  class( dataCollectionBase ), target :: dat
  

  integer ( kind = 8 ), parameter :: limit = 500

  real ( kind = 8 ) a
  real ( kind = 8 ) abserr
  real ( kind = 8 ) alist(limit)
  real ( kind = 8 ) b
  real ( kind = 8 ) blist(limit)
  real ( kind = 8 ) elist(limit)
  real ( kind = 8 ) c
  real ( kind = 8 ) epsabs
  real ( kind = 8 ) epsrel
  !real ( kind = 8 ), external :: f
  procedure (f_int_dat), intent(in), pointer :: f_ptr 
  integer ( kind = 8 ) ier
  integer ( kind = 8 ) iord(limit)
  integer ( kind = 8 ) last
  integer ( kind = 8 ) neval
  real ( kind = 8 ) result
  real ( kind = 8 ) rlist(limit)

  call qawce ( f_ptr, dat, a, b, c, epsabs, epsrel, limit, result, abserr, neval, ier, &
    alist, blist, rlist, elist, iord, last )

  return
end subroutine
subroutine qawce ( f_ptr, dat, a, b, c, epsabs, epsrel, limit, result, abserr, neval, &
  ier, alist, blist, rlist, elist, iord, last )

!*****************************************************************************80
!
!! QAWCE computes a Cauchy principal value.
!
!  Discussion:
!
!    The routine calculates an approximation RESULT to a Cauchy principal
!    value   
!      I = integral of F*W over (A,B),
!    with
!      W(X) = 1 / ( X - C ),
!    with C distinct from A and B, hopefully satisfying
!      | I - RESULT | <= max ( EPSABS, EPSREL * |I| ).
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real ( kind = 8 ) F, the name of the function routine, of the form
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!    which evaluates the integrand function.
!
!    Input, real ( kind = 8 ) A, B, the limits of integration.
!
!    Input, real ( kind = 8 ) C, a parameter in the weight function, which cannot be
!    equal to A or B.
!
!    Input, real ( kind = 8 ) EPSABS, EPSREL, the absolute and relative accuracy requested.
!
!    Input, integer ( kind = 8 ) LIMIT, the upper bound on the number of subintervals that
!    will be used in the partition of [A,B].  LIMIT is typically 500.
!
!    Output, real ( kind = 8 ) RESULT, the estimated value of the integral.
!
!    Output, real ( kind = 8 ) ABSERR, an estimate of || I - RESULT ||.
!
!    Output, integer ( kind = 8 ) NEVAL, the number of times the integral was evaluated.
!
!            ier    - integer ( kind = 8 )
!                     ier = 0 normal and reliable termination of the
!                             routine. it is assumed that the requested
!                             accuracy has been achieved.
!                     ier > 0 abnormal termination of the routine
!                             the estimates for integral and error are
!                             less reliable. it is assumed that the
!                             requested accuracy has not been achieved.
!                     ier = 1 maximum number of subdivisions allowed
!                             has been achieved. one can allow more sub-
!                             divisions by increasing the value of
!                             limit. however, if this yields no
!                             improvement it is advised to analyze the
!                             integrand, in order to determine the
!                             integration difficulties.  if the position
!                             of a local difficulty can be determined
!                             (e.g. singularity, discontinuity within
!                             the interval) one will probably gain
!                             from splitting up the interval at this
!                             point and calling appropriate integrators
!                             on the subranges.
!                         = 2 the occurrence of roundoff error is detec-
!                             ted, which prevents the requested
!                             tolerance from being achieved.
!                         = 3 extremely bad integrand behavior occurs
!                             at some interior points of the integration
!                             interval.
!                         = 6 the input is invalid, because
!                             c = a or c = b or
!                             epsabs < 0 and epsrel < 0,
!                             or limit < 1.
!                             result, abserr, neval, rlist(1), elist(1),
!                             iord(1) and last are set to zero.
!                             alist(1) and blist(1) are set to a and b
!                             respectively.
!
!    Workspace, real ( kind = 8 ) ALIST(LIMIT), BLIST(LIMIT), contains in entries 1 
!    through LAST the left and right ends of the partition subintervals.
!
!    Workspace, real ( kind = 8 ) RLIST(LIMIT), contains in entries 1 through LAST
!    the integral approximations on the subintervals.
!
!    Workspace, real ( kind = 8 ) ELIST(LIMIT), contains in entries 1 through LAST
!    the absolute error estimates on the subintervals.
!
!            iord    - integer ( kind = 8 )
!                      vector of dimension at least limit, the first k
!                      elements of which are pointers to the error
!                      estimates over the subintervals, so that
!                      elist(iord(1)), ...,  elist(iord(k)) with
!                      k = last if last <= (limit/2+2), and
!                      k = limit+1-last otherwise, form a decreasing
!                      sequence.
!
!            last    - integer ( kind = 8 )
!                      number of subintervals actually produced in
!                      the subdivision process
!
!  Local parameters:
!
!           alist     - list of left end points of all subintervals
!                       considered up to now
!           blist     - list of right end points of all subintervals
!                       considered up to now
!           rlist(i)  - approximation to the integral over
!                       (alist(i),blist(i))
!           elist(i)  - error estimate applying to rlist(i)
!           maxerr    - pointer to the interval with largest error
!                       estimate
!           errmax    - elist(maxerr)
!           area      - sum of the integrals over the subintervals
!           errsum    - sum of the errors over the subintervals
!           errbnd    - requested accuracy max(epsabs,epsrel*
!                       abs(result))
!           *****1    - variable for the left subinterval
!           *****2    - variable for the right subinterval
!           last      - index for subdivision
!
  implicit none

  integer ( kind = 8 ) limit

  real ( kind = 8 ) a
  real ( kind = 8 ) aa
  real ( kind = 8 ) abserr
  real ( kind = 8 ) alist(limit)
  real ( kind = 8 ) area
  real ( kind = 8 ) area1
  real ( kind = 8 ) area12
  real ( kind = 8 ) area2
  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  real ( kind = 8 ) b
  real ( kind = 8 ) bb
  real ( kind = 8 ) blist(limit)
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) c
  real ( kind = 8 ) elist(limit)
  real ( kind = 8 ) epsabs
  real ( kind = 8 ) epsrel
  real ( kind = 8 ) errbnd
  real ( kind = 8 ) errmax
  real ( kind = 8 ) error1
  real ( kind = 8 ) error2
  real ( kind = 8 ) erro12
  real ( kind = 8 ) errsum
  !real ( kind = 8 ), external :: f
  procedure (f_int_dat), intent(in), pointer :: f_ptr 
  integer ( kind = 8 ) ier
  integer ( kind = 8 ) iord(limit)
  integer ( kind = 8 ) iroff1
  integer ( kind = 8 ) iroff2
  integer ( kind = 8 ) krule
  integer ( kind = 8 ) last
  integer ( kind = 8 ) maxerr
  integer ( kind = 8 ) nev
  integer ( kind = 8 ) neval
  integer ( kind = 8 ) nrmax
  real ( kind = 8 ) result
  real ( kind = 8 ) rlist(limit)
  class( dataCollectionBase ), target :: dat
  
  
!
!  Test on validity of parameters.
!
  ier = 0
  neval = 0
  last = 0
  alist(1) = a
  blist(1) = b
  rlist(1) = 0.0E+00
  elist(1) = 0.0E+00
  iord(1) = 0
  result = 0.0E+00
  abserr = 0.0E+00

  if ( c == a  ) then
    ier = 6
    return
  else if ( c == b ) then
    ier = 6
    return
  else if ( epsabs < 0.0E+00 .and. epsrel < 0.0E+00 ) then
    ier = 6
    return
  end if
!
!  First approximation to the integral.
!
  if ( a <= b ) then
    aa = a
    bb = b
  else
    aa = b
    bb = a
  end if

  krule = 1
  call qc25c ( f_ptr, dat, aa, bb, c, result, abserr, krule, neval )
  last = 1
  rlist(1) = result
  elist(1) = abserr
  iord(1) = 1
  alist(1) = a
  blist(1) = b
!
!  Test on accuracy.
!
  errbnd = max ( epsabs, epsrel * abs(result) )

  if ( limit == 1 ) then
    ier = 1
    go to 70
  end if

  if ( abserr < min ( 1.0E-02 * abs(result), errbnd)  ) then
    go to 70
  end if
!
!  Initialization
!
  alist(1) = aa
  blist(1) = bb
  rlist(1) = result
  errmax = abserr
  maxerr = 1
  area = result
  errsum = abserr
  nrmax = 1
  iroff1 = 0
  iroff2 = 0

  do last = 2, limit
!
!  Bisect the subinterval with nrmax-th largest error estimate.
!
    a1 = alist(maxerr)
    b1 = 5.0E-01*(alist(maxerr)+blist(maxerr))
    b2 = blist(maxerr)

    if ( c <= b1 .and. a1 < c ) then
      b1 = 5.0E-01*(c+b2)
    end if

    if ( b1 < c .and. c < b2 ) then
      b1 = 5.0E-01 * ( a1 + c )
    end if

    a2 = b1
    krule = 2

    call qc25c ( f_ptr, dat, a1, b1, c, area1, error1, krule, nev )
    neval = neval+nev

    call qc25c ( f_ptr, dat, a2, b2, c, area2, error2, krule, nev )
    neval = neval+nev
!
!  Improve previous approximations to integral and error
!  and test for accuracy.
!
    area12 = area1 + area2
    erro12 = error1 + error2
    errsum = errsum + erro12 - errmax
    area = area + area12 - rlist(maxerr)

    if ( abs ( rlist(maxerr)-area12) < 1.0E-05 * abs(area12) &
      .and. erro12 >= 9.9E-01 * errmax .and. krule == 0 ) &
      iroff1 = iroff1+1

    if ( last > 10.and.erro12 > errmax .and. krule == 0 ) then
      iroff2 = iroff2+1
    end if

    rlist(maxerr) = area1
    rlist(last) = area2
    errbnd = max ( epsabs, epsrel * abs(area) )

    if ( errsum > errbnd ) then
!
!  Test for roundoff error and eventually set error flag.
!
      if ( iroff1 >= 6 .and. iroff2 > 20 ) then
        ier = 2
      end if
!
!  Set error flag in the case that number of interval
!  bisections exceeds limit.
!
      if ( last == limit ) then
        ier = 1
      end if
!
!  Set error flag in the case of bad integrand behavior at
!  a point of the integration range.
!
      if ( max ( abs(a1), abs(b2) ) <= ( 1.0E+00 + 1.0E+03 * epsilon ( a1 ) ) &
        *( abs(a2)+1.0E+03* tiny ( a2 ) )) then
        ier = 3
      end if

    end if
!
!  Append the newly-created intervals to the list.
!
    if ( error2 <= error1 ) then
      alist(last) = a2
      blist(maxerr) = b1
      blist(last) = b2
      elist(maxerr) = error1
      elist(last) = error2
    else
      alist(maxerr) = a2
      alist(last) = a1
      blist(last) = b1
      rlist(maxerr) = area2
      rlist(last) = area1
      elist(maxerr) = error2
      elist(last) = error1
    end if
!
!  Call QSORT to maintain the descending ordering
!  in the list of error estimates and select the subinterval
!  with NRMAX-th largest error estimate (to be bisected next).
!
    call qsort ( limit, last, maxerr, errmax, elist, iord, nrmax )

    if ( ier /= 0 .or. errsum <= errbnd ) then
      exit
    end if

  end do
!
!  Compute final result.
!
  result = sum ( rlist(1:last) )

  abserr = errsum

70 continue 

  if ( aa == b ) then
    result = - result
  end if

  return
end subroutine
subroutine qawf ( f_ptr, dat, a, omega, integr, epsabs, result, abserr, neval, ier )

!*****************************************************************************80
!
!! QAWF computes Fourier integrals over the interval [ A, +Infinity ).
!
!  Discussion:
!
!    The routine calculates an approximation RESULT to a definite integral  
! 
!      I = integral of F*COS(OMEGA*X) 
!    or 
!      I = integral of F*SIN(OMEGA*X) 
!
!    over the interval [A,+Infinity), hopefully satisfying
!
!      || I - RESULT || <= EPSABS.
!
!    If OMEGA = 0 and INTEGR = 1, the integral is calculated by means 
!    of QAGI, and IER has the meaning as described in the comments of QAGI.
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real ( kind = 8 ) F, the name of the function routine, of the form
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!    which evaluates the integrand function.
!
!    Input, real ( kind = 8 ) A, the lower limit of integration.
!
!    Input, real ( kind = 8 ) OMEGA, the parameter in the weight function.
!
!    Input, integer ( kind = 8 ) INTEGR, indicates which weight functions is used
!    = 1, w(x) = cos(omega*x)
!    = 2, w(x) = sin(omega*x)
!
!    Input, real ( kind = 8 ) EPSABS, the absolute accuracy requested.
!
!    Output, real ( kind = 8 ) RESULT, the estimated value of the integral.
!
!    Output, real ( kind = 8 ) ABSERR, an estimate of || I - RESULT ||.
!
!    Output, integer ( kind = 8 ) NEVAL, the number of times the integral was evaluated.
!
!            ier    - integer ( kind = 8 )
!                     ier = 0 normal and reliable termination of the
!                             routine. it is assumed that the
!                             requested accuracy has been achieved.
!                     ier > 0 abnormal termination of the routine.
!                             the estimates for integral and error are
!                             less reliable. it is assumed that the
!                             requested accuracy has not been achieved.
!                    if omega /= 0
!                     ier = 6 the input is invalid because
!                             (integr /= 1 and integr /= 2) or
!                              epsabs <= 0
!                              result, abserr, neval, lst are set to
!                              zero.
!                         = 7 abnormal termination of the computation
!                             of one or more subintegrals
!                         = 8 maximum number of cycles allowed
!                             has been achieved, i.e. of subintervals
!                             (a+(k-1)c,a+kc) where
!                             c = (2*int(abs(omega))+1)*pi/abs(omega),
!                             for k = 1, 2, ...
!                         = 9 the extrapolation table constructed for
!                             convergence acceleration of the series
!                             formed by the integral contributions
!                             over the cycles, does not converge to
!                             within the requested accuracy.
!
!  Local parameters:
!
!    integer ( kind = 8 ) LIMLST, gives an upper bound on the number of cycles, LIMLST >= 3.
!    if limlst < 3, the routine will end with ier = 6.
!
!    integer ( kind = 8 ) MAXP1, an upper bound on the number of Chebyshev moments which 
!    can be stored, i.e. for the intervals of lengths abs(b-a)*2**(-l), 
!    l = 0,1, ..., maxp1-2, maxp1 >= 1.  if maxp1 < 1, the routine will end
!    with ier = 6.
!
  implicit none

  integer ( kind = 8 ), parameter :: limit = 500
  integer ( kind = 8 ), parameter :: limlst = 50
  integer ( kind = 8 ), parameter :: maxp1 = 21

  real ( kind = 8 ) a
  real ( kind = 8 ) abserr
  real ( kind = 8 ) alist(limit)
  real ( kind = 8 ) blist(limit)
  real ( kind = 8 ) chebmo(maxp1,25)
  real ( kind = 8 ) elist(limit)
  real ( kind = 8 ) epsabs
  real ( kind = 8 ) erlst(limlst)
  !real ( kind = 8 ), external :: f
  procedure (f_int_dat), intent(in), pointer :: f_ptr 
  integer ( kind = 8 ) ier
  integer ( kind = 8 ) integr
  integer ( kind = 8 ) iord(limit)
  integer ( kind = 8 ) ierlst(limlst)
  integer ( kind = 8 ) lst
  integer ( kind = 8 ) neval
  integer ( kind = 8 ) nnlog(limit)
  real ( kind = 8 ) omega
  real ( kind = 8 ) result
  real ( kind = 8 ) rlist(limit)
  real ( kind = 8 ) rslst(limlst)
  class( dataCollectionBase), target :: dat
  ier = 6
  neval = 0
  result = 0.0E+00
  abserr = 0.0E+00

  if ( limlst < 3 .or. maxp1 < 1 ) then
    return
  end if

  call qawfe ( f_ptr, dat, a, omega, integr, epsabs, limlst, limit, maxp1, &
    result, abserr, neval, ier, rslst, erlst, ierlst, lst, alist, blist, &
    rlist, elist, iord, nnlog, chebmo )

  return
end subroutine
subroutine qawfe ( f_ptr, dat, a, omega, integr, epsabs, limlst, limit, maxp1, &
  result, abserr, neval, ier, rslst, erlst, ierlst, lst, alist, blist, &
  rlist, elist, iord, nnlog, chebmo )

!*****************************************************************************80
!
!! QAWFE computes Fourier integrals.
!
!  Discussion:
!
!    The routine calculates an approximation RESULT to a definite integral   
!      I = integral of F*COS(OMEGA*X) or F*SIN(OMEGA*X) over (A,+Infinity),
!    hopefully satisfying
!      || I - RESULT || <= EPSABS.
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real ( kind = 8 ) F, the name of the function routine, of the form
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!    which evaluates the integrand function.
!
!    Input, real ( kind = 8 ) A, the lower limit of integration.
!
!    Input, real ( kind = 8 ) OMEGA, the parameter in the weight function.
!
!    Input, integer ( kind = 8 ) INTEGR, indicates which weight function is used
!    = 1      w(x) = cos(omega*x)
!    = 2      w(x) = sin(omega*x)
!
!    Input, real ( kind = 8 ) EPSABS, the absolute accuracy requested.
!
!    Input, integer ( kind = 8 ) LIMLST, an upper bound on the number of cycles.
!    LIMLST must be at least 1.  In fact, if LIMLST < 3, the routine 
!    will end with IER= 6.
!
!    Input, integer ( kind = 8 ) LIMIT, an upper bound on the number of subintervals 
!    allowed in the partition of each cycle, limit >= 1.
!
!            maxp1  - integer ( kind = 8 )
!                     gives an upper bound on the number of
!                     Chebyshev moments which can be stored, i.e.
!                     for the intervals of lengths abs(b-a)*2**(-l),
!                     l=0,1, ..., maxp1-2, maxp1 >= 1
!
!    Output, real ( kind = 8 ) RESULT, the estimated value of the integral.
!
!    Output, real ( kind = 8 ) ABSERR, an estimate of || I - RESULT ||.
!
!    Output, integer ( kind = 8 ) NEVAL, the number of times the integral was evaluated.
!
!            ier    - ier = 0 normal and reliable termination of
!                             the routine. it is assumed that the
!                             requested accuracy has been achieved.
!                     ier > 0 abnormal termination of the routine
!                             the estimates for integral and error
!                             are less reliable. it is assumed that
!                             the requested accuracy has not been
!                             achieved.
!                    if omega /= 0
!                     ier = 6 the input is invalid because
!                             (integr /= 1 and integr /= 2) or
!                              epsabs <= 0 or limlst < 3.
!                              result, abserr, neval, lst are set
!                              to zero.
!                         = 7 bad integrand behavior occurs within
!                             one or more of the cycles. location
!                             and type of the difficulty involved
!                             can be determined from the vector ierlst.
!                             here lst is the number of cycles actually
!                             needed (see below).
!                             ierlst(k) = 1 the maximum number of
!                                           subdivisions (= limit)
!                                           has been achieved on the
!                                           k th cycle.
!                                       = 2 occurence of roundoff
!                                           error is detected and
!                                           prevents the tolerance
!                                           imposed on the k th cycle
!                                           from being acheived.
!                                       = 3 extremely bad integrand
!                                           behavior occurs at some
!                                           points of the k th cycle.
!                                       = 4 the integration procedure
!                                           over the k th cycle does
!                                           not converge (to within the
!                                           required accuracy) due to
!                                           roundoff in the
!                                           extrapolation procedure
!                                           invoked on this cycle. it
!                                           is assumed that the result
!                                           on this interval is the
!                                           best which can be obtained.
!                                       = 5 the integral over the k th
!                                           cycle is probably divergent
!                                           or slowly convergent. it
!                                           must be noted that
!                                           divergence can occur with
!                                           any other value of
!                                           ierlst(k).
!                         = 8 maximum number of  cycles  allowed
!                             has been achieved, i.e. of subintervals
!                             (a+(k-1)c,a+kc) where
!                             c = (2*int(abs(omega))+1)*pi/abs(omega),
!                             for k = 1, 2, ..., lst.
!                             one can allow more cycles by increasing
!                             the value of limlst (and taking the
!                             according dimension adjustments into
!                             account).
!                             examine the array iwork which contains
!                             the error flags over the cycles, in order
!                             to eventual look for local integration
!                             difficulties.
!                             if the position of a local difficulty can
!                             be determined (e.g. singularity,
!                             discontinuity within the interval)
!                             one will probably gain from splitting
!                             up the interval at this point and
!                             calling appopriate integrators on the
!                             subranges.
!                         = 9 the extrapolation table constructed for
!                             convergence acceleration of the series
!                             formed by the integral contributions
!                             over the cycles, does not converge to
!                             within the required accuracy.
!                             as in the case of ier = 8, it is advised
!                             to examine the array iwork which contains
!                             the error flags on the cycles.
!                    if omega = 0 and integr = 1,
!                    the integral is calculated by means of qagi
!                    and ier = ierlst(1) (with meaning as described
!                    for ierlst(k), k = 1).
!
!            rslst  - real
!                     vector of dimension at least limlst
!                     rslst(k) contains the integral contribution
!                     over the interval (a+(k-1)c,a+kc) where
!                     c = (2*int(abs(omega))+1)*pi/abs(omega),
!                     k = 1, 2, ..., lst.
!                     note that, if omega = 0, rslst(1) contains
!                     the value of the integral over (a,infinity).
!
!            erlst  - real
!                     vector of dimension at least limlst
!                     erlst(k) contains the error estimate
!                     corresponding with rslst(k).
!
!            ierlst - integer ( kind = 8 )
!                     vector of dimension at least limlst
!                     ierlst(k) contains the error flag corresponding
!                     with rslst(k). for the meaning of the local error
!                     flags see description of output parameter ier.
!
!            lst    - integer ( kind = 8 )
!                     number of subintervals needed for the integration
!                     if omega = 0 then lst is set to 1.
!
!            alist, blist, rlist, elist - real
!                     vector of dimension at least limit,
!
!            iord, nnlog - integer ( kind = 8 )
!                     vector of dimension at least limit, providing
!                     space for the quantities needed in the
!                     subdivision process of each cycle
!
!            chebmo - real
!                     array of dimension at least (maxp1,25),
!                     providing space for the Chebyshev moments
!                     needed within the cycles
!
!  Local parameters:
!
!           c1, c2    - end points of subinterval (of length
!                       cycle)
!           cycle     - (2*int(abs(omega))+1)*pi/abs(omega)
!           psum      - vector of dimension at least (limexp+2)
!                       (see routine qextr)
!                       psum contains the part of the epsilon table
!                       which is still needed for further computations.
!                       each element of psum is a partial sum of
!                       the series which should sum to the value of
!                       the integral.
!           errsum    - sum of error estimates over the
!                       subintervals, calculated cumulatively
!           epsa      - absolute tolerance requested over current
!                       subinterval
!           chebmo    - array containing the modified Chebyshev
!                       moments (see also routine qc25o)
!
  implicit none

  integer ( kind = 8 ) limit
  integer ( kind = 8 ) limlst
  integer ( kind = 8 ) maxp1

  real ( kind = 8 ) a
  real ( kind = 8 ) abseps
  real ( kind = 8 ) abserr
  real ( kind = 8 ) alist(limit)
  real ( kind = 8 ) blist(limit)
  real ( kind = 8 ) chebmo(maxp1,25)
  real ( kind = 8 ) correc
  real ( kind = 8 ) cycle
  real ( kind = 8 ) c1
  real ( kind = 8 ) c2
  real ( kind = 8 ) dl
! real ( kind = 8 ) dla
  real ( kind = 8 ) drl
  real ( kind = 8 ) elist(limit)
  real ( kind = 8 ) ep
  real ( kind = 8 ) eps
  real ( kind = 8 ) epsa
  real ( kind = 8 ) epsabs
  real ( kind = 8 ) erlst(limlst)
  real ( kind = 8 ) errsum
  !real ( kind = 8 ), external :: f
  procedure (f_int_dat), intent(in), pointer :: f_ptr 
  real ( kind = 8 ) fact
  integer ( kind = 8 ) ier
  integer ( kind = 8 ) ierlst(limlst)
  integer ( kind = 8 ) integr
  integer ( kind = 8 ) iord(limit)
  integer ( kind = 8 ) ktmin
  integer ( kind = 8 ) l
  integer ( kind = 8 ) ll
  integer ( kind = 8 ) lst
  integer ( kind = 8 ) momcom
  integer ( kind = 8 ) nev
  integer ( kind = 8 ) neval
  integer ( kind = 8 ) nnlog(limit)
  integer ( kind = 8 ) nres
  integer ( kind = 8 ) numrl2
  real ( kind = 8 ) omega
  real ( kind = 8 ), parameter :: p = 0.9E+00
  real ( kind = 8 ), parameter :: pi = 3.1415926535897932E+00
  real ( kind = 8 ) p1
  real ( kind = 8 ) psum(52)
  real ( kind = 8 ) reseps
  real ( kind = 8 ) result
  real ( kind = 8 ) res3la(3)
  real ( kind = 8 ) rlist(limit)
  real ( kind = 8 ) rslst(limlst)
  class( dataCollectionBase), target :: dat
!
!  The dimension of  psum  is determined by the value of
!  limexp in QEXTR (psum must be
!  of dimension (limexp+2) at least).
!
!  Test on validity of parameters.
!
  
  result = 0.0E+00
  abserr = 0.0E+00
  neval = 0
  lst = 0
  ier = 0

  if ( (integr /= 1 .and. integr /= 2 ) .or. &
    epsabs <= 0.0E+00 .or. &
    limlst < 3 ) then
    ier = 6
    return
  end if

  if ( omega == 0.0E+00 ) then

    if ( integr == 1 ) then
      call qagi ( f_ptr, dat, 0.0E+00, int8(1), epsabs, 0.0E+00, result, abserr, neval, ier )
    else
      result = 0.0E+00
      abserr = 0.0E+00
      neval = 0
      ier = 0
    end if

    rslst(1) = result
    erlst(1) = abserr
    ierlst(1) = ier
    lst = 1

    return
  end if
!
!  Initializations.
!
  l = int ( abs ( omega ) )
  dl = 2 * l + 1
  cycle = dl * pi / abs ( omega )
  ier = 0
  ktmin = 0
  neval = 0
  numrl2 = 0
  nres = 0
  c1 = a
  c2 = cycle+a
  p1 = 1.0E+00 - p
  eps = epsabs

  if ( epsabs > tiny ( epsabs ) / p1 ) then
    eps = epsabs * p1
  end if

  ep = eps
  fact = 1.0E+00
  correc = 0.0E+00
  abserr = 0.0E+00
  errsum = 0.0E+00

  do lst = 1, limlst
!
!  Integrate over current subinterval.
!
!   dla = lst
    epsa = eps * fact

    call qfour ( f_ptr, dat, c1, c2, omega, integr, epsa, 0.0E+00, limit, lst, maxp1, &
      rslst(lst), erlst(lst), nev, ierlst(lst), alist, blist, rlist, elist, &
      iord, nnlog, momcom, chebmo )

    neval = neval + nev
    fact = fact * p
    errsum = errsum + erlst(lst)
    drl = 5.0E+01 * abs(rslst(lst))
!
!  Test on accuracy with partial sum.
!
    if ((errsum+drl) <= epsabs.and.lst >= 6) then
      go to 80
    end if

    correc = max ( correc,erlst(lst))

    if ( ierlst(lst) /= 0 ) then
      eps = max ( ep,correc*p1)
      ier = 7
    end if

    if ( ier == 7 .and. (errsum+drl) <= correc*1.0E+01.and. lst > 5) go to 80

    numrl2 = numrl2+1

    if ( lst <= 1 ) then
      psum(1) = rslst(1)
      go to 40
    end if

    psum(numrl2) = psum(ll) + rslst(lst)

    if ( lst == 2 ) then
      go to 40
    end if
!
!  Test on maximum number of subintervals
!
    if ( lst == limlst ) then
      ier = 8
    end if
!
!  Perform new extrapolation
!
    call qextr ( numrl2, psum, reseps, abseps, res3la, nres )
!
!  Test whether extrapolated result is influenced by roundoff
!
    ktmin = ktmin + 1

    if ( ktmin >= 15 .and. abserr <= 1.0E-03 * (errsum+drl) ) then
      ier = 9
    end  if

    if ( abseps <= abserr .or. lst == 3 ) then

      abserr = abseps
      result = reseps
      ktmin = 0
!
!  If IER is not 0, check whether direct result (partial
!  sum) or extrapolated result yields the best integral
!  approximation
!
      if ( ( abserr + 1.0E+01 * correc ) <= epsabs ) then
        exit
      end if

      if ( abserr <= epsabs .and. 1.0E+01 * correc >= epsabs ) then
        exit
      end if

    end if

    if ( ier /= 0 .and. ier /= 7 ) then
      exit
    end if

40  continue

    ll = numrl2
    c1 = c2
    c2 = c2+cycle

  end do
!
!  Set final result and error estimate.
!
!60 continue

  abserr = abserr + 1.0E+01 * correc

  if ( ier == 0 ) then
    return
  end if

  if ( result /= 0.0E+00 .and. psum(numrl2) /= 0.0E+00) go to 70

  if ( abserr > errsum ) then
    go to 80
  end if

  if ( psum(numrl2) == 0.0E+00 ) then
    return
  end if

70 continue

  if ( abserr / abs(result) <= (errsum+drl)/abs(psum(numrl2)) ) then

    if ( ier >= 1 .and. ier /= 7 ) then
      abserr = abserr + drl
    end if

    return

  end if

80 continue

  result = psum(numrl2)
  abserr = errsum + drl

  return
end subroutine
subroutine qawo ( f_ptr, dat, a, b, omega, integr, epsabs, epsrel, result, abserr, &
  neval, ier )

!*****************************************************************************80
!
!! QAWO computes the integrals of oscillatory integrands.
!
!  Discussion:
!
!    The routine calculates an approximation RESULT to a given
!    definite integral
!      I = Integral ( A <= X <= B ) F(X) * cos ( OMEGA * X ) dx
!    or 
!      I = Integral ( A <= X <= B ) F(X) * sin ( OMEGA * X ) dx
!    hopefully satisfying following claim for accuracy
!      | I - RESULT | <= max ( epsabs, epsrel * |I| ).
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real ( kind = 8 ) F, the name of the function routine, of the form
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!    which evaluates the integrand function.
!
!    Input, real ( kind = 8 ) A, B, the limits of integration.
!
!    Input, real ( kind = 8 ) OMEGA, the parameter in the weight function.
!
!    Input, integer ( kind = 8 ) INTEGR, specifies the weight function:
!    1, W(X) = cos ( OMEGA * X )
!    2, W(X) = sin ( OMEGA * X )
!
!    Input, real ( kind = 8 ) EPSABS, EPSREL, the absolute and relative accuracy requested.
!
!    Output, real ( kind = 8 ) RESULT, the estimated value of the integral.
!
!    Output, real ( kind = 8 ) ABSERR, an estimate of || I - RESULT ||.
!
!    Output, integer ( kind = 8 ) NEVAL, the number of times the integral was evaluated.
!
!            ier    - integer ( kind = 8 )
!                     ier = 0 normal and reliable termination of the
!                             routine. it is assumed that the
!                             requested accuracy has been achieved.
!                   - ier > 0 abnormal termination of the routine.
!                             the estimates for integral and error are
!                             less reliable. it is assumed that the
!                             requested accuracy has not been achieved.
!                     ier = 1 maximum number of subdivisions allowed
!                             (= leniw/2) has been achieved. one can
!                             allow more subdivisions by increasing the
!                             value of leniw (and taking the according
!                             dimension adjustments into account).
!                             however, if this yields no improvement it
!                             is advised to analyze the integrand in
!                             order to determine the integration
!                             difficulties. if the position of a local
!                             difficulty can be determined (e.g.
!                             singularity, discontinuity within the
!                             interval) one will probably gain from
!                             splitting up the interval at this point
!                             and calling the integrator on the
!                             subranges. if possible, an appropriate
!                             special-purpose integrator should
!                             be used which is designed for handling
!                             the type of difficulty involved.
!                         = 2 the occurrence of roundoff error is
!                             detected, which prevents the requested
!                             tolerance from being achieved.
!                             the error may be under-estimated.
!                         = 3 extremely bad integrand behavior occurs
!                             at some interior points of the integration
!                             interval.
!                         = 4 the algorithm does not converge. roundoff
!                             error is detected in the extrapolation
!                             table. it is presumed that the requested
!                             tolerance cannot be achieved due to
!                             roundoff in the extrapolation table,
!                             and that the returned result is the best
!                             which can be obtained.
!                         = 5 the integral is probably divergent, or
!                             slowly convergent. it must be noted that
!                             divergence can occur with any other value
!                             of ier.
!                         = 6 the input is invalid, because
!                             epsabs < 0 and epsrel < 0,
!                             result, abserr, neval are set to zero.
!
!  Local parameters:
!
!    limit is the maximum number of subintervals allowed in the
!    subdivision process of QFOUR. take care that limit >= 1.
!
!    maxp1 gives an upper bound on the number of Chebyshev moments
!    which can be stored, i.e. for the intervals of lengths
!    abs(b-a)*2**(-l), l = 0, 1, ... , maxp1-2. take care that
!    maxp1 >= 1.

  implicit none

  integer ( kind = 8 ), parameter :: limit = 500
  integer ( kind = 8 ), parameter :: maxp1 = 21

  real ( kind = 8 ) a
  real ( kind = 8 ) abserr
  real ( kind = 8 ) alist(limit)
  real ( kind = 8 ) b
  real ( kind = 8 ) blist(limit)
  real ( kind = 8 ) chebmo(maxp1,25)
  real ( kind = 8 ) elist(limit)
  real ( kind = 8 ) epsabs
  real ( kind = 8 ) epsrel
  !real ( kind = 8 ), external :: f
  procedure (f_int_dat), intent(in), pointer :: f_ptr 
  integer ( kind = 8 ) ier
  integer ( kind = 8 ) integr
  integer ( kind = 8 ) iord(limit)
  integer ( kind = 8 ) momcom
  integer ( kind = 8 ) neval
  integer ( kind = 8 ) nnlog(limit)
  real ( kind = 8 ) omega
  real ( kind = 8 ) result
  real ( kind = 8 ) rlist(limit)
  class( dataCollectionBase), target :: dat

  call qfour ( f_ptr, dat, a, b, omega, integr, epsabs, epsrel, limit, int8(1), maxp1, &
    result, abserr, neval, ier, alist, blist, rlist, elist, iord, nnlog, &
    momcom, chebmo )

  return
end subroutine
subroutine qaws ( f_ptr, dat, a, b, alfa, beta, integr, epsabs, epsrel, result, &
  abserr, neval, ier )

!*****************************************************************************80
!
!! QAWS estimates integrals with algebraico-logarithmic endpoint singularities.
!
!  Discussion:
!
!    This routine calculates an approximation RESULT to a given
!    definite integral   
!      I = integral of f*w over (a,b) 
!    where w shows a singular behavior at the end points, see parameter
!    integr, hopefully satisfying following claim for accuracy
!      abs(i-result) <= max(epsabs,epsrel*abs(i)).
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real ( kind = 8 ) F, the name of the function routine, of the form
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!    which evaluates the integrand function.
!
!    Input, real ( kind = 8 ) A, B, the limits of integration.
!
!    Input, real ( kind = 8 ) ALFA, BETA, parameters used in the weight function.
!    ALFA and BETA should be greater than -1.
!
!    Input, integer ( kind = 8 ) INTEGR, indicates which weight function is to be used
!    = 1  (x-a)**alfa*(b-x)**beta
!    = 2  (x-a)**alfa*(b-x)**beta*log(x-a)
!    = 3  (x-a)**alfa*(b-x)**beta*log(b-x)
!    = 4  (x-a)**alfa*(b-x)**beta*log(x-a)*log(b-x)
!
!    Input, real ( kind = 8 ) EPSABS, EPSREL, the absolute and relative accuracy requested.
!
!    Output, real ( kind = 8 ) RESULT, the estimated value of the integral.
!
!    Output, real ( kind = 8 ) ABSERR, an estimate of || I - RESULT ||.
!
!    Output, integer ( kind = 8 ) NEVAL, the number of times the integral was evaluated.
!
!            ier    - integer ( kind = 8 )
!                     ier = 0 normal and reliable termination of the
!                             routine. it is assumed that the requested
!                             accuracy has been achieved.
!                     ier > 0 abnormal termination of the routine
!                             the estimates for the integral and error
!                             are less reliable. it is assumed that the
!                             requested accuracy has not been achieved.
!                     ier = 1 maximum number of subdivisions allowed
!                             has been achieved. one can allow more
!                             subdivisions by increasing the data value
!                             of limit in qaws (and taking the according
!                             dimension adjustments into account).
!                             however, if this yields no improvement it
!                             is advised to analyze the integrand, in
!                             order to determine the integration
!                             difficulties which prevent the requested
!                             tolerance from being achieved. in case of
!                             a jump discontinuity or a local
!                             singularity of algebraico-logarithmic type
!                             at one or more interior points of the
!                             integration range, one should proceed by
!                             splitting up the interval at these points
!                             and calling the integrator on the
!                             subranges.
!                         = 2 the occurrence of roundoff error is
!                             detected, which prevents the requested
!                             tolerance from being achieved.
!                         = 3 extremely bad integrand behavior occurs
!                             at some points of the integration
!                             interval.
!                         = 6 the input is invalid, because
!                             b <= a or alfa <= (-1) or beta <= (-1) or
!                             integr < 1 or integr > 4 or
!                             epsabs < 0 and epsrel < 0,
!                             result, abserr, neval are set to zero.
!
!  Local parameters:
!
!    LIMIT is the maximum number of subintervals allowed in the
!    subdivision process of qawse. take care that limit >= 2.
!
  implicit none

  integer ( kind = 8 ), parameter :: limit = 500

  real ( kind = 8 ) a
  real ( kind = 8 ) abserr
  real ( kind = 8 ) alfa
  real ( kind = 8 ) alist(limit)
  real ( kind = 8 ) b
  real ( kind = 8 ) blist(limit)
  real ( kind = 8 ) beta
  real ( kind = 8 ) elist(limit)
  real ( kind = 8 ) epsabs
  real ( kind = 8 ) epsrel
  !real ( kind = 8 ), external :: f
  procedure (f_int_dat), intent(in), pointer :: f_ptr 
  integer ( kind = 8 ) ier
  integer ( kind = 8 ) integr
  integer ( kind = 8 ) iord(limit)
  integer ( kind = 8 ) last
  integer ( kind = 8 ) neval
  real ( kind = 8 ) result
  real ( kind = 8 ) rlist(limit)
  class( dataCollectionBase), target :: dat
  call qawse ( f_ptr, dat, a, b, alfa, beta, integr, epsabs, epsrel, limit, result, &
    abserr, neval, ier, alist, blist, rlist, elist, iord, last )

  return
end subroutine
subroutine qawse ( f_ptr, dat, a, b, alfa, beta, integr, epsabs, epsrel, limit, &
  result, abserr, neval, ier, alist, blist, rlist, elist, iord, last )

!*****************************************************************************80
!
!! QAWSE estimates integrals with algebraico-logarithmic endpoint singularities.
!
!  Discussion:
!
!    This routine calculates an approximation RESULT to an integral
!      I = integral of F(X) * W(X) over (a,b), 
!    where W(X) shows a singular behavior at the endpoints, hopefully 
!    satisfying:
!      | I - RESULT | <= max ( epsabs, epsrel * |I| ).
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real ( kind = 8 ) F, the name of the function routine, of the form
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!    which evaluates the integrand function.
!
!    Input, real ( kind = 8 ) A, B, the limits of integration.
!
!    Input, real ( kind = 8 ) ALFA, BETA, parameters used in the weight function.
!    ALFA and BETA should be greater than -1.
!
!    Input, integer ( kind = 8 ) INTEGR, indicates which weight function is used:
!    = 1  (x-a)**alfa*(b-x)**beta
!    = 2  (x-a)**alfa*(b-x)**beta*log(x-a)
!    = 3  (x-a)**alfa*(b-x)**beta*log(b-x)
!    = 4  (x-a)**alfa*(b-x)**beta*log(x-a)*log(b-x)
!
!    Input, real ( kind = 8 ) EPSABS, EPSREL, the absolute and relative accuracy requested.
!
!    Input, integer ( kind = 8 ) LIMIT, an upper bound on the number of subintervals
!    in the partition of (A,B), LIMIT >= 2.  If LIMIT < 2, the routine 
!     will end with IER = 6.
!
!    Output, real ( kind = 8 ) RESULT, the estimated value of the integral.
!
!    Output, real ( kind = 8 ) ABSERR, an estimate of || I - RESULT ||.
!
!    Output, integer ( kind = 8 ) NEVAL, the number of times the integral was evaluated.
!
!            ier    - integer ( kind = 8 )
!                     ier = 0 normal and reliable termination of the
!                             routine. it is assumed that the requested
!                             accuracy has been achieved.
!                     ier > 0 abnormal termination of the routine
!                             the estimates for the integral and error
!                             are less reliable. it is assumed that the
!                             requested accuracy has not been achieved.
!                         = 1 maximum number of subdivisions allowed
!                             has been achieved. one can allow more
!                             subdivisions by increasing the value of
!                             limit. however, if this yields no
!                             improvement it is advised to analyze the
!                             integrand, in order to determine the
!                             integration difficulties which prevent
!                             the requested tolerance from being
!                             achieved. in case of a jump discontinuity
!                             or a local singularity of algebraico-
!                             logarithmic type at one or more interior
!                             points of the integration range, one
!                             should proceed by splitting up the
!                             interval at these points and calling the
!                             integrator on the subranges.
!                         = 2 the occurrence of roundoff error is
!                             detected, which prevents the requested
!                             tolerance from being achieved.
!                         = 3 extremely bad integrand behavior occurs
!                             at some points of the integration
!                             interval.
!                         = 6 the input is invalid, because
!                             b <= a or alfa <= (-1) or beta <= (-1) or
!                             integr < 1 or integr > 4, or
!                             epsabs < 0 and epsrel < 0,
!                             or limit < 2.
!                             result, abserr, neval, rlist(1), elist(1),
!                             iord(1) and last are set to zero.
!                             alist(1) and blist(1) are set to a and b
!                             respectively.
!
!    Workspace, real ( kind = 8 ) ALIST(LIMIT), BLIST(LIMIT), contains in entries 1 
!    through LAST the left and right ends of the partition subintervals.
!
!    Workspace, real ( kind = 8 ) RLIST(LIMIT), contains in entries 1 through LAST
!    the integral approximations on the subintervals.
!
!    Workspace, real ( kind = 8 ) ELIST(LIMIT), contains in entries 1 through LAST
!    the absolute error estimates on the subintervals.
!
!            iord   - integer ( kind = 8 )
!                     vector of dimension at least limit, the first k
!                     elements of which are pointers to the error
!                     estimates over the subintervals, so that
!                     elist(iord(1)), ..., elist(iord(k)) with k = last
!                     if last <= (limit/2+2), and k = limit+1-last
!                     otherwise, form a decreasing sequence.
!
!    Output, integer ( kind = 8 ) LAST, the number of subintervals actually produced in 
!    the subdivision process.
!
!  Local parameters:
!
!           alist     - list of left end points of all subintervals
!                       considered up to now
!           blist     - list of right end points of all subintervals
!                       considered up to now
!           rlist(i)  - approximation to the integral over
!                       (alist(i),blist(i))
!           elist(i)  - error estimate applying to rlist(i)
!           maxerr    - pointer to the interval with largest error
!                       estimate
!           errmax    - elist(maxerr)
!           area      - sum of the integrals over the subintervals
!           errsum    - sum of the errors over the subintervals
!           errbnd    - requested accuracy max(epsabs,epsrel*
!                       abs(result))
!           *****1    - variable for the left subinterval
!           *****2    - variable for the right subinterval
!           last      - index for subdivision
!
  implicit none

  integer ( kind = 8 ) limit

  real ( kind = 8 ) a
  real ( kind = 8 ) abserr
  real ( kind = 8 ) alfa
  real ( kind = 8 ) alist(limit)
  real ( kind = 8 ) area
  real ( kind = 8 ) area1
  real ( kind = 8 ) area12
  real ( kind = 8 ) area2
  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  real ( kind = 8 ) b
  real ( kind = 8 ) beta
  real ( kind = 8 ) blist(limit)
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) centre
  real ( kind = 8 ) elist(limit)
  real ( kind = 8 ) epsabs
  real ( kind = 8 ) epsrel
  real ( kind = 8 ) errbnd
  real ( kind = 8 ) errmax
  real ( kind = 8 ) error1
  real ( kind = 8 ) erro12
  real ( kind = 8 ) error2
  real ( kind = 8 ) errsum
  !real ( kind = 8 ), external :: f
  procedure (f_int_dat), intent(in), pointer :: f_ptr 
  integer ( kind = 8 ) ier
  integer ( kind = 8 ) integr
  integer ( kind = 8 ) iord(limit)
  integer ( kind = 8 ) iroff1
  integer ( kind = 8 ) iroff2
  integer ( kind = 8 ) last
  integer ( kind = 8 ) maxerr
  integer ( kind = 8 ) nev
  integer ( kind = 8 ) neval
  integer ( kind = 8 ) nrmax
  real ( kind = 8 ) resas1
  real ( kind = 8 ) resas2
  real ( kind = 8 ) result
  real ( kind = 8 ) rg(25)
  real ( kind = 8 ) rh(25)
  real ( kind = 8 ) ri(25)
  real ( kind = 8 ) rj(25)
  real ( kind = 8 ) rlist(limit)
  class( dataCollectionBase), target :: dat
!
!  Test on validity of parameters.
!
  ier = 0
  neval = 0
  last = 0
  rlist(1) = 0.0E+00
  elist(1) = 0.0E+00
  iord(1) = 0
  result = 0.0E+00
  abserr = 0.0E+00

  if ( b <= a .or. &
    (epsabs < 0.0E+00 .and. epsrel < 0.0E+00) .or. &
    alfa <= (-1.0E+00) .or. &
    beta <= (-1.0E+00) .or. &
    integr < 1 .or. &
    integr > 4 .or. &
    limit < 2 ) then
    ier = 6
    return
  end if
!
!  Compute the modified Chebyshev moments.
!
  call qmomo ( alfa, beta, ri, rj, rg, rh, integr )
!
!  Integrate over the intervals (a,(a+b)/2) and ((a+b)/2,b).
!
  centre = 5.0E-01 * ( b + a )

  call qc25s ( f_ptr, dat, a, b, a, centre, alfa, beta, ri, rj, rg, rh, area1, &
    error1, resas1, integr, nev )

  neval = nev

  call qc25s ( f_ptr, dat, a, b, centre, b, alfa, beta, ri, rj, rg, rh, area2, &
    error2, resas2, integr, nev )

  last = 2
  neval = neval+nev
  result = area1+area2
  abserr = error1+error2
!
!  Test on accuracy.
!
  errbnd = max ( epsabs, epsrel * abs ( result ) )
!
!  Initialization.
!
  if ( error2 <= error1 ) then
    alist(1) = a
    alist(2) = centre
    blist(1) = centre
    blist(2) = b
    rlist(1) = area1
    rlist(2) = area2
    elist(1) = error1
    elist(2) = error2
  else
    alist(1) = centre
    alist(2) = a
    blist(1) = b
    blist(2) = centre
    rlist(1) = area2
    rlist(2) = area1
    elist(1) = error2
    elist(2) = error1
  end if

  iord(1) = 1
  iord(2) = 2

  if ( limit == 2 ) then
    ier = 1
    return
  end if

  if ( abserr <= errbnd ) then
    return
  end if

  errmax = elist(1)
  maxerr = 1
  nrmax = 1
  area = result
  errsum = abserr
  iroff1 = 0
  iroff2 = 0

  do last = 3, limit
!
!  Bisect the subinterval with largest error estimate.
!
    a1 = alist(maxerr)
    b1 = 5.0E-01 * ( alist(maxerr) + blist(maxerr) )
    a2 = b1
    b2 = blist(maxerr)

    call qc25s ( f_ptr, dat, a, b, a1, b1, alfa, beta, ri, rj, rg, rh, area1, &
      error1, resas1, integr, nev )

    neval = neval + nev

    call qc25s ( f_ptr, dat, a, b, a2, b2, alfa, beta, ri, rj, rg, rh, area2, &
      error2, resas2, integr, nev )

    neval = neval + nev
!
!  Improve previous approximations integral and error and
!  test for accuracy.
!
    area12 = area1+area2
    erro12 = error1+error2
    errsum = errsum+erro12-errmax
    area = area+area12-rlist(maxerr)
!
!  Test for roundoff error.
!
    if ( a /= a1 .and. b /= b2 ) then

      if ( resas1 /= error1 .and. resas2 /= error2 ) then

        if ( abs ( rlist(maxerr) - area12 ) < 1.0E-05 * abs ( area12 ) &
          .and.erro12 >= 9.9E-01*errmax ) then
          iroff1 = iroff1 + 1
        end if

        if ( last > 10 .and. erro12 > errmax ) then
          iroff2 = iroff2 + 1
        end if

      end if

    end if

    rlist(maxerr) = area1
    rlist(last) = area2
!
!  Test on accuracy.
!
    errbnd = max ( epsabs, epsrel * abs ( area ) )

    if ( errsum > errbnd ) then
!
!  Set error flag in the case that the number of interval
!  bisections exceeds limit.
!
      if ( last == limit ) then
        ier = 1
      end if
!
!  Set error flag in the case of roundoff error.
!
      if ( iroff1 >= 6 .or. iroff2 >= 20 ) then
        ier = 2
     end if
!
!  Set error flag in the case of bad integrand behavior
!  at interior points of integration range.
!
      if ( max ( abs(a1),abs(b2)) <= (1.0E+00+1.0E+03* epsilon ( a1 ) )* &
        ( abs(a2) + 1.0E+03* tiny ( a2) )) then
        ier = 3
      end if

    end if
!
!  Append the newly-created intervals to the list.
!
    if ( error2 <= error1 ) then
      alist(last) = a2
      blist(maxerr) = b1
      blist(last) = b2
      elist(maxerr) = error1
      elist(last) = error2
    else
      alist(maxerr) = a2
      alist(last) = a1
      blist(last) = b1
      rlist(maxerr) = area2
      rlist(last) = area1
      elist(maxerr) = error2
      elist(last) = error1
    end if
!
!  Call QSORT to maintain the descending ordering
!  in the list of error estimates and select the subinterval
!  with largest error estimate (to be bisected next).
!
    call qsort ( limit, last, maxerr, errmax, elist, iord, nrmax )

    if ( ier /= 0 .or. errsum <= errbnd ) then
      exit
    end if

  end do
!
!  Compute final result.
!
  result = sum ( rlist(1:last) )

  abserr = errsum

  return
  end subroutine
  
  function qwgtc ( x, c, p2, p3, p4, kp )

!*****************************************************************************80
!
!! QWGTC defines the weight function used by QC25C.
!
!  Discussion:
!
!    The weight function has the form 1 / ( X - C ).
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which the weight function is evaluated.
!
!    Input, real ( kind = 8 ) C, the location of the singularity.
!
!    Input, real ( kind = 8 ) P2, P3, P4, parameters that are not used.
!
!    Input, integer ( kind = 8 ) KP, a parameter that is not used.
!
!    Output, real ( kind = 8 ) QWGTC, the value of the weight function at X.
!
  implicit none

  real ( kind = 8 ) c
  integer ( kind = 8 ) kp
  real ( kind = 8 ) p2
  real ( kind = 8 ) p3
  real ( kind = 8 ) p4
  real ( kind = 8 ) qwgtc
  real ( kind = 8 ) x

  qwgtc = 1.0E+00 / ( x - c )

  return
end function
  
subroutine qc25c ( f_ptr, dat, a, b, c, result, abserr, krul, neval )

!*****************************************************************************80
!
!! QC25C returns integration rules for Cauchy Principal Value integrals.
!
!  Discussion:
!
!    This routine estimates 
!      I = integral of F(X) * W(X) over (a,b) 
!    with error estimate, where 
!      w(x) = 1/(x-c)
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real ( kind = 8 ) F, the name of the function routine, of the form
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!    which evaluates the integrand function.
!
!    Input, real ( kind = 8 ) A, B, the limits of integration.
!
!    Input, real ( kind = 8 ) C, the parameter in the weight function.
!
!    Output, real ( kind = 8 ) RESULT, the estimated value of the integral.
!    RESULT is computed by using a generalized Clenshaw-Curtis method if
!    C lies within ten percent of the integration interval.  In the 
!    other case the 15-point Kronrod rule obtained by optimal addition
!    of abscissae to the 7-point Gauss rule, is applied.
!
!    Output, real ( kind = 8 ) ABSERR, an estimate of || I - RESULT ||.
!
!           krul   - integer ( kind = 8 )
!                    key which is decreased by 1 if the 15-point
!                    Gauss-Kronrod scheme has been used
!
!    Output, integer ( kind = 8 ) NEVAL, the number of times the integral was evaluated.
!
!  Local parameters:
!
!           fval   - value of the function f at the points
!                    cos(k*pi/24),  k = 0, ..., 24
!           cheb12 - Chebyshev series expansion coefficients, for the
!                    function f, of degree 12
!           cheb24 - Chebyshev series expansion coefficients, for the
!                    function f, of degree 24
!           res12  - approximation to the integral corresponding to the
!                    use of cheb12
!           res24  - approximation to the integral corresponding to the
!                    use of cheb24
!           qwgtc  - external function subprogram defining the weight
!                    function
!           hlgth  - half-length of the interval
!           centr  - mid point of the interval
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) abserr
  real ( kind = 8 ) ak22
  real ( kind = 8 ) amom0
  real ( kind = 8 ) amom1
  real ( kind = 8 ) amom2
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) cc
  real ( kind = 8 ) centr
  real ( kind = 8 ) cheb12(13)
  real ( kind = 8 ) cheb24(25)
  !real ( kind = 8 ), external :: f
  procedure (f_int_dat), intent(in), pointer :: f_ptr 
  real ( kind = 8 ) fval(25)
  real ( kind = 8 ) hlgth
  integer ( kind = 8 ) i
  integer ( kind = 8 ) isym
  integer ( kind = 8 ) k 
  integer ( kind = 8 ) kp
  integer ( kind = 8 ) krul
  integer ( kind = 8 ) neval
  real ( kind = 8 ) p2
  real ( kind = 8 ) p3
  real ( kind = 8 ) p4
  !real ( kind = 8 ), external :: qwgtc
  real ( kind = 8 ) resabs
  real ( kind = 8 ) resasc
  real ( kind = 8 ) result
  real ( kind = 8 ) res12
  real ( kind = 8 ) res24
  real ( kind = 8 ) u
  class( dataCollectionBase), target :: dat
  real ( kind = 8 ), parameter, dimension ( 11 ) :: x = (/ &
    9.914448613738104E-01, 9.659258262890683E-01, &
    9.238795325112868E-01, 8.660254037844386E-01, &
    7.933533402912352E-01, 7.071067811865475E-01, &
    6.087614290087206E-01, 5.000000000000000E-01, &
    3.826834323650898E-01, 2.588190451025208E-01, &
    1.305261922200516E-01 /)
 
!
!  Check the position of C.
!
  cc = ( 2.0E+00 * c - b - a ) / ( b - a )
!
!  Apply the 15-point Gauss-Kronrod scheme.
!

  if ( abs ( cc ) >= 1.1E+00 ) then
    krul = krul - 1
    call qk15w ( f_ptr, dat, qwgtc, c, p2, p3, p4, kp, a, b, result, abserr, &
      resabs, resasc )
    neval = 15
    if ( resasc == abserr ) then
      krul = krul+1
    end if
    return
  end if
!
!  Use the generalized Clenshaw-Curtis method.
!
  hlgth = 5.0E-01 * ( b - a )
  centr = 5.0E-01 * ( b + a )
  neval = 25
  fval(1) = 5.0E-01 * f_ptr(hlgth+centr,dat)
  fval(13) = f_ptr(centr,dat)
  fval(25) = 5.0E-01 * f_ptr(centr-hlgth,dat)

  do i = 2, 12
    u = hlgth * x(i-1)
    isym = 26 - i
    fval(i) = f_ptr(u+centr,dat)
    fval(isym) = f_ptr(centr-u,dat)
  end do
!
!  Compute the Chebyshev series expansion.
!
  call qcheb ( x, fval, cheb12, cheb24 )
!
!  The modified Chebyshev moments are computed by forward
!  recursion, using AMOM0 and AMOM1 as starting values.
!
  amom0 = log ( abs ( ( 1.0E+00 - cc ) / ( 1.0E+00 + cc ) ) )
  amom1 = 2.0E+00 + cc * amom0
  res12 = cheb12(1) * amom0 + cheb12(2) * amom1
  res24 = cheb24(1) * amom0 + cheb24(2) * amom1

  do k = 3, 13
    amom2 = 2.0E+00 * cc * amom1 - amom0
    ak22 = ( k - 2 ) * ( k - 2 )
    if ( ( k / 2 ) * 2 == k ) then
      amom2 = amom2 - 4.0E+00 / ( ak22 - 1.0E+00 )
    end if
    res12 = res12 + cheb12(k) * amom2
    res24 = res24 + cheb24(k) * amom2
    amom0 = amom1
    amom1 = amom2
  end do

  do k = 14, 25
    amom2 = 2.0E+00 * cc * amom1 - amom0
    ak22 = ( k - 2 ) * ( k - 2 )
    if ( ( k / 2 ) * 2 == k ) then
      amom2 = amom2 - 4.0E+00 / ( ak22 - 1.0E+00 )
    end if
    res24 = res24 + cheb24(k) * amom2
    amom0 = amom1
    amom1 = amom2
  end do

  result = res24
  abserr = abs ( res24 - res12 )

  return
end subroutine
subroutine qc25o ( f_ptr, dat, a, b, omega, integr, nrmom, maxp1, ksave, result, &
  abserr, neval, resabs, resasc, momcom, chebmo )

!*****************************************************************************80
!
!! QC25O returns integration rules for integrands with a COS or SIN factor.
!
!  Discussion:
!
!    This routine estimates the integral
!      I = integral of f(x) * w(x) over (a,b)
!    where
!      w(x) = cos(omega*x)
!    or 
!      w(x) = sin(omega*x),
!    and estimates
!      J = integral ( A <= X <= B ) |F(X)| dx.
!
!    For small values of OMEGA or small intervals (a,b) the 15-point
!    Gauss-Kronrod rule is used.  In all other cases a generalized
!    Clenshaw-Curtis method is used, that is, a truncated Chebyshev 
!    expansion of the function F is computed on (a,b), so that the 
!    integrand can be written as a sum of terms of the form W(X)*T(K,X), 
!    where T(K,X) is the Chebyshev polynomial of degree K.  The Chebyshev
!    moments are computed with use of a linear recurrence relation.
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real ( kind = 8 ) F, the name of the function routine, of the form
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!    which evaluates the integrand function.
!
!    Input, real ( kind = 8 ) A, B, the limits of integration.
!
!    Input, real ( kind = 8 ) OMEGA, the parameter in the weight function.
!
!    Input, integer ( kind = 8 ) INTEGR, indicates which weight function is to be used
!    = 1, w(x) = cos(omega*x)
!    = 2, w(x) = sin(omega*x)
!
!    ?, integer ( kind = 8 ) NRMOM, the length of interval (a,b) is equal to the length
!    of the original integration interval divided by
!    2**nrmom (we suppose that the routine is used in an
!    adaptive integration process, otherwise set
!    nrmom = 0).  nrmom must be zero at the first call.
!
!           maxp1  - integer ( kind = 8 )
!                    gives an upper bound on the number of Chebyshev
!                    moments which can be stored, i.e. for the intervals
!                    of lengths abs(bb-aa)*2**(-l), l = 0,1,2, ...,
!                    maxp1-2.
!
!           ksave  - integer ( kind = 8 )
!                    key which is one when the moments for the
!                    current interval have been computed
!
!    Output, real ( kind = 8 ) RESULT, the estimated value of the integral.
!
!           abserr - real
!                    estimate of the modulus of the absolute
!                    error, which should equal or exceed abs(i-result)
!
!    Output, integer ( kind = 8 ) NEVAL, the number of times the integral was evaluated.
!
!    Output, real ( kind = 8 ) RESABS, approximation to the integral J.
!
!    Output, real ( kind = 8 ) RESASC, approximation to the integral of abs(F-I/(B-A)).
!
!         on entry and return
!           momcom - integer ( kind = 8 )
!                    for each interval length we need to compute
!                    the Chebyshev moments. momcom counts the number
!                    of intervals for which these moments have already
!                    been computed. if nrmom < momcom or ksave = 1,
!                    the Chebyshev moments for the interval (a,b)
!                    have already been computed and stored, otherwise
!                    we compute them and we increase momcom.
!
!           chebmo - real
!                    array of dimension at least (maxp1,25) containing
!                    the modified Chebyshev moments for the first momcom
!                    interval lengths
!
!  Local parameters:
!
!    maxp1 gives an upper bound
!           on the number of Chebyshev moments which can be
!           computed, i.e. for the interval (bb-aa), ...,
!           (bb-aa)/2**(maxp1-2).
!           should this number be altered, the first dimension of
!           chebmo needs to be adapted.
!
!    x contains the values cos(k*pi/24)
!           k = 1, ...,11, to be used for the Chebyshev expansion of f
!
!           centr  - mid point of the integration interval
!           hlgth  - half length of the integration interval
!           fval   - value of the function f at the points
!                    (b-a)*0.5*cos(k*pi/12) + (b+a)*0.5
!                    k = 0, ...,24
!           cheb12 - coefficients of the Chebyshev series expansion
!                    of degree 12, for the function f, in the
!                    interval (a,b)
!           cheb24 - coefficients of the Chebyshev series expansion
!                    of degree 24, for the function f, in the
!                    interval (a,b)
!           resc12 - approximation to the integral of
!                    cos(0.5*(b-a)*omega*x)*f(0.5*(b-a)*x+0.5*(b+a))
!                    over (-1,+1), using the Chebyshev series
!                    expansion of degree 12
!           resc24 - approximation to the same integral, using the
!                    Chebyshev series expansion of degree 24
!           ress12 - the analogue of resc12 for the sine
!           ress24 - the analogue of resc24 for the sine
!
  implicit none

  integer ( kind = 8 ) maxp1

  real ( kind = 8 ) a
  real ( kind = 8 ) abserr
  real ( kind = 8 ) ac
  real ( kind = 8 ) an
  real ( kind = 8 ) an2
  real ( kind = 8 ) as
  real ( kind = 8 ) asap
  real ( kind = 8 ) ass
  real ( kind = 8 ) b
  real ( kind = 8 ) centr
  real ( kind = 8 ) chebmo(maxp1,25)
  real ( kind = 8 ) cheb12(13)
  real ( kind = 8 ) cheb24(25)
  real ( kind = 8 ) conc
  real ( kind = 8 ) cons
  real ( kind = 8 ) cospar
  real ( kind = 8 ) d(28)
  real ( kind = 8 ) d1(28)
  real ( kind = 8 ) d2(28)
  real ( kind = 8 ) d3(28)
  real ( kind = 8 ) estc
  real ( kind = 8 ) ests
!  real ( kind = 8 ), external :: f
  procedure (f_int_dat), intent(in), pointer :: f_ptr 
  real ( kind = 8 ) fval(25)
  real ( kind = 8 ) hlgth
  integer ( kind = 8 ) i
  integer ( kind = 8 ) integr
  integer ( kind = 8 ) isym
  integer ( kind = 8 ) j
  integer ( kind = 8 ) k
  integer ( kind = 8 ) ksave
  integer ( kind = 8 ) m
  integer ( kind = 8 ) momcom
  integer ( kind = 8 ) neval
  integer ( kind = 8 ), parameter :: nmac = 28
  integer ( kind = 8 ) noeq1
  integer ( kind = 8 ) noequ
  integer ( kind = 8 ) nrmom
  real ( kind = 8 ) omega
  real ( kind = 8 ) parint
  real ( kind = 8 ) par2
  real ( kind = 8 ) par22
  real ( kind = 8 ) p2
  real ( kind = 8 ) p3
  real ( kind = 8 ) p4
  !real ( kind = 8 ), external :: qwgto
  real ( kind = 8 ) resabs
  real ( kind = 8 ) resasc
  real ( kind = 8 ) resc12
  real ( kind = 8 ) resc24
  real ( kind = 8 ) ress12
  real ( kind = 8 ) ress24
  real ( kind = 8 ) result
  real ( kind = 8 ) sinpar
  real ( kind = 8 ) v(28)
  class( dataCollectionBase), target :: dat
  real ( kind = 8 ), dimension ( 11 ) :: x = (/ &
    9.914448613738104E-01,     9.659258262890683E-01, &
    9.238795325112868E-01,     8.660254037844386E-01, &
    7.933533402912352E-01,     7.071067811865475E-01, &
    6.087614290087206E-01,     5.000000000000000E-01, &
    3.826834323650898E-01,     2.588190451025208E-01, &
    1.305261922200516E-01 /)

  centr = 5.0E-01 * ( b + a )
  hlgth = 5.0E-01 * ( b - a )
  parint = omega * hlgth
  
!
!  Compute the integral using the 15-point Gauss-Kronrod
!  formula if the value of the parameter in the integrand
!  is small or if the length of the integration interval
!  is less than (bb-aa)/2**(maxp1-2), where (aa,bb) is the
!  original integration interval.
!
  if ( abs ( parint ) <= 2.0E+00 ) then

    call qk15w ( f_ptr, dat, qwgto, omega, p2, p3, p4, integr, a, b, result, &
      abserr, resabs, resasc )

    neval = 15
    return

  end if
!
!  Compute the integral using the generalized clenshaw-curtis method.
!
  conc = hlgth * cos(centr*omega)
  cons = hlgth * sin(centr*omega)
  resasc = huge ( resasc )
  neval = 25
!
!  Check whether the Chebyshev moments for this interval
!  have already been computed.
!
  if ( nrmom < momcom .or. ksave == 1 ) then
    go to 140
  end if
!
!  Compute a new set of Chebyshev moments.
!
  m = momcom + 1
  par2 = parint * parint
  par22 = par2 + 2.0E+00
  sinpar = sin(parint)
  cospar = cos(parint)
!
!  Compute the Chebyshev moments with respect to cosine.
!
  v(1) = 2.0E+00 * sinpar / parint
  v(2) = (8.0E+00*cospar+(par2+par2-8.0E+00)*sinpar/ parint)/par2
  v(3) = (3.2E+01*(par2-1.2E+01)*cospar+(2.0E+00* &
     ((par2-8.0E+01)*par2+1.92E+02)*sinpar)/ &
     parint)/(par2*par2)
  ac = 8.0E+00*cospar
  as = 2.4E+01*parint*sinpar

  if ( abs ( parint ) > 2.4E+01 ) then
    go to 70
  end if
!
!  Compute the Chebyshev moments as the solutions of a boundary value 
!  problem with one initial value (v(3)) and one end value computed
!  using an asymptotic formula.
!
  noequ = nmac-3
  noeq1 = noequ-1
  an = 6.0E+00

  do k = 1, noeq1
    an2 = an*an
    d(k) = -2.0E+00*(an2-4.0E+00) * (par22-an2-an2)
    d2(k) = (an-1.0E+00)*(an-2.0E+00) * par2
    d1(k) = (an+3.0E+00)*(an+4.0E+00) * par2
    v(k+3) = as-(an2-4.0E+00) * ac
    an = an+2.0E+00
  end do

  an2 = an*an
  d(noequ) = -2.0E+00*(an2-4.0E+00) * (par22-an2-an2)
  v(noequ+3) = as - ( an2 - 4.0E+00 ) * ac
  v(4) = v(4) - 5.6E+01 * par2 * v(3)
  ass = parint * sinpar
  asap = (((((2.10E+02*par2-1.0E+00)*cospar-(1.05E+02*par2 &
     -6.3E+01)*ass)/an2-(1.0E+00-1.5E+01*par2)*cospar &
     +1.5E+01*ass)/an2-cospar+3.0E+00*ass)/an2-cospar)/an2
  v(noequ+3) = v(noequ+3)-2.0E+00*asap*par2*(an-1.0E+00)* &
      (an-2.0E+00)
!
!  Solve the tridiagonal system by means of Gaussian
!  elimination with partial pivoting.
!
  d3(1:noequ) = 0.0E+00

  d2(noequ) = 0.0E+00

  do i = 1, noeq1

    if ( abs(d1(i)) > abs(d(i)) ) then
      an = d1(i)
      d1(i) = d(i)
      d(i) = an
      an = d2(i)
      d2(i) = d(i+1)
      d(i+1) = an
      d3(i) = d2(i+1)
      d2(i+1) = 0.0E+00
      an = v(i+4)
      v(i+4) = v(i+3)
      v(i+3) = an
    end if

    d(i+1) = d(i+1)-d2(i)*d1(i)/d(i)
    d2(i+1) = d2(i+1)-d3(i)*d1(i)/d(i)
    v(i+4) = v(i+4)-v(i+3)*d1(i)/d(i)

  end do

  v(noequ+3) = v(noequ+3) / d(noequ)
  v(noequ+2) = (v(noequ+2)-d2(noeq1)*v(noequ+3))/d(noeq1)

  do i = 2, noeq1
    k = noequ-i
    v(k+3) = (v(k+3)-d3(k)*v(k+5)-d2(k)*v(k+4))/d(k)
  end do

  go to 90
!
!  Compute the Chebyshev moments by means of forward recursion
!
70 continue

  an = 4.0E+00

  do i = 4, 13
    an2 = an*an
    v(i) = ((an2-4.0E+00)*(2.0E+00*(par22-an2-an2)*v(i-1)-ac) &
     +as-par2*(an+1.0E+00)*(an+2.0E+00)*v(i-2))/ &
     (par2*(an-1.0E+00)*(an-2.0E+00))
    an = an+2.0E+00
  end do

90 continue

  do j = 1, 13
    chebmo(m,2*j-1) = v(j)
  end do
!
!  Compute the Chebyshev moments with respect to sine.
!
  v(1) = 2.0E+00*(sinpar-parint*cospar)/par2
  v(2) = (1.8E+01-4.8E+01/par2)*sinpar/par2 &
     +(-2.0E+00+4.8E+01/par2)*cospar/parint
  ac = -2.4E+01*parint*cospar
  as = -8.0E+00*sinpar
  chebmo(m,2) = v(1)
  chebmo(m,4) = v(2)

  if ( abs(parint) <= 2.4E+01 ) then

    do k = 3, 12
      an = k
      chebmo(m,2*k) = -sinpar/(an*(2.0E+00*an-2.0E+00)) &
                 -2.5E-01*parint*(v(k+1)/an-v(k)/(an-1.0E+00))
    end do
!
!  Compute the Chebyshev moments by means of forward recursion.
!
  else

    an = 3.0E+00

    do i = 3, 12
      an2 = an*an
      v(i) = ((an2-4.0E+00)*(2.0E+00*(par22-an2-an2)*v(i-1)+as) &
       +ac-par2*(an+1.0E+00)*(an+2.0E+00)*v(i-2)) &
       /(par2*(an-1.0E+00)*(an-2.0E+00))
      an = an+2.0E+00
      chebmo(m,2*i) = v(i)
    end do

  end if

140 continue

  if ( nrmom < momcom ) then
    m = nrmom + 1
  end if

  if ( momcom < maxp1 - 1 .and. nrmom >= momcom ) then
    momcom = momcom + 1
  end if
!
!  Compute the coefficients of the Chebyshev expansions
!  of degrees 12 and 24 of the function F.
!
  fval(1) = 5.0E-01 * f_ptr(centr+hlgth,dat)
  fval(13) = f_ptr(centr,dat)
  fval(25) = 5.0E-01 * f_ptr(centr-hlgth,dat)

  do i = 2, 12
    isym = 26-i
    fval(i) = f_ptr(hlgth*x(i-1)+centr,dat)
    fval(isym) = f_ptr(centr-hlgth*x(i-1),dat)
  end do

  call qcheb ( x, fval, cheb12, cheb24 )
!
!  Compute the integral and error estimates.
!
  resc12 = cheb12(13) * chebmo(m,13)
  ress12 = 0.0E+00
  estc = abs ( cheb24(25)*chebmo(m,25))+abs((cheb12(13)- &
    cheb24(13))*chebmo(m,13) )
  ests = 0.0E+00
  k = 11

  do j = 1, 6
    resc12 = resc12+cheb12(k)*chebmo(m,k)
    ress12 = ress12+cheb12(k+1)*chebmo(m,k+1)
    estc = estc+abs((cheb12(k)-cheb24(k))*chebmo(m,k))
    ests = ests+abs((cheb12(k+1)-cheb24(k+1))*chebmo(m,k+1))
    k = k-2
  end do

  resc24 = cheb24(25)*chebmo(m,25)
  ress24 = 0.0E+00
  resabs = abs(cheb24(25))
  k = 23

  do j = 1, 12

    resc24 = resc24+cheb24(k)*chebmo(m,k)
    ress24 = ress24+cheb24(k+1)*chebmo(m,k+1)
    resabs = resabs+abs(cheb24(k))+abs(cheb24(k+1))

    if ( j <= 5 ) then
      estc = estc+abs(cheb24(k)*chebmo(m,k))
      ests = ests+abs(cheb24(k+1)*chebmo(m,k+1))
    end if

    k = k-2

  end do

  resabs = resabs * abs ( hlgth )

  if ( integr == 1 ) then
    result = conc * resc24-cons*ress24
    abserr = abs ( conc * estc ) + abs ( cons * ests )
  else
    result = conc*ress24+cons*resc24
    abserr = abs(conc*ests)+abs(cons*estc)
  end if

  return
end subroutine
subroutine qc25s ( f_ptr, dat, a, b, bl, br, alfa, beta, ri, rj, rg, rh, result, &
  abserr, resasc, integr, neval )

!*****************************************************************************80
!
!! QC25S returns rules for algebraico-logarithmic end point singularities.
!
!  Discussion:
!
!    This routine computes 
!      i = integral of F(X) * W(X) over (bl,br), 
!    with error estimate, where the weight function W(X) has a singular
!    behavior of algebraico-logarithmic type at the points
!    a and/or b. 
!
!    The interval (bl,br) is a subinterval of (a,b).
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real ( kind = 8 ) F, the name of the function routine, of the form
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!    which evaluates the integrand function.
!
!    Input, real ( kind = 8 ) A, B, the limits of integration.
!
!    Input, real ( kind = 8 ) BL, BR, the lower and upper limits of integration.
!    A <= BL < BR <= B.
!
!    Input, real ( kind = 8 ) ALFA, BETA, parameters in the weight function.
!
!    Input, real ( kind = 8 ) RI(25), RJ(25), RG(25), RH(25), modified Chebyshev moments 
!    for the application of the generalized Clenshaw-Curtis method,
!    computed in QMOMO.
!
!    Output, real ( kind = 8 ) RESULT, the estimated value of the integral, computed by 
!    using a generalized clenshaw-curtis method if b1 = a or br = b.
!    In all other cases the 15-point Kronrod rule is applied, obtained by
!    optimal addition of abscissae to the 7-point Gauss rule.
!
!    Output, real ( kind = 8 ) ABSERR, an estimate of || I - RESULT ||.
!
!    Output, real ( kind = 8 ) RESASC, approximation to the integral of abs(F*W-I/(B-A)).
!
!    Input, integer ( kind = 8 ) INTEGR,  determines the weight function
!    1, w(x) = (x-a)**alfa*(b-x)**beta
!    2, w(x) = (x-a)**alfa*(b-x)**beta*log(x-a)
!    3, w(x) = (x-a)**alfa*(b-x)**beta*log(b-x)
!    4, w(x) = (x-a)**alfa*(b-x)**beta*log(x-a)*log(b-x)
!
!    Output, integer ( kind = 8 ) NEVAL, the number of times the integral was evaluated.
!
!  Local Parameters:
!
!           fval   - value of the function f at the points
!                    (br-bl)*0.5*cos(k*pi/24)+(br+bl)*0.5
!                    k = 0, ..., 24
!           cheb12 - coefficients of the Chebyshev series expansion
!                    of degree 12, for the function f, in the interval
!                    (bl,br)
!           cheb24 - coefficients of the Chebyshev series expansion
!                    of degree 24, for the function f, in the interval
!                    (bl,br)
!           res12  - approximation to the integral obtained from cheb12
!           res24  - approximation to the integral obtained from cheb24
!           qwgts  - external function subprogram defining the four
!                    possible weight functions
!           hlgth  - half-length of the interval (bl,br)
!           centr  - mid point of the interval (bl,br)
!
!           the vector x contains the values cos(k*pi/24)
!           k = 1, ..., 11, to be used for the computation of the
!           Chebyshev series expansion of f.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) abserr
  real ( kind = 8 ) alfa
  real ( kind = 8 ) b
  real ( kind = 8 ) beta
  real ( kind = 8 ) bl
  real ( kind = 8 ) br
  real ( kind = 8 ) centr
  real ( kind = 8 ) cheb12(13)
  real ( kind = 8 ) cheb24(25)
  real ( kind = 8 ) dc
  !real ( kind = 8 ), external :: f
  procedure (f_int_dat), intent(in), pointer :: f_ptr 
  real ( kind = 8 ) factor
  real ( kind = 8 ) fix
  real ( kind = 8 ) fval(25)
  real ( kind = 8 ) hlgth
  integer ( kind = 8 ) i
  integer ( kind = 8 ) integr
  integer ( kind = 8 ) isym
  integer ( kind = 8 ) neval
  !real ( kind = 8 ), external :: qwgts
  real ( kind = 8 ) resabs
  real ( kind = 8 ) resasc
  real ( kind = 8 ) result
  real ( kind = 8 ) res12
  real ( kind = 8 ) res24
  real ( kind = 8 ) rg(25)
  real ( kind = 8 ) rh(25)
  real ( kind = 8 ) ri(25)
  real ( kind = 8 ) rj(25)
  real ( kind = 8 ) u
  class( dataCollectionBase), target :: dat
  real ( kind = 8 ), dimension ( 11 ) :: x = (/ &
    9.914448613738104E-01,     9.659258262890683E-01, &
    9.238795325112868E-01,     8.660254037844386E-01, &
    7.933533402912352E-01,     7.071067811865475E-01, &
    6.087614290087206E-01,     5.000000000000000E-01, &
    3.826834323650898E-01,     2.588190451025208E-01, &
    1.305261922200516E-01 /)

  neval = 25
  
  

  if ( bl == a .and. (alfa /= 0.0E+00 .or. integr == 2 .or. integr == 4)) then
    go to 10
  end if

  if ( br == b .and. (beta /= 0.0E+00 .or. integr == 3 .or. integr == 4)) &
    go to 140
!
!  If a > bl and b < br, apply the 15-point Gauss-Kronrod scheme.
!
  call qk15w ( f_ptr, dat, qwgts, a, b, alfa, beta, integr, bl, br, result, abserr, &
    resabs, resasc )

  neval = 15
  return
!
!  This part of the program is executed only if a = bl.
!
!  Compute the Chebyshev series expansion of the function
!  f1 = (0.5*(b+b-br-a)-0.5*(br-a)*x)**beta*f(0.5*(br-a)*x+0.5*(br+a))
!
10 continue

  hlgth = 5.0E-01*(br-bl)
  centr = 5.0E-01*(br+bl)
  fix = b-centr
  fval(1) = 5.0E-01*f_ptr(hlgth+centr,dat)*(fix-hlgth)**beta
  fval(13) = f_ptr(centr,dat)*(fix**beta)
  fval(25) = 5.0E-01*f_ptr(centr-hlgth,dat)*(fix+hlgth)**beta

  do i = 2, 12
    u = hlgth*x(i-1)
    isym = 26-i
    fval(i) = f_ptr(u+centr,dat)*(fix-u)**beta
    fval(isym) = f_ptr(centr-u,dat)*(fix+u)**beta
  end do

  factor = hlgth**(alfa+1.0E+00)
  result = 0.0E+00
  abserr = 0.0E+00
  res12 = 0.0E+00
  res24 = 0.0E+00

  if ( integr > 2 ) go to 70

  call qcheb ( x, fval, cheb12, cheb24 )
!
!  integr = 1  (or 2)
!
  do i = 1, 13
    res12 = res12+cheb12(i)*ri(i)
    res24 = res24+cheb24(i)*ri(i)
  end do

  do i = 14, 25
    res24 = res24 + cheb24(i) * ri(i)
  end do

  if ( integr == 1 ) go to 130
!
!  integr = 2
!
  dc = log ( br - bl )
  result = res24 * dc
  abserr = abs((res24-res12)*dc)
  res12 = 0.0E+00
  res24 = 0.0E+00

  do i = 1, 13
    res12 = res12+cheb12(i)*rg(i)
    res24 = res24+cheb24(i)*rg(i)
  end do

  do i = 14, 25
    res24 = res24+cheb24(i)*rg(i)
  end do

  go to 130
!
!  Compute the Chebyshev series expansion of the function
!  F4 = f1*log(0.5*(b+b-br-a)-0.5*(br-a)*x)
!
70 continue

  fval(1) = fval(1) * log ( fix - hlgth )
  fval(13) = fval(13) * log ( fix )
  fval(25) = fval(25) * log ( fix + hlgth )

  do i = 2, 12
    u = hlgth*x(i-1)
    isym = 26-i
    fval(i) = fval(i) * log ( fix - u )
    fval(isym) = fval(isym) * log ( fix + u )
  end do

  call qcheb ( x, fval, cheb12, cheb24 )
!
!  integr = 3  (or 4)
!
  do i = 1, 13
    res12 = res12+cheb12(i)*ri(i)
    res24 = res24+cheb24(i)*ri(i)
  end do

  do i = 14, 25
    res24 = res24+cheb24(i)*ri(i)
  end do

  if ( integr == 3 ) then
    go to 130
  end if
!
!  integr = 4
!
  dc = log ( br - bl )
  result = res24*dc
  abserr = abs((res24-res12)*dc)
  res12 = 0.0E+00
  res24 = 0.0E+00

  do i = 1, 13
    res12 = res12+cheb12(i)*rg(i)
    res24 = res24+cheb24(i)*rg(i)
  end do

  do i = 14, 25
    res24 = res24+cheb24(i)*rg(i)
  end do

130 continue

  result = (result+res24)*factor
  abserr = (abserr+abs(res24-res12))*factor
  go to 270
!
!  This part of the program is executed only if b = br.
!
!  Compute the Chebyshev series expansion of the function
!  f2 = (0.5*(b+bl-a-a)+0.5*(b-bl)*x)**alfa*f(0.5*(b-bl)*x+0.5*(b+bl))
!
140 continue

  hlgth = 5.0E-01*(br-bl)
  centr = 5.0E-01*(br+bl)
  fix = centr-a
  fval(1) = 5.0E-01*f_ptr(hlgth+centr,dat)*(fix+hlgth)**alfa
  fval(13) = f_ptr(centr,dat)*(fix**alfa)
  fval(25) = 5.0E-01*f_ptr(centr-hlgth,dat)*(fix-hlgth)**alfa

  do i = 2, 12
    u = hlgth*x(i-1)
    isym = 26-i
    fval(i) = f_ptr(u+centr,dat)*(fix+u)**alfa
    fval(isym) = f_ptr(centr-u,dat)*(fix-u)**alfa
  end do

  factor = hlgth**(beta+1.0E+00)
  result = 0.0E+00
  abserr = 0.0E+00
  res12 = 0.0E+00
  res24 = 0.0E+00

  if ( integr == 2 .or. integr == 4 ) then
    go to 200
  end if
!
!  integr = 1  (or 3)
!
  call qcheb ( x, fval, cheb12, cheb24 )

  do i = 1, 13
    res12 = res12+cheb12(i)*rj(i)
    res24 = res24+cheb24(i)*rj(i)
  end do

  do i = 14, 25
    res24 = res24+cheb24(i)*rj(i)
  end do

  if ( integr == 1 ) go to 260
!
!  integr = 3
!
  dc = log ( br - bl )
  result = res24*dc
  abserr = abs((res24-res12)*dc)
  res12 = 0.0E+00
  res24 = 0.0E+00

  do i = 1, 13
    res12 = res12+cheb12(i)*rh(i)
    res24 = res24+cheb24(i)*rh(i)
  end do

  do i = 14, 25
    res24 = res24+cheb24(i)*rh(i)
  end do

  go to 260
!
!  Compute the Chebyshev series expansion of the function
!  f3 = f2*log(0.5*(b-bl)*x+0.5*(b+bl-a-a))
!
200 continue

  fval(1) = fval(1) * log ( hlgth + fix )
  fval(13) = fval(13) * log ( fix )
  fval(25) = fval(25) * log ( fix - hlgth )

  do i = 2, 12
    u = hlgth*x(i-1)
    isym = 26-i
    fval(i) = fval(i) * log(u+fix)
    fval(isym) = fval(isym) * log(fix-u)
  end do

  call qcheb ( x, fval, cheb12, cheb24 )
!
!  integr = 2  (or 4)
!
  do i = 1, 13
    res12 = res12+cheb12(i)*rj(i)
    res24 = res24+cheb24(i)*rj(i)
  end do

  do i = 14, 25
    res24 = res24+cheb24(i)*rj(i)
  end do

  if ( integr == 2 ) go to 260

  dc = log(br-bl)
  result = res24*dc
  abserr = abs((res24-res12)*dc)
  res12 = 0.0E+00
  res24 = 0.0E+00
!
!  integr = 4
!
  do i = 1, 13
    res12 = res12+cheb12(i)*rh(i)
    res24 = res24+cheb24(i)*rh(i)
  end do

  do i = 14, 25
    res24 = res24+cheb24(i)*rh(i)
  end do

260 continue

  result = (result+res24)*factor
  abserr = (abserr+abs(res24-res12))*factor

270 continue

  return
end subroutine
subroutine qcheb ( x, fval, cheb12, cheb24 )

!*****************************************************************************80
!
!! QCHEB computes the Chebyshev series expansion.
!
!  Discussion:
!
!    This routine computes the Chebyshev series expansion
!    of degrees 12 and 24 of a function using a fast Fourier transform method
!
!      f(x) = sum(k=1, ...,13) (cheb12(k)*t(k-1,x)),
!      f(x) = sum(k=1, ...,25) (cheb24(k)*t(k-1,x)),
!
!    where T(K,X) is the Chebyshev polynomial of degree K.
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(11), contains the values of COS(K*PI/24), for K = 1 to 11.
!
!    Input/output, real ( kind = 8 ) FVAL(25), the function values at the points
!    (b+a+(b-a)*cos(k*pi/24))/2, k = 0, ...,24, where (a,b) is the 
!    approximation interval.  FVAL(1) and FVAL(25) are divided by two
!    These values are destroyed at output.
!
!    Output, real ( kind = 8 ) CHEB12(13), the Chebyshev coefficients for degree 12.
!
!    Output, real ( kind = 8 ) CHEB24(25), the Chebyshev coefficients for degree 24.
!
  implicit none

  real ( kind = 8 ) alam
  real ( kind = 8 ) alam1
  real ( kind = 8 ) alam2
  real ( kind = 8 ) cheb12(13)
  real ( kind = 8 ) cheb24(25)
  real ( kind = 8 ) fval(25)
  integer ( kind = 8 ) i
  integer ( kind = 8 ) j
  real ( kind = 8 ) part1
  real ( kind = 8 ) part2
  real ( kind = 8 ) part3
  real ( kind = 8 ) v(12)
  real ( kind = 8 ) x(11)

  do i = 1, 12
    j = 26-i
    v(i) = fval(i)-fval(j)
    fval(i) = fval(i)+fval(j)
  end do

  alam1 = v(1)-v(9)
  alam2 = x(6)*(v(3)-v(7)-v(11))
  cheb12(4) = alam1+alam2
  cheb12(10) = alam1-alam2
  alam1 = v(2)-v(8)-v(10)
  alam2 = v(4)-v(6)-v(12)
  alam = x(3)*alam1+x(9)*alam2
  cheb24(4) = cheb12(4)+alam
  cheb24(22) = cheb12(4)-alam
  alam = x(9)*alam1-x(3)*alam2
  cheb24(10) = cheb12(10)+alam
  cheb24(16) = cheb12(10)-alam
  part1 = x(4)*v(5)
  part2 = x(8)*v(9)
  part3 = x(6)*v(7)
  alam1 = v(1)+part1+part2
  alam2 = x(2)*v(3)+part3+x(10)*v(11)
  cheb12(2) = alam1+alam2
  cheb12(12) = alam1-alam2
  alam = x(1)*v(2)+x(3)*v(4)+x(5)*v(6)+x(7)*v(8) &
    +x(9)*v(10)+x(11)*v(12)
  cheb24(2) = cheb12(2)+alam
  cheb24(24) = cheb12(2)-alam
  alam = x(11)*v(2)-x(9)*v(4)+x(7)*v(6)-x(5)*v(8) &
    +x(3)*v(10)-x(1)*v(12)
  cheb24(12) = cheb12(12)+alam
  cheb24(14) = cheb12(12)-alam
  alam1 = v(1)-part1+part2
  alam2 = x(10)*v(3)-part3+x(2)*v(11)
  cheb12(6) = alam1+alam2
  cheb12(8) = alam1-alam2
  alam = x(5)*v(2)-x(9)*v(4)-x(1)*v(6) &
    -x(11)*v(8)+x(3)*v(10)+x(7)*v(12)
  cheb24(6) = cheb12(6)+alam
  cheb24(20) = cheb12(6)-alam
  alam = x(7)*v(2)-x(3)*v(4)-x(11)*v(6)+x(1)*v(8) &
    -x(9)*v(10)-x(5)*v(12)
  cheb24(8) = cheb12(8)+alam
  cheb24(18) = cheb12(8)-alam

  do i = 1, 6
    j = 14-i
    v(i) = fval(i)-fval(j)
    fval(i) = fval(i)+fval(j)
  end do

  alam1 = v(1)+x(8)*v(5)
  alam2 = x(4)*v(3)
  cheb12(3) = alam1+alam2
  cheb12(11) = alam1-alam2
  cheb12(7) = v(1)-v(5)
  alam = x(2)*v(2)+x(6)*v(4)+x(10)*v(6)
  cheb24(3) = cheb12(3)+alam
  cheb24(23) = cheb12(3)-alam
  alam = x(6)*(v(2)-v(4)-v(6))
  cheb24(7) = cheb12(7)+alam
  cheb24(19) = cheb12(7)-alam
  alam = x(10)*v(2)-x(6)*v(4)+x(2)*v(6)
  cheb24(11) = cheb12(11)+alam
  cheb24(15) = cheb12(11)-alam

  do i = 1, 3
    j = 8-i
    v(i) = fval(i)-fval(j)
    fval(i) = fval(i)+fval(j)
  end do

  cheb12(5) = v(1)+x(8)*v(3)
  cheb12(9) = fval(1)-x(8)*fval(3)
  alam = x(4)*v(2)
  cheb24(5) = cheb12(5)+alam
  cheb24(21) = cheb12(5)-alam
  alam = x(8)*fval(2)-fval(4)
  cheb24(9) = cheb12(9)+alam
  cheb24(17) = cheb12(9)-alam
  cheb12(1) = fval(1)+fval(3)
  alam = fval(2)+fval(4)
  cheb24(1) = cheb12(1)+alam
  cheb24(25) = cheb12(1)-alam
  cheb12(13) = v(1)-v(3)
  cheb24(13) = cheb12(13)
  alam = 1.0E+00/6.0E+00

  do i = 2, 12
    cheb12(i) = cheb12(i)*alam
  end do

  alam = 5.0E-01*alam
  cheb12(1) = cheb12(1)*alam
  cheb12(13) = cheb12(13)*alam

  do i = 2, 24
    cheb24(i) = cheb24(i)*alam
  end do

  cheb24(1) = 0.5E+00 * alam*cheb24(1)
  cheb24(25) = 0.5E+00 * alam*cheb24(25)

  return
end subroutine
subroutine qextr ( n, epstab, result, abserr, res3la, nres )

!*****************************************************************************80
!
!! QEXTR carries out the Epsilon extrapolation algorithm.
!
!  Discussion:
!
!    The routine determines the limit of a given sequence of approximations, 
!    by means of the epsilon algorithm of P. Wynn.  An estimate of the 
!    absolute error is also given.  The condensed epsilon table is computed.
!    Only those elements needed for the computation of the next diagonal
!    are preserved.
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, integer ( kind = 8 ) N, indicates the entry of EPSTAB which contains
!    the new element in the first column of the epsilon table.
!
!    Input/output, real ( kind = 8 ) EPSTAB(52), the two lower diagonals of the triangular
!    epsilon table.  The elements are numbered starting at the right-hand 
!    corner of the triangle.
!
!    Output, real ( kind = 8 ) RESULT, the estimated value of the integral.
!
!    Output, real ( kind = 8 ) ABSERR, estimate of the absolute error computed from
!    RESULT and the 3 previous results.
!
!    ?, real ( kind = 8 ) RES3LA(3), the last 3 results.
!
!    Input/output, integer ( kind = 8 ) NRES, the number of calls to the routine.  This
!    should be zero on the first call, and is automatically updated
!    before return.
!
!  Local Parameters:
!
!           e0     - the 4 elements on which the
!           e1       computation of a new element in
!           e2       the epsilon table is based
!           e3                 e0
!                        e3    e1    new
!                              e2
!           newelm - number of elements to be computed in the new
!                    diagonal
!           error  - error = abs(e1-e0)+abs(e2-e1)+abs(new-e2)
!           result - the element in the new diagonal with least value
!                    of error
!           limexp is the maximum number of elements the epsilon table
!           can contain. if this number is reached, the upper diagonal
!           of the epsilon table is deleted.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) delta1
  real ( kind = 8 ) delta2
  real ( kind = 8 ) delta3
  real ( kind = 8 ) epsinf
  real ( kind = 8 ) epstab(52)
  real ( kind = 8 ) error
  real ( kind = 8 ) err1
  real ( kind = 8 ) err2
  real ( kind = 8 ) err3
  real ( kind = 8 ) e0
  real ( kind = 8 ) e1
  real ( kind = 8 ) e1abs
  real ( kind = 8 ) e2
  real ( kind = 8 ) e3
  integer ( kind = 8 ) i
  integer ( kind = 8 ) ib
  integer ( kind = 8 ) ib2
  integer ( kind = 8 ) ie
  integer ( kind = 8 ) indx
  integer ( kind = 8 ) k1
  integer ( kind = 8 ) k2
  integer ( kind = 8 ) k3
  integer ( kind = 8 ) limexp
  integer ( kind = 8 ) n
  integer ( kind = 8 ) newelm
  integer ( kind = 8 ) nres
  integer ( kind = 8 ) num
  real ( kind = 8 ) res
  real ( kind = 8 ) result
  real ( kind = 8 ) res3la(3)
  real ( kind = 8 ) ss
  real ( kind = 8 ) tol1
  real ( kind = 8 ) tol2
  real ( kind = 8 ) tol3

  nres = nres+1
  abserr = huge ( abserr )
  result = epstab(n)

  if ( n < 3 ) then
    abserr = max ( abserr,0.5E+00* epsilon ( result ) *abs(result))
    return
  end if

  limexp = 50
  epstab(n+2) = epstab(n)
  newelm = (n-1)/2
  epstab(n) = huge ( epstab(n) )
  num = n
  k1 = n

  do i = 1, newelm

    k2 = k1-1
    k3 = k1-2
    res = epstab(k1+2)
    e0 = epstab(k3)
    e1 = epstab(k2)
    e2 = res
    e1abs = abs(e1)
    delta2 = e2-e1
    err2 = abs(delta2)
    tol2 = max ( abs(e2),e1abs)* epsilon ( e2 )
    delta3 = e1-e0
    err3 = abs(delta3)
    tol3 = max ( e1abs,abs(e0))* epsilon ( e0 )
!
!  If e0, e1 and e2 are equal to within machine accuracy, convergence 
!  is assumed.
!
    if ( err2 <= tol2 .and. err3 <= tol3 ) then
      result = res
      abserr = err2+err3
      abserr = max ( abserr,0.5E+00* epsilon ( result ) *abs(result))
      return
    end if

    e3 = epstab(k1)
    epstab(k1) = e1
    delta1 = e1-e3
    err1 = abs(delta1)
    tol1 = max ( e1abs,abs(e3))* epsilon ( e3 )
!
!  If two elements are very close to each other, omit a part
!  of the table by adjusting the value of N.
!
    if ( err1 <= tol1 .or. err2 <= tol2 .or. err3 <= tol3 ) go to 20

    ss = 1.0E+00/delta1+1.0E+00/delta2-1.0E+00/delta3
    epsinf = abs ( ss*e1 )
!
!  Test to detect irregular behavior in the table, and
!  eventually omit a part of the table adjusting the value of N.
!
    if ( epsinf > 1.0E-04 ) go to 30

20  continue

    n = i+i-1
    exit
!
!  Compute a new element and eventually adjust the value of RESULT.
!
30  continue

    res = e1+1.0E+00/ss
    epstab(k1) = res
    k1 = k1-2
    error = err2+abs(res-e2)+err3

    if ( error <= abserr ) then
      abserr = error
      result = res
    end if

  end do
!
!  Shift the table.
!
  if ( n == limexp ) then
    n = 2*(limexp/2)-1
  end if

  if ( (num/2)*2 == num ) then
    ib = 2
  else
    ib = 1
  end if

  ie = newelm+1

  do i = 1, ie
    ib2 = ib+2
    epstab(ib) = epstab(ib2)
    ib = ib2
  end do

  if ( num /= n ) then

    indx = num-n+1

    do i = 1, n
      epstab(i)= epstab(indx)
      indx = indx+1
    end do

  end if

  if ( nres < 4 ) then
    res3la(nres) = result
    abserr = huge ( abserr )
  else
    abserr = abs(result-res3la(3))+abs(result-res3la(2)) &
      +abs(result-res3la(1))
    res3la(1) = res3la(2)
    res3la(2) = res3la(3)
    res3la(3) = result
  end if

  abserr = max ( abserr,0.5E+00* epsilon ( result ) *abs(result))

  return
end subroutine
subroutine qfour ( f_ptr, dat, a, b, omega, integr, epsabs, epsrel, limit, icall, &
  maxp1, result, abserr, neval, ier, alist, blist, rlist, elist, iord, &
  nnlog, momcom, chebmo )

!*****************************************************************************80
!
!! QFOUR estimates the integrals of oscillatory functions.
!
!  Discussion:
!
!    This routine calculates an approximation RESULT to a definite integral
!      I = integral of F(X) * COS(OMEGA*X) 
!    or
!      I = integral of F(X) * SIN(OMEGA*X) 
!    over (A,B), hopefully satisfying:
!      | I - RESULT | <= max ( epsabs, epsrel * |I| ) ).
!
!    QFOUR is called by QAWO and QAWF.  It can also be called directly in 
!    a user-written program.  In the latter case it is possible for the 
!    user to determine the first dimension of array CHEBMO(MAXP1,25).
!    See also parameter description of MAXP1.  Additionally see
!    parameter description of ICALL for eventually re-using
!    Chebyshev moments computed during former call on subinterval
!    of equal length abs(B-A).
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real ( kind = 8 ) F, the name of the function routine, of the form
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!    which evaluates the integrand function.
!
!    Input, real ( kind = 8 ) A, B, the limits of integration.
!
!    Input, real ( kind = 8 ) OMEGA, the multiplier of X in the weight function.
!
!    Input, integer ( kind = 8 ) INTEGR, indicates the weight functions to be used.
!    = 1, w(x) = cos(omega*x)
!    = 2, w(x) = sin(omega*x)
!
!    Input, real ( kind = 8 ) EPSABS, EPSREL, the absolute and relative accuracy requested.
!
!    Input, integer ( kind = 8 ) LIMIT, the maximum number of subintervals of [A,B]
!    that can be generated.
!
!    icall  - integer ( kind = 8 )
!                     if qfour is to be used only once, ICALL must
!                     be set to 1.  assume that during this call, the
!                     Chebyshev moments (for clenshaw-curtis integration
!                     of degree 24) have been computed for intervals of
!                     lenghts (abs(b-a))*2**(-l), l=0,1,2,...momcom-1.
!                     the Chebyshev moments already computed can be
!                     re-used in subsequent calls, if qfour must be
!                     called twice or more times on intervals of the
!                     same length abs(b-a). from the second call on, one
!                     has to put then ICALL > 1.
!                     if ICALL < 1, the routine will end with ier = 6.
!
!            maxp1  - integer ( kind = 8 )
!                     gives an upper bound on the number of
!                     Chebyshev moments which can be stored, i.e.
!                     for the intervals of lenghts abs(b-a)*2**(-l),
!                     l=0,1, ..., maxp1-2, maxp1 >= 1.
!                     if maxp1 < 1, the routine will end with ier = 6.
!                     increasing (decreasing) the value of maxp1
!                     decreases (increases) the computational time but
!                     increases (decreases) the required memory space.
!
!    Output, real ( kind = 8 ) RESULT, the estimated value of the integral.
!
!    Output, real ( kind = 8 ) ABSERR, an estimate of || I - RESULT ||.
!
!    Output, integer ( kind = 8 ) NEVAL, the number of times the integral was evaluated.
!
!            ier    - integer ( kind = 8 )
!                     ier = 0 normal and reliable termination of the
!                             routine. it is assumed that the
!                             requested accuracy has been achieved.
!                   - ier > 0 abnormal termination of the routine.
!                             the estimates for integral and error are
!                             less reliable. it is assumed that the
!                             requested accuracy has not been achieved.
!                     ier = 1 maximum number of subdivisions allowed
!                             has been achieved. one can allow more
!                             subdivisions by increasing the value of
!                             limit (and taking according dimension
!                             adjustments into account). however, if
!                             this yields no improvement it is advised
!                             to analyze the integrand, in order to
!                             determine the integration difficulties.
!                             if the position of a local difficulty can
!                             be determined (e.g. singularity,
!                             discontinuity within the interval) one
!                             will probably gain from splitting up the
!                             interval at this point and calling the
!                             integrator on the subranges. if possible,
!                             an appropriate special-purpose integrator
!                             should be used which is designed for
!                             handling the type of difficulty involved.
!                         = 2 the occurrence of roundoff error is
!                             detected, which prevents the requested
!                             tolerance from being achieved.
!                             the error may be under-estimated.
!                         = 3 extremely bad integrand behavior occurs
!                             at some points of the integration
!                             interval.
!                         = 4 the algorithm does not converge. roundoff
!                             error is detected in the extrapolation
!                             table. it is presumed that the requested
!                             tolerance cannot be achieved due to
!                             roundoff in the extrapolation table, and
!                             that the returned result is the best which
!                             can be obtained.
!                         = 5 the integral is probably divergent, or
!                             slowly convergent. it must be noted that
!                             divergence can occur with any other value
!                             of ier > 0.
!                         = 6 the input is invalid, because
!                             epsabs < 0 and epsrel < 0,
!                             or (integr /= 1 and integr /= 2) or
!                             ICALL < 1 or maxp1 < 1.
!                             result, abserr, neval, last, rlist(1),
!                             elist(1), iord(1) and nnlog(1) are set to
!                             zero. alist(1) and blist(1) are set to a
!                             and b respectively.
!
!    Workspace, real ( kind = 8 ) ALIST(LIMIT), BLIST(LIMIT), contains in entries 1 
!    through LAST the left and right ends of the partition subintervals.
!
!    Workspace, real ( kind = 8 ) RLIST(LIMIT), contains in entries 1 through LAST
!    the integral approximations on the subintervals.
!
!    Workspace, real ( kind = 8 ) ELIST(LIMIT), contains in entries 1 through LAST
!    the absolute error estimates on the subintervals.
!
!            iord   - integer ( kind = 8 )
!                     vector of dimension at least limit, the first k
!                     elements of which are pointers to the error
!                     estimates over the subintervals, such that
!                     elist(iord(1)), ..., elist(iord(k)), form
!                     a decreasing sequence, with k = last
!                     if last <= (limit/2+2), and
!                     k = limit+1-last otherwise.
!
!            nnlog  - integer ( kind = 8 )
!                     vector of dimension at least limit, indicating the
!                     subdivision levels of the subintervals, i.e.
!                     iwork(i) = l means that the subinterval numbered
!                     i is of length abs(b-a)*2**(1-l)
!
!         on entry and return
!            momcom - integer ( kind = 8 )
!                     indicating that the Chebyshev moments have been
!                     computed for intervals of lengths
!                     (abs(b-a))*2**(-l), l=0,1,2, ..., momcom-1,
!                     momcom < maxp1
!
!            chebmo - real
!                     array of dimension (maxp1,25) containing the
!                     Chebyshev moments
!
!  Local Parameters:
!
!           alist     - list of left end points of all subintervals
!                       considered up to now
!           blist     - list of right end points of all subintervals
!                       considered up to now
!           rlist(i)  - approximation to the integral over
!                       (alist(i),blist(i))
!           rlist2    - array of dimension at least limexp+2 containing
!                       the part of the epsilon table which is still
!                       needed for further computations
!           elist(i)  - error estimate applying to rlist(i)
!           maxerr    - pointer to the interval with largest error
!                       estimate
!           errmax    - elist(maxerr)
!           erlast    - error on the interval currently subdivided
!           area      - sum of the integrals over the subintervals
!           errsum    - sum of the errors over the subintervals
!           errbnd    - requested accuracy max(epsabs,epsrel*
!                       abs(result))
!           *****1    - variable for the left subinterval
!           *****2    - variable for the right subinterval
!           last      - index for subdivision
!           nres      - number of calls to the extrapolation routine
!           numrl2    - number of elements in rlist2. if an appropriate
!                       approximation to the compounded integral has
!                       been obtained it is put in rlist2(numrl2) after
!                       numrl2 has been increased by one
!           small     - length of the smallest interval considered
!                       up to now, multiplied by 1.5
!           erlarg    - sum of the errors over the intervals larger
!                       than the smallest interval considered up to now
!           extrap    - logical variable denoting that the routine is
!                       attempting to perform extrapolation, i.e. before
!                       subdividing the smallest interval we try to
!                       decrease the value of erlarg
!           noext     - logical variable denoting that extrapolation
!                       is no longer allowed (true value)
!
  implicit none

  integer ( kind = 8 ) limit
  integer ( kind = 8 ) maxp1

  real ( kind = 8 ) a
  real ( kind = 8 ) abseps
  real ( kind = 8 ) abserr
  real ( kind = 8 ) alist(limit)
  real ( kind = 8 ) area
  real ( kind = 8 ) area1
  real ( kind = 8 ) area12
  real ( kind = 8 ) area2
  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  real ( kind = 8 ) b
  real ( kind = 8 ) blist(limit)
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) chebmo(maxp1,25)
  real ( kind = 8 ) correc
  real ( kind = 8 ) defab1
  real ( kind = 8 ) defab2
  real ( kind = 8 ) defabs
  real ( kind = 8 ) domega
  real ( kind = 8 ) dres
  real ( kind = 8 ) elist(limit)
  real ( kind = 8 ) epsabs
  real ( kind = 8 ) epsrel
  real ( kind = 8 ) erlarg
  real ( kind = 8 ) erlast
  real ( kind = 8 ) errbnd
  real ( kind = 8 ) errmax
  real ( kind = 8 ) error1
  real ( kind = 8 ) erro12
  real ( kind = 8 ) error2
  real ( kind = 8 ) errsum
  real ( kind = 8 ) ertest
  logical extall
  logical extrap
  !real ( kind = 8 ), external :: f
  procedure (f_int_dat), intent(in), pointer :: f_ptr 
  integer ( kind = 8 ) icall
  integer ( kind = 8 ) id
  integer ( kind = 8 ) ier
  integer ( kind = 8 ) ierro
  integer ( kind = 8 ) integr
  integer ( kind = 8 ) iord(limit)
  integer ( kind = 8 ) iroff1
  integer ( kind = 8 ) iroff2
  integer ( kind = 8 ) iroff3
  integer ( kind = 8 ) jupbnd
  integer ( kind = 8 ) k
  integer ( kind = 8 ) ksgn
  integer ( kind = 8 ) ktmin
  integer ( kind = 8 ) last
  integer ( kind = 8 ) maxerr
  integer ( kind = 8 ) momcom
  integer ( kind = 8 ) nev
  integer ( kind = 8 ) neval
  integer ( kind = 8 ) nnlog(limit)
  logical noext
  integer ( kind = 8 ) nres
  integer ( kind = 8 ) nrmax
  integer ( kind = 8 ) nrmom
  integer ( kind = 8 ) numrl2
  real ( kind = 8 ) omega
  real ( kind = 8 ) resabs
  real ( kind = 8 ) reseps
  real ( kind = 8 ) result
  real ( kind = 8 ) res3la(3)
  real ( kind = 8 ) rlist(limit)
  real ( kind = 8 ) rlist2(52)
  real ( kind = 8 ) small
  real ( kind = 8 ) width
  class( dataCollectionBase), target :: dat
!
!  the dimension of rlist2 is determined by  the value of
!  limexp in QEXTR (rlist2 should be of dimension
!  (limexp+2) at least).
!
!  Test on validity of parameters.
!
  ier = 0
  neval = 0
  last = 0
  result = 0.0E+00
  abserr = 0.0E+00
  alist(1) = a
  blist(1) = b
  rlist(1) = 0.0E+00
  elist(1) = 0.0E+00
  iord(1) = 0
  nnlog(1) = 0

  if ( (integr /= 1.and.integr /= 2) .or. (epsabs < 0.0E+00.and. &
    epsrel < 0.0E+00) .or. icall < 1 .or. maxp1 < 1 ) then
    ier = 6
    return
  end if
!
!  First approximation to the integral.
!
  domega = abs ( omega )
  nrmom = 0

  if ( icall <= 1 ) then
    momcom = 0
  end if

  call qc25o ( f_ptr, dat, a, b, domega, integr, nrmom, maxp1, int8(0), result, abserr, &
    neval, defabs, resabs, momcom, chebmo )
!
!  Test on accuracy.
!
  dres = abs(result)
  errbnd = max ( epsabs,epsrel*dres)
  rlist(1) = result
  elist(1) = abserr
  iord(1) = 1
  if ( abserr <= 1.0E+02* epsilon ( defabs ) *defabs .and. &
    abserr > errbnd ) ier = 2

  if ( limit == 1 ) then
    ier = 1
  end if

  if ( ier /= 0 .or. abserr <= errbnd ) then
    go to 200
  end if
!
!  Initializations
!
  errmax = abserr
  maxerr = 1
  area = result
  errsum = abserr
  abserr = huge ( abserr )
  nrmax = 1
  extrap = .false.
  noext = .false.
  ierro = 0
  iroff1 = 0
  iroff2 = 0
  iroff3 = 0
  ktmin = 0
  small = abs(b-a)*7.5E-01
  nres = 0
  numrl2 = 0
  extall = .false.

  if ( 5.0E-01*abs(b-a)*domega <= 2.0E+00) then
    numrl2 = 1
    extall = .true.
    rlist2(1) = result
  end if

  if ( 2.5E-01 * abs(b-a) * domega <= 2.0E+00 ) then
    extall = .true.
  end if

  if ( dres >= (1.0E+00-5.0E+01* epsilon ( defabs ) )*defabs ) then
    ksgn = 1
  else
    ksgn = -1
  end if
!
!  main do-loop
!
  do last = 2, limit
!
!  Bisect the subinterval with the nrmax-th largest error estimate.
!
    nrmom = nnlog(maxerr)+1
    a1 = alist(maxerr)
    b1 = 5.0E-01*(alist(maxerr)+blist(maxerr))
    a2 = b1
    b2 = blist(maxerr)
    erlast = errmax

    call qc25o ( f_ptr, dat, a1, b1, domega, integr, nrmom, maxp1, int8(0), area1, &
      error1, nev, resabs, defab1, momcom, chebmo )

    neval = neval+nev

    call qc25o ( f_ptr, dat, a2, b2, domega, integr, nrmom, maxp1, int8(1), area2, &
      error2, nev, resabs, defab2, momcom, chebmo )

    neval = neval+nev
!
!  Improve previous approximations to integral and error and
!  test for accuracy.
!
    area12 = area1+area2
    erro12 = error1+error2
    errsum = errsum+erro12-errmax
    area = area+area12-rlist(maxerr)
    if ( defab1 == error1 .or. defab2 == error2 ) go to 25
    if ( abs(rlist(maxerr)-area12) > 1.0E-05*abs(area12) &
    .or. erro12 < 9.9E-01*errmax ) go to 20
    if ( extrap ) iroff2 = iroff2+1

    if ( .not.extrap ) then
      iroff1 = iroff1+1
    end if

20  continue

    if ( last > 10.and.erro12 > errmax ) iroff3 = iroff3+1

25  continue

    rlist(maxerr) = area1
    rlist(last) = area2
    nnlog(maxerr) = nrmom
    nnlog(last) = nrmom
    errbnd = max ( epsabs,epsrel*abs(area))
!
!  Test for roundoff error and eventually set error flag
!
    if ( iroff1+iroff2 >= 10 .or. iroff3 >= 20 ) ier = 2

    if ( iroff2 >= 5) ierro = 3
!
!  Set error flag in the case that the number of subintervals
!  equals limit.
!
    if ( last == limit ) then
      ier = 1
    end if
!
!  Set error flag in the case of bad integrand behavior at
!  a point of the integration range.
!
    if ( max ( abs(a1),abs(b2)) <= (1.0E+00+1.0E+03* epsilon ( a1 ) ) &
    *(abs(a2)+1.0E+03* tiny ( a2 ) )) then
      ier = 4
    end if
!
!  Append the newly-created intervals to the list.
!
    if ( error2 <= error1 ) then
      alist(last) = a2
      blist(maxerr) = b1
      blist(last) = b2
      elist(maxerr) = error1
      elist(last) = error2
    else
      alist(maxerr) = a2
      alist(last) = a1
      blist(last) = b1
      rlist(maxerr) = area2
      rlist(last) = area1
      elist(maxerr) = error2
      elist(last) = error1
    end if
!
!  Call QSORT to maintain the descending ordering
!  in the list of error estimates and select the subinterval
!  with nrmax-th largest error estimate (to be bisected next).
!

    call qsort ( limit, last, maxerr, errmax, elist, iord, nrmax )

    if ( errsum <= errbnd ) then
      go to 170
    end if

    if ( ier /= 0 ) then
      exit
    end if

    if ( last == 2 .and. extall ) go to 120

    if ( noext ) then
      cycle
    end if

    if ( .not. extall ) go to 50
    erlarg = erlarg-erlast
    if ( abs(b1-a1) > small ) erlarg = erlarg+erro12
    if ( extrap ) go to 70
!
!  Test whether the interval to be bisected next is the
!  smallest interval.
!
50  continue

    width = abs(blist(maxerr)-alist(maxerr))

    if ( width > small ) then
      cycle
    end if

    if ( extall ) go to 60
!
!  Test whether we can start with the extrapolation procedure
!  (we do this if we integrate over the next interval with
!  use of a Gauss-Kronrod rule - see QC25O).
!
    small = small*5.0E-01

    if ( 2.5E-01*width*domega > 2.0E+00 ) then
      cycle
    end if

    extall = .true.
    go to 130

60  continue

    extrap = .true.
    nrmax = 2

70  continue

    if ( ierro == 3 .or. erlarg <= ertest ) go to 90
!
!  The smallest interval has the largest error.
!  Before bisecting decrease the sum of the errors over the
!  larger intervals (ERLARG) and perform extrapolation.
!
    jupbnd = last

    if ( last > (limit/2+2) ) then
      jupbnd = limit+3-last
    end if

    id = nrmax

    do k = id, jupbnd
      maxerr = iord(nrmax)
      errmax = elist(maxerr)
      if ( abs(blist(maxerr)-alist(maxerr)) > small ) go to 140
      nrmax = nrmax+1
    end do
!
!  Perform extrapolation.
!
90  continue

    numrl2 = numrl2+1
    rlist2(numrl2) = area

    if ( numrl2 < 3 ) go to 110

    call qextr ( numrl2, rlist2, reseps, abseps, res3la, nres )
    ktmin = ktmin+1

    if ( ktmin > 5.and.abserr < 1.0E-03*errsum ) then
      ier = 5
    end if

    if ( abseps >= abserr ) go to 100

    ktmin = 0
    abserr = abseps
    result = reseps
    correc = erlarg
    ertest = max ( epsabs, epsrel*abs(reseps))

    if ( abserr <= ertest ) then
      exit
    end if
!
!  Prepare bisection of the smallest interval.
!
100 continue

    if ( numrl2 == 1 ) then
      noext = .true.
    end if

    if ( ier == 5 ) then
      exit
    end if

110 continue

    maxerr = iord(1)
    errmax = elist(maxerr)
    nrmax = 1
    extrap = .false.
    small = small*5.0E-01
    erlarg = errsum
    cycle

120 continue

    small = small * 5.0E-01
    numrl2 = numrl2 + 1
    rlist2(numrl2) = area

130 continue

    ertest = errbnd
    erlarg = errsum

140 continue

  end do
!
!  set the final result.
!
  if ( abserr == huge ( abserr ) .or. nres == 0 ) then
    go to 170
  end if

  if ( ier+ierro == 0 ) go to 165
  if ( ierro == 3 ) abserr = abserr+correc
  if ( ier == 0 ) ier = 3
  if ( result /= 0.0E+00.and.area /= 0.0E+00 ) go to 160
  if ( abserr > errsum ) go to 170
  if ( area == 0.0E+00 ) go to 190
  go to 165

160 continue

  if ( abserr/abs(result) > errsum/abs(area) ) go to 170
!
!  Test on divergence.
!
  165 continue

  if ( ksgn == (-1) .and. max ( abs(result),abs(area)) <=  &
   defabs*1.0E-02 ) go to 190

  if ( 1.0E-02 > (result/area) .or. (result/area) > 1.0E+02 &
   .or. errsum >= abs(area) ) ier = 6

  go to 190
!
!  Compute global integral sum.
!
170 continue

  result = sum ( rlist(1:last) )

  abserr = errsum

190 continue

  if (ier > 2) ier=ier-1

200 continue

  if ( integr == 2 .and. omega < 0.0E+00 ) then
    result = -result
  end if

  return
end subroutine
subroutine qk15 ( f_ptr, dat, a, b, result, abserr, resabs, resasc )

!*****************************************************************************80
!
!! QK15 carries out a 15 point Gauss-Kronrod quadrature rule.
!
!  Discussion:
!
!    This routine approximates
!      I = integral ( A <= X <= B ) F(X) dx
!    with an error estimate, and
!      J = integral ( A <= X <= B ) | F(X) | dx
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real ( kind = 8 ) F, the name of the function routine, of the form
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!    which evaluates the integrand function.
!
!    Input, real ( kind = 8 ) A, B, the limits of integration.
!
!    Output, real ( kind = 8 ) RESULT, the estimated value of the integral.
!    RESULT is computed by applying the 15-point Kronrod rule (RESK) 
!    obtained by optimal addition of abscissae to the 7-point Gauss rule 
!    (RESG).
!
!    Output, real ( kind = 8 ) ABSERR, an estimate of | I - RESULT |.
!
!    Output, real ( kind = 8 ) RESABS, approximation to the integral of the absolute
!    value of F.
!
!    Output, real ( kind = 8 ) RESASC, approximation to the integral | F-I/(B-A) | 
!    over [A,B].
!
!  Local Parameters:
!
!           the abscissae and weights are given for the interval (-1,1).
!           because of symmetry only the positive abscissae and their
!           corresponding weights are given.
!
!           xgk    - abscissae of the 15-point Kronrod rule
!                    xgk(2), xgk(4), ...  abscissae of the 7-point
!                    Gauss rule
!                    xgk(1), xgk(3), ...  abscissae which are optimally
!                    added to the 7-point Gauss rule
!
!           wgk    - weights of the 15-point Kronrod rule
!
!           wg     - weights of the 7-point Gauss rule
!
!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc   - abscissa
!           fval*  - function value
!           resg   - result of the 7-point Gauss formula
!           resk   - result of the 15-point Kronrod formula
!           reskh  - approximation to the mean value of f over (a,b),
!                    i.e. to i/(b-a)
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) absc
  real ( kind = 8 ) abserr
  real ( kind = 8 ) b
  real ( kind = 8 ) centr
  real ( kind = 8 ) dhlgth
  !real ( kind = 8 ), external :: f
  procedure (f_int_dat), intent(in), pointer :: f_ptr 
  real ( kind = 8 ) fc
  real ( kind = 8 ) fsum
  real ( kind = 8 ) fval1
  real ( kind = 8 ) fval2
  real ( kind = 8 ) fv1(7)
  real ( kind = 8 ) fv2(7)
  real ( kind = 8 ) hlgth
  integer ( kind = 8 ) j
  integer ( kind = 8 ) jtw
  integer ( kind = 8 ) jtwm1
  real ( kind = 8 ) resabs
  real ( kind = 8 ) resasc
  real ( kind = 8 ) resg
  real ( kind = 8 ) resk
  real ( kind = 8 ) reskh
  real ( kind = 8 ) result
  real ( kind = 8 ) wg(4)
  real ( kind = 8 ) wgk(8)
  real ( kind = 8 ) xgk(8)
  class( dataCollectionBase), target :: dat

  data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8)/ &
       9.914553711208126E-01,   9.491079123427585E-01, &
       8.648644233597691E-01,   7.415311855993944E-01, &
       5.860872354676911E-01,   4.058451513773972E-01, &
       2.077849550078985E-01,   0.0E+00              /
  data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8)/ &
       2.293532201052922E-02,   6.309209262997855E-02, &
       1.047900103222502E-01,   1.406532597155259E-01, &
       1.690047266392679E-01,   1.903505780647854E-01, &
       2.044329400752989E-01,   2.094821410847278E-01/
  data wg(1),wg(2),wg(3),wg(4)/ &
       1.294849661688697E-01,   2.797053914892767E-01, &
       3.818300505051189E-01,   4.179591836734694E-01/
!
  centr = 5.0E-01*(a+b)
  hlgth = 5.0E-01*(b-a)
  dhlgth = abs(hlgth)
  
!
!  Compute the 15-point Kronrod approximation to the integral,
!  and estimate the absolute error.
!
  fc = f_ptr(centr,dat)
  resg = fc*wg(4)
  resk = fc*wgk(8)
  resabs = abs(resk)

  do j = 1, 3
    jtw = j*2
    absc = hlgth*xgk(jtw)
    fval1 = f_ptr(centr-absc,dat)
    fval2 = f_ptr(centr+absc,dat)
    fv1(jtw) = fval1
    fv2(jtw) = fval2
    fsum = fval1+fval2
    resg = resg+wg(j)*fsum
    resk = resk+wgk(jtw)*fsum
    resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
  end do

  do j = 1, 4
    jtwm1 = j*2-1
    absc = hlgth*xgk(jtwm1)
    fval1 = f_ptr(centr-absc,dat)
    fval2 = f_ptr(centr+absc,dat)
    fv1(jtwm1) = fval1
    fv2(jtwm1) = fval2
    fsum = fval1+fval2
    resk = resk+wgk(jtwm1)*fsum
    resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
  end do

  reskh = resk * 5.0E-01
  resasc = wgk(8)*abs(fc-reskh)

  do j = 1, 7
    resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
  end do

  result = resk*hlgth
  resabs = resabs*dhlgth
  resasc = resasc*dhlgth
  abserr = abs((resk-resg)*hlgth)

  if ( resasc /= 0.0E+00.and.abserr /= 0.0E+00 ) then
    abserr = resasc*min ( 1.0E+00,(2.0E+02*abserr/resasc)**1.5E+00)
  end if

  if ( resabs > tiny ( resabs ) / (5.0E+01* epsilon ( resabs ) ) ) then
    abserr = max (( epsilon ( resabs ) *5.0E+01)*resabs,abserr)
  end if

  return
end subroutine
subroutine qk15i ( f_ptr, dat, boun, inf, a, b, result, abserr, resabs, resasc )

!*****************************************************************************80
!
!! QK15I applies a 15 point Gauss-Kronrod quadrature on an infinite interval.
!
!  Discussion:
!
!    The original infinite integration range is mapped onto the interval 
!    (0,1) and (a,b) is a part of (0,1).  The routine then computes:
!
!    i = integral of transformed integrand over (a,b),
!    j = integral of abs(transformed integrand) over (a,b).
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real ( kind = 8 ) F, the name of the function routine, of the form
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!    which evaluates the integrand function.
!
!    Input, real ( kind = 8 ) BOUN, the finite bound of the original integration range,
!    or zero if INF is 2.
!
!    Input, integer ( kind = 8 ) INF, indicates the type of the interval.
!    -1: the original interval is (-infinity,BOUN),
!    +1, the original interval is (BOUN,+infinity),
!    +2, the original interval is (-infinity,+infinity) and
!    the integral is computed as the sum of two integrals, one 
!    over (-infinity,0) and one over (0,+infinity).
!
!    Input, real ( kind = 8 ) A, B, the limits of integration, over a subrange of [0,1].
!
!    Output, real ( kind = 8 ) RESULT, the estimated value of the integral.
!    RESULT is computed by applying the 15-point Kronrod rule (RESK) obtained 
!    by optimal addition of abscissae to the 7-point Gauss rule (RESG).
!
!    Output, real ( kind = 8 ) ABSERR, an estimate of | I - RESULT |.
!
!    Output, real ( kind = 8 ) RESABS, approximation to the integral of the absolute
!    value of F.
!
!    Output, real ( kind = 8 ) RESASC, approximation to the integral of the
!    transformated integrand | F-I/(B-A) | over [A,B].
!
!  Local Parameters:
!
!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc*  - abscissa
!           tabsc* - transformed abscissa
!           fval*  - function value
!           resg   - result of the 7-point Gauss formula
!           resk   - result of the 15-point Kronrod formula
!           reskh  - approximation to the mean value of the transformed
!                    integrand over (a,b), i.e. to i/(b-a)
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) absc
  real ( kind = 8 ) absc1
  real ( kind = 8 ) absc2
  real ( kind = 8 ) abserr
  real ( kind = 8 ) b
  real ( kind = 8 ) boun
  real ( kind = 8 ) centr
  real ( kind = 8 ) dinf
  !real ( kind = 8 ), external :: f
  procedure (f_int_dat), intent(in), pointer :: f_ptr 
  real ( kind = 8 ) fc
  real ( kind = 8 ) fsum
  real ( kind = 8 ) fval1
  real ( kind = 8 ) fval2
  real ( kind = 8 ) fv1(7)
  real ( kind = 8 ) fv2(7)
  real ( kind = 8 ) hlgth
  integer ( kind = 8 ) inf
  integer ( kind = 8 ) j
  real ( kind = 8 ) resabs
  real ( kind = 8 ) resasc
  real ( kind = 8 ) resg
  real ( kind = 8 ) resk
  real ( kind = 8 ) reskh
  real ( kind = 8 ) result
  real ( kind = 8 ) tabsc1
  real ( kind = 8 ) tabsc2
  real ( kind = 8 ) wg(8)
  real ( kind = 8 ) wgk(8)
  real ( kind = 8 ) xgk(8)
  class( dataCollectionBase), target :: dat
!
!  the abscissae and weights are supplied for the interval
!  (-1,1).  because of symmetry only the positive abscissae and
!  their corresponding weights are given.
!
!           xgk    - abscissae of the 15-point Kronrod rule
!                    xgk(2), xgk(4), ... abscissae of the 7-point Gauss
!                    rule
!                    xgk(1), xgk(3), ...  abscissae which are optimally
!                    added to the 7-point Gauss rule
!
!           wgk    - weights of the 15-point Kronrod rule
!
!           wg     - weights of the 7-point Gauss rule, corresponding
!                    to the abscissae xgk(2), xgk(4), ...
!                    wg(1), wg(3), ... are set to zero.
!
  data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8)/ &
       9.914553711208126E-01,     9.491079123427585E-01, &
       8.648644233597691E-01,     7.415311855993944E-01, &
       5.860872354676911E-01,     4.058451513773972E-01, &
       2.077849550078985E-01,     0.0000000000000000E+00/

  data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8)/ &
       2.293532201052922E-02,     6.309209262997855E-02, &
       1.047900103222502E-01,     1.406532597155259E-01, &
       1.690047266392679E-01,     1.903505780647854E-01, &
       2.044329400752989E-01,     2.094821410847278E-01/

  data wg(1),wg(2),wg(3),wg(4),wg(5),wg(6),wg(7),wg(8)/ &
       0.0000000000000000E+00,     1.294849661688697E-01, &
       0.0000000000000000E+00,     2.797053914892767E-01, &
       0.0000000000000000E+00,     3.818300505051189E-01, &
       0.0000000000000000E+00,     4.179591836734694E-01/

  dinf = min ( 1, inf )
  
  centr = 5.0E-01*(a+b)
  hlgth = 5.0E-01*(b-a)
  tabsc1 = boun+dinf*(1.0E+00-centr)/centr
  fval1 = f_ptr(tabsc1,dat)
  if ( inf == 2 ) fval1 = fval1+f_ptr(-tabsc1,dat)
  fc = (fval1/centr)/centr
!
!  Compute the 15-point Kronrod approximation to the integral,
!  and estimate the error.
!
  resg = wg(8)*fc
  resk = wgk(8)*fc
  resabs = abs(resk)

  do j = 1, 7

    absc = hlgth*xgk(j)
    absc1 = centr-absc
    absc2 = centr+absc
    tabsc1 = boun+dinf*(1.0E+00-absc1)/absc1
    tabsc2 = boun+dinf*(1.0E+00-absc2)/absc2
    fval1 = f_ptr(tabsc1,dat)
    fval2 = f_ptr(tabsc2,dat)

    if ( inf == 2 ) then
      fval1 = fval1+f_ptr(-tabsc1,dat)
      fval2 = fval2+f_ptr(-tabsc2,dat)
    end if

    fval1 = (fval1/absc1)/absc1
    fval2 = (fval2/absc2)/absc2
    fv1(j) = fval1
    fv2(j) = fval2
    fsum = fval1+fval2
    resg = resg+wg(j)*fsum
    resk = resk+wgk(j)*fsum
    resabs = resabs+wgk(j)*(abs(fval1)+abs(fval2))
  end do

  reskh = resk * 5.0E-01
  resasc = wgk(8) * abs(fc-reskh)

  do j = 1, 7
    resasc = resasc + wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
  end do

  result = resk * hlgth
  resasc = resasc * hlgth
  resabs = resabs * hlgth
  abserr = abs ( ( resk - resg ) * hlgth )

  if ( resasc /= 0.0E+00.and.abserr /= 0.0E+00) then
    abserr = resasc* min ( 1.0E+00,(2.0E+02*abserr/resasc)**1.5E+00)
  end if

  if ( resabs > tiny ( resabs ) / ( 5.0E+01 * epsilon ( resabs ) ) ) then
    abserr = max (( epsilon ( resabs ) *5.0E+01)*resabs,abserr)
  end if

  return
end subroutine
subroutine qk15w ( f_ptr, dat, w, p1, p2, p3, p4, kp, a, b, result, abserr, resabs, &
  resasc )

!*****************************************************************************80
!
!! QK15W applies a 15 point Gauss-Kronrod rule for a weighted integrand.
!
!  Discussion:
!
!    This routine approximates 
!      i = integral of f*w over (a,b), 
!    with error estimate, and
!      j = integral of abs(f*w) over (a,b)
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real ( kind = 8 ) F, the name of the function routine, of the form
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!    which evaluates the integrand function.
!
!              w      - real
!                       function subprogram defining the integrand
!                       weight function w(x). the actual name for w
!                       needs to be declared e x t e r n a l in the
!                       calling program.
!
!    ?, real ( kind = 8 ) P1, P2, P3, P4, parameters in the weight function
!
!    Input, integer ( kind = 8 ) KP, key for indicating the type of weight function
!
!    Input, real ( kind = 8 ) A, B, the limits of integration.
!
!    Output, real ( kind = 8 ) RESULT, the estimated value of the integral.
!    RESULT is computed by applying the 15-point Kronrod rule (RESK) obtained by
!    optimal addition of abscissae to the 7-point Gauss rule (RESG).
!
!    Output, real ( kind = 8 ) ABSERR, an estimate of | I - RESULT |.
!
!    Output, real ( kind = 8 ) RESABS, approximation to the integral of the absolute
!    value of F.
!
!    Output, real ( kind = 8 ) RESASC, approximation to the integral | F-I/(B-A) | 
!    over [A,B].
!
!  Local Parameters:
!
!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc*  - abscissa
!           fval*  - function value
!           resg   - result of the 7-point Gauss formula
!           resk   - result of the 15-point Kronrod formula
!           reskh  - approximation to the mean value of f*w over (a,b),
!                    i.e. to i/(b-a)
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) absc
  real ( kind = 8 ) absc1
  real ( kind = 8 ) absc2
  real ( kind = 8 ) abserr
  real ( kind = 8 ) b
  real ( kind = 8 ) centr
  real ( kind = 8 ) dhlgth
  !real ( kind = 8 ), external :: f
  procedure (f_int_dat), intent(in), pointer :: f_ptr 
  real ( kind = 8 ) fc
  real ( kind = 8 ) fsum
  real ( kind = 8 ) fval1
  real ( kind = 8 ) fval2
  real ( kind = 8 ) fv1(7)
  real ( kind = 8 ) fv2(7)
  real ( kind = 8 ) hlgth
  integer ( kind = 8 ) j
  integer ( kind = 8 ) jtw
  integer ( kind = 8 ) jtwm1
  integer ( kind = 8 ) kp
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  real ( kind = 8 ) p3
  real ( kind = 8 ) p4
  real ( kind = 8 ) resabs
  real ( kind = 8 ) resasc
  real ( kind = 8 ) resg
  real ( kind = 8 ) resk
  real ( kind = 8 ) reskh
  real ( kind = 8 ) result
  real ( kind = 8 ), external :: w
  real ( kind = 8 ), dimension ( 4 ) :: wg = (/ &
    1.294849661688697E-01,     2.797053914892767E-01, &
    3.818300505051889E-01,     4.179591836734694E-01 /)
  real ( kind = 8 ) wgk(8)
  real ( kind = 8 ) xgk(8)
  class( dataCollectionBase), target :: dat
!
!  the abscissae and weights are given for the interval (-1,1).
!  because of symmetry only the positive abscissae and their
!  corresponding weights are given.
!
!           xgk    - abscissae of the 15-point Gauss-Kronrod rule
!                    xgk(2), xgk(4), ... abscissae of the 7-point Gauss
!                    rule
!                    xgk(1), xgk(3), ... abscissae which are optimally
!                    added to the 7-point Gauss rule
!
!           wgk    - weights of the 15-point Gauss-Kronrod rule
!
!           wg     - weights of the 7-point Gauss rule
!
  data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8)/ &
       9.914553711208126E-01,     9.491079123427585E-01, &
       8.648644233597691E-01,     7.415311855993944E-01, &
       5.860872354676911E-01,     4.058451513773972E-01, &
       2.077849550789850E-01,     0.000000000000000E+00/

  data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8)/ &
       2.293532201052922E-02,     6.309209262997855E-02, &
       1.047900103222502E-01,     1.406532597155259E-01, &
       1.690047266392679E-01,     1.903505780647854E-01, &
       2.044329400752989E-01,     2.094821410847278E-01/
!
  centr = 5.0E-01*(a+b)
  hlgth = 5.0E-01*(b-a)
  dhlgth = abs(hlgth)
  
!
!  Compute the 15-point Kronrod approximation to the integral,
!  and estimate the error.
!
  fc = f_ptr(centr,dat)*w(centr,p1,p2,p3,p4,kp)
  resg = wg(4)*fc
  resk = wgk(8)*fc
  resabs = abs(resk)

  do j = 1, 3
    jtw = j*2
    absc = hlgth*xgk(jtw)
    absc1 = centr-absc
    absc2 = centr+absc
    fval1 = f_ptr(absc1,dat)*w(absc1,p1,p2,p3,p4,kp)
    fval2 = f_ptr(absc2,dat)*w(absc2,p1,p2,p3,p4,kp)
    fv1(jtw) = fval1
    fv2(jtw) = fval2
    fsum = fval1+fval2
    resg = resg+wg(j)*fsum
    resk = resk+wgk(jtw)*fsum
    resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
  end do

  do j = 1, 4
    jtwm1 = j*2-1
    absc = hlgth*xgk(jtwm1)
    absc1 = centr-absc
    absc2 = centr+absc
    fval1 = f_ptr(absc1,dat)*w(absc1,p1,p2,p3,p4,kp)
    fval2 = f_ptr(absc2,dat)*w(absc2,p1,p2,p3,p4,kp)
    fv1(jtwm1) = fval1
    fv2(jtwm1) = fval2
    fsum = fval1+fval2
    resk = resk+wgk(jtwm1)*fsum
    resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
  end do

  reskh = resk*5.0E-01
  resasc = wgk(8)*abs(fc-reskh)

  do j = 1, 7
    resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
  end do

  result = resk*hlgth
  resabs = resabs*dhlgth
  resasc = resasc*dhlgth
  abserr = abs((resk-resg)*hlgth)

  if ( resasc /= 0.0E+00.and.abserr /= 0.0E+00) then
    abserr = resasc*min ( 1.0E+00,(2.0E+02*abserr/resasc)**1.5E+00)
  end if

  if ( resabs > tiny ( resabs ) /(5.0E+01* epsilon ( resabs ) ) ) then
    abserr = max ( ( epsilon ( resabs ) * 5.0E+01)*resabs,abserr)
  end if

  return
end subroutine
subroutine qk21_x ( f_ptr, dat, a, b, result, abserr, resabs, resasc )

!*****************************************************************************80
!
!! QK21 carries out a 21 point Gauss-Kronrod quadrature rule.
!
!  Discussion:
!
!    This routine approximates
!      I = integral ( A <= X <= B ) F(X) dx
!    with an error estimate, and
!      J = integral ( A <= X <= B ) | F(X) | dx
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real ( kind = 8 ) F, the name of the function routine, of the form
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!    which evaluates the integrand function.
!
!    Input, real ( kind = 8 ) A, B, the limits of integration.
!
!    Output, real ( kind = 8 ) RESULT, the estimated value of the integral.
!    RESULT is computed by applying the 21-point Kronrod rule (resk) 
!    obtained by optimal addition of abscissae to the 10-point Gauss 
!    rule (resg).
!
!    Output, real ( kind = 8 ) ABSERR, an estimate of | I - RESULT |.
!
!    Output, real ( kind = 8 ) RESABS, approximation to the integral of the absolute
!    value of F.
!
!    Output, real ( kind = 8 ) RESASC, approximation to the integral | F-I/(B-A) | 
!    over [A,B].
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) absc
  real ( kind = 8 ) abserr
  real ( kind = 8 ) b
  real ( kind = 8 ) centr
  real ( kind = 8 ) dhlgth
  !real ( kind = 8 ), external :: f
  procedure (f_int_dat),intent(in), pointer :: f_ptr 
  real ( kind = 8 ) fc
  real ( kind = 8 ) fsum
  real ( kind = 8 ) fval1
  real ( kind = 8 ) fval2
  real ( kind = 8 ) fv1(10)
  real ( kind = 8 ) fv2(10)
  real ( kind = 8 ) hlgth
  integer ( kind = 8 ) j
  integer ( kind = 8 ) jtw
  integer ( kind = 8 ) jtwm1
  real ( kind = 8 ) resabs
  real ( kind = 8 ) resasc
  real ( kind = 8 ) resg
  real ( kind = 8 ) resk
  real ( kind = 8 ) reskh
  real ( kind = 8 ) result
  real ( kind = 8 ) wg(5)
  real ( kind = 8 ) wgk(11)
  real ( kind = 8 ) xgk(11)
  class(dataCollectionBase), target :: dat
  
 ! f_ptr => f
!
!           the abscissae and weights are given for the interval (-1,1).
!           because of symmetry only the positive abscissae and their
!           corresponding weights are given.
!
!           xgk    - abscissae of the 21-point Kronrod rule
!                    xgk(2), xgk(4), ...  abscissae of the 10-point
!                    Gauss rule
!                    xgk(1), xgk(3), ...  abscissae which are optimally
!                    added to the 10-point Gauss rule
!
!           wgk    - weights of the 21-point Kronrod rule
!
!           wg     - weights of the 10-point Gauss rule
!
  data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8), &
    xgk(9),xgk(10),xgk(11)/ &
       9.956571630258081E-01,     9.739065285171717E-01, &
       9.301574913557082E-01,     8.650633666889845E-01, &
       7.808177265864169E-01,     6.794095682990244E-01, &
       5.627571346686047E-01,     4.333953941292472E-01, &
       2.943928627014602E-01,     1.488743389816312E-01, &
       0.000000000000000E+00/
!
  data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8), &
    wgk(9),wgk(10),wgk(11)/ &
       1.169463886737187E-02,     3.255816230796473E-02, &
       5.475589657435200E-02,     7.503967481091995E-02, &
       9.312545458369761E-02,     1.093871588022976E-01, &
       1.234919762620659E-01,     1.347092173114733E-01, &
       1.427759385770601E-01,     1.477391049013385E-01, &
       1.494455540029169E-01/
!
  data wg(1),wg(2),wg(3),wg(4),wg(5)/ &
       6.667134430868814E-02,     1.494513491505806E-01, &
       2.190863625159820E-01,     2.692667193099964E-01, &
       2.955242247147529E-01/
!
!
!           list of major variables
!
!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc   - abscissa
!           fval*  - function value
!           resg   - result of the 10-point Gauss formula
!           resk   - result of the 21-point Kronrod formula
!           reskh  - approximation to the mean value of f over (a,b),
!                    i.e. to i/(b-a)
!
  centr = 5.0E-01*(a+b)
  hlgth = 5.0E-01*(b-a)
  dhlgth = abs(hlgth)
!
!  Compute the 21-point Kronrod approximation to the
!  integral, and estimate the absolute error.
!
  resg = 0.0E+00
  fc = f_ptr(centr, dat)
  resk = wgk(11)*fc
  resabs = abs(resk)

  do j = 1, 5
    jtw = 2*j
    absc = hlgth*xgk(jtw)
    fval1 = f_ptr(centr-absc,dat)
    fval2 = f_ptr(centr+absc,dat)
    fv1(jtw) = fval1
    fv2(jtw) = fval2
    fsum = fval1+fval2
    resg = resg+wg(j)*fsum
    resk = resk+wgk(jtw)*fsum
    resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
  end do

  do j = 1, 5
    jtwm1 = 2*j-1
    absc = hlgth*xgk(jtwm1)
    fval1 = f_ptr(centr-absc,dat)
    fval2 = f_ptr(centr+absc,dat)
    fv1(jtwm1) = fval1
    fv2(jtwm1) = fval2
    fsum = fval1+fval2
    resk = resk+wgk(jtwm1)*fsum
    resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
  end do

  reskh = resk*5.0E-01
  resasc = wgk(11)*abs(fc-reskh)

  do j = 1, 10
    resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
  end do

  result = resk*hlgth
  resabs = resabs*dhlgth
  resasc = resasc*dhlgth
  abserr = abs((resk-resg)*hlgth)

  if ( resasc /= 0.0E+00.and.abserr /= 0.0E+00) then
    abserr = resasc*min ( 1.0E+00,(2.0E+02*abserr/resasc)**1.5E+00)
  end if

  if ( resabs > tiny ( resabs ) /(5.0E+01* epsilon ( resabs ) )) then
    abserr = max (( epsilon ( resabs ) *5.0E+01)*resabs,abserr)
  end if

  return
end subroutine qk21_x

subroutine qk21_x_vec ( f_ptr, dat, a, b, result, abserr, resabs, resasc )

!*****************************************************************************80
!
!! QK21 carries out a 21 point Gauss-Kronrod quadrature rule.
!
!! Modified version by Kaspar K. Nielsen, 25-1-2018. The routine has been vectorized so that the function pointer is called with 20 points in one go
!  Discussion:
!
!    This routine approximates
!      I = integral ( A <= X <= B ) F(X) dx
!    with an error estimate, and
!      J = integral ( A <= X <= B ) | F(X) | dx
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real ( kind = 8 ) F, the name of the function routine, of the form
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!    which evaluates the integrand function.
!
!    Input, real ( kind = 8 ) A, B, the limits of integration.
!
!    Output, real ( kind = 8 ) RESULT, the estimated value of the integral.
!    RESULT is computed by applying the 21-point Kronrod rule (resk) 
!    obtained by optimal addition of abscissae to the 10-point Gauss 
!    rule (resg).
!
!    Output, real ( kind = 8 ) ABSERR, an estimate of | I - RESULT |.
!
!    Output, real ( kind = 8 ) RESABS, approximation to the integral of the absolute
!    value of F.
!
!    Output, real ( kind = 8 ) RESASC, approximation to the integral | F-I/(B-A) | 
!    over [A,B].
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) absc
  real ( kind = 8 ) abserr
  real ( kind = 8 ) b
  real ( kind = 8 ) centr
  real ( kind = 8 ) dhlgth
  !real ( kind = 8 ), external :: f
  procedure (f_int_dat_vec),intent(in), pointer :: f_ptr 
  real ( kind = 8 ) fc
  real ( kind = 8 ) fsum
  real ( kind = 8 ) fval1
  real ( kind = 8 ) fval2
  real ( kind = 8 ) fv1(10)
  real ( kind = 8 ) fv2(10)
  real ( kind = 8 ) fval(21)
  real ( kind = 8 ) xval(21)
  real ( kind = 8 ) hlgth
  integer ( kind = 8 ) j
  integer ( kind = 8 ) jtw
  integer ( kind = 8 ) jtwm1
  real ( kind = 8 ) resabs
  real ( kind = 8 ) resasc
  real ( kind = 8 ) resg
  real ( kind = 8 ) resk
  real ( kind = 8 ) reskh
  real ( kind = 8 ) result
  real ( kind = 8 ) wg(5)
  real ( kind = 8 ) wgk(11)
  real ( kind = 8 ) xgk(11)
  class(dataCollectionBase), target :: dat
  
 ! f_ptr => f
!
!           the abscissae and weights are given for the interval (-1,1).
!           because of symmetry only the positive abscissae and their
!           corresponding weights are given.
!
!           xgk    - abscissae of the 21-point Kronrod rule
!                    xgk(2), xgk(4), ...  abscissae of the 10-point
!                    Gauss rule
!                    xgk(1), xgk(3), ...  abscissae which are optimally
!                    added to the 10-point Gauss rule
!
!           wgk    - weights of the 21-point Kronrod rule
!
!           wg     - weights of the 10-point Gauss rule
!
  data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8), &
    xgk(9),xgk(10),xgk(11)/ &
       9.956571630258081E-01,     9.739065285171717E-01, &
       9.301574913557082E-01,     8.650633666889845E-01, &
       7.808177265864169E-01,     6.794095682990244E-01, &
       5.627571346686047E-01,     4.333953941292472E-01, &
       2.943928627014602E-01,     1.488743389816312E-01, &
       0.000000000000000E+00/
!
  data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8), &
    wgk(9),wgk(10),wgk(11)/ &
       1.169463886737187E-02,     3.255816230796473E-02, &
       5.475589657435200E-02,     7.503967481091995E-02, &
       9.312545458369761E-02,     1.093871588022976E-01, &
       1.234919762620659E-01,     1.347092173114733E-01, &
       1.427759385770601E-01,     1.477391049013385E-01, &
       1.494455540029169E-01/
!
  data wg(1),wg(2),wg(3),wg(4),wg(5)/ &
       6.667134430868814E-02,     1.494513491505806E-01, &
       2.190863625159820E-01,     2.692667193099964E-01, &
       2.955242247147529E-01/
!
!
!           list of major variables
!
!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc   - abscissa
!           fval*  - function value
!           resg   - result of the 10-point Gauss formula
!           resk   - result of the 21-point Kronrod formula
!           reskh  - approximation to the mean value of f over (a,b),
!                    i.e. to i/(b-a)
!
  centr = 5.0E-01*(a+b)
  hlgth = 5.0E-01*(b-a)
  dhlgth = abs(hlgth)
!
!  Compute the 21-point Kronrod approximation to the
!  integral, and estimate the absolute error.
!
  resg = 0.0E+00
  
  
  do j = 1, 5
    jtw = 2*j
    absc = hlgth*xgk(jtw)    
    xval(j) = centr-absc
    xval(j+5) = centr+absc
    
    jtwm1 = 2*j-1
    absc = hlgth*xgk(jtwm1)
    xval(j+10) = centr-absc
    xval(j+15) = centr+absc
  enddo
  !::Call the function with all points in one go
  xval(21) = centr
  fval = f_ptr( xval, dat, 21 )
  
  fc = fval(21)
  resk = wgk(11)*fc
  resabs = abs(resk)

  
  do j=1,5
    fv1(2*j) = fval(j)
    fv2(2*j) = fval(j+5)
    fsum = fval(j)+fval(j+5)
    resg = resg+wg(j)*fsum
    resk = resk+wgk(2*j)*fsum
    resabs = resabs+wgk(2*j)*(abs(fval(j))+abs(fval(j+5)))
    
    fv1(2*j-1) = fval(j+10)
    fv2(2*j-1) = fval(j+15)
    fsum = fval(j+10)+fval(j+15)
    resk = resk+wgk(2*j-1)*fsum
    resabs = resabs+wgk(2*j-1)*(abs(fval(j+10))+abs(fval(j+15)))
  enddo

  reskh = resk*5.0E-01
  resasc = wgk(11)*abs(fc-reskh)

  do j = 1, 10
    resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
  end do

  result = resk*hlgth
  resabs = resabs*dhlgth
  resasc = resasc*dhlgth
  abserr = abs((resk-resg)*hlgth)

  if ( resasc /= 0.0E+00.and.abserr /= 0.0E+00) then
    abserr = resasc*min ( 1.0E+00,(2.0E+02*abserr/resasc)**1.5E+00)
  end if

  if ( resabs > tiny ( resabs ) /(5.0E+01* epsilon ( resabs ) )) then
    abserr = max (( epsilon ( resabs ) *5.0E+01)*resabs,abserr)
  end if

  return
end subroutine qk21_x_vec

subroutine qk21_y ( f_ptr, dat, a, b, result, abserr, resabs, resasc )

!*****************************************************************************80
!
!! QK21 carries out a 21 point Gauss-Kronrod quadrature rule.
!
!  Discussion:
!
!    This routine approximates
!      I = integral ( A <= X <= B ) F(X) dx
!    with an error estimate, and
!      J = integral ( A <= X <= B ) | F(X) | dx
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real ( kind = 8 ) F, the name of the function routine, of the form
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!    which evaluates the integrand function.
!
!    Input, real ( kind = 8 ) A, B, the limits of integration.
!
!    Output, real ( kind = 8 ) RESULT, the estimated value of the integral.
!    RESULT is computed by applying the 21-point Kronrod rule (resk) 
!    obtained by optimal addition of abscissae to the 10-point Gauss 
!    rule (resg).
!
!    Output, real ( kind = 8 ) ABSERR, an estimate of | I - RESULT |.
!
!    Output, real ( kind = 8 ) RESABS, approximation to the integral of the absolute
!    value of F.
!
!    Output, real ( kind = 8 ) RESASC, approximation to the integral | F-I/(B-A) | 
!    over [A,B].
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) absc
  real ( kind = 8 ) abserr
  real ( kind = 8 ) b
  real ( kind = 8 ) centr
  real ( kind = 8 ) dhlgth
  !real ( kind = 8 ), external :: f
  procedure (f_int_dat),intent(in), pointer :: f_ptr 
  real ( kind = 8 ) fc
  real ( kind = 8 ) fsum
  real ( kind = 8 ) fval1
  real ( kind = 8 ) fval2
  real ( kind = 8 ) fv1(10)
  real ( kind = 8 ) fv2(10)
  real ( kind = 8 ) hlgth
  integer ( kind = 8 ) j
  integer ( kind = 8 ) jtw
  integer ( kind = 8 ) jtwm1
  real ( kind = 8 ) resabs
  real ( kind = 8 ) resasc
  real ( kind = 8 ) resg
  real ( kind = 8 ) resk
  real ( kind = 8 ) reskh
  real ( kind = 8 ) result
  real ( kind = 8 ) wg(5)
  real ( kind = 8 ) wgk(11)
  real ( kind = 8 ) xgk(11)
  class(dataCollectionBase), target :: dat
  
 ! f_ptr => f
!
!           the abscissae and weights are given for the interval (-1,1).
!           because of symmetry only the positive abscissae and their
!           corresponding weights are given.
!
!           xgk    - abscissae of the 21-point Kronrod rule
!                    xgk(2), xgk(4), ...  abscissae of the 10-point
!                    Gauss rule
!                    xgk(1), xgk(3), ...  abscissae which are optimally
!                    added to the 10-point Gauss rule
!
!           wgk    - weights of the 21-point Kronrod rule
!
!           wg     - weights of the 10-point Gauss rule
!
  data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8), &
    xgk(9),xgk(10),xgk(11)/ &
       9.956571630258081E-01,     9.739065285171717E-01, &
       9.301574913557082E-01,     8.650633666889845E-01, &
       7.808177265864169E-01,     6.794095682990244E-01, &
       5.627571346686047E-01,     4.333953941292472E-01, &
       2.943928627014602E-01,     1.488743389816312E-01, &
       0.000000000000000E+00/
!
  data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8), &
    wgk(9),wgk(10),wgk(11)/ &
       1.169463886737187E-02,     3.255816230796473E-02, &
       5.475589657435200E-02,     7.503967481091995E-02, &
       9.312545458369761E-02,     1.093871588022976E-01, &
       1.234919762620659E-01,     1.347092173114733E-01, &
       1.427759385770601E-01,     1.477391049013385E-01, &
       1.494455540029169E-01/
!
  data wg(1),wg(2),wg(3),wg(4),wg(5)/ &
       6.667134430868814E-02,     1.494513491505806E-01, &
       2.190863625159820E-01,     2.692667193099964E-01, &
       2.955242247147529E-01/
!
!
!           list of major variables
!
!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc   - abscissa
!           fval*  - function value
!           resg   - result of the 10-point Gauss formula
!           resk   - result of the 21-point Kronrod formula
!           reskh  - approximation to the mean value of f over (a,b),
!                    i.e. to i/(b-a)
!
  centr = 5.0E-01*(a+b)
  hlgth = 5.0E-01*(b-a)
  dhlgth = abs(hlgth)
!
!  Compute the 21-point Kronrod approximation to the
!  integral, and estimate the absolute error.
!
  resg = 0.0E+00
  fc = f_ptr(centr, dat)
  resk = wgk(11)*fc
  resabs = abs(resk)

  do j = 1, 5
    jtw = 2*j
    absc = hlgth*xgk(jtw)
    fval1 = f_ptr(centr-absc,dat)
    fval2 = f_ptr(centr+absc,dat)
    fv1(jtw) = fval1
    fv2(jtw) = fval2
    fsum = fval1+fval2
    resg = resg+wg(j)*fsum
    resk = resk+wgk(jtw)*fsum
    resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
  end do

  do j = 1, 5
    jtwm1 = 2*j-1
    absc = hlgth*xgk(jtwm1)
    fval1 = f_ptr(centr-absc,dat)
    fval2 = f_ptr(centr+absc,dat)
    fv1(jtwm1) = fval1
    fv2(jtwm1) = fval2
    fsum = fval1+fval2
    resk = resk+wgk(jtwm1)*fsum
    resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
  end do

  reskh = resk*5.0E-01
  resasc = wgk(11)*abs(fc-reskh)

  do j = 1, 10
    resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
  end do

  result = resk*hlgth
  resabs = resabs*dhlgth
  resasc = resasc*dhlgth
  abserr = abs((resk-resg)*hlgth)

  if ( resasc /= 0.0E+00.and.abserr /= 0.0E+00) then
    abserr = resasc*min ( 1.0E+00,(2.0E+02*abserr/resasc)**1.5E+00)
  end if

  if ( resabs > tiny ( resabs ) /(5.0E+01* epsilon ( resabs ) )) then
    abserr = max (( epsilon ( resabs ) *5.0E+01)*resabs,abserr)
  end if

  return
end subroutine qk21_y

subroutine qk21_y_vec ( f_ptr, dat, a, b, result, abserr, resabs, resasc )

!*****************************************************************************80
!
!! QK21 carries out a 21 point Gauss-Kronrod quadrature rule.
!
!! Modified version by Kaspar K. Nielsen, 25-1-2018. The routine has been vectorized so that the function pointer is called with 20 points in one go
!  Discussion:
!
!    This routine approximates
!      I = integral ( A <= X <= B ) F(X) dx
!    with an error estimate, and
!      J = integral ( A <= X <= B ) | F(X) | dx
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real ( kind = 8 ) F, the name of the function routine, of the form
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!    which evaluates the integrand function.
!
!    Input, real ( kind = 8 ) A, B, the limits of integration.
!
!    Output, real ( kind = 8 ) RESULT, the estimated value of the integral.
!    RESULT is computed by applying the 21-point Kronrod rule (resk) 
!    obtained by optimal addition of abscissae to the 10-point Gauss 
!    rule (resg).
!
!    Output, real ( kind = 8 ) ABSERR, an estimate of | I - RESULT |.
!
!    Output, real ( kind = 8 ) RESABS, approximation to the integral of the absolute
!    value of F.
!
!    Output, real ( kind = 8 ) RESASC, approximation to the integral | F-I/(B-A) | 
!    over [A,B].
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) absc
  real ( kind = 8 ) abserr
  real ( kind = 8 ) b
  real ( kind = 8 ) centr
  real ( kind = 8 ) dhlgth
  !real ( kind = 8 ), external :: f
  procedure (f_int_dat_vec),intent(in), pointer :: f_ptr 
  real ( kind = 8 ) fc
  real ( kind = 8 ) fsum
  real ( kind = 8 ) fval1
  real ( kind = 8 ) fval2
  real ( kind = 8 ) fv1(10)
  real ( kind = 8 ) fv2(10)
  real ( kind = 8 ) fval(21)
  real ( kind = 8 ) xval(21)
  real ( kind = 8 ) hlgth
  integer ( kind = 8 ) j
  integer ( kind = 8 ) jtw
  integer ( kind = 8 ) jtwm1
  real ( kind = 8 ) resabs
  real ( kind = 8 ) resasc
  real ( kind = 8 ) resg
  real ( kind = 8 ) resk
  real ( kind = 8 ) reskh
  real ( kind = 8 ) result
  real ( kind = 8 ) wg(5)
  real ( kind = 8 ) wgk(11)
  real ( kind = 8 ) xgk(11)
  class(dataCollectionBase), target :: dat
  
 ! f_ptr => f
!
!           the abscissae and weights are given for the interval (-1,1).
!           because of symmetry only the positive abscissae and their
!           corresponding weights are given.
!
!           xgk    - abscissae of the 21-point Kronrod rule
!                    xgk(2), xgk(4), ...  abscissae of the 10-point
!                    Gauss rule
!                    xgk(1), xgk(3), ...  abscissae which are optimally
!                    added to the 10-point Gauss rule
!
!           wgk    - weights of the 21-point Kronrod rule
!
!           wg     - weights of the 10-point Gauss rule
!
  data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8), &
    xgk(9),xgk(10),xgk(11)/ &
       9.956571630258081E-01,     9.739065285171717E-01, &
       9.301574913557082E-01,     8.650633666889845E-01, &
       7.808177265864169E-01,     6.794095682990244E-01, &
       5.627571346686047E-01,     4.333953941292472E-01, &
       2.943928627014602E-01,     1.488743389816312E-01, &
       0.000000000000000E+00/
!
  data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8), &
    wgk(9),wgk(10),wgk(11)/ &
       1.169463886737187E-02,     3.255816230796473E-02, &
       5.475589657435200E-02,     7.503967481091995E-02, &
       9.312545458369761E-02,     1.093871588022976E-01, &
       1.234919762620659E-01,     1.347092173114733E-01, &
       1.427759385770601E-01,     1.477391049013385E-01, &
       1.494455540029169E-01/
!
  data wg(1),wg(2),wg(3),wg(4),wg(5)/ &
       6.667134430868814E-02,     1.494513491505806E-01, &
       2.190863625159820E-01,     2.692667193099964E-01, &
       2.955242247147529E-01/
!
!
!           list of major variables
!
!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc   - abscissa
!           fval*  - function value
!           resg   - result of the 10-point Gauss formula
!           resk   - result of the 21-point Kronrod formula
!           reskh  - approximation to the mean value of f over (a,b),
!                    i.e. to i/(b-a)
!
  centr = 5.0E-01*(a+b)
  hlgth = 5.0E-01*(b-a)
  dhlgth = abs(hlgth)
!
!  Compute the 21-point Kronrod approximation to the
!  integral, and estimate the absolute error.
!
  resg = 0.0E+00
  
  
  do j = 1, 5
    jtw = 2*j
    absc = hlgth*xgk(jtw)    
    xval(j) = centr-absc
    xval(j+5) = centr+absc
    
    jtwm1 = 2*j-1
    absc = hlgth*xgk(jtwm1)
    xval(j+10) = centr-absc
    xval(j+15) = centr+absc
  enddo
  !::Call the function with all points in one go
  xval(21) = centr
  fval = f_ptr( xval, dat, 21 )
  
  fc = fval(21)
  resk = wgk(11)*fc
  resabs = abs(resk)

  
  do j=1,5
    fv1(2*j) = fval(j)
    fv2(2*j) = fval(j+5)
    fsum = fval(j)+fval(j+5)
    resg = resg+wg(j)*fsum
    resk = resk+wgk(2*j)*fsum
    resabs = resabs+wgk(2*j)*(abs(fval(j))+abs(fval(j+5)))
    
    fv1(2*j-1) = fval(j+10)
    fv2(2*j-1) = fval(j+15)
    fsum = fval(j+10)+fval(j+15)
    resk = resk+wgk(2*j-1)*fsum
    resabs = resabs+wgk(2*j-1)*(abs(fval(j+10))+abs(fval(j+15)))
  enddo

  reskh = resk*5.0E-01
  resasc = wgk(11)*abs(fc-reskh)

  do j = 1, 10
    resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
  end do

  result = resk*hlgth
  resabs = resabs*dhlgth
  resasc = resasc*dhlgth
  abserr = abs((resk-resg)*hlgth)

  if ( resasc /= 0.0E+00.and.abserr /= 0.0E+00) then
    abserr = resasc*min ( 1.0E+00,(2.0E+02*abserr/resasc)**1.5E+00)
  end if

  if ( resabs > tiny ( resabs ) /(5.0E+01* epsilon ( resabs ) )) then
    abserr = max (( epsilon ( resabs ) *5.0E+01)*resabs,abserr)
  end if

  return
end subroutine qk21_y_vec

subroutine qk31 ( f_ptr, dat, a, b, result, abserr, resabs, resasc )

!*****************************************************************************80
!
!! QK31 carries out a 31 point Gauss-Kronrod quadrature rule.
!
!  Discussion:
!
!    This routine approximates
!      I = integral ( A <= X <= B ) F(X) dx
!    with an error estimate, and
!      J = integral ( A <= X <= B ) | F(X) | dx
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real ( kind = 8 ) F, the name of the function routine, of the form
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!    which evaluates the integrand function.
!
!    Input, real ( kind = 8 ) A, B, the limits of integration.
!
!    Output, real ( kind = 8 ) RESULT, the estimated value of the integral.
!                       result is computed by applying the 31-point
!                       Gauss-Kronrod rule (resk), obtained by optimal
!                       addition of abscissae to the 15-point Gauss
!                       rule (resg).
!
!    Output, real ( kind = 8 ) ABSERR, an estimate of | I - RESULT |.
!
!    Output, real ( kind = 8 ) RESABS, approximation to the integral of the absolute
!    value of F.
!
!    Output, real ( kind = 8 ) RESASC, approximation to the integral | F-I/(B-A) | 
!    over [A,B].
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) absc
  real ( kind = 8 ) abserr
  real ( kind = 8 ) b
  real ( kind = 8 ) centr
  real ( kind = 8 ) dhlgth
  !real ( kind = 8 ), external :: f
  procedure (f_int_dat), intent(in), pointer :: f_ptr 
  real ( kind = 8 ) fc
  real ( kind = 8 ) fsum
  real ( kind = 8 ) fval1
  real ( kind = 8 ) fval2
  real ( kind = 8 ) fv1(15)
  real ( kind = 8 ) fv2(15)
  real ( kind = 8 ) hlgth
  integer ( kind = 8 ) j
  integer ( kind = 8 ) jtw
  integer ( kind = 8 ) jtwm1
  real ( kind = 8 ) resabs
  real ( kind = 8 ) resasc
  real ( kind = 8 ) resg
  real ( kind = 8 ) resk
  real ( kind = 8 ) reskh
  real ( kind = 8 ) result
  real ( kind = 8 ) wg(8)
  real ( kind = 8 ) wgk(16)
  real ( kind = 8 ) xgk(16)
  class( dataCollectionBase), target :: dat
!
!           the abscissae and weights are given for the interval (-1,1).
!           because of symmetry only the positive abscissae and their
!           corresponding weights are given.
!
!           xgk    - abscissae of the 31-point Kronrod rule
!                    xgk(2), xgk(4), ...  abscissae of the 15-point
!                    Gauss rule
!                    xgk(1), xgk(3), ...  abscissae which are optimally
!                    added to the 15-point Gauss rule
!
!           wgk    - weights of the 31-point Kronrod rule
!
!           wg     - weights of the 15-point Gauss rule
!
  data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8), &
    xgk(9),xgk(10),xgk(11),xgk(12),xgk(13),xgk(14),xgk(15),xgk(16)/ &
       9.980022986933971E-01,   9.879925180204854E-01, &
       9.677390756791391E-01,   9.372733924007059E-01, &
       8.972645323440819E-01,   8.482065834104272E-01, &
       7.904185014424659E-01,   7.244177313601700E-01, &
       6.509967412974170E-01,   5.709721726085388E-01, &
       4.850818636402397E-01,   3.941513470775634E-01, &
       2.991800071531688E-01,   2.011940939974345E-01, &
       1.011420669187175E-01,   0.0E+00               /
  data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8), &
    wgk(9),wgk(10),wgk(11),wgk(12),wgk(13),wgk(14),wgk(15),wgk(16)/ &
       5.377479872923349E-03,   1.500794732931612E-02, &
       2.546084732671532E-02,   3.534636079137585E-02, &
       4.458975132476488E-02,   5.348152469092809E-02, &
       6.200956780067064E-02,   6.985412131872826E-02, &
       7.684968075772038E-02,   8.308050282313302E-02, &
       8.856444305621177E-02,   9.312659817082532E-02, &
       9.664272698362368E-02,   9.917359872179196E-02, &
       1.007698455238756E-01,   1.013300070147915E-01/
  data wg(1),wg(2),wg(3),wg(4),wg(5),wg(6),wg(7),wg(8)/ &
       3.075324199611727E-02,   7.036604748810812E-02, &
       1.071592204671719E-01,   1.395706779261543E-01, &
       1.662692058169939E-01,   1.861610000155622E-01, &
       1.984314853271116E-01,   2.025782419255613E-01/
!
!
!           list of major variables
!
!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc   - abscissa
!           fval*  - function value
!           resg   - result of the 15-point Gauss formula
!           resk   - result of the 31-point Kronrod formula
!           reskh  - approximation to the mean value of f over (a,b),
!                    i.e. to i/(b-a)
!
  centr = 5.0E-01*(a+b)
  hlgth = 5.0E-01*(b-a)
  dhlgth = abs(hlgth)
  
!
!  Compute the 31-point Kronrod approximation to the integral,
!  and estimate the absolute error.
!
  fc = f_ptr(centr,dat)
  resg = wg(8)*fc
  resk = wgk(16)*fc
  resabs = abs(resk)

  do j = 1, 7
    jtw = j*2
    absc = hlgth*xgk(jtw)
    fval1 = f_ptr(centr-absc,dat)
    fval2 = f_ptr(centr+absc,dat)
    fv1(jtw) = fval1
    fv2(jtw) = fval2
    fsum = fval1+fval2
    resg = resg+wg(j)*fsum
    resk = resk+wgk(jtw)*fsum
    resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
  end do

  do j = 1, 8
    jtwm1 = j*2-1
    absc = hlgth*xgk(jtwm1)
    fval1 = f_ptr(centr-absc,dat)
    fval2 = f_ptr(centr+absc,dat)
    fv1(jtwm1) = fval1
    fv2(jtwm1) = fval2
    fsum = fval1+fval2
    resk = resk+wgk(jtwm1)*fsum
    resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
  end do

  reskh = resk*5.0E-01
  resasc = wgk(16)*abs(fc-reskh)

  do j = 1, 15
    resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
 end do

  result = resk*hlgth
  resabs = resabs*dhlgth
  resasc = resasc*dhlgth
  abserr = abs((resk-resg)*hlgth)

  if ( resasc /= 0.0E+00.and.abserr /= 0.0E+00) &
    abserr = resasc*min ( 1.0E+00,(2.0E+02*abserr/resasc)**1.5E+00)

  if ( resabs > tiny ( resabs ) /(5.0E+01* epsilon ( resabs ) )) then
    abserr = max (( epsilon ( resabs ) *5.0E+01)*resabs,abserr)
  end if

  return
end subroutine
subroutine qk41 ( f_ptr, dat, a, b, result, abserr, resabs, resasc )

!*****************************************************************************80
!
!! QK41 carries out a 41 point Gauss-Kronrod quadrature rule.
!
!  Discussion:
!
!    This routine approximates
!      I = integral ( A <= X <= B ) F(X) dx
!    with an error estimate, and
!      J = integral ( A <= X <= B ) | F(X) | dx
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real ( kind = 8 ) F, the name of the function routine, of the form
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!    which evaluates the integrand function.
!
!    Input, real ( kind = 8 ) A, B, the limits of integration.
!
!    Output, real ( kind = 8 ) RESULT, the estimated value of the integral.
!                       result is computed by applying the 41-point
!                       Gauss-Kronrod rule (resk) obtained by optimal
!                       addition of abscissae to the 20-point Gauss
!                       rule (resg).
!
!    Output, real ( kind = 8 ) ABSERR, an estimate of | I - RESULT |.
!
!    Output, real ( kind = 8 ) RESABS, approximation to the integral of the absolute
!    value of F.
!
!    Output, real ( kind = 8 ) RESASC, approximation to the integral | F-I/(B-A) | 
!    over [A,B].
!
!  Local Parameters:
!
!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc   - abscissa
!           fval*  - function value
!           resg   - result of the 20-point Gauss formula
!           resk   - result of the 41-point Kronrod formula
!           reskh  - approximation to mean value of f over (a,b), i.e.
!                    to i/(b-a)
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) absc
  real ( kind = 8 ) abserr
  real ( kind = 8 ) b
  real ( kind = 8 ) centr
  real ( kind = 8 ) dhlgth
  !real ( kind = 8 ), external :: f
  procedure (f_int_dat), intent(in), pointer :: f_ptr 
  real ( kind = 8 ) fc
  real ( kind = 8 ) fsum
  real ( kind = 8 ) fval1
  real ( kind = 8 ) fval2
  real ( kind = 8 ) fv1(20)
  real ( kind = 8 ) fv2(20)
  real ( kind = 8 ) hlgth
  integer ( kind = 8 ) j
  integer ( kind = 8 ) jtw
  integer ( kind = 8 ) jtwm1
  real ( kind = 8 ) resabs
  real ( kind = 8 ) resasc
  real ( kind = 8 ) resg
  real ( kind = 8 ) resk
  real ( kind = 8 ) reskh
  real ( kind = 8 ) result
  real ( kind = 8 ) wg(10)
  real ( kind = 8 ) wgk(21)
  real ( kind = 8 ) xgk(21)
  class( dataCollectionBase), target :: dat
!
!           the abscissae and weights are given for the interval (-1,1).
!           because of symmetry only the positive abscissae and their
!           corresponding weights are given.
!
!           xgk    - abscissae of the 41-point Gauss-Kronrod rule
!                    xgk(2), xgk(4), ...  abscissae of the 20-point
!                    Gauss rule
!                    xgk(1), xgk(3), ...  abscissae which are optimally
!                    added to the 20-point Gauss rule
!
!           wgk    - weights of the 41-point Gauss-Kronrod rule
!
!           wg     - weights of the 20-point Gauss rule
!
  data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8), &
    xgk(9),xgk(10),xgk(11),xgk(12),xgk(13),xgk(14),xgk(15),xgk(16), &
    xgk(17),xgk(18),xgk(19),xgk(20),xgk(21)/ &
       9.988590315882777E-01,   9.931285991850949E-01, &
       9.815078774502503E-01,   9.639719272779138E-01, &
       9.408226338317548E-01,   9.122344282513259E-01, &
       8.782768112522820E-01,   8.391169718222188E-01, &
       7.950414288375512E-01,   7.463319064601508E-01, &
       6.932376563347514E-01,   6.360536807265150E-01, &
       5.751404468197103E-01,   5.108670019508271E-01, &
       4.435931752387251E-01,   3.737060887154196E-01, &
       3.016278681149130E-01,   2.277858511416451E-01, &
       1.526054652409227E-01,   7.652652113349733E-02, &
       0.0E+00               /
  data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8), &
    wgk(9),wgk(10),wgk(11),wgk(12),wgk(13),wgk(14),wgk(15),wgk(16), &
    wgk(17),wgk(18),wgk(19),wgk(20),wgk(21)/ &
       3.073583718520532E-03,   8.600269855642942E-03, &
       1.462616925697125E-02,   2.038837346126652E-02, &
       2.588213360495116E-02,   3.128730677703280E-02, &
       3.660016975820080E-02,   4.166887332797369E-02, &
       4.643482186749767E-02,   5.094457392372869E-02, &
       5.519510534828599E-02,   5.911140088063957E-02, &
       6.265323755478117E-02,   6.583459713361842E-02, &
       6.864867292852162E-02,   7.105442355344407E-02, &
       7.303069033278667E-02,   7.458287540049919E-02, &
       7.570449768455667E-02,   7.637786767208074E-02, &
       7.660071191799966E-02/
  data wg(1),wg(2),wg(3),wg(4),wg(5),wg(6),wg(7),wg(8),wg(9),wg(10)/ &
       1.761400713915212E-02,   4.060142980038694E-02, &
       6.267204833410906E-02,   8.327674157670475E-02, &
       1.019301198172404E-01,   1.181945319615184E-01, &
       1.316886384491766E-01,   1.420961093183821E-01, &
       1.491729864726037E-01,   1.527533871307259E-01/
!
  centr = 5.0E-01*(a+b)
  hlgth = 5.0E-01*(b-a)
  dhlgth = abs(hlgth)
  
!
!  Compute 41-point Gauss-Kronrod approximation to the
!  the integral, and estimate the absolute error.
!
  resg = 0.0E+00
  fc = f_ptr(centr,dat)
  resk = wgk(21)*fc
  resabs = abs(resk)

  do j = 1, 10
    jtw = j*2
    absc = hlgth*xgk(jtw)
    fval1 = f_ptr(centr-absc,dat)
    fval2 = f_ptr(centr+absc,dat)
    fv1(jtw) = fval1
    fv2(jtw) = fval2
    fsum = fval1+fval2
    resg = resg+wg(j)*fsum
    resk = resk+wgk(jtw)*fsum
    resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
  end do

  do j = 1, 10
    jtwm1 = j*2-1
    absc = hlgth*xgk(jtwm1)
    fval1 = f_ptr(centr-absc,dat)
    fval2 = f_ptr(centr+absc,dat)
    fv1(jtwm1) = fval1
    fv2(jtwm1) = fval2
    fsum = fval1+fval2
    resk = resk+wgk(jtwm1)*fsum
    resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
  end do

  reskh = resk*5.0E-01
  resasc = wgk(21)*abs(fc-reskh)

  do j = 1, 20
    resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
  end do

  result = resk*hlgth
  resabs = resabs*dhlgth
  resasc = resasc*dhlgth
  abserr = abs((resk-resg)*hlgth)

  if ( resasc /= 0.0E+00.and.abserr /= 0.0E+00) &
    abserr = resasc*min ( 1.0E+00,(2.0E+02*abserr/resasc)**1.5E+00)

  if ( resabs > tiny ( resabs ) /(5.0E+01* epsilon ( resabs ) )) then
    abserr = max (( epsilon ( resabs ) *5.0E+01)*resabs,abserr)
  end if

  return
end subroutine
subroutine qk51 ( f_ptr, dat, a, b, result, abserr, resabs, resasc )

!*****************************************************************************80
!
!! QK51 carries out a 51 point Gauss-Kronrod quadrature rule.
!
!  Discussion:
!
!    This routine approximates
!      I = integral ( A <= X <= B ) F(X) dx
!    with an error estimate, and
!      J = integral ( A <= X <= B ) | F(X) | dx
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real ( kind = 8 ) F, the name of the function routine, of the form
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!    which evaluates the integrand function.
!
!    Input, real ( kind = 8 ) A, B, the limits of integration.
!
!    Output, real ( kind = 8 ) RESULT, the estimated value of the integral.
!                       result is computed by applying the 51-point
!                       Kronrod rule (resk) obtained by optimal addition
!                       of abscissae to the 25-point Gauss rule (resg).
!
!    Output, real ( kind = 8 ) ABSERR, an estimate of | I - RESULT |.
!
!    Output, real ( kind = 8 ) RESABS, approximation to the integral of the absolute
!    value of F.
!
!    Output, real ( kind = 8 ) RESASC, approximation to the integral | F-I/(B-A) | 
!    over [A,B].
!
!  Local Parameters:
!
!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc   - abscissa
!           fval*  - function value
!           resg   - result of the 25-point Gauss formula
!           resk   - result of the 51-point Kronrod formula
!           reskh  - approximation to the mean value of f over (a,b),
!                    i.e. to i/(b-a)
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) absc
  real ( kind = 8 ) abserr
  real ( kind = 8 ) b
  real ( kind = 8 ) centr
  real ( kind = 8 ) dhlgth
  !real ( kind = 8 ), external :: f
  procedure (f_int_dat), intent(in), pointer :: f_ptr 
  real ( kind = 8 ) fc
  real ( kind = 8 ) fsum
  real ( kind = 8 ) fval1
  real ( kind = 8 ) fval2
  real ( kind = 8 ) fv1(25)
  real ( kind = 8 ) fv2(25)
  real ( kind = 8 ) hlgth
  integer ( kind = 8 ) j
  integer ( kind = 8 ) jtw
  integer ( kind = 8 ) jtwm1
  real ( kind = 8 ) resabs
  real ( kind = 8 ) resasc
  real ( kind = 8 ) resg
  real ( kind = 8 ) resk
  real ( kind = 8 ) reskh
  real ( kind = 8 ) result
  real ( kind = 8 ) wg(13)
  real ( kind = 8 ) wgk(26)
  real ( kind = 8 ) xgk(26)
  class( dataCollectionBase), target :: dat
!
!           the abscissae and weights are given for the interval (-1,1).
!           because of symmetry only the positive abscissae and their
!           corresponding weights are given.
!
!           xgk    - abscissae of the 51-point Kronrod rule
!                    xgk(2), xgk(4), ...  abscissae of the 25-point
!                    Gauss rule
!                    xgk(1), xgk(3), ...  abscissae which are optimally
!                    added to the 25-point Gauss rule
!
!           wgk    - weights of the 51-point Kronrod rule
!
!           wg     - weights of the 25-point Gauss rule
!
  data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8), &
    xgk(9),xgk(10),xgk(11),xgk(12),xgk(13),xgk(14)/ &
       9.992621049926098E-01,   9.955569697904981E-01, &
       9.880357945340772E-01,   9.766639214595175E-01, &
       9.616149864258425E-01,   9.429745712289743E-01, &
       9.207471152817016E-01,   8.949919978782754E-01, &
       8.658470652932756E-01,   8.334426287608340E-01, &
       7.978737979985001E-01,   7.592592630373576E-01, &
       7.177664068130844E-01,   6.735663684734684E-01/
   data xgk(15),xgk(16),xgk(17),xgk(18),xgk(19),xgk(20),xgk(21), &
    xgk(22),xgk(23),xgk(24),xgk(25),xgk(26)/ &
       6.268100990103174E-01,   5.776629302412230E-01, &
       5.263252843347192E-01,   4.730027314457150E-01, &
       4.178853821930377E-01,   3.611723058093878E-01, &
       3.030895389311078E-01,   2.438668837209884E-01, &
       1.837189394210489E-01,   1.228646926107104E-01, &
       6.154448300568508E-02,   0.0E+00               /
  data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8), &
    wgk(9),wgk(10),wgk(11),wgk(12),wgk(13),wgk(14)/ &
       1.987383892330316E-03,   5.561932135356714E-03, &
       9.473973386174152E-03,   1.323622919557167E-02, &
       1.684781770912830E-02,   2.043537114588284E-02, &
       2.400994560695322E-02,   2.747531758785174E-02, &
       3.079230016738749E-02,   3.400213027432934E-02, &
       3.711627148341554E-02,   4.008382550403238E-02, &
       4.287284502017005E-02,   4.550291304992179E-02/
   data wgk(15),wgk(16),wgk(17),wgk(18),wgk(19),wgk(20),wgk(21), &
    wgk(22),wgk(23),wgk(24),wgk(25),wgk(26)/ &
       4.798253713883671E-02,   5.027767908071567E-02, &
       5.236288580640748E-02,   5.425112988854549E-02, &
       5.595081122041232E-02,   5.743711636156783E-02, &
       5.868968002239421E-02,   5.972034032417406E-02, &
       6.053945537604586E-02,   6.112850971705305E-02, &
       6.147118987142532E-02,   6.158081806783294E-02/
  data wg(1),wg(2),wg(3),wg(4),wg(5),wg(6),wg(7),wg(8),wg(9),wg(10), &
    wg(11),wg(12),wg(13)/ &
       1.139379850102629E-02,   2.635498661503214E-02, &
       4.093915670130631E-02,   5.490469597583519E-02, &
       6.803833381235692E-02,   8.014070033500102E-02, &
       9.102826198296365E-02,   1.005359490670506E-01, &
       1.085196244742637E-01,   1.148582591457116E-01, &
       1.194557635357848E-01,   1.222424429903100E-01, &
       1.231760537267155E-01/
!
  centr = 5.0E-01*(a+b)
  hlgth = 5.0E-01*(b-a)
  dhlgth = abs(hlgth)
  
!
!  Compute the 51-point Kronrod approximation to the integral,
!  and estimate the absolute error.
!
  fc = f_ptr(centr,dat)
  resg = wg(13)*fc
  resk = wgk(26)*fc
  resabs = abs(resk)

  do j = 1, 12
    jtw = j*2
    absc = hlgth*xgk(jtw)
    fval1 = f_ptr(centr-absc,dat)
    fval2 = f_ptr(centr+absc,dat)
    fv1(jtw) = fval1
    fv2(jtw) = fval2
    fsum = fval1+fval2
    resg = resg+wg(j)*fsum
    resk = resk+wgk(jtw)*fsum
    resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
  end do

  do j = 1, 13
    jtwm1 = j*2-1
    absc = hlgth*xgk(jtwm1)
    fval1 = f_ptr(centr-absc,dat)
    fval2 = f_ptr(centr+absc,dat)
    fv1(jtwm1) = fval1
    fv2(jtwm1) = fval2
    fsum = fval1+fval2
    resk = resk+wgk(jtwm1)*fsum
    resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
  end do

  reskh = resk*5.0E-01
  resasc = wgk(26)*abs(fc-reskh)

  do j = 1, 25
    resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
  end do

  result = resk*hlgth
  resabs = resabs*dhlgth
  resasc = resasc*dhlgth
  abserr = abs((resk-resg)*hlgth)

  if ( resasc /= 0.0E+00.and.abserr /= 0.0E+00) then
    abserr = resasc*min ( 1.0E+00,(2.0E+02*abserr/resasc)**1.5E+00)
  end if

  if ( resabs > tiny ( resabs ) / (5.0E+01* epsilon ( resabs ) ) ) then
    abserr = max (( epsilon ( resabs ) *5.0E+01)*resabs,abserr)
  end if

  return
end subroutine
subroutine qk61 ( f_ptr, dat, a, b, result, abserr, resabs, resasc ) 

!*****************************************************************************80
!
!! QK61 carries out a 61 point Gauss-Kronrod quadrature rule.
!
!  Discussion:
!
!    This routine approximates
!      I = integral ( A <= X <= B ) F(X) dx
!    with an error estimate, and
!      J = integral ( A <= X <= B ) | F(X) | dx
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real ( kind = 8 ) F, the name of the function routine, of the form
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!    which evaluates the integrand function.
!
!    Input, real ( kind = 8 ) A, B, the limits of integration.
!
!    Output, real ( kind = 8 ) RESULT, the estimated value of the integral.
!                    result is computed by applying the 61-point
!                    Kronrod rule (resk) obtained by optimal addition of
!                    abscissae to the 30-point Gauss rule (resg).
!
!    Output, real ( kind = 8 ) ABSERR, an estimate of | I - RESULT |.
!
!    Output, real ( kind = 8 ) RESABS, approximation to the integral of the absolute
!    value of F.
!
!    Output, real ( kind = 8 ) RESASC, approximation to the integral | F-I/(B-A) | 
!    over [A,B].
!
!  Local Parameters:
!
!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc   - abscissa
!           fval*  - function value
!           resg   - result of the 30-point Gauss rule
!           resk   - result of the 61-point Kronrod rule
!           reskh  - approximation to the mean value of f
!                    over (a,b), i.e. to i/(b-a)
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) absc
  real ( kind = 8 ) abserr
  real ( kind = 8 ) b
  real ( kind = 8 ) centr
  real ( kind = 8 ) dhlgth
  !real ( kind = 8 ), external :: f
  procedure (f_int_dat), intent(in), pointer :: f_ptr 
  real ( kind = 8 ) fc
  real ( kind = 8 ) fsum
  real ( kind = 8 ) fval1
  real ( kind = 8 ) fval2
  real ( kind = 8 ) fv1(30)
  real ( kind = 8 ) fv2(30)
  real ( kind = 8 ) hlgth
  integer ( kind = 8 ) j
  integer ( kind = 8 ) jtw
  integer ( kind = 8 ) jtwm1
  real ( kind = 8 ) resabs
  real ( kind = 8 ) resasc
  real ( kind = 8 ) resg
  real ( kind = 8 ) resk
  real ( kind = 8 ) reskh
  real ( kind = 8 ) result
  real ( kind = 8 ) wg(15)
  real ( kind = 8 ) wgk(31)
  real ( kind = 8 ) xgk(31)
  class( dataCollectionBase), target :: dat
!
!           the abscissae and weights are given for the
!           interval (-1,1). because of symmetry only the positive
!           abscissae and their corresponding weights are given.
!
!           xgk   - abscissae of the 61-point Kronrod rule
!                   xgk(2), xgk(4)  ... abscissae of the 30-point
!                   Gauss rule
!                   xgk(1), xgk(3)  ... optimally added abscissae
!                   to the 30-point Gauss rule
!
!           wgk   - weights of the 61-point Kronrod rule
!
!           wg    - weigths of the 30-point Gauss rule
!
  data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8), &
     xgk(9),xgk(10)/ &
       9.994844100504906E-01,     9.968934840746495E-01, &
       9.916309968704046E-01,     9.836681232797472E-01, &
       9.731163225011263E-01,     9.600218649683075E-01, &
       9.443744447485600E-01,     9.262000474292743E-01, &
       9.055733076999078E-01,     8.825605357920527E-01/
  data xgk(11),xgk(12),xgk(13),xgk(14),xgk(15),xgk(16),xgk(17), &
    xgk(18),xgk(19),xgk(20)/ &
       8.572052335460611E-01,     8.295657623827684E-01, &
       7.997278358218391E-01,     7.677774321048262E-01, &
       7.337900624532268E-01,     6.978504947933158E-01, &
       6.600610641266270E-01,     6.205261829892429E-01, &
       5.793452358263617E-01,     5.366241481420199E-01/
  data xgk(21),xgk(22),xgk(23),xgk(24),xgk(25),xgk(26),xgk(27), &
    xgk(28),xgk(29),xgk(30),xgk(31)/ &
       4.924804678617786E-01,     4.470337695380892E-01, &
       4.004012548303944E-01,     3.527047255308781E-01, &
       3.040732022736251E-01,     2.546369261678898E-01, &
       2.045251166823099E-01,     1.538699136085835E-01, &
       1.028069379667370E-01,     5.147184255531770E-02, &
       0.0E+00                   /
  data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8), &
    wgk(9),wgk(10)/ &
       1.389013698677008E-03,     3.890461127099884E-03, &
       6.630703915931292E-03,     9.273279659517763E-03, &
       1.182301525349634E-02,     1.436972950704580E-02, &
       1.692088918905327E-02,     1.941414119394238E-02, &
       2.182803582160919E-02,     2.419116207808060E-02/
  data wgk(11),wgk(12),wgk(13),wgk(14),wgk(15),wgk(16),wgk(17), &
    wgk(18),wgk(19),wgk(20)/ &
       2.650995488233310E-02,     2.875404876504129E-02, &
       3.090725756238776E-02,     3.298144705748373E-02, &
       3.497933802806002E-02,     3.688236465182123E-02, &
       3.867894562472759E-02,     4.037453895153596E-02, &
       4.196981021516425E-02,     4.345253970135607E-02/
  data wgk(21),wgk(22),wgk(23),wgk(24),wgk(25),wgk(26),wgk(27), &
    wgk(28),wgk(29),wgk(30),wgk(31)/ &
       4.481480013316266E-02,     4.605923827100699E-02, &
       4.718554656929915E-02,     4.818586175708713E-02, &
       4.905543455502978E-02,     4.979568342707421E-02, &
       5.040592140278235E-02,     5.088179589874961E-02, &
       5.122154784925877E-02,     5.142612853745903E-02, &
       5.149472942945157E-02/
  data wg(1),wg(2),wg(3),wg(4),wg(5),wg(6),wg(7),wg(8)/ &
       7.968192496166606E-03,     1.846646831109096E-02, &
       2.878470788332337E-02,     3.879919256962705E-02, &
       4.840267283059405E-02,     5.749315621761907E-02, &
       6.597422988218050E-02,     7.375597473770521E-02/
  data wg(9),wg(10),wg(11),wg(12),wg(13),wg(14),wg(15)/ &
       8.075589522942022E-02,     8.689978720108298E-02, &
       9.212252223778613E-02,     9.636873717464426E-02, &
       9.959342058679527E-02,     1.017623897484055E-01, &
       1.028526528935588E-01/

  centr = 5.0E-01*(b+a)
  hlgth = 5.0E-01*(b-a)
  dhlgth = abs(hlgth)
  
!
!  Compute the 61-point Kronrod approximation to the integral,
!  and estimate the absolute error.
!
  resg = 0.0E+00
  fc = f_ptr(centr,dat)
  resk = wgk(31)*fc
  resabs = abs(resk)

  do j = 1, 15
    jtw = j*2
    absc = hlgth*xgk(jtw)
    fval1 = f_ptr(centr-absc,dat)
    fval2 = f_ptr(centr+absc,dat)
    fv1(jtw) = fval1
    fv2(jtw) = fval2
    fsum = fval1+fval2
    resg = resg+wg(j)*fsum
    resk = resk+wgk(jtw)*fsum
    resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
  end do

  do j = 1, 15
    jtwm1 = j*2-1
    absc = hlgth*xgk(jtwm1)
    fval1 = f_ptr(centr-absc,dat)
    fval2 = f_ptr(centr+absc,dat)
    fv1(jtwm1) = fval1
    fv2(jtwm1) = fval2
    fsum = fval1+fval2
    resk = resk+wgk(jtwm1)*fsum
    resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
  end do

  reskh = resk * 5.0E-01
  resasc = wgk(31)*abs(fc-reskh)

  do j = 1, 30
    resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
  end do

  result = resk*hlgth
  resabs = resabs*dhlgth
  resasc = resasc*dhlgth
  abserr = abs((resk-resg)*hlgth)

  if ( resasc /= 0.0E+00 .and. abserr /= 0.0E+00) then
    abserr = resasc*min ( 1.0E+00,(2.0E+02*abserr/resasc)**1.5E+00)
  end if

  if ( resabs > tiny ( resabs ) / (5.0E+01* epsilon ( resabs ) )) then
    abserr = max ( ( epsilon ( resabs ) *5.0E+01)*resabs, abserr )
  end if

  return
end subroutine
subroutine qmomo ( alfa, beta, ri, rj, rg, rh, integr )

!*****************************************************************************80
!
!! QMOMO computes modified Chebyshev moments.
!
!  Discussion:
!
!    This routine computes modified Chebyshev moments.
!    The K-th modified Chebyshev moment is defined as the
!    integral over (-1,1) of W(X)*T(K,X), where T(K,X) is the
!    Chebyshev polynomial of degree K.
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALFA, a parameter in the weight function w(x), ALFA > -1.
!
!    Input, real ( kind = 8 ) BETA, a parameter in the weight function w(x), BETA > -1.
!
!           ri     - real
!                    vector of dimension 25
!                    ri(k) is the integral over (-1,1) of
!                    (1+x)**alfa*t(k-1,x), k = 1, ..., 25.
!
!           rj     - real
!                    vector of dimension 25
!                    rj(k) is the integral over (-1,1) of
!                    (1-x)**beta*t(k-1,x), k = 1, ..., 25.
!
!           rg     - real
!                    vector of dimension 25
!                    rg(k) is the integral over (-1,1) of
!                    (1+x)**alfa*log((1+x)/2)*t(k-1,x), k = 1, ...,25.
!
!           rh     - real
!                    vector of dimension 25
!                    rh(k) is the integral over (-1,1) of
!                    (1-x)**beta*log((1-x)/2)*t(k-1,x), k = 1, ..., 25.
!
!           integr - integer ( kind = 8 )
!                    input parameter indicating the modified moments
!                    to be computed
!                    integr = 1 compute ri, rj
!                           = 2 compute ri, rj, rg
!                           = 3 compute ri, rj, rh
!                           = 4 compute ri, rj, rg, rh
!
  implicit none

  real ( kind = 8 ) alfa
  real ( kind = 8 ) alfp1
  real ( kind = 8 ) alfp2
  real ( kind = 8 ) an
  real ( kind = 8 ) anm1
  real ( kind = 8 ) beta
  real ( kind = 8 ) betp1
  real ( kind = 8 ) betp2
  integer ( kind = 8 ) i
  integer ( kind = 8 ) im1
  integer ( kind = 8 ) integr
  real ( kind = 8 ) ralf
  real ( kind = 8 ) rbet
  real ( kind = 8 ) rg(25)
  real ( kind = 8 ) rh(25)
  real ( kind = 8 ) ri(25)
  real ( kind = 8 ) rj(25)
!
  alfp1 = alfa+1.0E+00
  betp1 = beta+1.0E+00
  alfp2 = alfa+2.0E+00
  betp2 = beta+2.0E+00
  ralf = 2.0E+00**alfp1
  rbet = 2.0E+00**betp1
!
!  Compute RI, RJ using a forward recurrence relation.
!
  ri(1) = ralf/alfp1
  rj(1) = rbet/betp1
  ri(2) = ri(1)*alfa/alfp2
  rj(2) = rj(1)*beta/betp2
  an = 2.0E+00
  anm1 = 1.0E+00

  do i = 3, 25
    ri(i) = -(ralf+an*(an-alfp2)*ri(i-1))/(anm1*(an+alfp1))
    rj(i) = -(rbet+an*(an-betp2)*rj(i-1))/(anm1*(an+betp1))
    anm1 = an
    an = an+1.0E+00
  end do

  if ( integr == 1 ) go to 70
  if ( integr == 3 ) go to 40
!
!  Compute RG using a forward recurrence relation.
!
  rg(1) = -ri(1)/alfp1
  rg(2) = -(ralf+ralf)/(alfp2*alfp2)-rg(1)
  an = 2.0E+00
  anm1 = 1.0E+00
  im1 = 2

  do i = 3, 25
    rg(i) = -(an*(an-alfp2)*rg(im1)-an*ri(im1)+anm1*ri(i))/ &
    (anm1*(an+alfp1))
    anm1 = an
    an = an+1.0E+00
    im1 = i
  end do

  if ( integr == 2 ) go to 70
!
!  Compute RH using a forward recurrence relation.
!
40 continue

  rh(1) = -rj(1) / betp1
  rh(2) = -(rbet+rbet)/(betp2*betp2)-rh(1)
  an = 2.0E+00
  anm1 = 1.0E+00
  im1 = 2

  do i = 3, 25
    rh(i) = -(an*(an-betp2)*rh(im1)-an*rj(im1)+ &
    anm1*rj(i))/(anm1*(an+betp1))
    anm1 = an
    an = an+1.0E+00
    im1 = i
  end do

  do i = 2, 25, 2
    rh(i) = -rh(i)
  end do

   70 continue

  do i = 2, 25, 2
    rj(i) = -rj(i)
  end do

!  90 continue

  return
end subroutine
subroutine qng ( f_ptr, dat, a, b, epsabs, epsrel, result, abserr, neval, ier )

!*****************************************************************************80
!
!! QNG estimates an integral, using non-adaptive integration.
!
!  Discussion:
!
!    The routine calculates an approximation RESULT to a definite integral   
!      I = integral of F over (A,B),
!    hopefully satisfying
!      || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
!
!    The routine is a simple non-adaptive automatic integrator, based on
!    a sequence of rules with increasing degree of algebraic
!    precision (Patterson, 1968).
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real ( kind = 8 ) F, the name of the function routine, of the form
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!    which evaluates the integrand function.
!
!    Input, real ( kind = 8 ) A, B, the limits of integration.
!
!    Input, real ( kind = 8 ) EPSABS, EPSREL, the absolute and relative accuracy requested.
!
!    Output, real ( kind = 8 ) RESULT, the estimated value of the integral.
!    RESULT is obtained by applying the 21-point Gauss-Kronrod rule (RES21)
!    obtained  by optimal addition of abscissae to the 10-point Gauss rule
!    (RES10), or by applying the 43-point rule (RES43) obtained by optimal
!    addition of abscissae to the 21-point Gauss-Kronrod rule, or by 
!    applying the 87-point rule (RES87) obtained by optimal addition of
!    abscissae to the 43-point rule.
!
!    Output, real ( kind = 8 ) ABSERR, an estimate of || I - RESULT ||.
!
!    Output, integer ( kind = 8 ) NEVAL, the number of times the integral was evaluated.
!
!           ier    - ier = 0 normal and reliable termination of the
!                            routine. it is assumed that the requested
!                            accuracy has been achieved.
!                    ier > 0 abnormal termination of the routine. it is
!                            assumed that the requested accuracy has
!                            not been achieved.
!                    ier = 1 the maximum number of steps has been
!                            executed. the integral is probably too
!                            difficult to be calculated by qng.
!                        = 6 the input is invalid, because
!                            epsabs < 0 and epsrel < 0,
!                            result, abserr and neval are set to zero.
!
!  Local Parameters:
!
!           centr  - mid point of the integration interval
!           hlgth  - half-length of the integration interval
!           fcentr - function value at mid point
!           absc   - abscissa
!           fval   - function value
!           savfun - array of function values which have already
!                    been computed
!           res10  - 10-point Gauss result
!           res21  - 21-point Kronrod result
!           res43  - 43-point result
!           res87  - 87-point result
!           resabs - approximation to the integral of abs(f)
!           resasc - approximation to the integral of abs(f-i/(b-a))
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) absc
  real ( kind = 8 ) abserr
  real ( kind = 8 ) b
  real ( kind = 8 ) centr
  real ( kind = 8 ) dhlgth
  real ( kind = 8 ) epsabs
  real ( kind = 8 ) epsrel
  !real ( kind = 8 ), external :: f
  procedure (f_int_dat), intent(in), pointer :: f_ptr 
  real ( kind = 8 ) fcentr
  real ( kind = 8 ) fval
  real ( kind = 8 ) fval1
  real ( kind = 8 ) fval2
  real ( kind = 8 ) fv1(5)
  real ( kind = 8 ) fv2(5)
  real ( kind = 8 ) fv3(5)
  real ( kind = 8 ) fv4(5)
  real ( kind = 8 ) hlgth
  integer ( kind = 8 ) ier
  integer ( kind = 8 ) ipx
  integer ( kind = 8 ) k
  integer ( kind = 8 ) l
  integer ( kind = 8 ) neval
  real ( kind = 8 ) result
  real ( kind = 8 ) res10
  real ( kind = 8 ) res21
  real ( kind = 8 ) res43
  real ( kind = 8 ) res87
  real ( kind = 8 ) resabs
  real ( kind = 8 ) resasc
  real ( kind = 8 ) reskh
  real ( kind = 8 ) savfun(21)
  real ( kind = 8 ) w10(5)
  real ( kind = 8 ) w21a(5)
  real ( kind = 8 ) w21b(6)
  real ( kind = 8 ) w43a(10)
  real ( kind = 8 ) w43b(12)
  real ( kind = 8 ) w87a(21)
  real ( kind = 8 ) w87b(23)
  real ( kind = 8 ) x1(5)
  real ( kind = 8 ) x2(5)
  real ( kind = 8 ) x3(11)
  real ( kind = 8 ) x4(22)
  class( dataCollectionBase), target :: dat
!
!           the following data statements contain the abscissae
!           and weights of the integration rules used.
!
!           x1      abscissae common to the 10-, 21-, 43- and 87-point
!                   rule
!           x2      abscissae common to the 21-, 43- and 87-point rule
!           x3      abscissae common to the 43- and 87-point rule
!           x4      abscissae of the 87-point rule
!           w10     weights of the 10-point formula
!           w21a    weights of the 21-point formula for abscissae x1
!           w21b    weights of the 21-point formula for abscissae x2
!           w43a    weights of the 43-point formula for absissae x1, x3
!           w43b    weights of the 43-point formula for abscissae x3
!           w87a    weights of the 87-point formula for abscissae x1,
!                   x2 and x3
!           w87b    weights of the 87-point formula for abscissae x4
!
  data x1(1),x1(2),x1(3),x1(4),x1(5)/ &
       9.739065285171717E-01,     8.650633666889845E-01, &
       6.794095682990244E-01,     4.333953941292472E-01, &
       1.488743389816312E-01/
  data x2(1),x2(2),x2(3),x2(4),x2(5)/ &
       9.956571630258081E-01,     9.301574913557082E-01, &
       7.808177265864169E-01,     5.627571346686047E-01, &
       2.943928627014602E-01/
  data x3(1),x3(2),x3(3),x3(4),x3(5),x3(6),x3(7),x3(8),x3(9),x3(10), &
    x3(11)/ &
       9.993333609019321E-01,     9.874334029080889E-01, &
       9.548079348142663E-01,     9.001486957483283E-01, &
       8.251983149831142E-01,     7.321483889893050E-01, &
       6.228479705377252E-01,     4.994795740710565E-01, &
       3.649016613465808E-01,     2.222549197766013E-01, &
       7.465061746138332E-02/
  data x4(1),x4(2),x4(3),x4(4),x4(5),x4(6),x4(7),x4(8),x4(9),x4(10), &
    x4(11),x4(12),x4(13),x4(14),x4(15),x4(16),x4(17),x4(18),x4(19), &
    x4(20),x4(21),x4(22)/         9.999029772627292E-01, &
       9.979898959866787E-01,     9.921754978606872E-01, &
       9.813581635727128E-01,     9.650576238583846E-01, &
       9.431676131336706E-01,     9.158064146855072E-01, &
       8.832216577713165E-01,     8.457107484624157E-01, &
       8.035576580352310E-01,     7.570057306854956E-01, &
       7.062732097873218E-01,     6.515894665011779E-01, &
       5.932233740579611E-01,     5.314936059708319E-01, &
       4.667636230420228E-01,     3.994248478592188E-01, &
       3.298748771061883E-01,     2.585035592021616E-01, &
       1.856953965683467E-01,     1.118422131799075E-01, &
       3.735212339461987E-02/
  data w10(1),w10(2),w10(3),w10(4),w10(5)/ &
       6.667134430868814E-02,     1.494513491505806E-01, &
       2.190863625159820E-01,     2.692667193099964E-01, &
       2.955242247147529E-01/
  data w21a(1),w21a(2),w21a(3),w21a(4),w21a(5)/ &
       3.255816230796473E-02,     7.503967481091995E-02, &
       1.093871588022976E-01,     1.347092173114733E-01, &
       1.477391049013385E-01/
  data w21b(1),w21b(2),w21b(3),w21b(4),w21b(5),w21b(6)/ &
       1.169463886737187E-02,     5.475589657435200E-02, &
       9.312545458369761E-02,     1.234919762620659E-01, &
       1.427759385770601E-01,     1.494455540029169E-01/
  data w43a(1),w43a(2),w43a(3),w43a(4),w43a(5),w43a(6),w43a(7), &
    w43a(8),w43a(9),w43a(10)/     1.629673428966656E-02, &
       3.752287612086950E-02,     5.469490205825544E-02, &
       6.735541460947809E-02,     7.387019963239395E-02, &
       5.768556059769796E-03,     2.737189059324884E-02, &
       4.656082691042883E-02,     6.174499520144256E-02, &
       7.138726726869340E-02/
  data w43b(1),w43b(2),w43b(3),w43b(4),w43b(5),w43b(6),w43b(7), &
    w43b(8),w43b(9),w43b(10),w43b(11),w43b(12)/ &
       1.844477640212414E-03,     1.079868958589165E-02, &
       2.189536386779543E-02,     3.259746397534569E-02, &
       4.216313793519181E-02,     5.074193960018458E-02, &
       5.837939554261925E-02,     6.474640495144589E-02, &
       6.956619791235648E-02,     7.282444147183321E-02, &
       7.450775101417512E-02,     7.472214751740301E-02/
  data w87a(1),w87a(2),w87a(3),w87a(4),w87a(5),w87a(6),w87a(7), &
    w87a(8),w87a(9),w87a(10),w87a(11),w87a(12),w87a(13),w87a(14), &
    w87a(15),w87a(16),w87a(17),w87a(18),w87a(19),w87a(20),w87a(21)/ &
       8.148377384149173E-03,     1.876143820156282E-02, &
       2.734745105005229E-02,     3.367770731163793E-02, &
       3.693509982042791E-02,     2.884872430211531E-03, &
       1.368594602271270E-02,     2.328041350288831E-02, &
       3.087249761171336E-02,     3.569363363941877E-02, &
       9.152833452022414E-04,     5.399280219300471E-03, &
       1.094767960111893E-02,     1.629873169678734E-02, &
       2.108156888920384E-02,     2.537096976925383E-02, &
       2.918969775647575E-02,     3.237320246720279E-02, &
       3.478309895036514E-02,     3.641222073135179E-02, &
       3.725387550304771E-02/
  data w87b(1),w87b(2),w87b(3),w87b(4),w87b(5),w87b(6),w87b(7), &
    w87b(8),w87b(9),w87b(10),w87b(11),w87b(12),w87b(13),w87b(14), &
    w87b(15),w87b(16),w87b(17),w87b(18),w87b(19),w87b(20),w87b(21), &
    w87b(22),w87b(23)/            2.741455637620724E-04, &
       1.807124155057943E-03,     4.096869282759165E-03, &
       6.758290051847379E-03,     9.549957672201647E-03, &
       1.232944765224485E-02,     1.501044734638895E-02, &
       1.754896798624319E-02,     1.993803778644089E-02, &
       2.219493596101229E-02,     2.433914712600081E-02, &
       2.637450541483921E-02,     2.828691078877120E-02, &
       3.005258112809270E-02,     3.164675137143993E-02, &
       3.305041341997850E-02,     3.425509970422606E-02, &
       3.526241266015668E-02,     3.607698962288870E-02, &
       3.669860449845609E-02,     3.712054926983258E-02, &
       3.733422875193504E-02,     3.736107376267902E-02/
!
!  Test on validity of parameters.
!
  result = 0.0E+00
  abserr = 0.0E+00
  neval = 0
  

  if ( epsabs < 0.0E+00 .and. epsrel < 0.0E+00 ) then
    ier = 6
    return
  end if

  hlgth = 5.0E-01 * ( b - a )
  dhlgth = abs ( hlgth )
  centr = 5.0E-01 * ( b + a )
  fcentr = f_ptr(centr,dat)
  neval = 21
  ier = 1
!
!  Compute the integral using the 10- and 21-point formula.
!
  do l = 1, 3

    if ( l == 1 ) then

      res10 = 0.0E+00
      res21 = w21b(6) * fcentr
      resabs = w21b(6) * abs(fcentr)

      do k = 1, 5
        absc = hlgth * x1(k)
        fval1 = f_ptr(centr+absc,dat)
        fval2 = f_ptr(centr-absc,dat)
        fval = fval1 + fval2
        res10 = res10 + w10(k)*fval
        res21 = res21 + w21a(k)*fval
        resabs = resabs + w21a(k)*(abs(fval1)+abs(fval2))
        savfun(k) = fval
        fv1(k) = fval1
        fv2(k) = fval2
      end do

      ipx = 5

      do k = 1, 5
        ipx = ipx + 1
        absc = hlgth * x2(k)
        fval1 = f_ptr(centr+absc,dat)
        fval2 = f_ptr(centr-absc,dat)
        fval = fval1 + fval2
        res21 = res21 + w21b(k) * fval
        resabs = resabs + w21b(k) * ( abs ( fval1 ) + abs ( fval2 ) )
        savfun(ipx) = fval
        fv3(k) = fval1
        fv4(k) = fval2
      end do
!
!  Test for convergence.
!
      result = res21 * hlgth
      resabs = resabs * dhlgth
      reskh = 5.0E-01 * res21
      resasc = w21b(6) * abs ( fcentr - reskh )

      do k = 1, 5
        resasc = resasc+w21a(k)*(abs(fv1(k)-reskh)+abs(fv2(k)-reskh)) &
                     +w21b(k)*(abs(fv3(k)-reskh)+abs(fv4(k)-reskh))
      end do

      abserr = abs ( ( res21 - res10 ) * hlgth )
      resasc = resasc * dhlgth
!
!  Compute the integral using the 43-point formula.
!
    else if ( l == 2 ) then

      res43 = w43b(12)*fcentr
      neval = 43

      do k = 1, 10
        res43 = res43 + savfun(k) * w43a(k)
      end do

      do k = 1, 11
        ipx = ipx + 1
        absc = hlgth * x3(k)
        fval = f_ptr(absc+centr,dat) + f_ptr(centr-absc,dat)
        res43 = res43 + fval * w43b(k)
        savfun(ipx) = fval
      end do
!
!  Test for convergence.
!
      result = res43 * hlgth
      abserr = abs((res43-res21)*hlgth)
!
!  Compute the integral using the 87-point formula.
!
    else if ( l == 3 ) then

      res87 = w87b(23) * fcentr
      neval = 87

      do k = 1, 21
        res87 = res87 + savfun(k) * w87a(k)
      end do

      do k = 1, 22
        absc = hlgth * x4(k)
        res87 = res87 + w87b(k) * ( f_ptr(absc+centr,dat) + f_ptr(centr-absc,dat) )
      end do

      result = res87 * hlgth
      abserr = abs ( ( res87 - res43) * hlgth )

    end if

    if ( resasc /= 0.0E+00.and.abserr /= 0.0E+00 ) then
      abserr = resasc * min ( 1.0E+00,(2.0E+02*abserr/resasc)**1.5E+00)
    end if

    if ( resabs > tiny ( resabs ) / ( 5.0E+01 * epsilon ( resabs ) ) ) then
      abserr = max (( epsilon ( resabs ) *5.0E+01) * resabs, abserr )
    end if

    if ( abserr <= max ( epsabs, epsrel*abs(result))) then
      ier = 0
    end if

    if ( ier == 0 ) then
      exit
    end if

  end do

  return
end subroutine
subroutine qsort ( limit, last, maxerr, ermax, elist, iord, nrmax )

!*****************************************************************************80
!
!! QSORT maintains the order of a list of local error estimates.
!
!  Discussion:
!
!    This routine maintains the descending ordering in the list of the 
!    local error estimates resulting from the interval subdivision process. 
!    At each call two error estimates are inserted using the sequential 
!    search top-down for the largest error estimate and bottom-up for the
!    smallest error estimate.
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, integer ( kind = 8 ) LIMIT, the maximum number of error estimates the list can
!    contain.
!
!    Input, integer ( kind = 8 ) LAST, the current number of error estimates.
!
!    Input/output, integer ( kind = 8 ) MAXERR, the index in the list of the NRMAX-th 
!    largest error.
!
!    Output, real ( kind = 8 ) ERMAX, the NRMAX-th largest error = ELIST(MAXERR).
!
!    Input, real ( kind = 8 ) ELIST(LIMIT), contains the error estimates.
!
!    Input/output, integer ( kind = 8 ) IORD(LAST).  The first K elements contain 
!    pointers to the error estimates such that ELIST(IORD(1)) through
!    ELIST(IORD(K)) form a decreasing sequence, with
!      K = LAST 
!    if 
!      LAST <= (LIMIT/2+2), 
!    and otherwise
!      K = LIMIT+1-LAST.
!
!    Input/output, integer ( kind = 8 ) NRMAX.
!
  implicit none

  integer ( kind = 8 ) last

  real ( kind = 8 ) elist(last)
  real ( kind = 8 ) ermax
  real ( kind = 8 ) errmax
  real ( kind = 8 ) errmin
  integer ( kind = 8 ) i
  integer ( kind = 8 ) ibeg
  integer ( kind = 8 ) iord(last)
  integer ( kind = 8 ) isucc
  integer ( kind = 8 ) j
  integer ( kind = 8 ) jbnd
  integer ( kind = 8 ) jupbn
  integer ( kind = 8 ) k
  integer ( kind = 8 ) limit
  integer ( kind = 8 ) maxerr
  integer ( kind = 8 ) nrmax
!
!  Check whether the list contains more than two error estimates.
!
  if ( last <= 2 ) then
    iord(1) = 1
    iord(2) = 2
    go to 90
  end if
!
!  This part of the routine is only executed if, due to a
!  difficult integrand, subdivision increased the error
!  estimate. in the normal case the insert procedure should
!  start after the nrmax-th largest error estimate.
!
  errmax = elist(maxerr)

  do i = 1, nrmax-1

    isucc = iord(nrmax-1)

    if ( errmax <= elist(isucc) ) then
      exit
    end if

    iord(nrmax) = isucc
    nrmax = nrmax-1

  end do
!
!  Compute the number of elements in the list to be maintained
!  in descending order.  This number depends on the number of
!  subdivisions still allowed.
!
  jupbn = last

  if ( (limit/2+2) < last ) then
    jupbn = limit+3-last
  end if

  errmin = elist(last)
!
!  Insert errmax by traversing the list top-down, starting
!  comparison from the element elist(iord(nrmax+1)).
!
  jbnd = jupbn-1
  ibeg = nrmax+1

  do i = ibeg, jbnd
    isucc = iord(i)
    if ( elist(isucc) <= errmax ) then
      go to 60
    end if
    iord(i-1) = isucc
  end do

  iord(jbnd) = maxerr
  iord(jupbn) = last
  go to 90
!
!  Insert errmin by traversing the list bottom-up.
!
60 continue

  iord(i-1) = maxerr
  k = jbnd

  do j = i, jbnd
    isucc = iord(k)
    if ( errmin < elist(isucc) ) then
      go to 80
    end if
    iord(k+1) = isucc
    k = k-1
  end do

  iord(i) = last
  go to 90

80 continue

  iord(k+1) = last
!
!  Set maxerr and ermax.
!
90 continue

  maxerr = iord(nrmax)
  ermax = elist(maxerr)

  return
end subroutine

function qwgto ( x, omega, p2, p3, p4, integr )

!*****************************************************************************80
!
!! QWGTO defines the weight functions used by QC25O.
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which the weight function is evaluated.
!
!    Input, real ( kind = 8 ) OMEGA, the factor multiplying X.
!
!    Input, real ( kind = 8 ) P2, P3, P4, parameters that are not used.
!
!    Input, integer ( kind = 8 ) INTEGR, specifies which weight function is used:
!    1. W(X) = cos ( OMEGA * X )
!    2, W(X) = sin ( OMEGA * X )
!
!    Output, real ( kind = 8 ) QWGTO, the value of the weight function at X.
!
  implicit none

  integer ( kind = 8 ) integr
  real ( kind = 8 ) omega
  real ( kind = 8 ) p2
  real ( kind = 8 ) p3
  real ( kind = 8 ) p4
  real ( kind = 8 ) qwgto
  real ( kind = 8 ) x

  if ( integr == 1 ) then
    qwgto = cos ( omega * x )
  else if ( integr == 2 ) then
    qwgto = sin ( omega * x )
  end if

  return
end function
function qwgts ( x, a, b, alfa, beta, integr )

!*****************************************************************************80
!
!! QWGTS defines the weight functions used by QC25S.
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which the weight function is evaluated.
!
!    Input, real ( kind = 8 ) A, B, the endpoints of the integration interval.
!
!    Input, real ( kind = 8 ) ALFA, BETA, exponents that occur in the weight function.
!
!    Input, integer ( kind = 8 ) INTEGR, specifies which weight function is used:
!    1. W(X) = (X-A)**ALFA * (B-X)**BETA
!    2, W(X) = (X-A)**ALFA * (B-X)**BETA * log (X-A)
!    3, W(X) = (X-A)**ALFA * (B-X)**BETA * log (B-X)
!    4, W(X) = (X-A)**ALFA * (B-X)**BETA * log (X-A) * log(B-X)
!
!    Output, real ( kind = 8 ) QWGTS, the value of the weight function at X.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) alfa
  real ( kind = 8 ) b
  real ( kind = 8 ) beta
  integer ( kind = 8 ) integr
  real ( kind = 8 ) qwgts
  real ( kind = 8 ) x

  if ( integr == 1 ) then
    qwgts = ( x - a )**alfa * ( b - x )**beta
  else if ( integr == 2 ) then
    qwgts = ( x - a )**alfa * ( b - x )**beta * log ( x - a )
  else if ( integr == 3 ) then
    qwgts = ( x - a )**alfa * ( b - x )**beta * log ( b - x )
  else if ( integr == 4 ) then
    qwgts = ( x - a )**alfa * ( b - x )**beta * log ( x - a ) * log ( b - x )
  end if

  return
end function
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 8 ) d
  integer ( kind = 8 ) h
  integer ( kind = 8 ) m
  integer ( kind = 8 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 8 ) n
  integer ( kind = 8 ) s
  integer ( kind = 8 ) values(8)
  integer ( kind = 8 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end subroutine
end module QUADPACK