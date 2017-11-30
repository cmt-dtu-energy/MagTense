module SPECIALFUNCTIONS
    implicit none
    contains 
    subroutine elit ( hk, phi, fe, ee )

!*****************************************************************************80
!
!! ELIT: complete and incomplete elliptic integrals F(k,phi) and E(k,phi).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    12 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) HK, the modulus, between 0 and 1.
!
!    Input, real ( kind = 8 ) PHI, the argument in degrees.
!
!    Output, real ( kind = 8 ) FE, EE, the values of F(k,phi) and E(k,phi).
!


  real ( kind = 8 ) a
  real ( kind = 8 ) a0
  real ( kind = 8 ) b
  real ( kind = 8 ) b0
  real ( kind = 8 ) c
  real ( kind = 8 ) ce
  real ( kind = 8 ) ck
  real ( kind = 8 ) d
  real ( kind = 8 ) d0
  real ( kind = 8 ) ee
  real ( kind = 8 ) fac
  real ( kind = 8 ) fe
  real ( kind = 8 ) g
  real ( kind = 8 ) hk
  integer ( kind = 4 ) n
  real ( kind = 8 ) phi
  real ( kind = 8 ) pi
  real ( kind = 8 ) r

  g = 0.0D+00
  pi = 3.14159265358979D+00
  a0 = 1.0D+00
  b0 = sqrt ( 1.0D+00 - hk * hk )
  d0 = ( pi / 180.0D+00 ) * phi
  r = hk * hk

  if ( hk == 1.0D+00 .and. phi == 90.0D+00 ) then

    fe = 1.0D+300
    ee = 1.0D+00

  else if ( hk == 1.0D+00 ) then

    fe = log ( ( 1.0D+00 + sin ( d0 ) ) / cos ( d0 ) )
    ee = sin ( d0 )

  else

    fac = 1.0D+00
    do n = 1, 40
      a = ( a0 + b0 ) /2.0D+00
      b = sqrt ( a0 * b0 )
      c = ( a0 - b0 ) / 2.0D+00
      fac = 2.0D+00 * fac
      r = r + fac * c * c
      if ( phi /= 90.0D+00 ) then
        d = d0 + atan ( ( b0 / a0 ) * tan ( d0 ) )
        g = g + c * sin( d )
        d0 = d + pi * int ( d / pi + 0.5D+00 )
      end if
      a0 = a
      b0 = b
      if ( c < 1.0D-07 ) then
        exit
      end if
    end do

    ck = pi / ( 2.0D+00 * a )
    ce = pi * ( 2.0D+00 - r ) / ( 4.0D+00 * a )
    if ( phi == 90.0D+00 ) then
      fe = ck
      ee = ce
    else
      fe = d / ( fac * a )
      ee = fe * ce / ck + g
    end if

  end if

  return
end
subroutine elit3 ( phi, hk, c, el3 )

!*****************************************************************************80
!
!! ELIT3 computes the elliptic integral of the third kind.
!
!  Discussion:
!
!    Gauss-Legendre quadrature is used.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    14 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) PHI, the argument in degrees.
!
!    Input, real ( kind = 8 ) HK, the modulus, between 0 and 1.
!
!    Input, real ( kind = 8 ) C, the parameter, between 0 and 1.
!
!    Output, real ( kind = 8 ) EL3, the value of the elliptic integral
!    of the third kind.
!
  

  real ( kind = 8 ) c
  real ( kind = 8 ) c0
  real ( kind = 8 ) c1
  real ( kind = 8 ) c2
  real ( kind = 8 ) el3
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) hk
  integer ( kind = 4 ) i 
  logical lb1
  logical lb2
  real ( kind = 8 ) phi
  real ( kind = 8 ), dimension ( 10 ), save :: t = (/ &
    0.9931285991850949D+00, 0.9639719272779138D+00, &
    0.9122344282513259D+00, 0.8391169718222188D+00, &
    0.7463319064601508D+00, 0.6360536807265150D+00, &
    0.5108670019508271D+00, 0.3737060887154195D+00, &
    0.2277858511416451D+00, 0.7652652113349734D-01 /)
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ), dimension ( 10 ), save :: w = (/ &
    0.1761400713915212D-01, 0.4060142980038694D-01, &
    0.6267204833410907D-01, 0.8327674157670475D-01, &
    0.1019301198172404D+00, 0.1181945319615184D+00, &
    0.1316886384491766D+00, 0.1420961093183820D+00, &
    0.1491729864726037D+00, 0.1527533871307258D+00 /)

  lb1 = ( hk == 1.0D+00 ) .and. ( abs ( phi - 90.0D+00 ) <= 1.0D-08 )

  lb2 = c == 1.0D+00 .and. abs ( phi - 90.0D+00 ) <= 1.0D-08

  if ( lb1 .or. lb2 ) then
    el3 = 1.0D+300
    return
  end if

  c1 = 0.87266462599716D-02 * phi
  c2 = c1

  el3 = 0.0D+00
  do i = 1, 10
    c0 = c2 * t(i)
    t1 = c1 + c0
    t2 = c1 - c0
    f1 = 1.0D+00 / ( ( 1.0D+00 - c * sin(t1) * sin(t1) ) &
      * sqrt ( 1.0D+00 - hk * hk * sin ( t1 ) * sin ( t1 ) ) )
    f2 = 1.0D+00 / ( ( 1.0D+00 - c * sin ( t2 ) * sin ( t2 ) ) &
      * sqrt( 1.0D+00 - hk * hk * sin ( t2 ) * sin ( t2 ) ) )
    el3 = el3 + w(i) * ( f1 + f2 )
  end do

  el3 = c1 * el3

  return
end

END MODULE SPECIALFUNCTIONS