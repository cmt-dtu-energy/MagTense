module SPECIALFUNCTIONS
    implicit none
    contains 
    
    !---------------------------------------------------------------------------    
    !> @author Kaspar K. Nielsen, kasparkn@gmail.com, DTU, 2019
    !> @brief
    !> Simple sorting function. Code implemented from Stack-overflow: https://stackoverflow.com/questions/54005339/sorting-an-array-from-lowest-to-greatest-fortran
    !> @param[in] arr the array to be sorted (n,1)
    !> @param[inout] arr_out sorted array (ascending order)
    !> @param[inout] ind indices such that arr_out = arr(ind)
    !---------------------------------------------------------------------------            
    subroutine simple_sort( arr, arr_out, ind )
        
    real,dimension(:),intent(in) :: arr
    real,dimension(:),intent(inout) :: arr_out
    integer,dimension(:),intent(inout) :: ind
    LOGICAL, DIMENSION(size(arr)) :: mk
    integer :: i,sz
    integer,dimension(1) :: tmp
    
    mk(:) = .true.
    
    
    sz = size(arr)
    do i = 1, sz
        arr_out(i) = MINVAL(arr,mk)
        tmp = MINLOC(arr,mk)
        ind(i) = tmp(1)
        mk(MINLOC(arr,mk)) = .FALSE.
    enddo
    
    
    end subroutine simple_sort
    
    
    !---------------------------------------------------------------------------    
    !> @author Rasmus Bjørk, rabj@dtu.dk, DTU, 2020
    !> @brief
    !> Simple function to find the unique elements in an array. Code implemented from Rosettacode: https://rosettacode.org/wiki/Remove_duplicate_elements#Fortran
    !> @param[in] arr the array to be sorted (n,1)
    !> @param[inout] arr_out sorted array (ascending order)
    !> @param[inout] ind indices such that arr_out = arr(ind)
    !---------------------------------------------------------------------------            
    subroutine simple_unique( arr, arr_out, k)
    
    real,dimension(:),intent(in) :: arr
    real,dimension(:),intent(inout) :: arr_out
    integer,intent(inout) :: k
    
    integer :: i
    
    arr_out(:) = 0;
    
    k = 1
    arr_out(1) = arr(1)
    do i=2,size(arr)
        ! if the number already exist in res check next
        if (any( arr_out == arr(i) )) cycle
        ! No match found so add it to the output
        k = k + 1
        arr_out(k) = arr(i)
    end do
    
    end subroutine simple_unique
    
    !---------------------------------------------------------------------------    
    !> @author Kaspar K. Nielsen, kasparkn@gmail.com, DTU, 2019
    !> @brief
    !> Calculates the cross product of a and b and returns it
    !---------------------------------------------------------------------------            
    FUNCTION cross(a, b)
      real, DIMENSION(3) :: cross
      real, DIMENSION(3), INTENT(IN) :: a, b

      cross(1) = a(2) * b(3) - a(3) * b(2)
      cross(2) = a(3) * b(1) - a(1) * b(3)
      cross(3) = a(1) * b(2) - a(2) * b(1)
    END FUNCTION cross

    
    subroutine elit ( hk, phi_, fe, ee )

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
!    Input, real ( kind = 8 ) PHI, the argument in degrees. (no, Kaspar changed it to radians on 13 November 2017)
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
  real ( kind = 8 ) phi,phi_
  real ( kind = 8 ) pi
  real ( kind = 8 ) r
  
  g = 0.0D+00
  pi = 3.14159265358979D+00
  phi = phi_ * 180./pi
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
end subroutine
subroutine elit3 ( phi_, hk, c, el3 )

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
!    Input, real ( kind = 8 ) PHI, the argument in degrees. (no, Kaspar changed it to radians on 13 November 2017)
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
  real ( kind = 8 ) phi,phi_
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

  !::convert from radians to degrees
  phi = phi_ * 180. / 3.14159265358979D+00
  
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
end subroutine

      FUNCTION ellf(phi,ak)
      REAL ellf,ak,phi
!    USES rf
      REAL s
      !s=sin(phi)
      !ellf=s*rf(cos(phi)**2,(1.-s*ak)*(1.+s*ak),1.)
      s=sin(phi)
      ellf=s*rf(cos(phi)**2,(1.-s**2*ak),1.)
      return
      END function


      FUNCTION elle(phi,ak)
      REAL elle,ak,phi
!     USES rd,rf
      REAL cc,q,s
      s=sin(phi)
      cc=cos(phi)**2
      !q=(1.-s*ak)*(1.+s*ak)
      q=(1.-s**2*ak)
      elle=s*(rf(cc,q,1.)-((s**2*ak))*rd(cc,q,1.)/3.)
      return
      END function


      FUNCTION ellpi(phi,en,ak)
      REAL ellpi,ak,en,phi
!     USES rf,rj
      REAL cc,enss,q,s
      s=sin(phi)
      enss=en*s*s
      cc=cos(phi)**2
      !q=(1.-s*ak)*(1.+s*ak)
      q=(1.-s**2*ak)
      ellpi=s*(rf(cc,q,1.)-enss*rj(cc,q,1.,1.+enss)/3.)
      return
      END function
      
      FUNCTION rf(x,y,z)
      REAL rf,x,y,z,ERRTOL,TINY,BIG,THIRD,C1,C2,C3,C4
      PARAMETER (ERRTOL=.08,TINY=1.e-308,BIG=1.E308,THIRD=1./3.,C1=1./24.,C2=.1,C3=3./44.,C4=1./14.)
      REAL alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt
      if(min(x,y,z).lt.0..or.min(x+y,x+z,y+z).lt.TINY.or.max(x,y,z).gt.BIG)pause 'invalid arguments in rf'
      xt=x
      yt=y
      zt=z
1     continue
        sqrtx=sqrt(xt)
        sqrty=sqrt(yt)
        sqrtz=sqrt(zt)
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
        xt=.25*(xt+alamb)
        yt=.25*(yt+alamb)
        zt=.25*(zt+alamb)
        ave=THIRD*(xt+yt+zt)
        delx=(ave-xt)/ave
        dely=(ave-yt)/ave
        delz=(ave-zt)/ave
      if(max(abs(delx),abs(dely),abs(delz)).gt.ERRTOL)goto 1
      e2=delx*dely-delz**2
      e3=delx*dely*delz
      rf=(1.+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave)
      return
      END FUNCTION

      FUNCTION rj(x,y,z,p)
      REAL rj,p,x,y,z,ERRTOL,TINY,BIG,C1,C2,C3,C4,C5,C6,C7,C8
      PARAMETER (ERRTOL=.05,TINY=2.5e-13,BIG=9.E11,C1=3./14.,C2=1./3.,C3=3./22.,C4=3./26.,C5=.75*C3,C6=1.5*C4,C7=.5*C2,C8=C3+C3)
!     USES rc,rf
      REAL a,alamb,alpha,ave,b,beta,delp,delx,dely,delz,ea,eb,ec,ed,ee,fac,pt,rcx,rho,sqrtx,sqrty,sqrtz,sum,tau,xt,yt,zt
      if(min(x,y,z).lt.0..or.min(x+y,x+z,y+z,abs(p)).lt.TINY.or.max(x,y,z,abs(p)).gt.BIG)pause 'invalid arguments in rj'
      sum=0.
      fac=1.
      if(p.gt.0.)then
        xt=x
        yt=y
        zt=z
        pt=p
      else
        xt=min(x,y,z)
        zt=max(x,y,z)
        yt=x+y+z-xt-zt
        a=1./(yt-p)
        b=a*(zt-yt)*(yt-xt)
        pt=yt+b
        rho=xt*zt/yt
        tau=p*pt/yt
        rcx=rc(rho,tau)
      endif
1     continue
        sqrtx=sqrt(xt)
        sqrty=sqrt(yt)
        sqrtz=sqrt(zt)
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
        alpha=(pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz)**2
        beta=pt*(pt+alamb)**2
        sum=sum+fac*rc(alpha,beta)
        fac=.25*fac
        xt=.25*(xt+alamb)
        yt=.25*(yt+alamb)
        zt=.25*(zt+alamb)
        pt=.25*(pt+alamb)
        ave=.2*(xt+yt+zt+pt+pt)
        delx=(ave-xt)/ave
        dely=(ave-yt)/ave
        delz=(ave-zt)/ave
        delp=(ave-pt)/ave
      if(max(abs(delx),abs(dely),abs(delz),abs(delp)).gt.ERRTOL)goto 1
      ea=delx*(dely+delz)+dely*delz
      eb=delx*dely*delz
      ec=delp**2
      ed=ea-3.*ec
      ee=eb+2.*delp*(ea-ec)
      rj=3.*sum+fac*(1.+ed*(-C1+C5*ed-C6*ee)+eb*(C7+delp*(-C8+delp*C4))+delp*ea*(C2-delp*C3)-C2*delp*ec)/(ave*sqrt(ave))
      if (p.le.0.) rj=a*(b*rj+3.*(rcx-rf(xt,yt,zt)))
      return
      END function

      FUNCTION rc(x,y)
      REAL rc,x,y,ERRTOL,TINY,SQRTNY,BIG,TNBG,COMP1,COMP2,THIRD,C1,C2,C3,C4
      PARAMETER (ERRTOL=.04,TINY=1.69e-38,SQRTNY=1.3e-19,BIG=3.E37,TNBG=TINY*BIG,COMP1=2.236/SQRTNY,COMP2=TNBG*TNBG/25.,THIRD=1./3.,C1=.3,C2=1./7.,C3=.375,C4=9./22.)
      REAL alamb,ave,s,w,xt,yt
      if(x.lt.0..or.y.eq.0..or.(x+abs(y)).lt.TINY.or.(x+abs(y)).gt.BIG.or.(y.lt.-COMP1.and.x.gt.0..and.x.lt.COMP2))pause 'invalid arguments in rc'
      if(y.gt.0.)then
        xt=x
        yt=y
        w=1.
      else
        xt=x-y
        yt=-y
        w=sqrt(x)/sqrt(xt)
      endif
1     continue
        alamb=2.*sqrt(xt)*sqrt(yt)+yt
        xt=.25*(xt+alamb)
        yt=.25*(yt+alamb)
        ave=THIRD*(xt+yt+yt)
        s=(yt-ave)/ave
      if(abs(s).gt.ERRTOL)goto 1
      rc=w*(1.+s*s*(C1+s*(C2+s*(C3+s*C4))))/sqrt(ave)
      return
      END function
      
      FUNCTION rd(x,y,z)
      REAL rd,x,y,z,ERRTOL,TINY,BIG,C1,C2,C3,C4,C5,C6
      PARAMETER (ERRTOL=.05,TINY=1.e-25,BIG=4.5E21,C1=3./14.,C2=1./6.,C3=9./22.,C4=3./26.,C5=.25*C3,C6=1.5*C4)
      REAL alamb,ave,delx,dely,delz,ea,eb,ec,ed,ee,fac,sqrtx,sqrty,sqrtz,sum,xt,yt,zt
      if(min(x,y).lt.0..or.min(x+y,z).lt.TINY.or.max(x,y,z).gt.BIG)pause 'invalid arguments in rd'
      xt=x
      yt=y
      zt=z
      sum=0.
      fac=1.
1     continue
        sqrtx=sqrt(xt)
        sqrty=sqrt(yt)
        sqrtz=sqrt(zt)
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
        sum=sum+fac/(sqrtz*(zt+alamb))
        fac=.25*fac
        xt=.25*(xt+alamb)
        yt=.25*(yt+alamb)
        zt=.25*(zt+alamb)
        ave=.2*(xt+yt+3.*zt)
        delx=(ave-xt)/ave
        dely=(ave-yt)/ave
        delz=(ave-zt)/ave
      if(max(abs(delx),abs(dely),abs(delz)).gt.ERRTOL)goto 1
      ea=delx*dely
      eb=delz*delz
      ec=ea-eb
      ed=ea-6.*eb
      ee=ed+ec+ec
      rd=3.*sum+fac*(1.+ed*(-C1+C5*ed-C6*delz*ee)+delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*sqrt(ave))
      return
      END function


end module
