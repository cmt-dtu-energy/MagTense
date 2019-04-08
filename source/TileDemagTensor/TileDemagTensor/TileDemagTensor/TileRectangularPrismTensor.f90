module TileRectanagularPrismTensor

    implicit none
    
    contains
    
    function getF_limit(a,b,c,x,y,z,func)
    real,intent(in) :: x,y,z,a,b,c
    real :: getF_limit
    real :: xh,yh,zh,xl,yl,zl,nom_l,nom_h,denom_l,denom_h,nom,denom
    real,parameter :: lim_scl_h=1.0001,lim_scl_l=0.9999
    real :: func
    external func
   
        !Find the limit
        xh = x * lim_scl_h
        yh = y * lim_scl_h
        zh = z * lim_scl_h
    
        xl = x * lim_scl_l
        yl = y * lim_scl_l
        zl = z * lim_scl_l
    
        
        nom_l = func(a,b,c,xl,yl,zl)  * func(-a,-b,c,xl,yl,zl) * func(a,-b,-c,xl,yl,zl) * func(-a,b,-c,xl,yl,zl)
        nom_h = func(a,b,c,xh,yh,zh)  * func(-a,-b,c,xh,yh,zh) * func(a,-b,-c,xh,yh,zh) * func(-a,b,-c,xh,yh,zh)
    
        denom_l = func(a,-b,c,xl,yl,zl) * func(-a,b,c,xl,yl,zl)  * func(a,b,-c,xl,yl,zl)  * func(-a,-b,-c,xl,yl,zl)
        denom_h = func(a,-b,c,xh,yh,zh) * func(-a,b,c,xh,yh,zh)  * func(a,b,-c,xh,yh,zh)  * func(-a,-b,-c,xh,yh,zh)
    
        nom = 0.5 * ( nom_l + nom_h )
    
        denom = 0.5 * ( denom_l + denom_h )
    
    
        getF_limit = nom / denom

    end function getF_limit

    !:: Dims is defined as the dimensions (a/2.,b/2.,c/2.) for each prism
    !:: As the code below assumes that dims is center +/- dims, we have to divide dims by two


    !::function f from eq. 11
    function f_3D( a, b, c, x, y, z )
    real,intent(in) :: a,b,c,x,y,z
    real :: f_3D, f_3D_l, f_3D_h, xh,xl
    real,parameter :: lim_scl_h=1.0001,lim_scl_l=0.9999

    !::Numerical problem: atan(inf) = pi/2 and when a/2 = x then f_3D returns NaN (inf). 
    !::So we need to check for this and return a large number to ensure numerical stability

    if ( a/2. - x .eq. 0 ) then
        !f_3D = sign(1e10, (b/2. - y) * (c/2. - z))
        xh = x * lim_scl_h
        xl = x * lim_scl_l
    
        f_3D_h = (b/2. - y) * (c/2. - z) / ( (a/2. - xh) * sqrt( (a/2. - xh)**2 + (b/2. - y)**2 + (c/2. - z)**2 ) )
        f_3D_l = (b/2. - y) * (c/2. - z) / ( (a/2. - xl) * sqrt( (a/2. - xl)**2 + (b/2. - y)**2 + (c/2. - z)**2 ) )
    
        f_3D = 0.5 * ( f_3D_h + f_3D_l )
    
    else    
        f_3D = (b/2. - y) * (c/2. - z) / ( (a/2. - x) * sqrt( (a/2. - x)**2 + (b/2. - y)**2 + (c/2. - z)**2 ) )
    endif



    return
    end function f_3D


    !::function g from eq. 11
    function g_3D( a, b, c, x, y, z )
    real,intent(in) :: a,b,c,x,y,z
    real :: g_3D, g_3D_l, g_3D_h, yh,yl
    real,parameter :: lim_scl_h=1.0001,lim_scl_l=0.9999

    !::Numerical problem: atan(inf) = pi/2 and when a/2 = x then f_3D returns NaN (inf). 
    !::So we need to check for this and return a large number to ensure numerical stability
    if ( b/2. - y .eq. 0 ) then
        yh = y * lim_scl_h
        yl = y * lim_scl_l
    
        g_3D_h = (a/2. - x) * (c/2. - z) / ( (b/2. - yh) * sqrt( (a/2. - x)**2 + (b/2. - yh)**2 + (c/2. - z)**2 ) )
        g_3D_l = (a/2. - x) * (c/2. - z) / ( (b/2. - yl) * sqrt( (a/2. - x)**2 + (b/2. - yl)**2 + (c/2. - z)**2 ) )
    
        g_3D = 0.5 * ( g_3D_h + g_3D_l )
    else    
        g_3D = (a/2. - x) * (c/2. - z) / ( (b/2. - y) * sqrt( (a/2. - x)**2 + (b/2. - y)**2 + (c/2. - z)**2 ) )
    endif
    return
    end function g_3D


    !::function h from eq. 11
    function h_3D( a, b, c, x, y, z )
    real,intent(in) :: a,b,c,x,y,z
    real :: h_3D,h_3D_l,h_3D_h,zh,zl
    real,parameter :: lim_scl_h=1.0001,lim_scl_l=0.9999
    !::Numerical problem: atan(inf) = pi/2 and when a/2 = x then f_3D returns NaN (inf). 
    !::So we need to check for this and return a large number to ensure numerical stability
    if ( c/2. - z .eq. 0 ) then
        zh = z * lim_scl_h
        zl = z * lim_scl_l
    
        h_3D_h = (a/2. - x) * (b/2. - y) / ( (c/2. - zh) * sqrt( (a/2. - x)**2 + (b/2. - y)**2 + (c/2. - zh)**2 ) )
        h_3D_l = (a/2. - x) * (b/2. - y) / ( (c/2. - zl) * sqrt( (a/2. - x)**2 + (b/2. - y)**2 + (c/2. - zl)**2 ) )
    
        h_3D = 0.5 * ( h_3D_l + h_3D_h )
    
    else    
        h_3D = (a/2. - x) * (b/2. - y) / ( (c/2. - z) * sqrt( (a/2. - x)**2 + (b/2. - y)**2 + (c/2. - z)**2 ) )
    endif
    return
    end function h_3D


    !::function F from eq. 11
    function FF_3D( a, b, c, x, y, z )
    real,intent(in) :: a,b,c,x,y,z
    real :: FF_3D

    FF_3D = (c/2. - z) + sqrt( (a/2. - x)**2 + (b/2. - y)**2 + (c/2. - z)**2 )

    return
    end function FF_3D


    !::function G from eq. 11
    function GG_3D( a, b, c, x, y, z )
    real,intent(in) :: a,b,c,x,y,z
    real :: GG_3D

    GG_3D = (a/2. - x) + sqrt( (a/2. - x)**2 + (b/2. - y)**2 + (c/2. - z)**2 )

    return
    end function GG_3D


    !::function H from eq. 11
    function HH_3D( a, b, c, x, y, z )
    real,intent(in) :: a,b,c,x,y,z
    real :: HH_3D

    HH_3D = (b/2. - y) + sqrt( (a/2. - x)**2 + (b/2. - y)**2 + (c/2. - z)**2 )

    return
    end function HH_3D
    
end module TileRectanagularPrismTensor
    