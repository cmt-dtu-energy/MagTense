module TileRectangularPrismAvgTensor

    implicit none
    
    contains

        ! Function to calculate the volume of a box, when given the corner coordinates (X1,X2,Y1,Y2,Z1,Z2)
        function volume(X1,X2,Y1,Y2,Z1,Z2) result(vol)
            real, intent(in) :: X1,X2,Y1,Y2,Z1,Z2
            real :: vol
            vol = abs((X1-X2)*(Y1-Y2)*(Z1-Z2))
        end function volume

        !Function to calculate euclidian distance from origo to a point (X,Y,Z)
        function dist(X,Y,Z) result(D)
            real, intent(in) :: X,Y,Z
            real :: D
            D = sqrt(X**2 + Y**2 + Z**2)
        end function dist

        ! Function F1 corresponding to eq. 15 in The Fukushima paper:
        ! Volume Average Demagnetizing Tensor of Rectangular Prisms, 
        ! Hiroshi Fukushima, Yoshinobu Nakatani, and Nobuo Hayashi, Member, ZEEE 
        function F1(X,Y,Z) result(res)
            real, intent(in) :: X,Y,Z
            real :: res
            real :: Xv,Yv,Zv,D
            real :: eps

            eps = 10d-8

            Xv = X
            Yv = Y
            Zv = Z
            D = dist(Xv,Yv,Zv)

            if(Xv == 0.0d0) Xv = eps
            if(Yv == 0.0d0) Yv = eps
            if(Zv == 0.0d0) Zv = eps
            if(D-Yv == 0.0d0) then
                !print *, 'D-Y=0'
                Yv = Yv + eps
            end if
            if(D-Zv == 0.0d0) then
                !print *, 'D-Z=0 D=', D, ' Z=', Zv, ' X=', Xv, ' Y=', Yv
                Zv = Zv + eps
            end if

            res = Xv*Yv*Zv*atan(Yv*Zv/(Xv*D)) &
                  +0.5d0*Yv*(Zv**2 - Xv**2)*log(abs(D-Yv)) &
                  +0.5d0*Zv*(Yv**2 - Xv**2)*log(abs(D-Zv)) &
                  + (1.0d0/6.0d0)*(Yv**2 + Zv**2 - 2.0d0*Xv**2)*D
        end function F1

        ! Function F2 corresponding to eq. 16 in Fukushima et al. 2002 [DOI: 10.1109/20.650225],
        ! Volume Average Demagnetizing Tensor of Rectangular Prisms
        function F2(X,Y,Z) result(res)
            real, intent(in) :: X,Y,Z
            real :: res
            real :: Xv,Yv,Zv,D,A,B,C
            real :: eps

            eps = 1.0d-8

            Xv = X
            Yv = Y
            Zv = Z
            
            if(Xv == 0.0d0) Xv = eps
            if(Yv == 0.0d0) Yv = eps
            if(Zv == 0.0d0) Zv = eps

            D = sqrt(Xv**2+Yv**2+Zv**2)

            if(D == 0.0d0) D = eps

            A = D + Xv
            B = D + Yv
            C = D + Zv

            if(abs(A) == 0.0d0) A = eps
            if(abs(B) == 0.0d0) B = eps
            if(abs(C) == 0.0d0) C = eps

            res = -Xv*Yv*Zv*log(abs(C)) &
                  + (1.0d0/6.0d0)*Yv*(Yv**2 - 3.0d0*Zv**2)*log(abs(A)) &
                  + (1.0d0/6.0d0)*Xv*(Xv**2 - 3.0d0*Zv**2)*log(abs(B)) &
                  + 0.5d0*Xv**2*Zv*atan(Yv*Zv/(Xv*D)) &
                  + 0.5d0*Yv**2*Zv*atan(Xv*Zv/(Yv*D)) &
                  + (1.0d0/6.0d0)*Zv**3*atan(Xv*Yv/(Zv*D)) &
                  + (1.0d0/3.0d0)*Xv*Yv*D

        end function F2
        
        ! Function AvgN_prism_definite_integral corresponding to eq. 13 in Fukushima et al. 2002 [DOI: 10.1109/20.650225],
        ! Volume Average Demagnetizing Tensor of Rectangular Prisms
        function AvgN_prism_definite_integral(func,X1,X2,Y1,Y2,Z1,Z2) result(res)
            real, intent(in) :: X1,X2,Y1,Y2,Z1,Z2
            interface
                function func(X,Y,Z) result(f)
                    real(8), intent(in) :: X,Y,Z
                    real(8) :: f
                end function func
            end interface
            real :: res

            res = func(X2,Y2,Z2) - func(X1,Y2,Z2) &
                - func(X2,Y1,Z2) + func(X1,Y1,Z2) &
                - func(X2,Y2,Z1) + func(X1,Y2,Z1) &
                + func(X2,Y1,Z1) - func(X1,Y1,Z1)
        end function AvgN_prism_definite_integral
    
    
end module TileRectangularPrismAvgTensor
    