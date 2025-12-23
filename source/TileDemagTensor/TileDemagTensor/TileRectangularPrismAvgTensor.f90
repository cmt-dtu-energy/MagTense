module TileRectangularPrismAvgTensor

    implicit none
    
    contains

        !Function to calculate the volume of a box, when given the corner coordinates (X1,X2,Y1,Y2,Z1,Z2)
        function volume(X1,X2,Y1,Y2,Z1,Z2) result(vol)
            real(8), intent(in) :: X1,X2,Y1,Y2,Z1,Z2
            real(8) :: vol
            vol = abs((X1-X2)*(Y1-Y2)*(Z1-Z2))
        end function volume

        !Function to calculate euclidian distance from origo to a point (X,Y,Z)
        function dist(X,Y,Z) result(D)
            real(8), intent(in) :: X,Y,Z
            real(8) :: D
            D = sqrt(X**2 + Y**2 + Z**2)
            !if(D <= 0.0d0) then
            !    print *, 'dist issue!'
            !end if
        end function dist

        !Function F1 corresponding to eq. 15 in The Fukushima paper:
        ! Volume Average Demagnetizing Tensor of Rectangular Prisms, 
        ! Hiroshi Fukushima, Yoshinobu Nakatani, and Nobuo Hayashi, Member, ZEEE 
        function F1(X,Y,Z) result(res)
            real(8), intent(in) :: X,Y,Z
            real(8) :: res
            real(8) :: Xv,Yv,Zv,D

            Xv = X
            Yv = Y
            Zv = Z
            D = dist(Xv,Yv,Zv)

            if(Xv == 0.0d0) Xv = 1.0d-5
            if(Yv == 0.0d0) Yv = 1.0d-5
            if(Zv == 0.0d0) Zv = 1.0d-5
            if(D-Yv == 0.0d0) then
                !print *, 'D-Y=0'
                Yv = Yv + 1.0d-5
            end if
            if(D-Zv == 0.0d0) then
                !print *, 'D-Z=0 D=', D, ' Z=', Zv, ' X=', Xv, ' Y=', Yv
                Zv = Zv + 1.0d-5
            end if

            res = Xv*Yv*Zv*atan(Yv*Zv/(Xv*D)) &
                  +0.5d0*Yv*(Zv**2 - Xv**2)*log(abs(D-Yv)) &
                  +0.5d0*Zv*(Yv**2 - Xv**2)*log(abs(D-Zv)) &
                  + (1.0d0/6.0d0)*(Yv**2 + Zv**2 - 2.0d0*Xv**2)*D
        end function F1

        !Function F2 corresponding to eq. 16 in The Fukushima paper:
        ! Volume Average Demagnetizing Tensor of Rectangular Prisms, 
        ! Hiroshi Fukushima, Yoshinobu Nakatani, and Nobuo Hayashi, Member, ZEEE 
        function F2(X,Y,Z) result(res)
            real(8), intent(in) :: X,Y,Z
            real(8) :: res
            real(8) :: Xv,Yv,Zv,D,A,B,C

            Xv = X
            Yv = Y
            Zv = Z
            
            if(Xv == 0.0d0) Xv = 1.0d-12
            if(Yv == 0.0d0) Yv = 1.0d-12
            if(Zv == 0.0d0) Zv = 1.0d-12

            !D = dist(Xv,Yv,Zv)
            D = sqrt(Xv**2+Yv**2+Zv**2)

            if(D == 0.0d0) D = 1.0d-12

            A = D + Xv
            B = D + Yv
            C = D + Zv

            if(abs(A) == 0.0d0) A = 1.0d-12
            if(abs(B) == 0.0d0) B = 1.0d-12
            if(abs(C) == 0.0d0) C = 1.0d-12

            res = -Xv*Yv*Zv*log(abs(C)) &
                  + (1.0d0/6.0d0)*Yv*(Yv**2 - 3.0d0*Zv**2)*log(abs(A)) &
                  + (1.0d0/6.0d0)*Xv*(Xv**2 - 3.0d0*Zv**2)*log(abs(B)) &
                  + 0.5d0*Xv**2*Zv*atan(Yv*Zv/(Xv*D)) &
                  + 0.5d0*Yv**2*Zv*atan(Xv*Zv/(Yv*D)) &
                  + (1.0d0/6.0d0)*Zv**3*atan(Xv*Yv/(Zv*D)) &
                  + (1.0d0/3.0d0)*Xv*Yv*D

            if (abs(res)>100) print *, res
        end function F2
        
        !Function definite_integral corresponding to eq. 13 in The Fukushima paper:
        ! Volume Average Demagnetizing Tensor of Rectangular Prisms, 
        ! Hiroshi Fukushima, Yoshinobu Nakatani, and Nobuo Hayashi, Member, ZEEE 
        function definite_integral(func,X1,X2,Y1,Y2,Z1,Z2) result(res)
            real(8), intent(in) :: X1,X2,Y1,Y2,Z1,Z2
            interface
                function func(X,Y,Z) result(f)
                    real(8), intent(in) :: X,Y,Z
                    real(8) :: f
                end function func
            end interface
            real(8) :: res

            res = func(X2,Y2,Z2) - func(X1,Y2,Z2) &
                - func(X2,Y1,Z2) + func(X1,Y1,Z2) &
                - func(X2,Y2,Z1) + func(X1,Y2,Z1) &
                + func(X2,Y1,Z1) - func(X1,Y1,Z1)
        end function definite_integral
    
    
end module TileRectangularPrismAvgTensor
    