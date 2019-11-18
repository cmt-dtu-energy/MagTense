

function [N,P] = getN_Triangle( v, r )
%% Calculates and returns the demag tensor, N, given the four vertices,
%% v (size 3,4) and at the positions r (size 3,n)) from a triangluar surface
%% defined by vertices 1-3. The last vertex defines the orientation of the normal to the triangular
%% surface such that the normal points away from this vertex.
%% the array v is the four vertices organized as column vectors, 
%% i.e. v(:,1) is the first vertex etc.


%test if the three vertices are colinear. If so, return without further
%calculation
numErr = 1e-15;
d1 = norm(v(:,1)-v(:,2));
d2 = norm(v(:,2)-v(:,3));
if abs( norm(v(:,1)-v(:,3)) - (d1+d2)) < numErr
    disp('Err. The vertices are collinear to within 1e-15');
    N = [];
    return
end

%ensure the vertices are ordered such that the largest angle is at the
%middle vertex
angles = zeros(3,1);
angles(1) = acos( dot(v(:,1)-v(:,2),v(:,1)-v(:,3)) / (norm(v(:,1)-v(:,2))*norm(v(:,1)-v(:,3))) );
angles(2) = acos( dot(v(:,2)-v(:,1),v(:,2)-v(:,3)) / (norm(v(:,2)-v(:,1))*norm(v(:,2)-v(:,3))) );
angles(3) = acos( dot(v(:,3)-v(:,2),v(:,3)-v(:,1)) / (norm(v(:,3)-v(:,2))*norm(v(:,3)-v(:,1))) );
%order the vertices after angle where the largest angle is in the middle
[angles,ind]=sort(angles);
v = v(:,[ind' 4]);
v = v(:,[1,3,2,4]);
angles = angles([1,3,2]);
%find the transformation matrix from the local to the global system
%first unit vector
e1 = v(:,1)-v(:,3);
e1 = e1/norm(e1);
%third unit vector
e3 = cross( e1, v(:,2)-v(:,3) );
e3 = e3 / norm(e3);


%ensure the fourth vertex is not in the triangluar plane

if abs( dot( (v(:,4)-v(:,1)), e3 ) ) <= numErr 
    disp('Vertices are all in the same plane');
    N = [];
    return
end

%check which way the surface normal should point
if dot( e3, v(:,4)-v(:,1) ) > 0
    %switch around v1 and v3
    v = v(:,[3,2,1,4]);
    angles = angles([3,2,1]);
    %find the two first unit vectors again
    e1 = v(:,1)-v(:,3);
    e1 = e1/norm(e1);
    %third unit vector
    e3 = cross( e1, v(:,2)-v(:,3) );
    e3 = e3 / norm(e3);
end

%second unit vector
e2 = cross( e3, e1 );

%p-matrix
P = [e1'; e2'; e3']';
%the inverse is the transpose as P is orthogonal
Pinv = P';

%find the point D in order to translate the triangle to coincide with the
%Origin
D = cos(angles(3)) * norm( v(:,2)-v(:,3) ) * e1 + v(:,3);

%transformed coordinates (from global to local)
vp = zeros(3,3);
for i=1:3;vp(:,i) = Pinv * ( v(:,i)-D );end


threshold = 1e-10;


%% find the tensor in the local system

%transform the point of interest to the local system
r_t = Pinv * ( r - D );


for i=1:3
    if abs(r_t(i)) < threshold
        sgn = sign(r_t(i));
        if sgn == 0
            sgn = 1;
        end
        r_t(i) = sgn * threshold;
    end
end


%right triangle in the first quadrant
N1_a = zeros(3,3);
N1_n = N1_a;
N1_a(1,3) = Nxz( r_t, vp(1,1), vp(2,2) );

[N1_a(2,3),N1_n(2,3)] = Nyz( r_t, vp(1,1), vp(2,2) );

[N1_a(3,3),N1_n(3,3)] = Nzz( r_t, vp(1,1), vp(2,2) );


%right triangle in the second quadrant
N2_a = zeros(3,3);
N2_n = N2_a;

N2_a(1,3) = -Nxz( r_t, vp(1,3), vp(2,2) );

N2_a(2,3) = -Nyz( r_t, vp(1,3), vp(2,2) );

N2_a(3,3) = -Nzz( r_t, vp(1,3), vp(2,2) );


%total N in the local coordinate system
Na = N1_a + N2_a;


%Apply change of basis
N = P * Na * Pinv;


end


function [N1] = Nxz( r, l, h )
%% Returns the Nxz tensor component in the local coordinate system
    N1 = -1/(4*pi) .* ( F_Nxz(r,h,l,h) - F_Nxz(r,0,l,h) - ( G_Nxz(r,h) - G_Nxz(r,0) ) );
    
    %N2 = N1;%-1/(4*pi) .* numint_xz( r, l, h );
end


function val = numint_xz( r, l, h )
        
    xmax = @(y) l.*(1-y./h);

    fct = @(y,x) dPdx(r,x,y);

    val = integral2( fct, 0, h, 0, xmax );
    
%     fct_out = @(y) int_loop( y, r, l, h );
%     
%     val2 = integral( fct_out, 0, h );
%     
%     fct_out = @(y) -1./sqrt((r(1)-l*(1-y/h)).^2+(r(2)-y).^2+r(3).^2 ) - (-1./sqrt( r(1).^2+(r(2)-y).^2+r(3).^2 ));
% 
%     val3 = integral( fct_out, 0, h );
end


function dP = dPdx(r,xp,yp)
    dP = -( r(1)-xp ) ./ ( (r(1)-xp).^2 + (r(2)-yp).^2 + r(3).^2 ).^(3./2.);
end

function F = F_Nxz( r, yp, l, h )
%% r is the position vectors, i.e. r(1,:) are the x-values, r(2,:) the y-values and r(3,:) the z-values
%RUBI solution (Int):
    arg = (l^2 - l*r(1,:) + h*r(2,:) - h*yp.*(1+l^2/h^2))./...
        ( sqrt(h^2+l^2).*sqrt( l^2 - 2*l*r(1,:) + r(1,:).^2 + r(2,:).^2 - 2*yp.*(l^2-l*r(1,:)+h*r(2,:))/h + yp^2*(1+l^2/h^2) + r(3,:).^2) );
    if abs(arg-1) < 1e-15 
        arg = 1-1e-15;
    end
    F = h / sqrt( h^2 + l^2 ) .* atanh( arg );

    %Mathematica solution (Integrate)
%    F2 = -h/sqrt(h^2+l^2) * log( h*l*r(1,:) + l^2*(yp-h) + h^2*(yp-r(2,:)) + ...
 %       h*sqrt(h^2+l^2).*sqrt( r(1,:).^2 + l^2*(h-yp)^2/h^2 + (r(2,:)-yp).^2 + 2*l*r(1,:)*(yp-h)/h + r(3,:).^2 ) ) ;
        
end




function G = G_Nxz( r, yp )
    G = atanh( ( r(2,:)-yp )./ (sqrt(r(1,:).^2+(r(2,:)-yp).^2+r(3,:).^2)) );
end

function [N1,N2] = Nyz( r, l, h )

   %% returns the Nyz component of the tensor in the local coordinate system
   N1 = -1/(4*pi) * ( K_Nyz(r,l,l,h) - K_Nyz(r,0,l,h) - ( L_Nyz(r,l) - L_Nyz(r,0) ) );
   
   N2 = N1;%-1/(4*pi) * numint_yz( r, l, h );
end

function val = numint_yz( r, l, h )
        
     ymax = @(x) h.*(1-x./l);

    fct = @(x,y) dPdy(r,x,y);

    val = integral2( fct, 0, l, 0, ymax );
end


function dP = dPdy(r,xp,yp)
    dP = -( r(2)-yp ) ./ ( (r(1)-xp).^2 + (r(2)-yp).^2 + r(3).^2 ).^(3./2.);
end

function K = K_Nyz( r, xp, l, h )

    K = l/sqrt(h^2+l^2) * atanh( (h^2-h*r(2,:)+l*r(1,:)-l*xp*(1+h^2/l^2) ) ./...
        (sqrt(h^2+l^2).*sqrt(h^2-2*h*r(2,:)+r(1,:).^2 + r(2,:).^2 - 2*xp*(h^2+l*r(1,:)-h*r(2,:))/l +xp^2*(1+h^2/l^2) + r(3,:).^2)) );

end


function L = L_Nyz( r, xp )
    L = atanh( (r(1,:) - xp) ./ (sqrt((r(1,:)-xp).^2+r(2,:).^2+r(3,:).^2)) );
end

function [N1,N2] = Nzz( r, l, h )

    N1 = -1/(4*pi) .* ( P_Nzz( r, l, l, h ) - P_Nzz( r, 0, l, h ) - ( Q_Nzz(r,l) - Q_Nzz(r,0) ) );
    N2 = N1;%-1/(4*pi) * numint_zz( r, l, h );
end
    
function val = numint_zz( r, l, h )
        
    xmax = @(y) l.*(1-y./h);

    fct = @(y,x) dPdz(r,x,y);

    val = integral2( fct, 0, h, 0, xmax );
end


function dP = dPdz(r,xp,yp)
    dP = -( r(3) ) ./ ( (r(1)-xp).^2 + (r(2)-yp).^2 + r(3).^2 ).^(3./2.);
end

function P = P_Nzz( r, xp, l, h )

    P = atan( ( r(1,:).*(h-r(2,:)) - xp*(h*(1-r(1,:)./l)-r(2,:)) - h*(r(1,:).^2+r(3,:).^2)/l ) ./...
        (r(3,:).*sqrt(h^2+r(1,:).^2 + xp^2*(1+h^2/l^2) - 2*h*r(2,:) + r(2,:).^2 - 2*xp*(h^2+l*r(1,:)-h*r(2,:))/l + r(3,:).^2 )) );
    %P = atan2( ( r(1,:).*(h-r(2,:)) - xp*(h*(1-r(1,:)./l)-r(2,:)) - h*(r(1,:).^2+r(3,:).^2)/l ),...
    %    (r(3,:).*sqrt(h^2+r(1,:).^2 + xp^2*(1+h^2/l^2) - 2*h*r(2,:) + r(2,:).^2 - 2*xp*(h^2+l*r(1,:)-h*r(2,:))/l + r(3,:).^2 )) );
end

function Q = Q_Nzz( r, xp )
    Q = -atan( (r(1,:)-xp).*r(2,:) ./ ( r(3,:) .* sqrt( (r(1,:)-xp).^2 + r(2,:).^2 + r(3,:).^2 ) ) );
    %Q = -atan2( (r(1,:)-xp).*r(2,:), ( r(3,:) .* sqrt( (r(1,:)-xp).^2 + r(2,:).^2 + r(3,:).^2 ) ) );
end

function Q = Q_nzz_arg( r, xp )
    Q = (r(1,:)-xp).*r(2,:) ./ ( r(3,:) .* sqrt( (r(1,:)-xp).^2 + r(2,:).^2 + r(3,:).^2 ) ) ;
end




