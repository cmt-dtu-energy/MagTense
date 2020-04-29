function [P] = plotCircPiece( cP )
    thg = hgtransform;
    col = cP.color;
    
    n = 10;
    
    
    R = cP.r0+cP.dr/2;
    
    h = 0;%cP.offset(1);
    k = 0;%cP.offset(2);
    
    th = cP.theta0;
    
    if cos(th) >= 0 && sin(th)>=0
        %first quadrant
        a = cP.theta0-cP.dtheta/2;
        b = cP.theta0+cP.dtheta/2;            
    elseif cos(th) <0 && sin(th) >=0
        %second quadrant
        a = cP.theta0+cP.dtheta/2;
        b = cP.theta0-cP.dtheta/2;
    elseif cos(th) <0 && sin(th) <0
        %third
        a = cP.theta0-cP.dtheta/2;
        b = cP.theta0+cP.dtheta/2;
    elseif cos(th) >=0 && sin(th) <0
        %fourth 
        a = cP.theta0+cP.dtheta/2;
        b = cP.theta0-cP.dtheta/2;
    end
    t = linspace(a,b,n);
    %define the corners
    x1 = R .* cos( a );
    y1 = R .* sin( a );

    x2 = R .* cos( b );
    y2 = R .* sin( b );

    x3 = x2;
    y3 = y1;
        
    %%plot the top and bottom
    x = R.*cos(t) + h;
    y = R.*sin(t) + k;
        
    x = [x, zeros(1,n)+x(end)];
    y = [y, linspace(y(end),y(1),n)];
    
    zup = ones(1,length(x)).*(cP.z0 + cP.dz/2);
    
    Pup = patch(x,y,zup,col,'parent',thg);
    
    zdn = ones(1,length(x)).*(cP.z0 - cP.dz/2);
    
    Pdn = patch(x,y,zdn,col,'parent',thg);
    
    %plot the outer arc piece    
    x = R.*cos(t) + h;
    y = R.*sin(t) + k;
    zup = ones(1,length(x)).*(cP.z0 + cP.dz/2);
    zdn = ones(1,length(x)).*(cP.z0 - cP.dz/2);
    z = [zup fliplr(zdn) zup(1)];
    x = [x fliplr(x) x(1)];
    y = [y fliplr(y) y(1)];
    Pout = patch(x,y,z,col,'parent',thg);
    
    %%plot the x-constant surface
    x = zeros(1,n) + x2 + h;
    y = linspace(y2,y3,n) + k;
    
    x = [x fliplr(x) x(1)];
    y = [y fliplr(y) y(1)];
    zup = ones(1,length(x)).*(cP.z0 + cP.dz/2);
    zdn = ones(1,length(x)).*(cP.z0 - cP.dz/2);
    Pleft = patch(x,y,z,col,'parent',thg);
    
    %%plot the y-constant piece
    x = linspace(x3,x1,n) + h;
    y = zeros(1,n) + y3 + k;
    
    x = [x fliplr(x) x(1)];
    y = [y fliplr(y) y(1)];
    zup = ones(1,length(x)).*(cP.z0 + cP.dz/2);
    zdn = ones(1,length(x)).*(cP.z0 - cP.dz/2);
    Pbot = patch(x,y,z,col,'parent',thg);
    
    P = struct();
    P.P1 = Pup;
    P.P2 = Pdn;
    P.P3 = Pout;
    P.P4 = Pleft;
    P.P5 = Pbot;
    
    Rx = makehgtform('xrotate',cP.rotAngles(1));
    Ry = makehgtform('yrotate',cP.rotAngles(2));
    Rz = makehgtform('zrotate',cP.rotAngles(3));
    
    %after rotation, translate to the offset
    T2 = makehgtform('translate', cP.offset);
    
    %Note the order of rotation is about z-axis first, then y and finally
    %x. This is when rotating to the global coordinate system. Also note
    %that we translate at the end
    thg.Matrix = T2 * Rx * Ry * Rz;
    
    
end