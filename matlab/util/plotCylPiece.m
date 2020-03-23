function [P] = plotCylPiece( cylP )
    thg = hgtransform;
    col = cylP.color;
    a = cylP.theta0-cylP.dtheta/2;
    b = cylP.theta0+cylP.dtheta/2;
    t = linspace(a,b,10);
    
    r2 = cylP.r0+cylP.dr/2;
    r1 = cylP.r0-cylP.dr/2;
    h = cylP.offset(1);
    k = cylP.offset(2);
    
    x2 = r2*cos(t) + h;
    y2 = r2*sin(t) + k;
    
    x1 = r1*cos(t) + h;
    y1 = r1*sin(t) + k;
    x = [x2 fliplr(x1) x1(1)];
    y = [y2 fliplr(y1) y1(1)];
    
    zup = ones(1,length(x)).*(cylP.z0 + cylP.offset(3)+cylP.dz/2);
    
    Pup = patch(x,y,zup,col,'parent',thg);
    
    zdn = ones(1,length(x)).*(cylP.z0 + cylP.offset(3)-cylP.dz/2);
    
    Pdn = patch(x,y,zdn,col,'parent',thg);
    
    x = [x2 fliplr(x2) x2(1)];
    y = [y2 fliplr(y2) y2(1)];
    zup = ones(1,length(x2)).*(cylP.z0 + cylP.offset(3)+cylP.dz/2);
    zdn = ones(1,length(x2)).*(cylP.z0 + cylP.offset(3)-cylP.dz/2);
    z = [zup fliplr(zdn) zup(1)];
    Pout = patch(x,y,z,col,'parent',thg);
    
    x = [x1 fliplr(x1) x1(1)];
    y = [y1 fliplr(y1) y1(1)];
    
    Pin = patch(x,y,z,col,'parent',thg);
    
    
    x1 = linspace(r1,r2,10) * cos(t(1)) + h;
    y1 = linspace(r1,r2,10) * sin(t(1)) + k;
    
    x = [x1 fliplr(x1) x1(1)];
    y = [y1 fliplr(y1) y1(1)];
    Plow = patch(x,y,z,col,'parent',thg);
    
    x1 = linspace(r1,r2,10) * cos(t(end)) + h;
    y1 = linspace(r1,r2,10) * sin(t(end)) + k;
    x = [x1 fliplr(x1) x1(1)];
    y = [y1 fliplr(y1) y1(1)];
    Phigh = patch(x,y,z,col,'parent',thg);
    
    P = struct();
    P.P1 = Pup;
    P.P2 = Pdn;
    P.P3 = Pout;
    P.P4 = Pin;
    P.P5 = Plow;
    P.P6 = Phigh;
    
    
    
end