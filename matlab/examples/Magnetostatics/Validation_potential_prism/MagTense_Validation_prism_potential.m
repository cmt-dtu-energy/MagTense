function MagTense_Validation_prism_potential
clearvars
close all

figure1= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig1 = axes('Parent',figure1,'Layer','top','FontSize',16); hold on; grid on; box on

%--- Prism dimensions
a = 1; %--- The half side length
b = 2; %--- The half side length
c = 3; %--- The half side length

%--- Evaluation point
xpoint_org = [10 0 0 10 10 0 8 -8 8 8 -8 -8 8 -8];
ypoint_org = [0 10 0 10 0 10 6 6 -6 6 -6 6 -6 -6];
zpoint_org = [0 0 10 0 10 10 9 9 9 -9 9 -9 -9 -9];

Mxp = [1 0 0 2 -2 2 2 -2 -2 2 -2];
Myp = [0 1 0 3 3 -3 3 -3 3 -3 -3];
Mzp = [0 0 1 4 4 4 -4 4 -4 -4 -4];

for m = 1:length(Mxp)
    figure2= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig2 = axes('Parent',figure2,'Layer','top','FontSize',16); hold on; grid on; box on

    for j = 1:length(xpoint_org)   
        xpoint = [(0:0.01:1)*xpoint_org(j)];
        ypoint = [(0:0.01:1)*ypoint_org(j)];
        zpoint = [(0:0.01:1)*zpoint_org(j)];
    
        for i = 1:length(xpoint)
            x = xpoint(i);
            y = ypoint(i);
            z = zpoint(i);
    
            %--- x'-faces
            % phi_sx1 = -1/(4*pi)*Mxp(m)*(-(F(x-a,y-b,z+c)-F(x-a,y+b,z+c))+(F(x-a,y-b,z-c)-F(x-a,y+b,z-c)));
            % phi_sx2 =  1/(4*pi)*Mxp(m)*(-(F(x+a,y-b,z+c)-F(x+a,y+b,z+c))+(F(x+a,y-b,z-c)-F(x+a,y+b,z-c)));
            % 
            % phi_sx(i) = phi_sx1+phi_sx2;
                    
            phi_sx1_2 = -1/(4*pi)*Mxp(m)*((F(x-a,y-b,z-c) - F(x-a,y-b,z+c)) - (F(x-a,y+b,z-c) - F(x-a,y+b,z+c)));
            phi_sx2_2 =  1/(4*pi)*Mxp(m)*((F(x+a,y-b,z-c) - F(x+a,y-b,z+c)) - (F(x+a,y+b,z-c) - F(x+a,y+b,z+c)));
                          
            phi_sx_2(i) = phi_sx1_2+phi_sx2_2;
    
            % phi_sx1_3 = -1/(4*pi)*Mxp(m)*((F_cart(+a,+b,+c,x,y,z)-F_cart(+a,+b,-c,x,y,z))-(F_cart(+a,-b,+c,x,y,z)-F_cart(+a,-b,-c,x,y,z)));
            % phi_sx2_3 =  1/(4*pi)*Mxp(m)*((F_cart(-a,+b,+c,x,y,z)-F_cart(-a,+b,-c,x,y,z))-(F_cart(-a,-b,+c,x,y,z)-F_cart(-a,-b,-c,x,y,z)));
            % 
            % phi_sx_3(i) = phi_sx1_3+phi_sx2_3;
    
            % % %--- Double integral
            % phi_sx_int1 =  1/(4*pi)*Mxp(m)*integral2(@(yp,zp) 1./sqrt((x-a).^2+(y-yp).^2+(z-zp).^2), -b, b, -c, c);
            % phi_sx_int2 = -1/(4*pi)*Mxp(m)*integral2(@(yp,zp) 1./sqrt((x+a).^2+(y-yp).^2+(z-zp).^2), -b, b, -c, c);
            % 
            % phi_sx_int(i) = phi_sx_int1+phi_sx_int2;
    
    
            %--- y'-faces              
            phi_sy1_2 = -1/(4*pi)*Myp(m)*((F(y-b,z-c,x-a) - F(y-b,z-c,x+a)) - (F(y-b,z+c,x-a) - F(y-b,z+c,x+a)));
            phi_sy2_2 =  1/(4*pi)*Myp(m)*((F(y+b,z-c,x-a) - F(y+b,z-c,x+a)) - (F(y+b,z+c,x-a) - F(y+b,z+c,x+a)));
                          
            phi_sy_2(i) = phi_sy1_2+phi_sy2_2;
    
            % % %--- Double integral
            % phi_sy_int1 =  1/(4*pi)*Myp(m)*integral2(@(xp,zp) 1./sqrt((x-xp).^2+(y-b).^2+(z-zp).^2), -a, a, -c, c);
            % phi_sy_int2 = -1/(4*pi)*Myp(m)*integral2(@(xp,zp) 1./sqrt((x+xp).^2+(y+b).^2+(z-zp).^2), -a, a, -c, c);
            % 
            % phi_sy_int(i) = phi_sy_int1+phi_sy_int2;
            
    
            %--- z'-faces              
            phi_sz1_2 = -1/(4*pi)*Mzp(m)*((F(z-c,x-a,y-b) - F(z-c,x-a,y+b)) - (F(z-c,x+a,y-b) - F(z-c,x+a,y+b)));
            phi_sz2_2 =  1/(4*pi)*Mzp(m)*((F(z+c,x-a,y-b) - F(z+c,x-a,y+b)) - (F(z+c,x+a,y-b) - F(z+c,x+a,y+b)));
                          
            phi_sz_2(i) = phi_sz1_2+phi_sz2_2;
    
            % % %--- Double integral
            % phi_sz_int1 =  1/(4*pi)*Mzp(m)*integral2(@(xp,yp) 1./sqrt((x-xp).^2+(y-yp).^2+(z-c).^2), -a, a, -b, b);
            % phi_sz_int2 = -1/(4*pi)*Mzp(m)*integral2(@(xp,yp) 1./sqrt((x+xp).^2+(y-yp).^2+(z+c).^2), -a, a, -b, b);
            % 
            % phi_sz_int(i) = phi_sz_int1+phi_sz_int2;
    
            phi_2(i) = phi_sx_2(i) + phi_sy_2(i) + phi_sz_2(i);
    
            % phi_int(i) = phi_sx_int(i) + phi_sy_int(i) + phi_sz_int(i);
    
        end
        disp('---------------')
        
        MagTense_path = '..\..\..\..\..\MagTense\documentation\examples_FEM_validation\Validation_potential_prism';
        data = load([MagTense_path '\Line' sprintf('%02.2i',j) '.txt']);
        plot(fig1,data(:,1),data(:,1+m),'-');
        
        % plot(sqrt(xpoint.^2+ypoint.^2+zpoint.^2),phi_sx_int,'d')
        plot(fig1,sqrt(xpoint.^2+ypoint.^2+zpoint.^2),phi_2,'o')
        
        data_interp = interp1(data(:,1),data(:,1+m),sqrt(xpoint.^2+ypoint.^2+zpoint.^2));
        
        plot(fig2,sqrt(xpoint.^2+ypoint.^2+zpoint.^2),(data_interp-phi_2)./phi_2*100,'.','markersize',20);
        
    end
end

end

function value = F(x,y,z) 
    r = sqrt(x.^2+y.^2+z.^2);

    if (x == 0)
        term1 = 0;
    else
        term1 = x.*atan(y.*z./(x.*r));
    end

    if (r == 0)
        value = 0;
    else
        value = term1-y.*log(z+r)-z.*log(y+r);
    end
end

function value = F_cart(xp,yp,zp,x,y,z)
    r_rp = sqrt((x-xp).^2+(y-yp).^2+(z-zp).^2);
    value = (x-xp)*atan(((y-yp)*(z-zp))/((x-xp)*r_rp)) - (y-yp)*log((z-zp) + r_rp) - (z-zp)*log((y-yp) + r_rp);
end