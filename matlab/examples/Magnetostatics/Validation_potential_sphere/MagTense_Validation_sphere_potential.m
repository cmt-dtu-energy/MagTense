function MagTense_Validation_sphere_potential
clearvars
close all

figure1= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig1 = axes('Parent',figure1,'Layer','top','FontSize',16); hold on; grid on; box on

%--- Sphere dimensions
Radius = 2.5;

%--- Evaluation points
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
    
            r = sqrt(x.^2+y.^2+z.^2);
    
            phi_2(i) =  1/(4*pi)*4/3*pi*(2*Radius/(r+Radius+abs(r-Radius)))^3*dot([x y z],[Mxp(m) Myp(m) Mzp(m)]);
            
        end
    disp('---------------')
    
    MagTense_path = '..\..\..\..\..\MagTense\documentation\examples_FEM_validation\Validation_potential_sphere';
    data = load([MagTense_path '\Line' sprintf('%02.2i',j) '.txt']);
    plot(fig1,data(:,1),data(:,1+m),'-');
    
    plot(fig1,sqrt(xpoint.^2+ypoint.^2+zpoint.^2),phi_2,'o')
    xlabel(fig1,'Distance [m]')
    ylabel(fig1,'Potential [A]')
    
    data_interp = interp1(data(:,1),data(:,1+m),sqrt(xpoint.^2+ypoint.^2+zpoint.^2));
    
    plot(fig2,sqrt(xpoint.^2+ypoint.^2+zpoint.^2),(data_interp-phi_2)./phi_2*100,'.','markersize',20);
    xlabel(fig2,'Distance [m]')
    ylabel(fig2,'Percentage error [%]')
    end
end

end

function value = F(x,y,z) 
    r = sqrt(x.^2+y.^2+z.^2);

    if (x == 0)
        term1 = 0;
    else
        term1 = -x.*atan(y.*z./(x.*r));
    end

    if (r == 0)
        value = 0;
    else
        value = term1 + y.*log(z+r) + z.*log(y+r);
    end
end