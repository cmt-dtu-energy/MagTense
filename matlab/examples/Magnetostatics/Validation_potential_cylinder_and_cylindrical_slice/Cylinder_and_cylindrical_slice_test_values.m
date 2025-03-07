function Cylinder_and_cylindrical_slice_test_values

clearvars
close all

figure1= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig1 = axes('Parent',figure1,'Layer','top','FontSize',16); hold on; grid on; box on

slice = 1;

info.use_matlab = 0;

%--- To be logical we define phi1 to always be smaller than phi2.
%--- If this is not the case, we should substract 2*pi from phi1.

%--- Cylinder information
if (slice == 0)
    Ro = 0.25;
    Height = 0.7;
    Ri = 0;
    phi1 = 0;
    phi2 = 2*pi;
end

if (slice == 1)
    % Ro = 0.35;
    % Ri = 0.25;
    % Height = 0.7;
    % phi1 = pi/7;
    % phi2 = pi/3;
    % file = '1';

    Ro = 0.35;
    Ri = 0.25;
    Height = 0.7;
    phi1 = pi/7;
    phi2 = 2*pi-pi/3;
    file = '2';

    % Ro = 0.35;
    % Ri = 0.25;
    % Height = 0.7;
    % phi1 = -pi/4;
    % phi2 = pi+pi/3;
    % file = '3';

    % Ro = 0.35;
    % Ri = 0.25;
    % Height = 0.7;
    % phi1 = -pi/32;
    % phi2 =  pi/32;
    % file = '4';

    % Ro = 0.35;
    % Ri = 0.25;
    % Height = 0.7;
    % phi1 = 0;
    % phi2 = 2*pi;
    % file = '5';

    % Ro = 0.35;
    % Ri = 0.25;
    % Height = 0.7;
    % phi1 = 0;
    % phi2 = pi;
    % file = '6';
end

z1 = -Height/2;
z2 =  Height/2;

%--- Evaluation points
xpoint_org_sta = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
xpoint_org_end = [2 0 0 2 2 0 1.5 -1.5 1.5 1.5 -1.5 -1.5 1.5 -1.5 -2 0 0];
ypoint_org_sta = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.3 0.2];
ypoint_org_end = [0 2 0 2 0 2 0.8 0.8 -0.8 0.8 -0.8 0.8 -0.8 -0.8 0 0.3 0.2];
zpoint_org_sta = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2 -2];
zpoint_org_end = [0 0 2 0 2 2 0.75 0.75 0.75 -0.75 0.75 -0.75 -0.75 -0.75 0 2 2];

%--- Magnetization directions
Mxp = [1 0 0 2 -2 2 2 -2 -2 2 -2];
Myp = [0 1 0 3 3 -3 3 -3 3 -3 -3];
Mzp = [0 0 1 4 4 4 -4 4 -4 -4 -4];

colorarr = turbo(length(xpoint_org_sta));

for m = 4%1:length(Mxp)

    
    for j = 7%1:length(xpoint_org_sta)
        xpoint = [linspace(xpoint_org_sta(j),xpoint_org_end(j),20)];
        ypoint = [linspace(ypoint_org_sta(j),ypoint_org_end(j),20)];
        zpoint = [linspace(zpoint_org_sta(j),zpoint_org_end(j),20)];

        if (j == 15)
            %--- Case 15
            xpoint = [linspace(xpoint_org_sta(j),xpoint_org_end(j),20) -0.35];
            ypoint = [linspace(ypoint_org_sta(j),ypoint_org_end(j),20) 0];
            zpoint = [linspace(zpoint_org_sta(j),zpoint_org_end(j),20) 0];
        end

        if (j == 16)
            %--- Case 16
            xpoint = [linspace(xpoint_org_sta(j),xpoint_org_end(j),20) 0];
            ypoint = [linspace(ypoint_org_sta(j),ypoint_org_end(j),20) 0.3];
            zpoint = [linspace(zpoint_org_sta(j),zpoint_org_end(j),20) Height/2];
        end
        
        if (j == 17)
            %--- Case 17
            xpoint = [linspace(xpoint_org_sta(j),xpoint_org_end(j),20) 0];
            ypoint = [linspace(ypoint_org_sta(j),ypoint_org_end(j),20) 0.2];
            zpoint = [linspace(zpoint_org_sta(j),zpoint_org_end(j),20) Height/2];
        end

        for i = 1:length(xpoint)
            [phi_arr(i),r_arr(i),z_arr(i)] = cart2pol(xpoint(i),ypoint(i),zpoint(i));
        end
        
        [Psi_tot, Psi_z_arr, Psi_r_arr, Psi_phi_arr] = MagTense_Validation_cylinder_and_clyndrical_slice_potential(Ro, Ri, phi1, phi2, z1, z2, [Mxp(m) Myp(m) Mzp(m)], r_arr, phi_arr, z_arr, slice, info);
    
    end
    disp('---------------')
    
    if (slice == 1)
        MagTense_path = '..\..\..\..\..\MagTense\documentation\examples_FEM_validation\Validation_potential_cylinder_slice';
        data = load([MagTense_path '\Line' sprintf('%02.2i',j) '_' file '.txt']);
    end
    if (slice == 0)
        MagTense_path = '..\..\..\..\..\MagTense\documentation\examples_FEM_validation\Validation_potential_cylinder';
        data = load([MagTense_path '\Line' sprintf('%02.2i',j) '.txt']);
    end
   
    plot(data(:,1),data(:,1+m),'.','markersize',3,'Color',colorarr(j,:));
    hold on;  
    plot(sqrt((xpoint-xpoint_org_sta(j)).^2+(ypoint-ypoint_org_sta(j)).^2+(zpoint-zpoint_org_sta(j)).^2),Psi_tot(:,1),'p','MarkerSize',10,'Color',colorarr(j,:))
    plot(sqrt((xpoint-xpoint_org_sta(j)).^2+(ypoint-ypoint_org_sta(j)).^2+(zpoint-zpoint_org_sta(j)).^2),Psi_tot(:,3),'o','MarkerSize',10,'Color',colorarr(j,:))
end
end