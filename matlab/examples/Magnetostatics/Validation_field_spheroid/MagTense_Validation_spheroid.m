function MagTense_Validation_spheroid()

clearvars
close all

%--- Choose a spheroid to compare against
use_existing_FEM_prolate = 1;
use_existing_FEM_oblate  = 0;
%--- If both variables are zero, a random spheroid will be computed using
%--- Comsol and compared with MagTensen
%--- This requires that Comsol is installed

%make sure to source the right path for the generic Matlab routines
addpath(genpath('../../../util/'));
addpath('../../../MEX_files/');

%% Geometric parameters
%%Get a default tile from MagTense
tile = getDefaultMagTile();

%ensure the tile is a permanent magnet
tile.magnetType = getMagnetType('hard');

%set the geometry to be a rectangular prism
tile.tileType = getMagTileType('spheroid');

%use existing FEM data or run a random model in Comsol to compare against
if use_existing_FEM_prolate
    % Prolate
    tile.abc = [0.00611061	 0.00305531	 0.00305531];
    tile.offset = [-0.05781690, 0.00804030, 0.01479210];
    
    %--- The axis specified as a vector, as done in Comsol Axis type: Cartesian
    axis_vector  = [0.101152 -0.0385283 -0.225353];
    [phi, theta] = cart2sph(axis_vector(1),axis_vector(2),axis_vector(3));
    %--- The corresponding angles in Comsol
    tile.ax   =  [pi/2-theta, phi, 0];
    
    tile.rot_axis = 'c';       
    
    do_Comsol_model = 0;

elseif use_existing_FEM_oblate
    % Oblate
    tile.abc = [0.00611061	 0.00611061	 0.00305531];
    tile.offset = [-0.05781690, 0.00804030, 0.01479210];

    %--- The axis specified as a vector, as done in Comsol Axis type: Cartesian
    axis_vector  = [0.101152 -0.0385283 -0.225353];
    [phi, theta] = cart2sph(axis_vector(1),axis_vector(2),axis_vector(3));
    %--- The corresponding angles in Comsol
    tile.ax   =  [pi/2-theta, phi, 0];

    tile.rot_axis = 'c';  
    
    do_Comsol_model = 0;
else
    %--- Do a random model and compute the corresponding Comsol model
    ab = rand(1,2);
    c = ab(randi(2));
    abc = [ab c];
    tile.abc = abc(randperm(3));
    
    tile.offset = -1 + 2*rand(1,3);
    tile.ax  = [1,1,1].*(-1 + 2*rand(1,3)) ;
        
    if (rand(1) < 0.5)
        tile.rot_axis = 'c';
    else
        tile.rot_axis = 'symm';
    end
    
    do_Comsol_model = 1;
end

%% Magnetic Parameters (not implemented in the COMSOL validation model)
%% All with respect to the global Cartesian coordinate system (x,y,z)
%% Therefore they have nothing to do with the geometrical orientation of the spheroid
%% can be implemented as for any other shape
%set the easy axis of the cube. This is expected to be with respect to the global
%coordinate system basis
mu0 = 4*pi*1e-7;

tile.u_ea = [1.1, 0.5, 0.3];
tile.u_ea = tile.u_ea ./ norm(tile.u_ea);

%ensure the two hard axes are perpendicular and normalized
tile.u_oa1 = [1.1, -0.5,  0.3];
tile.u_oa1 = tile.u_oa1 ./ norm(tile.u_oa1);

tile.u_oa2 = [-1.1,  0.5, -0.3 ];
tile.u_oa2 = tile.u_oa2 ./ norm(tile.u_oa2);

%set the relative permeability for the easy axis
tile.mu_r_ea = 1.00;
%and for the two hard axes
tile.mu_r_oa = 1.00;

%set the remanence of the magnet (1.2 T converted to A/m)
tile.Mrem = 1.2 / mu0;


%% Comsol Computation (only to compare the geometry for now)
if (do_Comsol_model == 1)
    [~] = getSpheroidComsolCompare( tile );
end

if use_existing_FEM_prolate
    data_FEM_coor = load('..\..\..\..\documentation\examples_FEM_validation\Validation_spheroid\Validation_spheroid_prolate_surface_coordinates.txt');
elseif use_existing_FEM_oblate
    data_FEM_coor = load('..\..\..\..\documentation\examples_FEM_validation\Validation_spheroid\Validation_spheroid_oblate_surface_coordinates.txt');
else
    data_FEM_coor = load('Comsol_spheroid_surface.txt');
end


%% Get the rotation angles of the specified axis 
[rotAngles] = getSpheroidRotAngles( tile );
tile.rotAngles = rotAngles;


%% Mag Tense computation
%%Let MagTense find the magnetization vector of the cube by iterating to a
%%self-consistent solution. The two empty arguments ( [] ) ensure that
%%default values are used and they are in any case not relevant for this
%%example. The 1e-6 argument is the maximum relative error (change in
%%magnetization from one iteration to the next) allowed before convergence
%%is reached. Finally, 100 reflects the max. no. of allowed iterations
tile = IterateMagnetization( tile, [], [], 1e-6, 100 );


%% Plot a visual representation of the spheroid
%--- MagTense
figure0 = figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig0 = axes('Parent',figure0,'Layer','top','FontSize',16); hold all; grid on; box on; axis equal; light('position',[1,1,1]); view(30,30);

%--- Find the field in a set of points around the spheroid
TheMaxLenght = sqrt(sum(tile.abc.^2)) ;    
Np = 101 ;
x = tile.offset(1) + (1.2*TheMaxLenght).*linspace(-1,+1,Np) ;
y = tile.offset(2) + (1.2*TheMaxLenght).*linspace(-1,+1,Np) ;
z = tile.offset(3) + (1.2*TheMaxLenght).*linspace(-1,+1,Np) ;

[X,Y,Z] = ndgrid(x,y,z);
pts = [reshape(X,numel(X),1,1) reshape(Y,numel(Y),1,1) reshape(Z,numel(Z),1,1)];

%--- Get the field
H = getHFromTiles_mex( tile, pts, int32( length(tile) ), int32( length(pts(:,1)) ) );

%--- Find the spheroid based on where the field is uniform, and thus equal
%--- to the mode of the field values
Hx = H(:,1) ; Hy = H(:,2) ; Hz = H(:,3) ;
HHx = reshape(Hx,size(X)) ; HHy = reshape(Hy,size(X)) ; HHz = reshape(Hz,size(X)) ;
MagTenseSpheroidSurface = abs(HHx-mode(Hx)) ;
[X2,Y2,Z2] = meshgrid(x,y,z);                                        % isosurface doesn't like ndgrid unfortunately
MagTenseSpheroidSurface = permute(MagTenseSpheroidSurface,[2,1,3]) ; % isosurface doesn't like ndgrid unfortunately
isosurface(X2,Y2,Z2,MagTenseSpheroidSurface,1e-4);

%--- Comsol FEM results
TheseIndexes = randperm(size(data_FEM_coor,1)) ;
TheseIndexes = TheseIndexes(1:round(size(data_FEM_coor,1)/3)) ;
plot3(data_FEM_coor(TheseIndexes,1), data_FEM_coor(TheseIndexes,2),data_FEM_coor(TheseIndexes,3),'.r') % Comsol

title(fig0,'Spheroid')
xlabel(fig0,'x');
ylabel(fig0,'y');
zlabel(fig0,'z');
h_l = legend(fig0,'MagTense spheroid','Comsol spheroid','Location','NorthEast');
set(h_l,'Fontsize',10);

%% Compare field with FEM data
%--- Compare the MagTense field Hx, Hy and Hz along the x-, y-, and z-axis
%--- passing though the center of the spheroid

%Make a figure
figure1 = figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig1 = axes('Parent',figure1,'Layer','top','FontSize',16); hold all; grid on; box on
figure2 = figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig2 = axes('Parent',figure2,'Layer','top','FontSize',16); hold all; grid on; box on
figure3 = figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig3 = axes('Parent',figure3,'Layer','top','FontSize',16); hold all; grid on; box on

%Load comparison data from FEM simulation
if use_existing_FEM_oblate
    FEM_str = 'oblate';
end
if use_existing_FEM_prolate
    FEM_str = 'prolate';
end
if (use_existing_FEM_oblate || use_existing_FEM_prolate)
    data_FEM_Hx_x = load(['..\..\..\..\documentation\examples_FEM_validation\Validation_spheroid\Validation_spheroid_' FEM_str '_Hx_x.txt']);
    data_FEM_Hx_y = load(['..\..\..\..\documentation\examples_FEM_validation\Validation_spheroid\Validation_spheroid_' FEM_str '_Hx_y.txt']);
    data_FEM_Hx_z = load(['..\..\..\..\documentation\examples_FEM_validation\Validation_spheroid\Validation_spheroid_' FEM_str '_Hx_z.txt']);
    data_FEM_Hy_x = load(['..\..\..\..\documentation\examples_FEM_validation\Validation_spheroid\Validation_spheroid_' FEM_str '_Hy_x.txt']);
    data_FEM_Hy_y = load(['..\..\..\..\documentation\examples_FEM_validation\Validation_spheroid\Validation_spheroid_' FEM_str '_Hy_y.txt']);
    data_FEM_Hy_z = load(['..\..\..\..\documentation\examples_FEM_validation\Validation_spheroid\Validation_spheroid_' FEM_str '_Hy_z.txt']);
    data_FEM_Hz_x = load(['..\..\..\..\documentation\examples_FEM_validation\Validation_spheroid\Validation_spheroid_' FEM_str '_Hz_x.txt']);
    data_FEM_Hz_y = load(['..\..\..\..\documentation\examples_FEM_validation\Validation_spheroid\Validation_spheroid_' FEM_str '_Hz_y.txt']);
    data_FEM_Hz_z = load(['..\..\..\..\documentation\examples_FEM_validation\Validation_spheroid\Validation_spheroid_' FEM_str '_Hz_z.txt']);
    xlim(fig1,[-0.09 -0.03])
    xlim(fig2,[-0.02 0.04])
    xlim(fig3,[-0.02 0.04])
else
    data_FEM_Hx_x = load(['Comsol_spheroid_Hx_x.txt']);
    data_FEM_Hx_y = load(['Comsol_spheroid_Hx_y.txt']);
    data_FEM_Hx_z = load(['Comsol_spheroid_Hx_z.txt']);
    data_FEM_Hy_x = load(['Comsol_spheroid_Hy_x.txt']);
    data_FEM_Hy_y = load(['Comsol_spheroid_Hy_y.txt']);
    data_FEM_Hy_z = load(['Comsol_spheroid_Hy_z.txt']);
    data_FEM_Hz_x = load(['Comsol_spheroid_Hz_x.txt']);
    data_FEM_Hz_y = load(['Comsol_spheroid_Hz_y.txt']);
    data_FEM_Hz_z = load(['Comsol_spheroid_Hz_z.txt']);
end

%--- Run a MagTense model
if (use_existing_FEM_prolate || use_existing_FEM_oblate)
    %--- Find the field in a set of points along each axis
    x = -0.13:0.0001:0.1;
    y = -0.13:0.0001:0.1;
    z = -0.13:0.0001:0.1;
else
    %--- Get the points to evaluate the field in from the FEM data
    x = data_FEM_Hx_x(:,1)';
    y = data_FEM_Hx_y(:,1)';
    z = data_FEM_Hx_z(:,1)';
end
pts = zeros( numel(x)+numel(y)+numel(z), 3 );
pts(1:numel(x),:) = [x; zeros(1,numel(x))+tile.offset(2); zeros(1,numel(x))+tile.offset(3)]';
pts((numel(x)+1):(numel(x)+numel(y)),:) = [zeros(1,numel(y))+tile.offset(1); y; zeros(1,numel(y))+tile.offset(3)]';
pts((numel(x)+1+numel(y)):(numel(x)+numel(y)+numel(z)),:) = [zeros(1,numel(z))+tile.offset(1); zeros(1,numel(z))+tile.offset(2); z]';

%get the field
H = getHFromTiles_mex( tile, pts, int32( length(tile) ), int32( length(pts(:,1)) ) );
         
%Find the norm of the field
Hnorm = squeeze( sqrt( sum(H.^2,2) ) );

%Plot the solution
plot(fig1,x,H(1:numel(x),1),'r.');
plot(fig1,x,H(1:numel(x),2),'g.');
plot(fig1,x,H(1:numel(x),3),'b.');

plot(fig2,y,H((numel(x)+1):(numel(x)+numel(y)),1),'r.');
plot(fig2,y,H((numel(x)+1):(numel(x)+numel(y)),2),'g.');
plot(fig2,y,H((numel(x)+1):(numel(x)+numel(y)),3),'b.');

plot(fig3,z,H((numel(x)+1+numel(y)):(numel(x)+numel(y)+numel(z)),1),'r.');
plot(fig3,z,H((numel(x)+1+numel(y)):(numel(x)+numel(y)+numel(z)),2),'g.');
plot(fig3,z,H((numel(x)+1+numel(y)):(numel(x)+numel(y)+numel(z)),3),'b.');


%--- Plot FEM data
plot(fig1,data_FEM_Hx_x(:,1),data_FEM_Hx_x(:,2),'ro')
plot(fig1,data_FEM_Hy_x(:,1),data_FEM_Hy_x(:,2),'go')
plot(fig1,data_FEM_Hz_x(:,1),data_FEM_Hz_x(:,2),'bo')

plot(fig2,data_FEM_Hx_y(:,1),data_FEM_Hx_y(:,2),'ro')
plot(fig2,data_FEM_Hy_y(:,1),data_FEM_Hy_y(:,2),'go')
plot(fig2,data_FEM_Hz_y(:,1),data_FEM_Hz_y(:,2),'bo')

plot(fig3,data_FEM_Hx_z(:,1),data_FEM_Hx_z(:,2),'ro')
plot(fig3,data_FEM_Hy_z(:,1),data_FEM_Hy_z(:,2),'go')
plot(fig3,data_FEM_Hz_z(:,1),data_FEM_Hz_z(:,2),'bo')

h_l = legend(fig1,'MagTense, H_x','MagTense, H_y','MagTense, H_z','FEM, H_x','FEM, H_y','FEM, H_z','Location','SouthWest');
set(h_l,'fontsize',10);
ylabel(fig1,'H [A/m]');
xlabel(fig1,'x [m]');
figure(figure1)

h_l = legend(fig2,'MagTense, H_x','MagTense, H_y','MagTense, H_z','FEM, H_x','FEM, H_y','FEM, H_z','Location','SouthEast');
set(h_l,'fontsize',10);
ylabel(fig2,'H [A/m]');
xlabel(fig2,'y [m]');
figure(figure2)

h_l = legend(fig3,'MagTense, H_x','MagTense, H_y','MagTense, H_z','FEM, H_x','FEM, H_y','FEM, H_z','Location','SouthWest');
set(h_l,'fontsize',10);
ylabel(fig3,'H [A/m]');
xlabel(fig3,'z [m]');
figure(figure3)

%--- Cleanup
if ~(use_existing_FEM_oblate || use_existing_FEM_prolate)
    delete('Comsol_spheroid_surface.txt');
    delete('Comsol_spheroid_Hx_x.txt');
    delete('Comsol_spheroid_Hx_y.txt');
    delete('Comsol_spheroid_Hx_z.txt');
    delete('Comsol_spheroid_Hy_x.txt');
    delete('Comsol_spheroid_Hy_y.txt');
    delete('Comsol_spheroid_Hy_z.txt');
    delete('Comsol_spheroid_Hz_x.txt');
    delete('Comsol_spheroid_Hz_y.txt');
    delete('Comsol_spheroid_Hz_z.txt');
end

end