%%This function compares MagTense to a FEM simulation for a full cylinder.
%%A cylindrical tile that spans the full circle (dtheta = 2*pi) with inner radius zero
%%is a full cylinder, and MagTense evaluates it with the closed-form field of
%%Caciagli et al., J. Magn. Magn. Mater. 456 (2018) 423, instead of the
%%cylinder-piece integrals. The cylinder has radius 1.1 m and height 0.75 m, is
%%centred at the origin and carries the magnetization (2, 3, 4) A/m, i.e. both an
%%axial and a transverse component. The field is compared with a COMSOL solution
%%along the line from the origin to (-1.5, -0.75, 1.25), which runs through the
%%magnet and leaves it through the top end surface.
function [rel_int_error] = MagTense_Validation_cylinder()

%make sure to source the right path for the generic Matlab routines
addpath(genpath('../../../util/'));
addpath('../../../MEX_files/');

radius = 1.1;
height = 0.75;
M = [2, 3, 4];

%%Get a default tile from MagTense
tile = DefaultMagTile();

%ensure the tile is a permanent magnet
tile = tile.setMagnetType('hard');

%set the geometry to be a cylindrical tile
tile = tile.setMagTileType('cylinder');

%set the dimensions: the full circle with inner radius r0 - dr/2 = 0 makes the
%tile a full cylinder
tile.r0 = radius/2;
tile.theta0 = 0;
tile.z0 = 0;

tile.dr = radius;
tile.dtheta = 2*pi;
tile.dz = height;

%set the center position of the cylinder (centered at Origo)
tile.offset = [0, 0, 0];

%no rotation (a full cylinder follows rotAngles like the other tile types)
tile.rotAngles = [0,0,0];

%set the easy axis along the magnetization. With mu_r = 1 the magnetization is
%the remanence along the easy axis, i.e. M
tile.u_ea = M ./ norm(M);

%ensure the two hard axes are perpendicular and normalized
tile.u_oa1 = cross(tile.u_ea, [1, 0, 0]);
tile.u_oa1 = tile.u_oa1 ./ norm(tile.u_oa1);
tile.u_oa2 = cross(tile.u_ea, tile.u_oa1);

%set the relative permeability for the easy axis
tile.mu_r_ea = 1.00;
%and for the two hard axes
tile.mu_r_oa = 1.00;

%set the remanence of the magnet in A/m
tile.Mrem = norm(M);

%%Let MagTense find the magnetization vector of the cylinder by iterating to a
%%self-consistent solution (trivial here, as mu_r = 1)
tile = IterateMagnetization( tile, [], [], 1e-6, 100 );

%%Load comparison data from the FEM simulation. Each file holds one component of
%%H in A/m against the x, y, z coordinates of the evaluation points
data_FEM{1} = load('../../../../documentation/examples_FEM_validation/Validation_cylinder/Validation_cylinder_full_Hx.txt');
data_FEM{2} = load('../../../../documentation/examples_FEM_validation/Validation_cylinder/Validation_cylinder_full_Hy.txt');
data_FEM{3} = load('../../../../documentation/examples_FEM_validation/Validation_cylinder/Validation_cylinder_full_Hz.txt');

%Make a figure
figure1 = figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto');
fig1 = axes('Parent',figure1,'Layer','top','FontSize',16);
hold all
grid on
box on
colors = {'r','g','b'};
comp = {'H_x','H_y','H_z'};

for i = 1:3
    %The export lists the nodes in mesh order, and the node on the top end
    %surface twice, with the value from each side of the surface. The points
    %are evaluated as given (MagTense takes the limit from outside on a
    %surface) and sorted along the line for the integrated error
    [FEM_dist, order] = sort(sqrt(sum(data_FEM{i}(:,1:3).^2,2)));
    pts = data_FEM{i}(order,1:3);
    H_FEM = data_FEM{i}(order,4);

    %%Now find the field in the points
    H = getHFromTiles_mex( tile, pts, int32( length(tile) ), int32( length(pts(:,1)) ) );

    %Plot the solution against the distance from the origin along the line
    plot(FEM_dist,H(:,i),[colors{i} '.'],'linewidth',2);
    plot(FEM_dist,H_FEM,[colors{i} 'o'],'linewidth',2);

    % Calculate the relative error in percent between the MagTense and the FEM
    % solution. The MagTense curve is interpolated to the FEM points, so the
    % duplicated node is passed once (its two MagTense values are identical)
    [~, iu] = unique(FEM_dist);
    rel_int_error(i) = calculate_relative_integral_error(FEM_dist,H_FEM,FEM_dist(iu),H(iu,i));
end

h_l = legend('MagTense, H_x','FEM, H_x','MagTense, H_y','FEM, H_y','MagTense, H_z','FEM, H_z','Location','NorthEast');
set(h_l,'fontsize',10);
ylabel('H_i [A m^{-1}]');
xlabel('Distance from the origin along the line to (-1.5,-0.75,1.25) [m]');
title('Full cylinder, M = (2, 3, 4) A/m');

end
