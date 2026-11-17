%%This function compares MagTense to a FEM simulation for a cone.
%%A cone tile is a right circular cone that stands on its base, and MagTense
%%evaluates it with the closed-form demagnetization tensor of R. Bjoerk, "The
%%magnetic field of a homogeneously magnetized cone". The cone has slant length
%%1.3 m and half-opening angle 0.6 rad, i.e. base radius 0.734 m and height
%%1.073 m, and carries the magnetization (0.32, 0.74, 0.89) A/m. In the FEM
%%model its apex is at the origin and its base in the plane z = height, so the
%%tile is turned upside down (rotation by pi about the x-axis) and the centre of
%%its base is moved to (0, 0, height). The field is compared with a COMSOL
%%solution along the line from (-0.6, -0.4, -0.1) to (1.1, 0.7, 1.3), which runs
%%through the magnet.
function [rel_int_error] = MagTense_Validation_cone()

%make sure to source the right path for the generic Matlab routines
addpath(genpath('../../../util/'));
addpath('../../../MEX_files/');

slant = 1.3;
half_angle = 0.6;
radius = slant*sin(half_angle);
height = slant*cos(half_angle);
M = [0.32, 0.74, 0.89];
line_start = [-0.6, -0.4, -0.1];

%%Get a default tile from MagTense
tile = DefaultMagTile();

%ensure the tile is a permanent magnet
tile = tile.setMagnetType('hard');

%set the geometry to be a cone
tile = tile.setMagTileType('cone');

%set the dimensions: abc = [base radius, 0, height]. The second entry is
%reserved for the top radius of a truncated cone and must be zero
tile.abc = [radius, 0, height];

%The tile has the centre of its base at offset and its apex at the distance
%height along the rotated z-axis. Turned upside down and moved up by its height,
%the apex is at the origin and the base in the plane z = height
tile.rotAngles = [pi, 0, 0];
tile.offset = [0, 0, height];

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

%%Let MagTense find the magnetization vector of the cone by iterating to a
%%self-consistent solution (trivial here, as mu_r = 1)
tile = IterateMagnetization( tile, [], [], 1e-6, 100 );

%%Load comparison data from the FEM simulation. Each file holds one component of
%%H in A/m against the x, y, z coordinates of the mesh nodes on the line
data_FEM{1} = load('../../../../documentation/examples_FEM_validation/Validation_cone/Validation_cone_Hx.txt');
data_FEM{2} = load('../../../../documentation/examples_FEM_validation/Validation_cone/Validation_cone_Hy.txt');
data_FEM{3} = load('../../../../documentation/examples_FEM_validation/Validation_cone/Validation_cone_Hz.txt');

%In the FEM model the air around the cone is a sphere of radius 1.5 m about the
%middle of the cone axis, enclosed by an infinite element domain. The line leaves
%the sphere 2.4506 m from its start, 11 mm before its end. The coordinates are
%scaled in the infinite element domain, so the nodes there cannot be compared
%and the comparison stops at 2.45 m
line_length = 2.45;

%Make a figure
figure1 = figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto');
fig1 = axes('Parent',figure1,'Layer','top','FontSize',16);
hold all
grid on
box on
colors = {'r','g','b'};

for i = 1:3
    %The export lists the nodes in mesh order, so they are sorted along the line
    [FEM_dist, order] = sort(sqrt(sum((data_FEM{i}(:,1:3) - line_start).^2,2)));
    order = order(FEM_dist < line_length);
    FEM_dist = FEM_dist(FEM_dist < line_length);
    pts = data_FEM{i}(order,1:3);
    H_FEM = data_FEM{i}(order,4);

    %%Now find the field in the points
    H = getHFromTiles_mex( tile, pts, int32( length(tile) ), int32( length(pts(:,1)) ) );

    %Plot the solution against the distance along the line
    plot(FEM_dist,H(:,i),[colors{i} '.'],'linewidth',2);
    plot(FEM_dist,H_FEM,[colors{i} 'o'],'linewidth',2);

    %A node on the surface of the cone is listed twice, with the field on each
    %side of the surface. The field jumps there, and the surface of the FEM mesh
    %lies about 0.2 mm inside the cone, so MagTense returns the field inside the
    %cone for both. These two nodes are left out of the error
    on_surface = [diff(FEM_dist) < 1e-12; false];
    on_surface = on_surface | [false; on_surface(1:end-1)];

    % Calculate the relative error in percent between the MagTense and the FEM
    % solution
    rel_int_error(i) = calculate_relative_integral_error(FEM_dist(~on_surface),H_FEM(~on_surface),FEM_dist(~on_surface),H(~on_surface,i));
end

h_l = legend('MagTense, H_x','FEM, H_x','MagTense, H_y','FEM, H_y','MagTense, H_z','FEM, H_z','Location','NorthEast');
set(h_l,'fontsize',10);
ylabel('H_i [A m^{-1}]');
xlabel('Distance along the line from (-0.6,-0.4,-0.1) to (1.1,0.7,1.3) [m]');
title('Cone, M = (0.32, 0.74, 0.89) A/m');

end
