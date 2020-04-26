
%%Two PM cubes and a single Fe cube
function [] = MagTense_Example003_PM_Fe_cubes()
%close all
%figure1= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig1 = axes('Parent',figure1,'Layer','top','FontSize',16); hold on; grid on; box on       
%figure2= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig2 = axes('Parent',figure2,'Layer','top','FontSize',16); hold on; grid on; box on
%figure4= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig4 = axes('Parent',figure4,'Layer','top','FontSize',16); hold on; grid on; box on
%figure5= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig5 = axes('Parent',figure5,'Layer','top','FontSize',16); hold on; grid on; box on

%make sure to source the right path for the generic Matlab routines
addpath(genpath('../util/'));
addpath('../../MEX_files/');
%define the vacuum permeability
mu0 = 4*pi*1e-7;

%%Get a default tile from MagTense
tile = getDefaultMagTile();
    
%ensure the tile is a permanent magnet
tile.magnetType = getMagnetType('hard');

%set the geometry to be a rectangular prism
tile.tileType = getMagTileType('prism');

%set the dimensions of the prism (5 mm on either side)
d = 0.005;
tile.abc = [d,d,d];

%set the center position of the prism (centered at Origin)
tile.offset = [0,0,0];

%set the easy axis of the cube. This is expected to be with respect to the global
%coordinate system basis
tile.u_ea = [0,1,0];
%ensure the two hard axes are perpendicular and normalized
tile.u_oa1 = [1,0,0];
tile.u_oa2 = [0,0,1];

%set the relative permeability for the easy axis (NdFeB magnet)
tile.mu_r_ea = 1.0;
%and for the two hard axes
tile.mu_r_oa = 1.0;

%set the remanence of the magnet (1.2 T converted to A/m)
tile.Mrem = 1.2 / mu0;

tiles(1) = tile;
tiles(2) = tile;
tiles(3) = tile;
tiles(3).offset(1) = 2*d;
tiles(3).offset(2) = 2*d;
tiles(3).u_ea = [1,0,0];
tiles(3).u_oa1 = [0,1,0];


tiles(2).magnetType = getMagnetType('soft_const_mur');
tiles(2).offset(2) = d;
tiles(2).Mrem = 0;
tiles(2).M = [0,0,0];
tiles(2).mu_r_ea = 1000;
tiles(2).mu_r_oa = 1000;
tiles(2).color = [0,0,1];
tiles(2).u_ea = [1,0,0];
tiles(2).u_oa1 = [0,1,0];
tiles(2).u_oa2 = [0,0,1];

tiles(4) = tiles(2);
tiles(4).offset(1) = d;
tiles(4).offset(2) = 0;

res = struct('nx',10,'ny',10,'nz',1);

tiles = [tiles(1) refineTiles(tiles(4),res)];

%%Let MagTense find the magnetization vector of the cube by iterating to a
%%self-consistent solution. The two empty arguments ( [] ) ensure that
%%default values are used and they are in any case not relevant for this
%%example. The 1e-6 argument is the maximum relative error (change in
%%magnetization from one iteration to the next) allowed before convergence
%%is reached. Finally, 100 reflects the max. no. of allowed iterations

mur = 200:10:200;
figure;
alpha 0.3
hold on
for i=1:length(tiles)
   rectangle('position',[tiles(i).offset(1:2)-tiles(i).abc(1:2)/2 tiles(i).abc(1:2)]); 
end
col = jet(length(mur));
for i=1:1%length(mur)
    tiles = IterateMagnetization( tiles, [], [], 1e-4, 1000 );
%     load('tiles_mur_200.mat')
%     if i==1
%         figure(2);
%         plotTiles(tiles,true);axis equal;
        off=reshape([tiles.offset],[3,length(tiles)]);
        x=off(1,:);y=off(2,:);
        
        M=reshape([tiles.M],[3,length(tiles)]);
        Mnorm = sqrt(sum(M.^2,1));
        M=M./Mnorm;
        
        colIndM = round((Mnorm-min(Mnorm))./(max(Mnorm)-min(Mnorm)) * ( length(tiles) - 1 ) + 1);
        colM = hot( length(colIndM) );
        for j=1:length(tiles)
            rectangle('position',[tiles(j).offset(1:2)-tiles(j).abc(1:2)/2 tiles(j).abc(1:2)],'facecolor',colM(colIndM(j),:)); 
        end
        colormap(hot);h=colorbar;set(h,'yticklabel',linspace(min(Mnorm),max(Mnorm),11));set(get(h,'title'),'string','Magnetization [A/m]');
        
        u=M(1,:);v=M(2,:);
        quiver(x,y,u,v);        
        disp(mur(i));
        drawnow;
%         alpha 0.3
%         figure(1);
%     end
%     plot( reshape( [tiles(2:end).M], [3,length(tiles(2:end))] ), 'color',col(i,:));
    mu = num2cell( [1 ones(1,length(tiles)-1)*mur(i)] );
    [tiles.mu_r_ea] = mu{:};
    [tiles.mu_r_oa] = mu{:};
end
% figure;
% plotTiles(tiles,true);axis equal;
% alpha 0.3
% disp('hej');
%%Now find the field in a set of points
%define a range of points spanning the xy plane at z=0
x = linspace( 0.5*d,2*d, 101);
y = linspace( -0.25*d,0.25*d, 101);
z = 0.0;

%use meshgrid to fill out the span
[X,Y,Z] = ndgrid(x,y,z);

%get the field
[H,Hnorm] = getHMagTense( tiles, X, Y, Z );
H = H .* mu0;
Hnorm = Hnorm .* mu0;
%figure;
%quiver(X,Y,H(:,:,1),H(:,:,2));
%Plot the solution in a contour map
surf_and_con( X,Y,Hnorm);
% % xlabel('x');
% % ylabel('y');
% 
% %%plot the magnetic field (H). Notice how the field opposes the
% %%magnetization inside the cube.
% % figure(figure2)
% % quiver(fig2,X,Y,H(:,:,1),H(:,:,2),2,'linewidth',2);
% % contour(fig2,X,Y,Hnorm,'linewidth',2);
% % plotTiles(tiles,true);
% % alpha 0.3;
% % axis equal;
% % xlabel('x');
% % ylabel('y');
% 
% %--- Compare with Comsol
% H_com = load('..\..\..\documentation\examples_FEM_validation\Validation_MagTense_Example003_PM_Fe_cubes\Validation_MagTense_Example003_PM_Fe_cubes_Hy_x.txt');
% plot(fig4,H_com(:,1),H_com(:,2),'o')
% H_mag = interpn(X,Y,1./mu0*H(:,:,2),H_com(:,1),0);
% plot(fig4,H_com(:,1),H_mag,'.')
% xlabel(fig4,'x');
% ylabel(fig4,'H_y [T]');
% 
% B_com = load('..\..\..\documentation\examples_FEM_validation\Validation_MagTense_Example003_PM_Fe_cubes\Validation_MagTense_Example003_PM_Fe_cubes_By_x.txt');
% plot(fig5,B_com(:,1),B_com(:,2),'o')
% B_mag = interpn(X,Y,H(:,:,2),B_com(:,1),0);
% indx = find(B_com(:,1) >= -d/2 & B_com(:,1) <= d/2);
% B_mag(indx) = B_mag(indx)+1.2;
% indx = find(B_com(:,1) >= -d/2 & B_com(:,1) <= d/2);
% B_mag(indx) = B_mag(indx)*tiles(2).mu_r_ea;
% plot(fig5,B_com(:,1),B_mag,'.')
% xlabel(fig5,'x');
% ylabel(fig5,'B_y [T]');
% 
% 
% % %define a range of points spanning the xy plane at z=0
% % x = linspace( -0.05,0.05, 20);
% % y = 0.025 + 5e-3; %height above the cube
% % z = linspace( -0.06,0.06, 21);
% % 
% % %use meshgrid to fill out the span
% % [X,Y,Z] = meshgrid(x,y,z);
% % 
% % 
% % %get the field
% % [H,Hnorm] = getHMagTense( tiles, X, Y, Z );
% % H = H .* mu0;
% % Hnorm = Hnorm .* mu0;
% % 
% % %Plot the solution in a contour map
% % surf_and_con( squeeze(X),squeeze(Z),Hnorm);
% % xlabel('x');
% % ylabel('z');
% % 
% % 
% % getFigure();
% % quiver(squeeze(X),squeeze(Z),H(:,:,1),H(:,:,3),2,'linewidth',2);
% % contour(squeeze(X),squeeze(Z),Hnorm,'linewidth',2);
% % 
% % alpha 0.3;
% % axis equal;
% % xlabel('x');
% % ylabel('z');
% % 
% % getFigure();contour(squeeze(X),squeeze(Z),H(:,:,1),'linewidth',2,'showtext','on');
% % xlabel('x');
% % ylabel('z');
end
