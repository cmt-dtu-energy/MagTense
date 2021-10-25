
%%Two PM cubes and a single Fe cube
function [] = MagTense_Example003_PM_Fe_cubes()
%close all
%figure1= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig1 = axes('Parent',figure1,'Layer','top','FontSize',16); hold on; grid on; box on       
%figure2= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig2 = axes('Parent',figure2,'Layer','top','FontSize',16); hold on; grid on; box on
%figure4= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig4 = axes('Parent',figure4,'Layer','top','FontSize',16); hold on; grid on; box on
%figure5= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig5 = axes('Parent',figure5,'Layer','top','FontSize',16); hold on; grid on; box on

fe_type = 2;
mur = 100;

%make sure to source the right path for the generic Matlab routines
addpath(genpath('../../../util/'));
addpath('../../../MEX_files/');
%define the vacuum permeability
mu0 = 4*pi*1e-7;
Ms = 20.1;

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
tiles(3).u_oa1 = [0,0,1];

if (fe_type == 2)
    tiles(2).magnetType = getMagnetType('soft_const_mur');
else
    tiles(2).magnetType = getMagnetType('soft');
end
tiles(2).offset(2) = 2*d;
tiles(2).Mrem = 0;
tiles(2).M = [0,0,0];
tiles(2).mu_r_ea = 100;
tiles(2).mu_r_oa = 100;
tiles(2).color = [0,0,1];
tiles(2).u_ea = [1,0,0];
tiles(2).u_oa1 = [0,1,0];
tiles(2).u_oa2 = [0,0,1];


tiles(4) = tiles(2);
tiles(4).offset(1) = 2*d;
tiles(4).offset(2) = 0;

res = struct('nx',10,'ny',10,'nz',1);

tile_ex = refineTiles(tiles(4),res);
tiles = [tiles(1) tiles(3) refineTiles([tiles(2) tiles(4)],res)];

%%Let MagTense find the magnetization vector of the cube by iterating to a
%%self-consistent solution. The two empty arguments ( [] ) ensure that
%%default values are used and they are in any case not relevant for this
%%example. The 1e-6 argument is the maximum relative error (change in
%%magnetization from one iteration to the next) allowed before convergence
%%is reached. Finally, 100 reflects the max. no. of allowed iterations

figure;
alpha 0.3
hold on
for i=1:length(tiles)
   rectangle('position',[tiles(i).offset(1:2)-tiles(i).abc(1:2)/2 tiles(i).abc(1:2)]); 
end
col = jet(length(mur));
for i=1:1%length(mur)
    
    if (fe_type == 2)
        tiles = IterateMagnetization( tiles, [], [], 1e-6, 1000 );
    else
        stFcn=MakeMH_Fe_const_mur( mur, Ms );
        tiles = IterateMagnetization( tiles, stFcn, 300, 1e-6, 1000 );
    end
    
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
%         disp(mur(i));
        drawnow;
end

%%Now find the field in a set of points define a range of points spanning the xy plane at z=0
off_ex = reshape([tile_ex.offset],[3,length(tile_ex)]);
x_e = unique(off_ex(1,:));
y_e = unique(off_ex(2,:));
z_e = 0.0;
[x_e,y_e,z_e] = ndgrid(x_e,y_e,0);
x_e = x_e(:)';
y_e = y_e(:)';
z_e = z_e(:)';

x2 = linspace( 1.25*d,2.75*d, 101);
y2 = linspace( -0.75*d,0.75*d, 101);
[x2,y2,z2] = ndgrid(x2,y2,0);

clear x_t y_t z_t 
k = 1;
for i = 1:length(x2)
    for j = 1:length(y2)
        if ~((x2(i,j) >= 1.5*d) && (x2(i,j) <= 2.5*d) && (y2(i,j) >= -0.5*d) && (y2(i,j) <= 0.5*d))
            x_t(k) = x2(i,j); 
            y_t(k) = y2(i,j);
            z_t(k) = 0;
            k = k+1;
        end
    end
end

%get the field
H = getHFromTiles_mex( tiles, [x_t y_t z_t], int32( length(tiles) ), int32( length(x_t) ) );
H2 = getHFromTiles_mex( tiles, [x_e y_e z_e], int32( length(tiles) ), int32( length(x_e) ) );
H = H .* mu0;
H2 = H2 .* mu0;

%Combine the two
H = [H; H2];
x = [x_t x_e];
y = [y_t y_e];
z = [z_t z_e];

%Plot the solution in a contour map
figure; quiver(x',y',H(:,1)./sqrt(sum(H(:,:).^2,2)),H(:,2)./sqrt(sum(H(:,:).^2,2)));

Hx_c = load(['Hx_' num2str(mur) '.txt']);
Hy_c = load(['Hy_' num2str(mur) '.txt']);

Hx_c(:,3) = mu0*Hx_c(:,3);
Hy_c(:,3) = mu0*Hy_c(:,3);

Fx_c = scatteredInterpolant(Hx_c(:,1),Hx_c(:,2),Hx_c(:,3));
Hx_c = Fx_c(x,y);

Fy_c = scatteredInterpolant(Hy_c(:,1),Hy_c(:,2),Hy_c(:,3));
Hy_c = Fy_c(x,y);

figure; quiver(x,y,Hx_c./sqrt(Hx_c.^2+Hy_c.^2),Hy_c./sqrt(Hx_c.^2+Hy_c.^2));

figure; plot3(x,y,H(:,1),'.'); hold all; plot3(x,y,Hx_c,'o');
figure; plot3(x,y,H(:,2),'.'); hold all; plot3(x,y,Hy_c,'o');
% figure; plot3(x,y,(H(:,1)-Hx_c')./Hx_c','.');

x2 = linspace( 1.30*d,2.70*d, 101);
y2 = linspace( -0.70*d,0.70*d, 101);
[x2,y2,z2] = meshgrid(x2,y2,0);

vqx = griddata(x,y,(H(:,1)-Hx_c')./Hx_c',x2,y2);
vqy = griddata(x,y,(H(:,2)-Hy_c')./Hy_c',x2,y2);

figure; surf(x2,y2,vqx,'linestyle','none'); shading interp;
figure; surf(x2,y2,vqy,'linestyle','none'); shading interp;
end
