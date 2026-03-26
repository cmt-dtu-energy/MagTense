
%defines a class for doing static evaluations of tiles
classdef MagTenseTilesPlot
    properties
       
    end
    
    methods (Static)
       
        function [] = plotTiles( tiles, arrows )
            if ~exist( 'arrows', 'var' )
                arrows = false;
            end
            
            for i=1:length(tiles)
                switch tiles(i).tileType
                    % case getMagTileType( 'cylinder' )
                    case MagTenseTilesUtil.getMagTileType( 'Cylinder' )
                        %cylinder piece
                        MagTenseTilesPlot.plotCylPiece( tiles(i));
                    case MagTenseTilesUtil.getMagTileType( 'Prism' )
                        %a cube
                        %cube_plot(tiles(i).offset,tiles(i).abc,tiles(i).rotAngles,3,tiles(i).color, tiles(i).graphRotxAng);
                        MagTenseTilesPlot.cube_plot( tiles(i), [1,2,3] );
                    case MagTenseTilesUtil.getMagTileType( 'Circpiece' )
                        MagTenseTilesPlot.plotCircPiece( tiles(i) );
                    case MagTenseTilesUtil.getMagTileType( 'Circpieceinv' )
                        MagTenseTilesPlot.plotCircPieceInv( tiles(i) );
                    case MagTenseTilesUtil.getMagTileType('Tetrahedron')
                        MagTenseTilesPlot.plotTetrahedron(tiles(i));
                end
            end
            
            if arrows
                plotFieldArrows( tiles );
            end
        end
        
       
        %%
        function [] = plotFieldArrows( tiles )
        t = hgtransform;
        
            for i=1:length(tiles)
        
                switch( tiles(i).tileType )
                    case {MagTenseTilesUtil.getMagTileType( 'Cylinder' ), MagTenseTilesUtil.getMagTileType( 'Circpiece' ) }
                        p0 = [tiles(i).r0 * cos( tiles(i).theta0 ), tiles(i).r0 * sin( tiles(i).theta0 ), tiles(i).z0+tiles(i).dz/2];
                        dl = tiles(i).dr/4;
                        
                    case (MagTenseTilesUtil.getMagTileType( 'Prism' ) )
                        p0 = [0,0,0];
                        dl = mean(tiles(i).abc);
                    
                end
                p0 = p0 + tiles(i).offset;
                
                un = tiles(i).M / sqrt(sum(tiles(i).M.^2));        
                p1 = p0 - dl * un;
                p2 = p0 + dl * un;
                MagTenseTilesPlot.mArrow3(p1,p2,tiles(i).rotAngles,1,'color','k','parent',t);
                
                %plot the easy axis as well
                if ((MagTenseTilesUtil.getMagTileType( 'Cylinder' )) && (MagTenseTilesUtil.getMagnetType( 'Hard' )))
                   un = tiles(i).u_ea;
                   p1 = p0 - dl * un;
                   p2 = p0 + dl * un;
        %            mArrow3(p1,p2,tiles(i).rotAngles,1,'color','b');
                end
                if isfield(tiles(i),'graphRotxAng')
                    Rx=makehgtform('xrotate',tiles(i).graphRotxAng);
                    t.Matrix=Rx;    
                end
            end
            
        end
       

        %%
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

        %%
        function [P] = plotCircPieceInv( cP )
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
        
            x3 = x1;
            y3 = y2;
                
            %%plot the top and bottom
            x = R.*cos(t) + h;
            y = R.*sin(t) + k;
                
            x = [x, linspace(x(end),x(1),n)];
            y = [y, zeros(1,n) + y(end)];
            
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
            x = zeros(1,n) + x1;
            y = linspace(y1,y2,n);
            
            x = [x fliplr(x) x(1)];
            y = [y fliplr(y) y(1)];
            zup = ones(1,length(x)).*(cP.z0 + cP.dz/2);
            zdn = ones(1,length(x)).*(cP.z0 - cP.dz/2);
            Pleft = patch(x,y,z,col,'parent',thg);
            
            %%plot the y-constant piece
            x = linspace(x3,x2,n);
            y = zeros(1,n) + y3;
            
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

        %%
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
       
        %%
        function plotTetrahedron(tile)
        
            inds = [[1,2,3];[1,3,4];[1,2,4];[2,3,4]];
            
            patch('faces',inds,'vertices',tile.vertices','facecolor',tile.color);
        
        end
       
        %%
        %function cube_plot(origin,dims,rot,order,color,grphRotAng)
        function cube_plot( tile, order )
            origin = tile.offset;
            dims = tile.abc;
            rot = tile.rotAngles;
            color = tile.color;
            if isfield(tile,'graphRotxAng')
                grphRotAng = tile.graphRotxAng;
            else
                grphRotAng = 0;
            end
            %addpath('../utils/');
            % CUBE_PLOT plots a cube with dimension of X, Y, Z.
            %
            % INPUTS:
            % origin = set origin point for the cube in the form of [x,y,z].
            % X      = cube length along x direction.
            % Y      = cube length along y direction.
            % Z      = cube length along z direction.
            % color  = STRING, the color patched for the cube.
            %         List of colors
            %         b blue
            %         g green
            %         r red
            %         c cyan
            %         m magenta
            %         y yellow
            %         k black
            %         w white
            % OUPUTS:
            % Plot a figure in the form of cubics.
            %
            % EXAMPLES
            % cube_plot(2,3,4,'red')
            %
            % ------------------------------Code Starts Here------------------------------ %
            % Define the vertexes of the unit cubic
            X = dims(1);
            Y = dims(2);
            Z = dims(3);
            
            ver = [1 1 0;
                0 1 0;
                0 1 1;
                1 1 1;
                0 0 1;
                1 0 1;
                1 0 0;
                0 0 0];
            %move the vertices so the cube is centered on origo
            for i=1:3
                ver(:,i) = ver(:,i) - 0.5;
            end
            %  Define the faces of the unit cubic
            fac = [1 2 3 4;
                4 3 5 6;
                6 7 8 5;
                1 2 8 7;
                6 7 1 4;
                2 3 5 8];
            %rotate the cube according to the angles specified in rot and the order of
            %rotation specified in order
            ver = MagTenseTilesPlot.rotateVertices( rot, ver, order );
            
            %move the vertices back to have lower left corner in (0,0,0)
            for i=1:3
                ver(:,i) = ver(:,i) + 0.5;
            end
            
            cube = [ver(:,1)*X-X/2+origin(1),ver(:,2)*Y-Y/2+origin(2),ver(:,3)*Z-Z/2+origin(3)];
            t = hgtransform;
            
            patch('Faces',fac,'Vertices',cube,'FaceColor',color,'parent',t);
            
            Rx=makehgtform('xrotate',grphRotAng);
            t.Matrix=Rx;
        end

        %%
        %%rot is an array with n elements, each representing an angle in radians
        %%ver is an array with vertices assumed to be [m,3] in dimensions
        %%order is an array with n elements each representing a rotation about an
        %%axis, i.e. the elements of order should be 1, 2 or 3. All other values
        %%are not allowed
        function [ver] = rotateVertices( rot, ver, order )
        
        
        RotOut = MagTenseTilesUtil.getRotationMatrices( rot );
        
        for i=1:length(order)
            %rotation about the axis defined in order    
            ver = ( RotOut{order(i)} * ver')';
        end
        
        
        end

        %%
        %%ang is an array with three elements (rotation angles about x-, y- and
        %%z-axes, respectively). Angles are assumed in radians
        %%Rerturns a cell array with three elements: the rotation matrices for
        %%rotation about the x-, y- and z-axes, respectively.
        function [RotOut] = getRotationMatrices( ang )
        
            RotOut = cell(3,1);
            
            %rotation about x-axis
            RotOut{1} = [ 1, 0, 0; 0, cos(ang(1)), -sin(ang(1)); 0, sin(ang(1)), cos(ang(1))];
            
            %rotation about y-axis
            RotOut{2} = [cos(ang(2)), 0, sin(ang(2)); 0, 1, 0; -sin(ang(2)), 0, cos(ang(2))];
            
            %rotation about z-axis
            RotOut{3} = [ cos(ang(3)), -sin(ang(3)), 0; sin(ang(3)), cos(ang(3)), 0; 0, 0, 1];
        
        end

        %%
        function h = mArrow3(p1,p2,rot,order,varargin)
            %mArrow3 - plot a 3D arrow as patch object (cylinder+cone)
            %
            % syntax:   h = mArrow3(p1,p2)
            %           h = mArrow3(p1,p2,'propertyName',propertyValue,...)
            %
            % with:     p1:         starting point
            %           p2:         end point
            %           properties: 'color':      color according to MATLAB specification
            %                                     (see MATLAB help item 'ColorSpec')
            %                       'stemWidth':  width of the line
            %                       'tipWidth':   width of the cone                       
            %
            %           Additionally, you can specify any patch object properties. (For
            %           example, you can make the arrow semitransparent by using
            %           'facealpha'.)
            %                       
            % example1: h = mArrow3([0 0 0],[1 1 1])
            %           (Draws an arrow from [0 0 0] to [1 1 1] with default properties.)
            %
            % example2: h = mArrow3([0 0 0],[1 1 1],'color','red','stemWidth',0.02,'facealpha',0.5)
            %           (Draws a red semitransparent arrow with a stem width of 0.02 units.)
            %
            % hint:     use light to achieve 3D impression
            %
            
            propertyNames = {'edgeColor'};
            propertyValues = {'none'};    
            
            %% evaluate property specifications
            for argno = 1:2:nargin-4
                switch varargin{argno}
                    case 'color'
                        propertyNames = {propertyNames{:},'facecolor'};
                        propertyValues = {propertyValues{:},varargin{argno+1}};
                    case 'stemWidth'
                        if isreal(varargin{argno+1})
                            stemWidth = varargin{argno+1};
                        else
                            warning('mArrow3:stemWidth','stemWidth must be a real number');
                        end
                    case 'tipWidth'
                        if isreal(varargin{argno+1})
                            tipWidth = varargin{argno+1};
                        else
                            warning('mArrow3:tipWidth','tipWidth must be a real number');
                        end
                    otherwise
                        propertyNames = {propertyNames{:},varargin{argno}};
                        propertyValues = {propertyValues{:},varargin{argno+1}};
                end
            end            
            
            %% default parameters
            if ~exist('stemWidth','var')
                ax = axis;
                if numel(ax)==4
                    stemWidth = norm(ax([2 4])-ax([1 3]))/300;
                elseif numel(ax)==6
                    stemWidth = norm(ax([2 4 6])-ax([1 3 5]))/300;
                end
            end
            if ~exist('tipWidth','var')
                tipWidth = 3*stemWidth;
            end
            tipAngle = 22.5/180*pi;
            tipLength = tipWidth/tan(tipAngle/2);
            ppsc = 50;  % (points per small circle)
            ppbc = 250; % (points per big circle)
            
            %% ensure column vectors
            p1 = p1(:);
            p2 = p2(:);
            
            %% basic lengths and vectors
            x = (p2-p1)/norm(p2-p1); % (unit vector in arrow direction)
            y = cross(x,[0;0;1]);    % (y and z are unit vectors orthogonal to arrow)
            if norm(y)<0.1
                y = cross(x,[0;1;0]);
            end
            y = y/norm(y);
            z = cross(x,y);
            z = z/norm(z);
            
            %% basic angles
            theta = 0:2*pi/ppsc:2*pi; % (list of angles from 0 to 2*pi for small circle)
            sintheta = sin(theta);
            costheta = cos(theta);
            upsilon = 0:2*pi/ppbc:2*pi; % (list of angles from 0 to 2*pi for big circle)
            sinupsilon = sin(upsilon);
            cosupsilon = cos(upsilon);
            
            %% initialize face matrix
            f = NaN([ppsc+ppbc+2 ppbc+1]);
            
            %% normal arrow
            if norm(p2-p1)>tipLength
                % vertices of the first stem circle
                for idx = 1:ppsc+1
                    v(idx,:) = p1 + stemWidth*(sintheta(idx)*y + costheta(idx)*z);
                end
                % vertices of the second stem circle
                p3 = p2-tipLength*x;
                for idx = 1:ppsc+1
                    v(ppsc+1+idx,:) = p3 + stemWidth*(sintheta(idx)*y + costheta(idx)*z);
                end
                % vertices of the tip circle
                for idx = 1:ppbc+1
                    v(2*ppsc+2+idx,:) = p3 + tipWidth*(sinupsilon(idx)*y + cosupsilon(idx)*z);
                end
                % vertex of the tiptip
                v(2*ppsc+ppbc+4,:) = p2;
            
                % face of the stem circle
                f(1,1:ppsc+1) = 1:ppsc+1;
                % faces of the stem cylinder
                for idx = 1:ppsc
                    f(1+idx,1:4) = [idx idx+1 ppsc+1+idx+1 ppsc+1+idx];
                end
                % face of the tip circle
                f(ppsc+2,:) = 2*ppsc+3:(2*ppsc+3)+ppbc;
                % faces of the tip cone
                for idx = 1:ppbc
                    f(ppsc+2+idx,1:3) = [2*ppsc+2+idx 2*ppsc+2+idx+1 2*ppsc+ppbc+4];
                end
            
            %% only cone v
            else
                tipWidth = 2*sin(tipAngle/2)*norm(p2-p1);
                % vertices of the tip circle
                for idx = 1:ppbc+1
                    v(idx,:) = p1 + tipWidth*(sinupsilon(idx)*y + cosupsilon(idx)*z);
                end
                % vertex of the tiptip
                v(ppbc+2,:) = p2;
                % face of the tip circle
                f(1,:) = 1:ppbc+1;
                % faces of the tip cone
                for idx = 1:ppbc
                    f(1+idx,1:3) = [idx idx+1 ppbc+2];
                end
            end
            
            %%Rotate the arrow (change made by Kaspar K. Nielsen, 4. May 2017),
            %%kaki@dtu.dk
            
            %rotate the cube according to the angles specified in rot and the order of
            %rotation specified in order
            %translate the arrow to the origo
            for i=1:length(p1)
                v(:,i) = v(:,i) - p1(i);
            end
            %rotate
            v = MagTenseTilesUtil.rotateVertices( rot, v, order );
            %move it back again
            for i=1:length(p1)
                v(:,i) = v(:,i) + p1(i);
            end
            
            %% draw
            fv.faces = f;
            fv.vertices = v;
            h = patch(fv);
            for propno = 1:numel(propertyNames)
                try
                    set(h,propertyNames{propno},propertyValues{propno});
                catch
                    disp(lasterr)
                end
            end
        end

   end
    
    
end