
function [] = plotFieldArrows( tiles )
t = hgtransform;

    for i=1:length(tiles)

        switch( tiles(i).tileType )
            case {getMagTileType( 'cylinder' ), getMagTileType( 'circpiece' ) }
                p0 = [tiles(i).r0 * cos( tiles(i).theta0 ), tiles(i).r0 * sin( tiles(i).theta0 ), tiles(i).z0+tiles(i).dz/2];
                dl = tiles(i).dr/4;
                
            case getMagTileType( 'prism' )
                p0 = [0,0,0];
                dl = mean(tiles(i).abc);
            
        end
        p0 = p0 + tiles(i).offset;
        
        un = tiles(i).M / sqrt(sum(tiles(i).M.^2));        
        p1 = p0 - dl * un;
        p2 = p0 + dl * un;
        mArrow3(p1,p2,tiles(i).rotAngles,1,'color','k','parent',t);
        
        %plot the easy axis as well
%         if tiles(i).tileType == getMagTileType( 'cylinder' ) && tiles(i).magnetType == getMagnetType('hard')
%            un = tiles(i).u_ea;
%            p1 = p0 - dl * un;
%            p2 = p0 + dl * un;
%            mArrow3(p1,p2,tiles(i).rotAngles,1,'color','b');
%         end
        if isfield(tiles(i),'graphRotxAng')
            Rx=makehgtform('xrotate',tiles(i).graphRotxAng);
            t.Matrix=Rx;    
        end
    end
    
end