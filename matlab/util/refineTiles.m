function [tiles_ref] = refineTiles( tiles, res )

    %refine each tile in tiles array given the specified resolution
    cnt = 1;
    
    for i=1:length(tiles)
        switch tiles(i).tileType
            case getMagTileType( 'cylinder' )
            dr = tiles(i).dr / res.nr;
            dtheta = tiles(i).dtheta / res.ntheta;
            dz = tiles(i).dz / res.nz;

            r1 = tiles(i).r0 - tiles(i).dr/2;
            theta1 = tiles(i).theta0 - tiles(i).dtheta/2;
            z1 = tiles(i).z0 - tiles(i).dz/2;
            for j=1:res.nr
               for k=1:res.ntheta
                   for l=1:res.nz

                       tiles_ref(cnt) = copyStruct( tiles(i) );
                       tiles_ref(cnt).dr = dr;
                       tiles_ref(cnt).dtheta = dtheta;
                       tiles_ref(cnt).dz = dz;
                       tiles_ref(cnt).r0 = r1 + dr * ( 0.5 + (j-1) );
                       tiles_ref(cnt).theta0 = theta1 + dtheta * ( 0.5 + (k-1) );
                       tiles_ref(cnt).z0 = z1 + dz * ( 0.5 + (l-1) );

                       cnt = cnt + 1;
                   end
               end
            end
            case getMagTileType( 'prism' )
            dx = tiles(i).abc(1) / res.nx;
            dy = tiles(i).abc(2) / res.ny;
            dz = tiles(i).abc(3) / res.nz;
            
            x1 = tiles(i).offset(1) - tiles(i).abc(1)/2;
            y1 = tiles(i).offset(2) - tiles(i).abc(2)/2;
            z1 = tiles(i).offset(3) - tiles(i).abc(3)/2;
            for j=1:res.nx
                for k=1:res.ny
                    for l=1:res.nz
                       tiles_ref(cnt) = copyStruct( tiles(i) );
                       tiles_ref(cnt).abc = [dx,dy,dz];
                       
                       tiles_ref(cnt).offset(1) = x1 + dx * ( 0.5 + (j-1) );
                       tiles_ref(cnt).offset(2) = y1 + dy * ( 0.5 + (k-1) );
                       tiles_ref(cnt).offset(3) = z1 + dz * ( 0.5 + (l-1) );

                       cnt = cnt + 1;
                    end
                end
            end
            case getMagTileType( 'circpieceinv') 
                %only refine the z-direction
                dz = tiles(i).dz / res.nz;
                z1 = tiles(i).z0 - tiles(i).dz/2;
                for l=1:res.nz
                    tiles_ref(cnt) = copyStruct( tiles(i) );
                    tiles_ref(cnt).dz = dz;
                    tiles_ref(cnt).z0 = z1 + dz * ( 0.5 + (l-1) );
                    cnt = cnt + 1;
                end
            case getMagTileType( 'circpiece') 
                %only refine the z-direction
                dz = tiles(i).dz / res.nz;
                z1 = tiles(i).z0 - tiles(i).dz/2;
                for l=1:res.nz
                    tiles_ref(cnt) = copyStruct( tiles(i) );
                    tiles_ref(cnt).dz = dz;
                    tiles_ref(cnt).z0 = z1 + dz * ( 0.5 + (l-1) );
                    cnt = cnt + 1;
                end
                
            otherwise
                tiles_ref(cnt) = copyStruct( tiles(i) );
                cnt = cnt + 1;
        end
    end
end
