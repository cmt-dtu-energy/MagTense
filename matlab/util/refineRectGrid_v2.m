function [pos_out,dims_out,children_out,nn] = refineRectGrid_v2( pos, dims, pts_central, ndim, L, criteria)
    
    if (~exist('criteria','var'))
        criteria.name = 'voronoi_dist';
    end
    
    if (isempty(criteria))
        criteria.name = 'voronoi_dist';
    end

    textprogressbar('init');
    disp('');
    disp('Making reference tree');
    ref_tree = getRefinedTree( [0,0,0], dims(1,:), ndim, L );
    
    pos_ = cell( length(pos(:,1)), 1 );
    dims_ = cell( length(pos(:,1)), 1 );
    children_ = cell( length(pos(:,1)), 1 );
    nn_ = cell( length(pos(:,1)), 1 );
    
    disp('Refine all cells');
    for i=1:length( pos(:,1) )
        cell_tree = getRefTree( pos(i,:), ref_tree );
        %disp(['Refining cell nr. ' num2str(i) ' out of ' num2str(length(pos(:,1)))]);
        cell_tree = refineCell2( cell_tree, pts_central, ndim, criteria );
        
        
        [pos_{i},dims_{i},children_{i}, nn_{i}] = structArrToDouble( cell_tree, ndim );
        %disp(['Number of refined elements ' num2str(length(pos_{i}(:,1)))]);
        textprogressbar(100*i/length(pos(:,1)));
    end
    disp('Converting to array');
    n = 0;
    for i=1:length(pos_)
       n = n + length(pos_{i}); 
    end
    pos_out = zeros( n, 3 );
    dims_out = zeros( n, 3 );
    children_out = zeros( n, 4 );
    
    cnt = 1;
    for i=1:length(pos_)
       
       pos_out( cnt:cnt+length(pos_{i}(:,1))-1, : ) = pos_{i};
       dims_out( cnt:cnt+length(pos_{i}(:,1))-1, : ) = dims_{i};
       children_out( cnt:cnt+length(pos_{i}(:,1))-1, 1:length(children_{i}(1,:)) ) = children_{i};       
       
       cnt = cnt + length(pos_{i}(:,1));
    end
    
    ind = find(pos_out(:,1)>0);
    pos_out = pos_out(ind,:);
    dims_out = dims_out(ind,:);
    children_out = children_out(ind,:);
    
%    if ndim == 2
%       dims_out(:,3) = dims_out(:,1); 
%    end
    
    nn = getNearestNeighbors( pos_out, pts_central, ndim );
    textprogressbar('Done');

end

function [nn] = getNearestNeighbors( pos, pts_central, ndim )
    
    nn = zeros( length(pos(:,1)), 1 );
    for i=1:length(pos(:,1))
        if ndim==2
            dist_sq = (pos(i,1)-pts_central(:,1)).^2 + (pos(i,2)-pts_central(:,2)).^2;% + (pos(i,3)-pts_central(:,3)).^2;
        else
            dist_sq = (pos(i,1)-pts_central(:,1)).^2 + (pos(i,2)-pts_central(:,2)).^2 + (pos(i,3)-pts_central(:,3)).^2;
        end

        nn(i) = find( dist_sq == min(dist_sq) );
    end
end

function [ref_tree] = getRefTree( pos, ref_tree )
    for i=1:length(ref_tree)
        for j=1:length(ref_tree{i}.pos(:,1))
            ref_tree{i}.pos(j,:) = ref_tree{i}.pos(j,:) + pos;
        end
    end
end

function [ref_tree] = refineCell2( ref_tree, pts_central, ndim, criteria )

    for i=length(ref_tree)-1:-1:1
        
        for j=1:length(ref_tree{i}.pos(:,1))
            %get the children of the current node
            chrd = ref_tree{i}.children(j,:);
            %positions of the current node's children
            pp = ref_tree{i+1}.pos( chrd, : );
            
            if sum( chrd ) == 0
                %if no children, then leave the node alone
                %disp('no children');
            else
                
                hasGrandChildren = false;
                for k=1:length(chrd)
                   if chrd(k)>0
                       g_chrd = ref_tree{i+1}.children(chrd(k),:);
                       if sum(g_chrd)>0
                           hasGrandChildren = true;
                           break;
                       end
                   end
                end
                if hasGrandChildren
                    %if children and grand children, leave the node alone
                    %disp('has grand children');
                else
                    %if the node has children but no grand children, check to see if it should
                    %be refined or not    
                    NN = zeros( length(chrd), 1 );
                    NN2 = zeros( length(chrd), 1 );
                    
                    for k=1:length(chrd)
                        if (strcmp(criteria.name,'voronoi_dist') || strcmp(criteria.name,'both'))
                            if ndim == 2
                                dist_sq = (pp(k,1)-pts_central(:,1)).^2 + (pp(k,2)-pts_central(:,2)).^2;
                            else
                                dist_sq = (pp(k,1)-pts_central(:,1)).^2 + (pp(k,2)-pts_central(:,2)).^2 + (pp(k,3)-pts_central(:,3)).^2;
                            end

                            NN(k) = find( dist_sq == min(dist_sq) );
                            
                            if (k == length(chrd))
                                NN = [0 diff(NN)'];
                            end
                        end
                        if (strcmp(criteria.name,'field_grad') || strcmp(criteria.name,'both'))
                            %--- Find current refinement level
                            ds = unique(abs(diff(pp)));
                            ds = ds(2); %--- The first value is always zero
%                             disp(num2str(k));
                            indx = find(abs(ds-1./criteria.ds) < 1e-8);
                            
                            %--- Find the current position in the provided
                            %--- field grid
                            if ndim == 2
                                indx_xy = find((abs(pp(k,1)-criteria.X{indx}) < 1e-8) & abs(pp(k,2)-criteria.Y{indx}) < 1e-8);
                            else
                                indx_xy = find((abs(pp(k,1)-criteria.X{indx}) < 1e-8) & abs(pp(k,2)-criteria.Y{indx}) < 1e-8 & abs(pp(k,3)-criteria.Z{indx}) < 1e-8);
                            end
                            
                            %--- Get the value of the field at the current position
                            g_norm_pos = criteria.field{indx}(indx_xy);
                                
                            %--- If the field contains NaN then refine
                            if (~isnan(g_norm_pos))
                                NN2(k) = max([0 (g_norm_pos > criteria.g_value)]);
                            else
                                NN2(k) = 1;
                            end
                            
                            %--- Merge the criteria for refinement or use
                            %--- only the field criteria
                            if (k == length(chrd))
                                if (strcmp(criteria.name,'both'))
                                    NN(k) = max([NN(k) NN2(k)]);
                                else
                                    NN = NN2;
                                end
                            end
                            
                        end
                    end
                    if (~any(NN))
                        %no reason to keep refined, delete the children
                        ref_tree{i+1}.pos(chrd,:) = 0;
                        ref_tree{i+1}.dims(chrd,:) = 0;
                        ref_tree{i+1}.children(chrd,:) = 0;
                        ref_tree{i+1}.nn(chrd) = 0;
                        
                        %delete the reference to the children in the
                        %current node
                        ref_tree{i}.children(j,:) = 0;
                    else
                        ref_tree{i+1}.nn(chrd) = NN;
                    end
                end
            end            
        end
        
    end

end



function [ pos, dims, children, nn ] = structArrToDouble( str_in, ndim )
    nOut = 0;
    
    for i=1:length(str_in)
        
        nOut = nOut + (2^ndim)^(i-1);        
    end
    
    pos = zeros( nOut, 3 );
    dims = zeros( nOut, 3 );
    children = zeros( nOut, 2^ndim );
    nn = zeros( nOut, 1 );
    
    pos(1,:) = str_in{1}.pos;
    dims(1,:) = str_in{1}.dims;
    children(1,:) = str_in{1}.children;
    nn(1) = str_in{1}.nn;
    
    cnt = 2;
    for i=2:length(str_in);
        L = length(str_in{i}.pos(:,1));
        
        pos( cnt:1:(L+cnt-1),: ) = str_in{i}.pos;
        dims( cnt:1:(L+cnt-1),: ) = str_in{i}.dims;
        children( cnt:1:(L+cnt-1),1:length(str_in{i}.children(1,:)) ) = str_in{i}.children;

        nn( cnt:1:(L+cnt-1), 1:length(str_in{i}.nn(1,:)) ) = str_in{i}.nn;
        
        cnt = cnt + L;
    end
    
  
    
    ind = find( dims(:,1) ~= 0 );
    pos = pos(ind,:);
    dims = dims(ind,:);
    children = children(ind,:);
    nn = nn(ind);
    
    k = sum( children, 2 );
    
    ind = find( k == 0 );
    
    pos = pos(ind,:);
    dims = dims(ind,:);
    children = children(ind,:);
    nn = nn(ind,:);
    
     
end

function [cell_out] = getRefinedTree( pos, dims, ndim, L )
    
    cell_out = cell( L, 1 );
        
    cell_out{1} = struct();
    cell_out{1}.pos = pos;
    cell_out{1}.dims = dims;
    cell_out{1}.nCh = 1;
    cell_out{1}.children = zeros( 1, 2^ndim );
    cell_out{1}.nn = 0;
    
    %loop over branch level
    for i=2:L
        
        nCh = (2^ndim)^(i-1);
        
        cell_out{i}.pos = zeros( nCh, 3 );
        cell_out{i}.dims = zeros( nCh, 3 );
        cell_out{i}.children = zeros( nCh, 2^ndim );
        cell_out{i}.nCh = nCh;
        cell_out{i}.nn = zeros( nCh, 1 );
        
        for j=1:cell_out{i-1}.nCh   
           if ndim==2
               [pos_tmp, dims_tmp] = getRefined2DRectangle( cell_out{i-1}.pos(j,:), cell_out{i-1}.dims(j,:) );
           else
               [pos_tmp, dims_tmp] = getRefined3DRectangle( cell_out{i-1}.pos(j,:), cell_out{i-1}.dims(j,:) );
           end

           indStart = (j-1) * 2^ndim + 1;
           indEnd = j * 2^ndim;
            
           cell_out{i-1}.children(j,:) = indStart:1:indEnd;
                      
           cell_out{i}.pos(indStart:1:indEnd,:) = pos_tmp;
           cell_out{i}.dims(indStart:1:indEnd,:) = dims_tmp;
           
        end

    end
    
    

end

function [pos_out,dims_out] = getRefined3DRectangle( pos, dims )
    [pos1,dims1] = getRefined2DRectangle( pos, dims );
    
    pos_out = zeros(8,3);
    dims_out = zeros(8,3);
    
    pos1(:,3) = pos(3) - dims(3)/4;
    dims1(:,3) = dims(3)/2;
    
    pos_out(1:4,:) = pos1;
    dims_out(1:4,:) = dims1;
    
    pos1(:,3) = pos(3) + dims(3)/4;
    
    pos_out(5:8,:) = pos1;
    dims_out(5:8,:) = dims1;
    
end

function [pos_out,dims_out] = getRefined2DRectangle( pos, dims)
    pos_out = zeros(4,3);
    dims_out = zeros(4,3);
    %%lower left rectangle
    pos_out(1,1) = pos(1) - dims(1)/4;
    pos_out(1,2) = pos(2) - dims(2)/4;
    pos_out(1,3) = pos(3);

    %upper left rectangle
    pos_out(2,1) = pos(1) - dims(1)/4;
    pos_out(2,2) = pos(2) + dims(2)/4;
    pos_out(2,3) = pos(3);

    %lower right rectangle
    pos_out(3,1) = pos(1) + dims(1)/4;
    pos_out(3,2) = pos(2) - dims(2)/4;
    pos_out(3,3) = pos(3);

    %upper right rectangle
    pos_out(4,1) = pos(1) + dims(1)/4;
    pos_out(4,2) = pos(2) + dims(2)/4;
    pos_out(4,3) = pos(3);
    
    dims_out(:,1) = dims(1)/2;
    dims_out(:,2) = dims(2)/2;
    dims_out(:,3) = dims(3);
end