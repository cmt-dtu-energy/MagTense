function [pos_out,dims_out,children_out,nn] = refineRectGrid_v4( pos, dims, pts_central, ndim, L, criteria, Lsize)
% N cells
% pos_out and dims_out are N x 3 array
% nn is a N x 1 array of integers (index of the corresponding grain/Voronoi cell)
% children_out is a N x (2^dim)    i.e. 4 for 2D and 8 for 3D


%     textprogressbar = @(s) disp(s) ;
        textprogressbar = @(s) s ;
    if (~exist('criteria','var'))
        criteria.name = 'voronoi_dist';
    end
    
    if (isempty(criteria))
        criteria.name = 'voronoi_dist';
    end

    textprogressbar('init');
    disp('');
    disp('Making reference tree');

    ref_tree = getRefinedTree( [0,0,0], dims(1,:), ndim, L ); % get the reference tree by refining a reference cell (with position pos = [0,0,0] and dimension dims)
    
    pos_ = cell( length(pos(:,1)), 1 );
    dims_ = cell( length(pos(:,1)), 1 );
    children_ = cell( length(pos(:,1)), 1 );
    nn_ = cell( length(pos(:,1)), 1 );
    
    disp('Refine all cells');
    for i=1:length( pos(:,1) ) % cycles over the cells of the coarse grid
        cell_tree = getRefTree( pos(i,:), ref_tree ); % Shift Reference tree to current cell
        
        % determine which nodes are relevant (set the rest to zero)
        cell_tree = refineCell2( cell_tree, pts_central, ndim, criteria ); 
        
        % extract relevant nodes from tree (i.e. not zero, and with no children) 
        % and assigns them to 1-D array (from first generation to last)
        [pos_{i},dims_{i},children_{i}, nn_{i}] = structArrToDouble( cell_tree, ndim ); 

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

    
    if 1 % what's the purpose of this code ?
%         ind = find(pos_out(:,1)>0); This is how it was
        ind = find(dims_out(:,1)>0);
        pos_out = pos_out(ind,:);
        dims_out = dims_out(ind,:);
        children_out = children_out(ind,:);
    end
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
    for i=1:length(ref_tree) % L
        for j=1:length(ref_tree{i}.pos(:,1))
            ref_tree{i}.pos(j,:) = ref_tree{i}.pos(j,:) + pos; % shifts the tree to "pos" 
        end
    end
end

function [ref_tree] = refineCell2( ref_tree, pts_central, ndim, criteria )
% analyzes the tree "ref_tree" and modifies it 
% by setting to zero the unnecessary nodes
%
% 
% PSEUDOCODE
% cycle (i) over level of the tree going backwards !
%    cycle (j) over all nodes on that level of the tree
%        get children of current node
%        if node has no children
%            then leave the node alone
%        else
%           cycle (k) over children of current node
%               get children of current children (i.e. grandchildren of current node)
%           end cycle (k)
%           if children and grand children
%                   leave the node alone
%           else
%              cycle (k) over children of current node
%                   VoronoiDistanceRefinement
%                   FieldGradientRefinement
%              end cycle (k)
%              if none of the children is special (no reason to keep refined)
%                   delete the children (set pos dims and children to zero)
%              else
%                   leave the node alone (store closest voronoi generator of its children)
%              end
%           end
%       end
%   end cycle(j)
% end cycle(i)

    for i=length(ref_tree)-1:-1:1 % cycle over level of the tree going backwards !

        for j=1:length(ref_tree{i}.pos(:,1)) % cycle over nodes on the same level of the tree
            %get the children of the current node
            chrd = ref_tree{i}.children(j,:);
            %positions of the current node's children
            pp = ref_tree{i+1}.pos( chrd, : );
            
            if sum( chrd ) == 0
                %if no children, then leave the node alone
                %disp('no children');
            else
                
                hasGrandChildren = false;
                for k=1:length(chrd) % cycle over children of current node
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
                            NN = VoronoiDistanceRefinement(pp,pts_central,ndim,chrd,NN,k) ;
                        end
                        if (strcmp(criteria.name,'field_grad') || strcmp(criteria.name,'both'))
                            NN = FieldGradientRefinement(pp,criteria,ndim,chrd,NN,NN2,k) ;     
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
% 1) Calculate the total number of nodes in the tree and create a 1-D array
% 2) Put the base node as first element of the 1-D array
% 3) Take all the other nodes in the tree and assign to the array
% 4) Take all the nodes with 1st dim = 0: remove them from the array
% 5) Take all the nodes with children and remove them from the array


    % 1)
    nOut = 0;

    for i=1:length(str_in)

        nOut = nOut + (2^ndim)^(i-1);
    end

    pos = zeros( nOut, 3 );
    dims = zeros( nOut, 3 );
    children = zeros( nOut, 2^ndim );
    nn = zeros( nOut, 1 );
    % 2)
    pos(1,:) = str_in{1}.pos;
    dims(1,:) = str_in{1}.dims;
    children(1,:) = str_in{1}.children;
    nn(1) = str_in{1}.nn;
    % 3)
    cnt = 2;
    for i=2:length(str_in);
        L = length(str_in{i}.pos(:,1));

        pos( cnt:1:(L+cnt-1),: ) = str_in{i}.pos;
        dims( cnt:1:(L+cnt-1),: ) = str_in{i}.dims;
        children( cnt:1:(L+cnt-1),1:length(str_in{i}.children(1,:)) ) = str_in{i}.children;

        nn( cnt:1:(L+cnt-1), 1:length(str_in{i}.nn(1,:)) ) = str_in{i}.nn;

        cnt = cnt + L;
    end
    % 4)


    ind = find( dims(:,1) ~= 0 );
    pos = pos(ind,:);
    dims = dims(ind,:);
    children = children(ind,:);
    nn = nn(ind);

    % 5)
    k = sum( children, 2 );

    ind = find( k == 0 );

    pos = pos(ind,:);
    dims = dims(ind,:);
    children = children(ind,:);
    nn = nn(ind,:);

     
end

function [cell_out] = getRefinedTree( pos, dims, ndim, L )
    
    cell_out = cell( L, 1 );
        
    % first level
    cell_out{1} = struct();
    cell_out{1}.pos = pos;
    cell_out{1}.dims = dims;
    cell_out{1}.nCh = 1; % number of nodes on the same level of the tree
    cell_out{1}.children = zeros( 1, 2^ndim );
    cell_out{1}.nn = 0;
    
    %loop over higher levels of the tree
    for i=2:L
        
        nCh = (2^ndim)^(i-1);  % number of nodes on the same level of the tree
        
        cell_out{i}.pos = zeros( nCh, 3 );
        cell_out{i}.dims = zeros( nCh, 3 );
        cell_out{i}.children = zeros( nCh, 2^ndim );
        cell_out{i}.nCh = nCh;
        cell_out{i}.nn = zeros( nCh, 1 );
        
        for j=1:cell_out{i-1}.nCh   % loop over nodes on the same level of the tree
            
            % compute position and dimensions of current node from those of its parent
           if ndim==2
               [pos_tmp, dims_tmp] = getRefined2DRectangle( cell_out{i-1}.pos(j,:), cell_out{i-1}.dims(j,:) );
           else
               [pos_tmp, dims_tmp] = getRefined3DRectangle( cell_out{i-1}.pos(j,:), cell_out{i-1}.dims(j,:) );
           end
           
           % compute indexes of children of current node

           indStart = (j-1) * 2^ndim + 1;
           indEnd = j * 2^ndim;
            
           cell_out{i-1}.children(j,:) = indStart:1:indEnd;
           
           % assign position and dimensions of current node
           cell_out{i}.pos(indStart:1:indEnd,:) = pos_tmp;
           cell_out{i}.dims(indStart:1:indEnd,:) = dims_tmp;
           
        end

    end
    
    

end

function [pos_out,dims_out] = getRefined3DRectangle( pos, dims )
% INPUT: one original block with 
% position pos (1x3 array)
% dimension dim (1x3 array)
% OUTPUT: 8 children blocks
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
function NN = VoronoiDistanceRefinement(pp,pts_central,ndim,chrd,NN,k)
% pp          : positions of the current node's children 
%                           (the current children corresponds to the index k)
% pts_central : positions of the Voronoi generators
    if ndim == 2
        dist_sq = (pp(k,1)-pts_central(:,1)).^2 + (pp(k,2)-pts_central(:,2)).^2;
    else
        dist_sq = (pp(k,1)-pts_central(:,1)).^2 + (pp(k,2)-pts_central(:,2)).^2 + (pp(k,3)-pts_central(:,3)).^2;
    end

    NN(k) = find( dist_sq == min(dist_sq) ); % find the closest generator

    if (k == length(chrd))
        NN = [0 diff(NN)'];
    end
end
function NN = FieldGradientRefinement(pp,criteria,ndim,chrd,NN,NN2,k)

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
        %                             else
        %                                 NN2(k) = 1;
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