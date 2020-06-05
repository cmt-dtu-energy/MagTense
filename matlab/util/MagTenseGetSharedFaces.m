
%By Kaspar K. Nielsen, kasparkn@gmail.com
%First version 12 May 2020
%find any shared faces between the two tiles and returns the relevant info
%in a struct containing the area of the shared face and the distances from
%the center of each tile to the center of the shared face
function face = MagTenseGetSharedFaces( tileA, tileB )
face = [];
%check the types of the two tiles and call the appropriate specific
%implementation for that tile pair

%many pairs need to be implemented here eventually!
if tileA.tileType==getMagTileType('prism')
   %all cases where the first tile is of type prism
   switch tileB.tileType
       case getMagTileType('prism')
           face = getSharedFacePrismPrism( tileA, tileB );
       case getMagTileType('cylinder')
       case getMagTileType('circpiece')
       case getMagTileType('tetrahedron')
           face = getSharedFacePrismTetrahedron( tileA, tileB );
   end
elseif tileA.tileType==getMagTileType('tetrahedron')
    switch tileB.tileType
        case getMagTileType('prism')
            face = getSharedFacePrismTetrahedron( tileB, tileA );
        case getMagTileType('tetrahedron')
            face = getSharedFaceTetraTetra( tileA, tileB );
            
    end
end

end

%Implements the face-finder for the case where tileA is a prism and tileB
%is a tetrahedron
function face = getSharedFacePrismTetrahedron( tileA, tileB )

face = [];

%wiggle room for the comparison
delta = 1e-10;

%permutations for getting the triangular faces of the tetrahedron
perm = [[1,2,3];[1,2,4];[1,3,4];[2,3,4]];

%all eight vertices for the prism
pv = zeros( 8, 3);

%pos x
pv(1:4,1) = tileA.offset(1) + tileA.abc(1)/2;
%neg x
pv(5:8,1) = tileA.offset(1) - tileA.abc(1)/2;

%pos y
pv([1,4,5,8],2) = tileA.offset(2) + tileA.abc(2)/2;
%neg y
pv([2,3,6,7],2) = tileA.offset(2) - tileA.abc(2)/2;

%pos z
pv([1,2,5,6],3) = tileA.offset(3) + tileA.abc(3)/2;
%neg z
pv([3,4,7,8],3) = tileA.offset(3) - tileA.abc(3)/2;

%indices of the vertices for each face
faceInds = [ [1,2,3,4];[5,6,7,8];[1,4,5,8];[2,3,6,7];[1,2,5,6];[3,4,7,8]; ];

inds = [1,1,2,2,3,3];

p_inds = [[2,3];[1,3];[1,2]];

%go through all six faces of the prism (tileA) and all four faces of the
%tetrahedron (tileB)
for i=1:6
    %current face vertices
    pvc = pv( faceInds(i,:), : );
    
    %index defining the coordinate in the current plane
    ind = inds(i);
    
    %get the value of this coordinate
    p = pvc(1,ind);
    
    for j=1:4
        %get the vertices for this particular triangle
       tv = tileB.vertices(:,perm(j,:));
       %get the coordinate comparable to the one for the prism
       q = tv(ind,:);
       if max(abs(p-q)) <= delta
           
           %find the clipping polygon
           cp = sutherlandHodgman( pvc(:,p_inds(ind,:)), tv(p_inds(ind,:),:)' );
           if ~isempty(cp)
               %get the area of the clipping polygon
               pa = polyarea(cp(:,1),cp(:,2));
               if pa>1e-10

                   %center of polygon
                   pc = zeros(3,1);
                   pc(ind) = p;
                   pc(p_inds(ind,1)) = mean(cp(:,1));
                   pc(p_inds(ind,2)) = mean(cp(:,2));

                   %center of tetrahedron
                   tc = mean( tileB.vertices, 2 );

                   face = struct('A',pa,'lA',sqrt(sum((pc'-tileA.offset).^2)),'lB',sqrt(sum((pc-tc).^2)), 'ID_A',i, 'ID_B',j);
               end
               break;
           end
       end
    end
end
end

%implements the face-finder for the case where tileA and tileB both are
%tetrahedra
function face = getSharedFaceTetraTetra( tileA, tileB )
face = [];

%error allowed for when checking whether the faces are in the same plane
delta = 1e-10;

%indices for the vertices in each face
inds = [[1,2,3];[1,3,4];[1,2,4];[2,3,4]];

%loop over each face in tile A
for i=1:4
    
    %find the orthonormal basis this face
    vA = tileA.vertices(:,inds(i,:));
    
    %find the unit vector denoted (1,0,0) in the face's local system
    un_A_1 = (vA(:,2)-vA(:,1))./sqrt(sum((vA(:,2)-vA(:,1)).^2));
    %and the unit vector denoted (0,0,1) in the local system    
    un_A_3 = cross( vA(:,2)-vA(:,1), vA(:,3)-vA(:,1) );
    un_A_3 = un_A_3 / sqrt(sum(un_A_3.^2));
    %and the unit vector denoted (0,1,0) in the local system
    un_A_2 = cross( un_A_3, un_A_1 );
    
    %find the vertices of the face from tile A in the local system
    %basis-change matrix from global to local is the unit vectors of the
    %local system written with respect to the global as column vectors and
    %then that matrix transposed (P = (un1,un2,un3), P^-1 = transp(P) since
    %P is orthonormal)
    Pinv = [un_A_1';un_A_2';un_A_3'];
    
    vA_loc = Pinv * vA;
    
    %and over each face in tile B
    for j=1:4
        
        %vertices of the current face
        vB = tileB.vertices(:,inds(j,:));
        
        %change coordinates of the vertices to the local system of the face
        %from til A
        vB_loc_A = Pinv * vB;
        
        %if the z-coordinates of the two faces' vertices are identical (in
        %the local system) then the faces must be coplanar and we should
        %try to find their intersection
        
        if sum( abs(vA_loc(3,:)-vB_loc_A(3,:)) ) < delta
           %find intersecting polygon
           cp = sutherlandHodgman( vA_loc(1:2,:)', vB_loc_A(1:2,:)' );
           if ~isempty(cp)
               %get the area of the clipping polygon
               pa = polyarea(cp(:,1),cp(:,2));
               if pa>1e-10

                   %center of polygon
                   pc = zeros(3,1);
                   pc(3) = vA_loc(3,1);
                   pc(1) = mean(cp(:,1));
                   pc(2) = mean(cp(:,2));

                   %center of tetrahedron A
                   tcA = mean( Pinv * tileA.vertices, 2 );
                   %center of tetrahedron B
                   tcB = mean( Pinv * tileB.vertices, 2 );

                   face = struct('A',pa,'lA',sqrt(sum((pc-tcA).^2)),'lB',sqrt(sum((pc-tcB).^2)), 'ID_A',i, 'ID_B',j);
                   break;
               end
           end
        end
    end
end
end

%Implements the face-finder for the case of two prisms
function face = getSharedFacePrismPrism( tileA, tileB )
%no rotation is supported yet!
%Go through each possible pair of mutual faces 
%(xA+ with xB-, xA- with xB+ and likewise for the y- and z-faces)
face = [];
pA = tileA.offset;
dA = tileA.abc;

pB = tileB.offset;
dB = tileB.abc;
%wiggle room for the comparison
delta = 0.001 * min([dA,dB]);

%indices for the dimensions that are not i
d=[[2,3];[1,3];[1,2]];
ind = 0;

%center of intersecting rectangle
c = zeros(3,1);
%Face ID's for the shared face for rectangle A and B, respectively.
% 1 = x+, 2 = x-, 3 = y+, 4 = y-, 5 = z+, 6 = z-
ID_A = 0;
ID_B = 0;
for i=1:3
    %check A+ with B-
    if pA(i)+dA(i)/2+delta > pB(i)-dB(i)/2 && pA(i)+dA(i)/2-delta < pB(i)-dB(i)/2
        ind = i;
        c(ind) = pA(i)+dA(i)/2;
        ID_A = 2*(i-1)+1;
        ID_B = 2*i;
        break;
    elseif pA(i)-dA(i)/2-delta < pB(i)+dB(i)/2 && pA(i)-dA(i)/2+delta > pB(i)+dB(i)/2
        %check A- with B+
        ind = i;
        c(ind) = pA(i)-dA(i)/2;
        ID_A = 2*i;
        ID_B = 2*(i-1)+1;
        break;
    end
end

if ind>0
   %make rectangles for the two prism faces in this plane
    r1 = [pA(d(ind,1))-dA(d(ind,1))/2,pA(d(ind,2))-dA(d(ind,2))/2,dA(d(ind,1)),dA(d(ind,2))];
    r2 = [pB(d(ind,1))-dB(d(ind,1))/2,pB(d(ind,2))-dB(d(ind,2))/2,dB(d(ind,1)),dB(d(ind,2))]; 
    %intersect the rects and get the surface area
    A = rectint(r1,r2);
    if A>1e-15 %has to be a finite number otherwise we can get very small artifical surface areas due to round off errors
        %then they intersect. Now, let us find the distance from the center
        %of prism A to the center of the common face and likewise from the center of
        %prism B to center of common face
        %lower left corner of intersecting rect
        clow = [max(r1(1),r2(1)),max(r1(2),r2(2))];
        %upper right corner
        cup = [min(r1(1)+r1(3),r2(1)+r2(3)),min(r1(2)+r1(4),r2(2)+r2(4))];
        
        %fill out the last two numbers for the center of the intersecting
        %prism
        c(d(ind,:)) = 0.5 * ( clow + cup );
        %return a struct with the surface area and distances from the
        %respective centers of the prisms to the center of the shared
        %rectangle.
        face = struct('A',A,'lA',sqrt(sum((c'-pA).^2)),'lB',sqrt(sum((c'-pB).^2)), 'ID_A',ID_A, 'ID_B',ID_B);
    end
end


end