
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