%The inputs are a table of x-y pairs for the verticies of the subject
%polygon and boundary polygon. (x values in column 1 and y values in column
%2) The output is a table of x-y pairs for the clipped version of the 
%subject polygon.
 %Modified by Kaspar K. Nielsen (kasparkn@gmail.com) to pre-sort the
 %vertices of the clipping polygon such that they are guaranteed to be
 %ordered counter-clockwise
function clippedPolygon = sutherlandHodgman(subjectPolygon,clipPolygon)
 
    %sorting function for ensuring the vertices of the clipping polygon are
    %ordered counter-clockwise
    function clipPolygon = orderVerticesCCW( polygon )
       %find the center of the polygon and use this as the Origin.
       %Then find the polar angle to each vertex and order these angles
       %ascending (CCW is in the positive direction of the increasing
       %angle)
       %then insert the vertices in this order in the out-polygon
       
       orig = [mean(polygon(:,1)),mean(polygon(:,2))];
       
       coords = [polygon(:,1)-orig(1),polygon(:,2)-orig(2)];
       
       %find the angles
       theta = atan2( coords(:,2), coords(:,1) );
       
       [srt,ind] = sort(theta);
       clipPolygon = [polygon(ind,1), polygon(ind,2)];
    end

%% Helper Functions
 
    %computerIntersection() assumes the two lines intersect
    function intersection = computeIntersection(line1,line2)
 
        %this is an implementation of
        %http://en.wikipedia.org/wiki/Line-line_intersection
 
        intersection = zeros(1,2);
 
        detL1 = det(line1);
        detL2 = det(line2);
 
        detL1x = det([line1(:,1),[1;1]]);
        detL1y = det([line1(:,2),[1;1]]);
 
        detL2x = det([line2(:,1),[1;1]]);
        detL2y = det([line2(:,2),[1;1]]);
 
        denominator = det([detL1x detL1y;detL2x detL2y]);
 
        intersection(1) = det([detL1 detL1x;detL2 detL2x]) / denominator;
        intersection(2) = det([detL1 detL1y;detL2 detL2y]) / denominator;
 
    end %computeIntersection
 
    %inside() assumes the boundary is oriented counter-clockwise
    function in = inside(point,boundary)
 
        pointPositionVector = [diff([point;boundary(1,:)]) 0];
        boundaryVector = [diff(boundary) 0];
        crossVector = cross(pointPositionVector,boundaryVector);
 
        if ( crossVector(3) <= 0 )
            in = true;
        else
            in = false;
        end
 
    end %inside
 
%% Sutherland-Hodgman Algorithm
    %ensure the vertices are ordered CCW
    clippedPolygon = subjectPolygon;
    clipPolygon = orderVerticesCCW( clipPolygon );
    numVerticies = size(clipPolygon,1);
    clipVertexPrevious = clipPolygon(end,:);
 
    for clipVertex = (1:numVerticies)
 
        clipBoundary = [clipPolygon(clipVertex,:) ; clipVertexPrevious];
 
        inputList = clippedPolygon;
 
        clippedPolygon = [];
        if ~isempty(inputList),
            previousVertex = inputList(end,:);
        end
 
        for subjectVertex = (1:size(inputList,1))
 
            if ( inside(inputList(subjectVertex,:),clipBoundary) )
 
                if( not(inside(previousVertex,clipBoundary)) )  
                    subjectLineSegment = [previousVertex;inputList(subjectVertex,:)];
                    clippedPolygon(end+1,1:2) = computeIntersection(clipBoundary,subjectLineSegment);
                end
 
                clippedPolygon(end+1,1:2) = inputList(subjectVertex,:);
 
            elseif( inside(previousVertex,clipBoundary) )
                    subjectLineSegment = [previousVertex;inputList(subjectVertex,:)];
                    clippedPolygon(end+1,1:2) = computeIntersection(clipBoundary,subjectLineSegment);                            
            end
 
            previousVertex = inputList(subjectVertex,:);
            clipVertexPrevious = clipPolygon(clipVertex,:);
 
        end %for subject verticies                
    end %for boundary verticies
end %sutherlandHodgman