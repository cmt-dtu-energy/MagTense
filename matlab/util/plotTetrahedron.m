

function plotTetrahedron(tile)

inds = [[1,2,3];[1,3,4];[1,2,4];[2,3,4]];

patch('faces',inds,'vertices',tile.vertices','facecolor',tile.color);

end