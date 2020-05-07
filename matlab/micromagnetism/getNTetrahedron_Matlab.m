function N_all = getNTetrahedron_Matlab( r, v )
%% calculates and returns the magnetic field from a tetrahedron defined in the vertices v
%% (v(:,1) = v1 etc) size (3,4) evaluated at the points r size (3,n) given M the magnetization vector
%% of the tetrahedron size (3,1)
 N_all = zeros(3*length(r(1,:)),3) ;
%     H = zeros(3,length(r(1,:)));
    for i=1:length(r(1,:))
        N = zeros(3,3);        
        for j=1:4
            k = circshift( v, j-1, 2 );
            [tmp,~] = getN_Triangle( k, r(:,i) );
            N = N + tmp;        
        end
%         H(:,i) = N * M;
N_all(((i-1)*3+1):(i*3),:) = N ;
    end
 
    
end