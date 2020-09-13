
%By Kaspar K. Nielsen, 20200615, kasparkn@gmail.com
%Method that finds the self-consistent solution to the local magnetic field
%in each tile given the tiles' temperature and state.
%Happ (input) is assumed an (n,3) matrix with the components of the externally applied field
%M_in (input) is the current magnetization (n,3)
%T_in (n,1) is the temperature
%p_in (n,1) is the pressure
%geom is an object of type MagTenseTransientGeometry containng central
%information about the geometry
%setts is an object of type MagTenseTransientSettings containing relevant
%settings for the simulation
%debug is a flag that indicates whether to produce debugging information
%returns H as the resulting field and M the corresponding magnetization.
%Both are (n,3) matrices
%The return value in iterationFlag indicates the status of the iteration
%process
function [H,M_out,iterationFlag] = MagTenseIterateMagneticField( Happ, M_in, T_in, p_in, hyst, geom, setts )

%no. of tiles
n = length(geom.tiles);

%check if the demag tensor has been calculated already. If not, then do it
if isempty( geom.N )
  
    %allocate the n tensors
    geom.N = cell(6,1);
    geom.N{1} = zeros(n,n); %Nxx
    geom.N{2} = zeros(n,n); %Nxy
    geom.N{3} = zeros(n,n); %Nxz
    geom.N{4} = zeros(n,n); %Nyy
    geom.N{5} = zeros(n,n); %Nyz
    geom.N{6} = zeros(n,n); %Nzz
    %loop over each tile to get the N tensor in that tile from all other
    %tiles including it self. Note that in order not to overwhelm the
    %memory, we loop over each tile rather than attempting to get all at
    %once.
    
    %tile centers
    C = MagTenseTilesUtil.getCenter( geom.tiles );
    for i=1:length(geom.tiles)
        
        N = getNTensor_experiment( geom.tiles, C(i,:) ); 
        %distribute N to the various components
        geom.N{1}(i,:) = N( 1, 1:3:(end-2) ); %Nxx
        geom.N{2}(i,:) = N( 1, 2:3:(end-1) ); %Nxy
        geom.N{3}(i,:) = N( 1, 3:3:end ); %Nxz
        
        geom.N{4}(i,:) = N( 2, 2:3:(end-1) ); %Nyy
        geom.N{5}(i,:) = N( 2, 3:3:end ); %Nyz
        
        geom.N{6}(i,:) = N( 3, 3:3:end ); %Nzz
    end
end

%now assume the tensor is there, iterate to find the magnetization
%output field
H = zeros(n,3);

%the field at the previous iteration
H_prev = H;
%dampening factor
lambda = 1.;
err = 2;
cnt = 0;
%initial guess of the magnetization
M_prev = M_in;
while err > 0.001 && cnt < 10
    
    %get the magnetic field
    H(:,1) = geom.N{1} * M_prev(:,1) + geom.N{2} * M_prev(:,2) + geom.N{3} * M_prev(:,3);
    H(:,2) = geom.N{2} * M_prev(:,1) + geom.N{4} * M_prev(:,2) + geom.N{5} * M_prev(:,3);
    H(:,3) = geom.N{3} * M_prev(:,1) + geom.N{5} * M_prev(:,2) + geom.N{6} * M_prev(:,3);
    
    %add the applied field
    H = H + Happ;
    
    %if the iteration is not easy, we can dampen it by a factor less than
    %one
    H = H_prev + lambda * ( H - H_prev );
    
    %get the new magnetization.
    M_out = setts.M( H, T_in, p_in, hyst );
    %calculate the error
    err = max( abs( (sqrt(sum(M_prev.^2,2))-sqrt(sum(M_out.^2,2))) ./ sqrt(sum(M_prev.^2,2)) ) );
    
    %Save the current iteration of the magnetization
    M_prev = M_out;
    
    %increment counter
    cnt = cnt + 1;
end
disp(['err = ' num2str(err) ', cnt = ' num2str(cnt)]);
iterationFlag = 0;
end