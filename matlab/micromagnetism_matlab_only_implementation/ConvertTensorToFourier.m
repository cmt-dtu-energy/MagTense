function InteractionMatrices = ConvertTensorToFourier(InteractionMatrices,N,dims)

%% Compute Change of Basis Matrix
NN = prod(N) ;
% Do3Dfft  = @(this1DVector)  reshape( fftn(reshape(this1DVector,N(1),N(2),N(3))),NN,1) ;
% Do3Difft = @(this1DVector)  reshape(ifftn(reshape(this1DVector,N(1),N(2),N(3))),NN,1) ;

Do3Dfft   = @(this1DVector)  DoFourier(this1DVector,N,dims,1) ; % direct
Do3Difft  = @(this1DVector)  DoFourier(this1DVector,N,dims,0) ; % inverse


ee = eye(NN) ;
UU = zeros(NN) ;
for n=1:NN
    UU(:,n) =  Do3Dfft(ee(:,n)) ;
end

UU1 = UU'./NN ;
InteractionMatrices.FFT.Do3Dfft = Do3Dfft ;
InteractionMatrices.FFT.Do3Difft = Do3Difft ;
%% Convert Demag Tensor to Fourier Space
 
InteractionMatrices.FFT.DemagTensor = ApplyFourierTransformToTensor(InteractionMatrices.DemagTensor,UU,UU1) ;

%% Filter tensor (and make sparse ?) 


if 0 
    [InteractionMatrices.FFT.DemagTensor,~] = FilterTensor(InteractionMatrices.FFT.DemagTensor,FractionOfEntries) ; % This will be done later in "CalculateInteractionMatrices.m"
    InteractionMatrices.FFT.DemagTensor = ConvertTensorToSparse(InteractionMatrices.FFT.DemagTensor) ;    
end