function thisVectorFFT = DoFourier(this1DVector,N,dims,DirectOrInverse)
% get column vector "this1DVector" as input
% reshape to a 3D vector using sizes specified by "N"
% apply Fourier transform along dimensions specified
% by the 1x3 array of booleans "dims"
% reshape to a column vector
% (if DirectOrInverse==0 apply inverse Fourier transform instead)

NN = prod(N) ;



if all(dims==1)
    if DirectOrInverse
        thisVectorFFT = reshape(fftn(reshape(this1DVector,N(1),N(2),N(3))),NN,1) ;
    else
        thisVectorFFT = reshape(ifftn(reshape(this1DVector,N(1),N(2),N(3))),NN,1) ;
    end
else
    if DirectOrInverse
        ThatFourier = @(this1DVectorReshaped,thedim) fft(this1DVector,[],thedim) ;
    else
        ThatFourier = @(this1DVectorReshaped,thedim) ifft(this1DVector,[],thedim) ;
    end
    thisReshape = reshape(this1DVector,N(1),N(2),N(3)) ;
    if dims(1)
        thisReshape = ThatFourier(thisReshape,1) ;
    end
    if dims(2)
        thisReshape = ThatFourier(thisReshape,2) ;
    end
    if dims(3)
        thisReshape = ThatFourier(thisReshape,3) ;
    end
    
    thisVectorFFT = reshape(thisReshape,NN,1) ;
end