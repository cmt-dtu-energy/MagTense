function val = getMicroMagDemagApproximation( mode )

    val = int32(1);
    
    switch mode
        case 'none'
            %no approximation applied to the demag tensor
            val = int32(1);
        case 'threshold'
            %cut-off all values below threshold and make the matrix sparse
            val = int32(2);
        case 'fft_thres'
            %apply the threshold in fourier-space through:
            % NP = FT * N * IFT
            % NP(NP<threshold) = 0
            % NP = sparse(NP)
            %H = IFT ( NP * FT * M )
            %with FT = fft( eye(n,n) ) and IFT = ifft( eye(n,n) )
            val = int32(3);
        case 'threshold_fraction'
            %cut-off all values below a certain fraction specified by threshold and make the matrix sparse
            val = int32(4);
        case 'fft_threshold_fraction'
            %cut-off all values below a certain fraction specified by threshold and make the matrix sparse
            val = int32(5);
    end
end