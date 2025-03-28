function [Mx,My,Mz,mx,my,mz] = computeMagneticMomentGeneralMesh(magnetizationsAll,Vols)

arguments
    magnetizationsAll   = [];     %--- The magnetizations
    Vols                = [];     %--- The volumes of theUse CUDA for the calculations
end

if (isempty(Vols))
    %--- The tiles all have the same volume, so we just return the mean
    Mx = mean(magnetizationsAll(:,:,1),2);
    My = mean(magnetizationsAll(:,:,2),2);
    Mz = mean(magnetizationsAll(:,:,3),2);

    mx = [];
    my = [];
    mz = [];
else
    %--- The tiles have different volumes 
    
    sumVol = sum(Vols) ;
    for k=1:size(magnetizationsAll,1) 
        Sigma = magnetizationsAll(k,:).' ;
        NN = round(numel(Sigma)/3) ;
    
        SigmaX = Sigma(0*NN+[1:NN]) ;
        SigmaY = Sigma(1*NN+[1:NN]) ;
        SigmaZ = Sigma(2*NN+[1:NN]) ;
   
        %-- The average magnetization of the system
        Mx(k) = sum(Vols.*SigmaX)./sumVol ;
        My(k) = sum(Vols.*SigmaY)./sumVol ;
        Mz(k) = sum(Vols.*SigmaZ)./sumVol ;
        
        %-- The individual relative magnetic moments of the tiles
        mx(k,:) = (Vols.*SigmaX)./sumVol ;
        my(k,:) = (Vols.*SigmaY)./sumVol ;
        mz(k,:) = (Vols.*SigmaZ)./sumVol ;
    end
    Mx = Mx';
    My = My';
    Mz = Mz';
end