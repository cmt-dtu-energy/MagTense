function [Mx,My,Mz,mx,my,mz] = computeMagneticMomentGeneralMesh(magnetizationsAll,Vols,Ms)

arguments
    magnetizationsAll   = [];     %--- The magnetizations
    Vols                = [];     %--- The volumes of the tiles used in the micromagnetism model
    Ms                  = 1 ;     %--- The saturation magnetization
end

if (isempty(Vols))
    %--- The tiles all have the same volume, so we just return the mean

    if (length(magnetizationsAll(:,1,1,1)) == 2)
        %--- We have a hysteresis problem
        for i = 1:length(magnetizationsAll(1,1,:,1))
            Mx(i) = mean(Ms'.*magnetizationsAll(end,:,i,1));
            My(i) = mean(Ms'.*magnetizationsAll(end,:,i,2));
            Mz(i) = mean(Ms'.*magnetizationsAll(end,:,i,3));
        end
    else
        %--- We have a time-varying problem
        for i = 1:length(magnetizationsAll(:,1,1,1))
            Mx(i) = mean(Ms'.*magnetizationsAll(i,:,1,1));
            My(i) = mean(Ms'.*magnetizationsAll(i,:,1,2));
            Mz(i) = mean(Ms'.*magnetizationsAll(i,:,1,3));
        end
    end

    mx = [];
    my = [];
    mz = [];
else
    %--- The tiles have different volumes 
    
    sumVol = sum(Vols) ;

    if (length(magnetizationsAll(:,1,1,1)) == 2)
        %--- We have a hysteresis problem
        for i = 1:length(magnetizationsAll(1,1,:,1))
            SigmaX(i,:) = Ms'.*magnetizationsAll(end,:,i,1);
            SigmaY(i,:) = Ms'.*magnetizationsAll(end,:,i,2);
            SigmaZ(i,:) = Ms'.*magnetizationsAll(end,:,i,3);
        end
    else
        %--- We have a time-varying problem
        for i = 1:length(magnetizationsAll(:,1,1,1))
            SigmaX(i,:) = Ms'.*magnetizationsAll(i,:,1,1);
            SigmaY(i,:) = Ms'.*magnetizationsAll(i,:,1,2);
            SigmaZ(i,:) = Ms'.*magnetizationsAll(i,:,1,3);
        end
    end

    for k=1:size(SigmaX,1) 
        %-- The average magnetization of the system
        Mx(k) = sum(Vols.*SigmaX(k,:)')./sumVol ;
        My(k) = sum(Vols.*SigmaY(k,:)')./sumVol ;
        Mz(k) = sum(Vols.*SigmaZ(k,:)')./sumVol ;
        
        %-- The individual relative magnetic moments of the tiles
        mx(k,:) = (Vols.*SigmaX(k,:)')./sumVol ;
        my(k,:) = (Vols.*SigmaY(k,:)')./sumVol ;
        mz(k,:) = (Vols.*SigmaZ(k,:)')./sumVol ;
    end
    Mx = Mx';
    My = My';
    Mz = Mz';
end