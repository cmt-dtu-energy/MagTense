function stateFunction = getDefaultMagStatStateFunction()

    stateFunction = struct();
    stateFunction.nT = int32(1);
    stateFunction.nH = int32(1);
    stateFunction.M = zeros(1,1);
    stateFunction.T = zeros(1,1);
    stateFunction.H = zeros(1,1);

end