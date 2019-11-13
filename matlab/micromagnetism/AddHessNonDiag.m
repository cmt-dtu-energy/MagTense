function HessH = AddHessNonDiag(HessH,HessND,admx,admy,admz,bdmx,bdmy,bdmz,NN)
    HessH(1:NN,1:NN) = HessH(1:NN,1:NN) +  HessND.*(admx*(bdmx.')) ;
    HessH(1:NN,NN+1:2*NN) = HessH(1:NN,NN+1:2*NN) +  HessND.*(admx*(bdmy.')) ;
    HessH(1:NN,2*NN+1:end) = HessH(1:NN,2*NN+1:end) +  HessND.*(admx*(bdmz.')) ;
    HessH(NN+1:2*NN,1:NN) = HessH(NN+1:2*NN,1:NN) +  HessND.*(admy*(bdmx.')) ;
    HessH(NN+1:2*NN,NN+1:2*NN) = HessH(NN+1:2*NN,NN+1:2*NN) +  HessND.*(admy*(bdmy.')) ;
    HessH(NN+1:2*NN,2*NN+1:end) = HessH(NN+1:2*NN,2*NN+1:end) +  HessND.*(admy*(bdmz.')) ;
    HessH(2*NN+1:end,1:NN) = HessH(2*NN+1:end,1:NN) +  HessND.*(admz*(bdmx.')) ;
    HessH(2*NN+1:end,NN+1:2*NN) = HessH(2*NN+1:end,NN+1:2*NN) +  HessND.*(admz*(bdmy.')) ;
    HessH(2*NN+1:end,2*NN+1:end) = HessH(2*NN+1:end,2*NN+1:end) +  HessND.*(admz*(bdmz.')) ;
end