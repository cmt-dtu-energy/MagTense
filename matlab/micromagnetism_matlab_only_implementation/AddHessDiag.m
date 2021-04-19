function HessH = AddHessDiag(HessH,HessDiag,dmx,dmy,dmz,NN)
    HessH(1:NN,1:NN) = HessH(1:NN,1:NN) +  HessDiag.*(dmx*(dmx.')) ;
    HessH(1:NN,NN+1:2*NN) = HessH(1:NN,NN+1:2*NN) +  HessDiag.*(dmx*(dmy.')) ;
    HessH(1:NN,2*NN+1:end) = HessH(1:NN,2*NN+1:end) +  HessDiag.*(dmx*(dmz.')) ;
    HessH(NN+1:2*NN,1:NN) = HessH(NN+1:2*NN,1:NN) +  HessDiag.*(dmy*(dmx.')) ;
    HessH(NN+1:2*NN,NN+1:2*NN) = HessH(NN+1:2*NN,NN+1:2*NN) +  HessDiag.*(dmy*(dmy.')) ;
    HessH(NN+1:2*NN,2*NN+1:end) = HessH(NN+1:2*NN,2*NN+1:end) +  HessDiag.*(dmy*(dmz.')) ;
    HessH(2*NN+1:end,1:NN) = HessH(2*NN+1:end,1:NN) +  HessDiag.*(dmz*(dmx.')) ;
    HessH(2*NN+1:end,NN+1:2*NN) = HessH(2*NN+1:end,NN+1:2*NN) +  HessDiag.*(dmz*(dmy.')) ;
    HessH(2*NN+1:end,2*NN+1:end) = HessH(2*NN+1:end,2*NN+1:end) +  HessDiag.*(dmz*(dmz.')) ;
end