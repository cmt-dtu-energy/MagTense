function HessH = AddHessAIP(HessH,dfxdmxdmx,dfydmxdmx,dfzdmxdmx,dfxdmxdmy,dfydmxdmy,dfzdmxdmy,dfxdmxdmz,dfydmxdmz,dfzdmxdmz,...
    dfxdmydmx,dfydmydmx,dfzdmydmx,dfxdmydmy,dfydmydmy,dfzdmydmy,dfxdmydmz,dfydmydmz,dfzdmydmz,...
    dfxdmzdmx,dfydmzdmx,dfzdmzdmx,dfxdmzdmy,dfydmzdmy,dfzdmzdmy,dfxdmzdmz,dfydmzdmz,dfzdmzdmz,Heff,NN)

for i=1:NN
    HessH(i,i) = HessH(i,i) + dfxdmxdmx(i)*Heff(i) + dfydmxdmx(i)*Heff(NN+i) + dfzdmxdmx(i)*Heff(2*NN+i);
    HessH(i,NN+i) = HessH(i,NN+i) + dfxdmxdmy(i)*Heff(i) + dfydmxdmy(i)*Heff(NN+i) + dfzdmxdmy(i)*Heff(2*NN+i) ;
    HessH(i,2*NN+i) = HessH(i,2*NN+i) + dfxdmxdmz(i)*Heff(i) + dfydmxdmz(i)*Heff(NN+i) + dfzdmxdmz(i)*Heff(2*NN+i) ;
    
    HessH(NN+i,i) = HessH(NN+i,i) + dfxdmydmx(i)*Heff(i) + dfydmydmx(i)*Heff(NN+i) + dfzdmydmx(i)*Heff(2*NN+i) ;
    HessH(NN+i,NN+i) = HessH(NN+i,NN+i) + dfxdmydmy(i)*Heff(i) + dfydmydmy(i)*Heff(NN+i) + dfzdmydmy(i)*Heff(2*NN+i) ;
    HessH(NN+i,2*NN+i) = HessH(NN+i,2*NN+i) + dfxdmydmz(i)*Heff(i) + dfydmydmz(i)*Heff(NN+i) + dfzdmydmz(i)*Heff(2*NN+i) ;
    
    HessH(2*NN+i,i) = HessH(2*NN+i,i) + dfxdmzdmx(i)*Heff(i) + dfydmzdmx(i)*Heff(NN+i) + dfzdmzdmx(i)*Heff(2*NN+i) ;
    HessH(2*NN+i,NN+i) = HessH(2*NN+i,NN+i) + dfxdmzdmy(i)*Heff(i) + dfydmzdmy(i)*Heff(NN+i) + dfzdmzdmy(i)*Heff(2*NN+i) ;
    HessH(2*NN+i,2*NN+i) = HessH(2*NN+i,2*NN+i) + dfxdmzdmz(i)*Heff(i) + dfydmzdmz(i)*Heff(NN+i) + dfzdmzdmz(i)*Heff(2*NN+i) ;
end
end