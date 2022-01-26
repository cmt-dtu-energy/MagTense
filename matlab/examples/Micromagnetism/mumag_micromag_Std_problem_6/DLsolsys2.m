function [Wkfull,Gk,Hk]=DLsolsys2(dxk,dyk,dzk,Wk,scale_bool)
    if nargin==5 && scale_bool
        scalex=10^(round(log10(mean(abs(dxk)))));
        scaley=10^(round(log10(mean(abs(dyk)))));
        scalez=10^(round(log10(mean(abs(dzk)))));
        scales=[scalex;scaley;scalez];
        dxk=dxk/scalex;
        dyk=dyk/scaley;
        dzk=dzk/scalez;
    end
    e=ones(length(dxk),1);
    Gkl1=[e,dxk,dyk,dzk];
    Gk=[Wk*Gkl1; Wk*(dxk.*Gkl1); Wk*(dyk.*Gkl1); Wk*(dzk.*Gkl1)];
    Hk=Wk.*Gkl1';
    Wkfull=Gk\Hk;
    if nargin==5 && scale_bool
        Wkfull=[Wkfull(1,:);Wkfull(2:end,:)./scales];
    end
end