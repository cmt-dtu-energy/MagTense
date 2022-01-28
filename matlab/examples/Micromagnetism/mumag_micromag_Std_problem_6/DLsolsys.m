function [Wkfull,Gk,Hk]=DLsolsys(dxk,dyk,dzk,Wk,scale_bool)
    if nargin==5 && scale_bool
        scale=10^(round(log10(mean(abs([dxk;dyk;dzk])))));
        dxk=dxk/scale;
        dyk=dyk/scale;
        dzk=dzk/scale;
    end
    e=ones(length(dxk),1);
    Gkl1=[e,dxk,dyk,dzk];
    Gk=[Wk*Gkl1; Wk*(dxk.*Gkl1); Wk*(dyk.*Gkl1); Wk*(dzk.*Gkl1)];
    Hk=Wk.*Gkl1';
    Wkfull=Gk\Hk;
    if nargin==5 && scale_bool
        Wkfull=[Wkfull(1,:);Wkfull(2:end,:)/scale];
    end
end