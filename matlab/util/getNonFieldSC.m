function [cNonField,sNonField] = getNonFieldSC( mat, T )

[J,g,Ns,muB,mu0,kB,rhos,Na,mol,Tdeb,gamma] = getMFTParams( mat );

%T = 1:0.5:350;

cNonField = zeros(length(T),1)';
sNonField = zeros(length(T),1)';
kB * Na/mol;

DebyeInt = @(x) x.^4.*exp(x)./(exp(x)-1.).^2;

for i=1:length(T)
   integ = integral( DebyeInt, 0.001, Tdeb./T(i) );
   cNonField(i) = 9 .* kB .* Na./mol .* ( T(i)./Tdeb ).^3 .* integ; 
end

%add the electronic contribution
cNonField = cNonField + gamma .* T;

%find the non-magnetic part of the entropy from the specific heat
DebyeInt = @(x) x.^3./(exp(x)-1);
for i=1:length(T)
    integ = integral( DebyeInt, 0.001, Tdeb/T(i) );
   sNonField(i) = -3*Na*kB .* ( log(1-exp(-Tdeb./T(i)) ) - 4./(Tdeb/T(i)).^3 * integ );
end
%add the electronic contribution
sNonField = sNonField + gamma .* T;

end

