function D2 = ComputeDifferentialOperatorsFromMesh_Neumann_FiniteDifference(GridInfo,Aexch)
%
% D2 = ComputeDifferentialOperatorsFromMesh_Neumann_FiniteDifference(GridInfo,Aexch)
%
% References: "Proposal for a micromagnetic standard problem: domain wall pinning at phase boundaries"
%
% Creates differential operators for a "un"structured mesh with varying exchange
% The mesh is composed by
% N elements (or tiles, or cells)
% K faces 
% 
% INPUT ARGUMENTS:
%
% Volumes,Xel,Yel,Zel, are N times 1 arrays with the corresponding properties 
% of the N elements (volumes, and centers)
%
% Signs is a N times K sparse matrix.
% The entry Signs(n,k) is equal to 1 if the normal to the k-th faces points outwards 
% with respect to the n-th element, and -1 if the normal points inwards.
% If the k-th face is not on the boundary of the n-th element Signs(n,k) is zero.
%
% OUTPUT ARGUMENTS:
% 
% D2 is an N times N sparse matrix performing the derivative div(A*grad(m))
%
%% Unpack the variables
Signs = GridInfo.TheSigns;
Xel = GridInfo.Xel;
Yel = GridInfo.Yel;
Zel = GridInfo.Zel;

%% Dimensions
N = size(Signs,1) ;

el2el=abs(Signs)*abs(Signs'); % compact (nn)

%% W
[n,k] = find(el2el) ;
if ~exist('Aexch','var') %normal method
    % either 1/dx^2, 1/dy^2 or 1/dz^2 since at all times only one of the
    % distances are nonzero.
    w =1./((Xel(n)-Xel(k))+(Yel(n)-Yel(k))+(Zel(n)-Zel(k))).^2;
else % varying Aexch method
    % either Afact/dx^2, Afact/dy^2 or Afact/dz^2 since at all times only 
    % one of the distances are nonzero.
    w =Aexch(k).*Aexch(n)*2./(((Xel(n)-Xel(k))+(Yel(n)-Yel(k))+(Zel(n)-Zel(k))).^2.*(Aexch(n)+Aexch(k)));
end

inds2=[find(diff(k));length(n)];
inds1=[1;inds2+1];
vw=zeros(length(w),1);

%% Main loop
for kk=1:N
    ind=inds1(kk):inds2(kk);
    ii=find(n(ind)==k(ind)); % find central element (element kk)
    for i=1:length(ind)
        if i==ii
            continue
        else
            vw(ind(ii))= vw(ind(ii))- w(ind(i)); % central element
            vw(ind(i)) = vw(ind(i)) + w(ind(i)); % neighbour element
        end
    end
end

D2=sparse(k,n,vw);
end