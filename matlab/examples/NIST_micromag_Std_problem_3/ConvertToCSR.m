function [v,c,rs,re] = ConvertToCSR(B)

% B = sparse(5,5);
% B(1,1) = 1;
% B(1,2) = -1;
% B(1,4) = -3;
% B(2,1) = -2;
% B(2,2) = 5;
% B(3,3) = 4;
% B(3,4) = 6;
% B(3,5) = 4;
% B(4,1) = -4;
% B(4,3) = 2;
% B(4,4) = 7;
% B(5,2) = 8;
% B(5,5) = -5;

[c,r,v] = find(B');
rs = [1; find(diff(r) > 0)+1];
re = [find(diff(r) > 0)+1; length(r)+1];

end