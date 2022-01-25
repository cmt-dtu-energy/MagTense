function [v,c,rs,re] = convertToCSR(B)

[c,r,v] = find(B');
rs = [1; find(diff(r) > 0)+1];
re = [find(diff(r) > 0)+1; length(r)+1];

end