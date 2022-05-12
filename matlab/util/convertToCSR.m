function [v,c,rs,re] = convertToCSR(B)

[c,r,v] = find(B');
dr = diff(r);
rtmp = repelem(find(dr > 0)+1,dr(dr>0));
rs = [1; rtmp];
re = [rtmp; length(r)+1];

end