

function [sph] = sph_field_fit( B, pts, L, R0 )

%convert coordinates so that z is the field direction
pts(:,:) = pts(:,[2,3,1]);

%get the spherical coordinates of the evaluation points
%[phi,theta,r] = cart2sph( pts(:,1), pts(:,2), pts(:,3) );
r = sqrt( sum(pts.^2,2) );
theta = acos( pts(:,3) ./ r );
phi = atan2(  pts(:,2), pts(:,1) );

%no. of points the field is evaluated at
np = length(B);

%no. of harmonics to decompose into
nh = 1 + 3;
for l=1:L
    nh = nh + l+1;
end

%%Construct the rectangular matrix of spherical harmonics
S = zeros( np, nh );
orders = zeros( 3, nh );
%loop over each order
ind_cnt = 1;
for l=0:L+3
    %associated Legendre polynomials
   Plm = legendre( l, cos( theta ) );
   
   
   S(:,ind_cnt) = legscale(l,0) .* (r./R0).^l .* Plm(1,:)';
  
   
   orders(:,ind_cnt) = [l,0,0];
   ind_cnt = ind_cnt + 1;
end

for l=1:L
   Plm = legendre( l, cos( theta ));
   for m=1:l
   
       S(:,ind_cnt) = legscale(l,m) .* (r./R0).^l .* Plm(m+1,:)' .* cos( m * phi );%exp(1i*m*phi);
       orders(:,ind_cnt) = [l,m,0];
       ind_cnt = ind_cnt + 1;
       
       
       S(:,ind_cnt) = legscale(l,m) .* (r./R0).^l .* Plm(m+1,:)' .* sin( m * phi );%exp(1i*-m*phi);
       orders(:,ind_cnt) = [l,m,1];
       ind_cnt = ind_cnt + 1;
   end
end
%S=real(S);
C = S\B;

order_labels = cell(length(C),1);
for i=1:length(C)
    order_labels{i} = char([num2str(orders(1,i)) ',' num2str(orders(2,i)) ',' num2str(orders(3,i))]);
end

sph = struct('C',C,'orders',orders);
sph.labels = order_labels;
end

function C = legscale(l,m)
   tmp = (l-m+1):2.0:(l+m);   
   %C = (-1.0)^m./ prod( tmp );
   C = 1./ prod( tmp );
end
   