function X = myinterp2(I, ratio)
%
% Bilinear interpolation that weights the output linearly proportional to
% the distance from the interpolated point. 
% 
% ratio should be 2
I = upsample(I', ratio);

I = I';

sz = size(I);


for m = 2:2:(sz(2)-1)
    for n = 2:1:(sz(1)-1)

          I(n,m) = .2072*(I(n,m-1) + I(n,m+1) )+... %+ I(n-1,m)+ I(n+1,m)
              .1464*(I(n-1,m-1) + I(n-1,m+1) + I(n+1,m-1) + I(n+1,m+1));
          
    end
end

X = I(2:sz(1)-1,1:sz(2)-1);
X = X'; 
