function Efinal = MRFmotion_detection(Ikm1,Ik,varargin)
%% Script:
% Description:
% Example:
% Required Functions:
% INPUTS ----------------------------------------------------------------
% Ik:   
% Ikm1:
% varargin - 'PropertyName','PropertyValue'
%   'thresh': [DEFUALT = 50]
%   'MtoSratio':  Moving to stationary priors ratio [DEFUALT = 5]
%   'costratio': Cost ratio [DEFUALT = 1]
%   'plot': Indicates whether to plot mask results or not [DEFUALT = 0]
%   'MRForder': Markov random field order [DEFUALT = 2 (8 pixels)]
% OUTPUTS ---------------------------------------------------------------
% Efinal: Final motion mask
%
%  Date           Author            E-mail                      Version
%  13 Feb  2013   Amanda Gaudreau   amanda.gaudreau@gmail.com     1

PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

if strmatch('thresh',PropertyNames)
  gamma = PropertyVal{strmatch('thresh',PropertyNames)};
else
  gamma = 50; % ** THIS PARAMETER REQUIRES TUNING BASED ON DEGREE OF MOTION! MORE MOTION => higher threshold. When frames 565 and 560 were tested, the best performing threshold was ~1500
end

if strmatch('MtoSratio',PropertyNames)
  sigMoverS = PropertyVal{strmatch('MtoSratio',PropertyNames)};
else
  sigMoverS = 5;
end

if strmatch('costratio',PropertyNames)
  costratio = PropertyVal{strmatch('costratio',PropertyNames)};
else
  costratio = 1;
end

if strmatch('MRForder',PropertyNames)
  MRForder = PropertyVal{strmatch('MRForder',PropertyNames)};
else
  MRForder = 2;
end

if strmatch('plot',PropertyNames)
  p = PropertyVal{strmatch('plot',PropertyNames)};
else
  p = 0;
end

RHSconst = log(costratio*sigMoverS);
psi = Ik - Ikm1;
psi_sq = psi.^2;
[M,N] = size(psi_sq);

%% Fixed-threshold hypothesis test
sigS_sq = gamma/(2*RHSconst);
MSmask = psi_sq > gamma;
multmask = double(MSmask);
multmask(MSmask == 0) = 0.5;
multmask(MSmask == 1) = 2;
multmaskimg = multmask.*Ik;
if p; 
  figure; 
  imagesc(MSmask); colormap(gray);
  %imagesc(multmaskimg); colormap(gray); 
  title(sprintf('Fixed threshold HT mask: \\gamma = %d, \\sigma_M/\\sigma_S = %d',gamma, sigMoverS));
end

Eupdate = zeros(M+2,N+2);
Eupdate(1:2:((M+2)*(N+2))) = 1;%(checkerboard(1,(M+2)/2,(N+2)/2) > 0.5);%zeros(M+2,N+2); %
Eupdate(2:end-1,2:end-1) = MSmask;
Eupdatei = ones(size(Eupdate));
Eupdatei(Eupdate == 1) = 0;
Epasti = Eupdatei;

if MRForder == 1
  Nhood = ones(3,3);
  Nhood(1:2:9) = 0;
  nhoodstr = '1-st';
elseif MRForder == 2
  Nhood = ones(3,3);
  Nhood(2,2) = 0;
  nhoodstr = '2-nd';
end

Ediff = 10;
T = 1;
alpha = 2*sigS_sq/T;
n = 1;
Ediffi = zeros(1,20).*1e5;
while Ediff > 5 && n < 20 %&& contflag
  Nsmat = conv2(Eupdatei,Nhood,'same');
  RHS = gamma + alpha.*(2.*Nsmat - sum(Nhood(:)));
  Eupdatei(2:end-1,2:end-1) = psi_sq < RHS(2:end-1,2:end-1);
  Ediffi(n) = sum(sum(abs(Eupdatei - Epasti)));
  Ediff = Ediffi(n);
  if n > 1 && abs(Ediffi(n) - Ediffi(n-1)) == 0
    for i = 2:M+1
      for j = 2:N+1
        Ns = sum(sum(Eupdatei(i-1:i+1,j-1:j+1).*Nhood));
        RHS = gamma + alpha.*(2.*Ns - sum(Nhood(:)));
        Eupdatei(i,j) = psi_sq(i-1,j-1) < RHS;
      end
    end
    Ediffi(n) = sum(sum(abs(Eupdatei - Epasti)));
    Ediff = Ediffi(n);
  end
  Epasti = Eupdatei;
  n = n+1;
end
Efinal = abs(Eupdatei(2:end-1,2:end-1)-1);
%imoverlay('',Efinal,'',Ik,'gray')
%imoverlay('',Efinal{2},'',Ik,'gray')
if p; 
  figure; 
  imagesc(Efinal); colormap(gray);
  %imagesc(multmaskimg); colormap(gray); 
  title(sprintf('MRF adaptive HT mask: %s order',nhoodstr));
end
end
