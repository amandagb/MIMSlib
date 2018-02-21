function drgb = intensity2rgb(dmat,colorstr,varargin)
%% function: intensity2rgb
% Description:
% Example:
%     >> 
%
%
% INPUTS ----------------------------------------------------------------
% dmat:    M x N intensity image
% colorstr: 
% linelabels:
% varargin - 'PropertyName','PropertyValue'
%   • 'thresh': 1x2 vector with the [lower threshold, upper threshold]
%   values used to limit the input data
%   • 'i2nan': value (or vector) specifying which value(s) should be set to
%   nan
%
% OUTPUTS ---------------------------------------------------------------
% Calstruct: Stucture with relevant output data
%   • Elnames:
%   • LineSlopeInt_Means: nel x 2 matrix with rows corresponding to the
%   [slope,intercept] values for the element corresponding to that row
%   • LineFitR2_Means = rsqmean;
%   • LineSlopeInt = p;
%   • LineFitR2 = rsq;
% dcellOut:
%   • d_div: L x 1 cell of structures where each structure is the
%        basic MIMS structure associated with the label
%   • line_types
%   • type_rows
%   • CalMasks
%   • dinMasks
%   • ppmVal
%   • CalLabels
%   • dmeans
%
%  Date           Author              E-mail                      Version
%  1  Nov  2015   Amanda Balderrama   amandagbalderrama@gmail.com     3.1


compdir = loaddirfun;

PropertyNames = lower(varargin(1:2:length(varargin)));
PropertyVal = varargin(2:2:length(varargin));

if strmatch('thresh',PropertyNames)
  thresh = PropertyVal{strmatch('thresh',PropertyNames)};
else
  thresh = [nanmin(dmat(:)),nanmax(dmat(:))];
end

if strmatch('i2nan',PropertyNames)
  i2nan = PropertyVal{strmatch('i2nan',PropertyNames)};
else
  i2nan = [];
end

phCmap = custom_colormaps(colorstr,256);
% chanord = 'rgbymcw';
% colr = strfind(chanord,colorstr);
% over_logical = [eye(3);1,1,0;1,0,1;0,1,1;1,1,1]; % R, G, B, Y, M, C, W

dmat(dmat < thresh(1)) = thresh(1);
dmat(dmat > thresh(2)) = thresh(2);
nanind = ismember(dmat,i2nan);
% phCmap = zeros(256,3);
% phCmap(:,find(over_logical(colr,:))) = repmat([0:255]'./255,1,sum(over_logical(colr,:)));%jet(256);
[Lm,Ln] = size(dmat);
phrng = nanmax(dmat(:)) - nanmin(dmat(:));
phnorm = (dmat - nanmin(dmat(:)))./phrng;
drgb = zeros(Lm,Ln,3);
mapind = round(phnorm(:).*(size(phCmap,1)-1)+1);
mapind(isnan(mapind)) = size(phCmap,1) + 1;
phCmap = [phCmap;0,0,0];
del3chv = nan(Lm*Ln,3);
del3chv(nanind(:) == 0,:) = phCmap(mapind(nanind == 0),:);
drgb(:,:,1) = reshape(del3chv(:,1),Lm,Ln);
drgb(:,:,2) = reshape(del3chv(:,2),Lm,Ln);
drgb(:,:,3) = reshape(del3chv(:,3),Lm,Ln);

end
