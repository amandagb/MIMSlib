function cmap256 = colormapFromPoints(basecolors)
%% Script:
% Description:
% Example:
% Required Functions:
% INPUTS ----------------------------------------------------------------
% map_name: string indicating name of the colormap that user will use
% n:
% varargin - 'PropertyName','PropertyValue'
%   ''
% OUTPUTS ---------------------------------------------------------------
%
%
%  Date            Author               E-mail                      Version
%   4  Mar  2015   Amanda GBalderrama   amandagbalderrama@gmail.com     1

nclrsgiven = size(basecolors,1);
npntsbtwnclrs = 256/(nclrsgiven-1);
filli = [1:floor(npntsbtwnclrs):(floor(npntsbtwnclrs)*(nclrsgiven-1)),256];
cmap256 = nan(256,3);
cmap256(filli,:) = basecolors;
for i = 1:(nclrsgiven-1)
  p1 = basecolors(i,:);
  p2 = basecolors(i+1,:);
  v = p2 - p1;
  nbtwnp = filli(i+1)-filli(i);
  t = [1:nbtwnp]./nbtwnp;
  pnts = repmat(p1,nbtwnp,1) + repmat(t(:),1,3).*repmat(v,nbtwnp,1);
  cmap256((filli(i)+1):(filli(i+1)-1),:) = pnts(1:end-1,:);
end