function centroidpos = computeImageCentroid(I,varargin)
%% Script:  centroidpos = computeImageCentroid(I,varargin)
% Description: 
% 
% INPUTS ----------------------------------------------------------------
% I           = any gray image (may be a binary mask if the centroid of
% the mask is the objective
% varargin - 'PropertyName','PropertyValue'
%   'mask' : string indicating mask type [DEFAULT = 'large']
%   'outorder': string indicating the order of the output [DEFAULT = 'yx']
%     > 'yx':   gives the centroid row and colum position
%     > 'xy':   gives the centroid column and row position
%   'boundingbox': typical bounding box parameters = (TL is top/left pixel)
%       [TL column, TL row, width, height]
%
% OUTPUTS ---------------------------------------------------------------
% centroidpos = [centroidRowY, controidColX] (default)
%
%--------------------------------------------------------------------------
%  Date           Author            E-mail                      Version
%  7  July 2016   Amanda Balderrama amanda.gaudreau@gmail.com     0
%==========================================================================

[M,N] = size(I);

startvarin = 1;
PropertyNames = lower(varargin(startvarin:2:length(varargin)));
PropertyVal = varargin(startvarin+1:2:length(varargin));

if strmatch('mask',PropertyNames)
  binMask = PropertyVal{strmatch('mask',PropertyNames)};
else
  binMask = true(M,N);
end

if strmatch('weighted',PropertyNames)
  weighted = PropertyVal{strmatch('weighted',PropertyNames)};
else
  weighted = 1;
end

if strmatch('outorder',PropertyNames)
  outorder = PropertyVal{strmatch('outorder',PropertyNames)};
else
  outorder = 'yx';
end

if strmatch('bounding',PropertyNames)
  roipos = PropertyVal{strmatch('bounding',PropertyNames)};
else
  roipos = [1,1,N,M];
end
rowsVec = roipos(2):(roipos(2)+roipos(4)-1);
colsVec = roipos(1):(roipos(1)+roipos(3)-1);

% finds the bounding box of ROI ---------------------------------
% maskRows = find(sum(binMask,1) > 0);
% maskCols = find(sum(binMask,2) > 0);
% regprop.BoundingBox = [maskRows(1),maskCols(1),maskRows(end)-maskRows(1), maskCols(end) - maskCols(1)];

% finds centroid (center of mass) of ROI ---------------------------------
if ~weighted && any(binMask(:) == 0)
  maskedI = binMask(rowsVec,colsVec);
else
  maskedI = I(rowsVec,colsVec).*binMask(rowsVec,colsVec);
end
totI = sum(maskedI(:));
[colMesh,rowMesh] = meshgrid(colsVec,rowsVec);
centroidColX = sum(sum(maskedI.*colMesh))/totI;
centroidRowY = sum(sum(maskedI.*rowMesh))/totI;
centroidpos = [centroidRowY,centroidColX];

if strmatch('xy',outorder)
  centroidpos = fliplr(centroidpos);
end

end
