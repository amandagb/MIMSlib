function [detectedXY,objDetails] = trackHSVkalmanFilter(I,varargin)
%% Script:  objDetails = detectHSVobject(I,varargin)
% Description: 
% 
% INPUTS ----------------------------------------------------------------
% I:  
% varargin - 'PropertyName','PropertyValue'
%   'seedLocation': column and row (x and y) position of the nearest
%     predicted location of the object. If threshold segmentation is
%     performaned and there is more than 1 component, the component closest
%     to the seed location will be selected
%   'BoundingBox': typical bounding box parameters = (TL is top/left pixel)
%       [TL column, TL row, width, height] 
%   'segmentMethod': string indicating method to use for identifying the mask
%     > 'CV': use CVcentroidPerFrame function
%     > 'thresh': use intensity threshold only
%     > 'kmeans': use kmeans function with two modes
%   'searchsize': integer indicating how many pixels to search surrounding
%     the specified location
%   'intensityThresh': value indicating the intensity threshold [DEFAULT =
%   mean]
%   'minCCsize': integer indcating the smallest allowed CC to be identified
%   as the object of interest  [ DEFAULT = 5]
%   
%
% OUTPUTS ---------------------------------------------------------------
% objDetails:   structure with relevant details regarding the location of
% the detected object
%
%--------------------------------------------------------------------------
%  Date           Author            E-mail                      Version
%  11 July 2016   Amanda Balderrama amanda.gaudreau@gmail.com     0
%==========================================================================

[M,N,d] = size(I);

param = getDefaultParameters;
if size(I,3) > 1
  I = rgb2gray(I);
end

startvarin = 1;
PropertyNames = lower(varargin(startvarin:2:length(varargin)));
PropertyVal = varargin(startvarin+1:2:length(varargin));

if strmatch('seed',PropertyNames)
  param.seedPosition = PropertyVal{strmatch('seed',PropertyNames)};
end

if strmatch('bounding',PropertyNames)
  param.BoundingBox = PropertyVal{strmatch('bounding',PropertyNames)};
end

if strmatch('search',PropertyNames)
  param.pxlborder = PropertyVal{strmatch('search',PropertyNames)};
end

if strmatch('seg',PropertyNames)
  param.segMethod = PropertyVal{strmatch('seg',PropertyNames)};
end

if strmatch('int',PropertyNames)
  param.Ith = PropertyVal{strmatch('int',PropertyNames)};
end

if strmatch('minCC',PropertyNames)
  param.minCCsize = PropertyVal{strmatch('minCC',PropertyNames)};
end

if ~isempty(param.seedPosition) && param.pxlborder
  maskArea = [floor(param.seedPosition(1)), floor(param.seedPosition(2)),round(param.minCCsize/2),round(param.minCCsize/2)] ...
    + param.pxlborder.*[-1,-1,2,2];
elseif param.BoundingBox(3) ~= N || param.BoundingBox(4) ~= M
  maskArea = floor(param.BoundingBox) + param.pxlborder.*[-1,-1,2,2];
else
  maskArea = floor(param.BoundingBox);
  rowsVec = maskArea(2):(maskArea(2)+maskArea(4)-1);
  rowsVec = rowsVec(rowsVec >= 1 & rowsVec <= M);
  colsVec = maskArea(1):(maskArea(1)+maskArea(3)-1);
  colsVec = colsVec(colsVec >= 1 & colsVec <= N);
  
  M0 = false(M,N);
  M0(rowsVec,colsVec) = I(rowsVec,colsVec) >= param.Ith;
  ccM0 = bwconncomp(M0);
  rpM0 = regionprops(M0);
  ccM0.Area = cellfun(@(x) length(x),ccM0.PixelIdxList);
  validind = find(ccM0.Area > param.minCCsize);
  [~,compoi] = min(ccM0.Area(validind));
  compoi = validind(compoi);
  roipos = floor(rpM0(compoi).BoundingBox);
  bbwidth = round(roipos(3)*1.5);
  bbheight = round(roipos(4)*1.5);
  
  maskArea = roipos + [-bbwidth,-bbheight,2*bbwidth,2*bbheight] ...
    + param.pxlborder.*[-1,-1,2,2];
end
rowsVec = maskArea(2):(maskArea(2)+maskArea(4)-1);
rowsVec = rowsVec(rowsVec >= 1 & rowsVec <= M);
colsVec = maskArea(1):(maskArea(1)+maskArea(3)-1);
colsVec = colsVec(colsVec >= 1 & colsVec <= N);

switch param.segMethod
  case 'thresh'
    M0 = false(M,N);
    M0(rowsVec,colsVec) = I(rowsVec,colsVec) >= param.Ith;
    objDetails.foregroundMask = M0;
  case 'kmeans'
    imgpxls = reshape(I(rowsVec,colsVec),[],1);
    klabels = kmeans(imgpxls,2,'Replicates',10);
    kmRegionMask = reshape(klabels,maskArea(4),maskArea(3)) - 1;
    if kmRegionMask(1,1) == 1
      kmRegionMask = abs(kmRegionMask - 1);
    end
    % Debugging images
    %figure('position',[680,680,1000,300]); subplot(1,3,1); imagesc(I(rowsVec,colsVec)); subplot(1,3,2); imagesc(kmRegionMask); subplot(1,3,3); hist(reshape(I(rowsVec,colsVec),[],1),50);
    kmMask = zeros(M,N);
    kmMask(rowsVec,colsVec) = kmRegionMask;
    objDetails.foregroundMask = logical( kmMask );
  otherwise
    CVout = CVcentroidPerFrame(I,maskArea);
    objDetails.foregroundMask = logical( CVout.foregroundMask );
end

maskCC = bwconncomp(objDetails.foregroundMask);
if maskCC.NumObjects > 1
  [~,compoi] = max(cellfun(@(x) length(x),maskCC.PixelIdxList));
  M0 = false(M,N);
  M0(maskCC.PixelIdxList{compoi}) = true;
  objDetails.foregroundMask = M0;
end
objDetails.foregroundMinIntensity = min(I(objDetails.foregroundMask == 1));
objDetails.BoundingBox = getfield(regionprops(objDetails.foregroundMask,'BoundingBox'),'BoundingBox');

detectedXY = computeImageCentroid(I,'mask',objDetails.foregroundMask,...
  'outorder','xy');

  function param = getDefaultParameters
    param.seedPosition          = [];
    param.pxlborder             = 0;
    param.BoundingBox           = [1,1,N,M];
    param.segMethod             = 'kmeans';
    lowth = quantile(I(:),0.99);
    upth = max(I(:));
    param.Ith                   = mean([lowth,upth]);
    param.minCCsize = 5;
    
    param.segmentationThreshold = 0.05;
    param.movementThreshold = 1;
    
    param.frameRange = [500,3e3]; % starts tracking procedure after 5 ms - 30 ms. entire duration of HSV experiments is typically 250 msec or 25000 frames
  end
end
