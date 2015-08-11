function [Imask,ExpParam] = makeTissueMask(I,varargin)
%% function: makeTissueMask
% Description:
%
% INPUTS ----------------------------------------------------------------
% I:    Input image (M x N x 3 uintXX number) OR a string indicating the
%       file name (should include either entire path name or file name with
%       .jpg/.tif/.png file type at the end)
% varargin - 'PropertyName','PropertyValue'
%   • 'Imasks':   structure (see outputs) => if indicates, F/B segmentation
%   procedure will not be executed and the masks in the structure will be
%   used instead
%   • 'dir':      string indicating the parent directory of the image file.
%       should end with a "\". If unspecified compdir.EBfldr is used
%   • 'dsfact':   N x 1 vector indicating the downsampling factors that
%       should be used [DEFUALT = [8,1]']
%   • 'iter':     N x 1 vector indicating the number of iterations that
%       should be performed in the active contours method using the DS-ed
%       image [DEFAULT = [1500,1000]']
%   • 'smoothing':N x 1 vector indicating the smoothing factor to use in
%       the active contours algorithm [DEFAULT = [2,2]']
%   • 'plotimgs': logical indicating whether to plot images or not [DEFAULT = 0]
%   • 'saveimgs': logical indicating whether to save images or not [DEFAULT = 0]
%   • 'graythonly': logical indicating wehther to construct mask using only
%       grayscale threshold technique and now active contours function
%       [DEFAULT = 0]
%
% OUTPUTS ---------------------------------------------------------------
% Imasks: structure with fields for each animal number (as indicated by
%         animals.id)
%   >> an#: structue whose fields are M x N matrices containing the mask
%   for a particular anaomical direction (D/L/R/V) and class
%   (fgndbkgnd,bld,EB,hemo,tsu)
%    - ie: Imasks.an2.D_fgndbkgnd = M x N binary image corresponding to
%    final F/B seg mask
%    - If field is empty, it means no mask was available or produced
%
%  Date           Author              E-mail                      Version
%  14 Jan  2015   Amanda Balderrama   amandagbalderrama@gmail.com     1
%     Based off of masking lines in EB_brainimg_processing function
%  15 Jan  2015   Amanda Balderrama   amandagbalderrama@gmail.com     2
%     Removes background pixels from the image to be processed by the AC
%     algorithm
%  23 Feb  2015   Amanda Balderrama   amandagbalderrama@gmail.com     2.1
%  23 Mar  2015   Amanda Balderrama   amandagbalderrama@gmail.com     2.1

compdir = loaddirfun;

PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

if strmatch('Imasks',PropertyNames)
  Imasks = PropertyVal{strmatch('Imasks',PropertyNames)};
  findmasks = 0;
else
  findmasks = 1;
end

if strmatch('graythonly',PropertyNames)
  gthonly = PropertyVal{strmatch('graythonly',PropertyNames)};
else
  gthonly = 0;
end

if strmatch('dir',PropertyNames)
  imgdir = PropertyVal{strmatch('dir',PropertyNames)};
else
  imgdir = compdir.EBfldr;
end
if strmatch('\',imgdir(end))
  imgdir = strcat(imgdir,'\');
end

if isstr(I)
  if ~strmatch('.',I(end-3))
    Io = imread(imgdir,I,'.tif');
  else
    Io = imread(strcat(imgdir,I));
  end
else
  Io = I;
end

if isa(Io,'double')
  requint = nextpow2(max(Io(:)));
  if requint <= 8
    Io = uint8(Io);
  elseif requint <= 16
    Io = uint16(Io);
  elseif requint <=32
    Io = uint32(Io);
  end
end

if strmatch('plot',PropertyNames)
  plotimgs = PropertyVal{strmatch('plot',PropertyNames)};
else
  plotimgs = 0;
end

if strmatch('dsfact',PropertyNames)
  dsfactv = PropertyVal{strmatch('dsfact',PropertyNames)};
else
  dsfactv = [8,1];
end

if strmatch('iter',PropertyNames)
  iter = PropertyVal{strmatch('iter',PropertyNames)};
else
  iter = [1500,1000];
end

if strmatch('smoothing',PropertyNames)
  smthfact = PropertyVal{strmatch('smoothing',PropertyNames)};
else
  nlevels = length(iter);
  smthfact = repmat(2,1,nlevels);
end

if strmatch('actype',PropertyNames)
  ACalgostr = PropertyVal{strmatch('actype',PropertyNames)};
else
  ACalgostr = 'Edge';
end

Id = double(Io);
[M,N,d] = size(Id);
% ---------------------------------------------------------------------
% • Define input images
% ---------------------------------------------------------------------
ILeq = max(Id,[],3);%Iods;%mean(Id,3);%Id(:,:,3);%double(rgb2gray(Iods));%
[Gx,Gy] = imgradientxy(ILeq,'prewitt'); %figure; imshowpair(Gx,Gy,'montage')
[Gm,Gd] = imgradient(Gx,Gy); %figure; imshowpair(Gm.*IoL_black,Gd.*IoL_black,'montage')

% ---------------------------------------------------------------------
% • Create preliminary mask from intensity threshold of image
% ---------------------------------------------------------------------
% >> fgnd/bkgnd level, gray level effectiveness metric
[FBlvl glvlEM] = graythresh(Io);
% >> initial binary mask from thresholding channels
if gthonly
  scaleFBlvl = 1;
else
  scaleFBlvl = 0.75;
end
M0th = im2bw(Io,FBlvl*scaleFBlvl);
% >> code which looks at the "1" label and redefines the mask so that
% there is only one connected component
M0th_keep = keeplargestCC(M0th);

tinyseR = round(max(M,N)*0.005);
se = strel('rectangle',[tinyseR,tinyseR]);
M0th_keep = imdilate(M0th_keep,se);
M0th_keep = imerode(M0th_keep,se);

% >> Ensure that there is only 1 cc fgnd and one cc bkgnd --- this is
% background validation
M0th_inv = abs(M0th_keep-1);
cc = bwconncomp(M0th_inv);
Lcc = cellfun(@(x) length(x),cc.PixelIdxList);
[leni,keepi] = max(Lcc);
M0th_inv = true(size(M0th_inv));
M0th_inv(cc.PixelIdxList{keepi}) = false;

if ~gthonly
  % ---------------------------------------------------------------------
  % • dilates the cc intensity mask so that the initial contour is
  % guaranteed to be OUTSIDE the tissue boundary (requirment for AC
  % edge algorithm)
  % ---------------------------------------------------------------------
  obj_width = sum(M0th_inv,2);
  obj_height = sum(M0th_inv,1);
  seR = round(max([obj_width(:);obj_height(:)])*0.03);
  se = strel('disk',seR);
  
  Minit = imdilate(M0th_inv,se);
  keeprows = find(sum(Minit,2));%) = 1;
  keeprows = [max([keeprows(1) - 1,1]); keeprows; min([keeprows(end) + 1,M])];
  keepcols = find(sum(Minit,1));%) = 1;
  keepcols = [max([keepcols(1) - 1,1]), keepcols, min([keepcols(end) + 1,N])];
  % keeprows = 1:M;
  % keepcols = 1:N;
  nrows = length(keeprows);
  ncols = length(keepcols);
  padc = round(length(keeprows)/100);
  padr = round(length(keepcols)/100);
  Minit = padarray(Minit(keeprows,keepcols),[padr,padc],'replicate');
  Itest = padarray(Gm(keeprows,keepcols),[padr,padc]);
  
  % ---------------------------------------------------------------------
  % • Perform multi-level active contour segementation on gradient
  % magnitude image using threshold mask for initialization. Parameters
  % of down-sample factor (dsfactv) and iterations could be adjusted
  % for increased speed segementation accuracy
  % ---------------------------------------------------------------------
  dsfactv = dsfactv(:)';
  iter = iter(:)';%[2000,1000,500,200]./2;
  smthfact = smthfact(:)';
  for ell = 1:nlevels
    if ell == 1
      timev = zeros(nlevels,1);
      ACmasks = cell(nlevels,1);
    end
    Iac = imresize(Itest,1/dsfactv(ell));
    if ell == 1
      M0th = imresize(Minit,1/dsfactv(ell),'nearest');
      seR = round(max(size(M0th))*0.05);
      se = strel('disk',seR);
    else
      M0th = imresize(ACmasks{ell-1},size(Iac),'nearest');
      M0th = imdilate(M0th,se);
    end
    
    tic;
    ACmasks{ell} = activecontour(Iac,M0th,iter(ell),ACalgostr,smthfact(ell));
    timev(ell) = toc;
    orig_rows = ceil((padr + (1:dsfactv(ell):nrows))./dsfactv(ell));
    orig_cols = ceil((padc + (1:dsfactv(ell):ncols))./dsfactv(ell));
    M0thnopad = M0th(orig_rows,orig_cols);
    ACnopad = ACmasks{ell}(orig_rows,orig_cols);
  end
  ACnopad = imresize(ACmasks{ell}(orig_rows,orig_cols),[nrows,ncols],'nearest');
  Imask = zeros(M,N);
  Imask(keeprows,keepcols) = ACnopad;
  ExpParam.dsfact = dsfactv;
  ExpParam.iter = iter;
  ExpParam.ACsmoothing = smthfact;
  ExpParam.ACalgo = ACalgostr;
  ExpParam.ACmasks = ACmasks;
  ExpParam.iterTime = timev;
else
  Imask = M0th_inv;
  ExpParam.dsfact = [];
  ExpParam.iter = [];
  ExpParam.ACsmoothing = [];
  ExpParam.ACalgo = 'NONE';
  ExpParam.ACmasks = [];
  ExpParam.iterTime = [];
end

end
