function [reginfo, fixed_reg, moving_reg] = register_image(photo,MIMS,varargin)
%% Script: [reginfo, fixed_reg ,moving_reg] = config_register(photo,MIMS,varargin)
% [reginfo, fixed_reg ,moving_reg] = EC720reg([],1,'filtered',0,'element','Cu3273',...
%     'fix',1,'resize',1,'plot',1)
% Description:
% Example:
% Required Functions: 
% INPUTS ----------------------------------------------------------------
% photo: either 2 or 3-D matrix corresponding to the photographic image. It
% could also be a string corresponding to the path & location of the image.
% MIMS: metallomic image
% varargin - 'PropertyName','PropertyValue'
%   'filtered': binary 0 or 1 to apply 'fix_sngl_pxl'
%   'element': string associated with the element that will be used
%   [DEFAULT = 'Cu3273']
%   'fix': binary indicating whether to fix photo (0) or MSI (1) [DEFAULT = 1]
%   'resize': binary indicating whether to resize images to match
%   dimensions of photo (0) or MSI (1) [DEFAULT = 0]
%   'plot': binary indicating whether to plot (1) or not (0) [DEFAULT = 1]
%   'optimization_mode': string indicate whether optimizer should be set up
%   for monomodal or multimodal registrtion [DEFAULT = 'monomodal']. See
%   MATLAB Help, imregconfig, for more information
%   'lpf': string indicating whether to lpf the MIMS images ('mi'), the
%   photo ('ph') or both ('any string') with a gaussian lpf [DEFAULT = NO
%   FILTERING]
%   'lpf_size': one or two element vector indicating window size of the 2D
%   filter [DEFAULT = [20,20]]
%   'lpf_sigma': integer indicating standard deviation of the filter
%   [DEFAULT = 5]
% OUTPUTS ---------------------------------------------------------------
%
%
%  Date           Author            E-mail                      Version
%  09 July 2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     1
%     Adapted from EC720reg

PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));
notMSIflg = 0;
if ischar(MIMS)
  d = readfilesOES('comp','pc','dat',MIMS);
  
  if strmatch('filt',PropertyNames)
    filt01 = PropertyVal{strmatch('filtered',PropertyNames)};
    if filt01
      d = fix_sngl_pxl(d,[]);
    end
  end
else
  d.Cu3273 = MIMS;
  d.dataset = 'USRimg';
  notMSIflg = 1;
  MIMS = 0;
end

if strmatch('element',PropertyNames)
  elem = PropertyVal{strmatch('element',PropertyNames)};
else
  elem = 'Cu3273';
end

MI = getfield(d,elem);

if ischar(photo)
  PI = double(rgb2gray(imread(photo)));
  slashes = strfind(photo,'\');
  PIname = photo(slashes(end) + 1:end);
elseif sum(size(photo)) > 2
  PI = photo;
  PIname = 'USRimg';
end

if strmatch('lpf',PropertyNames) % default is to resize MIMS to photo
  lpf_control = PropertyVal{strmatch('lpf',PropertyNames)}; 
  if strmatch('lpf_size',PropertyNames)
    lpf_size = PropertyVal{strmatch('lpf_size',PropertyNames)};
    if sum(size(lpf_size)) == 2
      lpf_size = [lpf_size,lpf_size];
    end
  else
    lpf_size = [20,20];
  end
  
  if strmatch('lpf_sigma',PropertyNames)
    lpf_sigma = PropertyVal{strmatch('lpf_sigma',PropertyNames)};
  else
    lpf_sigma = 5;
  end
  
  gfilt = fspecial('gaussian',lpf_size,lpf_sigma);
  if strmatch('ph',lpf_control)
    PI = filter2(gfilt,PI);
  elseif strmatch('mi',lpf_control)
    MI = filter2(gfilt,MI);
  else
    PI = filter2(gfilt,PI);
    MI = filter2(gfilt,MI);
  end
end

reginfo.MIMS.size = size(MI);
reginfo.MIMS.original = MI;
reginfo.MIMS.fname  = d.dataset;

reginfo.photo.fname = PIname;
reginfo.photo.original = PI;
reginfo.photo.size = size(PI);

if strmatch('plot',PropertyNames)
  p = PropertyVal{strmatch('plot',PropertyNames)};
else
  p = 1;
end

if strmatch('resize',PropertyNames) % default is to resize MIMS to photo
  if ~PropertyVal{strmatch('resize',PropertyNames)};
    MI = imresize(MI,size(PI),'bicubic');
  else
    PI = imresize(PI,size(MI),'bicubic');
  end
else
  MI = imresize(MI,size(PI),'bicubic');
end
reginfo.scaled2 = size(MI);
fixed = MI;
moving = PI;
reginfo.fixed = 'Metallomic';%d.dataset;
reginfo.moving = 'Photo';%PIname;

if strmatch('fix',PropertyNames)
  if ~PropertyVal{strmatch('fix',PropertyNames)};
    fixed = PI;
    moving = MI;
    reginfo.fixed = 'Photo';%PIname;
    reginfo.moving = 'Metallomic';%d.dataset;
  end
end

if strmatch('opt',PropertyNames)
  opt_mode = PropertyVal{strmatch('opt',PropertyNames)};
else
  opt_mode = 'monomodal';%'multimodal';
end

[optimizer, metric]  = imregconfig(opt_mode);
optimizer.MaximumIterations = 1000;
movingN = (moving-min(moving(:)))./(max(moving(:)) - min(moving(:))); % normalizes image to values in 0-1
fixedN = (fixed-min(fixed(:)))./(max(fixed(:)) - min(fixed(:))); % normalizes image to values in 0-1

tic;
[moving_reg,T] = imregister(fixed,moving,'affine',optimizer,metric,'DisplayOptimization',0); % Im_reg is the product of applying t to Im using imtransform
endt = toc; 

[ysize,xsize] = size(moving);
x_center = xsize/2;
y_center = ysize/2;
t = reshape(T(1:6),1,6);
iT = inv(T); iT(7:9) = [0;0;1];
it = reshape(iT(1:6),1,6);

%% Recalcuate geometric parameters based on optimal transform determined in imregister
theta = atan(t(4)/t(5));
tc = cos(theta); ts = sin(theta);
sx = t(5)/tc;
k = (t(2)*tc + t(1)*ts)/(t(1)*tc - t(2) * ts);
sy = t(1)/(tc + k*ts);
tx = (2*t(3) - xsize + xsize*sy*(tc + k*ts) - ysize*sy*(ts - k*tc))/2;
ty = (2*t(6) - ysize + xsize*sx*ts + ysize*sx*tc)/2;
alpha_est = [tx ty theta sx sy k];
[Aest,Mest,iMest] = alpha2tmat(alpha_est,ysize,xsize);

moving_reg = transform_image(moving,Mest,'transtype','trans');%,'rotpnt','tl');
movingN_reg = transform_image(movingN,Mest,'transtype','trans');%,'rotpnt','tl');

fixed_reg = transform_image(fixed,iMest,'transtype','trans');%,'rotpnt','tl');
fixedN_reg = transform_image(fixedN,iMest,'transtype','trans');%,'rotpnt','tl');

MIb4 = mi(fixed,moving);
MIafter = mi(fixed_reg,moving);
CCb4 = corr2(fixed,moving);
CCafter = corr2(fixed_reg,moving);
disp(sprintf('Mutual Information: before = %1.3f; after = %1.3f \nCross Correlation: before = %1.3f; after = %1.3f',MIb4, MIafter, CCb4, CCafter));
ae = alpha_est; ae(3) = alpha_est(3)*180/pi; 
disp(sprintf('tx = %1.2f, ty = %1.2f, angle = %1.2f, sx = %1.2f, sy = %1.2f, k = %1.2f',ae));

reginfo.fixed_img = fixed;
reginfo.moving_reg = moving_reg;
reginfo.reg_time = endt;
reginfo.alpha_est = alpha_est;
reginfo.MIbefore = MIb4;
reginfo.MIafter = MIafter;
reginfo.NCCbefore = CCb4;
reginfo.NCCafter = CCafter;

if p
  if strmatch('Metallomic',reginfo.fixed)
    th = auto_thresh(fixed,'auto',[]);
    fixed(fixed >= th(2)) = th(2); fixed(fixed <= th(1)) = th(1);
  else
    th = auto_thresh(moving,'auto',[]);
    moving(moving >= th(2)) = th(2); moving(moving <= th(1)) = th(1);
  end
  imoverlay1x4(fixed, reginfo.fixed(1), moving, moving_reg, reginfo.moving(1))
  set(gcf,'PaperPositionMode','auto'); 
end

end
