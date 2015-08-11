function [seg,varargout] = chanvese(I,varargin)
%% Script:  [seg,varargout] = chanvese(I,varargin)
% Description: This code implements the paper: "Active Contours Without
% Edges" by Chan and Vese for method 'chen', the paper:"Active Contours Without
% Edges for vector image" by Chan and Vese for method 'vector', and the paper
% "A Multiphase Level Set Framework for Image Segmentation Using the
% Mumford and Shah Model" by Chan and Vese.
%   Active contour with Chen-Vese Method
%   for image segementation
%
%   Adaptation of code implemented by Yue Wu (yue.wu@tufts.edu)
%   http://sites.google.com/site/rexstribeofimageprocessing/
% Required Functions: interp_data_elemrng
% Internal Sub-functions: phi0, Heaviside
% INPUTS ----------------------------------------------------------------
% I           = any gray/double/RGB input image
% varargin - 'PropertyName','PropertyValue'
%   'mask' : string indicating mask type [DEFAULT = 'large']
%   'method': string indicating method [DEFUALT = 'chen']
%   'iter': Number of iterations [DEFAULT = 5000]
%   'lambda_1': Weights applied to importance of mean inside the contour
%   [DEFUALT = 1]
%   'lambda_2': Weights applied to importance of mean outside the contour
%   [DEFUALT = 1]
%   'nu': Weighs the importance of area [DEFAULT = 1]
%   'mu': Weighs impact of length of contour [DEFAULT = 0.2]
%   'stepsize': Step size (dt) [DEFAULT = 0.5]
%   'epsilon': value of epsilon for Heaviside function [DEFAULT = 1e-5]
%   'H_type': type of heaviside function to use [DEFUALT = 2]
%   'plot_title': string which would contain description of data being
%   plotted [DEFAULT = '']
%   'verbose': binary indicating whether to give real-time output plot (1)
%   or not (0) [DEFAULT = 1]
%   'singlecomp': binary indicating whether user would like to process
%   segmentation map such that there is only a single component (1) or to
%   output the raw CV outcome (0). The component with the largest number of
%   pixels will be selected as the component of interest [DEFAULT = 0]
%
% LOCATION OF VARIABLES IN ACTIVE CONTOURS WITHOUT EDGES, CHAN, VESE:
% lambda_1, lambda_2, nu, mu : pg 268, F(c1,c2,C)
% Heaviside function: pg 270, H_{1,eps} and H_{2_eps}
%
%   Types of built-in mask functions
%       'small'     = create a small circular mask
%       'medium'    = create a medium circular mask
%       'large'     = create a large circular mask
%       'whole'     = create a mask with holes around
%       'whole+small' = create a two layer mask with one layer small
%                       circular mask and the other layer with holes around
%                       (only work for method 'multiphase')
%       'circles'   = 
%   Types of methods
%       'chen'      = general CV method
%       'vector'    = CV method for vector image
%       'multiphase'= CV method for multiphase (2 phases applied here)
%
% OUTPUTS ---------------------------------------------------------------
% seg        = binary mask where foreground = 1 (original image high
% intesnity, phi0 <= 0) and background = 0 (original image low intensity,
% phi0 > 0) pixels
% varargout - 
%   {1} phi0 = level set corresponding to the final image segmentation
%   {2} regprop = structure which gives the properties of the final
%   segmentation
%           .BoundingBox: bounding box of the segmentation mask
%           .IntCentroid: intensity centroid [row cent, col cent]
%           .BinCentroid: binary centroid [row cent, col cent]
%           .niter: number of iterations to find segementation
%           .t: time to find segmentation
%
%--------------------------------------------------------------------------
%  Date           Author            E-mail                      Version
%  25 Sept 2011   Amanda Gaudreau   amanda.gaudreau@gmail.com     1
%  11 Oct  2011   Amanda Gaudreau   amanda.gaudreau@gmail.com     2
%     Added MIMS specific naming functionality to figures
%  16 Feb  2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     3
%     Added line 222-224 so that dynamic range of the data has a maximum
%     value of 255 based on experiements conduced in Feb 9 notes
%  03 May  2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     4
%     changed mask creation to use initlevelset function & allowed for
%     larger maximum image dimensions
%==========================================================================

%% Initialization of Parameters
tic;
startvarin = 1;
if ~ischar(varargin{1}); 
  phi0 = varargin{1}; 
  startvarin = 2; 
end;

PropertyNames = lower(varargin(startvarin:2:length(varargin)));
PropertyVal = varargin(startvarin+1:2:length(varargin));

% Algorithm Specifictions
if strmatch('iter',PropertyNames)
  num_iter = PropertyVal{strmatch('iter',PropertyNames)};
else
  num_iter = 5000;
end

if strmatch('epsilon',PropertyNames)
  epsilon = PropertyVal(strmatch('epsilon',PropertyNames));
else
  epsilon = 1e-5;
end

if strmatch('h_type',PropertyNames)
  Htype = PropertyVal{strmatch('h_type',PropertyNames)};
else
  Htype =2;
end

% Parameter initializations
if strmatch('mu',PropertyNames)
  mu = PropertyVal{strmatch('mu',PropertyNames)};
else
  mu=0.2;
end

if strmatch('nu',PropertyNames)
  nu = PropertyVal(strmatch('nu',PropertyNames));
else
  nu=0;
end

if strmatch('lambda_1',PropertyNames)
  lambda_1 = PropertyVal{strmatch('lambda_1',PropertyNames)};
else
  lambda_1 = 1;
end

if strmatch('lambda_2',PropertyNames)
  lambda_2 = PropertyVal{strmatch('lambda_2',PropertyNames)};
else
  lambda_2 = 1;
end

if strmatch('stepsize',PropertyNames)
  dt = PropertyVal{strmatch('stepsize',PropertyNames)};
else
  dt = 0.5;
end

if strmatch('singlecomp',PropertyNames)
  onecc = PropertyVal{strmatch('singlecomp',PropertyNames)};
else
  onecc = 0;
end

% Plotting Notes
if strmatch('verbose',PropertyNames)
  v = PropertyVal{strmatch('verbose',PropertyNames)};
else
  v = 1;
end

if strmatch('plot_title',PropertyNames)
  plot_title = PropertyVal{strmatch('plot_title',PropertyNames)};
else
  plot_title = '';
end
fntsz = 14;

title_str2 = sprintf('$$\\lambda_1 = %g, \\lambda_2 = %g, \\mu = %g, \\nu = %g,  dt = %g, H method = %g$$',...
  lambda_1, lambda_2, mu, nu, dt, Htype);
% title_str3 = sprintf('H method = %g, dt = %g',Htype, dt);
%-- End Initialization

%% Initialize input image I and mask
%  resize original image -- Not sure why this is necessary... the largest
%  the minimum size of the image can be is 200 -- note left 9/26
% s = 200./min(size(I,1),size(I,2)); % resize scale
% if s<1
%   I = imresize(I,s);
% end

%% Initialize Mask
tic;
if ~exist('phi0')
  if strmatch('type',PropertyNames)
    type = PropertyVal{strmatch('type',PropertyNames)};
  else
    type = 'circles';
  end
  
  phi0 = initlevelset(I,varargin{:}); 
end
% * in v4 deleted cases for 'chan' (gray) and 'multiphase' images...one
% case now

layer = size(I,3);
P = double(I);

for l = 1:layer
  P(:,:,l) = (255/max(max(P(:,:,l)))).*P(:,:,l);
end

%-- End Initializations on input image I and mask

%--   Core function
%-- Get the distance map of the initial mask
force = eps; % initial force, set to eps to avoid division by zeros
indicator = 0;
%-- End Initialization
phih = [];
if v
  phih = figure('color','w');
  subplot(1,2,1);
end
%-- Main loop
n = 1;
while ~indicator && n <= num_iter
  inidx = find(phi0>=0); % frontground index==> equivalent statement could be find(mask == 0). outside contour is considered to be mask == 0
  outidx = find(phi0<0); % background index ==> equivalent statement could be find(mask == 1). outside contour is considered to be mask == 1
  force_image = 0; % initial image force for each layer
  H = Heaviside(phi0,epsilon,Htype);
  c1 = sum(sum(P.*repmat(H,[1,1,layer]),1),2)/(sum(sum(H))+eps);%(length(inidx)+eps);
  c2 = sum(sum(P.*repmat((1-H),[1,1,layer]),1),2)/(sum(sum(1-H))+eps);%(length(outidx)+eps);
  force_image =  sum(-lambda_1.*(P-repmat(c1,[size(P,1),size(P,2),1])).^2 ...
    + lambda_2.*(P-repmat(c2,[size(P,1),size(P,2),1])).^2,3);
  
  % calculate the external force of the image
  KG = kappa(phi0);
  force = mu*KG./max(max(abs(KG)))+1/layer.*force_image;
  
  % normalized the force
  force = force./max(abs(force(:)));
  
  % get parameters for checking whether to stop
  old = phi0;
  phi0 = phi0+dt.*force;
  new = phi0;
  indicator = checkstop(old,new,dt);
  
  % intermediate output
  if v
    if(mod(n,20) == 0)
      showphi(I,phi0,n);
    end;
  end
  n = n+1;
end;

telaps = toc;
seg = phi0<=0; %-- Get mask from levelset

if v
  showphi(I,phi0,n);
  title_str1 = sprintf('%s Segmentation (%g seconds)',plot_title,telaps);
  hold on;
  %contour(mask0(:,:,1), [0 0], 'g','LineWidth',1);
  op = get(gca,'outerposition');
  set(gca,'OuterPosition',[op(1),op(2),op(3),0.9]);
  %make mask from SDF
  subplot(1,2,2); imagesc(seg);colormap(gray); title('Global Region-Based Segmentation');
  op = get(gca,'outerposition');
  set(gca,'OuterPosition',[op(1),op(2),op(3),0.9]);
  annotation('textbox', [0 0.9 1 0.1], ...%'String', {[title_str1],[title_str2],[title_str3]}, ...
    'String', title_str1, 'EdgeColor', 'none','HorizontalAlignment', 'center',...
    'fontsize',fntsz,'interpreter','none');
  annotation('textbox', [0 0.85 1 0.1], ...%'String', {[title_str1],[title_str2],[title_str3]}, ...
    'String', title_str2, 'EdgeColor', 'none','HorizontalAlignment', 'center',...
    'fontsize',fntsz,'interpreter','latex');
  %annotation('textbox', [0 0.8 1 0.1], ...%'String', {[title_str1],[title_str2],[title_str3]}, ...
  %'String', title_str3, 'EdgeColor', 'none','HorizontalAlignment', 'center',...
  %'fontsize',fntsz,'interpreter','latex');
  %suptitle({[title_str1],[title_str2],[title_str3]},'fontsize',fntsz,'interpreter','latex')
  %subplot(2,2,4); imshow(seg); title('Global Region-Based Segmentation');
end

if onecc
  CC = bwconncomp(seg);
  if CC.NumObjects > 1
    lenCC = cellfun(@(x) length(x),CC.PixelIdxList);
    [v,i] = max(lenCC);
    seg = zeros(size(seg));
    seg(CC.PixelIdxList{i}(:)) = 1;
    phi0 = initlevelset(I,'type',seg);
  end
end

varargout{1} = phi0;
regprop = struct;
% finds the bounding box of ROI ---------------------------------
seg_rows = find(sum(seg,1) > 0);
seg_cols = find(sum(seg,2) > 0);
regprop.BoundingBox = [seg_rows(1),seg_cols(1),seg_rows(end)-seg_rows(1), seg_cols(end) - seg_cols(1)];
% finds centroid (center of mass) of ROI ---------------------------------
maskedI = I.*seg;
totI = sum(maskedI(:));
[c1mesh,r1mesh] = meshgrid(1:size(I,2),1:size(I,1));
meanx = sum(sum(maskedI.*c1mesh))/totI;
meany = sum(sum(maskedI.*r1mesh))/totI;
regprop.IntCentroid = [meany, meanx];

totI = sum(seg(:));
meanx = sum(sum(seg.*c1mesh))/totI;
meany = sum(sum(seg.*r1mesh))/totI;
regprop.BinCentroid = [meany, meanx];
regprop.niter = n;
regprop.t = telaps;
regprop.phih = phih;
varargout{2} = regprop;
end

% ------------------------------------------------------------------
%%   SUBFUNCTIONS
% -------------------------------------------------------------------

function H = Heaviside(z,epsilon,H_type);
%epsilon = 1e-5;
if H_type == 1
  H = zeros(size(z));
  H(z > epsilon) = 1;
  idx = find(z <= epsilon & z >= -epsilon);
  if ~isempty(idx)
    H(idx) = 0.5.*(1 +(z(idx)./epsilon)+sin((pi.*z(idx))./epsilon)./pi);
  end
else
  H = 0.5.*(1+(2/pi).*atan(z./epsilon));
end

end

function KG = kappa(I)
% get curvature information of input image
% input: 2D image I
% output: curvature matrix KG

% Copyright (c) 2009,
% Yue Wu @ ECE Department, Tufts University
% All Rights Reserved

I = double(I);
[m,n] = size(I);
P = padarray(I,[1,1],1,'pre');
P = padarray(P,[1,1],1,'post');

% central difference
fy = P(3:end,2:n+1)-P(1:m,2:n+1);
fx = P(2:m+1,3:end)-P(2:m+1,1:n);
fyy = P(3:end,2:n+1)+P(1:m,2:n+1)-2*I;
fxx = P(2:m+1,3:end)+P(2:m+1,1:n)-2*I;
fxy = 0.25.*(P(3:end,3:end)-P(1:m,3:end)+P(3:end,1:n)-P(1:m,1:n));
G = (fx.^2+fy.^2).^(0.5);
K = (fxx.*fy.^2-2*fxy.*fx.*fy+fyy.*fx.^2)./((fx.^2+fy.^2+eps).^(1.5));
KG = K.*G;
KG(1,:) = eps;
KG(end,:) = eps;
KG(:,1) = eps;
KG(:,end) = eps;
KG = KG./max(max(abs(KG)));
end

function indicator = checkstop(old,new,dt)
% indicate whether we should performance further iteraions or stop

% Copyright (c) 2009,
% Yue Wu @ ECE Department, Tufts University
% All Rights Reserved

layer = size(new,3);

for i = 1:layer
  old_{i} = old(:,:,i);
  new_{i} = new(:,:,i);
end

if layer
  ind = find(abs(new)<=.5);
  M = length(ind);
  Q = sum(abs(new(ind)-old(ind)))./M;
  if Q<=dt*.18^2
    indicator = 1;
  else
    indicator = 0;
  end
else
  ind1 = find(abs(old_{1})<1);
  ind2 = find(abs(old_{2})<1);
  M1 = length(ind1);
  M2 = length(ind2);
  Q1 = sum(abs(new_{1}(ind1)-old_{1}(ind1)))./M1;
  Q2 = sum(abs(new_{2}(ind2)-old_{2}(ind2)))./M2;
  if Q1<=dt*.18^2 && Q2<=dt*.18^2
    indicator = 1;
  else
    indicator = 0;
  end
end
return
end

function D = reinitialization(D,dt)
% reinitialize the distance map for active contour

% Copyright (c) 2009,
% Yue Wu @ ECE Department, Tufts University
% All Rights Reserved

T = padarray(D,[1,1],0,'post');
T = padarray(T,[1,1],0,'pre');
% differences on all directions
a = D-T(1:end-2,2:end-1);
b = T(3:end,2:end-1)-D;
c = D-T(2:end-1,1:end-2);
d = T(2:end-1,3:end)-D;

a_p = max(a,0);
a_m = min(a,0);
b_p = max(b,0);
b_m = min(b,0);
c_p = max(c,0);
c_m = min(c,0);
d_p = max(d,0);
d_m = min(d,0);

G = zeros(size(D));
ind_plus = find(D>0);
ind_minus = find(D<0);
G(ind_plus) = sqrt(max(a_p(ind_plus).^2,b_m(ind_plus).^2)+max(c_p(ind_plus).^2,d_m(ind_plus).^2))-1;
G(ind_minus) = sqrt(max(a_m(ind_minus).^2,b_p(ind_minus).^2)+max(c_m(ind_minus).^2,d_p(ind_minus).^2))-1;

sign_D = D./sqrt(D.^2+1);
D = D-dt.*sign_D.*G;
end

function phih = showphi(I, phi, i)
% show curve evolution of phi

% Copyright (c) 2009,
% Yue Wu @ ECE Department, Tufts University
% All Rights Reserved

for j = 1:size(phi,3)
  phi_{j} = phi(:,:,j);
end
% I(I < 0) = 0;
% I(I > 1) = 1;
% I = sum(I,3);
imagesc(I); colormap(gray);
% imshow(I,'initialmagnification','fit','displayrange',[0 255]);
% imshow(uint8(I))
hold on;

if size(phi,3) == 1
  contour(phi_{1}, [0 0], 'r','LineWidth',2);
  %contour(phi_{1}, [0 0], 'g','LineWidth',1.3);
else
  contour(phi_{1}, [0 0], 'r','LineWidth',2);
  %contour(phi_{1}, [0 0], 'x','LineWidth',1.3);
  %contour(phi_{2}, [0 0], 'g','LineWidth',4);
  %contour(phi_{2}, [0 0], 'x','LineWidth',1.3);
end
hold off;
title([num2str(i) ' Iterations']);
drawnow;
end
