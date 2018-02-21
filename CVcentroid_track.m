function [segarea] = CVcentroid_track(ana_frame);
%% Script:  CVcentroid_track --- not a function in v1, v3
% Description: <E-mailed to Chad, Andy, Noel 2013 Feb 12 04:00>
% CVcentroid_track will prompt for user interaction when necessary (only at
% the very beginning). INSTRUCTIONS FOR USER INTERACTION ARE GIVEN IN THE
% MATLAB COMMAND WINDOW. While a message box may have been more obvious,
% It's easier for now to display in command window.
%
% If you would like live feedback, set db = 1. This will plot everything
% for you (output from CV algo, real time tracking result). If you would
% like to save the image sequence with the centroid pixel indicated, set
% save_vid = 1 and define your savepath. Note that some variables are
% artifact variables which I used for a previous tracking algorithm and
% which may be used again in the future (ie enhancefact, fadefact,
% upsampfact, ect). Please refer to comments for more information on other
% variables.
%
% INPUTS ----------------------------------------------------------------
% varargin - 'PropertyName','PropertyValue'
%   'path':  Path of parent directory of tif images to be analyzed. 
%
% OUTPUTS ---------------------------------------------------------------
% 
%--------------------------------------------------------------------------
%  Date           Author            E-mail                      Version
%  12 Feb  2013   Amanda Gaudreau   amanda.gaudreau@gmail.com     1
%  14 Feb  2013   Amanda Gaudreau   amanda.gaudreau@gmail.com     2
%  20 July 2013   Amanda Gaudreau   amanda.gaudreau@gmail.com     3
%   Version used by Andy Fisher for tracking blast experiment videos. Most
%   notable change (perhaps only) is the file read in method from .tif
%   series to multipage tif method. 
%  21 July 2013   Amanda Gaudreau   amanda.gaudreau@gmail.com     3.1
%   Changed conditional statement on line 172 so difference in centroid
%   mean must be within +/- 1stdev rather than +/- 0.5stdev
%==========================================================================

set(0,'DefaultAxesFontSize',14,'defaultfigurecolor','w',...%'DefaultAxesFontWeight','bold',...
  'DefaultLineLineWidth',2, 'Defaultaxesposition','remove',...
  'defaultAxesFontName','Arial');

% find image headers and file properties ---------------------------------
if exist('fpanel'); clear('fpanel'); end
path1 = 'C:\Users\ADGB\Documents\Dropbox\MADLab Research\Data\Impacter Nose Tracking\20130605_0910_HSV69_B1_C1_convert\20130605_0910_HSV69_B1_C1_convert';
if ~exist('path1') || isempty(path1)
%   path1 =
%   uigetdir('C:\Users\Andy\Documents\Research\MatLab\01.25.13\','Select the home directly of the blast images you would like to analyze');
  path1 = uigetdir('C:\Users\ADGB\Documents\Dropbox\MADLab Research\Data\Impacter Nose Tracking\20130605_0910_HSV69_B1_C1_convert','Select the home directly of the blast images you would like to analyze');
  path1 = strcat(path1,'\');
end

fcount = 0;
info = imfinfo([path1 '.tif']);
num_images = numel(info);

% I1orig = double(rgb2gray(imread([path1 '.tif'],1,'Info',info)));
I1orig = double(imread([path1 '.tif'],1,'Info',info));
[imgR,imgC] = size(I1orig);

% define parameters & initialize variables -------------------------------
upsampfact = 1; % upsampling factor determines the resolution of the velocity vectors
lam1 = 1;
lam2 = 1;
mu = 0.01*255^2;
nCViter = 1000;
enhancefact = 1.1; % enhancement factor for foreground labled pixels
fadefact = 0.9; % fade factor for background labeled pixels
roi_border = 4; % number of pixels (in original image dimensions) around roi
maxf2fmove = 5; % maximum frame to frame movement
nframe = 0;
db = 0;
save_vid = 0;
savepath = 'C:\Users\ADGB\Documents\Dropbox\MADLab Research\Data\Impacter Nose Tracking\20130605_0910_HSV69_B1_C1_convert';
t1frame = 1; %ana_frame(1);
f = 1;
seg2 = 'small';
centpxl = zeros(length(25501),2); % gives the row and column (which is the y and x coordinate)
segarea = zeros(length(25501),3); % [1: total area of mask, 2: mean of pixels in mask, 3: standard dev of pixels in mask]

% User should select ---------------------------------
% roifig = figure; set(roifig,'toolbar','figure');
% h1 = imagesc(I1orig);colormap(gray);
% h = uicontrol('String','Continue',...
%   'Callback','uiresume(gcbf)');
% disp('Zoom into region of image where the nose is & press Continue');
% uiwait(gcf);
% h = uicontrol('String','Continue',...
%   'Callback','uiresume(gcbf)');
% disp('Using cross-hairs, draw a rectangle around the nose, leaving some bordering dark pixels & press Continue');
% rh = imrect;
% uiwait(gcf);
% roipos = getPosition(rh);
% close(roifig); clear h
roipos = [ 197    80    18    14];%round(roipos);
roipos(3:4) = roipos(3:4) + 1.*rem(roipos(3:4),2); % <= makes even; 0.^(rem(roipos(3:4),2)); <= makes odd
roi_rows = round(roipos(2)):round((roipos(2)+roipos(4)));
roi_cols = round(roipos(1)):round((roipos(1)+roipos(3)));

% define image block and segment ---------------------------------
I1roi = I1orig(roi_rows,roi_cols);
[seg1,phi_seg1,prop1] = chanvese(I1roi, 'type','small', 'singlecomp',1,'verbose',0,...
  'lambda_1', lam1, 'lambda_2', lam2, 'mu', mu, 'iter', 1000);
bkgndind = (seg1 == 0);
forgndind = (seg1 == 1);
centpxl(1,:) = [roi_rows(floor(prop1.IntCentroid(1)))+rem(prop1.IntCentroid(1),1), ... %[roi_rows(round(prop1.IntCentroid(1))), ...
  roi_cols(floor(prop1.IntCentroid(2)))+rem(prop1.IntCentroid(2),1)]; %roi_cols(round(prop1.IntCentroid(2)))];
bndbx1 = prop1.BoundingBox;
% redefine roi rows and cols based on center of mass --------------------
roi_rows = max(1,roi_rows(bndbx1(2)) - roi_border) + ...
  [0:(bndbx1(4) + 2*roi_border)];
roi_cols = max(1,roi_cols(bndbx1(1)) - roi_border) + ...
  [0:(bndbx1(3) + 2*roi_border)];
borderoffset = round(max(roipos(3:4))/2) + maxf2fmove;

if db
  trackfig = figure;
  imagesc(I1orig); colormap(gray); hold on;
  plot(centpxl(1:nframe+1,2),centpxl(1:nframe+1,1),'.');
  title(sprintf('Frame %04d, centroid= %1.2f x, %1.2f y',t1frame,centpxl(1,2),centpxl(1,1)));
end

% for f = ana_frame(2:end)
for fcount = 0:0 %2
for k = 1:num_images
        
  nframe = nframe + 1;
  %video{nframe} = I1oirg;
  t2frame = f;
  border_rows = max((roi_rows(1) - borderoffset),1):...
    min((roi_rows(end) + borderoffset),imgR);
  border_cols = max((roi_cols(1) - borderoffset),1):...
    min((roi_cols(end) + borderoffset),imgC);
  UProi_rows = borderoffset*upsampfact+[1:length(roi_rows)*upsampfact];
  UProi_cols = borderoffset*upsampfact+[1:length(roi_cols)*upsampfact];
  I1b = I1orig(border_rows,border_cols);
  I2orig = double(imread([path1 '.tif'],k,'Info',info));
  I2b = I2orig(border_rows,border_cols);
  
  if upsampfact ~= 1
    I1border = imresize(I1b,upsampfact);
    I2border = imresize(I2b,upsampfact);
  else
    I1border = I1b;
    I2border = I2b;
  end
  I1 = I1border(UProi_rows,UProi_cols);
  I2 = I2border(UProi_rows,UProi_cols);
  
  % MASK METHOD 3: Chan Vese
  [seg1,phi_seg1,prop1] = chanvese(I1, 'type',seg2, 'singlecomp',1,'verbose',db,...
    'lambda_1', lam1, 'lambda_2', lam2, 'mu', mu, 'iter', nCViter);
  segarea(nframe,1) = sum(seg1(:));
  segarea(nframe,2:end) = [mean(I1(seg1==1)),std(I1(seg1==1))];
  bkgndind = (seg1 == 0);
  forgndind = (seg1 == 1);
  I1(bkgndind) = I1(bkgndind)*fadefact;
  I1(forgndind) = I1(forgndind)*enhancefact;
  I1border = I1border.*fadefact;
  I1border(UProi_rows,UProi_cols) = I1;
  
  [seg2,phi_seg2,prop2] = chanvese(I2, 'type',seg1, 'singlecomp',1,'verbose',db,...
    'lambda_1', lam1, 'lambda_2', lam2, 'mu', mu, 'iter', nCViter);
  segarea(nframe+1,1) = sum(seg2(:));
  bkgndind = (seg2 == 0);
  forgndind = (seg2 == 1);
  I2(bkgndind) = I2(bkgndind)*fadefact;
  I2(forgndind) = I2(forgndind)*enhancefact;
  I2border = I2border.*fadefact;
  I2border(UProi_rows,UProi_cols) = I2;
  centpxl(nframe+1,:) = [border_rows(UProi_rows(floor(prop2.IntCentroid(1)))) + rem(prop2.IntCentroid(1),1),...
    border_cols(UProi_cols(floor(prop2.IntCentroid(2)))) + rem(prop2.IntCentroid(2),1)];
  
  if nframe > 2
    if segarea(nframe,2) > segarea(nframe-1,2) + segarea(nframe-1,3) ...
        || segarea(nframe,2) < segarea(nframe-1,2) - segarea(nframe-1,3)
      disp(sprintf('Detected segementation error at frame %d',f))
      return;
    end
  end
  
  if db;
    % Real time tracking output ---------------------------------
    set(0,'CurrentFigure',trackfig);
    hold off;
    imagesc(I2orig); colormap(gray); hold on;
    plot(centpxl(1:nframe+1,2),centpxl(1:nframe+1,1),'.');
    title(sprintf('Frame %04d, centroid= %1.2f x, %1.2f y',...
      t2frame,centpxl(nframe+1,2),centpxl(nframe+1,1)));
    
    if exist('fpanel');
      close(fpanel)
    end
    close(prop1.phih)
    close(prop2.phih)
    fpanel = figure('position',[129    67   801   459]);% 408   345   801   459]);
    subplot(2,3,1); imagesc(I1,[0,255]); colormap(gray); hold on;
    title(sprintf('I1blk, Frame %04d',t1frame));
    
    subplot(2,3,2); imagesc(I2,[0,255]); colormap(gray);
    title(sprintf('I2blk, Frame %04d',t2frame));
  elseif (mod(nframe,500) == 0)
    if ~exist('trackfig');
        trackfig = figure;
    end
    set(0,'CurrentFigure',trackfig);
    hold off;
    imagesc(I2orig); colormap(gray); hold on;
    plot(centpxl(1:nframe+1,2),centpxl(1:nframe+1,1),'.');
    title(sprintf('Frame %04d, centroid= %1.2f x, %1.2f y',...
      t2frame,centpxl(nframe+1,2),centpxl(nframe+1,1)));
    pause(1)
  end
  
  if (mod(nframe,100) == 0)
    disp(nframe)
  end
  
  roi_rows = max(1,roi_rows(prop2.BoundingBox(2)) - roi_border) + ...
    [0:(prop2.BoundingBox(4) + 2*roi_border)];
  roi_cols = max(1,roi_cols(prop2.BoundingBox(1)) - roi_border) + ...
    [0:(prop2.BoundingBox(3) + 2*roi_border)];
  seg2mask = seg2(prop2.BoundingBox(2)+[0:prop2.BoundingBox(4)],...
    prop2.BoundingBox(1)+[0:prop2.BoundingBox(3)]);
  s2 = zeros(length(roi_rows),length(roi_cols));%zeros(size(seg2mask,1)+2*roi_border,size(seg2mask,2)+2*roi_border);
  s2([roi_border+[1:size(seg2mask,1)]],[roi_border+[1:size(seg2mask,2)]]) = seg2mask;
  seg2 = s2;
  nCViter = 200;
  I1 = I2orig(roi_rows,roi_cols);
  I1orig = I2orig;
  t1frame = t2frame;
  
  if save_vid
    Itst = repmat(I2orig./256,[1,1,3]);
    for i = 1:nframe+1
      Itst(round(centpxl(i,1)),round(centpxl(i,2)),1) = 1;
    end
    imwrite(Itst,sprintf('%svideo%04d.tif',savepath,f),...
      'tif');
  end
f = f + 1;  
end
fcount = fcount + 1;
info = imfinfo([path1 num2str(fcount) '.tif']);
num_images = numel(info);
end
set(0,'CurrentFigure',trackfig);
hold off;
imagesc(I2orig); colormap(gray); hold on;
plot(centpxl(:,2),centpxl(:,1),'.');
title(sprintf('Frame %04d, centroid= %1.2f x, %1.2f y',...
  t2frame,centpxl(nframe+1,2),centpxl(nframe+1,1)));