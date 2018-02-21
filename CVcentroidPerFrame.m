function detectedInfo = CVcentroidPerFrame(grayImage,roipos,pastMask);
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
%  28 June 2016   Amanda Balderrama amandagbalderrama@gmail.com   0
%   Based on v3.1 of CVcentroid_track function
%==========================================================================

if size(grayImage,3) == 3
  grayImage = rgb2gray(grayImage);
elseif size(grayImage,3) > 1
  grayImage = grayImage(:,:,1);
end
[imgR,imgC] = size(grayImage);

if ~exist('pastMask')
  initMask = 'small';
else
  initMask = pastMask;
end

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
t1frame = 1; %ana_frame(1);
f = 1;
centpxl = zeros(1,2); % gives the row and column (which is the y and x coordinate)
segarea = zeros(1,3); % [1: total area of mask, 2: mean of pixels in mask, 3: standard dev of pixels in mask]

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

roipos(3:4) = roipos(3:4) + 1.*rem(roipos(3:4),2); % <= makes even; 0.^(rem(roipos(3:4),2)); <= makes odd
roi_rows = round(roipos(2)):round((roipos(2)+roipos(4)));
roi_cols = round(roipos(1)):round((roipos(1)+roipos(3)));

% define image block and segment ---------------------------------
I1roi = grayImage(roi_rows,roi_cols);
[seg1,phi_seg1,prop1] = chanvese(I1roi, 'type',initMask, 'singlecomp',1,'verbose',0,...
  'lambda_1', lam1, 'lambda_2', lam2, 'mu', mu, 'iter', 1000);
bkgndind = (seg1 == 0);
forgndind = (seg1 == 1);

detectedInfo.foregroundMask = zeros(imgR,imgC);
detectedInfo.foregroundMask(roi_rows,roi_cols) = seg1;
detectedInfo.levelSet = ones(imgR,imgC).*max(phi_seg1(:));
detectedInfo.levelSet(roi_rows,roi_cols) = phi_seg1;
detectedInfo.maskArea = sum(seg1(:));
detectedInfo.maskMeanIntensity = mean(I1roi(seg1 == 1));
detectedInfo.maskSTDIntensity = std(I1roi(seg1 == 1));
detectedInfo.weightedCentroid = [roi_cols(floor(prop1.IntCentroid(2)))+rem(prop1.IntCentroid(2),1),...
  roi_rows(floor(prop1.IntCentroid(1)))+rem(prop1.IntCentroid(1),1)]; %roi_cols(round(prop1.IntCentroid(2)))];

% % redefine roi rows and cols based on center of mass --------------------
% roi_rows = max(1,roi_rows(bndbx1(2)) - roi_border) + ...
%   [0:(bndbx1(4) + 2*roi_border)];
% roi_cols = max(1,roi_cols(bndbx1(1)) - roi_border) + ...
%   [0:(bndbx1(3) + 2*roi_border)];
% borderoffset = round(max(roipos(3:4))/2) + maxf2fmove;
% 
% if db
%   trackfig = figure;
%   imagesc(grayImage); colormap(gray); hold on;
%   plot(centpxl(1:nframe+1,2),centpxl(1:nframe+1,1),'.');
%   title(sprintf('Frame %04d, centroid= %1.2f x, %1.2f y',t1frame,centpxl(1,2),centpxl(1,1)));
% end

end