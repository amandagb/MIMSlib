function grayI_CmapOver(grayI,colIover)
%% Script: grayI_CmapOver
% Description:  The script (which can easily be made into a function) takes
% in two single channeled images => grayI and colIover The grayI is the
% image that will be displayed as the "bottom" grayscale image. The
% colIover image is the image that you would like to display using a
% colormap. You must also ensure that the following variables are defined
% appropriately for your application:
%     minColVal, maxColVal, figname, cmapstr, tranfact
% Example:
% Required Functions: auto_thresh, export_fig, custom_colormaps
% -------------------------------------------------------------------------
% INPUTS: Define in workspace
% -------------------------------------------------------------------------
%   grayI:    M x N x 1 image that will serve as the base grayscale image
%   colIover: K x L x 1 image that will serve as the overlayed image in the
%       indicated colormap
%   
% -------------------------------------------------------------------------
% OUTPUTS
% -------------------------------------------------------------------------
%
% Date           Author                   E-mail                      Version
% 26 Sept 2016   A Gaudreau-Balderrama    amandagbalderrama@gmail.com     0

%--- IMPORTANT: make sure these variables are defined before running the
%script!!
if ~exist('grayI')
  errordlg('Variable grayI is not defined in the workspace')
end

if ~exist('colIover')
  errordlg('Variable colIover is not defined in the workspace')
end

flipmap = 0;
cmapstr = 'parula';  % see help custom_colormaps for colormap choices 
figname = 'imageOverlay';   % sting to use for naming the figure
tranfact = 0.7;  % larger gives a less transparent colormap overlayed image
savefig = 0;
if savefig
  close all;
end
minColVal = min(colIover(:)) ; maxColVal = max(colIover(:));  % maximum and minimum values for colorbar

% Condition gray image using autothresh function. This step is not
% necessary unless dynamic range of the data need to be adjusted (as in
% MIMS or a low contrast image for instance)
[~,~,grayI] = auto_thresh(grayI,'auto',[]);

% Condition image which will be the colormap overlay image
I2 = colIover;
I2(I2 == 0) = NaN;
Il2up = imresize(colIover,size(grayI));
Il2th = Il2up;
Il2th(Il2th < minColVal) = NaN;
Il2th(Il2th > maxColVal) = maxColVal;

nonNANi = find(isnan(Il2th) == 0);
%     [~,~,Il2th] = auto_thresh(Il2th,'auto',[]);

%---- Create RGB image of grayscale image => grayI3chv
grayCmap = repmat([0:255]'./255,1,3);
[M,N] = size(grayI);
phrng = max(grayI(:)) - min(grayI(:));
phnorm = (grayI - min(grayI(:)))./phrng;
phnorm(isnan(phnorm)) = 0;
grayI3ch = zeros(M,N,3);
mapind = round(phnorm(:).*(size(grayCmap,1)-1)+1);
gray3chv = grayCmap(mapind,:);
grayI3ch(:,:,1) = reshape(gray3chv(:,1),M,N);
grayI3ch(:,:,2) = reshape(gray3chv(:,2),M,N);
grayI3ch(:,:,3) = reshape(gray3chv(:,3),M,N);
%imwrite(P3ch,strcat(cdir.mtngpath,pad_IVISinfo.ivisfldr.ivisSUfldr{i},'img.png'),'png')

% Plot grayscale image
% hPh = figure('position',[1,1,round(M/3),round(N/3)]);
% image(grayI3ch); axis image;
% set(gca,'position',[0,0,1,1]);
% if savefig
%   export_fig(strcat(cdir.mtngpath,pad_IVISinfo.ivisfldr.ivisSUfldr{i},'_photo'),...
%     '-png');
%   close(hPh);
% end

%---- Create RGB image of colormap overlay image with indicated colormap =>
% L3chv
if flipmap
  lCmap = flipud(custom_colormaps(cmapstr));
else
  lCmap = custom_colormaps(cmapstr);
end
[M,N] = size(Il2th);
lrng = max(Il2th(:)) - min(Il2th(:));
lnorm = (Il2th - min(Il2th(:)))./lrng;
L3chv = zeros(M*N,3);
L3ch = nan(M,N,3);
mapind = round(lnorm(:).*(size(lCmap,1)-1)+1);
L3chv(nonNANi,:) = lCmap(mapind(nonNANi),:);
L3ch(:,:,1) = reshape(L3chv(:,1),M,N);
L3ch(:,:,2) = reshape(L3chv(:,2),M,N);
L3ch(:,:,3) = reshape(L3chv(:,3),M,N);

% Plot colormap overlay image
% hL = figure('position',[1,1,round(M/3),round(N/3)]);
% imagesc(Il2th,[minColVal,maxColVal]); %(pad_Ilum2{i})
% axis image; colormap(lCmap); cbh = colorbar;
% set(gca,'xticklabel','','yticklabel','','position',[0,0,1,1]);
% set(cbh,'position',[0.8,0.05,0.05,0.9],'YColor','w','ylim',[minColVal,maxColVal])
% if savefig
%   export_fig(strcat(cdir.mtngpath,pad_IVISinfo.ivisfldr.ivisSUfldr{i},'_umx2_',cmapstr),...
%     '-png');
%   close(hL);
% end

%--- Plot overlay of two images. im_ghost is the unlying image which
%controls the colorbar properties
hall = figure;%('position',[1,1,round(M/3),round(N/3)]);
im_ghost = imagesc(Il2th,[minColVal,maxColVal]);
colormap(lCmap); hold on;
axh = gca;
im2 = image(L3ch,'cdatamapping','scaled');
im1 = image(grayI3ch,'Cdatamapping','scaled');
set(im1,'AlphaData',ones(size(im1,1),size(im1,2)).*tranfact)
cbh = colorbar('peer',axh,'yLim',[minColVal,maxColVal]);
set(gcf,'name',figname);
axis image;
set(gca,'xticklabel','','yticklabel','','position',[0,0,1,1]);
set(cbh,'position',[0.8,0.05,0.05,0.9],'YColor','w','ylim',[minColVal,maxColVal])
xv = xlim; yv = ylim;

%text(yv(1),xv(1),{[sprintf('IVIS %s',pad_IVISinfo.animalseval{i})];...
%  [injstr];...
%  [pad_IVISinfo.ivisXL{hdroffset+i,XLhcols.injectioncol}]},...
%  'fontsize',14,'HorizontalAlignment','left',...
%  'VerticalAlignment','top','Color','w','FontName','ariel');
%if savefig
%  export_fig(strcat(cdir.mtngpath,pad_IVISinfo.ivisfldr.ivisSUfldr{i},'_umx2overlay_',cmapstr,xtrastr),...
%    '-png');
%  close(hall);
%end

% dockf on all;
% save_open_figures(cdir.mtngpath,[],[],'umx2overlay');
end