function imoverlay(fixed, Fstr, moving, moving_r, Mstr,varargin)
%% Script: imoverlay(I1str, I1, Fnum, I2)
% Example: imoverlay('Fixed',fixed,'Moving',moving_reg)
% imoverlay(F,Fstr,M,Mr,Mstr,'displaymode','rgbgray','colormap','fire','trans',0.4);
% imoverlay(F,Fstr,M,Mr,Mstr,'displaymode','checker','colormap','jet','blocksz',6);
% imoverlay(F,Fstr,M,Mr,Mstr,'displaymode','rgbchan','1channel',[1,1,0]);
% imoverlay(F,Fstr,M,Mr,Mstr,'displaymode','edgeoverlay','colormap','gray','edgecolor','r','edgeth',20);
% 
% 
% Description: Based off of information provided at:
% http://blogs.mathworks.com/steve/2009/02/18/image-overlay-using-transparency/
% http://www.mathworks.com/support/solutions/en/data/1-1AK7N/
% Example:
% Required Functions:
% -------------------------------------------------------------------------
% INPUTS
% -------------------------------------------------------------------------
% fixed: fr x fc matrix which represents the fixed image
% Fstr: string used to identify the fixed image (in title with text
%   "\mathcal{F}_{Fstr}"
% moving: mr x mc matrix which represents the original moving image
% moving_r: M x N matrix which represents the moving image with applied
%   transformation
% Mstr: string used to identify the moving image (in title with text
%   "\mathcal{M}^r_{Mstr}"
% varargin - 'PropertyName','PropertyValue'
%   'nbins': number of colorbins that should be used [DEFAULT = 256]
%   'displaymode': string indicating display mode for overlayed images. 
%     [DEFAULT = 'rgbgray']
%     Valid strings for 'PropertyValue' include
%       • 'rgbgray': fixed image is pseudocolored into jet colormap and
%         moving_r is overlayed in grayscale (if fixed is fr x fc x 3, then
%         fixed will be plotted as an rgb colored image)
%         -> 'transparency': integer (0,1) [DEFAULT = 0.7]
%         -> 'colormap': fixed image colormap string ('gray','jet','fire')
%         [DEFAULT = 'jet']
%       • 'checker': image with 4x4 checkerboard pattern where "black" (top
%         left) blocks are the fixed image and "white" blocks are moving_r
%         -> 'blocksz': integer [DEFAULT = 4]
%         -> 'colormap': both image colormap string ('gray','jet') 
%         [DEFAULT = 'gray]
%       • 'rgbchan': creates a 3 channel image where R & B channels are the
%         fixed image and the G channel is moving_r
%         -> 'colormap': 1x3 "channel indicator" vector for fixed image
%         [DEFAULT = [1,0,1]]
%       • 'edgeoverlay': plots the fixed image using colormap indicated and
%         the moving_r image
%         -> 'colormap': fixed image colormap string ('gray','jet','fire)
%         [DEFAULT = 'jet']
%         -> 'edgeth': integer indicating the edge detector threshold for moving_r.
%         Suitable values are 20 for 255 image [DEFAULT = 0.1*image_range]
%         -> 'edgecolor': string indicating color of the edges
%         ('r','g','b','k','c','y','w') [DEFAULT = 'k']
%       • 'indiv': plots the two images individual. This variable overrides
%         the npanels variable
%       **NOTE: For all PropertyValues with 'colormap' as a property, if
%       the fixed image is 3D, then channels will be treated as RGB and
%       this will be the default
%   'npanels': integer coding for type of image display [DEFAULT = 4]
%     'PropertyValue' options
%       1: Single image which overlays the moving_r and fixed
%       2: plots same as 1 AND an image which overlays moving and fixed 
%       4: plots in a subplot(1,4,x) fashion (i) the fixed image alone,
%         (ii) the moving image alone, (iii) the transformed moving image
%         alone, (iv) the same as 1
%     
% -------------------------------------------------------------------------
% OUTPUTS
% -------------------------------------------------------------------------
%
% Date           Author            E-mail                      Version
% 10 July 2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     1
% 10  May 2013   Amanda Gaudreau   amanda.gaudreau@gmail.com     2

%% Image Overlay
set(0,'DefaultAxesFontSize',18,'defaultfigurecolor','w',...%'DefaultAxesFontWeight','bold',...
  'DefaultLineLineWidth',2, 'Defaultaxesposition','remove');%, 'defaultfontsize', 14);

% Determine user options
PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

if strmatch('displaymode',PropertyNames)
  dspmode = PropertyVal{strmatch('displaymode',PropertyNames)};
else
  dspmode = 'rgbgray';
end

if strmatch('npanels',PropertyNames)
  npanels = PropertyVal{strmatch('npanels',PropertyNames)};
else
  npanels = 4;
end

if strmatch('nbins',PropertyNames)
  ncol = PropertyVal{strmatch('nbins',PropertyNames)};
else
  ncol = 256;
end

if strmatch('colormap',PropertyNames)
  colmapstr = PropertyVal{strmatch('colormap',PropertyNames)};
else
  colmapstr = 'jet';
end

if strmatch('colorbar',PropertyNames)
  cbar = PropertyVal{strmatch('colorbar',PropertyNames)};
else
  cbar = 0;
end

if strmatch('rgbchan',dspmode)
  if strmatch('colormap',PropertyNames)
    colmapstr = PropertyVal{strmatch('colormap',PropertyNames)};
  else
    colmapstr = [1,0,1];
  end
end

if strmatch('figh',PropertyNames)
  h = PropertyVal{strmatch('figh',PropertyNames)};
else
  h = [];
end

fchorder = 1:2;

Cmap = custom_colormaps(colmapstr,ncol);

% Prepare image variables
fr = size(fixed,1); fc = size(fixed,2); 
mr = size(moving,1); mc = size(moving,2);
[M,N] = size(moving_r);
TxMrng = max(moving_r(:)) - min(moving_r(:));
TxMnorm = ( moving_r - min(moving_r(:)) )./TxMrng;
Mrng = max(moving(:)) - min(moving(:));
Mnorm = ( moving - min(moving(:)) )./Mrng;

if length(size(fixed)) > 2
  fixedCmap = [];
  if size(fixed,3) >= 3
    Frgb = fixed(:,:,1:3);
  else
    Frgb = zeros(fr,fc,3);
    Frgb(:,:,fchorder) = fixed;
  end
else
  fixedCmap = Cmap;
  Frng = max(fixed(:)) - min(fixed(:));
  Fnorm = (fixed - min(fixed(:)))./Frng;
  Frgb = zeros(M,N,3);
  mapind = round(Fnorm(:).*(ncol-1)+1);
  CmapF = Cmap(mapind,:);
  Frgb(:,:,1) = reshape(CmapF(:,1),M,N);
  Frgb(:,:,2) = reshape(CmapF(:,2),M,N);
  Frgb(:,:,3) = reshape(CmapF(:,3),M,N);
  colmapstr = 'jet';
end
Frgb = (Frgb - min(Frgb(:)))./(max(Frgb(:))- min(Frgb(:)));
switch dspmode
%--------------------------------------------------------------------------
  case 'rgbgray' % Plots two separate images: Frgb and moving in gray
    if strmatch('trans',PropertyNames)
      trans_factor = PropertyVal{strmatch('trans',PropertyNames)};
    else
      trans_factor = 0.7;
    end
    movingCmap = custom_colormaps(colmapstr,ncol);%
%-------------------------------------------------------------------------- 
  case 'checker'
    if strmatch('blocksz',PropertyNames)
      blocksz = PropertyVal{strmatch('blocksz',PropertyNames)};
    else
      blocksz = 4;
    end
    movingCmap = custom_colormaps('gray',ncol);%'gray';
    mask = (checkerboard(blocksz,ceil(mr/(2*blocksz)),ceil(mc/(2*blocksz))) > 0.5);
    mask = repmat(mask(1:mr,1:mc),[1,1,3]);
    Mrgb = zeros(M,N,3);
    if length(size(TxMnorm)) > 2
      Mrgb = TxMnorm;
    else
      mapind = round(TxMnorm(:).*(ncol-1)+1);
      CmapM = movingCmap(mapind,:);
      Mrgb(:,:,1) = reshape(CmapM(:,1),M,N);
      Mrgb(:,:,2) = reshape(CmapM(:,2),M,N);
      Mrgb(:,:,3) = reshape(CmapM(:,3),M,N);
    end
    Icomb = Frgb;
    Icomb(mask == 1) = Mrgb(mask == 1);
    
    
    if npanels == 2
      Mrgb = zeros(M,N,3);
      mapind = round(Mnorm(:).*(ncol-1)+1);
      CmapM = Cmap(mapind,:);
      Mrgb(:,:,1) = reshape(CmapM(:,1),M,N);
      Mrgb(:,:,2) = reshape(CmapM(:,2),M,N);
      Mrgb(:,:,3) = reshape(CmapM(:,3),M,N);
      Icomb0 = Frgb;
      Icomb0(mask == 1) = Mrgb(mask == 1);
    end
    
%--------------------------------------------------------------------------
  case 'rgbchan'
    Mind = abs(ones(1,3) - colmapstr);
    Mindnon0 = find(Mind == 1);
    Icomb = Frgb;
    for i = Mindnon0(:)'
      Icomb(:,:,i) = TxMnorm;
    end
    movingCmap = custom_colormaps(Mind,ncol);%'gray';
    
    if npanels == 2
      Icomb0 = Frgb;
      for i = Mindnon0(:)'
        Icomb0(:,:,i) = Mnorm;
      end
    end
    
%--------------------------------------------------------------------------  
  case 'edgeoverlay'
    if strmatch('edgeth',PropertyNames)
      edgeth = PropertyVal{strmatch('edgeth',PropertyNames)};
    else
      edgeth = 0.1*TxMrng;
    end
    
    if strmatch('edgecolor',PropertyNames)
      edgecol = PropertyVal{strmatch('edgecolor',PropertyNames)};
    else
      edgecol = 'k';
    end
    
    switch edgecol
      case 'r'; edgebin = [1,0,0];
      case 'g'; edgebin = [0,1,0];
      case 'b'; edgebin = [0,0,1];
      case 'k'; edgebin = [0,0,0];
      case 'w'; edgebin = [1,1,1];
      case 'c'; edgebin = [0,1,1];
    end
    movingCmap = 'gray';%custom_colormaps(colmapstr,ncol);%
    Mtedge = edge(moving_r,edgeth);
    Icomb = Frgb;
    ch1 = Icomb(:,:,1); ch2 = Icomb(:,:,2); ch3 = Icomb(:,:,3);
    ch1(find(Mtedge == 1)) = edgebin(1); ch2(find(Mtedge == 1)) = edgebin(2); ch3(find(Mtedge == 1)) = edgebin(3);
    Icomb(:,:,1) = ch1; Icomb(:,:,2) = ch2; Icomb(:,:,3) = ch3;
    
    if npanels == 2
      Mtedge = edge(moving,edgeth);
      Icomb0 = Frgb;
      ch1 = Icomb0(:,:,1); ch2 = Icomb0(:,:,2); ch3 = Icomb0(:,:,3);
      ch1(find(Mtedge == 1)) = edgebin(1); ch2(find(Mtedge == 1)) = edgebin(2); ch3(find(Mtedge == 1)) = edgebin(3);
      Icomb0(:,:,1) = ch1; Icomb0(:,:,2) = ch2; Icomb0(:,:,3) = ch3;
    end
%--------------------------------------------------------------------------
  case 'indiv'
    
end

if isempty(h)
  figure;
end

switch npanels
  case 1
    if ~exist('Icomb')
      im2 = image(Frgb); set(im2,'cdatamapping','scaled');
      hold on;
      if size(TxMnorm,3) > 1
        im1 = image(TxMnorm);
        set(im1,'AlphaData',ones(size(im1,1),size(im1,2),size(im1,3)).*trans_factor)
      else
        im1 = image(TxMnorm(:,:,1)); set(im1,'Cdatamapping','scaled'); colormap(movingCmap);
        set(im1,'AlphaData',ones(size(im1,1),size(im1,2)).*trans_factor)
      end
      %set(gca,'xticklabel','','yticklabel','');
      hold off;
    else
      im3 = image(Icomb); set(im3,'cdatamapping','scaled');
    end
    title(sprintf('$$\\mathcal{M}^r_{%s}, \\mathcal{F}_{%s}$$ Overlayed',...
      Mstr,Fstr),'interpreter','latex')
    set(gcf,'position',[50   338   731   645])
    if cbar
      colorbar(im1);
    end
  case 2
    subplot(1,2,1)
    if ~exist('Icomb')
      im0 = image(Frgb); set(im0,'cdatamapping','scaled');
      hold on;
      im1 = image(moving); set(im1,'Cdatamapping','scaled'); colormap(movingCmap);
      %set(gca,'xticklabel','','yticklabel','');
      set(im1,'AlphaData',ones(size(im1,1),size(im1,2)).*trans_factor)
      hold off;
    else
      im3 = image(Icomb0); set(im3,'cdatamapping','scaled');
    end
    title(sprintf('$$\\mathcal{M}_{%s}, \\mathcal{F}_{%s}$$ Overlayed',...
      Mstr,Fstr),'interpreter','latex')
    
    subplot(1,2,2)
    if ~exist('Icomb')
      im2 = image(Frgb); set(im2,'cdatamapping','scaled');
      hold on;
      im1 = image(moving_r); set(im1,'Cdatamapping','scaled'); colormap(movingCmap);
      %set(gca,'xticklabel','','yticklabel','');
      set(im1,'AlphaData',ones(size(im1,1),size(im1,2)).*trans_factor)
      hold off;
    else
      im3 = image(Icomb); set(im3,'cdatamapping','scaled');
    end
    title(sprintf('$$\\mathcal{M}^r_{%s}, \\mathcal{F}_{%s}$$ Overlayed',...
      Mstr,Fstr),'interpreter','latex')
    set(gcf,'position',[50   338   731   645])
  case 4
    subplot(1,4,1);
    im0 = imagesc(moving); set(im0,'Cdatamapping','scaled'); colormap(movingCmap);
    title(sprintf('$$\\mathcal{M}_{%s}(\\vec{q})$$',Mstr),'interpreter','latex')
  
    subplot(1,4,2);
    im1 = imagesc(moving_r); set(im1,'Cdatamapping','scaled'); colormap(movingCmap);
    title(sprintf('$$\\mathcal{M}_{%s}^r(\\hat{\\vec{p}})$$',...
      Mstr),'interpreter','latex');
    
    subplot(1,4,3);
    im2 = image(Frgb); set(im2,'cdatamapping','scaled');
    %set(gca,'xticklabel','','yticklabel','');
    title(sprintf('$$\\mathcal{F}_{%s}(\\vec{p})$$',Fstr),'interpreter','latex')
    
    subplot(1,4,4);
    if ~exist('Icomb')
      im2 = image(Frgb); set(im2,'cdatamapping','scaled');
      hold on;
      im1 = image(moving_r); set(im1,'Cdatamapping','scaled'); colormap(movingCmap);
      %set(gca,'xticklabel','','yticklabel','');
      set(im1,'AlphaData',ones(size(im1,1),size(im1,2)).*trans_factor)
      hold off;
    else
      im3 = image(Icomb); set(im3,'cdatamapping','scaled');
    end
    title(sprintf('$$\\mathcal{M}^r_{%s}, \\mathcal{F}_{%s}$$ Overlayed',...
        Mstr,Fstr),'interpreter','latex')
    set(gcf,'position',[50   338   731   645])
end
end













