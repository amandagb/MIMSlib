function imoverlay1x4(fixed, Fstr, moving, moving_r, Mstr)
%% Script: imoverlay(I1str, I1, Fnum, I2)
% Example: imoverlay('Fixed',fixed,'Moving',moving_reg)
% Description: Based off of information provided at:
% http://blogs.mathworks.com/steve/2009/02/18/image-overlay-using-transparency/
% http://www.mathworks.com/support/solutions/en/data/1-1AK7N/
% Example:
% Required Functions:
% INPUTS  ----------------------------------------------------------------
%
% OUTPUTS ---------------------------------------------------------------
%
%  Date           Author            E-mail                      Version
%  10 July 2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     1

%% Image Overlay
% Create the background This example uses a blend of colors from left to
% right, converted to a TrueColor image Use repmat to replicate the pattern
% in the matrix Use the "jet" colormap to specify the color space
set(0,'DefaultAxesFontSize',18,'defaultfigurecolor','w',...%'DefaultAxesFontWeight','bold',...
  'DefaultLineLineWidth',2, 'Defaultaxesposition','remove');

trans_factor = 0.7;
if ~exist('colind'); colind = 'gray'; end;
fntsz = 14;
[M,N] = size(moving_r);
jetmap = jet(256);
moving = round( ( ( moving - min(moving(:)) )./( max(moving(:)) - min(moving(:)) ) ).*255 )./255;
moving_r = round( ( ( moving_r - min(moving_r(:)) )./( max(moving_r(:)) - min(moving_r(:)) ) ).*255 )./255;
I1rgb = zeros(M,N,3);
% Nice configurations include I1 in R & B (fuscia) channels OR in B & G
% (teal) OR R & G (yellow)
I1rgb(:,:,1) = moving_r;%I1;
I1rgb(:,:,2) = 0;
I1rgb(:,:,3) = moving_r;

[M,N] = size(fixed);
fixed = round(((fixed - min(fixed(:)))./(max(fixed(:)) - min(fixed(:)))).*255)./255;
I2rgb = zeros(M,N,3);
jetI2 = jetmap(fixed(:).*255+1,:);
I2rgb(:,:,1) = reshape(jetI2(:,1),M,N);
I2rgb(:,:,2) = reshape(jetI2(:,2),M,N);
I2rgb(:,:,3) = reshape(jetI2(:,3),M,N);


figure;
subplot(1,4,1);
im0 = imagesc(moving,[0,1]); set(im0,'Cdatamapping','scaled');
colormap(gray);
%set(gca,'xticklabel','','yticklabel','');
title(sprintf('$$\\mathcal{M}_{%s}(\\vec{q})$$',Mstr),'interpreter','latex')

subplot(1,4,2);
im1 = imagesc(moving_r,[0,1]); set(im1,'Cdatamapping','scaled');
colormap(gray);
%set(gca,'xticklabel','','yticklabel','');
title(sprintf('$$\\mathcal{M}_{%s}^r(\\hat{\\vec{p}})$$',...
  Mstr),'interpreter','latex')

subplot(1,4,3);
im2 = image(I2rgb); set(im2,'cdatamapping','scaled');
%set(gca,'xticklabel','','yticklabel','');
title(sprintf('$$\\mathcal{F}_{%s}(\\vec{p})$$',Fstr),'interpreter','latex')

subplot(1,4,4);
im2 = image(I2rgb);
set(im2,'cdatamapping','scaled');
hold on;
im1 = image(moving_r); set(im1,'Cdatamapping','scaled');
colormap(gray);
%set(gca,'xticklabel','','yticklabel','');
title(sprintf('$$\\mathcal{M}^r_{%s}, \\mathcal{F}_{%s}$$ Overlayed',...
  Mstr,Fstr),'interpreter','latex')
set(im1,'AlphaData',ones(size(im1,1),size(im1,2)).*trans_factor)
hold off;

set(gcf,'position',[  50         688        1408         290])
end













