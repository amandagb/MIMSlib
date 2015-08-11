function h = overlayimgs(data,elem_ana_range)
%% Script:
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
%  25 Oct 2011    Amanda Gaudreau   amanda.gaudreau@gmail.com     1
%% Image Overlay
% Create the background This example uses a blend of colors from left to
% right, converted to a TrueColor image Use repmat to replicate the pattern
% in the matrix Use the "jet" colormap to specify the color space
fntsz = 14;
scale = 4;
[elem_ana_range,fldnm] = interp_data_elemrng(data,elem_ana_range);
d = zeros(size(getfield(data,fldnm{elem_ana_range(1)}),1),...
  size(getfield(data,fldnm{elem_ana_range(1)}),2),length(elem_ana_range));
i = 0;
for f = elem_ana_range
  i=i+1;
  if isstruct(data)
    d(:,:,i) = getfield(data,fldnm{f});
  else
    d(:,:,i) = data(:,:,i);
  end
end
d = d([1:52*scale],:,:);

d1 = (round((d(:,:,1)./max(max(d(:,:,1)))).*255)./255);
d1(d1 < 0) = 0;
d1rgb = zeros(size(d1,1),size(d1,2),3);
d1rgb(:,:,1) = d1;
d1rgb(:,:,3) = d1;%0;%
d1rgb(:,:,2) = d1;%0;%
subplot(i,3,1);
image(d1rgb);
title(sprintf('%s Original',fldnm{elem_ana_range(1)}),'fontsize',fntsz)
set(gca,'xticklabel','','yticklabel','');

% Create an image and its corresponding transparency data This example uses
% a random set of pixels to create a TrueColor image
for j = 2:i
  subplot(i,3,1+3*(j-1))
  im = d(:,:,j);
  im(im<0) = 0;
  %d1rgb(:,:,3) = (round((im./max(max(im))).*255)./255);
  %im = d1rgb; im(:,:,1) = 0;
  k(i-1) = image(im);
  set(k(i-1),'Cdatamapping','scaled');
  set(gca,'xticklabel','','yticklabel','');
  title(sprintf('%s Original',fldnm{elem_ana_range(j)}),'fontsize',fntsz)
end

subplot(i,3,[2,3,5,6]);
image(d1rgb); hold on; set(gca,'xticklabel','','yticklabel','');
for j = 2:i
  h(j-1) = image(im);
  set(h(j-1),'AlphaData',ones(size(im,1),size(im,2)).*0.3, 'cdatamapping','scaled');
end
title(sprintf('Overlayed images: %s gray, %s colored',...
  fldnm{elem_ana_range(1)},fldnm{elem_ana_range(2)}),'fontsize',fntsz);

end