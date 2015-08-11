function [phi0,mask] = initlevelset(I,varargin)
%% Script: m = initlevelset(I,type)
% m = initCVmask(I,'type','radius',[#,#],'density',[#,#]);
% Description: Initializes a mask for the Chan-Vese algorithm
% Example:
% Required Functions:
% INPUTS ----------------------------------------------------------------
% I: matrix containing image or image dimensions ([# rows, # cols])
% varargin - 'PropertyName','PropertyValue'
%   'type': string indicating type of mask
%     'small': small oval in the center of the page (radius = 10% of
%     pixel width and height)
%     'medium': medium oval in the center of the page (radius = 25% of
%     pixel width and height)
%     'large': large oval in the center of the page (radius = 50% of
%     pixel width and height)
%     'diamonds': 
%     'circles': default
%     M x N matirx containing mask -- output will be levelset
%   'radius': one or two element vector containing the integer(s)
%   indicating radius of the circle(s) in the mask
%   'density': [1 x 2] matrix of the density of mask elements one would
%   like in the rows (1) and columns (2) -- used in 'circles'
% OUTPUTS ---------------------------------------------------------------
%
%
%  Date           Author            E-mail                      Version
%  02 May 2012    Amanda Gaudreau   amanda.gaudreau@gmail.com     1

PropertyNames = lower(varargin(1:2:length(varargin)));
PropertyVal = varargin(2:2:length(varargin));

if strmatch('type',PropertyNames)
  type = PropertyVal{strmatch('type',PropertyNames)};
else
  type = 'circles';
end

%% Function based off of Yue Wu's 2009 Tuft's CV algorithm (maskcircle2.m)
if ischar(type)
  if sum(size(I)) == 2
    [M,N] = I;
  elseif size(I,3)~=3
    temp = double(I(:,:,1));
    [M,N] = size(temp);
  else
    temp = double(rgb2gray(I));
    [M,N] = size(temp);
  end
  
  % Algorithm Specifictions
  if strmatch('radius',PropertyNames)
    r = PropertyVal{strmatch('radius',PropertyNames)};
    if length(size(r)) == 2
      rx = r; ry = r;
    elseif length(size(r)) > 2
      rx = r(1); ry = r(2);
    end
  else
    switch lower (type)
      case 'small' % Ellipse with radius == 10% of picture H and W
        rx = M*0.1;
        ry = N*0.1;
      case 'medium' % Ellipse with radius == 25% of picture H and W
        rx = M*0.25;
        ry = N*0.25;
      case 'large' % Ellipse with radius == 50% of picture H and W
        rx = M*0.5;
        ry = N*0.5;
      case 'diamonds'
        rx = 20;
      case 'circles'
        rx = 20;
    end
  end
  
  if strmatch('density',PropertyNames)
    d = PropertyVal{strmatch('density',PropertyNames)};
    dx = d(1); dy = d(2);
  else
    dx = 5; dy = 5;
  end
  
  cx = M/2; %centers phi0 at center of x axis
  cy = N/2; %centers phi0 at center of y axis
  x = repmat([1:M]',1,N); %creates a matrix that is szx by szy has all columns == row number i => aij = i
  y = repmat([1:N],M,1); %creates a matrix that is szx by szy has all columns == column number j => aij = j
  
  switch lower (type)
    case 'diamonds'
      mask = zeros(round(ceil(max(M,N)/(2*(rx+1)))*3*(rx+1)));
      siz = size(mask,1);
      sx = round(siz/2);
      i = 1:round(siz/(2*(rx+1)));
      j = 1:round(0.9*(siz/(2*(rx+1))));
      j = j-round(median(j));
      mask(sx+2*j*(rx+1),(2*i-1)*(rx+1)) = 1;
      se = strel('disk',rx);
      mask = imdilate(mask,se);
      mask = mask(round(siz/2-M/2-6):round(siz/2-M/2-6)+M-1,round(siz/2-N/2-6):round(siz/2-N/2-6)+N-1);
    case 'circles'
      phi0 = Inf * ones(M,N);
      mask = zeros(M,N);
      for ch = 1:dx
        for cv = 1:dy
          h = 1 * (sqrt((y - (ch*N / (dx+1))).^2 + (x - (cv*M / (dy+1))).^2) - rx);
          phi0 = min(phi0, h);
        end
      end
  end
  
  if ~exist('mask','var')
    mask = zeros(size(temp));
    %m((x-cx).^2+(y-cy).^2<min(rx,ry).^2) = 1; %Gives circular boundary
    mask(ry^2.*(x-cx).^2+rx^2.*(y-cy).^2<(rx^2*ry^2)) = 1; %Gives elliptical boundary
    tem(:,:,1) = mask;
    tem(:,:,2) = mask;
    mask = tem;
  else
    tem(:,:,1) = mask;
    M = padarray(mask,[floor(2/3*rx),floor(2/3*rx)],0,'post');
    tem(:,:,2) = M(floor(2/3*rx)+1:end,floor(2/3*rx)+1:end);
    mask = tem;
  end
else
  mask = type;
end

if ~exist('phi0');
  phi0 = bwdist(mask(:,:,1))-bwdist(1-mask(:,:,1))+im2double(mask(:,:,1))-.5;
end;
mask = zeros(size(phi0));
mask(phi0 >= 0) = 1;
end