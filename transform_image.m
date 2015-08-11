function [Iuv,A,Mmatlab,iMmatlab] = transform_image(Ixy,varargin)
%% Script: new_image = transform_image(moving_image, alpha,extrapval,varargin)
% transform_image(moving_image, [x_shift, y_shift, rotation (radians), x_scale, y_scale, skew],'transtype','physicalORtrans',...
%   'extrapval',#,'fixedstr',{'tx'},'fixedvals',[#],...
%   'varstr',{'theta'},'varval',[5*pi/180])
% Description: rotation is applied about the center of the
% image. The parameters can be arbitrary floating point values. Default
% parameter values are as follows: [tx=0; ty=0; theta=0; sx=1; sy=1; k=0]
% Example:
%   im = imread('cameraman.tif');
%   imt = transform_image(im,[-20 -20 pi/8 1 1 0]);
%   clf; subplot(2,1,1); imagesc(im); title('Original');
%   subplot(2,1,2); imagesc(imt); title('Transformed');
%   colormap gray;
% Required Functions:
% INPUTS ----------------------------------------------------------------
% Ixy: Image on which affine transform will be performed according
%         to parameters of alpha. Only works on grayscale images
% alpha:  1x6 vector with the transformation parameters OR
%         3x3 matrix transformation matrix
%         [x_shift, y_shift, rotation (radians), x_scale, y_scale, skew]
%         NOTE: Positive displacement is to the left, up, and counterclockwise,
%           respectively.
% varargin - 'PropertyName','PropertyValue'
%   • 'extrapval': uses EXTRAPVAL to fill in any missing data. [DEFAULT = nan]
%   • 'interpmethod': string indiciating the interpolation method (see
%   interp2 function for possible inputs) [DEFAULT = 'bicubic']
%   • 'transtype': string indicating whether the user is inputting physical
%   parameters (alpha, n=6) -> 'physical'; or transform parameters (a,n=6)
%   • 'trans' [DEFAULT = 'physical']
%   • 'fixedstr': 1 x (6-n) cell with strings indicating the fixed parameters
%   • 'fixedvals': 1 x (6-n) vector of doubles indicate the value of the
%   parameter specified in 'fixedstr' value in the same order
%   • 'varstr': 1 x n cell with strings indicating the parameters that will
%   be varied
%   • 'varvals': 1 x n vector of doubles indicate the value of the parameter
%   specified in 'varstr' value in the same order
%   • 'rotpnt': string which specifies the point of rotation, either
%   'tl' (top-left) or 'center' [DEFAULT = 'center']
%   • 'shape': string or matrix indicating the desired size of the output
%   image [DEFAULT = 'same']
%     – 'same': Output image has the same dimensions as the input image
%     - 'valid': Output image is the space required to represent the
%     transformed image entirely
%     - [ROWS,COLS]: integers specifying the
%
% OUTPUTS ---------------------------------------------------------------
% Iuv
% A
% Mmatlab
% iMmatlab
%
%  Date           Author            E-mail                      Version
%  25 Apr 2011    Amanda Gaudreau   amanda.gaudreau@gmail.com     1
%  25 Apr 2011    Amanda Gaudreau   amanda.gaudreau@gmail.com     2
%  26 Apr 2011    Amanda Gaudreau   amanda.gaudreau@gmail.com     3
%  04 Apr 2013    Amanda Gaudreau   amanda.gaudreau@gmail.com     4
%  21 Aug 2013    Amanda Gaudreau   amanda.gaudreau@gmail.com     5
%   Adjusted for multidimensional input images (v5.1 saved 3 Oct 2013)
%  03 Oct 2013    Amanda Gaudreau   amanda.gaudreau@gmail.com     6
%   Fixing output image size option in varagin and resolving
%   inconsistencies between this function and imresize. Final solution
%   involved first setting a 0.5 offset to the x/y offset variables and
%   having the input space map to [0.5, OUT + 0.5]. A final hack was the
%   addition of 1/scale in dimensions where it was mapped from odd to odd
%   or even to even.
%  14 Oct 2013    Amanda Gaudreau   amanda.gaudreau@gmail.com     6.1

if ~isfloat(Ixy), Ixy = double(Ixy); end;
if ~isstr(varargin{1})
  alpha = varargin{1};
  start_argin = 2;
else
  alpha = [0, 0, 0, 1, 1,0];
  start_argin = 1;
end

PropertyNames = varargin(start_argin:2:length(varargin));
PropertyVal = varargin((start_argin + 1):2:length(varargin));

if strmatch('extrapval',PropertyNames)
  extrapval = PropertyVal{strmatch('extrapval',PropertyNames)};
else
  extrapval = nan;
end

if strmatch('interp',PropertyNames)
  interpstr = PropertyVal{strmatch('interp',PropertyNames)};
else
  interpstr = 'cubic';
end

if strmatch('rotpnt',PropertyNames)
  rotpnt = PropertyVal{strmatch('rotpnt',PropertyNames)};
else
  rotpnt = 'center';
end

[M, N, d] = size(Ixy);

if strmatch('tl',rotpnt)
  x_center = 0.5;
  y_center = 0.5;
else
  x_center = (N+1) / 2; % Changed to be xsize+1 since this give true center
  y_center = (M+1) / 2;
end

if strmatch('transtype',PropertyNames)
  ttype = PropertyVal{strmatch('transtype',PropertyNames)};
else
  ttype = 'physical';
end

if strmatch('shape',PropertyNames)
  outimsz = PropertyVal{strmatch('shape',PropertyNames)};
else
  outimsz = 'same';
end

% ------------------------------------------------------------------------
% Construct transformation matrix
% ------------------------------------------------------------------------
if sum(size(alpha) == [3,3]) == 2 || ~isempty(strmatch('trans',ttype))
  if size(alpha) == [3,3]
    if sum(alpha(3,:) == [0,0,1]) == 1  % input matrix is of the form M
      Mmatlab = alpha;
      A = Mmatlab';
      
      A(1,3) = Mmatlab(3) + Mmatlab(1)*x_center + Mmatlab(2)*y_center - x_center;
      A(2,3) = Mmatlab(6) + Mmatlab(4)*x_center + Mmatlab(5)*y_center - y_center;
      iMmatlab = inv(Mmatlab);
      iMmatlab(7:9) = [0;0;1];
    else
      A = alpha;
      Mmatlab = A';
      Mmatlab(3) = Mmatlab(3) - Mmatlab(1)*x_center - Mmatlab(2)*y_center + x_center;
      Mmatlab(6) = Mmatlab(6) - Mmatlab(4)*x_center - Mmatlab(5)*y_center + y_center;
      
      iMmatlab = inv(Mmatlab);
      iMmatlab(7:9) = [0;0;1];
    end
  else % This case is executed if the elements of the transform vector are entered as a 6 element vector
    Mmatlab = reshape(alpha,3,2);
    Mmatlab = [Mmatlab,[0;0;1]];
    A = Mmatlab';
    A(1,3) = Mmatlab(3) + Mmatlab(1)*x_center + Mmatlab(2)*y_center - x_center;
    A(2,3) = Mmatlab(6) + Mmatlab(4)*x_center + Mmatlab(5)*y_center - y_center;
    iMmatlab = inv(Mmatlab);
    iMmatlab(7:9) = [0;0;1];
  end
else
  if strmatch('varstr',PropertyNames)
    varstr = PropertyVal{strmatch('varstr',PropertyNames)};
    if strmatch('varvals',PropertyNames)
      varvals = PropertyVal{strmatch('varvals',PropertyNames)};
      if length(varvals) == length(varstr)
        alpha = change_alpha(alpha,varstr,varvals);
      else
        error('Your "varstr" variables do not match your "varvals"')
      end
    else
      error('You have not specified variable values in Property "varvals"')
    end
  else
    varstr = '';
    if strmatch('varvals',PropertyNames)
      error('You have not specified the variable strings and order in Property "varstr"')
    end
  end
  
  if strmatch('fixedstr',PropertyNames)
    fixedstr = PropertyVal{strmatch('fixedstr',PropertyNames)};
    if strmatch('fixedvals',PropertyNames)
      fixedvals = PropertyVal{strmatch('fixedvals',PropertyNames)};
      if length(fixedvals) == length(fixedstr)
        alpha = change_alpha(alpha,fixedstr,fixedvals);
      else
        error('Your "fixedstr" variables do not match your "fixedvals"')
      end
    else
      error('You have not specified variable values in Property "fixedvals"')
    end
  else
    fixedstr = '';
    if strmatch('fixedvals',PropertyNames)
      error('You have not specified the variable strings and order in Property "fixedstr"')
    end
  end
  
  [A,Mmatlab,iMmatlab] = alpha2tmat(alpha,M+1,N+1);
end


% ------------------------------------------------------------------------
% Apply transformation to the image
% ------------------------------------------------------------------------
if sum(sum(A ~= eye(3))) ~= 0 % Checks if the input transform was an identity transform
  Ap = A';
  a = reshape(Ap(1:6),1,6);
  x_offset = a(3)+x_center;%x_center + a(3) + 0.5;
  y_offset = a(6)+y_center;%y_center + a(6) + 0.5;

  T = [a(1),a(2),x_offset; a(4),a(5),y_offset;0,0,1];
  %   Mout = round(M/a(5));
  %   voffset = 0^xor(rem(M,2),rem(Mout,2));
  %   Nout = round(N/a(1));
  %   uoffset = 0^xor(rem(N,2),rem(Nout,2));
  Tinv = inv(T);
  
  if strmatch('same',outimsz)
    xstart = 1; xend = N;
    ystart = 1; yend = M;
    ustart = 1; uend = N;
    vstart = 1; vend = M;
    Iuv = zeros(size(Ixy));
  elseif strmatch('valid',outimsz)
    urng = [0.5,N+0.5];
    vrng = [0.5,M+0.5];
    XYmap2UV_start = Tinv*[urng(1);vrng(1);1] + T(:,3) + 1 - 1e-6.*ones(3,1);
    if  XYmap2UV_start(1) < 0
      ustart = 1;
      xstart = -1*abs(XYmap2UV_start(1))-0.5;%round(abs(XYmap2UV_start(1)));
    else 
      ustart = round(XYmap2UV_start(1));%XYmap2UV_start(1);%
      xstart = 1;
    end
    
    if  XYmap2UV_start(2) < 0
      vstart = 1;
      ystart = -1*abs(XYmap2UV_start(2))-0.5;%round(abs(XYmap2UV_start(2)));
    else 
      vstart = round(XYmap2UV_start(2));%XYmap2UV_start(2);%
      ystart = 1;
    end
    
    XYmap2UV_end = Tinv*[urng(2);vrng(2);1] + T(:,3) + 1 - 1e-6.*ones(3,1);
    if  XYmap2UV_end(1) > N
      xend = XYmap2UV_end(1)-1;%round(XYmap2UV_end(1))-1;
      uend = xend + abs(xstart) + 0.5;
    else 
      uend = round(XYmap2UV_end(1))-1;%XYmap2UV_end(1) - 1;%
      xend = N;
    end
    
    if  XYmap2UV_end(2) > M
      yend = XYmap2UV_end(2) - 1;%round(XYmap2UV_end(2))-1;
      vend = yend + abs(ystart) + 0.5;
    else 
      vend = round(XYmap2UV_end(2))-1;%XYmap2UV_end(2) - 1;%
      yend = M;
    end
  end
  
  [x, y] = meshgrid(xstart:xend, ystart:yend);
  u = a(1) * (x - x_center) + a(2) * (y - y_center) + x_offset;
  %u = u + uoffset*(u(1,2) - u(1,1))/2;
  v = a(4) * (x - x_center) + a(5) * (y - y_center) + y_offset;
  %v = v + voffset*(v(2,1) - v(1,1))/2;
  
  for k = 1:d
    Iuv(:,:,k) = interp2(1:N, 1:M, Ixy(:,:,k), u, v,interpstr,extrapval);
  end
  
  Iuv = Iuv(vstart:vend,ustart:uend,:);
else % If identity transform in, send the identical image out
  Iuv = Ixy;
end


if strmatch('outputims',PropertyNames)
  sz = PropertyVal{strmatch('outputims',PropertyNames)};
  outy = sz(1); outx = sz(2);
  reszM = ones(outy,outx,d).*extrapval;
  outx_center = (outx+1) / 2; % Changed to be xsize+1 since this give true center
  outy_center = (outy+1) / 2;
  ycentdiff = round(y_center - outy_center);%0;%
  xcentdiff = round(x_center - outx_center);%0;%
  if outy == M
    reszMyind = 1:M;
    Myind = 1:M;
  elseif outy < M
    reszMyind = 1:outy;
    Myind = (1+ycentdiff):(outy +ycentdiff);
  else
    Myind = 1:M;
    reszMyind = floor(outy_center - y_center + 1):(M + floor(outy_center - y_center + 1)-1);
  end
  
  if outx == N
    reszMxind = 1:N;
    Mxind = 1:N;
  elseif outx < N
    reszMxind = 1:outx;
    Mxind = (1+xcentdiff):(outx +xcentdiff);
  else
    Mxind = 1:N;
    reszMxind = floor(outx_center - x_center + 1):(N + floor(outx_center - x_center + 1)-1);
  end
  reszM(reszMyind,reszMxind,:) = Iuv(Myind,Mxind,:);
  Iuv = reszM;
  Iuv(isnan(Iuv)) = 0;
end

  function alpha_out = change_alpha(alpha_in,varstr,varvals);
    alpha_out = alpha_in;
    paramstr = {'tx','ty','theta','sx','sy','gx','gy'};
    for s = 1:length(varstr)
      i = strmatch(varstr{s},paramstr);
      alpha_out(i) = varvals(s);
    end
  end
end
