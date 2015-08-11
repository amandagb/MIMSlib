%% Script:   boundary = tissue_bkgnd_seg(data,elem_ana_range,varargin);
% Description: Chan-Vese Implementation
% By Rami C., Technion - Israel Institute of Technology
% email: rc@tx.technion.ac.il, website: tx.technion.ac.il/~rc
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
%   Types of methods
%       'chen'      = general CV method
%       'vector'    = CV method for vector image
%       'multiphase'= CV method for multiphase (2 phases applied here)
%
% OUTPUTS ---------------------------------------------------------------
% phi0        = updated level set function
%
%--------------------------------------------------------------------------
%  Date           Author            E-mail                      Version
%  11 Nov 2011    Amanda Gaudreau   amanda.gaudreau@gmail.com     1
%  13 Nov 2011    Amanda Gaudreau   amanda.gaudreau@gmail.com     2
%==========================================================================

function seg = CV(I,varargin)
% For academic use only. Please don't distribute this code without a
% permission from the author.
% read: http://tx.technion.ac.il/~rc/Copyright.htm
%% Initialization of Parameters
tic;
PropertyNames = lower(varargin(1:2:length(varargin)));
PropertyVal = varargin(2:2:length(varargin));

% More parameters, as defined in my report.
% Default values: epsilon=1, mu=0.5, ni=0, dt=0.1, lambda1=1, lambda2=1,
%                 h=1, p=1
h=1;
epsilon=1;
p=1;

% Algorithm Specifictions
if strmatch('method',PropertyNames)
  method = PropertyVal{strmatch('method',PropertyNames)};
else
  method = 'chan';
end

if strmatch('mask',PropertyNames)
  mask = PropertyVal(strmatch('mask',PropertyNames));
else
  mask = 'large';
end

if strmatch('iter',PropertyNames)
  iter_max = PropertyVal{strmatch('iter',PropertyNames)};
else
  iter_max = 5000;
end

% Increase inner_iter if the image is noisy.
% This may lead to a longer segmentation process,
% so you should take care of iter_max
inner_iter=1;

if strmatch('epsilon',PropertyNames)
  epsilon = PropertyVal(strmatch('epsilon',PropertyNames));
else
  epsilon = 1e-5;
end

if strmatch('h_type',PropertyNames)
  Htype = PropertyVal(strmatch('h_type',PropertyNames));
else
  Htype =2;
end

if strmatch('mu',PropertyNames)
  mu = PropertyVal{strmatch('mu',PropertyNames)};
else
  mu=0.5;
end

if strmatch('nu',PropertyNames)
  nu = PropertyVal(strmatch('nu',PropertyNames));
else
  nu=0;
end

if strmatch('lambda1',PropertyNames)
  lambda1 = PropertyVal{strmatch('lambda1',PropertyNames)};
else
  lambda1 = 1;
end

if strmatch('lambda2',PropertyNames)
  lambda2 = PropertyVal{strmatch('lambda2',PropertyNames)};
else
  lambda2 = 1;
end

if strmatch('stepsize',PropertyNames)
  dt = PropertyVal{strmatch('stepsize',PropertyNames)};
else
  dt = 0.1;
end

% Plotting Notes
if strmatch('plot_title',PropertyNames)
  plot_title = PropertyVal{strmatch('plot_title',PropertyNames)};
else
  plot_title = '';
end
fntsz = 14;

title_str2 = sprintf('$$\\lambda_1 = %g, \\lambda_2 = %g, \\mu = %g, \\nu = %g,  dt = %g, H method = %g$$',...
  lambda1, lambda2, mu, nu, dt, Htype);
%-- End Initialization

%% Initializing mask
if strmatch('mask',PropertyNames)
  mask = PropertyVal{strmatch('mask',PropertyNames)};
else
  mask = 'large';
end

switch lower (mask)
  case 'small'
    mask = init_mask(I,'small');
  case 'medium'
    mask = init_mask(I,'medium');
  case 'large'
    mask = init_mask(I,'large');
  case 'whole'
    mask = init_mask(I,'whole');
  case 'whole+small'
    m1 = init_mask(I,'whole');
    m2 = init_mask(I,'small');
    mask = zeros(size(I,1),size(I,2),2);
    mask(:,:,1) = m1(:,:,1);
    mask(:,:,2) = m2(:,:,2);
  otherwise
    error('unrecognized mask shape name (MASK).');
end
phi = mask;

rect_I = I;
rect_I = sum(rect_I,3);
tic; % initializing time counter

%% Defining the initial contour - a circle in the center of the image
[rows,cols]=size(rect_I(:,:,1));
cent_rows=floor(rows/2);
cent_cols=floor(cols/2);
radius=cent_rows/3;

[j,i] = meshgrid(1:cols,1:rows);
phi= repmat(radius,rows,cols)-((i-cent_rows).^2+(j-cent_cols).^2).^0.5;

%% figures
figure('name','Segmentation process');
imagesc(rect_I);  colormap(gray)
title('Initial state');
hold on;
contour(phi,[0 0],'g');

pause(0.5);
figure(2);

%% Main loop
fin=0;
for i=1:iter_max
  
  phi_old = phi;
  h_side_phi = 0.5*(1+(2/pi)*atan(phi/epsilon));;
  
  c1=sum(sum(h_side_phi.*rect_I))/sum(sum(h_side_phi));
  c2=sum(sum((1-h_side_phi).*rect_I))/sum(sum(1-h_side_phi));
  delta = (1/pi)*epsilon./(epsilon^2+phi.^2);
  m=dt*delta;
  [phi_x,phi_y]=gradient(phi);
  abs_grad=(phi_x.^2 +phi_y.^2).^(0.5);
  delta = (1/pi)*h./(h^2+phi.^2);
  len_phi=sum(sum(delta.*abs_grad));
  
  % Each C_i is defined accoring to my report, eps (matlab const.)
  % is used to prevent division by 0
  
  for q=1:inner_iter
    C_1 = 1./((phi(:,[2:cols,cols])-phi+eps).^2 ...
      + 0.25*(phi([2:rows,rows],:)-phi([1,1:rows-1],:)+eps).^2).^(0.5);
    C_2 = 1./((phi-phi(:,[1,1:cols-1])+eps).^2 ...
      +0.25*(phi([2:rows,rows],[1,1:cols-1])...
      -phi([1,1:rows-1],[1,1:cols-1])+eps).^2).^(0.5);
    C_3 = 1./((phi([2:rows,rows],:)-phi+eps).^2 ...
      +0.25*(phi(:,[2:cols,cols])-phi(:,[1,1:cols-1])+eps).^2).^(0.5);
    C_4 = 1./((phi-phi([1,1:rows-1],:)+eps).^2 ...
      +0.25*(phi([1,1:rows-1],[2:cols,cols]) ...
      -phi([1,1:rows-1],[1,1:cols-1])+eps).^2).^(0.5);
    
    
    C = 1+p*mu*m*len_phi^(p-1).*(C_1+C_2+C_3+C_4);
    phi = (phi+p*mu*m*len_phi^(p-1).*(C_1.*phi(:,[2:cols,cols]) ...
      + C_2.*phi(:,[1,1:cols-1]) + C_3.*phi([2:rows,rows],:) ...
      + C_4.*phi([1,1:rows-1],:) ) ...
      + m.*(nu+lambda1*(rect_I-c2).^2-lambda2*(rect_I-c1).^2))./C;
    
  end
  
  % Reinitialization of phi
  phi=reinit(phi,dt,h);
  phi_new=phi;
  
  %    Displaying current state after each 20 iterations
  if (mod(i,20)==0)
    figure(2);
    imagesc(rect_I);
    contour(phi,[0 0],'g');
    title('Current segmentation');
    xlabel(['Time: ' num2str(toc) ' seconds, ' num2str(i) ' iterations']);
    pause(0.1);
    
    [Q,fin]=cur_diff(phi_old,phi_new,dt,h);
    
    if fin && i>59
      break
    end
    
  end
  
  
end

%% Final segmentation results
close 2;
figure(2);
suptitle('Final segmentation');
subplot(211);
imagesc(uint8(rect_I));
hold on;
contour(phi,[0 0],'g');

seg_image=rect_I;
seg_image(phi>0)=c1;
seg_image(phi<0)=c2;
subplot(212);
imagesc(seg_image,[]);
title('Segments');
xlabel(['Time elpased: ' num2str(toc) ' seconds, ' num2str(i) ' iterations']);
end
%% ------------------------------------------------------------------
%   SUBFUNCTIONS
% -------------------------------------------------------------------
function phi = reinit(phi,dt,h)
mask=zeros(size(phi));
sign_mat=zeros(size(phi));
max_iter=25;

Q=10e6;
for i=1:10
  %i=i+1;
  phi_old=phi;
  padded_phi = padarray(phi,[1 1]);
  
  % the constants below are defined in section 4.1 of my report
  
  a = (phi-padded_phi(1:end-2,2:end-1))/h;
  b = (padded_phi(3:end,2:end-1)-phi)/h;
  c = (phi-padded_phi(2:end-1,1:end-2))/h;
  d = (padded_phi(2:end-1,3:end)-phi)/h;
  
  a_p = max(a,0);a_m = min(a,0);
  b_p = max(b,0);b_m = min(b,0);
  c_p = max(c,0);c_m = min(c,0);
  d_p = max(d,0);d_m = min(d,0);
  
  G = zeros(size(phi));
  ind_plus = find(phi>0);
  ind_minus = find(phi<0);
  G(ind_plus) = sqrt(max(a_p(ind_plus).^2,b_m(ind_plus).^2)+max(c_p(ind_plus).^2,d_m(ind_plus).^2))-1;
  G(ind_minus) = sqrt(max(a_m(ind_minus).^2,b_p(ind_minus).^2)+max(c_m(ind_minus).^2,d_p(ind_minus).^2))-1;
  
  sign_mat=sign(phi_old);
  phi = phi-dt*sign_mat.*G;
  phi_new=phi;
  % mask=(phi_old<h); M=sum(sum(mask));
  % Q=sum(sum(mask.*abs(phi_new-phi_old)))/M;
  %
  % if i>max_iter
  %     break
  % end
end
end

function [Q,res]=cur_diff(phi_old,phi_new,dt,h)

mask=(phi_old<h); M=sum(sum(mask));
Q=sum(sum(mask.*abs(phi_new-phi_old)))/M;

if Q>dt*h^2
  res=0;
else
  res=1;
end
end

function m = init_mask(I,type)
% auto pick a mask for image I
% built-in mask creation function
% Input: I   : input image
%        type: mask shape keywords
% Output: m  : mask image
%
%  Date           Author            E-mail                      Version
%  25 Sept 2011   Amanda Gaudreau   amanda.gaudreau@gmail.com     1

temp = double(I(:,:,1));

[M,N] = size(temp);
cx = M/2; %centers phi0 at center of x axis
cy = N/2; %centers phi0 at center of y axis
x = repmat([1:M]',1,N); %creates a matrix that is szx by szy has all columns == row number i => aij = i
y = repmat([1:N],M,1); %creates a matrix that is szx by szy has all columns == column number j => aij = j

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
  case 'whole'
    r = 9;
    m = zeros(round(ceil(max(M,N)/(2*(r+1)))*3*(r+1)));
    siz = size(m,1);
    sx = round(siz/2);
    i = 1:round(siz/(2*(r+1)));
    j = 1:round(0.9*(siz/(2*(r+1))));
    j = j-round(median(j));
    m(sx+2*j*(r+1),(2*i-1)*(r+1)) = 1;
    se = strel('disk',r);
    m = imdilate(m,se);
    m = m(round(siz/2-M/2-6):round(siz/2-M/2-6)+M-1,round(siz/2-N/2-6):round(siz/2-N/2-6)+N-1);
end
if ~exist('m','var')
  m = zeros(size(temp));
  %m((x-cx).^2+(y-cy).^2<min(rx,ry).^2) = 1; %Gives circular boundary
  m(ry^2.*(x-cx).^2+rx^2.*(y-cy).^2<(rx^2*ry^2)) = 1; %Gives elliptical boundary
  tem(:,:,1) = m;
  tem(:,:,2) = m;
  m = tem;
else
  tem(:,:,1) = m;
  M = padarray(m,[floor(2/3*r),floor(2/3*r)],0,'post');
  tem(:,:,2) = M(floor(2/3*r)+1:end,floor(2/3*r)+1:end);
  m = tem;
end
end