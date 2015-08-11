function [kval,varargout] = kernel_val(kern_str,x1,x2,bw,varargin)
%% Script: [kval] = kernel_val(kern_str,x,bw)
% transform_image(moving_image, [x_shift, y_shift, rotation (radians), x_scale, y_scale, skew],'transtype','physicalORtrans',...
%   'extrapval',#,'fixedstr',{'tx'},'fixedvals',[#],...
%   'varstr',{'theta'},'varval',[5*pi/180])
% Description:
% Example:
% Required Functions:
% INPUTS ----------------------------------------------------------------
% kern_str: Image on which affine transform will be performed according
%         to parameters of alpha. Only works on grayscale images
% x:  1x6 vector with the transformation parameters OR
%         3x3 matrix transformation matrix
%         [x_shift, y_shift, rotation (radians), x_scale, y_scale, skew]
%         NOTE: Positive displacement is to the left, up, and counterclockwise,
%           respectively.
% varargin - 'PropertyName','PropertyValue'
%   'extrapval': uses EXTRAPVAL to fill in any missing data. [DEFAULT = 0]
%
% OUTPUTS ---------------------------------------------------------------
%
%
%  Date           Author            E-mail                      Version
%  25 Apr 2011    Amanda Gaudreau   amanda.gaudreau@gmail.com     1
%  25 Apr 2011    Amanda Gaudreau   amanda.gaudreau@gmail.com     2
%  26 Apr 2011    Amanda Gaudreau   amanda.gaudreau@gmail.com     3
%  04 Aug 2013    Amanda Gaudreau   amanda.gaudreau@gmail.com     4
%  27 Jan 2014    Amanda Gaudreau   amanda.gaudreau@gmail.com     5
%  29 Jan 2014    Amanda Gaudreau   amanda.gaudreau@gmail.com     6
%   6 Feb 2014    Amanda Gaudreau   amanda.gaudreau@gmail.com     6.1
%     Changed the useind parameter so that only one of the columns needs to
%     be nonzero to be used

PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

if strmatch('extrapval',PropertyNames)
  extrapval = PropertyVal{strmatch('extrapval',PropertyNames)};
else
  extrapval = 0;
end

if strmatch('comp_prod',PropertyNames)
  compprod = PropertyVal{strmatch('comp_prod',PropertyNames)};
else
  compprod = 0;
end

L = size(x1,1);
d = size(x1,2);
kval = zeros(L,1);
kvald = zeros(L,d);
arg = (x1 - x2)./repmat(bw(:)',L,1);
useindLOGIC = abs(arg) < 1;
useind = find(sum(useindLOGIC,2)>0);

% tol=1000*eps;
default=0;
% K=struct('type',char(0),'name',char(0),'coef',[],'support',[-1,1],'nu',NaN,'k',NaN,'mu',NaN,'var',NaN,'beta',NaN);
% if nargin==0 default=1;
% elseif strncmpi(kern_str,'defa',4) | strncmpi(kern_str,'epan',4)
%     default=1;
% end
if default
  K.type='opt';
  K.coef=[-0.75,0,0.75];
  K.nu=0;
  K.k=2;
  K.mu=0;
  p1=K.coef;
  p2=conv(p1,p1);
  p3=polyint(p2);
  K.var=polyval(p3,1)-polyval(p3,-1);
  p2=conv(p1,[1 0 0]);
  p3=polyint(p2);
  K.beta=polyval(p3,1)-polyval(p3,-1);
elseif strncmpi(kern_str,'rect',4)
  kvald(useind,:) = 0.5;
elseif strncmpi(kern_str,'norm',4)
  kvald(useind,:) = exp(-arg.^2./2)/sqrt(2*pi);
elseif strncmpi(kern_str,'epan',4)
  kvald(useind,:) = 0.75.*(1-arg(useind,:).^2);
elseif strncmpi(kern_str,'quar',4)
  K.type='opt';
  K.coef=conv([-1,0,1],[-1,0,1])*15/16;
  K.nu=0;
  K.k=2;
  K.mu=1;
  p1=K.coef;
  p2=conv(p1,p1);
  p3=polyint(p2);
  K.var=polyval(p3,1)-polyval(p3,-1);
  p2=conv(p1,[1 0 0]);
  p3=polyint(p2);
  K.beta=polyval(p3,1)-polyval(p3,-1);
elseif strcmpi(kern_str,'opt')
  K.type='opt';
  nu=par1;
  k=par2;
  mu=par3;
  K.nu=nu;
  K.k=k;
  K.mu=mu;
  K.coef=K_opt(nu,k,mu);
  p1=K.coef;
  p2=conv(p1,p1);
  p3=polyint(p2);
  K.var=polyval(p3,1)-polyval(p3,-1);
  p2=conv(p1,[1 zeros(1,k)]);
  p3=polyint(p2);
  K.beta=polyval(p3,1)-polyval(p3,-1);
elseif strcmpi(kern_str,'pol')
  K.type='pol';
  P=par1;
  a=par2;
  b=par3;
  nu=0;
  P0=P;
  P1=polyint(P0);
  I=polyval(P1,b)-polyval(P1,a);
  while abs(I)<tol
    nu=nu+1;
    P0=conv(P0,[1,0]);
    P1=polyint(P0);
    I=polyval(P1,b)-polyval(P1,a);
  end
  P=(-1)^nu*prod(1:nu)/I*P;
  K.nu=nu;
  k=nu+1;
  P0=conv(P,[1,zeros(1,k)]);
  P1=polyint(P0);
  I=polyval(P1,b)-polyval(P1,a);
  while abs(I)<tol
    k=k+1;
    P0=conv(P0,[1,0]);
    P1=polyint(P0);
    I=polyval(P1,b)-polyval(P1,a);
  end
  K.k=k;
  K.beta=I;
  P0=conv(P,P);
  P1=polyint(P0);
  I=polyval(P1,b)-polyval(P1,a);
  K.var=I;
  mu=-1;
  P0=P;
  while abs(polyval(P0,a)) < tol & abs(polyval(P0,b)) < tol
    mu=mu+1;
    P0=polyder(P0);
  end
  K.mu=mu;
  K.coef=P;
  K.support=[a,b];
elseif strcmpi(kern_str,'str')
  K.type='str';
  fs=strrep(par1,'.*','*');
  fs=strrep(fs,'./','/');
  fs=strrep(fs,'.^','^');
  K.name=fs;
  K.support=[par2,par3];
  if exist('sym')
    a=par2;
    b=par3;
    nu=0;
    f=sym(fs);
    f0=f;
    I=int(f0,a,b);
    while abs(eval(I))<tol
      nu=nu+1;
      f0=f0*sym('x');
      I=int(f0,a,b);
    end
    f=(-1)^nu*prod(1:nu)/I*f;
    fs=char(f);
    K.name=fs;
    K.nu=nu;
    k=nu+1;
    f0=f*sym(['x^',num2str(k)]);
    I=int(f0,a,b);
    while abs(eval(I))<tol
      k=k+1;
      f0=f0*sym('x');
      I=int(f0,a,b);
    end
    K.k=k;
    K.beta=eval(I);
    f0=f*f;
    I=int(f0,a,b);
    K.var=eval(I);
    mu=-1;
    f0=f;
    fa=limit(f0,a);
    fb=limit(f0,b);
    while (abs(eval(fa)) < tol & abs(eval(fb)) < tol) & mu<20
      mu=mu+1;
      f0=diff(f0);
      fa=limit(f0,a);
      fb=limit(f0,b);
    end
    if mu>=20 mu=inf; end
    K.mu=mu;
  end
elseif strcmpi(kern_str,'fun')
  K.type='fun';
  K.name=par1;
  K.support=[par2,par3];
elseif strncmpi(kern_str,'gaus',4) % Gauss kernel
  K.type='gau';
  %    K.name='1/sqrt(2*pi)*exp(-x^2/2)';
  if nargin==1
    K.beta=1;
  else
    K.beta=par1^2;
  end
  s=sqrt(K.beta);
  K.support=[-inf,inf];
  K.nu=0;
  K.k=2;
  K.mu=inf;
  K.var=1/(2*s*sqrt(pi));
  K.coef=1; % used for derivatives of gaussian kernel
  %    K.beta=s^2;
elseif strncmpi(kern_str,'tria',4) % triangular kernel
  K.type='tri';
  K.support=[-1,1];
  K.nu=0;
  K.k=2;
  K.mu=0;
  K.var=2/3;
  K.beta=1/6;
else
  error('Unrecognized parameter')
end

if compprod
  kval(useind) = prod(kvald(useind),2);
else 
  kval = kvald;
end
varargout{1}.Kd = kvald;
varargout{1}.useind = useind;
varargout{1}.useindLOGIC = useindLOGIC;
