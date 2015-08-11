function [bw,varargout] = bandwidth_selector(X,Kstr, method,varargin)
%% Script: [kval] = kernel_val(kern_str,x,bw)
% transform_image(moving_image, [x_shift, y_shift, rotation (radians), x_scale, y_scale, skew],'transtype','physicalORtrans',...
%   'extrapval',#,'fixedstr',{'tx'},'fixedvals',[#],...
%   'varstr',{'theta'},'varval',[5*pi/180])
% Description:
% Example:
% Required Functions:
% INPUTS ----------------------------------------------------------------
% • X:
% • 'kernel': {'string',bandwidth,dimensions} cell identifying and
%   parameterizing the kernel that should be used to construct the pdf
%   estimates. All kernels can be specified or entered as 1-d separable
%   kernels. All kernels are separable second order kernels. The
%   (dA+dB)-dim kernel will be constructed as the product of the separable
%   1-d kernel. Kernels have a functional form K(u) where
%             f_hat(x) = 1/(bw*nsamp)*sum(K((SAMPLES - x)/bw))
%     - 'string': specify the kernel, K(u) => 'norm', 'rect', 'epan', 'b3'
%     [DEFAULT = 'norm'] - bandwidth: integer specifying the bandwidth (see
%     f_hat above) [DEFAULT = Silverman's Rule-of-Thumb =
%     sig_hat*C(K)*n^(-1/5)] - dimensions: indicates the dimension the
%     spline should be applied to
%   Multiple kernels can be chosen by letting the varaible have multiple
%   rows:
%             {'kernel1',[bw_d1,bw_d2,...],[d1,d2,...];
%             'kernel2',[bw_d1,bw_d2,...],[d1,d2,...]}
% varargin - 'PropertyName','PropertyValue'
%   'extrapval': uses EXTRAPVAL to fill in any missing data. [DEFAULT = 0]
%
% OUTPUTS ---------------------------------------------------------------
%
%
%  Date           Author            E-mail                      Version
%  11 Feb 2014    Amanda Gaudreau   amanda.gaudreau@gmail.com     0

PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

if strmatch('extrapval',PropertyNames)
  extrapval = PropertyVal{strmatch('extrapval',PropertyNames)};
else
  extrapval = 0;
end
method = 'srot';
switch lower(method)
  case 'srot'
    nu = 2;
    switch lower(Kstr)
      case 'norm'
        RK = 1/(2*sqrt(pi));
        kap2 = 1;
        K = @(u) exp(-u.^2./2)./sqrt(2*pi);
        %Ksupp(:,Kdim) = [-1;1];
      case 'rect'
        RK = 1/2;
        kap2 = 1/3;
        K = @(u) 0.5.*(abs(u) <= 1);
        Ksupp = [-1;1];
      case 'epan'
        RK = 3/5;
        kap2 = 1/5;
        K = @(u) 0.75.*(1-u.^2).*(abs(u) <= 1);
        Ksupp = [-1;1];
      case 'b3'
        RK = 1/3;
        kap2 = 151/315;
        K = @(u) (1/6).*((4-6*u.^2+3.*abs(u)).*(abs(u) < 1) + (2-abs(u)).^3.*(abs(u) >=1 & abs(u) < 2));
        Ksupp = [-2;2];
    end
    
    sig_hat = std(X);
    CK = 2*((pi^(1/2)*factorial(nu)^3*RK)/(2*nu*factorial(2*nu)*kap2^2))^(1/(2*nu+1));
    bw = sig_hat.*CK.*size(X,1)^(-1/(2*nu+1));
end