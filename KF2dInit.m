function param = KF2dInit(varargin)

%% Script: KF2dInit
% Description: Adaptation of kalmanFilterForTracking MATLAB function for
%   high speed video tracking
% The following variable definitions are used to describe the Kalman Filter
% I. System
%   x+  = A * x + B * u + w
%   y   = C * x + v
% 
%   x   : state vector          -> P  : cov. matrix
%   u   : control value
%   w   : perturbation vector   -> Q  : process cov. matrix
%   y   : measurment vector 
%   v   : measurment noise      -> R  : meas. cov. matrix
% 
%   A  : state transition matrix
%   B  : control matrix
%   C  : measurement matrix
% 
% Example:
% 
% INPUTS ----------------------------------------------------------------
%
% varargin - 'PropertyName','PropertyValue'
%   > stateTrans: state transition matrix     => A
%   > procCov: process covariance matrix      => Q = weighted by kinematic
%   eq
%   > procSig: process standard deviation     => siga = sqrt(1e-5)
%   > obsMtx: matrix indicating states that are observed    => C
%   > obsCov: observation covariance matrix  => R
%   > obsSig: observation standard deviation => sigm = sqrt(1)
%   > highestDir: integeter indicating the highest order derivative of x
%   and y position that should be included in the state space
%                                            => ndir = 1
%   > dt
%   > pixel scale
%   > controlVal: control value              => u = 0
% 
%
% OUTPUTS ---------------------------------------------------------------
%   param:  strucutre with fields containing necessary information for a 2D
%     kalman filtering model
%
%  Date           Author              E-mail                      Version
% 21 Sept 2016   Amanda Balderrama  amandagbalderrama@gmail.com     0
%

startvarin = 1;
PropertyNames = lower(varargin(startvarin:2:length(varargin)));
PropertyVal = varargin(startvarin+1:2:length(varargin));
Ain = []; Bin = []; Qin = []; siga = sqrt(1e-5); 
Cin = []; Rin = []; sigm = sqrt(1);
ndir = 1; dt = 1; pixelscale = 1; u = 0;

if strmatch('statetrans',PropertyNames)
  Ain = PropertyVal{strmatch('statetrans',PropertyNames)};
end

if strmatch('procc',PropertyNames) % process covariance matrix
  Qin = PropertyVal{strmatch('procc',PropertyNames)};
end

if strmatch('procs',PropertyNames) % process covariance matrix
  siga = PropertyVal{strmatch('procs',PropertyNames)};
end

if strmatch('obsm',PropertyNames) % measurement covariance matrix
  Cin = PropertyVal{strmatch('obsm',PropertyNames)};
end

if strmatch('obsc',PropertyNames) % measurement covariance matrix
  Rin = PropertyVal{strmatch('obsc',PropertyNames)};
end

if strmatch('obss',PropertyNames) % process covariance matrix
  sigm = PropertyVal{strmatch('obss',PropertyNames)};
end

if strmatch('high',PropertyNames) % process covariance matrix
  ndir = PropertyVal{strmatch('high',PropertyNames)};
end

if strmatch('dt',PropertyNames) % process covariance matrix
  dt = PropertyVal{strmatch('dt',PropertyNames)};
end

if strmatch('pix',PropertyNames) % process covariance matrix
  pixelscale = PropertyVal{strmatch('pix',PropertyNames)};
end

if strmatch('cont',PropertyNames) % process covariance matrix
  u = PropertyVal{strmatch('cont',PropertyNames)};
end

param.motionModel           = 'ConstantAcceleration'; % -- used for
param.initialLocation       = 'Same as first detection';
param.initialEstimateError  = 1E5 * ones(1, 3);
param.motionNoise           = [2500, 100, 10];
param.measurementNoise      = 0.5;
param.frameRange = [200,3e3]; % starts tracking procedure after 5 ms - 30 ms. entire duration of HSV experiments is typically 250 msec or 25000 frames
param.detectionMethod       = 'kmeans';
param.frameRate = 1e5; % frame rate of the video in frames/sec
param.metersppxl = 3.7007e-4;
estr = 'xy';
icnt = 0;
for i = 1:length(estr)
  icnt = icnt + 1;
  eval(sprintf('param.%si = icnt;',estr(i)));
  for d = 1:ndir
    icnt = icnt + 1;
    eval(sprintf('param.%s%si = icnt;',repmat('d',1,d),estr(i)));
  end
end

%--- Set up state transition model parameters
param.dt = dt;
param.pixelscale = pixelscale;
param.controlval = u;
param.nx = (ndir+1)*2;
param.ny = 2;
A = zeros(param.nx,param.nx);
Q = zeros(param.nx,param.nx);
B = zeros(param.nx,1);
C = zeros(param.ny,param.nx);
for i = 1:length(estr)
  for d = 1:(ndir+1)
    for c = d:(ndir+1)
      eval(sprintf('A(param.%s%si,param.%s%si) = dt^(c-d)*(1/factorial(c-d));',repmat('d',1,d-1),estr(i),repmat('d',1,c-1),estr(i)))
    end
  end
   
  for d = 1:(ndir+1)
    for c = 1:(ndir+1)
      eval(sprintf('Q(param.%s%si,param.%s%si) = (dt^(ndir+1-(d-1))/factorial(ndir+1-(d-1))) * (dt^(ndir+1-(c-1))/factorial(ndir+1-(c-1)));',repmat('d',1,d-1),estr(i),repmat('d',1,c-1),estr(i)));
      eval(sprintf('Q(param.%s%si,param.%s%si) = (dt^(ndir+1-(d-1))/factorial(ndir+1-(d-1))) * (dt^(ndir+1-(c-1))/factorial(ndir+1-(c-1)));',repmat('d',1,d-1),estr(i),repmat('d',1,c-1),estr(i)));
    end
  end
  
  for d = 1:(ndir+1)
    eval(sprintf('B(param.%s%si) = (dt^(ndir+1-(d-1))/factorial(ndir+1-(d-1)));',repmat('d',1,d-1),estr(i)));
  end
  
  eval(sprintf('C(%d,param.%si) = 1;',i,estr(i)));
end

if ~isempty(Bin)
  B = Bin;p.C
end
param.B = B;

param.siga = siga; % This can be considered the amount of variability in the acceleration
if ~isempty(Qin)
  Q = Qin;
end
param.Q = Q.*param.siga^2;

if ~isempty(Ain)
  A = Ain;
end
param.A = A;

param.Pinit = eye(param.nx).*1e5;
param.sys = @(k, xkm1, uk) ...
  A * xkm1 + B * param.controlval + ...
  ones(param.nx,1) * uk; % attempt at creating a "constant jerk" model with mean intensity tracking and side of block tracking

%--- Set up observation model parameters
param.sigm = sigm;
param.C = C;
if ~isempty(Cin)
  param.C = Cin;
else

if ~isempty(Rin)
  param.R = Rin;
else
  param.R = eye(param.ny).*param.sigm^2;
end
param.obs = @(k,xk,vk) ...
  C * xk + ones(param.ny,1)*vk;
end
