function [param, utilities] = HSVtrackingKF(filename,varargin)

%% Script: HSVtracking
% Description: Adaptation of kalmanFilterForTracking MATLAB function for
%   high speed video tracking
% Example:
% INPUTS ----------------------------------------------------------------
%
% varargin - 'PropertyName','PropertyValue'
%   > objectdetions:    N x 2 vector with the x, y pairs in the rows of the
%       centroid position of the object. This eliminates the need to read
%       frames from the video file and 
%   > utilities:      structure of previous runs of HSVtracking function
%   > track:    logical indicating whether to run filtering algorithm for
%   tracking [DEFAULT = 1]
%   > dframerange:   1 x 2 vector indicating the starting and ending frame
%   to conduct detection of object 
%   > tframerange:   1 x 2 vector indicating the starting and ending frame
%   to conduct tracking anlsysi on
%
% OUTPUTS ---------------------------------------------------------------
%
%  Date           Author              E-mail                      Version
% 18 Sept 2016   Amanda Balderrama  amandagbalderrama@gmail.com     0
%   adapted from MATLAB function and HSVtrackingPF function
% 20 Sept 2016   Amanda Balderrama  amandagbalderrama@gmail.com     1
%   More basic version of function with variables and parameters presented
%   up front
  
cdir = loaddirfun;
vidstr = '20130605_0910_HSV69_B1_C1_convert'; % 'HSI_12_C1H1_convert'; %
% filename = strcat(cdir.fpath,'Impacter Nose Tracking\HSV69_B1_C1 Comparison\blastHSV69_B1_C1.avi');
db = 0;
verbose = 0;
if ~exist('filename')
  filename = strcat(cdir.dboxpath(1:end-7),'Desktop\HSVtests\',vidstr,'.avi'); 
elseif isempty(filename)
  filename = strcat(cdir.dboxpath(1:end-7),'Desktop\HSVtests\',vidstr,'.avi'); 
end

%% Initialize variables
startvarin = 1;
PropertyNames = lower(varargin(startvarin:2:length(varargin)));
PropertyVal = varargin(startvarin+1:2:length(varargin));

if strmatch('detect',PropertyNames)
  detect = PropertyVal{strmatch('detect',PropertyNames)};
else 
  detect = 1;
end

if strmatch('track',PropertyNames)
  track = PropertyVal{strmatch('track',PropertyNames)};
else 
  track = 1;
end

if strmatch('dfr',PropertyNames)
  utilities.detfr = PropertyVal{strmatch('dfr',PropertyNames)};
else 
  utilities.detfr = [-1,Inf];
end

if strmatch('tfr',PropertyNames)
  utilities.tfr = PropertyVal{strmatch('tfr',PropertyNames)};
else 
  utilities.tfr =  [];
end

%-------------- initialize filter parameters
if strmatch('param',PropertyNames)
  param = PropertyVal{strmatch('param',PropertyNames)};
else 
  param = KF2dInit(varargin{:});  % get Kalman configuration that works well
end

%-------------- initialize recordering parameters
if strmatch('utilities',PropertyNames)
  utilities = PropertyVal{strmatch('utilities',PropertyNames)};
  if any(isnumeric(utilities.accumulatedDetections(:)))
    detect = 0;
  end
else
  utilities.frameidx = 0;
  utilities.videoReader = vision.VideoFileReader(filename);
  vI = info(utilities.videoReader);
  utilities.accumulatedImage      = 0;
  utilities.trackind = 0;
  utilities.accumulatedTrackings  = nan(3e4, 2);
  utilities.detectind = false(3e4,2);% param.frameRange(1)-1;
  utilities.accumulatedDetections = nan(3e4,2);
  utilities.foregroundMinMeanIntensity = nan(3e4,2);
  utilities.matchIntensity = 0;
  utilities.matchPos = 0;
  utilities.perframeMask = cell(3e4,1);
  Iind = 2;
  utilities.states = nan(3e4,param.nx);
  utilities.covar = nan(3e4,param.nx^2);
  utilities.MLstates = nan(3e4,param.nx);
  utilities.MLcovar = nan(3e4,param.nx^2);
end

if strmatch('objectdet',PropertyNames)
  detections = PropertyVal{strmatch('objectdet',PropertyNames)};
  utilities.accumulatedDetections = detections;
  nonNans = find(sum(isnan(utilities.accumulatedDetections),2) == 0);
  utilities.detectind = nonNans(end);
else
  detections = [];
end

%----------------- Separate memory space
y = zeros(param.ny,3e4);
if detect
  fcnt = 0;
  while ~isDone(utilities.videoReader) && fcnt < utilities.detfr(2)
    fcnt = fcnt + 1;
    frame = step(utilities.videoReader);
    utilities.frameidx = fcnt;
    
    if mod(utilities.frameidx,50) == 0 && verbose % debugging sequence to show when frames have been read and processed
      disp(utilities.frameidx);
    end
    
    if fcnt >= utilities.detfr(1)
      if ~isempty(detections)
        xyPos = detections(fcnt,:);
        xyInfo = 1-sum(isnan(xyPos));
        utilities.foregroundMask = zeros(size(frame,1), size(frame,2));
      elseif ~exist('xyInfo')
        [xyPos, xyInfo] = detectHSVobject(frame,'past',utilities,'init',1);
      else
        [xyPos, xyInfo] = detectHSVobject(frame,'past',utilities,...
          'segmethod',param.detectionMethod,...
          'intensity',utilities.matchIntensity,...
          'seed',utilities.accumulatedDetections(fcnt-1,:),...
          'search',15);
      end
      
      utilities.accumulatedImage = max(utilities.accumulatedImage, frame);
      utilities.foregroundMask = xyInfo.foregroundMask;
      utilities.foregroundMinMeanIntensity(fcnt,:) = [xyInfo.foregroundMinIntensity xyInfo.foregroundMeanIntensity];
      utilities.matchIntensity = utilities.foregroundMinMeanIntensity(fcnt,Iind);
      utilities.matchPos = xyPos;
      utilities.perframeMask{fcnt} = find(xyInfo.foregroundMask == 1);
      utilities.detectind(fcnt) = true;
      utilities.accumulatedDetections(fcnt,:) = xyPos;
      
      %----------------- Fill observations vector
      %y(:,fcnt) = ...
      %  [utilities.accumulatedDetections(fcnt,:)'];
      %utilities.foregroundMinIntensity;...
      %sum(utilities.foregroundMask(:))];
    end
  end
  utilities.detfr(2) = fcnt;
end

if isempty(utilities.tfr)
  utilities.tfr = utilities.detfr;
end
nonNans = find(sum(isnan(utilities.accumulatedDetections),2) == 0);
utilities.detfr = [1,nonNans(end)];
T = nonNans(end);  %# non-nan frames steps


if track
  k = 1;
  x = zeros(param.nx,1);
  x([param.xi,param.yi]) = utilities.accumulatedDetections(k,:)';
  P = param.Pinit;
  utilities.accumulatedTrackings(k,:) = x([param.xi,param.yi],k);
  utilities.states(k,:) = x(:)';
  utilities.covar(k,:) = P(:)';
  
  kalmanFilter = configureKalmanFilter(param.motionModel, ...
    utilities.accumulatedDetections(k,:)', param.initialEstimateError, ...
    param.motionNoise, param.measurementNoise);
  trackedLocation = correct(kalmanFilter, utilities.accumulatedDetections(k,:)');
  utilities.MLstates(k,:) = kalmanFilter.State([1,2,4,5]);
  utilities.MLcovar(k,:) = kalmanFilter.StateCovariance([1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23]);
  
  for k = 2:T    
    xmeas = utilities.accumulatedDetections(k,:)';
    % Prediction
    x = param.A * x + param.B * param.controlval;
    P = param.A * P * param.A' + param.Q;
    
    %Update
    K = P*param.C'*inv(param.C*P*param.C' + param.R);
    if ~any(isnan(xmeas))
      x = x + K*(xmeas - param.C*x);
    end
    P = (eye(param.nx) - K*param.C)*P;
    
    utilities.accumulatedTrackings(k,:) = x([param.xi,param.yi]);
    utilities.states(k,:) = x(:)';
    utilities.covar(k,:) = P(:)';
    
    % Use the Kalman filter to track the ball.
    if ~any(isnan(utilities.accumulatedDetections(k,:)'))%isObjectDetected % The ball was detected.
      % Reduce the measurement noise by calling predict followed by
      % correct.
      %predict(kalmanFilter)
      trackedLocation = correct(kalmanFilter, utilities.accumulatedDetections(k,:)');
      label = 'Corrected';
    else % The ball was missing.
      % Predict the ball's location.
      trackedLocation = predict(kalmanFilter);
      label = 'Predicted';
    end
    utilities.MLstates(k,:) = kalmanFilter.State([1,2,4,5]);
    utilities.MLcovar(k,:) = kalmanFilter.StateCovariance([1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23]);
    
    if mod(k,200) == 0 && db
      figure('position',[117,26,1036,901],'name',sprintf('Iter%d',k)); subplot(3,3,2);
      plot(y(1,1:k-1),'.');  hold on; plot(utilities.states(1:k-1,param.xi),'.');  plot(utilities.MLstates(1:k-1,param.xi),'.');
      title('x Pos')
      subplot(3,3,5); plot([nan(1),diff(y(1,1:k-1))],'.');  hold on; plot(utilities.states(1:k-1,param.dxi),'.'); plot(utilities.MLstates(1:k-1,param.dxi),'.');
      title('x Vel')
      subplot(3,3,8); plot([nan(1,2),diff(diff(y(1,1:k-1)))],'.');  hold on; plot([nan(1); diff(utilities.states(1:k-1,param.dxi))],'.');
      title('x Acc')
      
      subplot(3,3,3);
      plot(y(2,1:k-1),'.');  hold on; plot(utilities.states(1:k-1,param.yi),'.');  plot(utilities.MLstates(1:k-1,param.yi),'.');
      title('y Pos')
      subplot(3,3,6); plot([nan(1),diff(y(2,1:k-1))],'.');  hold on; plot(utilities.states(1:k-1,param.dyi),'.'); plot(utilities.MLstates(1:k-1,param.dyi),'.');
      title('y Vel')
      subplot(3,3,9); plot([nan(1,2),diff(diff(y(2,1:k-1)))],'.');  hold on; plot([nan(1); diff(utilities.states(1:k-1,param.dyi))],'.');
      title('y Acc')
      
      subplot(3,3,1); plot(y(1,1:k-1),y(2,1:k-1),'.',...
        utilities.states(1:k-1,param.xi),utilities.states(1:k-1,param.yi),'.',...
        utilities.MLstates(1:k-1,param.xi),utilities.MLstates(1:k-1,param.yi),'.');
      xlabel('x Pos'); ylabel('y Pos');
    end
  end
end


% showTrajectory();
%%
%{.
detectedind = find(sum(isnan(utilities.accumulatedDetections),2) == 0);
trackedind = find(sum(isnan(utilities.accumulatedTrackings),2) == 0);
fs = 100000;
fc = 2000;
[b,a] = butter(2,fc/(fs/2));
pstr = 'xy';
rowstr = {'Raw','Butter Filt','Kalman Filt'};
colstr = {'Pos','Vel','Acc','Jerk'};
lcol = lines(3);

for p = 1:2
  alldat{p} = zeros(length(utilities.accumulatedDetections(detectedind,p)),12);
  adi = 1;        alldat{p}(:,adi) = utilities.accumulatedDetections(detectedind,p);
  adi = adi + 1;  alldat{p}(:,adi) = [nan(1);diff(alldat{p}(:,adi-1))];
  adi = adi + 1;  alldat{p}(:,adi) = [nan(1);diff(alldat{p}(:,adi-1))];
  adi = adi + 1;  alldat{p}(:,adi) = [nan(1);diff(alldat{p}(:,adi-1))];
  
  adi = adi + 1;  alldat{p}(:,adi) = filtfilt(b,a,alldat{p}(:,1));
  adi = adi + 1;  alldat{p}(:,adi) = [nan(1);diff(alldat{p}(:,adi-1))];
  adi = adi + 1;  alldat{p}(:,adi) = [nan(1);diff(alldat{p}(:,adi-1))];
  adi = adi + 1;  alldat{p}(:,adi) = [nan(1);diff(alldat{p}(:,adi-1))];
  
  adi = adi + 1;  eval(sprintf('alldat{%d}(:,adi) = utilities.states(detectedind,param.%si);',p,pstr(p)));
  adi = adi + 1;  eval(sprintf('alldat{%d}(:,adi) = utilities.states(detectedind,param.d%si);',p,pstr(p)));
  adi = adi + 1;  alldat{p}(:,adi) = [nan(1);diff(alldat{p}(:,adi-1))];
  adi = adi + 1;  alldat{p}(:,adi) = [nan(1);diff(alldat{p}(:,adi-1))];
end
tv = ((utilities.detfr(1)-1):(utilities.detfr(2)-1));%./(param.frameRate*1e-3);
  
%=============================
% plots x and y on separate figures and shows raw, butter and KF in their
% own plots
% close all
for p = 1%:2
  %dl = alldat{p}(:,9);
  figure('position',[220, 220, 1545, 784],'name','HSVxder');
  i = 0;
  yl = zeros(4,2);
  xl = zeros(4,2);
  for r = 1:3
    for j = 1:4
      i = i + 1;
      subplot(3,4,i); plot(tv,alldat{p}(:,i),'.','color',lcol(r,:)); axis tight; hold on;
      %if r == 3 && j > 1
      %  dl = [nan;diff(dl)];
      %  plot(dl,'.');
      %end
      axis tight; title(sprintf('%s %s-%s',rowstr{r},pstr(p),colstr{j}));
      switch r
        case 2
          yl(j,:) = ylim;
          xl(j,:) = xlim;
        case 3
          ylim(yl(j,:));
          xlim(xl(j,:));
          xlabel('Frame #')
      end
      
      switch j
        case 1
          ylabel('Pixel Pos.')
        case 2
          ylabel('Pixels/frame');
        case 3
          ylabel('Pixels/frame^2')
        case 4
          ylabel('Pixels/frame^3')
      end
    end
  end
end
% save_open_figures('C:\Users\ADGB\Dropbox\MADLab Data\ThesisFig',[],[],'','fmt','png');

%=============================
% % plots all together
% figure('position',[220          220        1545         784]);
% i = 4;
% for p = 1:2
%   for j = 1:4
%     i = i + 1;
%     subplot(3,4,i); plot(alldat{p}(:,[j]),'.'); hold on;
%     plot(alldat{p}(:,[j+4,j+8]),'.'); hold on;
%     axis tight;
%     axis tight; title(sprintf('%s %s',pstr(p),colstr{j}));
%   end
% end
% subplot(3,4,[1,2,3]); plot(alldat{1}(:,1),alldat{2}(:,1)); hold on;
% plot(alldat{1}(:,5),alldat{2}(:,5),'.',...
%   alldat{1}(:,9),alldat{2}(:,9),'.');
% legend(rowstr,'location','northeastoutside');

% close all
figure('name','HSVdetTrackRaw','position',[680   724   437   374]);%[680   853   560   245]); 
plot(alldat{1}(:,1),alldat{2}(:,1)); hold on;
plot(alldat{1}(:,5),alldat{2}(:,5),'o','markersize',1);
plot(alldat{1}(:,9),alldat{2}(:,9),'o','markersize',1); grid on;
axis tight;
axis ij; %axis image;
xlabel('x pixel location (column)');
ylabel('y pixel location (row)');
legend(rowstr,'location','northeast');
% save_open_figures('C:\Users\ADGB\Dropbox\MADLab Data\ThesisFig',[],[],'','fmt','png');


% figure; plot3(alldat{1}(:,1),alldat{2}(:,1),tv); hold on;
% plot(alldat{1}(:,5),alldat{2}(:,5),'.',...
%   alldat{1}(:,9),alldat{2}(:,9),'.'); grid on;
% legend(rowstr,'location','northwest');
%}

%%
% Show trajectory of the ball by overlaying all video frames on top of
% each other.
  function showTrajectory
    % Close the window which was used to show individual video frame.
    uiscopes.close('All');
    
    % Create a figure to show the processing results for all video frames.
    figure; imshow(utilities.accumulatedImage/2+0.5); hold on;
    plot(utilities.accumulatedDetections(:,1), ...
      utilities.accumulatedDetections(:,2), 'r.');
    
    if ~isempty(utilities.accumulatedTrackings)
      plot(utilities.accumulatedTrackings(:,1), ...
        utilities.accumulatedTrackings(:,2), 'g-o');
      legend('Detection', 'Tracking');
    end
    set(gca,'position',[0,0,1,1])
  end


%% 

  function pf = setDefaultPF(nFrames,param)
    %----------------- Transition prior PDF p(x[k] | x[k-1])
    % (under the suposition of additive process noise)
    % p_xk_given_xkm1 = @(k, xk, xkm1) p_sys_noise(xk - sys(k, xkm1, 0));
    
    pf.k               = 1;                   % initial iteration number
    pf.Ns              = 200;                 % number of particles
    pf.w               = zeros(pf.Ns, nFrames);     % weights
    pf.particles       = zeros(param.nx, pf.Ns, nFrames); % particles
    
    %----------------- Initial PDF
    %Q0 = sqrt(10);
    pf.gen_x0          = @(x) normrnd(0, Q0); % function for sampling from initial pdf p_x0
    %pf.gen_x0 = @(x) normrnd(x,sqrt([10;2;1;1;10;2;1;1;0.001;0.5]));
    
    %----------------- PDF of observation noise and noise generator function
    nv = param.ny;                                           % size of the vector of observation noise
    sigma_v = sqrt(0.1);
    p_obs_noise   = @(v) mvnpdf(v, zeros(nv,1) , eye(nv)*sigma_v);
    pf.sigma_v = sigma_v;
    pf.p_obs_noise = p_obs_noise;
    
    %----------------- Observation likelihood PDF p(y[k] | x[k])
    % (under the suposition of additive process noise)
    pf.p_yk_given_xk   = @(k, yk, xk) p_obs_noise(yk - param.obs(k, xk, 0));       % function of the observation likelihood PDF p(y[k] | x[k])
    
    % ------ Process noise definition
    sigma_u = sqrt(1); % simplified noise model where Q is diagonal with all variances equal
    pf.sigma_u = sigma_u;
    pf.gen_sys_noise   = @(u) normrnd(0, sigma_u);         
                % function for generating system noise => sample from p_sys_noise (returns column vector)

    %pf.p_x0 = p_x0;                          % initial prior PDF p(x[0])
    %pf.p_xk_given_ xkm1 = p_xk_given_xkm1;   % transition prior PDF p(x[k] | x[k-1])
  end

end
