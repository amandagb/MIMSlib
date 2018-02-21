function [param, utilities] = HSVtrackingPF(filename,varargin)

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
%   > stateTrans:   matrix 
%
% OUTPUTS ---------------------------------------------------------------
%
%  Date           Author              E-mail                      Version
%  8 July 2016   Amanda Balderrama  amandagbalderrama@gmail.com     0
%    from HSVtrackingDebeug v1. Example code deleted
% 30 Aug  2016   Amanda Balderrama  amandagbalderrama@gmail.com     1
%    Cleaner version with particle filter parameters in a new subfunction

cdir = loaddirfun;
% filename = strcat(cdir.fpath,'Impacter Nose Tracking\HSV69_B1_C1 Comparison\blastHSV69_B1_C1.avi');

if ~exist('filename')
  filename = strcat(cdir.dboxpath(1:end-7),'Desktop\HSVtests\20130605_0910_HSV69_B1_C1_convert.avi');
elseif isempty(filename)
  filename = strcat(cdir.dboxpath(1:end-7),'Desktop\HSVtests\20130605_0910_HSV69_B1_C1_convert.avi');
end

%% Initialize variables
startvarin = 1;
PropertyNames = lower(varargin(startvarin:2:length(varargin)));
PropertyVal = varargin(startvarin+1:2:length(varargin));

if strmatch('objectdet',PropertyNames)
  detections = PropertyVal{strmatch('objectdet',PropertyNames)};
  utilities.accumulatedDetections = detections;
  nonNans = find(sum(isnan(utilities.accumulatedDetections),2) == 0);
  utilities.detectind = nonNans(end);
else
  detections = [];
end

if strmatch('track',PropertyNames)
  track = PropertyVal{strmatch('track',PropertyNames)};
else 
  track = 1;
end

if strmatch('stateTrans',PropertyNames)
  F = PropertyVal{strmatch('stateTrans',PropertyNames)};
else 
  F = [];
end

%----------------- initialize variables
frame            = [];  % A video frame
detectedLocation = [];  % The detected location
trackedLocation  = [];  % The tracked location
label            = '';  % Label for the ball
utilities = [];
param = getDefaultParameters(F);  % get Kalman configuration that works well

%----------------- initialize the utilities variable
utilities.frameidx = 0;
utilities.videoReader = vision.VideoFileReader(filename);%,'ImageColorSpace','Intensity');
% utilities.videoPlayer = vision.VideoPlayer('Position', [100,100,500,400]);
% utilities.foregroundDetector = vision.ForegroundDetector(...
%   'NumTrainingFrames', 10, 'InitialVariance', param.segmentationThreshold);
% utilities.blobAnalyzer = vision.BlobAnalysis('AreaOutputPort', false, ...
%   'MinimumBlobArea', param.minBlobSize, 'CentroidOutputPort', true);

utilities.accumulatedImage      = 0;
utilities.trackind = 0;
utilities.accumulatedTrackings  = nan(param.frameRange(2), 2);
utilities.detectind = 0;% param.frameRange(1)-1;
utilities.accumulatedDetections = nan(param.frameRange(2),2);
utilities.foregroundMinIntensity = Inf;

if strmatch('utilities',PropertyNames)
  utilities = PropertyVal{strmatch('utilities',PropertyNames)};
end

T = (param.frameRange(2) - param.frameRange(1) + 1);  %# time steps
pf = setDefaultPF(T,param);

%----------------- Separate memory space
y = zeros(param.ny,T);

%----------------- Separate memory
xh = zeros(param.nx, T); %xh(:,1) = xh0;
yh = zeros(param.ny, T); %yh(:,1) = obs(1, xh0, 0);

isTrackInitialized = false;
dind = 0;
while ~isDone(utilities.videoReader)
  frame = readFrame();
  
  if utilities.frameidx >= param.frameRange(1) && utilities.frameidx <= param.frameRange(2)
    if ~isempty(detections)
      dind = dind + 1;
      detectedLocation = detections(dind,:);
      isObjectDetected = 1-sum(isnan(detectedLocation));
      utilities.foregroundMask = zeros(size(frame,1), size(frame,2));
    else
      [detectedLocation, isObjectDetected] = detectObject(frame);
    end
    
    annotateTrackedObject();
    %----------------- Fill observations vector
    y(:,utilities.detectind) = ...
          [utilities.accumulatedDetections(utilities.detectind,:)';...
          utilities.foregroundMinIntensity;...
          sum(utilities.foregroundMask(:))];
  end
end

if ~isfield(utilities,'observations')
  utilities.observations = y;
else
  y = utilities.observations;
end

if track
  k = 1;
  xh0 = [y(1,k),0,0,0,y(2,k),0,0,0,y(3,k),y(4,k)]';
  xh(:,k) = xh0;
  yh(:,k) = param.obs(k,xh0,0);
  utilities.accumulatedTrackings(k,:) = yh(1:2,k);
  pf.gen_x0 = @(x) normrnd(xh0,sqrt([10;2;1;1;10;2;1;1;0.001;0.5]));
  
  kalmanFilter = configureKalmanFilter(param.motionModel, ...
    y(:,k), param.initialEstimateError, ...
    param.motionNoise, param.measurementNoise);
  
  for k = 2:T
    pf.k = k;
    [xh(:,k), pf] = particle_filter(param.sys, y(:,k), pf, 'systematic_resampling');  
    yh(:,k) = param.obs(k,xh(:,k),0);
    utilities.accumulatedTrackings(k,:) = yh(1:2,k);
    
    if ~isTrackInitialized
     if ~any(isnan(y(:,k)))
        % Initialize a track by creating a Kalman filter when the ball is
        % detected for the first time.
        %kalmanFilter = configureKalmanFilter(param.motionModel, ...
        %  y(:,k), param.initialEstimateError, ...
        %  param.motionNoise, param.measurementNoise);
       
       isTrackInitialized = true;
       trackedLocation = correct(kalmanFilter, y(:,k));
       label = 'Initial';
    else
      trackedLocation = [];
       label = '';
     end
    
    else
     % Use the Kalman filter to track the ball.
     if ~any(isnan(y(:,k)))%isObjectDetected % The ball was detected.
       % Reduce the measurement noise by calling predict followed by
       % correct.
       %predict(kalmanFilter)
       trackedLocation = correct(kalmanFilter, y(:,k));
       label = 'Corrected';
     else % The ball was missing.
       % Predict the ball's location.
       trackedLocation = predict(kalmanFilter);
       label = 'Predicted';
     end
    end
  end
  utilities.states = xh;
end

showTrajectory();

%%
% Get default parameters for creating Kalman filter and for segmenting the
% ball.
  function param = getDefaultParameters(F)
    param.motionModel           = 'ConstantAcceleration';
    param.initialLocation       = 'Same as first detection';
    param.initialEstimateError  = 1E5 * ones(1, 3);
    param.motionNoise           = [2500, 100, 10];
    param.measurementNoise      = 0.5;
    param.segmentationThreshold = 0.05;
    param.movementThreshold = 1;
    param.minBlobSize = 5;
    param.frameRange = [500,3e3]; % starts tracking procedure after 5 ms - 30 ms. entire duration of HSV experiments is typically 250 msec or 25000 frames
    param.detectionMethod       = 'kmeans';
    param.frameRate = 1e5; % frame rate of the video in frames/sec
    param.nx = 10;
    if isempty(F)
      F = [1,1,0,0,0,0,0,0,0,0;...[1,1,1/2,1/6,0,0,0,0,0,0;...
        0,1,1,0,0,0,0,0,0,0;...0,1,1,1/2,0,0,0,0,0,0;...
        0,0,1,1,0,0,0,0,0,0;...0,0,1,1,0,0,0,0,0,0;...
        0,0,0,1,0,0,0,0,0,0;...
        0,0,0,0,1,1,0,0,0,0;...0,0,0,0,1,1,1/2,1/6,0,0;...
        0,0,0,0,0,1,1,0,0,0;...0,0,0,0,0,1,1,1/2,0,0;...
        0,0,0,0,0,0,1,1,0,0;...
        0,0,0,0,0,0,0,1,0,0;...
        0,0,0,0,0,0,0,0,1,0;...
        0,0,0,0,0,0,0,0,0,1];
    end
    param.stateMTX = F;
    param.sys = @(k, xkm1, uk) ...
      F * xkm1 + ...
      ones(param.nx,1) * uk; % attempt at creating a "constant jerk" model with mean intensity tracking and side of block tracking 
     
    param.ny = 4;
    H = [1,0,0,0,0,0,0,0,0,0;...
       0,0,0,0,1,0,0,0,0,0;...
       0,0,0,0,0,0,0,0,1,0;...
       0,0,0,0,0,0,0,0,0,1];
    param.obsMTX = H;
    param.obs = @(k,xk,vk) ...
      H * xk + ones(param.ny,1)*vk;
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


%%
% Read the next video frame from the video file.
  function frame = readFrame()
    frame = step(utilities.videoReader);
    utilities.frameidx = utilities.frameidx + 1;
  end

%% Detect and annotate the ball in the video.
  function showDetections()
    param = getDefaultParameters();
    utilities = createUtilities(param);
    trackedLocation = [];
    
    while ~isDone(utilities.videoReader)
      frame = readFrame();
      if utilities.frameidx >= param.frameRange(1) && utilities.frameidx <= param.frameRange(2)
        detectedLocation = detectObject(frame);
        
        % -------- plotting tools useful for debugging
        %figure; plot(utilities.accumulatedDetections(1:utilities.detectind,1),utilities.accumulatedDetections(1:utilities.detectind,2),'.')
        %figure; plot3(utilities.accumulatedDetections(:,1),utilities.accumulatedDetections(:,2),1:3000,'.')
        %unique(utilities.accumulatedDetections(1:utilities.detectind,1))
        
        % Show the detection result for the current video frame.
        annotateTrackedObject();
        
        % To highlight the effects of the measurement noise, show the detection
        % results for the 40th frame in a separate figure.
        if 0 %mod(utilities.frameidx,100) == 0
          combinedImage = max(repmat(utilities.foregroundMask, [1,1,size(frame,3)]), frame);
          figure, imshow(combinedImage);
        end
      end
    end % while
    
    % Close the window which was used to show individual video frame.
    uiscopes.close('All');
  end

%%
% Detect the ball in the current video frame.
  function [detection, isObjectDetected] = detectObject(frame)
    if ~utilities.detectind
      [initXY,initOut] = detectHSVobject(frame,'segmethod','thresh');
      [detection,objOut] = detectHSVobject(frame,'segmethod',param.detectionMethod,...
        'searchsize',10,'bounding',initOut.BoundingBox);
      hasObjectMoved = false;
    else
      SA = 15;
      [detection,objOut] = detectHSVobject(frame,'segmethod',param.detectionMethod,...
        'intensity',utilities.foregroundMinIntensity,...
        'seed',utilities.accumulatedDetections(utilities.detectind,:),...
        'search',SA);
      maskAreaDiff = sum(objOut.foregroundMask(:)) - sum(utilities.foregroundMask(:));
      
      while (maskAreaDiff > sum(utilities.foregroundMask(:))*.5) || ...
          objOut.foregroundMinIntensity < utilities.foregroundMinIntensity/2
        SA = SA - 2;
        [detection,objOut] = detectHSVobject(frame,'segmethod','CV',...
          'intensity',utilities.foregroundMinIntensity,...
          'seed',utilities.accumulatedDetections(utilities.detectind,:),...
          'search',SA);
        maskAreaDiff = sum(objOut.foregroundMask(:)) - sum(utilities.foregroundMask(:));
      end
    end
    utilities.foregroundMask = objOut.foregroundMask;
    utilities.foregroundMinIntensity = objOut.foregroundMeanIntensity;
    if isempty(detection)
      isObjectDetected = false;
    else
      % To simplify the tracking process, only use the first detected object.
      detection = detection(1, :);
      isObjectDetected = true;
    end
  end

%%
% Show the current detection and tracking results.
  function annotateTrackedObject()
    accumulateResults();
    % Combine the foreground mask with the current video frame in order to
    % show the detection result.
    combinedImage = max(repmat(utilities.foregroundMask, [1,1,size(frame,3)]), frame);
    
    if ~isempty(trackedLocation)
      shape = 'circle';
      region = trackedLocation;
      region(:, 3) = 5;
      combinedImage = insertObjectAnnotation(combinedImage, shape, ...
        region, {label}, 'Color', 'red');
    end
    %step(utilities.videoPlayer, combinedImage);
  end

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
  end

%%
% Accumulate video frames, detected locations, and tracked locations to
% show the trajectory of the ball.
  function accumulateResults()
    utilities.accumulatedImage      = max(utilities.accumulatedImage, frame);
    if ~isempty(detectedLocation)
      utilities.detectind = utilities.detectind + 1;
      utilities.accumulatedDetections(utilities.detectind,:) = detectedLocation;
    end
    
    if ~isempty(trackedLocation)
      utilities.trackind = utilities.trackind + 1;
      utilities.accumulatedTrackings(utilities.trackind,:) = trackedLocation;
    end
  end

%%
% Create utilities for reading video, detecting moving objects, and
% displaying the results.
  function utilities = createUtilities(param)
    % Create System objects for reading video, displaying video, extracting
    % foreground, and analyzing connected components.
    utilities.frameidx = 0;
    utilities.videoReader = vision.VideoFileReader(filename);%,'ImageColorSpace','Intensity');
    utilities.videoPlayer = vision.VideoPlayer('Position', [100,100,500,400]);
    utilities.foregroundDetector = vision.ForegroundDetector(...
      'NumTrainingFrames', 10, 'InitialVariance', param.segmentationThreshold);
    utilities.blobAnalyzer = vision.BlobAnalysis('AreaOutputPort', false, ...
      'MinimumBlobArea', param.minBlobSize, 'CentroidOutputPort', true);
    
    utilities.accumulatedImage      = 0;
    utilities.trackind = 0;
    utilities.accumulatedTrackings  = nan(param.frameRange(2), 2);
    utilities.detectind = 0;% param.frameRange(1)-1;
    utilities.accumulatedDetections = nan(param.frameRange(2),2);
    utilities.foregroundMinIntensity = Inf;
  end

end
