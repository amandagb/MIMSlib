function [param, utilities] = HSVtracking(filename,varargin)

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
%
% OUTPUTS ---------------------------------------------------------------
%
%  Date           Author              E-mail                      Version
%  8 July 2016   Amanda Balderrama  amandagbalderrama@gmail.com     0
%    from HSVtrackingDebeug v1. Example code deleted

cdir = loaddirfun;
% filename = strcat(cdir.fpath,'Impacter Nose Tracking\HSV69_B1_C1 Comparison\blastHSV69_B1_C1.avi');

%% Initialize variables
% The |trackSingleObject| function includes nested helper functions. The
% following top-level variables are used to transfer the data between the
% nested functions.
if ~exist('filename')
  filename = strcat(cdir.dboxpath(1:end-7),'Desktop\HSVtests\20130605_0910_HSV69_B1_C1_convert.avi');
elseif isempty(filename)
  filename = strcat(cdir.dboxpath(1:end-7),'Desktop\HSVtests\20130605_0910_HSV69_B1_C1_convert.avi');
end

%----------------- initialize variables
frame            = [];  % A video frame
detectedLocation = [];  % The detected location
trackedLocation  = [];  % The tracked location
label            = '';  % Label for the ball
utilities        = [];  % Utilities used to process the video
param = getDefaultParameters();  % get Kalman configuration that works well

%----------------- initialize the utilities variable
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

%----------------- check for other specified varaibles
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


% showDetections();
% showTrajectory();

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
    
    if ~isTrackInitialized
      if isObjectDetected
        % Initialize a track by creating a Kalman filter when the ball is
        % detected for the first time.
        kalmanFilter = configureKalmanFilter(param.motionModel, ...
          detectedLocation, param.initialEstimateError, ...
          param.motionNoise, param.measurementNoise);
        
        isTrackInitialized = true;
        trackedLocation = correct(kalmanFilter, detectedLocation);
        label = 'Initial';
      else
        trackedLocation = [];
        label = '';
      end
      
    else
      % Use the Kalman filter to track the ball.
      if isObjectDetected % The ball was detected.
        % Reduce the measurement noise by calling predict followed by
        % correct.
        %predict(kalmanFilter)
        trackedLocation = correct(kalmanFilter, detectedLocation);
        label = 'Corrected';
      else % The ball was missing.
        % Predict the ball's location.
        trackedLocation = predict(kalmanFilter);
        label = 'Predicted';
      end
    end
    
    annotateTrackedObject();
    
    %figure; plot(utilities.accumulatedDetections(1:utilities.detectind,1),utilities.accumulatedDetections(1:utilities.detectind,2),'.')
    %hold on;
    %plot(utilities.accumulatedTrackings(:,1),utilities.accumulatedTrackings(:,2),'.r')
    
  end % per frame evaluation (detection & tracking)
end % while
showTrajectory();

%%
% Get default parameters for creating Kalman filter and for segmenting the
% ball.
  function param = getDefaultParameters
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
    utilities.foregroundMinIntensity = objOut.foregroundMinIntensity;
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
