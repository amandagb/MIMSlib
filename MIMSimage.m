function [mask,varargout] = MIMSimage(data,elem_ana_range,varargin)
%% Script:  [mask,varargout] = MIMSimage(data,elem_ana_range,varargin)
% Description:
% Example:
% INPUTS ----------------------------------------------------------------
%
% OUTPUTS ---------------------------------------------------------------
% mask : l x s logical tells if each pixel is part of the background ( ==
% 1) or the tissue (== 0)
%
%  Date           Author            E-mail                      Version
%                 Tim Connelly      timpconnelly@gmail.com        0
%  25 July 2011   Amanda Gaudreau   amanda.gaudreau@gmail.com     1

dir_name = 'D:\My Documents\My Dropbox\MADLab Research\Data\12-07-20 TgBrain';
contig_pnts = 100;
% Isoptope to process.
elem = 'Zn66';

% Character string to be used in titles and filenames
dir_div = strfind(dir_name,'\');
set_name = dir_name(dir_div(end)+1:end);

% Isotope to use for mask. Cu or Zn work well.
mask_iso = 'Zn66';

% Tail end probability to threshold range. Lower -> less threshold.
%                                          High  -> more threshold.
p = 0.005;

% Angle to rotate image. You will want to use multiples of 90. Non-90
% degree multples can be used but will require heavy interpolation.
angle = 0;

% Turn Image output On/Off. 0 -> no image ouput.
print_control = 0;

% Controls whether to shift rows, not generally needed now.
shift_control = 0;

% Control Drift Adjustment
drift_adjust = 0;

% change 0/1 to Flip final image Left-Right and Up-Down
LRflip = 0;
UDflip = 0;

% Segmentation parameter. 1 is default. Lower gets more background and
% higher captures less. Very sensitive to small changes!
eta = 1;

% Ablation Parameters
diameter = 25;  %um
line_separation = 2*diameter; %um
rep_time = 0.258; %seconds
scan_speed = 100; %um/seconds
deltax = rep_time*scan_speed;
deltay = line_separation;
y_ratio = round(deltay/deltax);

%% Read in raw data and display in log_10
if isstruct(data)
  fldnm = fieldnames(data);
  v = elem_ana_range;
  elem_ana_range = [];
  if ~isnumeric(v)
    if iscell(v)
      for i = 1:length(v)
        if isempty(strmatch(lower(v{i}),lower(fldnm)))
          disp(['The element you have indicated, (', v{i},'), does not exist in the dataset']);
        else
          elem_ana_range = [elem_ana_range,strmatch(lower(v{i}),lower(fldnm))];
        end
      end
    else
      elem_ana_range = [elem_ana_range,strmatch(lower(v),lower(fldnm))];
    end
  elseif isempty(v)
    elem_ana_range = 1:length(fldnm);
  end
else
  fldnm{1} = v;
  elem_ana_range = 1;
  d = data;
end

for f = elem_ana_range
  elem = fldnm{f};
  if ~exist('d','var')
    d = getfield(data,fldnm{f});
  end
end
time = data.Time;
% A = circshift(A,[-133 -46]);  %uncomment this line and fiddle with
%                       numbers if there is a large offset in the image.\
[l,s] = size(d);

%% Normalization and mean subtraction
% Calculates instrument drift and corrects by either subtracting off
% difference from mean or normalizing while keeping range the same.

% if(drift_adjust)
%     sub = 1;
%     div = 2;
%     range = [1 25];
%     [A, track] = bg_adjust(A,range,sub);
%
%     figure
%     plot(track)
%     title('Instrument Drift')
%     xlabel('row')
%     box off
% end


%% Segmentation and Threshold: Assumes both background and overall image
% distributations can be approximated by gaussians
bg_lines = 1:25;
bg_mean = mean2(d(bg_lines,:));      % Estimate BG mean from 1-25 rows
bg_var = std2(d(bg_lines,:));   % Estimate BG std deviation

im_mean = mean(d(:));       % Image/FG mean
im_sdev = std(d(:));        % Image/FG std dev.

[N, bins] = hist(d(:),1e4); % histogram data
N = N/sum(N);               %normalize histogram to a pdf
PDF = cumsum(N);            % pdf to PDF

[Y,I] = sort(PDF);          % sort PDF and find where to threshold
clip_ind = find(Y>(1-p));

if (~isempty(clip_ind))     %threshold image
  clip_val = bins(I(clip_ind(1)));
  locs = find(d>clip_val);
  d(locs) = clip_val;
end

% Segment image / set min value

% segmentation parameter
gamma = ((bg_mean+im_mean)/2)+((im_sdev^2)*log(eta))/(im_mean - bg_mean);

param.bg_mean = bg_mean;
param.bg_var = bg_var;
param.im_mean = im_mean;
param.im_var = im_sdev;
param.thresh = clip_val;
param.contig_pnts = contig_pnts;
param.gamma = gamma;
varargout{1} = param;

% use saved mask or
% if strcmp(elem, mask_iso)
mask = d < gamma;
temp = bwareaopen(~mask, contig_pnts);
mask = ~temp;
low_locs = find(mask);

% Display background-foreground image and output
%   figure
%   imagesc(mask);
%   colormap(gray)
%   title([set_name ' Mask']);
%   print( '-dpng',  [set_name '_mask.png'])
% else
%   mask1 = d < gamma;
%   load mask
%   low_locs = find(mask|mask1);
%
%   % Display background-foreground image
%   figure
%   imagesc(mask);
%   colormap(gray)
%   title([set_name ' Mask']);
% end

d(low_locs) = gamma;


%% Resize and Rotate

di = myinterp2(d,y_ratio); % myinterp2() is a custom function
% Ai = Ai';   %re-orient matrix to its initial form

di = imrotate(di, angle, 'nearest', 'loose');
d = imrotate(d, angle, 'nearest', 'loose');


%% alignment
%
% If alignment is on a window asking to crop the image will appear.
% Select the desired area and proceed.

if (shift_control)
  aligned = align_image(di,30,gamma);
  figure
  
  ds = imcrop((aligned/max(aligned(:))));
  figure
  imshow(ds,[])
  colormap jet
  di = ds;
end

%% Image display and print

fname1 = strcat(elem,set_name, '_log.tif');
fname3 = strcat(elem,set_name, '.tif');

% Get element symbol and isotope
if (strcmp(elem, 'P31'))
  esymbol = elem(1);
  isoval = elem(2:length(elem));
else
  esymbol = elem(1:2);
  isoval = elem(3:length(elem));
end

% Flip image as needed.

if (LRflip)
  di = fliplr(di);
  d = fliplr(d);
end

if (UDflip)
  di = flipud(di);
  d = flipud(d);
end

% figure('Name',elem,'color','w','position',[360,502,560,420]);
% imagesc(log10(d));
% title([set_name '^{' isoval '}' esymbol ' log_{10} Scale']);
% colormap(jet)
% colorbar
% if (print_control)
%   print( '-dtiffnocompression', fname1)
% end


fig_handle = figure('Name',elem,'color','w','position',[360,502,560,420]);
imagesc(d);
title([set_name ' ^{' isoval '}' esymbol ' Linear Scale'],'fontsize',12);
xlabel('Sample #','fontsize',12);
ylabel('Line #','fontsize',12);
colormap jet
colorbar
if(print_control)
  print( '-dtiffnocompression',  fname3)
end
dockf on all
end