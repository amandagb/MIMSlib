function outputs = mfile_header(inputs,varargin)
%% Script:
% Description:
% Example:
% Required Functions:
% -------------------------------------------------------------------------
% INPUTS
% -------------------------------------------------------------------------
% data:  data structure which is output from 'readfiles.m' script. Should
%       contain fields which have <# lines> x <# samples> matrix for each
%       isotype analyzed
% varargin - 'PropertyName','PropertyValue'
%   'color_res':  numeric specifying color resolution of the heatmap
%   (maximum 256) [DEFAULT = 256]
%   'thresh_data':  1 x 2 matrix with minimum threshold value as first
%   element and maximum threshold value as second
%   'lines': matrix containing line numbers user would like to plot
%   [DEFAULT: ALL LINES]
%   'samples': matrix containing samples (S) user would like to plot [DEFAULT:
%   ALL SAMPLES]
%   'plot_type': string - 'log','linear','both' [DEFUALT: BOTH]
%   'elim_blanks': 1 or 0 - eliminates any lines which have NAN data
% -------------------------------------------------------------------------
% OUTPUTS
% -------------------------------------------------------------------------% 
%
%  Date           Author            E-mail                      Version
%  19 July 2011   Amanda Gaudreau   amanda.gaudreau@gmail.com     1


PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

if strmatch('nbins',PropertyNames)
  nbins = PropertyVal{strmatch('nbins',PropertyNames)};
else
  nbins = 1e3;
end

[elem_ana_range,fldnm] = interp_data_elemrng(data,elem_ana_range);
for f = elem_ana_range
  if isstruct(data)
    d = getfield(data,fldnm{f});
  else
    d = data;
  end
  if isempty(strmatch(fldnm{f},skip_fields)) && ...
      f <= length(fldnm) %isempty(strmatch(fldnm{f},'Time')) && isempty(strmatch(fldnm{f},'line_datevec')) && f <= length(fldnm)
    %isempty(strmatch(fldnm{f},{'line_datevec','dataset'}))
    % function
  end
end

%% Adding data to an existing structure
dout = setfield(dout,fldnm{f},mu_dout);

%% Initializing file
data = readfiles('comp','pc','data','120710Brain');

%% LATEX STRINGS IN LEGEND
h = legend('holder1','holder2');
h1 = findobj(get(h,'Children'),'String','holder1');
set(h1,'String','$frac{d\mu_{MA}}{ds}$','Interpreter','latex')
set(h1,'String','$\frac{d\mu_{MA}}{ds}$','Interpreter','latex')
h2 = findobj(get(h,'Children'),'String','holder2');
set(h2,'String','$\left(\displaystyle\frac{d\mu_{MA}}{ds}\right)_{MA}$','Interpreter','latex')

%% Saving all open figures

h = get(0,'children');
fign = get(h,'Name');

for i=h'
  saveas(h(i), strcat('D:\My Documents\Dropbox\MADLab Research\Data\012712_OEScalibration\',...
    fign{i}(1:6),' thresh LINheatmap '), 'jpg');%fldnm{i},' LINheatmap'),'jpg');
  saveas(h(i), strcat('D:\My Documents\Dropbox\MADLab Research\Meeting Notes\LatexNotes\',...
    fign{i}(1:6),' thresh LINheatmap '), 'jpg');%fldnm{i},' LINheatmap'),'jpg');
end

h = get(0,'children');
set(h,'PaperPositionMode','auto');
fign = get(h,'Name');
j = 1;
for i=h'
  %set(h,'colormap',gray);
  a = get(h(j),'currentaxes');
  set(a,'xlim',[0,30]);
  imgnm = sprintf('D:\\My Documents\\MADLab Data\\Laser Ablation\\%s\\Images\\%s spatial %s',...
    d.dataset, d.dataset, fign{j}(1:6));
  saveas(h(j), imgnm, 'jpg');
  j = j + 1;
end

%% Systematically removes rows and/or columns from a structure
d = rmfield(d,'dataset');
d = structfun(@(x) (x(nonzeroind,:)),d,'uniformoutput',false);
d = structfun(@(x) (rot90(rot90(x))),d,'uniformoutput',false);
d.dataset = din.dataset;

%% General rules for using fir1 filtering
n = 10; % filter order, must be even
%freqz(data); % Allows you to visualize the spectrum of the data
data = [data(1:n/2-1);data;data(end-n/2+1:end)];

% Low pass filtering
filttype = 'lpf'; filtparam = 0.07; % bandwidth of the filter
h = fir1(n,filtparam);
filtdata = filter(h,1,data);
filtdata = filtdata(n:end);

%% Filling an empty cell
klog_cell = cell(Abins+1,1);
[klog_cell{1:end}] = deal(zeros(size(A,1),size(A,2)));

%% Plotting tools
set(gca,'LineStyleOrder',{'-*',':','o'})
markers = {'o','*','s','d','^','x','h'};

