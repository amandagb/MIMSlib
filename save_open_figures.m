function save_open_figures(d,xlimits,ylimits,figure_tag,varargin)
%% Script: save_open_figures(d,xlimits,ylimits,figure_tag,varargin)
%         save_open_figures(mtngpath,[],[],'','fmt','png');
% Description:
% Example:
% Required Functions:
% INPUTS ----------------------------------------------------------------
% d:  this variable is used to define "datastr" and "folderpath". It can be
% of the following data types:
%   • structure with a field 'dataset' containing the parent folder of the
%   dataset corresponding to the image
%     => folderpath = metallomic image parent directory, datastr = d.dataset;
%   • string of a parent directory folder where the images should be saved
%   (isdir must be true for the images to be saved into this path)
%     => folderpath = d, datastr = '';
%   • describptive string identifying the the data that are plotted
%     => folderpath = current matlab directory, datastr = d;
%   • if none, the data is identifyied by an emptry string ('')
%     => folderpath = current matlab directory, datastr = '';
% xlimits: [MIN,MAX] values which will be used to adjust the x-limits on
%   the axis of each of the plots (leave empty for no adjustments)
% ylimits: [MIN,MAX] values which will be used to adjust the x-limits on
%   the axis of each of the plots (leave empty for no adjustments)
% figure_tag: descriptive string to add to the file name
% varargin - 'PropertyName','PropertyValue'
%   'fcurrent': binary indicating whether to save only the current figure
%   (1) or all the figures (0) [DEFAULT = 0]
%   'fmt':  string indicating the formate of the image to be saved (see
%   HELP saveas for options) [DEFAULT = 'tif']
%   'parent': string indicating parent directory data corresponds to
%       • 'la': 'Laser Albation\' folder
%       • 'ivis': 'EB_Experiments\IVIS_EBFluor\' folder
%       • 'mtng': '\Latex Citation Summary\' folder
%   'folder': string indicating the folder within the parent directory that
%       the image should be saved into
%       • 'noborder': logical indicating whether to remove white space
%       around the image or not [DEFAULT = 0]
% OUTPUTS ---------------------------------------------------------------
% figure saved with the following path designation:
%   imgnm = sprintf('%s\\%s %s %s.%s',...
%       folderpath, datastr, figure_tag, fign{j}(1:cutoff{j}-2),fmt);
%
%  Date           Author            E-mail                      Version
%  17 Sept 2013   Amanda Gaudreau   amanda.gaudreau@gmail.com     1

PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

cpaths = loaddirfun;

%% Initialize variables that deal with scaling the image (converting pixels --> cm)
if strmatch('fcurrent',PropertyNames)
  cf = PropertyVal{strmatch('fcurrent',PropertyNames)};
else
  cf = 0; %microns
end

if strmatch('fmt',PropertyNames)
  fmt = PropertyVal{strmatch('fmt',PropertyNames)};
else
  fmt = 'png'; %CHANGED TO PNG from 'tif' 23 Feb 2015
end

if strmatch('parent',PropertyNames)
  pdir = PropertyVal{strmatch('parent',PropertyNames)};
else
  pdir = 'la'; %microns
end

if strmatch('folder',PropertyNames)
  fldstr = PropertyVal{strmatch('folder',PropertyNames)};
else
  fldstr = ''; %
end

if strmatch('noborder',PropertyNames)
  rmborder = PropertyVal{strmatch('noborder',PropertyNames)};
else
  rmborder = 0; %
end

if ~exist('figure_tag')
  figure_tag = '';
end

%dockf off all
if isstruct(d)
  if isfield(d,'dataset')
    datastr = d.dataset;
  else
    datastr = '';
  end
  
  if strmatch(pdir,'tex')
    folderpath = strcat(cpaths.mtngpath,fldstr);
    fmt = 'png';
  elseif isfield(d,'origseq') || ~isempty(strmatch(pdir,'ivis'))
    if ~isempty(fldstr)
      fpath = strcat(cpaths.ivisfldr,datastr,'\',fldstr);
    else
      fldstr = 'Images';
      fpath = strcat(cpaths.ivisfldr,datastr,'\Images');
    end
    
    if ~isdir(fpath)
      mkdir(cpaths.ivisfldr,fldstr)
    end
    
    folderpath = fpath;
  elseif  strmatch(pdir,'la')
    folderpath = strcat(cpaths.lafldr,datastr,'\Images');
    if ~isdir(folderpath)
      mkdir(strcat(cpaths.lafldr,datastr),'Images')
    end
  end
  %   folderpath = strcat('D:\My Documents\MADLab Data\Laser Ablation\',datastr,'\Images');
elseif isdir(d)
  datastr = '';
  if strmatch('\',d(end))
    d = d(1:end-1);
  end
  folderpath = d;
elseif isstr(d)
  datastr = d;
  folderpath = pwd;
else
  datastr = '';
  folderpath = pwd;
end

if cf
  h = gcf;
else
  h = get(0,'children');
end

set(h,'PaperPositionMode','auto');
fign = get(h,'Name');
if ~iscell(fign)
  fign = {fign};
end
cutoff = strfind(fign,'Heat');

if all(cellfun(@(x) (isempty(x)),cutoff))
  cutoff = cell(length(cutoff),1);
  cutoff = cellfun(@(x) (length(x)+2),fign);
  cutoff = num2cell(cutoff);
end
% empty_cutoff = cellfun(@(x) isempty(x),cutoff);
% cutoff{emtpy_cutoff} = length(fign{empty_cutoff});
j = 1;
for i=h'
  %set(h(j),'position',[1     1   970   400]);
  a = get(h(j),'currentaxes');
  if exist('xlimits') && ~isempty(xlimits)
    set(a,'xlim',xlimits);
  end
  
  if exist('ylimits') && ~isempty(ylimits)
    set(a,'ylim',ylimits);
  end
  
  if isempty(figure_tag) && isempty(datastr)
    if isempty(fign{j}(1:cutoff{j}-2))
      imgnm = sprintf('%s\\Img%d',...
        folderpath,i);
    else
      imgnm = sprintf('%s\\%s',...
        folderpath, fign{j}(1:cutoff{j}-2));
    end
  elseif isempty(figure_tag)
    if isempty(fign{j}(1:cutoff{j}-2))
      imgnm = sprintf('%s\\%s',...
        folderpath, datastr);
    else
      imgnm = sprintf('%s\\%s_%s',...
        folderpath, datastr, fign{j}(1:cutoff{j}-2));
    end
  elseif isempty(datastr)
    if isempty(fign{j}(1:cutoff{j}-2))
      imgnm = sprintf('%s\\%s',...
        folderpath, figure_tag);
    else
      imgnm = sprintf('%s\\%s_%s',...
        folderpath, figure_tag, fign{j}(1:cutoff{j}-2));
    end
  else
    if isempty(fign{j}(1:cutoff{j}-2))
      imgnm = sprintf('%s\\%s_%s',...
        folderpath, datastr, figure_tag);
    else
      imgnm = sprintf('%s\\%s_%s_%s',...
        folderpath, datastr, figure_tag, fign{j}(1:cutoff{j}-2));
    end
  end
  
  img = getframe(h(j));
  if strmatch(fmt,'fig')
    saveas(h(j), imgnm, fmt);
  elseif rmborder
    export_fig(imgnm,strcat('-',fmt));
  else
    imwrite(img.cdata, [imgnm, '.',fmt]);
  end
  
  j = j + 1;
end