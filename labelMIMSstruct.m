function [datastruct,hdrtxt,d_div,line_types,type_rows,labelMask,varargout] = ...
  labelMIMSstruct(datastruct,hdrtxt,varargin)
%% function: labelMIMSstruct
% Description:
% Example:
%     [d,hdrtxt,d_div,line_types,type_rows,labelMask,laONmask,tsumask,dchmask] = ...
%         labelMIMSstruct(d,hdrtxt,'ncc',1);
%
% INPUTS ----------------------------------------------------------------
% I:    Input image (M x N)
% varargin - 'PropertyName','PropertyValue'
%   • 'Imasks':   structure (see outputs) => if indicates, F/B segmentation
%   procedure will not be executed and the masks in the structure will be
%   used instead
%   • 'Nmedfilter': integer (best if odd) indicating the width of the median
%   filter used to smooth the data [DEFAULT = 9]
%   • 'nbins':  integer indicating how many histogram bins to use for
%   binning the data. IF nbins == 0, then unique values will be used
%   instead of a number of histogram bins [DEFAULT = 50]
%   • 'plotimgs': logical indicating whether to plot images or not [DEFAULT = 0]
%   • 'saveimgs': logical indicating whether to save images or not [DEFAULT = 0]
%   • 'dataname': string indicating the data's label
%   • 'sqmask': logical indicating whether to square the mask or not
%   [DEFAULT = 1]
%
% OUTPUTS ---------------------------------------------------------------
% Imask:
%
%  Date           Author              E-mail                      Version
%   1 Apr  2015   Amanda Balderrama   amandagbalderrama@gmail.com     1

if ~exist('hdrtxt')
  hdrtxt = [];
end

[datastruct,hdrtxt,d_div,line_types,type_rows] = divMIMSstruct(datastruct,hdrtxt);
ntype = length(line_types);
elnames = fieldnames(rmfield(datastruct,skip_fields));
compdir = loaddirfun;

PropertyNames = lower(varargin(1:2:length(varargin)));
PropertyVal = varargin(2:2:length(varargin));

if strmatch('plotimgs',PropertyNames)
  plotimgs = PropertyVal{strmatch('plotimgs',PropertyNames)};
else
  plotimgs = 0;
end

if strmatch('saveimgs',PropertyNames)
  saveimgs = PropertyVal{strmatch('saveimgs',PropertyNames)};
else
  saveimgs = 0;
end

if strmatch('laoffel',PropertyNames) %% FUTURE IDEA: attempt to find the best element for masking
  laoffel = PropertyVal{strmatch('laoffel',PropertyNames)};
else
  nistel = 'Gd';
  otherel = 'C_13';
  laoffel = cell(ntype,1);
  nistline = cellfun(@(x) ~isempty(strfind(x,'NIST')),line_types);
  laoffel(nistline) = {nistel};
  laoffel(nistline == 0) = {otherel};
end

if strmatch('tsuel',PropertyNames) %% FUTURE IDEA: attempt to find the best element for masking
  tsuel = PropertyVal{strmatch('tsuel',PropertyNames)};
else
  tsuel = 'S_34';
end

if ~ismember(tsuel,elnames)
  tsuel = 'S_33';
end

if strmatch('label',PropertyNames) %% FUTURE IDEA: attempt to find the best element for masking
  labelmethod = PropertyVal{strmatch('label',PropertyNames)};
else
  labelmethod = 'ch';
end

if strmatch('ncc',PropertyNames)
  ncc = PropertyVal{strmatch('ncc',PropertyNames)};
else
  ncc = [];
end

laONmask = cellfun(@(x) false(size(x.Time)),d_div,'uniformoutput',0);
if plotimgs
  spcols = ceil( ntype/2 );
  fh = figure('name',strcat(datastruct.dataset,' Laser ON Mask'),...
    'position',[50,50,350*spcols,475*2]);
else 
  fh = [];
end
for c = 1:ntype
  if iscell(laoffel)
    claoff = elnames{interp_data_elemrng(datastruct,laoffel{c})};
  elseif isstr(laoffel)
    claoff = laoffel;
  else
    claoff = elnames{laoffel(c)};
  end
  if ~isempty(strfind(line_types{c},'HIPP'))
    nmodes = 3;
  else 
    nmodes = 2;
  end
  dmaskraw = getfield(d_div{c},claoff);
  [~,~,d4mask] = auto_thresh(dmaskraw,'auto',[]);
  if plotimgs
    sprow = mod(c,spcols);
    h = subplot(4,spcols,c + floor((c-1)/spcols)*spcols);
  else
    h = [];
  end
  laONmask{c} = bimodalHistMask(d4mask,'Nmedfilter',1,'quantfrac',0,...
    'nbins',50,'plotimgs',plotimgs,'ncc',ncc,'axh',h,'figh',fh,'nmodes',nmodes);
  if plotimgs
    title(sprintf('pixel count dist %s',strrep(line_types{c},'_','\_')))
    subplot(4,spcols,(floor((c-1)/spcols)+1)*spcols+c);
    imagesc(d4mask); hold on;
    tM3 = abs(double(repmat(laONmask{c},[1,1,3]))-1);
    %tM3(:,:,2) = tM3(:,:,2).*0.7;
    im1 = image(tM3);
    set(im1,'AlphaData',tM3(:,:,1).*0.7);
    %imagesc(laONmask{c});
  end
end

if saveimgs
  save_open_figures(datastruct,[],[],'','fcurrent',1);
end


tissueMask = cellfun(@(x) false(size(x.Time)),d_div,'uniformoutput',0);
tsutypes = cellfun(@(x) ~isempty(strfind(x,'Calib')),line_types) + cellfun(@(x) ~isempty(strfind(x,'NIST')),line_types);
if plotimgs
  spcols = sum(tsutypes == 0);
  colcnt = 1;
  fh = figure('name',strcat(datastruct.dataset,' Tissue Mask'),...
    'position',[50,50,350*spcols,475]);
else
  fh = [];
end

for c = find(tsutypes == 0)'
  if iscell(tsuel)
    claoff = elnames{interp_data_elemrng(datastruct,tsuel{c})};
  elseif isstr(tsuel)
    claoff = tsuel;
  else
    claoff = elnames{tsuel(c)};
  end
  if ~isempty(strfind('cugrid',lower(line_types{c})))
    dmaskraw = getfield(d_div{c},'Cu63');
    [~,~,d4mask] = auto_thresh(dmaskraw,'auto',[]);
    if plotimgs
      h = subplot(2,spcols,colcnt);
    else
      h = [];
    end
    tissueMask{c} = bimodalHistMask(d4mask,'Nmedfilter',1,'sqmask',0,...
      'quantfrac',0,'nbins',50,'plotimgs',plotimgs,'fill',0,'axh',h,'figh',fh);
  else
    dmaskraw = getfield(d_div{c},tsuel);
    [~,~,d4mask] = auto_thresh(dmaskraw,'auto',[]);
    if plotimgs
      h = subplot(2,spcols,colcnt);
    else
      h = [];
    end
    tissueMask{c} = bimodalHistMask(d4mask,'Nmedfilter',3,'sqmask',0,...
      'quantfrac',0.01,'nbins',50,'plotimgs',plotimgs,'fill',0,'close',1,...
      'thupperI',0,'ncc',ncc,'axh',h,'figh',fh);
  end
  if plotimgs
    title(sprintf('pixel count dist %s',strrep(line_types{c},'_','\_')))
    subplot(2,spcols,colcnt+spcols);
    imagesc(d4mask); hold on;
    tM3 = abs(double(repmat(tissueMask{c},[1,1,3]))-1);
    %tM3(:,:,2) = tM3(:,:,2).*0.7;
    im1 = image(tM3);
    set(im1,'AlphaData',tM3(:,:,1).*0.7);
    colcnt = colcnt + 1;
  end
end

if saveimgs
  save_open_figures(datastruct,[],[],'','fcurrent',1);
end

labelMask = cellfun(@(x,y) double(x) + double(y).*2,laONmask,tissueMask,'uniformoutput',0);

ncol = 4;
[M,N] = size(datastruct.Time);
chlabelMask = nan(M,N,ncol);
for c = 1:ntype
  %dlaONall(type_rows{c},1:typesize(c,2)) = labelMask{c}; % combines all maps together
  [tM,tN] = size(laONmask{c});
  if ~any(laONmask{c}(:))
    chlabelMask(type_rows{c},1:tN,1) = ones(tM,tN);
  else
    halfcol  = round(tN/2);
    LA0beg = nan(tM,tN);
    i0beg = find(laONmask{c}(:,1:halfcol)==0);
    LA0beg(i0beg) = 1;
    LA0end = nan(tM,tN);
    iend = find(laONmask{c}(:,(halfcol+1):end)==0) + halfcol*tM;
    LA0end(iend) = 1;
    chlabelMask(type_rows{c},1:tN,1) = LA0beg;
    chlabelMask(type_rows{c},1:tN,2) = LA0end;
  end
  
  bkgndind = double(labelMask{c} == 1);
  bkgndind(bkgndind == 0) = nan;
  chlabelMask(type_rows{c},1:tN,3) = bkgndind;
  
  tsuind = double(labelMask{c} == 3);
  tsuind(tsuind == 0) = nan;
  chlabelMask(type_rows{c},1:tN,4) = tsuind;
end


varargout{1} = laONmask;
varargout{2} = tissueMask;
varargout{3} = chlabelMask;
end
