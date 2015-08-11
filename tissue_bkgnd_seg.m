function [boundary,varargout] = tissue_bkgnd_seg(data,elem_ana_range,varargin);

%% Script:   boundary = tissue_bkgnd_seg(data,elem_ana_range,varargin);
% Description:
% Required Functions: interp_data_elemrng
% INPUTS ----------------------------------------------------------------
% data:  data structure which is output from 'readfiles.m' script. Should
%       contain fields which have <# lines> x <# samples> matrix for each
%       isotope analyzed
% varargin - 'PropertyName','PropertyValue'
%   'imgdim': image-dimension => scalar indicating whether to find boundary for each
%   isotope (1), or whether to combine isotopes specified and find boundary
%   using N-isotope data (2) [DEFAULT = 2]
% OUTPUTS ---------------------------------------------------------------
% datafxed: structure with same fields as data which contains data with
% single pixel anomalies  fixed.
%
%  Date           Author            E-mail                      Version
%  25 Sept 2011   Amanda Gaudreau   amanda.gaudreau@gmail.com     1

boundary = [];
fntsz = 14;
[elem_ana_range,fldnm] = interp_data_elemrng(data,elem_ana_range);

PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

if strmatch('imgdim',lower(PropertyNames))
  imD = PropertyVal{strmatch('imgdim',lower(PropertyNames))};
else
  imD = 2;
end

%i_la0 = find_laoff_ind(data,ele_ana_range,varargin{:});
%regions = setfield(regions,'i_la0',i_la0);
if imD == 1
  for f = elem_ana_range
    if isstruct(data)
      d = getfield(data,fldnm{f});
    else
      d = data;
    end
    if isempty(strmatch(fldnm{f},'Time')) && isempty(strmatch(fldnm{f},'line_datevec')) ...
        && isempty(strmatch(fldnm{f},'dataset')) && f <= length(fldnm)
      t = auto_thresh(d,[0.1,2],[]);
      d(d >= t(2)) = t(2);
      d(d <= t(1)) = t(1);
      d = d./(max(max(d)));
      elem_str = fldnm{f};
      if isfield(data,'dataset')
        data_str = getfield(data,'dataset');
      else
        data_str = 'test image';
      end
      
      seg = chanvese(d,'method','vector',varargin{:},...
        'plot_title',sprintf('%s %s ',data_str,elem_str));
      boundary = setfield(boundary,fldnm{f},seg);
    end
  end
  
elseif imD == 2
  if isstruct(data)
    dummy = getfield(data,fldnm{2});
    if ~isempty(find(isnan(dummy) == 1))
      nan_r = find(any(isnan(dummy),2)==1);
      included_rows = find(any(isnan(dummy),2)==0);
      d_non0 = dummy(included_rows,:);
      if isempty(d_non0)
        d_non0 = dummy(:,any(isnan(dummy),1));
      end
    else
      d_non0 = dummy;
    end
    dummy = d_non0;
  end
  
  d = zeros(size(dummy,1),size(dummy,2),length(elem_ana_range));
  elem_str = '';
  data_str = '';
  if isfield(data,'dataset')
    data_str = getfield(data,'dataset');
  end
  i =1;
  for f = elem_ana_range
    if isempty(strmatch(fldnm{f},'Time')) && isempty(strmatch(fldnm{f},'line_datevec')) ...
        && isempty(strmatch(fldnm{f},'dataset')) && f <= length(fldnm)
      if isstruct(data)
        dummy = getfield(data,fldnm{f});
        if ~isempty(find(isnan(dummy) == 1))
          nan_r = find(any(isnan(dummy),2)==1);
          included_rows = find(any(isnan(dummy),2)==0);
          d_non0 = dummy(included_rows,:);
          if isempty(d_non0)
            d_non0 = dummy(:,any(isnan(dummy),1));
          end
        else
          d_non0 = dummy;
        end
        dummy = d_non0;
        t = auto_thresh(dummy,[0.1,2],[]);
        dummy(dummy >= t(2)) = t(2);
        dummy(dummy <= t(1)) = t(1);
        dummy = ((dummy-min(dummy(:)))/(max(dummy(:))-min(dummy(:))));%.*255;
        d(:,:,i) = dummy; %maps data so that entire range falls between 0 and 1
        if ~isempty(elem_str)
          elem_str = sprintf('%s, %s',elem_str,fldnm{f});
        else
          elem_str = fldnm{f};
        end
      else
        d(:,:,i) = (data-min(data(:)))/(max(data(:))-min(data(:)));
        elem_str = fldnm{f};
      end
      i = i+1;
    end
  end
  
  [boundary,phi] = chanvese(d,'method','vector',varargin{:},...
    'plot_title',sprintf('%s %s ',data_str,elem_str));
  dockf on all
  varargout{1} = phi;
end
end

