function regions = id_MIMS_regions(data,elem_ana_range,varargin);

%% Script:  regions = id_MIMS_regions(data,elem_ana_range,varargin);
% Description:
% Required Functions: interp_data_elemrng
% INPUTS ----------------------------------------------------------------
% data:  data structure which is output from 'readfiles.m' script. Should
%       contain fields which have <# lines> x <# samples> matrix for each
%       isotope analyzed
% varargin - 'PropertyName','PropertyValue'
%   'lines': row matrix indicating specific lines user would like to
%   analyze [DEFAULT = 'all']
%   'plot': binary indicating whether to show plots (1) or not (0) [DEFAULT
%   = 0]
% OUTPUTS ---------------------------------------------------------------
% datafxed: structure with same fields as data which contains data with
% single pixel anomalies  fixed.
%
%  Date           Author            E-mail                      Version
%  13 Sept 2011   Amanda Gaudreau   amanda.gaudreau@gmail.com     1
%  15 Sept 2011   Amanda Gaudreau   amanda.gaudreau@gmail.com     3

regions = [];
fntsz = 14;
[elem_ana_range,fldnm] = interp_data_elemrng(data,elem_ana_range);

PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

i_la0 = find_laoff_ind(data,ele_ana_range,varargin{:});
regions = setfield(regions,'i_la0',i_la0);

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
i =1;
for f = elem_ana_range
  if isempty(strmatch(fldnm{f},'Time')) && isempty(strmatch(fldnm{f},'line_datevec')) && f <= length(fldnm)
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
      d(:,:,i) = (dummy-min(dummy(:)))/(max(dummy(:))-min(dummy(:))); %maps data so that entire range falls between 0 and 1
    else
      d(:,:,i) = (data-min(data(:)))/(max(data(:))-min(data(:)));
    end
    i = i+1;
  end
end

%% Finds tissue edges on a line-per-line basis


end

