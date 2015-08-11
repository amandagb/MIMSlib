function [datastruct,hdrtxt,d_div,line_types,type_rows,varargout] = divMIMSstruct(datastruct,hdrtxt,varargin)
%% function: divMIMSstruct
% Description:
% Example:
%     
%
% INPUTS ----------------------------------------------------------------
% datastruct:    MIMS structure OR a string indiciating the dataset the
%     user would like to use. 
% hdrtxt: nrows + 1 x 86 cell which contains the header data from the
%     instrument CSV files. Columns 84 - 86 have element information while
%     all other columns have per line data
% 
%
% OUTPUTS ---------------------------------------------------------------
%
%  Date           Author              E-mail                      Version
%  24 Mar  2015   Amanda Balderrama   amandagbalderrama@gmail.com     1

compdir = loaddirfun;

PropertyNames = lower(varargin(1:2:length(varargin)));
PropertyVal = varargin(2:2:length(varargin));

if strmatch('plotimgs',PropertyNames)
  plotimgs = PropertyVal{strmatch('plotimgs',PropertyNames)};
else
  plotimgs = 0;
end

%% Reads raw data
if isstr(datastruct) || ~exist('hdrtxt')
  [datastruct,hdrtxt] = readfiles(datastruct);
end

if isempty(hdrtxt)
  [~,hdrtxt] = readfiles(datastruct.dataset);
end
linelabels = hdrtxt(2:end,1);
uselines = find(cellfun(@(x) ischar(x),linelabels));
linelabels = linelabels(uselines);
delonly = rmfield(datastruct,skip_fields);
elnames = sort(fieldnames( delonly ));

%% Automatically divides data by line names indicated in header files
lineCO = strfind(lower(linelabels),'_line');
relhtxt = cellfun(@(x,y) x(1:(y-1)),linelabels,lineCO,'uniformoutput',0);
reptxtind = find(cellfun(@(x) isempty(x),relhtxt));
relhtxt(reptxtind) = linelabels(reptxtind);
[line_types,~,linesind] = unique(relhtxt,'stable');
dnodataset = rmfield(datastruct,'dataset');
dnodataset = rmfield(dnodataset,'line_datevec');
dnodataset.Time(isnan(datastruct.Time)) = 0;
d_div = cell(length(line_types),1);
type_rows = cell(length(line_types),1);
for i = 1:length(d_div)
  type_rows{i} = uselines(linesind == i);
  dtype = structfun(@(x) x(type_rows{i},:),dnodataset,'uniformoutput',0);
  dtype.line_datevec = datevec(hdrtxt(type_rows{i}+1,2));
  dtype.dataset = strcat(datastruct.dataset,'_',line_types{i});
  dtype_elonly = rmfield(dtype,skip_fields);
  keepcols = find(sum(dtype.Time,1) > 0);
  dtype.Time = dtype.Time(:,keepcols);
  dtype_elonly = structfun(@(x) x(:,keepcols),dtype_elonly,'uniformoutput',0);
  for el = 1:length(elnames)
   dtype = setfield(dtype,elnames{el},getfield(dtype_elonly,elnames{el}));
  end
  d_div{i} = dtype;
  %plot_heatmap(dtype,{'cu63','zn66','fe57'},'thresh','auto');
end

end
