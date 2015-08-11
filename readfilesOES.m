function data_struct = readfilesOES(varargin);
%% Script:  data_struct = readfilesOES(varargin);
%
% Description:  Reads optical emissions spectrometer data from a text file.
% Data file(s) should be in a subfolder 'RawData' by themselves. When
% prompted, user should select the folder above RawData. Data from each
% isotope will be saved into a .csv file in the selected directory. If the
% csv files already exist for a selected directory, the data structure will
% read in data from these files (since it is about 100x faster). Script was
% tested on '012712_OEScalibration' data (34 MB txt file with 113 x 535 dim
% data) and took 21.13 seconds to run (which includes generating csv files
% for each element with data) and 0.30 seconds to read data from csv files.
% Example:  data_struct = readfilesOES('comp','pc','dat','012712_OEScalibration','file','1-27-2012 114 of 145 lines');
% INPUTS ----------------------------------------------------------------
% varargin - 'PropertyName','PropertyValue'
%   'computer':  string indicating 'pc' or 'pho' computer [DEFAULT = cd]
%   'data_dir':  string indicating the data directory [DEFAULT = user selected]
%   'file_name': string indicate the name of the text files the data is in
%   [DEFAULT = user selected]
%   'readcsv':   binary indicating whether to read csv files if they exist
%   for the data [DEFAULT = 1]
%   'type': string indicating 'i' for intensity data or 'c' for cps data.
%   'c' is only possible if there are 'csv' files with ELEMENT_ppm.csv
%   'save_data': binary indicating whether to save the raw data read in (1)
%   or not (0) [DEFAULT = 1]
% OUTPUTS ---------------------------------------------------------------
% data_struct:  a structure contining a field with the name of each
%   element for which data was collected. Each field contains a <# lines> x
%   <# time samples> double matrix.
% csv files: # elements + 2 *.csv files are written to the data directory
%   (the parent directory of the 'RawData' folder) one csv file is written
%   for each element and is organized so that each row is a line and each
%   column is a sample (acquired ponit). NOTE: the names of the *.csv files
%   for the elements all include a suffix letter. This letter corresponds
%   to the column of data from the txt file (L: left bkgnd, C: center peak,
%   R: right bkgn). For example, data from the central peak of Carbon
%   (lambda = 247.8nm) has a csv file name of 'CC_2478.csv'. There are NO
%   HEADERS. The Time and line_datevec *.csv files give the relative sample
%   time (indicated by the first isotope in the data file) and the
%   line_datevec contains a <#lines> x 6 matrix of data where the columns
%   correspond to YYYY MM DD HH MM SS. These can be interpretered as
%   strings using the datestr matlab function.
%
%  Date           Author            E-mail                      Version
%  30 Jan  2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     0
%  06 Feb  2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     1
%     This version is the same as the readfilesOES_analysis.m file
%  28 Feb  2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     2
%  29 Feb  2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     3
%     Saves data without preceeding signal string tags and does not read in
%     L and R background data
%  21 May  2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     4

%% ---- Read data from text files, assumes comma delimiter ----
signal_point = 1;
%signal_str = {'L','C','R'};

PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

if strmatch('comp',PropertyNames)
  v = PropertyVal{strmatch('comp',PropertyNames)};
  if strmatch(v,'pc')
    current_dir = 'D:\My Documents\MADLab Data\Laser Ablation';%'C:\Users\Amanda\Documents\MADLab Data\Laser Ablation\';%'C:\Users\MADLAB\Desktop\Laser Ablation';%cd; %
  elseif strmatch(v,'pho')
    current_dir = 'C:\Users\MADLAB\Desktop\Laser Ablation\';
  elseif ~isempty(v)
    current_dir = v;
  else
    current_dir = cd;
  end
else
  current_dir = 'D:\My Documents\MADLab Data\Laser Ablation';%'C:\Users\Amanda\Documents\MADLab Data\Laser Ablation\';
end

if strmatch('dat',PropertyNames)
  v = PropertyVal{strmatch('dat',PropertyNames)};
  dir_name = strcat(current_dir,v);
else
  dir_name = uigetdir(current_dir,'Select the directory containing the data you wish to analyze'); % Select directory with .Fin2 files
end
while ~isdir(dir_name)
  dir_name = uigetdir(current_dir,'Select the directory containing the data you wish to analyze'); % Select directory with .Fin2 files
end

if strmatch('readcsv',PropertyNames)
  readcsv = PropertyVal{strmatch('readcsv',PropertyNames)};
else
  readcsv = 1;
end

if strmatch('type',PropertyNames)
  datatype = PropertyVal{strmatch('type',PropertyNames)};
else
  datatype = 'i';
end

if strmatch('c',datatype)
  data_str = 'ppm';
else
  data_str = '';
end

if strmatch('save_data',PropertyNames)
  svdat = PropertyVal{strmatch('save_data',PropertyNames)};
else
  svdat = 1;
end

cd (dir_name); % got to that directory
allfiles = ls;
files = cellstr(allfiles);
csvfiles = strfind(files,strcat(data_str,'.csv'));
f_indx = find(cellfun(@isempty,csvfiles) == 0);
files = files(f_indx);
if strmatch('i',datatype)
  ppmfiles = strfind(files,'ppm');
  f_indx = find(cellfun(@isempty,ppmfiles) == 1);
  files = files(f_indx);
else
  ppmfiles = strfind(files,'ppm');
  readcsv = sum(cellfun(@isempty,ppmfiles) == 0) == length(ppmfiles);
  if readcsv == 0; disp('There is no concentration csv file available, Raw Data will be used for output'); end
end
OESfiles = strfind(files,'Time');
f_indx = find(cellfun(@isempty,OESfiles) == 0);
date_vec = [];
exit_while = 0;
if ~isempty(f_indx) && readcsv
  el_names = cellfun(@(x) x(1:end-4),files,'uniformoutput',false);
  ppm_part = strfind(el_names,'_ppm');
  if sum(cellfun(@isempty,ppm_part)) ~= length(ppm_part)
    el_names = cellfun(@(x,y) x(1:y-1),el_names,ppm_part,'uniformoutput',false);
  end
  data_struct = [];
  for i = 1:length(el_names)
    d = csvread(files{i});
    data_struct = setfield(data_struct,el_names{i},d);
  end
  d = csvread('Time.csv');
  data_struct = setfield(data_struct,'Time',d);
  i = strfind(dir_name,'\');
  data_struct = setfield(data_struct,'dataset',dir_name(i(end)+1:end));
  d = csvread('line_datevec.csv');
  data_struct = setfield(data_struct,'line_datevec',d);
else
  data_struct = [];
  line_number = 0;
  l = cell(0);
  cd('RawData');
  dir_name = cd;
  if strmatch('file',PropertyNames)
    fname = PropertyVal{strmatch('file',PropertyNames)};
  elseif ~exist('fname')
    [fname,dir_name,fid] = uigetfile('*.txt','Select the OES TXT data file');
  end
  fname = strcat(dir_name,'\',fname);
  fid = fopen(fname);
  validrows_perelement = cell(13,15);
  while iscell(l) || ischar(l) && ~exit_while
    if ~ischar(l)
      l = fgetl(fid);
      if l == -1
        exit_while = 1;
        break
      end
    end
    if ~exist('next_line_info'); next_line_info = ''; end
    while isempty(strmatch('Element',l))
      if ~isempty(strfind(l,'/'))
        l = strcat(next_line_info,l);
        line_number = line_number + 1;%   str2num( l(1) );
        slashes = strfind(l,'/');
        colons = strfind(l,':');
        date_str = l(slashes(1)-2:colons(end)+2);
        commas = strfind(date_str,',');
        if ~isempty(commas)
          date_str = date_str(commas(1)+1:end);
          date_str = strcat('0',date_str);
        end
        if  ~ismember(datevec(date_str,'mm/dd/yyyy HH:MM:SS'),date_vec,'rows')
          date_vec(line_number,:) = datevec(date_str,'mm/dd/yyyy HH:MM:SS');
        else
          exit_while = 1;
          line_number = line_number - 1;
        end
      end
      l = fgetl(fid);
    end
    if ~exit_while
      validrows_perelement{line_number+1,1} = line_number;
      next_line_info = '';
      headers = l;
      l = textscan(fid,'%s %f %d  %d %d',1e6,'delimiter',',');
      cell_len = cellfun(@length,l);
      min_cell_len = min(cell_len);
      uniq_el = unique(l{1}(1:min_cell_len));
      
      for k = 1:length(l)
        if length(l{k}) == max(cell_len)
          if k == 1
            next_line_info = strcat(next_line_info,l{k}{max(cell_len)},',');
          else
            next_line_info = strcat(next_line_info,num2str(l{k}(max(cell_len))));
          end
        end
      end
      if line_number == 1
        len_rows = zeros(length(uniq_el),1);
        for i = 1:length(uniq_el)
          rows = strfind(l{1},uniq_el{i});
          rows = cellfun(@isempty,rows);
          len_rows(i) = sum(rows == 0);
        end
        [v,usetime] = max(len_rows);
        Time = zeros(300,v);
        d = cell(1,length(uniq_el));
        [d{:}] = deal(zeros(300,v));
      end
      for i = 1:length(uniq_el)
        rows = strfind(l{1},uniq_el{i});
        rows = cellfun(@isempty,rows);
        valid_rows = find(rows == 0);
        Ndata_pnts = length(valid_rows) - 1; % The first row is ignored since for all elements, the instrument begins the count at 0. This sets the baseline count and is irrelevant to record as a data point
        validrows_perelement{line_number+1,i+1} = Ndata_pnts;
        el_line_time = l{2}(valid_rows);
        %el_data_lineL = diff(l{3}(valid_rows));
        %el_data_lineR = diff(l{5}(valid_rows));
        el_data_lineC = diff(l{4}(valid_rows));
        if line_number == 1
          validrows_perelement{1,i+1} = uniq_el(i);
        end
        if i == usetime; Time(line_number,1:length(valid_rows)) = el_line_time'; end
        d{i}(line_number,1:Ndata_pnts) = double(el_data_lineC');
        %d{2,i}(line_number,1:Ndata_pnts) = el_data_lineC';
        %d{3,i}(line_number,1:Ndata_pnts) = el_data_lineR';
      end
    end
  end
  
  d = cellfun(@(x) x(1:line_number,:),d,'uniformoutput',false);
  data_struct = [];
  Time = Time(1:line_number,2:end);
  dir_sep = strfind(fname,'\');
  dir_LA = strfind(fname,'Laser Ablation');
  if ~isempty(dir_LA)
    start_ind = dir_sep(find(dir_sep == dir_LA-1)+1)+1;
    end_ind = dir_sep(find(dir_sep == dir_LA - 1) + 2) - 1;
    data_struct.dataset = fname(start_ind:end_ind);
  end
  home_dir = fname(1:end_ind+1);
  cd(home_dir);
  files = ls;
  
  
  for i = 1:length(uniq_el)
    data_struct = setfield(data_struct,uniq_el{i},d{signal_point,i});
    csvname = strcat(home_dir,uniq_el{i},'.csv');
    existing_files = strfind(cellstr(files),strcat(uniq_el{i},'.csv'));
    if svdat %isempty(cell2mat(existing_files)) || ~cell2mat(existing_files)
      csvwrite(csvname,d{signal_point,i});
    end
  end
  data_struct = setfield(data_struct,'Time',Time);
  csvname = 'Time.csv';
  if svdat %isempty(cell2mat(existing_files)) || ~cell2mat(existing_files)
    csvwrite(csvname,Time);
  end
  
  data_struct = setfield(data_struct,'line_datevec',date_vec);
  csvname = 'line_datevec.csv';
  if svdat %isempty(cell2mat(existing_files)) || ~cell2mat(existing_files)
    csvwrite(csvname,date_vec);
  end
  
  fclose(fid);
end


