function [data_struct,varargout] = readfilesNC(varargin);
%% Script:  data_struct = readfiles;
% Description:  Read Mass Spec Data in TXT or FIN2 files. Data files should be in a
% subfolder 'RawData' by themselves. When prompted, user should select the
% folder above RawData. Data from each isotope will be saved into a .csv
% file in the selected directory. Script was tested on 12-11-10 Tg 18
% Brain' data and took 7 seconds to run w/o writing csv files and 10.8
% seconds when writing the 10 csv files.
% Example:  d = readfiles('dat','20150102_icapqbrain');
%           d = readfiles('dat','20150102_icapqbrain','save_data',0);
%           [d,hdrtxt] = readfiles('20150218_1534_MIMSTest__ICAPQ_S1');
% d = readfiles('comp','obsolete','data_dir','parentdir', 'readcsv',1,...
%       'type','i','save_data',1,'line_offset',#);
% Required Function: loaddirfun, skip_fields, auto_thresh, save32bitTIF
% INPUTS ----------------------------------------------------------------
% varargin - 'PropertyName','PropertyValue'
%   • 'computer':  string indicating 'pc' or 'pho' computer [DEFAULT = ADGB'S PC]
%   • 'data_dir':  string indicating the data directory [DEFAULT = user selected]
%   • 'readcsv':   binary indicating whether to read csv files if they exist
%   for the data [DEFAULT = 1]
%   • 'type': string indicating 'i' for intensity data (cps) or 'c' for
%   analytical concentration data (in ppm)
%     NOTE: 'c' is only possible if there are 'csv' files with ELEMENT_ppm.csv
%   • 'save_data': binary indicating whether to save the raw data read in (1)
%   or not (0) [DEFAULT = 1]
%   • 'line_offset': integer indicating the offset of the line (ie if first
%   line is line 19, line offset should equal 18) [DEFAULT = 0]
%   • 'thresh': 1 x 2 matrix with [MIN%,MAX%] values (must be less than 10)
%   indicating the percentages of data contained within values between min
%   and lower thresh AND max and upper thresh. If one number, it will be
%   treated as a logical -- 1 to save thresholded data (with 'auto' input
%   into auto_thresh function) and 0 not to [DEFAULT = 0]
%   • 'save32': logical indicating whether to save a 32-bit tif file or not
%   in the 'Images' folder
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
%  Date           Author            E-mail                      Version
%  20 July 2011   Amanda Gaudreau   amanda.gaudreau@gmail.com     0
%  11 Oct  2011   Amanda Gaudreau   amanda.gaudreau@gmail.com     1
%     Redefining data structure so that it includes data set information
%     (i.e. parent directory file name of RawData folder)
%  28 June 2012   Amanda Gaudreau   amanda.gaudreau@gmail.com     2
%     Combined readfilesOES and readfiles (originally created exclusively
%     to read ICP-MS data). The data acquisition instrument is
%     automatically detected based on file type and header structure.
%   6 Jan  2015   Amanda Balderrama amanda.gaudreau@gmail.com     3
%     Accommodations for iCAP-Q *.csv files
%  10 Feb  2015   Amanda Balderrama amanda.gaudreau@gmail.com     3.1
%     Added a conditional statement at line 344-346 which differentiates
%     the letter "i" from numbers
%  11 Feb  2015   Amanda Balderrama amanda.gaudreau@gmail.com     4
%     Added ability to save thresholded *.csv files in addition to raw data
%     *.csv files and 32-bit tifs
%  11 Feb  2015   Amanda Balderrama amanda.gaudreau@gmail.com     5
%  13 Feb  2015   Amanda Balderrama amanda.gaudreau@gmail.com     5.1
%  17 Feb  2015   Amanda Balderrama amanda.gaudreau@gmail.com     6
%     Changed folders data are saved into per Noel's request (e-mail 17
%     Feb)
%  21 Feb  2015   Amanda Balderrama amanda.gaudreau@gmail.com     7
%     Enabled ability to parse header information

%% ---- Read data from text files, assumes comma delimiter ----
signal_point = 1;
compdir = loaddirfun;

if nargin == 1
  PropertyNames = {'dat'};
  PropertyVal = varargin;
else
  PropertyNames = varargin(1:2:length(varargin));
  PropertyVal = varargin(2:2:length(varargin));
end

% ---------------------------------------------------------------------
% • Laser Ablation data parent directory
% ---------------------------------------------------------------------
if strmatch('comp',PropertyNames)
  v = PropertyVal{strmatch('comp',PropertyNames)};
  if strmatch(v,'pc')
    current_dir = compdir.lafldr; % 'C:\Users\ADGB\Documents\MADLab Data\Laser Ablation\';
  elseif strmatch(v,'pho')
    current_dir = 'C:\Users\MADLAB\Desktop\Laser Ablation\';
  elseif ~isempty(v)
    current_dir = v;
  else
    current_dir = cd;
  end
else
  current_dir = compdir.lafldr; % 'C:\Users\ADGB\Documents\MADLab Data\Laser Ablation\';
end

if ~isdir(current_dir)
  current_dir = cd;
end

% ---------------------------------------------------------------------
% • Specific dataset folder
% ---------------------------------------------------------------------
if strmatch('dat',PropertyNames)
  v = PropertyVal{strmatch('dat',PropertyNames)};
  if strmatch('\',v)
    datahomedir = v;
  else
    datahomedir = strcat(current_dir,v);
  end
else
  datahomedir = uigetdir(current_dir,'Select the directory containing the data you wish to analyze'); % Select directory with .Fin2 files
end
while ~isdir(datahomedir)
  datahomedir = uigetdir(current_dir,'Select the directory containing the data you wish to analyze'); % Select directory with .Fin2 files
end

if strmatch('readcsv',PropertyNames)
  readcsv = PropertyVal{strmatch('readcsv',PropertyNames)};
else
  readcsv = 1;
end

if strmatch('thresh',PropertyNames)
  th10 = PropertyVal{strmatch('thresh',PropertyNames)};
else
  th10 = 1;
end

if strmatch('type',PropertyNames)
  datatype = PropertyVal{strmatch('type',PropertyNames)};
else
  datatype = 'i';
end

if strmatch('c',datatype)
  data_str = '_ppm';
else
  data_str = '';
end

if th10
  data_str = strcat('_th',data_str);
end

if strmatch('save_data',PropertyNames)
  svdat = PropertyVal{strmatch('save_data',PropertyNames)};
elseif strmatch('savedata',PropertyNames)
  svdat = PropertyVal{strmatch('savedata',PropertyNames)};
else
  svdat = 1;
end

if strmatch('save32',PropertyNames)
  svtif = PropertyVal{strmatch('save32',PropertyNames)};
else
  svtif = 1;
end

if strmatch('csv2folder',PropertyNames)
  csv2folder = PropertyVal{strmatch('csv2folder',PropertyNames)};
else
  csv2folder = 1;
end

if strmatch('writexls',PropertyNames)
  multisheetxl = PropertyVal{strmatch('writexls',PropertyNames)};
else
  multisheetxl = 1;
end

if strmatch('line',PropertyNames)
  line_offset = PropertyVal{strmatch('line',PropertyNames)};
else
  line_offset = 0;
end


if strmatch('elorder',PropertyNames)
  elordercell = PropertyVal{strmatch('elorder',PropertyNames)};
else
  elordercell = {'Fe','Zn','Cu','Mn','C_','P_','S_','Li','K_','Mo','Ce','Gd','Au'};
end

patterncases = {'_th_ppm','_ppm','_th',''};
nmodifiers = length(patterncases);
patternnum = find(ismember(patterncases,data_str));
date_vec = [];
exit_while = 0;

switch patternnum
  case 3
    disp('Output data structure represents thresholded raw data')
  case 4
    disp('Output data structure represents raw data')
end

% Checks for existing csv files so that information does not need to be
% re-read
cd (datahomedir); % go to that directory

allfiles = ls;
files_cell = cellstr(allfiles);
csvindx = strfind(files_cell,'.csv');
f_indx = find(cellfun(@isempty,csvindx) == 0);
allcsvfiles = files_cell(f_indx);
csvfiles = allcsvfiles;

if ~isempty(csvfiles) && readcsv
  nothfiles = 0;
  thppmstr = strfind(csvfiles,strcat(patterncases{1},'.csv'));
  thppmlogical = cellfun(@(x) ~isempty(x),thppmstr);
  thppmfiles = csvfiles(thppmlogical);
  csvfiles = csvfiles(~thppmlogical);
  
  ppmstr = strfind(csvfiles,strcat(patterncases{2},'.csv'));
  ppmlogical = cellfun(@(x) ~isempty(x),ppmstr);
  ppmfiles = csvfiles(ppmlogical);
  csvfiles = csvfiles(~ppmlogical);
  
  thstr = strfind(csvfiles,strcat(patterncases{3},'.csv'));
  thlogical = cellfun(@(x) ~isempty(x),thstr);
  thfiles = csvfiles(thlogical);
  csvfiles = csvfiles(~thlogical);
  
  switch patternnum
    case 1 %_th_ppm
      if isempty(thppmfiles)
        nothfiles = 1;
        if isempty(ppmfiles)
          disp('There are no concentration csv files available, Raw Data will be used')
          readcsvfiles = csvfiles;
          fileending = '.csv';
        else
          disp('Output data structure represents thresholded ppm concentration data')
          readcsvfiles = ppmfiles;
          fileending = '_ppm.csv';
        end
      else
        disp('Output data structure represents thresholded ppm concentration data')
        readcsvfiles = thppmfiles;
        fileending = strcat(patterncases{patternnum},'.csv');
      end
    case 2 %_ppm
      if isempty(ppmfile)
        disp('There are no concentration csv files available, Raw Data will be used')
        readcsvfiles = csvfiles;
        fileending = '.csv';
      else
        disp('Output data structure represents ppm concentration data')
        readcsvfiles = ppmfiles;
        fileending = strcat(patterncases{patternnum},'.csv');
      end
    case 3
      if isempty(thfiles)
        nothfiles = 1;
        readcsvfiles = csvfiles;
        fileending = '.csv';
      else
        readcsvfiles = thfiles;
        fileending = strcat(patterncases{patternnum},'.csv');
      end
    case 4
      readcsvfiles = csvfiles;
      fileending = '.csv';
  end
  
  el_names = cellfun(@(x) strrep(x,fileending,''),readcsvfiles,'uniformoutput',0);
  data_struct = [];
  for i = 1:length(el_names)
    d = csvread(readcsvfiles{i});
    if th10 && nothfiles
      if isempty(strmatch(el_names{i},skip_fields))
        [~,~,d] = auto_thresh(d,'auto',[]);
      end
      data_struct = setfield(data_struct,el_names{i},d);
      csvwrite(strcat(el_names{i},'_th',fileending),d);
    else
      data_struct = setfield(data_struct,el_names{i},d);
    end
  end
  d = csvread('Time.csv');
  data_struct = setfield(data_struct,'Time',d);
  i = strfind(datahomedir,'\');
  data_struct = setfield(data_struct,'dataset',datahomedir(i(end)+1:end));
  d = csvread('line_datevec.csv');
  data_struct = setfield(data_struct,'line_datevec',d);
else % this is executed if the raw data files need to be read
  % Need to determine the instrument that data is being read from. If it's
  % from the ICP-MS, the header of the files as well as the naming and
  % number of files will differ. Both checks will be performed.
  cd('RawData');
  rawdatadir = cd;
  % 1) determine the format of the files in the RawData folder
  rawdatafiles = ls;
  rawdatafiles_cell = cellstr(rawdatafiles);
  if strmatch(lower(rawdatafiles_cell{3}(end-2:end)),'txt')
    fmt = '.TXT';
    format = 't';
    hlines = 6;
  elseif strmatch(lower(rawdatafiles_cell{3}(end-2:end)),'csv')
    fmt = '.CSV';
    format = 'c';
    hlines = 15;
  elseif strmatch(lower(rawdatafiles_cell{3}(end-3:end)),'fin2')
    fmt = '.FIN2';
    format = 'f';
    hlines = 8;
  end
  
  if strmatch('.TXT',fmt)
    ftest = strcat(rawdatadir,'\',rawdatafiles_cell{3});
    fidtest = fopen(ftest);
    for i = 1:20
      line = fgetl(fidtest);
    end
    if length(strtok(line,',')) == 6 % this assumes that if the text file is from the OES instrument (=0) then first element of a line with data will be a 6 letter string specifying the element that the data is for.
      ms1_oes0 = 0; % FROM OES
    else % this assumes that the first element of a text file from the ICP-MS instrument will be a number that is longer than 6 digits (usually ~16 digits)
      ms1_oes0 = 1; % FROM MASS SPEC
    end
    fclose(fidtest);
  elseif strmatch('.FIN2',fmt)
    ms1_oes0 = 1;
  elseif strmatch('.CSV',fmt)
    ms1_oes0 = 1;
  end
  
  if svdat
    cd(datahomedir)
    if ~isdir('CSV')
      mkdir('CSV')
    end
    cd(rawdatadir)
  end
  
  if ms1_oes0 == 1 % data came from ICP-MS
    %% Script read out ICP-MS data files
    keepfiles = cellfun(@(x) ~isempty(x),strfind(lower(rawdatafiles_cell),lower(fmt)));
    rawdatafiles_cell = rawdatafiles_cell(keepfiles);
    file_num_cell = strrep(lower(rawdatafiles_cell),lower(fmt),'');
    uniq_el = cell(0,1);
    %The following code within the "if isnan" loop attempts to recognize the
    %naming pattern of the data files if there is one. The pattern must be
    %PATTERN###.FMT where ### indicates the line number, FMT indicates file
    %format and PATTERN is the header pattern on each file.
    file_num = str2double(file_num_cell) - line_offset; % vector containing numeric values of .txt files
    if isnan(file_num(1))
      end_pattern_ind = 0;
      len_fnames = length(file_num_cell{1});
      pat = zeros(len_fnames,1);
      if length(file_num_cell) > 1
        c = len_fnames;
        while ~end_pattern_ind && c > 0
          pat(c) = length(strmatch(file_num_cell{1}(1:c),file_num_cell(2:end)));
          if pat(c) == length(file_num_cell)-1
            end_pattern_ind = c;
          end
          c = c - 1;
        end
        if ~end_pattern_ind
          end_pattern_ind = max(find(max(pat) == pat));
        end
      else
        end_pattern_ind = 0;
      end
      rawdatafiles_cell = cellfun(@(x) x(end_pattern_ind+1:end),file_num_cell,'uniformoutput',0);
      file_num = str2double(rawdatafiles_cell)- line_offset; % vector containing numeric values of .txt files
    end
    
    max_line = max(file_num);
    if length(file_num) >= max_line
      line_rng = 1:max_line;
    elseif isnan(file_num)
      line_rng = 1;
      max_line = 1;
      file_num = 1;
    else
      line_rng = 1:length(file_num);
    end
    
    A = importdata(strcat(file_num_cell{1},fmt), ',', hlines); %
    [r,c] = size(A.data); % r: number of time samples per line; c: # elements + 1 (one column reserved for time)
    data = nan(r,max_line,c);
    line_datevec = cell(max_line,1);
    hdrcol = cell(1,83);
    hdrdata = cell(max_line,83);
    hdrcol{1} = 'Line Name';
    hdrcol{2} = 'Time Stamp';
    
    for k = line_rng
      A = importdata(strcat(file_num_cell{k},fmt), ',', hlines);
      if strmatch(lower(format),'f')
        date_str = A.textdata{2,1};
        if ~isnan(file_num(k))
          line_datevec{file_num(k),1} = datevec(date_str,'dddd, mmmm dd, yyyy HH:MM:SS');
        end
      elseif strmatch(lower(format),'c')
        linetxt = A.textdata{1,1};
        colpos = strfind(linetxt,':');
        if strfind(linetxt,'PM')
          hhstr = str2num(linetxt(colpos(2)-2:colpos(2)-1));
          hhstr = hhstr + 12 - floor(hhstr/12)*12;
          linetxt(colpos(2)-2:colpos(2)-1) = num2str(hhstr);
        end
        date_str = linetxt(colpos+1:end-3);
        if ~isnan(file_num(k))
          line_datevec{file_num(k),1} = datevec(date_str,'mm/dd/yyyy HH:MM:SS');
        end
        
        hdrline = 1;
        [htxt,linetxt] = strtok(A.textdata{hdrline},':');
        linetxt = linetxt(2:end-1);
        hdrdata{file_num(k),1} = htxt;
        hdrdata{file_num(k),2} = linetxt;
        fillcol = 3;
        
        hdrline = 2;
        while ~isempty(A.textdata{hdrline})
          [htxt,linetxt] = strtok(A.textdata{hdrline},':');
          linetxt = linetxt(2:end);
          linefields = strsplit(linetxt,';');
          parsetxt = cellfun(@(x) ~isempty(x),linefields);
          for p = 1:sum(parsetxt)
            fldandval = strsplit(linefields{p},'=');
            if k == 1
              hdrcol{fillcol} = fldandval{1};
            end
            hdrdata{file_num(k),fillcol} = fldandval{2};
            fillcol = fillcol + 1;
          end
          hdrline = hdrline + 1;
        end
      else
        line_datevec = nan;
      end
      if size(A.data,1) ~= r % error display if the size of data isn't the same as the first line's size
        %disp(sprintf('The file %s%s has %d samples, differing from the standard %d samples expected.',file_num_cell{k},fmt, size(A.data,1),r));
      end
      %disp(file_num_cell{k});
      data(1:size(A.data,1),file_num(k),:) = A.data; %data(:,k,:) = A.data; %<time samples> x <# lines> x <# elements (+1 for time)> data matrix
    end
    cd ..
    hdrcell = [hdrcol;hdrdata];
    if svdat %isempty(cell2mat(existing_files)) || ~cell2mat(existing_files)
     xlswrite('headerinfo.xls',hdrcell);
    end
    varargout{1} = hdrcell;
    
    % Parse header text
    if strmatch(lower(format),'t')
      names =  A.textdata(1,1);
      n = strvcat(names(1,:));
      elem = textscan(n, '%s' ,'delimiter', ',') ;
      na = cell2struct(elem,'isos',1);
      na.isos{1} = 'Time';% [sec]';
      vn = na.isos;%strtok(na.isos, '(LR)');
    elseif strmatch(lower(format),'f')
      names =  A.colheaders;
      n = strvcat(names(1,:));
      elem{1} = cellstr(n);
      na = cell2struct(elem,'isos',1);
      na.isos{1} = 'Time';% [sec]';
      vn = na.isos;%(2:end);%strtok(na.isos, '(LR)');
    elseif strmatch(lower(format),'c')
      names =  A.textdata{14,1};
      elem = textscan(names, '%s' ,'delimiter', ',');
      na = cell2struct(elem,'isos',1);
      na.isos{1} = 'Time';
      vn = na.isos;
    end
    
    % Reshape data
    sz = size(data);
    col = sz(1);
    d = reshape(data(:),[sz(1)*sz(2),sz(3)]);
    %rr = floor(sz(1)*sz(2)/col);
    %dd = reshape(d(1:rr*col,:),col,rr,c); %<time samples> x <# lines> x <# elements (+1 for time)> data matrix
    
    % Save isotope data: Data is saved so that the rows contain line data
    % and the columns represent the intensity value at a given sample time
    for k = 1:sz(3)
      temp = reshape(d(:,k),[sz(1),sz(2)])';%transpose(dd(:,:,k)); % <# lines> x <# time samples>
      % temp(isnan(temp)) = min(temp(:)); % UNCOMMENT if you'd like only
      % numeric (non-NaN) data in the matrix
      elem_name = vn{k};
      if ~isempty(strfind(elem_name,'>'))
        elem_name = elem_name(1:strfind(elem_name,'>')-1);
      elseif ~isempty(strfind(elem_name,'('))
        elem_name = elem_name(1:strfind(elem_name,'(')-1);
      elseif ~isempty(str2num(elem_name(1))) % indicates ##EL format (ex: 12C)
        elparsed = strsplit(elem_name,'.');
        elem_name = '';
        for e = 1:length(elparsed)
          nel = length(elparsed{e});
          letters = [];
          while isempty(letters)
            letters = str2num(elparsed{e}(nel));
            nel = nel - 1;
            if ~isreal(letters)
              letters = [];
            end
          end
          nel = nel+1;
          if (nel+1) == length(elparsed{e})
            elparsed{e} = strcat(elparsed{e}(nel+1:end),'_',elparsed{e}(1:nel));
          else
            elparsed{e} = strcat(elparsed{e}(nel+1:end),elparsed{e}(1:nel));
          end
          elem_name = strcat(elem_name,elparsed{e});
        end
      end
      if k == 1
        data_struct = cell2struct({temp},'Time',1);
      else
        data_struct = setfield(data_struct,elem_name,temp);%cell2struct({temp},elem_name,1);
        uniq_el = [uniq_el,elem_name];
      end
      csvname = strcat(elem_name,'.csv');
      if svdat %isempty(cell2mat(existing_files)) || ~cell2mat(existing_files)
        csvwrite(csvname,temp);
        if csv2folder && ~ismember(elem_name,skip_fields)
          copyfile(csvname,strcat('CSV\',csvname));
        end
      end
    end
    
    csvname = 'Time.csv';%strcat(datahomedir,'\CSV\','Time.csv');
    if svdat %isempty(cell2mat(existing_files)) || ~cell2mat(existing_files)
      csvwrite(csvname,data_struct.Time);
    end
    
    data_struct.line_datevec = line_datevec;
    csvname = 'line_datevec.csv';%strcat(datahomedir,'\CSV\','line_datevec.csv');
    if svdat %isempty(cell2mat(existing_files)) || ~cell2mat(existing_files)
      csvwrite(csvname,data_struct.line_datevec);
    end
    
    dir_sep = strfind(datahomedir(1:end-1),'\');
    dir_LA = [];%strfind(datahomedir,'Laser Ablation');
    if ~isempty(dir_LA)
      start_ind = dir_sep(find(dir_sep == dir_LA-1)+1)+1;
      end_ind = dir_sep(find(dir_sep == dir_LA - 1) + 2) - 1;
      data_struct.dataset = rawdatadir(start_ind:end_ind);
    else
      start_ind = dir_sep(end) + 1;
      end_ind = length(datahomedir);
      data_struct.dataset = datahomedir(start_ind:end_ind);
    end
    
    if svdat 
      elnames = fieldnames(rmfield(data_struct,skip_fields));
      xlname = strcat('CSV\',data_struct.dataset,' Element Data.xls')
      for el = 1:length(elnames)
        csvname = strcat(elnames{el},'.csv');
        csvwrite(csvname,getfield(data_struct,elnames{el}));
        
        if csv2folder
          copyfile(csvname,strcat('CSV\',csvname));
        end
        
        if multisheetxl
          xlswrite(xlname,getfield(data_struct,elnames{el}),elnames{el});
        end
        
      end
    end
    %========================================================================
  elseif ms1_oes0 == 0
    %% Script read out OES data files
    data_struct = [];
    line_number = 0;
    l = cell(0);
    
    if strmatch('file',PropertyNames)
      fname = PropertyVal{strmatch('file',PropertyNames)};
    elseif ~exist('fname')
      [fname,rawdatadir,fid] = uigetfile('*.txt','Select the OES TXT data file');
    end
    fname = strcat(rawdatadir,'\',fname);
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
          if  isempty(date_vec) || ~ismember(datevec(date_str,'mm/dd/yyyy HH:MM:SS'),date_vec,'rows')
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
    usecol = cellfun(@(x) find(sum(x) ~= 0),d,'uniformoutput',false);
    userow = cellfun(@(x) find(sum(x,2) ~= 0),d,'uniformoutput',false);
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
    csvfiles = ls;
    
    
    for i = 1:length(uniq_el)
      data_struct = setfield(data_struct,uniq_el{i},d{signal_point,i}(userow{i},usecol{i}));
      csvname = strcat(home_dir,uniq_el{i},'.csv');
      existing_files = strfind(cellstr(csvfiles),strcat(uniq_el{i},'.csv'));
      if svdat %isempty(cell2mat(existing_files)) || ~cell2mat(existing_files)
        csvwrite(csvname,d{signal_point,i}(userow{i},usecol{i}));
        if csv2folder
          copyfile(csvname,strcat('\CSV\',csvname));
        end
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
  
  if svtif; save32bitTIF(data_struct,'','saveto','32Tiff','dataheader',0); end
  
  if th10
    if svdat
      cd(datahomedir)
      if ~isdir(strcat('CSV',data_str))
        mkdir(strcat('CSV',data_str))
      end
    end
    
    dth = data_struct;
    skpfldnm = skip_fields;
    dth = rmfield(dth,skpfldnm);
    evalelnms = fieldnames(dth);
    dthcell = struct2cell(dth);
    for d = 1:length(dthcell)
      [~,~,dthcell{d}] = auto_thresh(dthcell{d},'auto',[]);
    end
    dth = cell2struct(dthcell,evalelnms);
    for s = 1:length(skpfldnm)
      dth = setfield(dth,skpfldnm{s},getfield(data_struct,skpfldnm{s}));
    end
    
    data_struct = dth;
    if svdat %isempty(cell2mat(existing_files)) || ~cell2mat(existing_files)
      fldnm = fieldnames(data_struct);
      evalflds = find(cellfun(@(x) isempty(x),strfind(fldnm,'dataset')));
      for f = evalflds(:)'
        csvname = strcat(fldnm{f},data_str,'.csv');
        csvwrite(csvname,getfield(data_struct,fldnm{f}));
        if csv2folder && ~ismember(fldnm{f},skip_fields)
          copyfile(csvname,strcat('CSV',data_str,'\',csvname));
        end
      end
      
      if multisheetxl
        elnames = fieldnames(rmfield(data_struct,skip_fields));
        xlname = strcat('CSV',data_str,'\',data_struct.dataset,' Element Data.xls');
        for el = 1:length(elnames)
          xlswrite(xlname,getfield(data_struct,elnames{el}),elnames{el});
        end
      end
      
    end
  end
end

if svtif; save32bitTIF(data_struct,data_str,'saveto',strcat('32Tiff',data_str),'dataheader',0); end

if exist('headerinfo.xls') == 2
  [~,~,varargout{1}] = xlsread('headerinfo.xls');
else
  disp('No header info file is available, please rerun the function with readcsv parameter set to 0');
  varargout{1} = [];
end

