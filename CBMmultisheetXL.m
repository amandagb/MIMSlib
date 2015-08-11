function CBMmultisheetXL(xlpath,d,elendnum,hdrcell);
%% Script:  CBMmultisheetXL(xlname,xlpath,d,elendnum,hdrcell);
% Description:  
% Example:  d = readfiles('dat','20150102_icapqbrain');
%           d = readfiles('dat','20150102_icapqbrain','save_data',0);
%           [d,hdrtxt] = readfiles('20150218_1534_MIMSTest__ICAPQ_S1');
% d = readfiles('comp','obsolete','data_dir','parentdir', 'readcsv',1,...
%       'type','i','save_data',1,'line_offset',#);
% Required Function: loaddirfun, skip_fields, auto_thresh, save32bitTIF
% INPUTS ----------------------------------------------------------------
% OUTPUTS ---------------------------------------------------------------
%  Date           Author            E-mail                      Version
%  29 Mar  2015   Amanda Balderrama amanda.gaudreau@gmail.com     0

% xlname = strcat(datasetfdr,'_DataMatrix.xlsx');
% xlpath = strcat(datahomedir,'\',xlname);

elnumfirst = cellfun(@(x) strrep(strcat(x(3:end),x(1:2)),'_',''),elendnum,'uniformoutput',0);

nel = length(elendnum);
M = size(d.Time,1)+1;
N = size(d.Time,2) + 2;
basecell = cell(M, N);
basecell(1:end,1) = hdrcell(1:M,2);
basecell{1,1} = 'Date';
basecell(1:end,2) = hdrcell(1:M,1);
basecell(1,3:end) = num2cell(d.Time(1,:));
elorder = hdrcell(2:nel+1,strmatch('Element Names',hdrcell(1,:)):end);
% Excel = actxserver ('Excel.Application');
% File=xlpath;
% if ~exist(File,'file')
%     ExcelWorkbook = Excel.workbooks.Add;
%     ExcelWorkbook.SaveAs(File,1);
%     ExcelWorkbook.Close(false);
% end
% invoke(Excel.Workbooks,'Open',File);

for i = 1:nel
  elnum1str = elorder{i};
  eldata = getfield(d,elendnum{ismember(elnumfirst,elnum1str)});
  dcell = basecell;
  dcell(2:end,3:end) = num2cell(eldata);
  xlswrite(xlpath,dcell,elnum1str);
end
eldata = getfield(d,'Time');
dcell = basecell;
dcell(2:end,3:end) = num2cell(eldata);
xlswrite(xlpath,dcell,'TimeMatrix');

dcell = hdrcell(1:nel+1,strmatch('Element Names',hdrcell(1,:)):end);
xlswrite(xlpath,dcell,'Element List');

dcell = hdrcell(:,1:(strmatch('Element Names',hdrcell(1,:))-1));
xlswrite(xlpath,dcell,'HeaderFile Info');

sheetName = 'Sheet'; % EN: Sheet, DE: Tabelle, etc. (Lang. dependent)
Excel = actxserver('Excel.Application');
Excel.Workbooks.Open(fullfile(xlpath)); % Full path is necessary!
try
  Excel.ActiveWorkbook.Worksheets.Item([sheetName '1']).Delete;
catch
      ; % Do nothing.
end
% Save, close and clean up.
Excel.ActiveWorkbook.Save;
Excel.ActiveWorkbook.Close;
Excel.Quit;
Excel.delete;
