function [Iph,Lcell,Lcorrcell,bucsblk,IVISinfo] = umxIVISdata(varargin)
%% function: [Iph,Ilum1,Ilum2,bucsblk] = umxIVISdata(querypairs)
%
% Description: 
% Example:
% Required Functions: loaddirfun, IVISinfo.xlsx
%
% INPUTS ----------------------------------------------------------------
% 
% OUTPUTS ---------------------------------------------------------------
% colstruct
%
%  Date           Author            E-mail                      Version
%  25 Feb  2015   A.G. Balderrama   amanda.gaudreau@gmail.com      0
%  12 Mar  2015   A.G. Balderrama   amanda.gaudreau@gmail.com      1
%     Uses float corrected if available & nests all lum data into a single
%     cell

compdir = loaddirfun;

PropertyNames = varargin(1:2:length(varargin));
PropertyVal = varargin(2:2:length(varargin));

if strmatch('unmixed',PropertyNames)
  umxonly = PropertyVal{strmatch('unmixed',PropertyNames)};
else
  umxonly = 1;
end

[ivisfldr,ivisXL] = getfilelist(varargin{:});
XLhcols = ivisXLcol(ivisXL(1,:));
hdroffset = 2;
aninq = cellfun(@(x) num2str(x),ivisXL(hdroffset+1:end,XLhcols.animalnumcol),'uniformoutput',0);
evalf = 1:length(ivisfldr.ivisSUfldr);
nfldrs = length(evalf);

Iph = cell(nfldrs,1);
Lcell = cell(nfldrs,1);
Lcorrcell = cell(nfldrs,1);
for i = 1:nfldrs
  fprintf('Processing %s (%d of %d)......................................\n',ivisfldr.ivisSUfldr{evalf(i)},i,nfldrs);
  Ist = readIVISfiles('dat',ivisfldr.ivisSUfldr{evalf(i)},varargin{:});
  if isfield(Ist,'photograph') && isfield(Ist,'luminescent')
    clkannumrow = find(cellfun(@(x) ~isempty(x),strfind(Ist.ClickTable(:,2),'Animal Number')));
    clkannuminfo = Ist.ClickTable(clkannumrow,3:end);
    XLannum = num2str(ivisXL{i + hdroffset,XLhcols.animalnumcol});
    correctfldr = find(cellfun(@(x) ~isempty(x),strfind(clkannuminfo,XLannum)));
    
    %if ~isempty(correctfldr)
      Iph{i} = double(Ist.photograph.image{1});
      Lcell{i} = cellfun(@(x) double(x),Ist.luminescent.image,'uniformoutput',0);
      Lcorrcell{i} = cellfun(@(x) double(x),Ist.luminescent.imageCorrected,'uniformoutput',0);
    %else
    %  correctfldr = find(cellfun(@(x) ~isempty(x),clkannuminfo));
    %  fprintf('Orig sequence %s seems to be incorrectly relabeled as %s. Seq indicates %s while XL files indicates %s.\n',...
    %    Ist.origseq,Ist.dataset,clkannuminfo{correctfldr(1)},XLannum)
    %end
  else
    fprintf('Spectral Unmixing has not been conducted for %s (original sequence %s)\n',...
      ivisfldr.ivisSUfldr{evalf(i)}, ivisXL{evalf(i)+hdroffset,XLhcols.origivissuSEQ})
  end
end

qbucs = ivisXL(evalf+hdroffset,XLhcols.bucscols);
usebucs = cellfun(@(x) sum(~isnan(x)),qbucs) > 0;
evalrows = find(sum(usebucs,2));
evalcols = find(sum(usebucs,1));
bucsblk = cell(nfldrs,1);
bucs_sum = nan(nfldrs,1);
for i = evalrows(:)'
  bucsblk{i} = cell2mat(cellfun(@(x) str2num(x),qbucs(i,evalcols),'uniformoutput',0)');
  bucs_sum(i) = sum(bucsblk{i}(:));
end

IVISinfo.ivisfldr = ivisfldr;
IVISinfo.ivisXL = ivisXL;
IVISinfo.animalseval = aninq;
end