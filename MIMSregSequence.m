function MIMSregSequence(varargin)
%% Script: MIMSregSequence
% Description:
% Example:
% Required Functions: twoimgMIkde
% INPUTS ----------------------------------------------------------------
% varargin - 'PropertyName','PropertyValue'
%   • 'regdata':    string indicating the registration datafile name to use
%   for analysis [DEFAULT: none]
%       >>> MIMSregSequence('regd','MIMS30-1ch-P_ABA83-ni_Rightv2-20161117.mat')
%   • 'usemu0':  string of the registration datafile name to use the best
%   data from OR six-element vector [DEFAULT = uses mu0 from LvsR reg
%   experiment]
%   • 'expParam': structure containing the fields the use would like to
%   specify for a given experiment. If a field is not present, one is
%   created with the default values given below
%     DEFAULT FIELD NAMES & VALUES
%       exptag: ''
%       bndstr: 'large'
%    regminDim: []
%          eli: 1
%     regtoABA: 'nissl'
%       sumdim: 1
%     normcost: 1
%      normval: 5
%      tracksa: 1
%     tempInit: 5
%     tempFact: 0.9900
%     reannInt: 250
%      maxIter: 5000
%        nbins: 40
%      kernstr: 'norm'
%     costfstr: 'mi'
%   • 'MIMSstr': MIMS number or string associated with the dataset to be
%   analyzed
%       For initial experiment, the following options are recommended:
%         '20150923_MIMS30_iTBI'    '20150818_MIMS22_sham'
%         '20150324_MIMS12_sham'    '20150413_MIMS35_control'
%         '20150415_MIMS27_iTBI'    '20150818_MIMS14_iTBI'
%         '20150819_MIMS26_iTBI'
%   • 'back': logical indicating whether to count the background in the MI
%   cost (1) or not (0) [DEFAULT = 0]
%   • 'evalside': string indicating the side of the brain to perform
%   registration on (either 'l' or 'r')
%
%
% OUTPUTS ---------------------------------------------------------------
%
%
%  Date           Author            E-mail                      Version
%  16 Nov 2016   AGBalderrama      amandagbalderrama@gmail.com     0
%  21 Nov 2016   AGBalderrama      amandagbalderrama@gmail.com     1
%  23 Nov 2016   AGBalderrama      amandagbalderrama@gmail.com     1.1
%     Added coniditional statement which uses nearest neighbors for
%     resizing if atlas image is used 
%  26 Jan 2017   AGBalderrama      amandagbalderrama@gmail.com     1.2
%     Corrected part of code, didn't save old version


if exist('varargin')
  PropertyNames = lower(varargin(1:2:length(varargin)));
  PropertyVal = varargin(2:2:length(varargin));
else
  PropertyNames = [];
  PropertyVal = [];
end

%.....................................................................
%   simanneal_MI varargin
%.....................................................................

if strmatch('regd',PropertyNames)
  regDataName = PropertyVal{strmatch('regd',PropertyNames)};
else
  regDataName = '';
end

if strmatch('usemu0',PropertyNames)
  mu0v = PropertyVal{strmatch('usemu0',PropertyNames)};
else
  mu0v = [];
end

if strmatch('expp',PropertyNames)
  exParam = PropertyVal{strmatch('expp',PropertyNames)};
else
  exParam = [];
end

if strmatch('mims',PropertyNames)
  dstr = PropertyVal{strmatch('mims',PropertyNames)};
else
  dstr = 'MIMS30';
end


if strmatch('evalside',PropertyNames)
  evalside = lower(PropertyVal{strmatch('evalside',PropertyNames)});
else
  evalside = 'lr';
end

if strmatch('back',PropertyNames)
  cntbkgnd = PropertyVal{strmatch('back',PropertyNames)};
else
  cntbkgnd = 0;
end

if length(evalside) > 1
  mvals = 0:1;
elseif isequal('l',evalside)
  mvals = 0;
elseif isequal('r',evalside)
  mvals = 1;
else
  mvals = 0:1;
end


set(0,'DefaultAxesFontSize',14,'defaultfigurecolor','w',...%'DefaultAxesFontWeight','bold',...
  'DefaultLineLineWidth',2, 'Defaultaxesposition','remove');
lnmarkers = 'ospd^h*ospd^h*ospd^h*ospd^h*ospd^h*ospd^h*';
lnsty = {'-','-.',':','--','-','-.',':','--','-','-.',':','--','-','-.',':','--','-','-.',':','--','-','-.',':','--','-','-.',':','--'};
lncol = colorcube(64);%parula(64);%
% lncol = custom_colormaps('lines',30);
cdir = loaddirfun;

mimsInfo = MIMSsummary;
anvec = [84,46,6,15,17,45,47,37,36]; %[3,7,44];
mimslabels = cellfun(@(x) x(10:end),mimsInfo.dsetInfo(anvec,1),'uniformoutput',0);

    %LEFT ABAimg            %RIGHT ABAimg
abastrs = {{'ABA_coronal84_334','ABA_coronal83_330'};... %MIMS30
  {'ABA_coronal84_334', 'ABA_coronal83_330'};... %MIMS22
  {'ABA_coronal84_334', 'ABA_coronal83_330'};...%MIMS12
  {'ABA_coronal86_342', 'ABA_coronal84_334'};... %MIMS35
  {'ABA_coronal84_334', 'ABA_coronal84_334'};... %MIMS27
  {'ABA_coronal84_334', 'ABA_coronal82_326'};... %MIMS14
  {'ABA_coronal84_334', 'ABA_coronal84_334'};... %MIMS26
  {'ABA_coronal72_285', 'ABA_coronal72_285'};...%ALON10
  {'ABA_coronal75_297', 'ABA_coronal74_293'}}; %ALON12

for i = 1:length(mimslabels)
  mimsABAmatch.(mimslabels{i}) = abastrs{i};
end


if ~isempty(strfind(dstr,'ALON'))
  atlasColSeed = [125, 208, 75;... HPF
    103, 167, 59;...pyramidal cell layer of HPF
    29, 166, 153;... RETctx
    24, 128, 101;... SSctx
    1, 159, 171;... SENSctx
    13, 185, 176;... RHIctx
    108, 203, 186;... PIRctx
    170, 234, 211;... AMYctx
    255, 255, 75;... ctxsp
    127, 192, 228;...STRI
    205, 205, 204;... WMT
    150, 150, 150;... VENT
    255, 126, 132;... TH
    251, 75, 63]; % HypoThal
  chorder = {'Fe57','Mn','Zn64','Cu63','P_','Fe57','Gd'};
else
  atlasColSeed = [125, 207, 78;... HPF
    103, 167, 59;...pyramidal cell layer of HPF
    75, 125, 50;... SUB of HPF
    29, 166, 153;... RETctx
    8, 132, 140;... VISctx
    1, 159, 171;... SENSctx
    13, 185, 176;... RHIctx
    170, 234, 211;... AMYctx
    255, 255, 75;... ctxsp
    205, 205, 204;... WMT
    255, 126, 132;... TH
    254, 144, 255;... MidB
    251, 75, 63]; % HypoThal
  chorder = {'P_','Zn64','Cu63','La','Fe57','Gd'};
end
nseg = size(atlasColSeed,1);

if isempty(exParam) || ~isfield(exParam,'exptag')
  exParam.exptag = '';
end

if isempty(exParam) || ~isfield(exParam,'bndstr')
  exParam.bndstr = 'large';% 'regular';% 'precise';% 'rngFind';%
end

if isempty(exParam) || ~isfield(exParam,'eli')
  exParam.eli = [1];
end

if isempty(exParam) || ~isfield(exParam,'regminDim')
  regminDim = [];% 100;%   empty for full MIMS size, otherwise integer value of smallest dimension
else
  regminDim = exParam.regminDim;
end

if isempty(exParam) || ~isfield(exParam,'regtoABA')
  exParam.regtoABA = 'nissl';%'atlas';%
end

dfilesnum = mimsInfo.dfilesindx;
rotvec = mimsInfo.rotvec;
fileINX = 1:length(dfilesnum);
usefveci = find(cellfun(@(x) ~isempty(x),strfind(mimsInfo.dsetInfo(:,1),dstr)));
fvec = dfilesnum(usefveci);
evalf = mimsInfo.dsetInfo{usefveci,1}; % dfiles(fvec);%{'20150415_MIMS27_iTBI';'20150324_MIMS12_sham'};
an = strmatch(evalf(10:end),mimslabels);

% for an = 1%:length(mimslabels)
matname = [cdir.saMIMSbrain evalf(10:end) '.mat'];
if exist(matname,'file')
  load(matname); % loads dp, usedat, Jmat, chanInd, and I0 into the workspace
else
  [dp,usedat,Jmat,chanInd,I0] = loadMIMSdata(evalf,chorder,'IYmve');
  usedat.rotval = mimsInfo.rotvec(anvec(an));
  save(matname,'dp','usedat','Jmat','chanInd','I0');
end

Iuse = I0(:,:,exParam.eli);
LRstrs = {'Left','Right'};
switch exParam.bndstr
  case 'rngFind'
    exprng = repmat([10,10,6*pi/180,0.2,0.2,0.1],[2,1]);
  case 'large'
    exprng = repmat([5,5,4*pi/180,0.15,0.15,0.08],[2,1]);
  case 'regular'
    exprng = repmat([3,3,3*pi/180,0.08,0.08,0.06],[2,1]);
  case 'precise'
    exprng = repmat([1.5,1.5,2*pi/180,0.08,0.08,0.06],[2,1]);
  otherwise
    exprng = repmat([3,3,3*pi/180,0.08,0.08,0.06],[2,1]);
end

for m = mvals
  abaimgstr = mimsABAmatch.(dp.mveparams.dataset(10:end)){m+1};
  
  %------------- Read in appropriate atlas image
  % CREATES: ABAimg sized -atlas half mask, -atlas mask, -atlas integer
  % labeled image
  Iat = imread([cdir.fpath 'Brain Images\' abaimgstr 'AtlasOnly.tif']);
  Iatg = double(rgb2gray(Iat));
  Mat = Iatg <= 243;
  ccM0aba = bwconncomp(Mat);
  npxlcc = cellfun(@(x) length(x),ccM0aba.PixelIdxList);
  [~,cci] = max(npxlcc);
  Mat(:) = false;
  Mat(ccM0aba.PixelIdxList{cci}) = true;
  Mat = imclose(Mat,strel('square',5));
  Mat = imfill(Mat,'holes');
  rpMat = regionprops(Mat,'boundingbox','PixelIdxList');
  Mhalf0 = zeros(size(Mat));
  Mhalf0(:,floor(rpMat.BoundingBox(1)):end) = 1;
  
  Iatv = reshape(double(Iat(:)),[numel(Mat),3]);
  if cntbkgnd
    IatLabel = zeros(size(Mat));
  else
    IatLabel = nan(size(Mat));
  end
  IatLabel(Mat(:)) = kmeans(Iatv(Mat(:),:),nseg,'start',atlasColSeed);
  
  %------------- Read in appropriate nissl image
  Ini = imread([cdir.fpath 'Brain Images\' abaimgstr 'Nissl.jpg']);
  Inig = double(rgb2gray(Ini));
  Mabatsu =  Inig <= 243;
  ccMabatsu = bwconncomp(Mabatsu);
  npxlcc = cellfun(@(x) length(x),ccMabatsu.PixelIdxList);
  [~,cci] = max(npxlcc);
  Mabatsu(:) = false;
  Mabatsu(ccMabatsu.PixelIdxList{cci}) = true;
  Mabatsu = imclose(Mabatsu,strel('square',5));
  Mabatsu = imfill(Mabatsu,'holes');
  rpMni = regionprops(Mabatsu,'boundingbox','PixelIdxList');
  abaBBr2chalf = rpMat.BoundingBox(4)/rpMat.BoundingBox(3);
  abaBBr2cfull = rpMni.BoundingBox(4)/rpMni.BoundingBox(3);
  
  %------------- Crop all necessary images to tight region around tissue
  Macrp = Mabatsu([1:rpMni.BoundingBox(4)]+(rpMni.BoundingBox(2)+0.5),...
    [1:rpMni.BoundingBox(3)]+(rpMni.BoundingBox(1)+0.5));
  Mhalf0 = Mhalf0([1:rpMni.BoundingBox(4)]+(rpMni.BoundingBox(2)+0.5),...
    [1:rpMni.BoundingBox(3)]+(rpMni.BoundingBox(1)+0.5));
  Mat = Mat([1:rpMni.BoundingBox(4)]+(rpMni.BoundingBox(2)+0.5),...
    [1:rpMni.BoundingBox(3)]+(rpMni.BoundingBox(1)+0.5));
  
  if strmatch(exParam.regtoABA,'nissl')
    mifname = [abaimgstr 'Nissl.jpg'];
    IabaReg = Inig([1:rpMni.BoundingBox(4)]+(rpMni.BoundingBox(2)+0.5),...
      [1:rpMni.BoundingBox(3)]+(rpMni.BoundingBox(1)+0.5));
    binPos = [];
  elseif strmatch(exParam.regtoABA,'atlas')
    mifname = [abaimgstr 'AtlasOnly.tif'];
    IabaReg = IatLabel([1:rpMni.BoundingBox(4)]+(rpMni.BoundingBox(2)+0.5),...
      [1:rpMni.BoundingBox(3)]+(rpMni.BoundingBox(1)+0.5));
    binPos = 1:nseg;
  end
  
  if cntbkgnd
    Mat = ones(size(Mat,1),size(Mat,2));
  end
  IabaReg(Mat == 0) = nan;
  abaAR = size(Macrp,1)/size(Macrp,2);
  
  if m == 0
    IabaReg = fliplr(IabaReg);
    IatLabel = fliplr(IatLabel);
    Macrp = fliplr(Macrp);
    Mhalf0 = fliplr(Mhalf0);
    Mat = fliplr(Mat);
    rpMni = regionprops(fliplr(Mabatsu),'boundingbox','PixelIdxList');
  end
  
  %------------- Create MIMS mask
  Mmimstsu = imdilate(imfill(dp.mveparams.tsuMask,'holes'),strel('disk',3));
  rpMm = regionprops(Mmimstsu,'boundingbox','PixelIdxList');
  Mmcrp0 = Mmimstsu([1:rpMm.BoundingBox(4)]+(rpMm.BoundingBox(2)+0.5),...
    [1:rpMm.BoundingBox(3)]+(rpMm.BoundingBox(1)+0.5));
  mimsAR = size(Mmcrp0,1)/size(Mmcrp0,2);
  mimsRawSize = size(Mmcrp0);
  
  if mimsAR < abaAR*0.75
    rsize = [size(Mmcrp0,1),size(Mmcrp0,1)/abaAR];
  else
    rsize = size(Mmcrp0);
  end
  Macrp = imresize(double(Macrp),rsize,'nearest');
  Mmcrp = imresize(double(Mmcrp0),rsize,'nearest');
  
  siderng = repmat([25,25,25*pi/180,0.25,0.25,0.4],[2,1]); 
  mu0 = [0,0,0,1,1,0];
  bnd1 = repmat(mu0,[2,1]) + [ones(1,length(mu0));-1.*ones(1,length(mu0))].*siderng;
  
  %------------- Run registration function & aggregate variables
  sideMname = [matname(1:end-4) '_' LRstrs{m+1} '.mat'];
  
  if exist(sideMname,'file')
    load(sideMname);
  else
    clk1 = datenum(datetime('now'));
    [mu_tst,exp_outputs,phiMI] = simannealGeneral(Mmcrp,Macrp,'bounds',bnd1,......'Fmask',FMask,'Mmask',MMask,
      'costf','max1','pval',1,'tempfact',0.98,'t0',10,'mu0',mu0,...
      'reann',150,'maxiter',5000,...
      'normcost',1,'normto',10,'track_accept',0);%,'Bbins',{binPos}); %);%

%       'reann',50,'maxiter',1500,...
%       'normcost',1,'normto',150,'track_accept',0);%,'Bbins',{binPos}); %);%
    clk2 = datenum(datetime('now'));
    
    disp(sprintf('%s %s side mask found, Minutes Elapsed = %1.2f',mimslabels{an},LRstrs{m+1},(clk2-clk1)*24*60))
    
    mu_tst = change_alpha_dim(mu_tst,rsize,mimsRawSize);
    
    Mhalf = imresize(double(Mhalf0),mimsRawSize,'nearest');
    Mhalf = padarray(Mhalf,[50,50],'replicate');% zeros(size(Mhalf,1)+20,size(Mhalf,2)+20);
    Mhalftx = transform_image(Mhalf,mu_tst,'interp','nearest','extrapval',1);
    Mhalftx = Mhalftx(51:end-50,51:end-50);
    m_sideMask = and(Mhalftx,Mmcrp0);
    m_sideMask = m_sideMask + Mmcrp0;
    Mmimstsu = double(Mmimstsu);
    Mmimstsu(:) = 0;
    Mmimstsu([1:rpMm.BoundingBox(4)]+(rpMm.BoundingBox(2)+0.5),...
      [1:rpMm.BoundingBox(3)]+(rpMm.BoundingBox(1)+0.5)) = m_sideMask;
    save(sideMname,'Mmimstsu','mu_tst');
  end
  
  exParam = initRegParam(exParam);
  %------------- Prepared fixed image data
  Mmims = Mmimstsu == 2;
  ccuse = bwconncomp(Mmims);
  [~,icc] = max(cellfun(@(x) length(x),ccuse.PixelIdxList));
  Mmims = false(size(Mmims));
  Mmims(ccuse.PixelIdxList{icc}) = true;
  rpuse = regionprops(Mmims,'BoundingBox','PixelIdxList');
  Imims0 = Iuse([1:rpMm.BoundingBox(4)]+(rpMm.BoundingBox(2)+0.5),...
    [1:rpMm.BoundingBox(3)]+(rpMm.BoundingBox(1)+0.5),:);
  Mmims = Mmims([1:rpMm.BoundingBox(4)]+(rpMm.BoundingBox(2)+0.5),...
    [1:rpMm.BoundingBox(3)]+(rpMm.BoundingBox(1)+0.5));
  MmimsNAN = double(Mmims);
  MmimsNAN(Mmims == 0) = nan;
  if cntbkgnd
    MmimsNAN = ones(size(MmimsNAN,1),size(MmimsNAN,2));
  end
  Imims0 = Imims0.*repmat(MmimsNAN,[1,1,size(Imims0,3)]);
  
  %------------- Determine final registration parameters (image sizes,
  %transformation bounds)
  if isempty(regminDim)
    RegSize = size(Mmims);
  elseif length(regminDim) == 1 && abaBBr2chalf > 1
    RegSize = round([abaBBr2chalf,1].*regminDim);%size(FMask); %
    RegSize = round([RegSize(1),RegSize(1)/abaBBr2cfull]);
  elseif length(regminDim) == 1 && abaBBr2chalf < 1
    RegSize = round([1,1/abaBBr2chalf].*regminDim);%size(FMask); %
    RegSize = round([RegSize(1),RegSize(1)/abaBBr2cfull]);
  elseif length(regminDim) == 2
    RegSize = regminDim;
  end
  
  ImimsReg = zeros([RegSize size(Imims0,3)]);
  for ch = 1:size(Imims0,3)
    ImimsReg(:,:,ch) = imresize(Imims0(:,:,ch),RegSize);
  end
  
  if strmatch(exParam.regtoABA,'nissl')
    IabaReg = imresize(IabaReg,RegSize);%,'nearest');
  else
    IabaReg = imresize(IabaReg,RegSize,'nearest');
  end
  
  %   figure; subplot(2,2,1); imagesc(FImg); colormap(gray);
  %   subplot(2,2,2); imagesc(FMask); colormap(gray);
  %   subplot(2,2,3); imagesc(Muse); colormap(gray);
  %   subplot(2,2,4); imagesc(MMaskuse); colormap(gray);
  
  if isstr(mu0v) && exist([cdir.saMIMSbrain mu0v])
    if m == 0
      mu0v = strrep(mu0v,'Right','Left');
    else
      mu0v = strrep(mu0v,'Left','Right');
    end
    abai = strfind(mu0v,'ABA');
    mu0v(abai:abai+4) = ['ABA' abaimgstr(12:13)];
    Rmu0 = load([cdir.saMIMSbrain mu0v]);
    mu0 = change_alpha_dim(Rmu0.RegStruct.FdimBestmu,mimsRawSize,RegSize);
    clear Rmu0
  elseif ~isempty(mu0v) && isnumeric(mu0v)
    mu0 = mu0v;
  else
    mu0 = change_alpha_dim(mu_tst,mimsRawSize,RegSize);
  end
  
  bnd1 = repmat(mu0,[2,1]) + [ones(1,length(mu0));-1.*ones(1,length(mu0))].*exprng;
  
  %------------- Run registration function & aggregate variables
  channelsUsed = {dp.mveparams.elnames{chanInd(exParam.eli)}};
  elstrs = cell2mat(cellfun(@(x) strrep(x(1:2),'_',''),channelsUsed,'uniformoutput',0));
  if isempty(regDataName)
    expName = sprintf('%s-%dch-%s_ABA%s-%s_%s%s',...
      strtok(evalf(10:end),'_'),size(ImimsReg,3),elstrs,abaimgstr(12:13),exParam.regtoABA(1:2),...
      LRstrs{m+1},exParam.exptag);
    baseFigName = [expName '-' cdir.dateyyyymmddstr];
    RegMname = [baseFigName '.mat'];
  else
    RegMname = regDataName;
    baseFigName = RegMname(1:end-4);
  end
  RegMname = [cdir.saMIMSbrain RegMname];
  
  if exist(RegMname,'file')
    load(RegMname);
  else
    clk1 = datenum(datetime('now'));
    [mu_tst,exp_outputs,phiMI] = simannealGeneral(ImimsReg,IabaReg,...
      'normcost',exParam.normcost,'normto',exParam.normval,'bounds',bnd1,......'Fmask',FMask,'Mmask',MMask,
      'costf',exParam.costfstr,'tempfact',exParam.tempFact,'t0',exParam.tempInit,...
      'mu0',mu0,'reann',exParam.reannInt,'maxiter',exParam.maxIter,...
      'track_accept',exParam.tracksa,'nbins',exParam.nbins,...
      'kernel',{exParam.kernstr},'sumdim',exParam.sumdim,'Bbins',{binPos}); %);%
    clk2 = datenum(datetime('now'));
    
    expDetails.costf = exParam.costfstr;
    expDetails.sumdim = exParam.sumdim;
    expDetails.tempfact = exParam.tempFact;
    expDetails.initTemp = exParam.tempInit;
    expDetails.bnds = bnd1;
    expDetails.reannInt = exParam.reannInt;
    expDetails.maxIter = sum(exp_outputs.fval ~= 0);
    expDetails.nbins = exParam.nbins;
    expDetails.kernelstr = exParam.kernstr;
    expDetails.regSize = RegSize;
    expDetails.date = cdir.dateyyyymmddstr;
    expDetails.runTimeMin = (clk2 - clk1)*24*60;
    
    FimgInfo.fname = dp.mveparams.dataset;
    FimgInfo.channelsUsed = channelsUsed;
    FimgInfo.dataType = usedat.type;
    FimgInfo.ImgUsed = ImimsReg;
    FimgInfo.RawImg = Iuse;
    FimgInfo.RawImgSize = dp.mveparams.dsize(1:2);
    FimgInfo.Mask = Mmims;
    FimgInfo.AxSymData = [];%axparam; --- Old fields
    FimgInfo.AxSymSlope = [];%bestSlope0; --- Old fields
    FimgInfo.AxSymInt = [];%bestInt0; --- Old fields
    FimgInfo.imgInd = rpMm.PixelIdxList;
    FimgInfo.imgRows = [1:rpMm.BoundingBox(4)]+(rpMm.BoundingBox(2)+0.5);
    FimgInfo.imgCols = [1:rpMm.BoundingBox(3)]+(rpMm.BoundingBox(1)+0.5);
    
    MimgInfo.fname = mifname;
    MimgInfo.channelsUsed = 'rgb2gray';
    MimgInfo.fliplrUsed = abs(m-1);
    MimgInfo.RawImgSize = size(Inig);
    MimgInfo.Mask = Mat;
    MimgInfo.LabelImg = IatLabel;
    MimgInfo.imgInd = rpMat.PixelIdxList;
    MimgInfo.imgRows = [1:rpMni.BoundingBox(4)]+(rpMni.BoundingBox(2)+0.5);
    MimgInfo.imgCols = [1:rpMni.BoundingBox(3)]+(rpMni.BoundingBox(1)+0.5);
    
    mup = change_alpha_dim(mu_tst,expDetails.regSize,[length(FimgInfo.imgRows),length(FimgInfo.imgCols)]);
    RegStruct.expName = expName;
    RegStruct.bestmu = mu_tst;
    RegStruct.FdimBestmu = mup;
    RegStruct.outputs = exp_outputs;
    RegStruct.expDetails = expDetails;
    RegStruct.FimgInfo = FimgInfo;
    RegStruct.MimgInfo = MimgInfo;
    
    movefile(strcat(cdir.sapath,cdir.dateyyyymmddstr,'_nvar6Paccepth1_exp1.csv'),...
      strcat(cdir.sapath,sprintf('%s.csv',baseFigName)))
    
    set(gcf,'position',[169, 70, 1467, 894],...[56   191   934   772],...
      'name',[baseFigName '_SAresults'])
    
    save(RegMname,'RegStruct');
  end
  
  plotOverlays(RegStruct,baseFigName);
  
  %------------- conduct regional analysis of MIMS images
  if ~isfield(RegStruct,'perRegionData')
    anaSize = [length(RegStruct.FimgInfo.imgRows),length(RegStruct.FimgInfo.imgCols)];
    muRescale = anaSize./RegStruct.expDetails.regSize;
    mu_tst = RegStruct.outputs.bestx(:,RegStruct.expDetails.maxIter)'./sastatescalefactors;
    mu_tst(1) = mu_tst(1)*muRescale(2);
    mu_tst(2) = mu_tst(2)*muRescale(1);
    Ilabel = RegStruct.MimgInfo.LabelImg(RegStruct.MimgInfo.imgRows,RegStruct.MimgInfo.imgCols);
    Ilabel = imresize(Ilabel,anaSize,'nearest');
    Ilabel = transform_image(Ilabel,mu_tst,'interp','nearest');
    
    dataTypes = {'Iraw','IYmve','IXshift'};
    labelvals = unique(Ilabel(~isnan(Ilabel)))';
    if ~isempty(strfind(dstr,'ALON'))
      tabvars = {'HIPP' 'HIPPpy' 'RETctx' 'SSctx' 'SENSctx' 'RHIctx' 'PIRctx' 'AMYctx' 'ctxsp' 'STRI' 'WMT' 'VENT' 'THAL' 'HYTH'};% 'VariableNames',{'Gender' 'Age' 'State' 'Vote'}
    else
      tabvars = {'HIPP' 'HIPPpy' 'SUB' 'RETctx' 'VISctx' 'SENSctx' 'RHIctx' 'AMYctx' 'ctxsp' 'WMT' 'THAL' 'MID' 'HYTH'};% 'VariableNames',{'Gender' 'Age' 'State' 'Vote'}
    end
    evalEl = dp.mveparams.elnames;
    perRegionData.ElementNames = evalEl;
    perRegionData.LabelNames = tabvars;
    perRegionData.LabelValues = labelvals;
    for t = 1:3
      meantmat = nan(numel(evalEl),numel(labelvals));
      medtmat = nan(numel(evalEl),numel(labelvals));
      stdtmat = nan(numel(evalEl),numel(labelvals));
      for el = 1:length(evalEl)% relevantEli
        for l = labelvals(:)'
          elI = dp.(dataTypes{t}).(dp.mveparams.elnames{el});
          elI(RegStruct.FimgInfo.Mask == 2) = nan;
          elI = elI(RegStruct.FimgInfo.imgRows,RegStruct.FimgInfo.imgCols);
          meantmat(el,find(labelvals == l)) = nanmean(elI(Ilabel == l));
          medtmat(el,find(labelvals == l)) = nanmedian(elI(Ilabel == l));
          stdtmat(el,find(labelvals == l)) = nanstd(elI(Ilabel == l));
        end
      end
      perRegionData.(dataTypes{t}).mean = meantmat;
      perRegionData.(dataTypes{t}).median = medtmat;
      perRegionData.(dataTypes{t}).std = stdtmat;
    end
    RegStruct.perRegionData = perRegionData;
  end
  
  %----------------- Plots results -------
  plotData = 'IXshift';
  plotSignal = 'Median';
  nel = length(RegStruct.perRegionData.ElementNames);
  nlab = length(RegStruct.perRegionData.LabelValues);
  
  figure('name',[baseFigName '_LabelvMed'],'position',[194         392        1011         527]);
  elLines = plot(RegStruct.perRegionData.LabelValues,RegStruct.perRegionData.(plotData).(lower(plotSignal))','-','marker','o');
  legend(strrep(RegStruct.perRegionData.ElementNames,'_','\_'));
  set(gca,'xtick',RegStruct.perRegionData.LabelValues,...
    'xticklabel',RegStruct.perRegionData.LabelNames,'fontsize',13);
  xlim([min(RegStruct.perRegionData.LabelValues),max(RegStruct.perRegionData.LabelValues)])
  set(elLines,{'color'},mat2cell(lncol(1:nel,:),ones(1,nel),3),...
    {'markeredgecolor'},mat2cell(lncol(1:nel,:),ones(1,nel),3),...
    {'markerfacecolor'},mat2cell(lncol(1:nel,:),ones(1,nel),3));
  xlabel('Isotope'); ylabel([plotData ' level'])
  title(sprintf('%s amount of isotopes (lines) within brain regions',plotSignal))
  grid on;
  
  figure('name',[baseFigName '_ElvMed'],'position',[194         392        1011         527]);
  typeLines = plot(RegStruct.perRegionData.(plotData).(lower(plotSignal)),'-','marker','o');
  legend(RegStruct.perRegionData.LabelNames);
  set(gca,'xtick',1:nel,'xticklabel',strrep(RegStruct.perRegionData.ElementNames,'_',''));
  set(typeLines,{'color'},mat2cell(lncol(1:nlab,:),ones(1,nlab),3),...
    {'markeredgecolor'},mat2cell(lncol(1:nlab,:),ones(1,nlab),3),...
    {'markerfacecolor'},mat2cell(lncol(1:nlab,:),ones(1,nlab),3));
  xlabel('Isotope'); ylabel([plotData ' level'])
  title(sprintf('%s of isotopes within indicated brain region (lines)',plotSignal))
  grid on;
  
  save(RegMname,'RegStruct');
  
  cSA = sastatescalefactors;
  bMu = RegStruct.outputs.bestx;
  bMui = find(sum(diff(bMu'),2) ~= 0);
  param_abb = {'t_x','t_y', '\theta', 's_x', 's_y', 'k'};
  figure('name',[baseFigName '_' upper(exParam.costfstr) 'perState'],'position',[94, 153, 1440, 722]);
  linecolrs = lines(6);
  gcol = gray(length(bMui)+6);
  for i = 1:6
    subplot(2,3,i);
    
    plot(RegStruct.outputs.x(i,RegStruct.outputs.fval ~= 0)./cSA(i),abs(RegStruct.outputs.fval(RegStruct.outputs.fval ~= 0)),'.',...
      'markeredgecolor',linecolrs(i,:),'markerfacecolor',linecolrs(i,:));
    axis tight; grid on;
    hold on;
    
    for g = 1:length(bMui)
      plot(RegStruct.outputs.bestx(i,bMui(g))./cSA(i),abs(RegStruct.outputs.fval(bMui(g))),'s',...
        'markeredgecolor',gcol(g,:),'markerfacecolor',gcol(g,:));
    end
    %set(pk,{'markeredgecolor'},mat2cell(gray(length(bMui)),ones(1,length(bMui)),3),...
    %  {'markerfacecolor'},mat2cell(gray(length(bMui)),ones(1,length(bMui)),3));
    title(param_abb{i})
  end
  
  save_open_figures(cdir.sapath,[],[],''); close all;
  save([cdir.saMIMSbrain  baseFigName '.mat'],'RegStruct');
end
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function  [dmve,usedat,Imat,relevantEli,Inorm] = loadMIMSdata(dstr,relevantEl,mvetype);
    useElstr = [];%{'P_','Fe','Cu','C_13','Zn64','Gd'};%
    anatomystr = 'HIPP';
    anatomystr2 = 'MIMS';
    
    [df,hdrf] = readfiles('dat',dstr);
    [df,hdrf,d_div,line_types,type_rows] = divMIMSstruct(df,hdrf);
    usetype = find(cellfun(@(x) ~isempty(x),strfind(line_types,anatomystr)));
    usetype2 = find(cellfun(@(x) ~isempty(x),strfind(line_types,anatomystr2)));
    if isempty(usetype)
      usetype = usetype2;
    end
    df = d_div{usetype};
    elCh = fieldnames(rmfield(df,skip_fields));
    for eldf = 1:length(elCh)
      nanmask = isnan(df.(elCh{eldf}));
      df.(elCh{eldf})(nanmask) = quantile(df.(elCh{eldf})(:),0.1);
    end
    hdrf = hdrf([1;type_rows{usetype}+1],:);
    df.dataset = dstr;
    usedat.MF = 0;
    dmve = MIMSmve(df,useElstr,'hdrtxt',hdrf,'outtag','rXY','medfilter',usedat.MF,'repback',0);%,'rotdat',0);
    
    usedat.type = mvetype;
    usedat.rotreq = 0;
    if usedat.MF
      usedat.datadiscrip = strcat('_',anatomystr,usedat.type(1:2),'medfilt');
    else
      usedat.datadiscrip = strcat('_',anatomystr,usedat.type);
    end
    
    switch usedat.type
      case 'Iraw'
        usedat.datarng = 'auto';
      case 'IYmve'
        usedat.datarng = [-3,3];
      case 'IXshift'
        usedat.datarng = [0,255];
    end
    
    relevantEli = interp_data_elemrng(dmve.(usedat.type),relevantEl);  %{'Fe','mn','cu63'});%[1,2,3,4,6];
    % % c2mJ = cell2mat(cellfun(@(x) dmve.(usedat.type).(x),dmve.mveparams.elnames(useEli),'uniformoutput',false)');
    Jcell = cellfun(@(x) dmve.(usedat.type).(x),dmve.mveparams.elnames(relevantEli),'uniformoutput',false)';
    Imat = reshape(cell2mat(Jcell),...
      [dmve.mveparams.dsize(1:2),length(relevantEli)]);
    
    Inorm = Imat;
    Inorm(isnan(Inorm)) = 0;
    Inorm(Inorm < usedat.datarng(1)) = usedat.datarng(1);
    Inorm(Inorm > usedat.datarng(2)) = usedat.datarng(2);
    Inorm = Inorm + abs(usedat.datarng(1)); % ensures image is all positive
    Inorm = Inorm./max(Inorm(:)); %normalizes image so between 0 and 1
  end


  function regParam = initRegParam(regParam)
    %------------- Initialize all parameters -----
    if ~isfield(regParam,'exptag'); regParam.exptag = ''; end
    if ~isfield(regParam,'sumdim'); regParam.sumdim = 1; end
    if ~isfield(regParam,'normcost'); regParam.normcost = 1; end
    if ~isfield(regParam,'normval'); regParam.normval = 5; end
    if ~isfield(regParam,'tracksa'); regParam.tracksa = 1; end
    if ~isfield(regParam,'tempInit'); regParam.tempInit = 5; end
    if ~isfield(regParam,'tempFact'); regParam.tempFact = 0.99; end
    if ~isfield(regParam,'reannInt'); regParam.reannInt = 250; end
    if ~isfield(regParam,'maxIter'); regParam.maxIter = 5e3; end
    if ~isfield(regParam,'nbins'); regParam.nbins = 40; end
    if ~isfield(regParam,'kernstr'); regParam.kernstr = 'norm'; end
    if ~isfield(regParam,'costfstr'); regParam.costfstr = 'mi'; end%'max1';
  end

  function plotOverlays(RS,figstr)
    Ik = RS.MimgInfo.LabelImg(RS.MimgInfo.imgRows,RS.MimgInfo.imgCols);
    Ikp = imresize(Ik,[RS.expDetails.regSize],'nearest');
    Ikpr = transform_image(Ikp,RS.bestmu,'interp','nearest');
    Ikpr(isnan(Ikpr)) = 0;
    
    Fimg = RS.FimgInfo.ImgUsed;
    nanI = sum(isnan(Fimg),3) < size(Fimg,3);
    plotCols = find(sum(nanI) > 0);
    plotRows = find(sum(nanI') > 0);
    Fimg(isnan(Fimg)) = 0;
    Fimg = Fimg(plotRows,plotCols,:);
    Ikpr = Ikpr(plotRows,plotCols);
    
    imoverlay(Fimg,'',Ikpr,Ikpr,'','displaymode','edgeoverlay','edgecolor','c','npanels',1,'edgeth',0,'colormap','gray')
    % imoverlay(Fimg,'',Inipr,Inipr,'','displaymode','edgeoverlay','edgecolor','c','npanels',1,'edgeth',10,'colormap','gray')
    % imoverlay(Ikpr,'',Fimg,Fimg,'','displaymode','rgbgray','edgecolor','r','npanels',1,'edgeth',0,'colormap','parula')
    % imoverlay(Fimg,'',Ikpr,Ikpr,'','displaymode','rgbgray','edgecolor','r','npanels',1,'edgeth',0,'colormap','parula','transparency',0.5)
    %     imoverlay(repmat(Fimg,[1,1,3]),'',Ikpr,Ikpr,'','displaymode','rgbgray','edgecolor','r','npanels',1,'edgeth',0,'colormap','parula','transparency',0.5)
    title(sprintf('%s image with ABA%s edges: Reg Dim [%d %d]',...
      strtok(RS.FimgInfo.fname(10:end),'_'),RS.MimgInfo.fname(12:13),...
      RS.expDetails.regSize))
    set(gcf,'name',[figstr 'OverlayRegDim'])
   
    showFsize = [length(RS.FimgInfo.imgRows), length(RS.FimgInfo.imgCols)];
    %     Fsc = showFsize./RS.expDetails.regSize;
    %     muBest = RS.bestmu;
    %     % [~,AmuBest] = transform_image(checkerboard(100),muBest);
    %     muBest(1) = muBest(1)*Fsc(2);
    %     muBest(2) = muBest(2)*Fsc(1);
    %     %muBest(4) = muBest(4)*Fsc(2);
    %     %muBest(5) = muBest(5)*Fsc(1);
    Fimg = RS.FimgInfo.RawImg(RS.FimgInfo.imgRows,RS.FimgInfo.imgCols,:);
    Ikp = imresize(Ik,[showFsize],'nearest');
    Ikpr = transform_image(Ikp,RS.FdimBestmu','interp','nearest');
    Ikpr(isnan(Ikpr)) = 0;
    plotCols = find(sum(RS.FimgInfo.Mask) > 0);
    plotRows = find(sum(RS.FimgInfo.Mask') > 0);
    Fimg = Fimg(plotRows,plotCols,:);
    Ikpr = Ikpr(plotRows,plotCols);
    
    imoverlay(Fimg,'',Ikpr,Ikpr,'','displaymode','edgeoverlay','edgecolor','c','npanels',1,'edgeth',0,'colormap','gray')
    title(sprintf('%s image with ABA%s edges: MIMS Dim [%d %d]',...
      strtok(RS.FimgInfo.fname(10:end),'_'),RS.MimgInfo.fname(12:13),...
      length(RS.FimgInfo.imgRows), length(RS.FimgInfo.imgCols)))
    set(gcf,'name',[figstr 'OverlayMIMSDim'])
  end
end