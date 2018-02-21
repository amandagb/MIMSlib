function [Ipdf, pdfstruct] = EB_brainimg_classifyimgs(pibarc_struct,animals,Imasks,varargin)
%% Script: EB_brainimg_classifyimgs
% Description:
% INPUTS ----------------------------------------------------------------
% pibarc_struct: use either pdfkde or pdfhist from EB_brainimag_processing
% varargin - 'PropertyName','PropertyValue'
%   • 'prior':   string indicating what class probability to use [DEFAULT = 'uni']
%       >'uni': uniform
%       >'guess': informed guess about clas probabilities
%   • 'plotimgs': logical indicating whether to plot image classifications
%       (1) or not (0) [DEFAULT = 1]
%   • 'saveimgs': logical indicating whether to save classification images
%       (1) or not (0) in the Latex path [DEFAULT = 0]
%   • 'saveclassI': string indicating folder path where the *.mat files of
%       the binary masks should be saved [DEFAULT = 0] -- if this
%       PropertyName is detected at all, it's assumed that the user would
%       like to save class binary images and an empty str values indicates
%       to save in the Latex path
%   • 'usepdf': string indicating which pdf should be used to determine
%        classification. Whatever pdf is indicated, the classification is
%        done by finding the maximum value per row and assigning the class
%        column corresponding to the max
%       >'p(i|c)': likelihood
%       >'p(i,c)': joint probability [DEFAULT]
%       >'p(c|i)': posterior -- most reasonable to use
%   • 'evalanimal': vector indicating animal identification number to
%   analyze [DEFAULT = ALL in animals.id]
%   • 'evalview': string indicating views to create classification maps for
%   [DEFAULT = 'DLRV']
%
% OUTPUTS ---------------------------------------------------------------
%
%
%  Date           Author              E-mail                      Version
%  17 Oct  2014   Amanda Balderrama   amanda.gaudreau@gmail.com     0
%     Reads in data generated on 20141016 and creates a probability image
%     for the entire image (uses 5 classes: bkgnd, bld, EB, hemo, tsu
%  19 Oct  2014   Amanda Balderrama   amanda.gaudreau@gmail.com     1
%  30 Oct  2014   Amanda Balderrama   amanda.gaudreau@gmail.com     2
%   9 Nov  2014   Amanda Balderrama   amanda.gaudreau@gmail.com     3
%   4 Dec  2014   Amanda Balderrama   amanda.gaudreau@gmail.com     4
%     Determined that using the p(i) built from previous experiments leads
%     to erroneous results since pi ~= sum over c {p(i,c)}, so p_imarg is
%     not used

infocusimgs = {'20140819_pad2D_n08IPEB4','20140819_pad_IVIS2L_EB_pol0','20140819_pad2R_n16IPEB4','20140819_pad2V_n14IPEB4';...
  '20140819_pad3D_n08IPEB4','20140819_pad_IVIS3L_EB_pol0','20140819_pad3R_n12IPEB4','20140819_pad3V_n07IPEB4';...
  '20140819_pad4D_n06IPEB4','20140819_pad_IVIS4L_EB_pol0','20140819_pad4R_n07IPEB4','20140819_pad4V_n02IPEB4';...
  '20140819_nopad5D_n02IPEB4','20140819_nopad_IVIS5L_EB_pol0','20140819_nopad5R_n07IPEB4','20140819_nopad5V_n08IPEB4';...
  '20140819_nopad6D_n07IPEB4','20140819_nopad_IVIS6L_EB_pol0','20140819_nopad6R_n05IPEB4','20140819_nopad_IVIS6V_EB_pol0';...
  '20140821_nopad8D_n08IPEB4','20140821_nopad_IVIS8L_EB_pol0','20140821_nopad8R_n09IPEB4','20140821_nopad8V_n01IPEB4';...
  '20140821_nopad9D_n06IPEB4','20140821_nopad_IVIS9L_EB_pol0','20140821_nopad9R_n08IPEB4','20140821_nopad9V_n04IPEB4';...
  '20140821_nopad10D_n03IPEB4','20140821_nopad_IVIS10L_EB_pol0','20140821_nopad10R_n06IPEB4','20140821_nopad10V_n02IPEB4';...
  '20140821_pad12D_n03IPEB4','20140821_pad_IVIS12L_EB_pol0','20140821_pad_IVIS12R_EB_pol0','20140821_pad12V_n02IPEB4';...
  '','20141118_pad_IVIS12L_EBAu_pol2','',''};

%=========================================================================
% Initialize folder & File Variables
%=========================================================================
comppaths = loaddirfun;
if ~exist('pibarc_struct') || isempty(pibarc_struct)
  load(strcat(comppaths.dbEBexppath,'20141109_workspace.mat'),'pdfkde','animals','Imasks');
  pibarc_struct = pdfkde;
  clear pdfkde
end
startvin = 1;
animals.imgnames = infocusimgs;

PropertyNames = varargin(startvin:2:length(varargin));
PropertyVal = varargin((startvin+1):2:length(varargin));

if strmatch('prior',PropertyNames)
  pcstr = PropertyVal{strmatch('prior',PropertyNames)};
else
  pcstr = 'uni';
end

if strmatch('plot',PropertyNames)
  plot_imgs = PropertyVal{strmatch('plot',PropertyNames)};
else
  plot_imgs = 1;
end

if strmatch('saveimgs',PropertyNames)
  save_imgs = PropertyVal{strmatch('saveimgs',PropertyNames)};
else
  save_imgs = 0;
end

if strmatch('saveclass',PropertyNames)
  saveMpath = PropertyVal{strmatch('saveclass',PropertyNames)};
  if isempty(saveMpath)
    saveMpath = strcat(comppaths.mtngpath,'\EBprocImgs');
  end
else
  saveMpath = 0;
end

if strmatch('usepdf',PropertyNames)
  usepdf = PropertyVal{strmatch('usepdf',PropertyNames)};
else
  usepdf = 'p(i,c)';
end

if strmatch('evalanimal',PropertyNames)
  animalnumvec = PropertyVal{strmatch('evalanimal',PropertyNames)};
else
  animalnumvec = animals.id;
end

if strmatch('evalview',PropertyNames)
  anatomy_str = PropertyVal{strmatch('evalview',PropertyNames)};
else
  anatomy_str = 'L';%animals.anatomy_str;
end

class_str = fieldnames(pibarc_struct);
nclasses = length(class_str);
full_class_str = {'Blood','EB','Contusion','Tissue'};
class_color = [1,0,0;0,0,1;1,2/3,0 ;0,1,0];
tstdata = getfield(pibarc_struct,class_str{1});
lenx = length(tstdata(:));
p_ibarc = zeros(lenx,nclasses);
for c = 1:nclasses;
  cdata = getfield(pibarc_struct,class_str{c});
  p_ibarc(:,c) = cdata(:);
end

% >>> Set the probability of each class {blood, evans blue, contusion, tissue}
switch pcstr
  case 'uni'
    p_c = 1/nclasses .* ones(1,nclasses); % probability of classes
  case 'guess'
    p_c = [0.14,0.22,0.14,0.5];
end
p_c = p_c./sum(p_c);

% >>> Construct the empirical p(i)
p_ic = p_ibarc.*repmat(p_c(:)',lenx,1);
if strmatch('pi',PropertyNames)
  p_iemp = PropertyVal{strmatch('pi',PropertyNames)};
else
  p_iemp = sum(p_ic,2);
end

% >>> required if kde is not used 
non0p_ibarc = find(full(any(p_ibarc,2)) == 1);
non0rgbv = zeros(length(non0p_ibarc),3);
[non0rgbv(:,1),non0rgbv(:,2),non0rgbv(:,3)] = ind2sub([256,256,256],non0p_ibarc);

infocusimgs = animals.imgnames;

plotimgs = 1;

for an = 5%1:length(animalnumvec) % for each animal number
  ana_animalnum = animalnumvec(an);
  i = find(animals.id == ana_animalnum);
  animal_masks = getfield(Imasks,sprintf('an%d',ana_animalnum));
  
  dateyyyymmdd = datevec(now);
  disp(sprintf('*** Begin animal %d @ %02d:%02d ***',ana_animalnum,dateyyyymmdd(4:5)))
  
  for j = 1:length(anatomy_str) % for each anatomical plane
    % >>> read in the image
    useview = anatomy_str(j);
    viewi = strfind(animals.anatomy_str,useview);
    imgname = infocusimgs{i,viewi};
    pad_str = strfind(imgname,'_');
    Io = imread(strcat(comppaths.EBfldr,'\',imgname,'.tif'));
    [M,N,d] = size(Io);
    figpos = [680,186,size(Io,2)/4,size(Io,1)/4];
    
    % >>> Find image mask & erode to classify only clear foreground pixels
    maskname = sprintf('%s_fgndbkgnd',useview);
    seR = round(max(M,N)*0.01);
    se = strel('disk',seR);
    FBmask = getfield(animal_masks,maskname);
    FBmask = imerode(FBmask,se);
    
    analyzeind = find(FBmask == 1);
    bkgndind = find(FBmask == 0);
    
    Idv = reshape(double(Io(:)),M*N,d); % stores each channel (RGB) column-wise
    % Changes the RGB index (+1) triplet into scalar index
    Idind = sub2ind([256,256,256],Idv(analyzeind,1)+1,Idv(analyzeind,2)+1,Idv(analyzeind,3)+1);
    nfgndpxl = length(Idind);
    Ivecfgndpibarc = zeros(nfgndpxl,nclasses); % columns contain the p_ibarc value for the given pixel index
    Ivecfgndpcbari = zeros(nfgndpxl,nclasses); % columns contain the p_cbari value for the given pixel index
    [Uind,ia,ic] = unique(Idind); % only poll unique RGB value -- rather than visiting each pixel, visit each unique RGB value and populate relevant indices
    nU = length(Uind);
    indprobdef = ismember(Uind,non0p_ibarc); %ones(nU,1); % logical where 1 = probability is defined for those values of RGB and 0 zero probability for all classes
    Uindpibarc = zeros(nU,nclasses); %
    Uindpibarc(indprobdef==1,:) = p_ibarc(Uind(indprobdef==1),:);
    
    % For KDE type implementation, this should not be necessary since 
    rgbnondef = zeros(sum(indprobdef==0),3);
    if ~isempty(rgbnondef)
      [rgbnondef(:,1),rgbnondef(:,2),rgbnondef(:,3)] = ind2sub([256,256,256],Uind(indprobdef==0));
      [vi,d] =nearestNeighbor(DT,rgbnondef);
      DTclasspdfind = non0p_ibarc(vi);
      Uindpibarc(indprobdef == 0,:) = p_ibarc(non0p_ibarc(vi),:);
    end
    Uindpic = p_ic(Uind,:);
    Uindpcbari = Uindpic./repmat(p_iemp(Uind),1,nclasses);
    Ivecfgndpic = Uindpic(ic,:); % Gives p(i|c) column-wise for each class
    Ivecfgndpibarc = Uindpibarc(ic,:); % Gives p(i|c) column-wise for each class
    Ivecfgndpcbari = Uindpcbari(ic,:); % Gives p(c|i) (1 col/class) for the i @ that pixel position
    
    Ipic = zeros(M*N,nclasses);
    Ipic(analyzeind,:) = Ivecfgndpic;
    
    Ipibarc = zeros(M*N,nclasses);
    Ipibarc(analyzeind,:) = Ivecfgndpibarc;
    
    Ipcbari = zeros(M*N,nclasses);
    Ipcbari(analyzeind,:) = Ivecfgndpcbari;
    
    switch usepdf
      case 'p(i,c)'
        Idist = Ipic;
        pmaxstr = 'jnt';
      case 'p(c|i)'
        Idist = Ipcbari;
        pmaxstr = 'cbari';
      case 'p(i|c)'
        Idist = Ipibarc;
        pmaxstr = 'ibarc';
    end
    
    [iv,ind] = max(Idist,[],2);
    bldind = find(ind(analyzeind) == 1);
    EBind = find(ind(analyzeind) == 2);
    hemoind = find(ind(analyzeind) == 3);
    tsuind = find(ind(analyzeind) == 4);
    
    IclassM = zeros(M,N,nclasses);
    IclassM(analyzeind(bldind)) = 1;
    IclassM(analyzeind(EBind) + M*N) = 1;
    IclassM(analyzeind(hemoind) + 2*M*N) = 1;
    IclassM(analyzeind(tsuind) + 3*M*N) = 1;
    
    Ipcbari = reshape(Ipcbari,[M,N,nclasses]);
    Ipibarc = reshape(Ipibarc,[M,N,nclasses]);
    Ipic = reshape(Ipic,[M,N,nclasses]);
    switch usepdf
      case 'p(i,c)'
        Idist = Ipic;
      case 'p(c|i)'
        Idist = Ipcbari;
      case 'p(i|c)'
        Idist = Ipibarc;
    end
    
    IdxImask = Idist.*IclassM;
    
    close all;
    if plot_imgs
      % plots all classes with contour style
      %{
    close all;
    if plot_imgs
      figure('position',[100,500,400,420],...%figure('position',[104         497        1762         420]);
        'name',sprintf('%s_allclasses_%spc_EBfcnt4p3',...
        imgname(1:pad_str(2)-1),p_cstr));
      im2 = imagesc(Io); axis image
      hold on;
      for c = 1:length(class_str)
        contour(Ianmasks(:,:,c),[0,0],class_color{c})
      end
      set(gca,'position',[0,0,1,1])
      set(gcf,'position',[680,498,size(Io,2)/4,size(Io,1)/4])
      legend({'Blood', 'Evans Blue', 'Complex','Tissue'},...
        'location','northwest','color',[0.7,0.7,0.7])
    end
      %}
      
      % plots all classes with image overlay style
      %{.
      Iclass3 = nan(M*N,3);
      Iclass3(bkgndind,[1,2]) = 0;%1;        % Background yellow
      Iclass3(analyzeind(bldind),1) = 1;  % Blood red
      Iclass3(analyzeind(EBind),3) = 1;   % EB blue
      %FOR PURPLE HEMO CLASS
      %Iclass3(analyzeind(hemoind),3) = 0.75;
      %Iclass3(analyzeind(hemoind),1) = 0.75;
      %FOR GREEN HEMO CLASS
      Iclass3(analyzeind(hemoind),1) = 1; % complex green
      Iclass3(analyzeind(hemoind),2) = 2/3; % complex green
      Iclass3(analyzeind(tsuind),2) = 1; % complex green
      
      Iclass3 = reshape(Iclass3,[M,N,3]);
      
      figure('position',figpos,...%figure('position',[104         497        1762         420]);
        'name',['EBclass' imgname(pad_str(2)+1:pad_str(3)-1) 'raw']) %sprintf('%s_%sallclasses_%spc_EBfcnt1030', imgname(1:pad_str(2)-1),pmaxstr,pcstr));
      trans_factor = 0.15;%1;%
      im2 = imagesc(Io); set(im2,'cdatamapping','scaled'); axis image
      set(gca,'position',[0,0,1,1],'xticklabel','','yticklabel','');
      save_open_figures(comppaths.thesisfig,[],[],'','fmt','png');

      hold on;
      set(gcf,'name',['EBclass' imgname(pad_str(2)+1:pad_str(3)-1) 'class'])
      im1 = imagesc(Iclass3); set(im1,'Cdatamapping','scaled'); axis image;
      set(im1,'AlphaData',ones(size(im1,1),size(im1,2)).*trans_factor)
      hold off;
            %text(1,1,{[strrep(sprintf('%s, pdf: %s, p(c): %s',imgname(1:pad_str(2)-1),usepdf,pcstr),'_','\_')]},...
      %  'VerticalAlignment','top','color','w','fontsize',14)
      for c = 1:nclasses
        if ~isempty(class_color)
          text(N/10,c*(M/10),{[strrep(sprintf('%s',full_class_str{c}),'_','\_')]},...
            'VerticalAlignment','top','color',class_color(c,:),'fontsize',20,...
            'fontweight','bold')
        end
      end
      save_open_figures(comppaths.thesisfig,[],[],'','fmt','png');
      %}
      
%       for c = 1:nclasses
%         figure('position',[100,500,400,420],...%figure('position',[104         497        1762         420]);
%           'name',sprintf('%s_%s%s_%spc_EBfcnt1030',...
%           imgname(1:pad_str(2)-1), pmaxstr,class_str{c},pcstr));
%         imagesc(Idist(:,:,c)); axis image;
%         hold on;
%         contour(IclassM(:,:,c),[0,0],'w');
%         text(1,1,{[strrep(sprintf('%s',imgname(1:pad_str(2)-1)),'_','\_')],...
%           [sprintf('%s',full_class_str{c})],...
%           [sprintf('colormap: %s',usepdf)]},...
%           'VerticalAlignment','top','color','w','fontsize',14)
%         set(gca,'position',[0,0,1,1])
%         set(gcf,'position',figpos)
%       end
      if save_imgs; save_open_figures(strcat(comppaths.mtngpath,'\EBprocImgs'),[],[],'','fmt','png'); end
    end % end plot_imgs
    
    if saveMpath ~= 0; 
      save(strcat(saveMpath,...
        sprintf('%s_classmask',imgname(1:pad_str(2)-1))...
        ),'IclassM');
    end
    
    maskstr = sprintf('%s_pcbari', useview);
    if ~exist('anpdfI')
      anpdfI = struct(maskstr,Ipcbari);
    else
      anpdfI = setfield(anpdfI,maskstr,Ipcbari);
    end
    
    anpdfI = setfield(anpdfI,sprintf('%s_pibarc', useview),Ipibarc);
    anpdfI = setfield(anpdfI,sprintf('%s_pic', useview),Ipic);
    anpdfI = setfield(anpdfI,sprintf('%s_classes', useview),IclassM);
  end % anatomical directions for
  
  if ~exist('Ipdf')
    Ipdf = struct(sprintf('an%d',ana_animalnum),anpdfI);
  else
    Ipdf = setfield(Ipdf,sprintf('an%d',ana_animalnum),anpdfI);
  end
end % animals for
Ipdf.classes = class_str;


% >>> Construct the empirical p(c|i)
% p_cbari = p_ic ./ repmat(p_iemp,1,nclasses);
% p_cbari = p_cbari./repmat(sum(p_cbari,1),1,nclasses); %ROW WISE pdf

% >>> Structure for saving
pdfstruct.pkde_ibarc = pibarc_struct;
% pdfstruct.
pdfstruct.p_ibarc = p_ibarc;
pdfstruct.p_c = p_c;

end