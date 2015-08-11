%% Script: EB_brainimg_processing
% Description: 
% Let i = {r,g,b} and c = {bkgnd, blood, EB, complex, tissue}
% All input images are used to estimate p(i) and p(i,c = bkgnd). Only
% images with manually defined class masks will be used to build p(i|c = CLASS)
% while any tissue with no manually labeled pathology will be used to
% define p(i|c=tissue).
% INPUTS ----------------------------------------------------------------
%
%
% OUTPUTS ---------------------------------------------------------------
% • pdfhist: structure with fields related to class occurence matrix
% formation. The fields are as follows
%   >> all uint32, 256x256x256 matrices <<: 7 fields containing histogram
%   information for training data
%   FILED     DESCRIPTION                 REQUIRES LABEL?
%   - rgb:    all pixels in the image         N
%   - frgb:   all pixels in foreground        N (implicit via F/B seg)
%   - bkgnd:  all pixels in background        N         "
%   - bld:    red colored pathology           Y
%   - EB:     blue colored pathology          Y
%   - hemo:   black colored pathology         Y
%   - tsu:    all not pathological values
%   >> used4___: logicals indicating whether view (4, columns: D, L, R, V)
%   and animal (9, rows ordered according to animals.id) are used to
%   contribute to the class distribution
% • Imasks: structure with fields for each animal number (as indicated by
% animals.id)
%   >> an#: structue whose fields are M x N matrices containing the mask
%   for a particular anaomical direction (D/L/R/V) and class
%   (fgndbkgnd,bld,EB,hemo,tsu)
%    - ie: Imasks.an2.D_fgndbkgnd = M x N binary image corresponding to
%    final F/B seg mask
%    - If field is empty, it means no mask was available or produced
%
%  Date           Author              E-mail                      Version
%  29 Sept 2014   Amanda Balderrama   amanda.gaudreau@gmail.com     0
%     Draft version
%  29 Sept 2014   Amanda Balderrama   amanda.gaudreau@gmail.com     1
%     Conducts segmentation based on CV algorithm
%   5 Oct  2014   Amanda Balderrama   amanda.gaudreau@gmail.com     1_1
%     Uses CV algorithm and histogram manipulations to improve
%     segementation. Works well on animal 6, anatomical direction 4, but
%     seems to fail on more "simple" cases (directions 2,3)
%   6 Oct  2014   Amanda Balderrama   amanda.gaudreau@gmail.com     3
%     Segmentation via gradient magnitude edge filling
%   6 Oct  2014   Amanda Balderrama   amanda.gaudreau@gmail.com     4
%     Continued development of active contours and other off the shelf
%     matlab functions
%  14 Oct  2014   Amanda Balderrama   amanda.gaudreau@gmail.com     4.1
%     Uses 2-level (ds x 8, full image) edge contour segementation, then
%     uses files from L images with pathology to define class histograms
%  24 Oct  2014   Amanda Balderrama   amanda.gaudreau@gmail.com     4.2
%     1) Creates a small buffer around Fgnd/Bkgnd to form distributions -
%     does not use pixels near tissue/bkgnd boundary; 2) also build an rgb
%     distribution
%  30 Oct  2014   Amanda Balderrama   amanda.gaudreau@gmail.com     4.3
%     Create p(i) for just non-background pixels as well

loaddir
plotimgs = 1;
EBprocessingversion = 4.3;

%=========================================================================
% Initialize folder & File Variables
%=========================================================================
EBfldr = strcat(fpath,'EB_brainimg');
padno_and_logic = [...
  2,3,4,5,6,8,9,10,12;...1: animal idenifier
  1,1,1,0,0,0,0,0 ,1;...2: pad or no pad
  1,0,0,0,0,0,0,0 ,0;...3: blood visible on D surface
  1,0,1,1,1,0,1,1 ,0;...4: blood visible on L surface
  1,0,0,0,0,0,0,0 ,0;...5: blood visible on R surface
  1,0,1,0,1,0,0,0 ,0;...6: blood visible on V surface
  0,0,0,0,0,0,0,0 ,0;...7: EB visible on D surface
  0,1,1,1,1,1,1,1 ,0;...8: EB visible on L surface
  0,0,0,0,0,0,0,0 ,0;...9: EB visible on R surface
  0,1,1,1,1,1,1,1 ,0;...10: EB visible on V surface
  0,0,0,0,0,0,0,0 ,0;...11: Hemo visible on D surface
  1,0,1,1,1,1,1,1 ,0;...12: Hemo visible on L surface
  0,0,0,0,0,0,0,0 ,0;...13: Hemo visible on R surface
  1,0,1,1,1,1,1,1 ,0;...14: Hemo visible on V surface
  ]';
anatomy_str = 'DLRV';
class_str = {'bld','EB','hemo','tsu'};
full_class_str = {'Blood','EB','Hemorrhage','Tissue'};
class_color = {'r','b','k','g'};

animals.id = padno_and_logic(:,1);
animals.anatomy_str = anatomy_str;
animals.pad = padno_and_logic(:,2);
animals.bld = padno_and_logic(:,3:6);
animals.EB = padno_and_logic(:,7:10);
animals.hemo = padno_and_logic(:,11:14);
infocusimgs = {'20140819_pad2D_n08IPEB4','20140819_pad2L_n13IPEB4','20140819_pad2R_n16IPEB4','20140819_pad2V_n14IPEB4';...
  '20140819_pad3D_n08IPEB4','20140819_pad3L_n09IPEB4','20140819_pad3R_n12IPEB4','20140819_pad3V_n07IPEB4';...
  '20140819_pad4D_n06IPEB4','20140819_pad4L_n04IPEB4','20140819_pad4R_n07IPEB4','20140819_pad4V_n02IPEB4';...
  '20140819_nopad5D_n02IPEB4','20140819_nopad5L_n06IPEB4','20140819_nopad5R_n07IPEB4','20140819_nopad5V_n08IPEB4';...
  '20140819_nopad6D_n07IPEB4','20140819_nopad6L_n12IPEB4','20140819_nopad6R_n05IPEB4','20140819_nopad6V_n01IPEB4';...
  '20140821_nopad8D_n08IPEB4','20140821_nopad8L_n06IPEB4','20140821_nopad8R_n09IPEB4','20140821_nopad8V_n01IPEB4';...
  '20140821_nopad9D_n06IPEB4','20140821_nopad9L_n07IPEB4','20140821_nopad9R_n08IPEB4','20140821_nopad9V_n04IPEB4';...
  '20140821_nopad10D_n03IPEB4','20140821_nopad10L_n04IPEB4','20140821_nopad10R_n06IPEB4','20140821_nopad10V_n02IPEB4';...
  '20140821_pad12D_n03IPEB4','20140821_pad12L_n04IPEB4','20140821_pad12R_n07IPEB4','20140821_pad12V_n02IPEB4'};
animals.imgnames = infocusimgs;
animalnumvec = animals.id;%[animals.id(1)];%

pdfhist.rgb = uint32(zeros([256,256,256]));
pdfhist.frgb = uint32(zeros([256,256,256]));
pdfhist.bkgnd = uint32(zeros([256,256,256]));
pdfhist.bld = uint32(zeros([256,256,256]));
pdfhist.used4bld = zeros(length(animalnumvec),4);
pdfhist.EB = uint32(zeros([256,256,256]));
pdfhist.used4EB = zeros(length(animalnumvec),4);
pdfhist.hemo = uint32(zeros([256,256,256]));
pdfhist.used4hemo = zeros(length(animalnumvec),4);
pdfhist.tsu = uint32(zeros([256,256,256]));
pdfhist.used4tsu = zeros(length(animalnumvec),4);


for an = 1:length(animalnumvec) % for each animal number
  ana_animalnum = animalnumvec(an);
  i = find(animals.id == ana_animalnum);
  
  %   pdfhistANIMAL.rgb = uint32(zeros([256,256,256]));
  %   pdfhistANIMAL.bkgnd = uint32(zeros([256,256,256]));
  %   for c = 1:length(class_str)
  %     pdfhistANIMAL = setfield(pdfhistANIMAL,class_str{c},uint32(zeros([256,256,256])));
  %   end
  
  pdfhistANIMAL.rgb = uint32(zeros([256,256,256]));
  pdfhistANIMAL.frgb = uint32(zeros([256,256,256]));
  pdfhistANIMAL.bkgnd = uint32(zeros([256,256,256]));
  pdfhistANIMAL.bld = uint32(zeros([256,256,256]));
  pdfhistANIMAL.EB = uint32(zeros([256,256,256]));
  pdfhistANIMAL.hemo = uint32(zeros([256,256,256]));
  pdfhistANIMAL.tsu = uint32(zeros([256,256,256]));
  
  dateyyyymmdd = datevec(now);
  disp(sprintf('*** Begin animal %d @ %d:%02d ***',ana_animalnum,dateyyyymmdd(4:5)))
  %figure('name',sprintf('EBthreshmask_animal%d',ana_animalnum),'position',[680   263   578   715]);
  for j = 1:length(anatomy_str) % for each anatomical plane
    
    imgname = infocusimgs{i,j};
    pad_str = strfind(imgname,'_');
    Io = imread(strcat(EBfldr,'\',imgname,'.tif'));
    [M,N,d] = size(Io);
    f1 = sprintf('%s_fgndbkgnd',anatomy_str(j));
    
    F = ls(EBfldr);  % list all contents in the EB folder
    Fcell = cellstr(F);
    relFi = cellfun(@(x) isempty(x) == 0,strfind(Fcell,imgname));
    
    F = F(relFi,:);
    Fcell = cellstr(F);
    if ~exist('anstruct') || ~isfield(anstruct,f1) ||...
        isfield(anstruct,f1) && isempty(getfield(anstruct,f1))
      
      %% #2) Foreground/background segmentation using activecontour
      %{.
      %close all
      title_str = sprintf('%s\\%s\\%s',imgname(1:pad_str(1)-1),...
        imgname(pad_str(1):pad_str(2)-1), imgname(pad_str(2):end));
      
      % ---------------------------------------------------------------------
      % • Define input images
      % ---------------------------------------------------------------------
      Id = double(Io);
      ILeq = max(Id,[],3);%Iods;%mean(Id,3);%Id(:,:,3);%double(rgb2gray(Iods));%
      [M,N] = size(ILeq);
      [Gx,Gy] = imgradientxy(ILeq,'prewitt'); %figure; imshowpair(Gx,Gy,'montage')
      [Gm,Gd] = imgradient(Gx,Gy); %figure; imshowpair(Gm.*IoL_black,Gd.*IoL_black,'montage')
      
      % ---------------------------------------------------------------------
      % • Create preliminary mask from intensity threshold of image
      % ---------------------------------------------------------------------
      % >> fgnd/bkgnd level, gray level effectiveness metric
      [FBlvl glvlEM] = graythresh(Io);
      % >> initial binary mask from thresholding channels
      M0th = im2bw(Io,FBlvl);
      % >> code which looks at the "1" label and redefines the mask so that
      % there is only one connected component
      cc = bwconncomp(M0th);
      Lcc = cellfun(@(x) length(x),cc.PixelIdxList);
      [leni,keepi] = max(Lcc);
      M0th_keep = false(size(M0th));
      M0th_keep(cc.PixelIdxList{keepi}) = true;
      tinyseR = round(max(M,N)*0.005);
      se = strel('rectangle',[tinyseR,tinyseR]);
      M0th_keep = imdilate(M0th_keep,se);
      M0th_keep = imerode(M0th_keep,se);
      
      % >> Ensure that there is only 1 cc fgnd and one cc bkgnd --- this is
      % background validation
      M0th_inv = abs(M0th_keep-1);
      cc = bwconncomp(M0th_inv);
      Lcc = cellfun(@(x) length(x),cc.PixelIdxList);
      [leni,keepi] = max(Lcc);
      M0th_inv = true(size(M0th_inv));
      M0th_inv(cc.PixelIdxList{keepi}) = false;
      
      %subplot(2,2,j); imagesc(Io); axis image; hold on;
      %contour(M0th_keep,[0,0],'r');
      %ylabel(title_str)
      %set(gca,'xticklabel','','yticklabel','')
      
      % ---------------------------------------------------------------------
      % • dilates the cc intensity mask so that the initial contour is
      % guaranteed to be OUTSIDE the tissue boundary (requirment for AC
      % edge algorithm)
      % ---------------------------------------------------------------------
      seR = round(max(M,N)*0.06);
      se = strel('disk',seR);
      Minit = imdilate(M0th_inv,se);
      padc = round(N/100); 
      padr = round(M/100);
      Minit = padarray(Minit,[padr,padc],'replicate');
      Itest = padarray(Gm,[padr,padc]);
      
      %figure; imagesc(Io); hold on;
      %contour(Minit,[0,0],'w')
      
      % ---------------------------------------------------------------------
      % • Perform multi-level active contour segementation on gradient
      % magnitude image using threshold mask for initialization. Parameters
      % of down-sample factor (dsfactv) and iterations could be adjusted
      % for increased speed segementation accuracy
      % ---------------------------------------------------------------------
      dsfactv = [8,1];%[8,4,2,1];
      iter = [1500,1000];%[2000,1000,500,200]./2;
      nlevels = length(iter);
      smthfact = repmat(2,1,nlevels);
      for ell = 1:nlevels
        if ell == 1
          timev = zeros(nlevels,1);
          ACmasks = cell(nlevels,1);
        end
        Iac = imresize(Itest,1/dsfactv(ell));
        if ell == 1
          pnum = 1;
          M0th = imresize(Minit,1/dsfactv(ell),'nearest');
          seR = round(max(size(M0th))*0.05);
          se = strel('disk',seR);
        else
          M0th = imresize(ACmasks{ell-1},size(Iac),'nearest');
          M0th = imdilate(M0th,se);
        end
        
        tic;
        ACmasks{ell} = activecontour(Iac,M0th,iter(ell),'Edge',smthfact(ell));
        orig_rows = ceil((padr + (1:dsfactv(ell):M))./dsfactv(ell));
        orig_cols = ceil((padc + (1:dsfactv(ell):N))./dsfactv(ell));
        M0thnopad = M0th(orig_rows,orig_cols);
        ACnopad = ACmasks{ell}(orig_rows,orig_cols);
        timev(ell) = toc;
        
        if plotimgs;
          nsubfc = 1; nsubfr = 1;
          figure('position',[46,47,400*nsubfc,420*nsubfr],...
            'name',sprintf('%s_bkgndmaskdsf%d_EBfcnt4p2',...
            imgname(1:pad_str(2)-1),dsfactv(ell)));
          trans_factor = 0.5;
          im2 = imagesc(imresize(Io,1/dsfactv(ell))); axis image
          hold on;
          %im1 = imagesc(Iac); set(im1,'Cdatamapping','scaled'); axis image;
          set(gca,'xticklabel','','yticklabel','');
          %set(im1,'AlphaData',ones(size(im1,1),size(im1,2)).*trans_factor)
          contour(M0thnopad,[0,0],'w')
          contour(ACnopad,[0,0],'edgecolor','r');%,'linewidth',2);
          set(gca,'position',[0,0,1,1])
          set(gcf,'position',[680,498,size(Io,2)/4,size(Io,1)/4])
          text(1,1,{sprintf('DSfactor = %d',dsfactv(ell)),...
            sprintf('# iter = %d',iter(ell)),...
            sprintf('time = %dm%ds',floor(timev(ell)/60),round(rem(timev(ell),60)))},...
            'VerticalAlignment','top','color','w');
        end
      end
      save_open_figures(EBfldr,[],[],'','fmt','png');
          close all;
     
      %set(gca,'position',[0,0,1,1])
      %set(gcf,'position',[680,560,size(Io,2)/8,size(Io,1)/8])
      
      fbmask = ACnopad;
      if ~exist('anstruct')
        anstruct = struct(f1,logical(fbmask));
      else
        anstruct = setfield(anstruct,f1,logical(fbmask));
      end
    end % segement of code is skipped if anstruct.X_fgndbkgnd already exists
    
    seR = round(max(M,N)*0.005);
    se = strel('disk',seR);
    % ---------------------------------------------------------------------
    % • Incrimentally build histogram of rgb values
    % ---------------------------------------------------------------------
    new_hist = build_rgbpdf(reshape(Id(:),M*N,3));
    pdfhistANIMAL.rgb = pdfhistANIMAL.rgb + new_hist;
    pdfhist.rgb = pdfhist.rgb + new_hist;
    
    % ---------------------------------------------------------------------
    % • Incrimentally build histogram of rgb values using only foreground
    % pixels
    % ---------------------------------------------------------------------
    fbmask_e = imerode(fbmask,se);
    pxls = find(repmat(fbmask_e,[1,1,3]) == 1);
    rgb_fgnd = reshape(Id(pxls),length(pxls)/3,3);
    new_fhist = build_rgbpdf(rgb_fgnd);
    pdfhistANIMAL.frgb = pdfhistANIMAL.frgb + new_fhist;
    pdfhist.frgb = pdfhist.frgb + new_fhist;
    
    % ---------------------------------------------------------------------
    % • Incrimentally build histogram of rgb values belonging to background
    % (mask == 0)
    % ---------------------------------------------------------------------
    fbmask_d = imdilate(fbmask,se);
    bkgndpxls = find(repmat(fbmask_d,[1,1,3]) == 0);
    rgb_bkgnd = reshape(Id(bkgndpxls),length(bkgndpxls)/3,3);
    new_hist = build_rgbpdf(rgb_bkgnd);
    pdfhistANIMAL.bkgnd = pdfhistANIMAL.bkgnd + new_hist;
    pdfhist.bkgnd = pdfhist.bkgnd + new_hist;
    disp(sprintf('=> DONE: Fgnd/Bkgnd seg for animal %d%s (total time %1.2f min)',ana_animalnum,anatomy_str(j),sum(timev)/60))
    
    %% #2.3) Loads & processes Photoshop mask information
    %{.
    %close all; figure; image(Io); axis image
    clogicvec = zeros(1,3);
    for c = 1:length(class_str)
      fc = sprintf('%s_%s',anatomy_str(j),class_str{c});
      if isempty(strmatch('tsu',class_str{c}))
        classlog = getfield(animals,class_str{c});
        clogicvec(c) = classlog(i,j);
        relFi = cellfun(@(x) isempty(x) == 0,strfind(Fcell,strcat('mask',class_str{c})));
        if clogicvec(c) == 0 && sum(relFi) == 0
          disp(sprintf('Class %s NOT EXPECTED for image %s, no file found',full_class_str{c},imgname));
        elseif clogicvec(c) == 1 && sum(relFi) == 0
          disp(sprintf('Missing File: Class %s IS EXPECTED for image %s',full_class_str{c},imgname));
        elseif clogicvec(c) == 0 && sum(relFi) == 1
          disp(sprintf('Incorrectly definited class: Class %s NOT EXPECTED for image %s, but file exists',full_class_str{c},imgname));
        end
        % If any files for the anatomical direction and and class being
        % considered are found, the masking file is read in and added to
        % the animal structure anstruct
        if any(relFi) && clogicvec(c) == 1% relevant mask file found
          if ~isfield(anstruct,fc)
            classan10 = getfield(pdfhist,sprintf('used4%s',class_str{c}));
            classan10(an,j) = 1;
            pdfhist = setfield(pdfhist,sprintf('used4%s',class_str{c}),classan10);
            
            classimgname = strcat(EBfldr,'\',Fcell{relFi});
            Ic = imread(classimgname); % read in mask file
            bld_ind = find(Ic ~= max(Ic(:))); % identify pixels that are not WHITE
            Icmask3 = zeros(size(Ic)); % M x N x 3
            Icmask3(bld_ind) = 1;
            anstruct = setfield(anstruct,fc,logical(Icmask3(:,:,1)));
            
            % ---------------------------------------------------------------------
            % • Incrimentally build histogram of rgb values which describe
            % class_str{c}
            % ---------------------------------------------------------------------
            pxls = find(Icmask3 == 1);
            rgb_pxls = reshape(Id(pxls),length(pxls)/3,3);
            new_hist = build_rgbpdf(rgb_pxls);
            pdfhistANIMAL = setfield(pdfhistANIMAL,class_str{c},getfield(pdfhistANIMAL,class_str{c}) + new_hist);
            pdfhist = setfield(pdfhist,class_str{c},getfield(pdfhist,class_str{c}) + new_hist);
          end
        else % no mask image file exists for the class being analyzed
          anstruct = setfield(anstruct,fc,[]);
        end
      else % tissue class definition
        if sum(clogicvec) == 0 % anatomical direction being considered no pathology so F/B mask is fully part of the tissue class
          anstruct = setfield(anstruct,fc,logical(fbmask));
          classan10 = getfield(pdfhist,sprintf('used4%s',class_str{c}));
          classan10(an,j) = 1;
          pdfhist = setfield(pdfhist,sprintf('used4%s',class_str{c}),classan10);
          % ---------------------------------------------------------------------
          % • Incrimentally build histogram of rgb values which describe
          % tissue
          % ---------------------------------------------------------------------
          pdfhistANIMAL = setfield(pdfhistANIMAL,class_str{c},...
            getfield(pdfhistANIMAL,class_str{c}) + new_fhist);
          pdfhist = setfield(pdfhist,class_str{c},...
            getfield(pdfhist,class_str{c}) + new_fhist);
        else
          anstruct = setfield(anstruct,fc,[]);
        end
      end % pathology vs tissue if
    end % image classes for
    disp(sprintf('----- FINISHED animal %d%s (total time %1.2f min)',ana_animalnum,anatomy_str(j),toc/60))
    save(strcat(EBfldr,'\',sprintf('%s_EBclasshistmed',dateyyyymmddstr)))
  end % anatomical directions for
  anstruct.file_header = sprintf('%s',imgname(1:pad_str(2)-2));
  save(strcat(EBfldr,'\',sprintf('%s_EBclasshistmed_animal%d',dateyyyymmddstr,ana_animalnum)),'pdfhistANIMAL')
  save(strcat(EBfldr,'\',sprintf('%s_EBclassmasksmed_animal%d',dateyyyymmddstr,ana_animalnum)),'anstruct')
  if ~exist('Imasks')
    Imasks = struct(sprintf('an%d',ana_animalnum),anstruct);
  else
    Imasks = setfield(Imasks,sprintf('an%d',ana_animalnum),anstruct);
  end
  save(strcat(EBfldr,'\',sprintf('%s_EBclasshistmed',dateyyyymmddstr)))
  clear anstruct
end % animals for