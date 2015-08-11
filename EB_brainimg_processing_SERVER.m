%% Script: EB_brainimg_processing
% Description:
% INPUTS ----------------------------------------------------------------
%
%
% OUTPUTS ---------------------------------------------------------------
%
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

loaddir
plotimgs = 1;

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
animalnumvec = 6;%animals.id;

for an = 1:length(animalnumvec) % for each animal number
  ana_animalnum = animalnumvec(an);
  i = find(animals.id == ana_animalnum);
  for j = 1:length(anatomy_str) % for each anatomical plane
    
    imgname = infocusimgs{i,j};
    pad_str = strfind(imgname,'_');
    Io = imread(strcat(EBfldr,'\',imgname,'.tif'));
    Iraw = imread(strcat(EBfldr,'\',imgname,'.tif'));
    [M,N,d] = size(Io);
    dsfact = 1;
    Io = imresize(Io,1/dsfact);
    f1 = sprintf('%s_fgndbkgnd',anatomy_str(j));
    %anstruct = struct([]);
    
    F = ls(EBfldr);
    Fcell = cellstr(F);
    relFi = cellfun(@(x) isempty(x) == 0,strfind(Fcell,imgname));
    
    F = F(relFi,:);
    Fcell = cellstr(F);
    if ~exist('anstruct') || ~isfield(anstruct,f1) ||...
        isfield(anstruct,f1) && isempty(getfield(anstruct,f1))
      
      %% #2) Foreground/background segmentation using activecontour
      %{.
      close all
      title_str = sprintf('%s\\%s\\%s',imgname(1:pad_str(1)-1),...
        imgname(pad_str(1):pad_str(2)-1), imgname(pad_str(2):end));
      
      Id = double(Io);
      ILeq = max(Id,[],3);%Iods;%mean(Id,3);%Id(:,:,3);%double(rgb2gray(Iods));%
      [M,N] = size(ILeq);
      [Gx,Gy] = imgradientxy(ILeq,'prewitt'); %figure; imshowpair(Gx,Gy,'montage')
      [Gm,Gd] = imgradient(Gx,Gy); %figure; imshowpair(Gm.*IoL_black,Gd.*IoL_black,'montage')
      % Mask: dilated intensity mask of the image
      [FBlvl glvlEM] = graythresh(Iods); % fgnd/bkgnd level, gray level effectiveness metric
      M0th = im2bw(Iods,FBlvl); % M0th = im2bw(ILeq./max(ILeq(:)),FBlvl);
      seR = round(max(M,N)*0.1);
      se = strel('disk',seR);
      M0th = imdilate(M0th,se);
      Itest = Gm;

      dsfactv = [8,1];%[8,4,2,1];
      iter = [2000,750];%[2000,1000,500,200]./2;
      nlevels = length(iter);
      smthfact = repmat(2,1,nlevels);
      nsm = length(smthfact);
      for i = 1:nlevels
        if i == 1
          timev = zeros(nlevels,1);
          ACmasks = cell(nlevels,1);
        end
        Iac = imresize(Gm,1/dsfactv(i));
        if i == 1
          pnum = 1;
          M0th = imresize(M0th,1/dsfactv(i),'nearest');
          seR = round(max(size(M0th))*0.05);
          se = strel('disk',seR);
        else
          M0th = imresize(ACmasks{i-1},size(Iac),'nearest');
          M0th = imdilate(M0th,se);
        end
        
        tic;
        ACmasks{i} = activecontour(Iac,M0th,iter(i),acalgostr{a},smthfact(i));
        timev(i) = toc;
      end
      
      Iot = zeros(size(Io));
      Iomask3 = repmat(ACmasks{i},[1,1,3]);
      Iot = Iomask3.*Id;
      
      if plotimgs; figure; image(Iomask3); axis image; end
      if plotimgs; figure; image(uint8(Iot)); axis image; end
      if plotimgs; figure; image(uint8(Io)); axis image; end
      %set(gca,'position',[0,0,1,1])
      %set(gcf,'position',[680,560,size(Io,2)/8,size(Io,1)/8])
      fv1 = Iomask3(:,:,1);
      if ~exist('anstruct')
        anstruct = struct(f1,logical(fv1));
      else
        anstruct = setfield(anstruct,f1,logical(fv1));
      end
    end
    %% #2.3) Loads & processes Photoshop mask information
    %{.
    close all
    % figure; image(Io); axis image
    clogicvec = zeros(1,3);
    for c = 1:length(class_str)
      fc = sprintf('%s_%s',anatomy_str(j),class_str{c});
      if isempty(strmatch('tsu',class_str{c}))
        classlog = getfield(animals,class_str{c});
        clogicvec(c) = classlog(i,j);
        relFi = cellfun(@(x) isempty(x) == 0,strfind(Fcell,strcat('mask',class_str{c})));
        if clogicvec(c) ~= sum(relFi); disp(sprintf('File not found for %s class, but no class information expected from image %s',full_class_str{c},imgname)); end
        % If any files for the anatomical direction and and class being
        % considered are found, the masking file is read in and added to
        % the animal structure anstruct
        if any(relFi) % relevant mask file found
          if ~isfield(anstruct,fc)
            classimgname = strcat(EBfldr,'\',Fcell{relFi});
            Ic = imread(classimgname); % read in mask file
            bld_ind = find(Ic ~= max(Ic(:))); % identify pixels that are not WHITE
            Icmask3 = zeros(size(Ic)); % M x N x 3
            Icmask3(bld_ind) = 1;
            anstruct = setfield(anstruct,fc,logical(Icmask3(:,:,1)));
            
            if plotimgs
              Ib_v = reshape(Ic(bld_ind),[length(bld_ind)/3,3]);
              figure('position',[403,251,1067,415]); subplot(1,2,2);
              image(Ic); hold on; axis image;
              title(sprintf('%s %s',title_str,full_class_str{c}));
              
              subplot(1,2,1); plot3(Ib_v(:,1),Ib_v(:,2),Ib_v(:,3),sprintf('%s.',class_color{c}));
              grid on;
              xlabel('R'); ylabel('G'); zlabel('B');
              xlim([0,255]);ylim([0,255]);zlim([0,255]);
              title(sprintf('RGB data of %s',full_class_str{c}));
              hold on;
              %   plot3(IEB_v(:,1),IEB_v(:,2),IEB_v(:,3),'b.'); hold on;
              %   plot3(Itsu_v(:,1),Itsu_v(:,2),Itsu_v(:,3),'g.');
              %   plot3(Ibld_v(:,1),Ibld_v(:,2),Ibld_v(:,3),'r.');
              %   plot3(Ihmg_v(:,1),Ihmg_v(:,2),Ihmg_v(:,3),'k.');
            end
          end
        else % no mask image file exists for the class being analyzed
          anstruct = setfield(anstruct,fc,[]);
        end
      else % tissue class definition
        if sum(clogicvec) == 0 % anatomical direction being considered no pathology so F/B mask is fully part of the tissue class
          anstruct = setfield(anstruct,fc,logical(fv1));
        else
          anstruct = setfield(anstruct,fc,[]);
        end
      end % pathology vs tissue if
    end % image classes for
  end % anatomical directions for
  anstruct.file_header = sprintf('%s',imgname(1:pad_str(2)-1));
  if ~exist('Imasks')
    Imasks = struct(sprintf('an%d',ana_animalnum),anstruct);
  else
    Imasks = setfield(Imasks,sprintf('an%d',ana_animalnum),anstruct);
  end
end % animals for