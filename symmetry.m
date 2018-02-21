%% This Algorithm determines whether a skin lesion has 'similarity' lines.
% Copyrighted by Tyler L. Coye (2015)
% Publication or commercial use of any part of this script is not
% authorized unless approved by the Author.
%
% I say that this script determines similarity lines because it searches
% for lines that split the lesion into two sides which are similar and not
% necessarily completely symmetric. I believe that for practical purposes,
% lines that yeild >90% similarity can be considered symmetry lines in
% medical images of lesions.
%
% This script first segments a skin lesion using an algorithm I
% developed previously. Then, this segmented image is then rotated by 1
% degree for 360 degrees. At each degree rotation, a Jaccard Index is
% calculated and the Jaccard Distance is calculated. If the Jaccard distace
% is less than .072 ( which is equivalent to 92.8% similarity) 1 is added to a
% count. I had made the assumption beforehand that each side of the mole does
% not need to be identical (jaccard distance = 0). This script gives you a
% basic idea of how one can count lines of 'similarity' in a skin lesion.
%
% The output of this script is the number of similatiry lines as well as their
% loation
%
% Note: You can change the Jaccard distance to whichever value you like, for
% looser or tighter similarities. Just alter its value in the code below.
% clear all
%{
%% ***PART I-SEGMENTATION***
% Read image
I = imread('test3.jpg');
im = im2double(I);
% im = Jmat(:,:,1);
% Convert RGB to Gray via PCA
lab = rgb2lab(im);
f = 0;
wlab = reshape(bsxfun(@times,cat(3,1-f,f/2,f/2),lab),[],3);
[C,S] = pca(wlab);
S = reshape(S,size(lab));
S = S(:,:,1);
imgray = (S-min(S(:)))./(max(S(:))-min(S(:)));
% Morphological Closing
se = strel('disk',1);
imgclosed = imclose(imgray,se);
% Complement Image
K= imcomplement(imgclosed);
% 2-D wavelet Decomposition using B-Spline
[cA,cH,cV,cD] = dwt2(K,'bior1.1');
%% Otsu thresholding on each of the 4 wavelet outputs
thresh1 = multithresh(cA);
thresh2 = multithresh(cH);
thresh3 = multithresh(cV);
thresh4 = multithresh(cD);
% Calculating new threshold from sum of the 4 otsu thresholds and dividing by 2
level = (thresh1 + thresh2 + thresh3 + thresh4)/2;
% single level inverse discrete 2-D wavelet transform
X = idwt2(cA,cH,cV,cD,'bior1.1');
% Black and White segmentation
BW=imquantize(X,level);
% Iterative Canny Edge (Novel Method)
BW1 = edge(edge(BW,'canny'), 'canny');
% Post-Processing
BW3 = imclearborder(BW1);
CC = bwconncomp(BW3);
S = regionprops(CC, 'Area');
L = labelmatrix(CC);
BW4 = ismember(L, find([S.Area] >= 100));
%}
BW4 = tMask;%dp.mveparams.tsuMask;
origimgsize = size(BW4);
sqimsize = origimgsize;%min(origimgsize).*[1,1];%
BW4 = imresize(BW4,sqimsize);
BW51 = imfill(BW4,'holes');
BW5 = imcomplement(BW51);
imgbw = bwlabel(BW51);
ss_mask = bwlabel(imgbw);
stats = regionprops(ss_mask, 'BoundingBox', 'Area','centroid');
Astats = [stats.Area];
idx = find(Astats == max(Astats));
ss_mask = ismember(ss_mask, idx);
out_img = imcrop(BW51, stats(idx).BoundingBox);
Imask = imresize( dp.(usedat.type).Cu63,sqimsize);
Imask = imcrop(Imask, stats(idx).BoundingBox);

%% Search for Number Symmety Lines using Jaccard Similarity Coefficient
count = 0;
close all
degspace = 2;
evaldeg = -60:degspace:60;%-10:-0.5:-15;%
numerator = zeros(1,length(evaldeg));
denomenator = zeros(1,length(evaldeg));
jaccardDistance = zeros(1,length(evaldeg));
overlapU = zeros(1,length(evaldeg));
overlapI = zeros(1,length(evaldeg));
i = 0;
for k=evaldeg(:)'
  i = i + 1;
  B = imrotate(out_img,k);
  C = fliplr(B);
  intBC = B & C;
  uBC = B | C;
  overlapU(i) = sum(xor(uBC(:),B(:)));
  overlapI(i) = sum(xor(intBC(:),B(:)));
  numerator(i) = sum(intBC(:));
  denomenator(i) = sum(uBC(:));
  jaccardIndex = numerator(i)/denomenator(i);
  jaccardDistance(i) = 1 - jaccardIndex;
  if jaccardDistance(i) < 0.19%.072 % you can change this value ( The closer to zero, the more similar)
    count = count + 1;
    if plotimgs
      figure('name',num2str(k));
      subplot(3,2,1); imagesc(B);
      subplot(3,2,2); imagesc(C);
      subplot(3,2,3); imagesc(intBC);
      subplot(3,2,4); imagesc(uBC);
      colormap(gray);
      subplot(3,2,5);
      D = imresize(B, 2);
      [m n] = size(D);
      imagesc(D); colormap(gray); hold on;
      line([n/2,n/2],[0,m],'Color','r','LineWidth',2)
      text(0,-12,strcat('\color{blue}\fontsize{16}Degrees:',num2str(k)))
      hold off
    end
  end
end
final_count = count/2; % Generally you will get 2 k's that correspond to the same line, so you divide by 2.
[ppos,pind] = findpeaks(numerator./denomenator,'MINPEAKDISTANCE',10/degspace);
if plotimgs
  figure; plot(evaldeg, [numerator./sum(B(:)); denomenator./sum(B(:)); numerator./denomenator],'.');%jaccardDistance,'.'); %
  hold on; plot(evaldeg(pind),ppos,'xk');
  legend({'\Sigma B \cap C','\Sigma B \cup C','1-dist_{jac}'})
  xlabel('angle (\circ)'); ylabel('\Sigma B \cap C/\Sigma B \cup C');
end
% angDeg = -12.5;
% rot0 = angDeg;
% tx0 = 0;%maskcc{Mi}.Centroid(1) - maskcc{Fi}.Centroid(1);
% ty0 = 0;%maskcc{Mi}.Centroid(2) - maskcc{Fi}.Centroid(2);
% sx0 = 1;%maskcc{Mi}.MinorAxisLength/maskcc{Fi}.MinorAxisLength;
% sy0 = 1;% maskcc{Mi}.MajorAxisLength/maskcc{Fi}.MajorAxisLength;
% mu0 = [tx0,ty0,rot0*pi/180,sx0,sy0,0];
% [TxM,Amask] = transform_image(double(out_img),mu0,'extrapval',0);
% figure; imagesc(TxM); colormap(gray); colorbar;
%


%% Procedure to select the most relevant peak by testing which has the minimum
%   difference between L and R side
%[jval,jind ] = sort(jaccardDistance);
jind = pind;
npks = length(pind);
nel = length(dp.mveparams.elnames);
elL = zeros(nel,npks);
elR = zeros(nel,npks);
nL = zeros(npks,1);
nR = zeros(npks,1);
elmean = zeros(nel,1);
for j = 1:npks%find(find(evaldeg == 21) == jind);%1; %
  k = evaldeg(jind(j));
  B = imrotate(out_img,k);
  C = fliplr(B);
  intBC = B & C;
  uBC = B | C;
  [m n] = size(B);
  if plotimgs
    figure('name',num2str(k),'position',[680   285   547   693]);
    subplot(3,2,1); imagesc(B); title('B: original rotated image');
    subplot(3,2,2); imagesc(C); title('C: mirrored rotated image');
    subplot(3,2,3); imagesc(intBC);  title('B \cap C');
    subplot(3,2,4); imagesc(uBC); title('B \cup C');
    colormap(gray);
    subplot(3,2,5);
    D = B;
    imagesc(D); colormap(gray); hold on;
    line([n/2,n/2],[0,m],'Color','r','LineWidth',2)
    text(0,-12,strcat('\color{blue}\fontsize{16}Degrees:',num2str(k)))
    contour(out_img,[0,0],':r');
    hold off
    Imgrot = imrotate(Imask,k);
    subplot(3,2,6); imagesc(Imgrot,usedat.datarng); colormap(gray); hold on;
    line([n/2,n/2],[0,m],'Color','r','LineWidth',2)
  end
  
  nL(j) = sum(sum(B(:,1:floor(n/2))));
  nR(j) = sum(sum(B(:,ceil(n/2):end)));
  for el = 1:nel
    Iel = imresize( dp.(usedat.type).(dp.mveparams.elnames{el}),sqimsize);
    Iel = imrotate(imcrop(Iel, stats(idx).BoundingBox),k);
    Iel(~B) = 0;
    elL(el,j) = sum(sum(Iel(:,1:floor(n/2))));
    elR(el,j) = sum(sum(Iel(:,ceil(n/2):end)));
    if j == 1
      elmean(el) = nanmean(Iel(B(:)));
    end
  end
end

% 2015-09-09 Note: Method for selecting optimal peak to search around has
%     many potential weaknesses and failure modes. More robust methods for
%     selection should be explored
dLR = abs(elL./repmat(nL(:)',nel,1) - elR./repmat(nR(:)',nel,1));%./repmat(elmean.*sum(out_img(:)),1,npks);
[mdLRval, dLRi] = min(sum(dLR));

%% Searches a smaller region around the best matching angle
ang_opt = evaldeg(pind(dLRi));
degspace2 = 0.5;
evaldeg2 = (ang_opt - 5):degspace2:(ang_opt + 5);
numerator2 = zeros(1,length(evaldeg2));
denomenator2 = zeros(1,length(evaldeg2));
jaccardDistance2 = zeros(1,length(evaldeg2));
elL = zeros(nel,length(evaldeg2));
elR = zeros(nel,length(evaldeg2));
nL = zeros(length(evaldeg2),1);
nR = zeros(length(evaldeg2),1);
i = 0;
for k=evaldeg2(:)'
  i = i + 1;
  B = imrotate(out_img,k);
  C = fliplr(B);
  intBC = B & C;
  uBC = B | C;
  numerator2(i) = sum(intBC(:));
  denomenator2(i) = sum(uBC(:));
  jaccardIndex = numerator2(i)/denomenator2(i);
  jaccardDistance2(i) = 1 - jaccardIndex;
  if jaccardDistance2(i) < 0.19%.072 % you can change this value ( The closer to zero, the more similar)
    count = count + 1;
    if plotimgs
      figure('name',num2str(k));
      subplot(3,2,1); imagesc(B);
      subplot(3,2,2); imagesc(C);
      subplot(3,2,3); imagesc(intBC);
      subplot(3,2,4); imagesc(uBC);
      colormap(gray);
      subplot(3,2,5);
      D = imresize(B, 2);
      [m n] = size(D);
      imagesc(D); colormap(gray); hold on;
      line([n/2,n/2],[0,m],'Color','r','LineWidth',2)
      text(0,-12,strcat('\color{blue}\fontsize{16}Degrees:',num2str(k)))
      hold off
    end
  end
  
  nL(i) = sum(sum(B(:,1:floor(n/2))));
  nR(i) = sum(sum(B(:,ceil(n/2):end)));
  for el = 1:nel
    Iel = imresize( dp.(usedat.type).(dp.mveparams.elnames{el}),sqimsize);
    Iel = imrotate(imcrop(Iel, stats(idx).BoundingBox),k);
    Iel(~out_img) = 0;
    elL(el,i) = sum(sum(Iel(:,1:floor(n/2))));
    elR(el,i) = sum(sum(Iel(:,ceil(n/2):end)));
  end
end
[ppos2,pind2] = findpeaks(numerator2./denomenator2,'MINPEAKDISTANCE',5/degspace2);
if plotimgs
  figure; plot(evaldeg2, numerator2./denomenator2,'.');%[numerator2./sum(B(:)); denomenator2./sum(B(:)); numerator2./denomenator2],'.');%jaccardDistance,'.'); %
  hold on; plot(evaldeg2(pind2),ppos2,'xk');
  legend({'1-dist_{jac}'})
  xlabel('angle (\circ)'); ylabel('\Sigma B \cap C/\Sigma B \cup C');
end

k = (evaldeg2(pind2(1)));%-10;%19;%evaldeg2(dLRi);%
if plotimgs%1% 
  B = imrotate(out_img,k);
  C = fliplr(B);
  intBC = B & C;
  uBC = B | C;
  figure('name',num2str(k),'position',[680   285   547   693]);
  subplot(3,2,1); imagesc(B); title('B: original rotated image');
  subplot(3,2,2); imagesc(C); title('C: mirrored rotated image');
  subplot(3,2,3); imagesc(intBC);  title('B \cap C');
  subplot(3,2,4); imagesc(uBC); title('B \cup C');
  colormap(gray);
  subplot(3,2,5);
  D = B;
  [m n] = size(D);
  imagesc(D); colormap(gray); hold on;
  line([n/2,n/2],[0,m],'Color','r','LineWidth',2)
  text(0,-12,strcat('\color{blue}\fontsize{16}Degrees:',num2str(k)))
  contour(out_img,[0,0],':r');
  hold off
  Imgrot = imrotate(Imask,k);
  subplot(3,2,6); imagesc(Imgrot,usedat.datarng); colormap(gray); hold on;
  line([n/2,n/2],[0,m],'Color','r','LineWidth',2)
  disp(sprintf('Final angle: %1.2f',k))
end

% figure; plot(elL./elR)
% hold on; plot(elR,':');
% dLR = abs(elL - elR);
% figure; plot(dLR);
% [mdLRval, dLRi] = min(sum(dLR));
% evaldeg2(dLRi)

