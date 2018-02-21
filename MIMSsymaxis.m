% Copyright © 2014 New York University.
% See notice at the end of this file.

function axOut = MIMSsymaxis(I)
plotimgs = 0;
%=========================================================================
%_________________________________________________________________________

%% compute wavelet coefficients
% --------------------------------------------------
nangs = 32;
stretch = 1;
scale = 2;
hopsize = 5;
halfwindowsize = 1;
magthreshold = 0.01;
ignoredirections = 1;
[m,a,x,y] = coefficientslist(ignoredirections,I,nangs,stretch,scale,hopsize,halfwindowsize,magthreshold);

%% triangles
% --------------------------------------------------
disp('Listing triangles...')
% tic
ts = triangleslist(m,a,x,y,@force_sym,size(I,2),size(I,1));
% toc

%% projection
% --------------------------------------------------
m = 0.5+0.5*m;
p4 = [1 0 0 0]; % wmp, mp, wpp, pp
accwidth = 2*ceil(sqrt(size(I,1)^2+size(I,2)^2))+1; % displacement
accheight = 360; % angle
disp('Projecting...')
% tic
A = project2(ts,@location_sym,@force_sym,accwidth,accheight,size(I,2),size(I,1),p4);
% toc

%% finding local maxima in accumulator
% --------------------------------------------------
hsize = 10;
halfwindow = 10;
mindistbetcent = 10;
lowerbound = 0.8;
minarea = 0.1;
nlocmax = 10; % just looking for the main symmetry axis
[locs,Out1,Out2] = localmaxima(A,hsize,halfwindow,lowerbound,minarea,mindistbetcent,nlocmax);
% return

axOut.gamma = (locs(1,:)-1)/accheight*pi;
axOut.gammaDEG = axOut.gamma.*180/pi;
axOut.displacement = locs(2,:) - (accwidth-1)/2;
axOut.v = [cos(axOut.gamma);sin(axOut.gamma)];
axOut.closestpnt = repmat(axOut.displacement,2,1).*axOut.v;
axOut.locs = locs;
axOut.vp = nan(size(axOut.v));
axOut.AccMean = nan(1,size(locs,2));
axOut.LineSlope = nan(1,size(locs,2));
axOut.LineInt = nan(1,size(locs,2));
for i = 1:size(locs,2)
  axOut.AccMean(i) = mean(mean(Out1(max([1,locs(1,i)-5]):min([accheight,locs(1,i)+5]),...
    max([1,locs(2,i)-5]):min([accwidth,locs(2,i)+5]))));
  if axOut.gamma(i) <= pi/2
    axOut.vp(:,i) = [-axOut.v(2,i);axOut.v(1,i)];
  else
    axOut.vp(:,i) = [axOut.v(2,i);-axOut.v(1,i)];
  end
  y1 = 0;
  k = ((y1-axOut.closestpnt(1,i))/(axOut.vp(1,i)));
  x1 = axOut.closestpnt(2,i) + k*axOut.vp(2,i);
  x2 = 0;
  k = ((x2-axOut.closestpnt(2,i))/(axOut.vp(2,i)));
  y2 = axOut.closestpnt(1,i) + k*axOut.vp(1,i);
  axOut.LineSlope(i) = (y1 - y2)/(x1 - x2);
  axOut.LineInt(i) = (y2*x1 - y1*x2)/(x1-x2);
end
axOut.imgsize = size(I);


% --------------------------------------------------
%% draw results
% --------------------------------------------------
% [G,J,lineparams] = paintlines(Out1,Image,accwidth,accheight,locs(:,:));
% G = imresize(G,[size(J,1) size(J,2)]);
if plotimgs
  lncol = custom_colormaps('lines',30);
  figure; subplot(2,3,1); imagesc(A); colorbar; title('Raw Accumulator (A), project2');
  subplot(2,3,2); imagesc(Out1); colormap(gray); colorbar; title('Smoothed A (Out1), localmaxima');
  subplot(2,3,3); imagesc(Out2); hold on; plot(locs(2,:),locs(1,:),'.r'); colormap(gray); colorbar;
  title('Thresh. Accumulator (Out2), localmaxima');
  subplot(2,3,4); imagesc(Out1); hold on;%imagesc(G);%figure; imshow([G J])
  for i = 1:size(axOut.locs,2)
    plot(axOut.locs(2,i),axOut.locs(1,i),'s','markerfacecolor',lncol(i,:),'markeredgecolor',lncol(i,:),'markersize',5);
  end
  title('G, paintlines')
  subplot(2,3,5); imshow(I); hold on;
  axJ = gca;
  %   imagesc(Image); colormap(gray); hold on;
  for i = 1:length(axOut.displacement)
    y1 = 1;
    k = (y1-lineparams.closestpnt(1,i))/(lineparams.vp(1,i));
    x1 = lineparams.closestpnt(2,i) + k*lineparams.vp(2,i);
    y2 = isz(1);
    k = (y2-lineparams.closestpnt(1,i))/(lineparams.vp(1,i));
    x2 = lineparams.closestpnt(2,i) + k*lineparams.vp(2,i);

    axes(axJ);
    plot([x1,x2],[y1,y2],'color',lncol(i,:))
%     LRmask = zeros(isz.*2);
%     LRmask(:,(isz(2)/2 + 1):end) = 1;
%     LRtx = transform_image(LRmask,[0,0,atan(1/m),1,1,0],'interpmethod','nearest');
%     axAng = atand(1/m);
%     if m < 0
%       axAng = axAng + pi/2;
%     end
%     
%     orthBasis = [1;-1/m];
%     tMaskisz = imresize(tMask,isz,'nearest');
%     tMprop = regionprops(tMaskisz,'centroid');
%     [tMr,tMc] = find(tMaskisz);
%     tst = ones(isz).*2;
%     dval = dot([tMc';tMr'-b],repmat(orthBasis,1,length(tMr)));
%     dval10 = dval > 0;
%     tst(sub2ind(isz,tMr,tMc)) = dval10;
%     figure; imagesc(tst); hold on; plot([x1,x2],[y1,y2],'color',lncol(i,:))
%     IvalLR = [mean(Ih(tst == 0)),mean(Ih(tst == 1))]
%     npxlLR = [sum(tst(:) == 0),sum(tst(:) == 1)]
  end
  
end
end
% figure; imshow(J);

% Copyright © 2014 New York University.
%
% All Rights Reserved. A license to use and copy this software and its documentation
% solely for your internal research and evaluation purposes, without fee and without a signed licensing agreement,
% is hereby granted upon your download of the software, through which you agree to the following:
% 1) the above copyright notice, this paragraph and the following paragraphs
% will prominently appear in all internal copies;
% 2) no rights to sublicense or further distribute this software are granted;
% 3) no rights to modify this software are granted; and
% 4) no rights to assign this license are granted.
% Please Contact The Office of Industrial Liaison, New York University, One Park Avenue, 6th Floor,
% New York, NY 10016 (212) 263-8178, for commercial licensing opportunities,
% or for further distribution, modification or license rights.
%
% Created by Marcelo Cicconet.
%
% IN NO EVENT SHALL NYU, OR ITS EMPLOYEES, OFFICERS, AGENTS OR TRUSTEES (?COLLECTIVELY ?NYU PARTIES?)
% BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES OF ANY KIND ,
% INCLUDING LOST PROFITS, ARISING OUT OF ANY CLAIM RESULTING FROM YOUR USE OF THIS SOFTWARE AND ITS DOCUMENTATION,
% EVEN IF ANY OF NYU PARTIES HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH CLAIM OR DAMAGE.
%
% NYU SPECIFICALLY DISCLAIMS ANY WARRANTIES OF ANY KIND REGARDING THE SOFTWARE,
% INCLUDING, BUT NOT LIMITED TO, NON-INFRINGEMENT, THE IMPLIED WARRANTIES OF  MERCHANTABILITY
% AND FITNESS FOR A PARTICULAR PURPOSE, OR THE ACCURACY OR USEFULNESS,
% OR COMPLETENESS OF THE SOFTWARE. THE SOFTWARE AND ACCOMPANYING DOCUMENTATION,
% IF ANY, PROVIDED HEREUNDER IS PROVIDED COMPLETELY "AS IS".
% NYU HAS NO OBLIGATION TO PROVIDE FURTHER DOCUMENTATION, MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
%
% Please cite the following reference if you use this software in your research:
%
% Marcelo Cicconet, Davi Geiger, Kristin Gunsalus, and Michael Werman.
% Mirror Symmetry Histograms for Capturing Geometric Properties in Images.
% IEEE Conference on Computer Vision and Pattern Recognition. Columbus, Ohio. 2014.