%% Script: manualregGUI
% Description:
%   CODE ADAPTED FROM: C:\Users\Amanda\Documents\My Box Files\MIT Courses\Senior Year, S2009\6.555J\Labs\Lab 4\Lab4DH\showManualComposite.m
%   Overlays the moving and fixed images in the workspace and allows user
%   to manually change the 6 affine registration parameters. The title
%   reflects the mutual information value associated with the
%   transformation parameters specified by the sliders.
% Example:
% Required Functions:
% 
% -------------------------------------------------------------------------
% WORKSPACE VARIALBE REQUIERED
% -------------------------------------------------------------------------
% F:  fixed image
% Fstr: string indicating the label associated with the fixed image
% M:  moving image
% Mstr: string indicating the label associated with the moving image
% mu0: 1x6 vector with initial transformation parameter values
%   (DEFAULT = [xtran = 0, ytran = 0, deg_ang = 0.*-pi/180, xsc = 1, ysc = 1, sk = 0])
%
%  Date           Author            E-mail                      Version
%  ORIGINAL AUTHOR LILLA ZOLLEI for 6.555
%  25 Apr 2013   Amanda Gaudreau   amanda.gaudreau@gmail.com     1
%  20 Oct 2013   Amanda Gaudreau   amanda.gaudreau@gmail.com     2
%  2  Apr 2017   Amanda Gaudreau   amanda.gaudreau@gmail.com     3

% F = double(rgb2gray(PItimms));%SIhipp(:,:,[1,3]);
% M = double(rgb2gray(PItimms));%PIatlas(:,:,4);
% Fstr = 'Fe/Cu/Zn MSI';
% Mstr = 'Atlas';

if ~exist('F')
  F = double(dicomread('knee1.dcm'));
end

if ~exist('M')
  M = double(dicomread('knee2.dcm'));
end

if ~exist('mu0')
  mu0 = [0,0,0,1,1,0]; %tx, ty, theta, sx, sy, sk;
end

if ~exist('Fstr')
  Fstr = '';
end

if ~exist('Mstr')
  Mstr = '';
end

% Load the data and set some of the parameters
T0 = maketform('affine',[1 0; 0 1; 0 0]);
orientation = 3;

% addpath /afs/athena.mit.edu/course/6/6.555/data/reg
maxData = max( max(F(:)), max(M(:)) );
F = F./max(F(:)); %F / maxData;
M = M./max(M(:)); %M / maxData;

sz = size(F);
midslices = round(sz/2);

% Setting GUI parameters
min_dx = mu0(1)-round(size(F,2)*0.2);%-20;%
max_dx = mu0(1)+round(size(F,2)*0.2);%20;%
def_dx = mu0(1);
range_dx = max_dx - min_dx;

min_dy = mu0(2)-round(size(F,1)*0.2);%-20;%
max_dy = mu0(2)+round(size(F,1)*0.2);%20;%
def_dy = mu0(2);
range_dy = max_dy - min_dy;

min_rot = mu0(3)*180/pi-20;
max_rot = mu0(3)*180/pi+20;
def_rot = mu0(3)*180/pi;
range_rot = max_rot - min_rot;

min_sx = 0.6;
max_sx = 1.5;
def_sx = mu0(4);
range_sx = max_sx - min_sx;

min_sy = 0.6;
max_sy = 1.5;
def_sy = mu0(5);
range_sy = max_sy - min_sy;

min_sk = -0.2;
max_sk = 0.2;
def_sk = mu0(6);
range_sk = max_sk - min_sk;

Medgth = .05;%.15;
usiF = [];
reF = (F - min(F(:)))./(max(F(:)) - min(F(:)));
if length(size(F)) > 2
  ch1 = reF(:,:,1); ch2 = reF(:,:,2);
  [Fm,Fn] = size(F(:,:,1));
  if size(F,3) >= 3
    ch3 = reF(:,:,3);
  else
    ch3 = zeros(Fm,Fn);
  end
  jetF = [ch1(:),ch2(:),ch3(:)];
  if any(ch1(:)>min(ch1(:)))
    usiF = [usiF,1];
  end
  if any(ch2(:)>min(ch2(:)))
    usiF = [usiF,2];
  end
  if any(ch3(:) > min(ch3(:)))
    usiF = [usiF,3];
  end
else
  [Fm,Fn] = size(F);
  jetmap = parula(256);
  reF = round(((F - min(F(:)))./(max(F(:)) - min(F(:)))).*255)./255;
  jetF = jetmap(reF(:).*255+1,:);
  usiF = 1;
  %ch1 = reF; ch2 = reF; ch3 = reF;
end
% ch1(find(Mtedge == 1)) = 1;
% ch2(find(Mtedge == 1)) = 0;
% ch3(find(Mtedge == 1)) = 0;
% Fimg3 = repmat(ch1,[1,1,3]);
% Fimg3(:,:,1) = ch1;
% Fimg3(:,:,2) = ch2;
% Fimg3(:,:,3) = ch3;
% image(Fimg3);
% 
% jetmap = jet(256);
% reF = round(((F - min(F(:)))./(max(F(:)) - min(F(:)))).*255)./255;
% jetF = jetmap(reF(:).*255+1,:);
[MIorig,info] = twoimgMIkde(F(:,:,usiF(1)),M(:,:,1),'mu',[def_dx, def_dy, def_rot*pi/180, def_sx, def_sy, def_sk],...
  'extrapval',nan,'kernel',{'epan'},'nbins',20,'sumdim',1);
titlestr2 = sprintf('Original MI: %1.4f, Original mu = [%1.1f,%1.1f,%1.1f,%1.1f,%1.1f,%1.1f]',...
  MIorig, def_dx, def_dy, def_rot, def_sx, def_sy, def_sk);
titlestr1 = '';

hman = figure;
callbackString = [ 'T = [ get( dx_slider, ''value'' ),  get( dy_slider, ''value'' ),'...
  'get( rot_slider, ''value'' )*pi/180, get( sx_slider, ''value'' ),'...
  'get( sy_slider, ''value'' ), get( sk_slider, ''value'' )]; ' ...
  'moved = transform_image(M, T,''extrapval'',nan); ' ...
  'Emoved = transform_image(M, T,''extrapval'',0); ' ...
  'if length(size(Emoved)) > 2; Emoved = sum(Emoved.*1/size(Emoved,3),3); end;' ...
  'im1 = F; moved_im2 = moved;'...
  'set(dx_str, ''string'', [''dx = '' num2str(T(1), 3)] );'...def_dx = T(1);' ...
  'set(dy_str, ''string'', [''dy = '' num2str(T(2), 3)] );'...def_dx = T(1);' ...
  'set(rot_str, ''string'', [''rot = '' num2str(T(3)*180/pi, 3)] );'...def_dx = T(1);' ...
  'set(sx_str, ''string'', [''sx = '' num2str(T(4), 3)] );'...def_dx = T(1);' ...
  'set(sy_str, ''string'', [''sy = '' num2str(T(5), 3)] );'...
  'set(sk_str, ''string'', [''sk = '' num2str(T(6), 3)] );'...
  'imoverlay(im1,'''',moved_im2,moved_im2,'''',''displaymode'',''checker'',''blocksz'',20,''npanels'',1,''figh'',hman); axis image;'...
  'set(gca,''xticklabel'','''',''yticklabel'','''');'...
  'if get(compMI,''value''); titlestr1 = sprintf(''I_F: %s, I_M: %s, MI = %1.4f, T = [%1.1f,%1.1f,%1.1f,%1.1f,%1.1f,%1.1f]'''...
  ',Fstr, Mstr,twoimgMIkde(F(:,:,1),moved(:,:,1),''kernel'',{''epan''},''nbins'',20,''sumdim'',1),T); end;'...
  'title({[titlestr1],[titlestr2]});'];
% 'imoverlay(im1,'''',moved_im2,moved_im2,'''',''displaymode'',''rgbgray'',''npanels'',1,''transparency'',0.8,''figh'',hman);'...
%   ',Fstr, Mstr,twoimgMIkde(F,moved,''kernel'',{''epan''},''nbins'',30),T)];'...
%   '[sprintf(''Original MI: %1.4f, Original mu = [%1.1f,%1.1f,%1.1f,%1.1f,%1.1f,%1.1f]'','...
%   'MIorig, def_dx, def_dy, def_rot, def_sx, def_sy, def_sk)]}); end' ];

compMI = uicontrol('style','pushbutton','string','Compute MI','Position', [460 20 75 20],...
  'callback',callbackString);

% Displacement sliders (dx, dy, rot)
dx_slider = uicontrol( 'style', 'slider', ...
  'min', min_dx, 'max', max_dx, 'value', def_dx, ...
  'sliderstep', [.1/range_dx 1/range_dx], ...
  'position', [120 32 250 14], ...
  'callback', callbackString );
dy_slider = uicontrol( 'style', 'slider', ...
  'min', min_dy, 'max', max_dy, 'value', def_dy, ...
  'sliderstep', [.1/range_dy 1/range_dy], ...
  'position', [120 16 250 14], ...
  'callback', callbackString );
rot_slider = uicontrol( 'style', 'slider', ...
  'min', min_rot, 'max', max_rot, 'value', def_rot, ...
  'sliderstep', [.1/range_rot 1/range_rot], ...
  'position', [120 1 250 14], ...
  'callback', callbackString );

% angles_str   = uicontrol( 'style', 'text', ...
%     'string', 'Angles:' , ...
%     'position', [25  330 60 18], ...
%                'callback', callbackString );
sx_slider = uicontrol( 'style', 'slider', 'String', 'Xscale', ...
  'min', min_sx, 'max', max_sx, 'value', def_sx, ...
  'sliderstep', [0.01 0.1], ...
  'position', [30 78 15 250], ...
  'callback', callbackString );
sy_slider = uicontrol( 'style', 'slider', 'String', 'Yscale', ...
  'min', min_sy, 'max', max_sy, 'value', def_sy, ...
  'sliderstep', [0.01 0.1], ...
  'position', [50 78 15 250], ...
  'callback', callbackString );
sk_slider = uicontrol( 'style', 'slider', 'String', 'skew', ...
  'min', min_sk, 'max', max_sk, 'value', def_sk, ...
  'sliderstep', [0.001 0.01], ...
  'position', [70 78 15 250], ...
  'callback', callbackString );

% GUI strings
dx_str   = uicontrol( 'style', 'text', ...
  'string', ['dx = ' num2str(def_dx)] , ...
  'position', [375 32 80 14], ...
  'callback', callbackString );
dy_str   = uicontrol( 'style', 'text', ...
  'string', ['dy = ' num2str(def_dy)] , ...
  'position', [375 16 80 14], ...
  'callback', callbackString );
rot_str   = uicontrol( 'style', 'text', ...
  'string', ['rot = ' num2str(def_rot)] , ...
  'position', [375 1 80 14], ...
  'callback', callbackString );

sx_str   = uicontrol( 'style', 'text', ...
  'string', ['sx = ' num2str(def_sx)] , ...
  'position', [17 60 80 15], ...
  'callback', callbackString );
sy_str   = uicontrol( 'style', 'text', ...
  'string', ['sy = ' num2str(def_sy)] , ...
  'position', [17 44 80 15], ...
  'callback', callbackString );
sk_str   = uicontrol( 'style', 'text', ...
  'string', ['sk = ' num2str(def_sk)] , ...
  'position', [17 28 80 15], ...
  'callback', callbackString );

% 3D volume orientation menu
% orient_menu   = uicontrol( 'style', 'popup', ...
%   'string', 'coronal|sagittal|axial' , ...
%   'value', orientation , ...
%   'position', [20 370 80 20], ...
%   'callback', callbackString );

eval( callbackString );





