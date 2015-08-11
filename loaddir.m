%% Script: loaddir
% Description: Loads appropriate research directories for the computer
% INPUTS ----------------------------------------------------------------
%
% OUTPUTS ---------------------------------------------------------------
% comp    Tries to find which computer is begin used based on the parent
%         directly of the "My Documents" folder. This assumes that the
%         default path upon statrtup of MATLAB is somewhere in the "My
%         Documents" folder
% dboxpath  Path of the dropbox folder (DOES NOT end in \)
% fpath     Path of the Data folder (ends in \)
% mtngpath  Latex folder -- needed for saving pngs to use in tex docs
% img_path  path of the Cu grid image used in ground truth registration
%           experiments
% cugridpath ''
% datapath  Mutual Information experimental data path (ends in \)
% sapath    Simmulated annealing experimental data path (ends in \)
% lafldr    laser ablation data folder (ends in \)
% mlabpath  MATLAB path where research scripts live (ends in \)
% antspath  Path of ANTS bin folder (DOES NOT end in \)
% EBfldr    Path of EB RGB brain images (DOES NOT end in \)
% dateyyyymmdd 6x1 double with datainfo
% dateyyyymmddstr string with data info
%
%  Date           Author            E-mail                      Version
%  11 Mar  2014   Amanda Gaudreau   amanda.gaudreau@gmail.com     0
%  19 Nov  2014   Amanda Balderrama amanda.gaudreau@gmail.com     1
%     Changed EBfldr: moved images from dboxpath to mdatapath

set(0,'defaultfigurecolor','w','DefaultAxesFontSize',12,...'DefaultAxesFontWeight','bold','DefaultAxesLineStyleOrder','-|--|:|-.',...
  'DefaultLineLineWidth',2,... 'Defaultaxesposition','remove',...
  'defaultAxesFontName','Arial');

% comp = 'TENN1_000'; % MED comp %'ADGB'; % ThinkPad; %
currentfldr = cd;
Dpos = strfind(currentfldr,'Documents');
comp = currentfldr(strfind(currentfldr,'Users\')+6:Dpos-2);%'TENN1_000'; % MED comp %'ADGB'; % ThinkPad; %
dboxpath = strcat(currentfldr(1:Dpos-1),'Documents\Dropbox');
fpath = strcat(dboxpath,'\MADLab Research\Data\');
mtngpath = strcat(dboxpath,'\MADLab Research\Literature\Latex Citation Summary\');
img_path = strcat(dboxpath,'\MADLab Research\Data\Cugrid Images\Photos\12-3-27 G50 Copper Grid Black background3 1212.jpg');
cugridpath = strcat(dboxpath,'\MADLab Research\Data\Cugrid Images\Photos\12-3-27 G50 Copper Grid Black background3 1212.jpg');
datapath = strcat(dboxpath,'\MADLab Research\Data\Mutual Information\');
sapath = strcat(dboxpath,'\MADLab Research\Data\simanneal Experiments\');
lafldr = strcat(currentfldr(1:Dpos-1),'Documents\MADLab Data\LaserAblation\');
mdatapath = strcat(currentfldr(1:Dpos-1),'Documents\MADLab Data\');
dbEBexppath = strcat(fpath,'EB Experiments\');
EBexppath = strcat(mdatapath,'EB_Experiments\');
ivisfldr = strcat(EBexppath,'IVIS_EBFluor\');
mlabpath = strcat(dboxpath,'\MADLab Research\MATLAB\MIMSlib');
antspath = 'C:\Program Files (x86)\ANTS\bin';
EBfldr = strcat(EBexppath,'RGBmacroscopic');%strcat(fpath,'EB_brainimg');
rgbEBfldr = strcat(EBexppath,'RGBmacroscopic');%strcat(fpath,'EB_brainimg');
draperdata = strcat(mlabpath,'draper\data\');
dateyyyymmdd = datevec(now);
dateyyyymmddstr = sprintf('%d%02d%02d',dateyyyymmdd(1:3));
