function varargout = sa_analysisgui(varargin)
% SA_ANALYSISGUI MATLAB code for sa_analysisgui.fig
%      SA_ANALYSISGUI, by itself, creates a new SA_ANALYSISGUI or raises the existing
%      singleton*.
%
%      H = SA_ANALYSISGUI returns the handle to a new SA_ANALYSISGUI or the handle to
%      the existing singleton*.
%
%      SA_ANALYSISGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SA_ANALYSISGUI.M with the given input arguments.
%
%      SA_ANALYSISGUI('Property','Value',...) creates a new SA_ANALYSISGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before sa_analysisgui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to sa_analysisgui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help sa_analysisgui

% Last Modified by GUIDE v2.5 11-Mar-2013 18:06:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @sa_analysisgui_OpeningFcn, ...
                   'gui_OutputFcn',  @sa_analysisgui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before sa_analysisgui is made visible.
function sa_analysisgui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to sa_analysisgui (see VARARGIN)

% Choose default command line output for sa_analysisgui
handles.output = hObject;
fpath = 'D:\My Documents\Dropbox\MADLab Research\Data\simanneal Experiments';%'C:\Users\Amanda\Documents\Dropbox\MADLab Research\Data\simanneal Experiments';
set(handles.pathtxt,'String',strcat(fpath,'\*.csv'));
cbinfo.CallbackCode = 'init';
updatefilelist(handles,cbinfo)

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes sa_analysisgui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = sa_analysisgui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
% --- Executes on selection change in mmddpopup.
function mmddpopup_Callback(hObject, eventdata, handles)
cbinfo.CallbackCode = 'mmdd';
cbinfo = getpopupinfo(handles,cbinfo);
updatefilelist(handles,cbinfo)
% --- Executes during object creation, after setting all properties.
function mmddpopup_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on selection change in nvarpopup.
function nvarpopup_Callback(hObject, eventdata, handles)
cbinfo.CallbackCode = 'nvar';
cbinfo = getpopupinfo(handles,cbinfo);
updatefilelist(handles,cbinfo)
% --- Executes during object creation, after setting all properties.
function nvarpopup_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on selection change in T0popup.
function T0popup_Callback(hObject, eventdata, handles)
cbinfo.CallbackCode = 'temp';
cbinfo = getpopupinfo(handles,cbinfo);
updatefilelist(handles,cbinfo)
% --- Executes during object creation, after setting all properties.
function T0popup_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on selection change in varbinarypopup.
function varbinarypopup_Callback(hObject, eventdata, handles)
cbinfo.CallbackCode = 'varbin';
cbinfo = getpopupinfo(handles,cbinfo);
updatefilelist(handles,cbinfo)
% --- Executes during object creation, after setting all properties.
function varbinarypopup_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on button press in filebrowsebtn.
function filebrowsebtn_Callback(hObject, eventdata, handles)
% hObject    handle to filebrowsebtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fpath = get(handles.pathtxt,'String');
fname = uigetfile(strcat(fpath,'\*.csv'));
set(handles.pathtxt,'String',strcat(fpath,'\',fname));
cbinfo.CallbackCode = 'browse';
updatefilelist(handles,cbinfo)
% --- Executes on selection change in varbinarypopup.
function pathtxt_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function pathtxt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function cbinfo = getpopupinfo(handles,cbinfo)
cbinfo.varbini = get(handles.varbinarypopup,'value');
cbinfo.varbinstr = get(handles.varbinarypopup,'string');
cbinfo.varbinsel = cbinfo.varbinstr{cbinfo.varbini};
cbinfo.T0i = get(handles.T0popup,'value');
cbinfo.T0str = get(handles.T0popup,'string');
cbinfo.T0sel = cbinfo.T0str{cbinfo.T0i};
cbinfo.nvari = get(handles.nvarpopup,'value');
cbinfo.nvarstr = get(handles.nvarpopup,'string');
cbinfo.nvarsel = cbinfo.nvarstr{cbinfo.nvari};
cbinfo.mmddi = get(handles.mmddpopup,'value');
cbinfo.mmddstr = get(handles.mmddpopup,'string');
cbinfo.mmddsel = cbinfo.mmddstr{cbinfo.mmddi};

% --- Executes during object creation, after setting all properties.
function updatefilelist(handles,cbinfo)
mmddslashi = 6; % end index for MM_DD string at beginning of file names
modstr = '';
% Parse out file path from existing handles
currentdir = 'D:\My Documents\Dropbox\MADLab Research\Data\simanneal Experiments';%'C:\Users\Amanda\Documents\Dropbox\MADLab Research\Data\simanneal Experiments\';%uigetdir(); get(handles.pathtxt,'string');
bkslashi = strfind(currentdir,'\');
currentdir = currentdir(1:bkslashi(end)-1);
if isdir(currentdir)
  fpath = currentdir;
else
  fpath = 'D:\My Documents\Dropbox\MADLab Research\Data\simanneal Experiments';%'C:\Users\Amanda\Documents\Dropbox\MADLab Research\Data\simanneal Experiments';%uigetdir();
end

switch cbinfo.CallbackCode
  case 'init'
    cbinfo.nvarsel = '';
    cbinfo.mmddsel = '';
    cbinfo.T0sel = '';
    cbinfo.varbinsel = '';
  case 'browse'
    
  case 'mmdd'
    modstr = cbinfo.mmddsel;
  case 'temp'
    modstr = cbinfo.T0sel;
    if rem(str2num(modstr),1)
      modstr = regexprep(modstr,'.','p');
    end
    modstr = sprintf('_T%s',modstr);
  case 'nvar'
    modstr = strcat('nvar',cbinfo.nvarsel);
  case 'varbin'
    modstr = cbinfo.varbinsel;
end

F = ls(fpath);
Fcell = cellstr(F);
% Isolate .csv files
Fcsvi = strfind(Fcell,'.csv');
relF = cellfun(@(x) isempty(x) == 0,Fcsvi);
F = F(relF,:);
Fcell = Fcell(relF);
Fcsvi = Fcsvi(relF);
% Isolate modifier files selected by the user
if ~isempty(modstr)
  Fmodi = strfind(Fcell,modstr);
  relF = cellfun(@(x) isempty(x) == 0,Fmodi);
  F = F(relF,:);
  Fcell = Fcell(relF);
  Fcsvi = Fcsvi(relF);
end

% Fill up the popup menus so that they reflect the user's input
Fnvari = strfind(Fcell,'nvar');
relF = cellfun(@(x) isempty(x) == 0,Fnvari);
F = F(relF,:);
Fcell = Fcell(relF);
Fnvari = Fnvari(relF);
Fcsvi = Fcsvi(relF);
Fnvari = cell2mat(Fnvari) + 4;
nvall = F(sub2ind(size(F),[1:size(F,1)]',Fnvari));
nvunq = cellstr(unique(nvall,'rows'));
nvunq = [{''}; nvunq];
if length(nvunq) == 2
  set(handles.nvarpopup,'String',nvunq,'value',2);
elseif isempty(cbinfo.nvarsel)
  set(handles.nvarpopup,'String',nvunq,'value',1);
else
  nvarval = strfind(nvunq,cbinfo.nvarsel);
  nvarval = find(cellfun(@(x) ~isempty(x),nvarval));
  set(handles.nvarpopup,'String',nvunq,'value',nvarval);
end

Fslashi = strfind(Fcell,'_');
mmddall = F(:,1:mmddslashi-1);
mmddunq = cellstr(unique(mmddall,'rows'));
mmddunq = [{''}; mmddunq];
if length(mmddunq) == 2
  set(handles.mmddpopup,'String',mmddunq,'value',2);
elseif isempty(cbinfo.mmddsel)
  set(handles.mmddpopup,'String',mmddunq,'value',1);
else
  mmddval = strfind(mmddunq,cbinfo.mmddsel);
  mmddval = find(cellfun(@(x) ~isempty(x),mmddval));
  set(handles.mmddpopup,'String',mmddunq,'value',mmddval);
end

Ftempi = strfind(Fcell,'_T');
tempall = cell(length(Ftempi),1);
bincodeTF = zeros(length(Ftempi),1);
bincodeall = cell(length(Ftempi),1);
for t = 1:length(Ftempi)
  eqslashi = find(Fslashi{t} == Ftempi{t});
  if eqslashi == length(Fslashi{t})
    tempstr = F(t,(Ftempi{t}+2):(Fcsvi{t}-1));
    tempall{t}= regexprep(regexprep(tempstr,'p','.'),'i','');
    bincodeTF(t) = 0;
    bincodeall{t} = 'none';
  else
    tempstr = F(t,(Ftempi{t}+2):(Fslashi{t}(eqslashi+1)-1))
    tempall{t} = regexprep(regexprep(tempstr,'p','.'),'i','');
    bincodeTF(t) = 1;
    bincodeall{t} = F(t,Fslashi{t}(eqslashi+1)+(1:6));
  end
end
tempunq = unique(tempall);
[tsort,tind] = sort(cellfun(@(x) str2num(x),tempunq));
tempunq = tempunq(tind);
tempunq = [{''}; tempunq];
if length(tempunq) == 2
  set(handles.T0popup,'String',tempunq,'value',2);
elseif isempty(cbinfo.T0sel)
  set(handles.T0popup,'String',tempunq,'value',1);
else
  T0val = strfind(tempunq,cbinfo.T0sel);
  T0val = find(cellfun(@(x) ~isempty(x),T0val));
  set(handles.T0popup,'String',tempunq,'value',T0val);
end

bincodeunq = unique(bincodeall);
bincodeunq = [{''}; bincodeunq];
if length(bincodeunq) == 2
  set(handles.varbinarypopup,'String',bincodeunq,'value',2);
elseif isempty(cbinfo.varbinsel)
  set(handles.varbinarypopup,'String',bincodeunq,'value',1);
else
  vbinval = strfind(bincodeunq,cbinfo.varbinsel);
  vbinval = find(cellfun(@(x) ~isempty(x),vbinval));
  set(handles.varbinarypopup,'String',bincodeunq,'value',vbinval);
end

if size(F,1) == 1
  set(handles.pathtxt,'String',F,'value',1);
else
  Fpossible = [{''}; cellstr(F)];
  set(handles.pathtxt,'String',Fpossible,'value',1);
end


% --- Executes on button press in plotslice.
function plotslice_Callback(hObject, eventdata, handles)
% hObject    handle to plotslice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotslice


% --- Executes on button press in plot3d.
function plot3d_Callback(hObject, eventdata, handles)
% hObject    handle to plot3d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plot3d
