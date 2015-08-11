function varargout = ablation(varargin)
% ABLATION M-file for ablation.fig
%      ABLATION, by itself, creates a new ABLATION or raises the existing
%      singleton*.
%
%      H = ABLATION returns the handle to a new ABLATION or the handle to
%      the existing singleton*.
%
%      ABLATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ABLATION.M with the given input arguments.
%
%      ABLATION('Property','Value',...) creates a new ABLATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ablation_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ablation_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help ablation

% Last Modified by GUIDE v2.5 06-Apr-2010 23:59:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ablation_OpeningFcn, ...
                   'gui_OutputFcn',  @ablation_OutputFcn, ...
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


% --- Executes just before ablation is made visible.
function ablation_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ablation (see VARARGIN)

% Choose default command line output for ablation
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ablation wait for user response (see UIRESUME)
% uiwait(handles.base_element);
types = {'Empty Contour','Filled Contour','Pixel Image'};
set(handles.type,'String',types);
global schemes;
schemes = {'gray','hot','jet','hsv','copper'};
set(handles.scheme,'String',schemes);
plot_types = {'Linear Intensity','Logarithmic Intensity'};
set(handles.plot_type,'String',plot_types);
global massdata elements base_value;

% --- Outputs from this function are returned to the command line.
function varargout = ablation_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function start_Callback(hObject, eventdata, handles)
% hObject    handle to start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of start as text
%        str2double(get(hObject,'String')) returns contents of start as a double


% --- Executes during object creation, after setting all properties.
function start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function base_element_Callback(hObject, eventdata, handles)
% hObject    handle to base_element (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of base_element as text
%        str2double(get(hObject,'String')) returns contents of base_element as a double
global massdata elements base_value;
base_value = get(handles.elements,'Value');
base = elements{base_value};
    

    

% --- Executes during object creation, after setting all properties.
function base_element_CreateFcn(hObject, eventdata, handles)
% hObject    handle to base_element (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on selection change in elements.
function elements_Callback(hObject, eventdata, handles)
% hObject    handle to elements (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns elements contents as cell array
%        contents{get(hObject,'Value')} returns selected item from elements
global massdata elements base_value;
element_choice = get(handles.elements,'Value');
hold(handles.raw_plot,'off');
plot(handles.raw_plot,massdata(:,1),massdata(:,element_choice+1));
drawnow;
% --- Executes during object creation, after setting all properties.
function elements_CreateFcn(hObject, eventdata, handles)
% hObject    handle to elements (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in scheme.
function scheme_Callback(hObject, eventdata, handles)
% hObject    handle to scheme (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns scheme contents as cell array
%        contents{get(hObject,'Value')} returns selected item from scheme


% --- Executes during object creation, after setting all properties.
function scheme_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scheme (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function width_Callback(hObject, eventdata, handles)
% hObject    handle to width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of width as text
%        str2double(get(hObject,'String')) returns contents of width as a double


% --- Executes during object creation, after setting all properties.
function width_CreateFcn(hObject, eventdata, handles)
% hObject    handle to width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function num_Callback(hObject, eventdata, handles)
% hObject    handle to num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num as text
%        str2double(get(hObject,'String')) returns contents of num as a double


% --- Executes during object creation, after setting all properties.
function num_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function inter_Callback(hObject, eventdata, handles)
% hObject    handle to inter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of inter as text
%        str2double(get(hObject,'String')) returns contents of inter as a double


% --- Executes during object creation, after setting all properties.
function inter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in type.
function type_Callback(hObject, eventdata, handles)
% hObject    handle to type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from type


% --- Executes during object creation, after setting all properties.
function type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in clear.
function clear_Callback(hObject, eventdata, handles)
% hObject    handle to clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla(handles.map);



% --- Executes on button press in load.
function load_Callback(hObject, eventdata, handles)
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global massdata elements base_value;
current_location = cd;
set(handles.status,'String','Select one Data file.');
drawnow;
[filename,pathname] = uigetfile;
cd(pathname);
set(handles.status,'String','Merging Files...');
drawnow;
allfiles = struct2cell(dir);
t = 0;
y = 1;
[sizesize sizez] = size(allfiles);
for i = 1:sizez
    fid = fopen(allfiles{1,i});
    if(fid>0)
        set(handles.status,'String',strcat('Merging Files...',num2str(y)));
        drawnow;
        y = y+1;
        data1 = textscan(fid,'%s');
        data2 = data1{1};
        for j = 1:length(data2)
            if(~isempty(strfind(data2{j},'Time,')))
                data3 = data2{j};
                elements_c = data3(6:length(data3));
                begin = j+1;
                break;
            end
        end
        for k = begin:length(data2)
            t = t+1;
            datas{t} = data2{k};
        end
    end
end
set(handles.status,'String','Loading Data...');
drawnow;
cd(current_location);
fid = fopen('massdata_c.txt','w+');
for i = 1:length(datas)
    fprintf(fid,'%s\n',datas{i});
end
massdata = csvread('massdata_c.txt');
time_base = massdata(10,1)-massdata(9,1);
time = 0:time_base:length(massdata)*time_base-time_base;
massdata(:,1) = time;
save massdata.txt -ASCII -TABS massdata;
t = 1;
commas = strfind(elements_c,',');
for i = 1:length(commas)+1
    if(i==1)
        elements{t} = elements_c(i:commas(i)-1);
        t = t + 1;
    elseif(i==length(commas)+1)
        elements{t} = elements_c(commas(i-1)+1:length(elements_c));    
    else
        elements{t} = elements_c(commas(i-1)+1:commas(i)-1);
        t = t + 1;
    end
end
fid2 = fopen('elements_c.txt','w+');
for i = 1:length(elements)
    fprintf(fid2,'%s\t',elements{i});
end
fclose('all');
fid3 = fopen('elements_c.txt','r+');
set(handles.elements,'String',elements(1:length(elements)));
set(handles.status,'String','Data Ready.');




% --- Executes on button press in base.
function base_Callback(hObject, eventdata, handles)
% hObject    handle to base (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)





function t_max_Callback(hObject, eventdata, handles)
% hObject    handle to t_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t_max as text
%        str2double(get(hObject,'String')) returns contents of t_max as a double


% --- Executes during object creation, after setting all properties.
function t_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function t_min_Callback(hObject, eventdata, handles)
% hObject    handle to t_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t_min as text
%        str2double(get(hObject,'String')) returns contents of t_min as a double


% --- Executes during object creation, after setting all properties.
function t_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in set_t_axis.
function set_t_axis_Callback(hObject, eventdata, handles)
% hObject    handle to set_t_axis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmin = str2double(get(handles.t_min,'String'));
tmax = str2double(get(handles.t_max,'String'));
set(handles.raw_plot,'XLim',[tmin tmax]);



% --- Executes on button press in map_button.
function map_button_Callback(hObject, eventdata, handles)
% hObject    handle to map_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global massdata elements base_value schemes d;
set(handles.status,'String','Preparing...');
massdata_orig = massdata;
drawnow;
i = get(handles.elements,'Value');
i = i+1;
time_base = massdata(11,1)-massdata(10,1);
start = round(str2double(get(handles.start,'String'))/time_base);
width = round(str2double(get(handles.width,'String'))/time_base);
inter = round(str2double(get(handles.inter,'String'))/time_base);
num = str2double(get(handles.num,'String'));
[a,b] = size(massdata);
avg = mean(massdata(:,i));
%separate the data into lines

    done = 0;
    k = 1;
    g = 1;
    j = start+1;
    set(handles.status,'String','Filtering...');
    drawnow;
    if(get(handles.filter,'Value'))
        massdata(massdata>str2double(get(handles.thresh,'String'))*avg) = 0;
    end
    set(handles.status,'String','Calculating...');
    while(done == 0)
        while(k<=width)
            line(k,g,i-1) = massdata(j,i);
            line_index(k,g,i-1) = j;
            k = k + 1;
            j = j + 1;
        end
        g = g+1;
        if(g>num)
            done = 1;
        end
        j = j + inter;
        k = 1;
    end

set(handles.status,'String','Plotting...');
drawnow;
d = line(:,:,i-1);
dd = line_index(:,:,i-1);
d = d';
dd = dd';
%file_save_data = strcat(dir,'\',folder2,'\',elements{1}{i-1},'_data.txt');
ddd = d';
%save(file_save_data,'-ASCII','-TABS','ddd');
plot_type = get(handles.type,'Value');
if(plot_type == 1)
    axes(handles.map);
    figure(get(handles.map,'parent'));
    if(get(handles.plot_type,'Value')==1)
        contour(d);
    else
        contour(log(d));
    end
    colormap(schemes{get(handles.scheme,'Value')});
    drawnow;
elseif(plot_type == 2)
    axes(handles.map);
    figure(get(handles.map,'parent'));
    if(get(handles.plot_type,'Value')==1)
        contourf(d);
    else
        contourf(log(d));
    end
    colormap(schemes{get(handles.scheme,'Value')});
    drawnow;
else
    axes(handles.map);
    figure(get(handles.map,'parent'))
    if(get(handles.plot_type,'Value')==1)
        imagesc(d);
    else
        imagesc(log(d));
    end
    colormap(schemes{get(handles.scheme,'Value')});
    drawnow;
    set(handles.map,'YDir','normal');
    save dd.txt -ASCII -TABS d;
end
if(get(handles.plot_type,'Value')==2)
    title_log = strcat(elements{i-1},'-Logarithmic Intensity');
    title(title_log);
else
    title(elements{i-1});
end
colorbar;
%file_save_image = strcat(dir,'\',folder1,'\',elements{1}{i-1},'_contourplot.jpg');
%saveas(gcf,file_save_image);
%delete(gcf);

[aa,bb] = size(d);
for i = 1:aa
    starts(i) = dd(i,1)*time_base;
    stops(i) = dd(i,bb)*time_base;
end
axes(handles.raw_plot);
plot(handles.raw_plot,massdata(:,1),massdata(:,base_value+1));
hold(handles.raw_plot,'on');
reference = str2double(get(handles.reference,'String'));
if(isnan(reference))
    reference = 1;
end
plot(starts,max(massdata(:,base_value+1))*reference,'go');
plot(stops,max(massdata(:,base_value+1))*reference,'ro');
drawnow;
set(handles.status,'String','Data Ready.');
hold off;
massdata = massdata_orig;


% --- Executes on button press in filter.
function filter_Callback(hObject, eventdata, handles)
% hObject    handle to filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of filter



function thresh_Callback(hObject, eventdata, handles)
% hObject    handle to thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thresh as text
%        str2double(get(hObject,'String')) returns contents of thresh as a double


% --- Executes during object creation, after setting all properties.
function thresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in ref.
function ref_Callback(hObject, eventdata, handles)
% hObject    handle to ref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global massdata elements base_value;
time_base = massdata(11,1)-massdata(10,1);
start = round(str2double(get(handles.start,'String'))/time_base);
width = round(str2double(get(handles.width,'String'))/time_base);
inter = round(str2double(get(handles.inter,'String'))/time_base);
num = str2double(get(handles.num,'String'));
[a,b] = size(massdata);
for i = 2:b
    done = 0;
    k = 1;
    g = 1;
    j = start + 1;
    while(done == 0)
        while(k<=width)
            line(k,g,i-1) = massdata(j,i);
            line_index(k,g,i-1) = j;
            k = k + 1;
            j = j + 1;
        end
        g = g+1;
        if(g>num)
            done = 1;
        end
        j = j + inter;
        k = 1;
    end
end

for i = 2:b
    d = line(:,:,i-1);
    dd = line_index(:,:,i-1);
    d = d';
    dd = dd';
    end
[aa,bb] = size(d);
for i = 1:aa
    starts(i) = dd(i,1)*time_base;
    stops(i) = dd(i,bb)*time_base;
end
axes(handles.raw_plot);
plot(handles.raw_plot,massdata(:,1),massdata(:,base_value+1));
hold(handles.raw_plot,'on');
reference = str2double(get(handles.reference,'String'))/100;
if(isnan(reference))
    reference = .01;
end
plot(starts,max(massdata(:,base_value+1))*reference,'go');
plot(stops,max(massdata(:,base_value+1))*reference,'ro');
drawnow;
hold off;


% --- Executes on button press in plot_all.
function plot_all_Callback(hObject, eventdata, handles)
% hObject    handle to plot_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global massdata elements base_value schemes;
set(handles.status,'String','Plotting...');
drawnow;
time_base = massdata(11,1)-massdata(10,1);
start = round(str2double(get(handles.start,'String'))/time_base);
width = round(str2double(get(handles.width,'String'))/time_base);
inter = round(str2double(get(handles.inter,'String'))/time_base);
num = str2double(get(handles.num,'String'));
dir = cd;
[a,b] = size(massdata);
for i = 2:b
    done = 0;
    k = 1;
    g = 1;
    j = start+1;
    while(done == 0)
        while(k<=width)
            if(get(handles.filter,'Value'))
                if(massdata(j,i)>str2double(get(handles.thresh,'String'))*avg(i))
                    massdata(j,i) = 0;
                end
            end
            line(k,g,i-1) = massdata(j,i);
            line_index(k,g,i-1) = j;
            k = k + 1;
            j = j + 1;
        end
        g = g+1;
        if(g>num)
            done = 1;
        end
        j = j + inter;
        k = 1;
    end
end
try_again = 1;
for i = 2:b
    d = line(:,:,i-1);
    dd = line_index(:,:,i-1);
    d = d';
    dd = dd';
    folder2 = get(handles.folder2,'String');
    if(isempty(folder2))
        set(handles.status,'String','Enter File Name.');
        drawnow;
        try_again = 0;
        break
    end
    mkdir(folder2);
    warning off all;
    file_save_data = strcat(dir,'\',folder2,'\',elements{i-1},'_data.txt');
    ddd = d';
    plot_type = get(handles.type,'Value');
    if(plot_type == 1)
        figure;
        if(get(handles.plot_type,'Value')==1)
            contour(d);
        else
            contour(log(d));
        end
        colormap(schemes{get(handles.scheme,'Value')});
        drawnow;
        colorbar;
    elseif(plot_type == 2)
        figure;
        if(get(handles.plot_type,'Value')==1)
            contourf(d);
        else
            contourf(log(d));
        end
        colormap(schemes{get(handles.scheme,'Value')});
        drawnow;
        colorbar;
    else
        figure;
        if(get(handles.plot_type,'Value')==1)
            imagesc(d);
        else
            imagesc(log(d));
        end
        colormap(schemes{get(handles.scheme,'Value')});
        colorbar;
        set(gca,'YDir','normal');
        drawnow;
    end
    if(get(handles.plot_type,'Value')==2)
        title_log = strcat(elements{i-1},'-Logarithmic Intensity');
        title(title_log);
    else
        title(elements{i-1});
    end
    save(file_save_data,'-ASCII','-TABS','d');
    file_save_image = strcat(dir,'\',folder2,'\',elements{i-1},'_plot.jpg');
    saveas(gcf,file_save_image);
    delete(gcf);
end
if(try_again)
    [aa,bb] = size(d);
    for i = 1:aa
        starts(i) = dd(i,1);
        stops(i) = dd(i,bb);
    end
    figure;
    plot(massdata(:,7));
    hold on;
    reference = str2double(get(handles.reference,'String'))/100;
    if(isnan(reference))
        reference = 1;
    end
    plot(starts,max(massdata(:,base_value+1))*reference,'go');
    plot(stops,max(massdata(:,base_value+1))*reference,'ro');
    file_save_image = strcat(dir,'\',folder2,'\','reference_points.jpg');
    title('Reference Points: Green/Red circles define the Data Considered');
    saveas(gcf,file_save_image);
    hold off;
    delete(gcf);
end
set(handles.status,'String','Data Ready.');



function folder2_Callback(hObject, eventdata, handles)
% hObject    handle to folder2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of folder2 as text
%        str2double(get(hObject,'String')) returns contents of folder2 as a double


% --- Executes during object creation, after setting all properties.
function folder2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to folder2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in save_plot.
function save_plot_Callback(hObject, eventdata, handles)
% hObject    handle to save_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global elements d schemes;
set(handles.status,'String','Saving Current Plot...');
drawnow;
dir = cd;
plot_type = get(handles.type,'Value');
element_selection = get(handles.elements,'value');
if(plot_type == 1)
    figure;
    if(get(handles.plot_type,'Value')==1)
        contour(d);
    else
        contour(log(d));
    end
    colormap(schemes{get(handles.scheme,'Value')});
    drawnow;
elseif(plot_type == 2)
    figure;
    if(get(handles.plot_type,'Value')==1)
        contourf(d);
    else
        contourf(log(d));
    end
    colormap(schemes{get(handles.scheme,'Value')});
    drawnow;
else
    figure;
    if(get(handles.plot_type,'Value')==1)
        imagesc(d);
    else
        imagesc(log(d));
    end
    colormap(schemes{get(handles.scheme,'Value')});
    drawnow;
    set(handles.map,'YDir','normal');
end
if(get(handles.plot_type,'Value')==2)
    title_log = strcat(elements{element_selection},'-Logarithmic Intensity');
    title(title_log);
else
    title(elements{element_selection});
end
colorbar;
ddd = d';
current_time = clock;
file_save_data = strcat(elements(element_selection),'_data',num2str(current_time(4)),num2str(current_time(5)));
save(file_save_data{1},'-ASCII','-TABS','ddd');
file_save_image = strcat(elements(element_selection),'_plot',num2str(current_time(4)),num2str(current_time(5)));
saveas(gcf,file_save_image{1},'jpg');
set(handles.status,'String','Data Ready.');
delete(gcf);




function reference_Callback(hObject, eventdata, handles)
% hObject    handle to reference (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of reference as text
%        str2double(get(hObject,'String')) returns contents of reference as a double


% --- Executes during object creation, after setting all properties.
function reference_CreateFcn(hObject, eventdata, handles)
% hObject    handle to reference (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in plot_type.
function plot_type_Callback(hObject, eventdata, handles)
% hObject    handle to plot_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns plot_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plot_type


% --- Executes during object creation, after setting all properties.
function plot_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


