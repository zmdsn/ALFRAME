function varargout = figscan(varargin)
% FIGSCAN MATLAB code for figscan.fig
%      FIGSCAN, by itself, creates a new FIGSCAN or raises the existing
%      singleton*.
%
%      H = FIGSCAN returns the handle to a new FIGSCAN or the handle to
%      the existing singleton*.
%
%      FIGSCAN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FIGSCAN.M with the given input arguments.
%
%      FIGSCAN('Property','Value',...) creates a new FIGSCAN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before figscan_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to figscan_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help figscan

% Last Modified by GUIDE v2.5 24-Oct-2016 17:29:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @figscan_OpeningFcn, ...
                   'gui_OutputFcn',  @figscan_OutputFcn, ...
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


% --- Executes just before figscan is made visible.
function figscan_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to figscan (see VARARGIN)

% Choose default command line output for figscan
handles.output = hObject;
handles.dirpath = './figure';
handles.h = 0;
handles.value = 0;
handles.filename = '';

set(handles.foldname,'String',handles.dirpath);
handles.filenames = get_fig_files(handles.dirpath);
set(handles.listbox1,'String',handles.filenames);


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes figscan wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = figscan_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in selfold.
function selfold_Callback(hObject, eventdata, handles)
% hObject    handle to selfold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.dirpath = uigetdir();
if isequal(handles.dirpath,0)
    disp('Users Selected Canceled');
else
    set(handles.foldname,'String',handles.dirpath);
    handles.filenames = get_fig_files(handles.dirpath);
    set(handles.listbox1,'String',handles.filenames);
end
guidata(hObject, handles);

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    if (handles.value == 0)
        handles.filename = get(handles.listbox1,'string');
        handles.value = get(handles.listbox1,'value');
    end
    nlang = length(handles.filename);
    if handles.value < nlang
        handles.value = handles.value + 1;
        str = char(strcat(handles.dirpath ,'\', handles.filename(handles.value)));
        set(handles.listbox1,'value',handles.value);
        if handles.h ~= 0 
            close(handles.h);
        end
        handles.h = open(str);
        guidata(hObject, handles);
    else
        msgbox('the end');
    end
    
% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    if (handles.value == 0)
        handles.filename = get(handles.listbox1,'string');
        handles.value = get(handles.listbox1,'value');
    end
    nlang = length(handles.filename);
    if handles.value > 1
        handles.value = handles.value - 1;
        str = char(strcat(handles.dirpath ,'\', handles.filename(handles.value)));
        set(handles.listbox1,'value',handles.value);
        if handles.h ~= 0 
            close(handles.h);
        end
        handles.h = open(str);
        guidata(hObject, handles);
    else
        msgbox('the end');
    end
    

% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1
    handles.filename = get(handles.listbox1,'string');
    handles.value = get(handles.listbox1,'value');
    str = char(strcat(handles.dirpath ,'\', handles.filename(handles.value)));
    if handles.h ~= 0 
        close(handles.h);
    end
    handles.h = open(str);
    guidata(hObject, handles);
    



% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function [ fileNames ] = get_fig_files( f_road )
% 获取文件夹的路径    f_road
    fileFolder = fullfile( f_road );
    dirOutput = dir(fullfile(fileFolder,'*.fig'));
    fileNames = { dirOutput.name };
    for j = 1:length(fileNames)
        str = fileNames(j);
        str{1} = str{1}(1:end);
        fileNames(j) =  [str];
    end
  
  


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    filenames = get_fig_files(handles.dirpath);
    n = length(filenames);
    con_fold = char(strcat(handles.dirpath ,'/img'));
    if ~exist('con_fold','dir')
        mkdir(con_fold);
    end
    for i = 1:n
        fig = char(strcat(handles.dirpath ,'/', filenames(i)));
        name = char(strcat(handles.dirpath ,'/img/',filenames{i}(1:end-3),'jpg'))
        h = open(fig);
        saveas(gcf,name);
        close(h);
    end
    msgbox('success ! ');
