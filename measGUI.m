function varargout = measGUI(varargin)
% MEASGUI MATLAB code for measGUI.fig
%      MEASGUI, by itself, creates a new MEASGUI or raises the existing
%      singleton*.
%
%      H = MEASGUI returns the handle to a new MEASGUI or the handle to
%      the existing singleton*.
%
%      MEASGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MEASGUI.M with the given input arguments.
%
%      MEASGUI('Property','Value',...) creates a new MEASGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before measGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to measGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help measGUI

% Last Modified by GUIDE v2.5 31-Aug-2015 17:08:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @measGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @measGUI_OutputFcn, ...
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


% --- Executes just before measGUI is made visible.
function measGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to measGUI (see VARARGIN)

% Choose default command line output for measGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% default date is current date
today = datestr(date, 29);
dashes = strfind(today, '-');
currYear = today(3:4);
currMonth = num2str(str2num(today(dashes(1)+1:dashes(2)-1)));
currDay = num2str(str2num(today(dashes(2)+1:end)));

setappdata(handles.GUI, 'currDate', [currMonth '-' currDay '-' currYear]);
setappdata(handles.GUI, 'dateU', [currMonth '_' currDay '_' currYear]);
setappdata(handles.GUI, 'type', 'Mouse Embryo');
setappdata(handles.GUI, 'embryoNum', '1');
setappdata(handles.GUI, 'pipSize', 128);
setappdata(handles.GUI, 'manualCorner', 0);
setappdata(handles.GUI, 'manualMeasure', 0);
setappdata(handles.GUI, 'filePathRaw', '');
setappdata(handles.GUI, 'filePathProc', '');
setappdata(handles.GUI, 'videoLoaded', 0);
setappdata(handles.GUI, 'frames', []);
setappdata(handles.GUI, 't', []);
setappdata(handles.GUI, 'params', []);
setappdata(handles.GUI, 'displayFig', NaN);

% UIWAIT makes measGUI wait for user response (see UIRESUME)
% uiwait(handles.GUI);


% --- Outputs from this function are returned to the command line.
function varargout = measGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in TypeSelector.
function TypeSelector_Callback(hObject, eventdata, handles)
% hObject    handle to TypeSelector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns TypeSelector contents as cell array
%        contents{get(hObject,'Value')} returns selected item from TypeSelector

contents = cellstr(get(hObject,'String'));
type = contents{get(hObject,'Value')};
setappdata(handles.GUI, 'type', type);



% --- Executes during object creation, after setting all properties.
function TypeSelector_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TypeSelector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function DateTextEdit_Callback(hObject, eventdata, handles)
% hObject    handle to DateTextEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DateTextEdit as text
%        str2double(get(hObject,'String')) returns contents of DateTextEdit as a double

% update currDate and dateU
currDate = get(hObject, 'String');
dateU = currDate;
dateUI = strfind(currDate, '-');
dateU(dateUI) = '_';
setappdata(handles.GUI, 'currDate', currDate);
setappdata(handles.GUI, 'dateU', dateU);


% --- Executes during object creation, after setting all properties.
function DateTextEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DateTextEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PipSizeEdit_Callback(hObject, eventdata, handles)
% hObject    handle to PipSizeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PipSizeEdit as text
%        str2double(get(hObject,'String')) returns contents of PipSizeEdit as a double
pipSize = str2double(get(hObject,'String'))
setappdata(handles.GUI, 'pipSize', pipSize);


% --- Executes during object creation, after setting all properties.
function PipSizeEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PipSizeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in MeasPipBtn.
function MeasPipBtn_Callback(hObject, eventdata, handles)
% hObject    handle to MeasPipBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% display a frame in the video, tell user to click top and bottom of
% pipette

displayFig = getappdata(handles.GUI, 'displayFig');

if ~ishandle(displayFig)
    
    if isnan(displayFig)
        fig = figure(1);
        setappdata(handles.GUI, 'displayFig', fig);
    else
        figure(displayFig);
    end
    
    clf;
    frames = getappdata(handles.GUI, 'frames');
    imshow(frames(:,:,1));
    
else
    figure(displayFig);
end

[x y] = ginput(2)
pipSize = round(abs(y(2) - y(1)));
setappdata(handles.GUI, 'pipSize', pipSize);
set(handles.PipSizeEdit, 'String', num2str(pipSize));


function EmbryoNameEdit_Callback(hObject, eventdata, handles)
% hObject    handle to EmbryoNameEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EmbryoNameEdit as text
%        str2double(get(hObject,'String')) returns contents of EmbryoNameEdit as a double
embryoNum = get(hObject, 'String');
setappdata(handles.GUI, 'embryoNum', embryoNum);


% --- Executes during object creation, after setting all properties.
function EmbryoNameEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EmbryoNameEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ManualMeasure.
function ManualMeasure_Callback(hObject, eventdata, handles)
% hObject    handle to ManualMeasure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ManualMeasure
if get(hObject,'Value')
    manualMeasure = 1;
else
    manualMeasure = 0;
end

setappdata(handles.GUI, 'manualMeasure', manualMeasure);



% --- Executes on button press in ManualCorner.
function ManualCorner_Callback(hObject, eventdata, handles)
% hObject    handle to ManualCorner (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ManualCorner
if get(hObject,'Value')
    manualCorner = 1;
else
    manualCorner = 0;
end

setappdata(handles.GUI, 'manualCorner', manualCorner);


% --- Executes on button press in EmbryoMeasBtn.
function EmbryoMeasBtn_Callback(hObject, eventdata, handles)
% hObject    handle to EmbryoMeasBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in LoadVideoBtn.
function LoadVideoBtn_Callback(hObject, eventdata, handles)
% hObject    handle to LoadVideoBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currDate = getappdata(handles.GUI, 'currDate')
dateU = getappdata(handles.GUI, 'dateU')
type = getappdata(handles.GUI, 'type')
embryoNum = getappdata(handles.GUI, 'embryoNum')
pipSize = getappdata(handles.GUI, 'pipSize')
manualCorner = getappdata(handles.GUI, 'manualCorner')
manualMeasure = getappdata(handles.GUI, 'manualMeasure')
filePathRaw = getappdata(handles.GUI, 'filePathRaw')
filePathProc = getappdata(handles.GUI, 'filePathProc')


[frames, t, params] = LoadVideo(type, currDate, pipSize, filePathRaw, ...
    embryoNum, handles);
setappdata(handles.GUI, 'frames', frames);
setappdata(handles.GUI, 't', t);
setappdata(handles.GUI, 'params', params);

setappdata(handles.GUI, 'videoLoaded', 1);
set(handles.MeasPipBtn, 'Enable', 'on');
set(handles.EmbryoMeasBtn, 'Enable', 'on');
set(handles.NewDisplayBtn, 'Enable', 'on');


% --- Executes on button press in RawDataBtn.
function RawDataBtn_Callback(hObject, eventdata, handles)
% hObject    handle to RawDataBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = cellstr(get(handles.TypeSelector, 'String'));
type = contents{get(handles.TypeSelector,'Value')};
filePathRaw = uigetdir(['C:\Users\Livia\Desktop\IVF\Raw Data\Videos\' type], ...
    'Select Raw Data Folder');
setappdata(handles.GUI, 'filePathRaw', filePathRaw);
set(handles.RawDataLabel, 'String', filePathRaw);


% --- Executes on button press in ProcDataBtn.
function ProcDataBtn_Callback(hObject, eventdata, handles)
% hObject    handle to ProcDataBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = cellstr(get(handles.TypeSelector, 'String'));
type = contents{get(handles.TypeSelector,'Value')};
filePathProc = uigetdir(['C:\Users\Livia\Desktop\IVF\Processed Data\' type], ...
    'Select Processed Data Folder');
setappdata(handles.GUI, 'filePathProc', filePathProc);
set(handles.ProcDataLabel, 'String', filePathProc);


% --- Executes during object creation, after setting all properties.
function GUI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in NewDisplayBtn.
function NewDisplayBtn_Callback(hObject, eventdata, handles)
% hObject    handle to NewDisplayBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


frames = getappdata(handles.GUI, 'frames');
displayFig = getappdata(handles.GUI, 'displayFig');

if ~isempty(frames)
    
    if isnan(displayFig)
        fig = figure(1);
        setappdata(handles.GUI, 'displayFig', fig);
    else
        figure(displayFig);
    end
    
    clf;
    imshow(frames(:,:,1));
else
    errordlg('Error: no frames loaded');
end
