function varargout = gui2(varargin)
% GUI2 MATLAB code for gui2.fig
%      GUI2, by itself, creates a new GUI2 or raises the existing
%      singleton*.
%
%      H = GUI2 returns the handle to a new GUI2 or the handle to
%      the existing singleton*.
%
%      GUI2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI2.M with the given input arguments.
%
%      GUI2('Property','Value',...) creates a new GUI2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui2

% Last Modified by GUIDE v2.5 15-Jul-2019 15:32:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui2_OpeningFcn, ...
                   'gui_OutputFcn',  @gui2_OutputFcn, ...
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


% --- Executes just before gui2 is made visible.
function gui2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui2 (see VARARGIN)

% Choose default command line output for gui2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);
clearvars -except hObject eventdata handles
update(hObject, eventdata, handles);


% --- Outputs from this function are returned to the command line.
function varargout = gui2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%%%%%%%%%%
% update %
%%%%%%%%%%
function update(hObject, eventdata, handles)
global y band rect lpf_out map_out active_partial mod_out out;
global fs t channels nsamples a b;
global pulses_per_sec_per_channel start_freq log_gap_freq; %for compressed analogue 
in_filename = get(handles.inputfile_name, 'String');
[x,fs] = audioread(in_filename); %p4.mp3
y = x(:, 1);
nsamples = size(y);
filter_order = 2;
channels = str2num(get(handles.channels, 'String'));
ts = 1/fs;
t = 0 : ts : (nsamples-1)*ts; %xaxis showing seconds
plot(handles.input_plot, t, y);
title(handles.input_plot, 'Input Signal'); xlabel(handles.input_plot, 'time'); ylabel(handles.input_plot, 'amplitude'); 

%% BPF
start_freq = str2num(get(handles.cutoff1, 'String'));
end_freq = str2num(get(handles.cutoff2, 'String'));
log_gap_freq = ( log10(end_freq) - log10(start_freq) ) / channels; %(start freq + end freq)/channels
logf_prev = log10(start_freq);
logf_next = log10(start_freq) + log_gap_freq;

band = zeros(nsamples(1), channels);
for c = 1:channels
    w_prev = ( (10^logf_prev) / (fs/2) );
    w_next = ( (10^logf_next) / (fs/2) );
    Wn1 = [w_prev, w_next];
    [b(c,:), a(c,:)] = butter(filter_order, Wn1, 'bandpass');
    band(:, c) = filter(b(c,:), a(c,:), y);
    logf_prev = logf_next;
    logf_next = logf_next + log_gap_freq;
end

%% RECTIFICATION
rect = abs(band);
 
%% LPF - ENVELOPE
Wn_l = (200/(fs/2));
[b_l, a_l] = butter(2, Wn_l, 'low');
lpf_out = zeros(nsamples(1), channels);
for c = 1:channels
    lpf_out(:, c) = filter(b_l, a_l, rect(:, c) );
end

%% NON-LINEAR MAPPING. DOWNSAMPLE
map_out = zeros(nsamples(1), channels);
amp_nlevels = str2num(get(handles.amplitudes, 'String'));
amp_levels = zeros(amp_nlevels,1); % thresholds
amp_levels_out = zeros(amp_nlevels,1); % levels output of the mappging
min_amp = str2num(get(handles.min_amp, 'String'));
max_amp = str2num(get(handles.max_amp, 'String'));

amp_prev = 0;
logamp_next = log10(min_amp);
log_gap_amp = ( log10(max_amp) - log10(min_amp) ) / (amp_nlevels-1);
for i = 1:amp_nlevels
    amp_levels(i,1) = 10^logamp_next;
    amp_levels_out(i,1) = (amp_prev+10^logamp_next)/2;
    amp_prev = 10^logamp_next;
    logamp_next = logamp_next + log_gap_amp;
end

%DO BINARY SEARCH!!
for c = 1:channels
    for r = 1:nsamples
        i=1;
        while lpf_out(r,c)>amp_levels(i,1)
            i=i+1;
        end
        map_out(r,c) = amp_levels_out(i,1);
    end
end

%% ONE ACTIVE CHANNEL + TIME SPREAD
%one active channel
pulses_per_sec_per_channel = str2num(get(handles.pps_channel, 'String'));
pulses_per_sec = pulses_per_sec_per_channel * channels;
t_pulse = 1 / pulses_per_sec;
n = round( t_pulse / ts );
active_partial = zeros(nsamples(1), channels);
for r = 1:nsamples
    active_ch = mod( (floor((r-1)/n)), channels) + 1;
    active_partial(r, active_ch) = map_out(r, active_ch);
end

%time spread (lpf)
if get(handles.time_spread, 'Value')
    timecutoff = str2num(get(handles.time_cutoff, 'String'));
    Wn_l = (timecutoff/(fs/2));
    [b_l, a_l] = butter(2, Wn_l, 'low');
    for c = 1:channels
        active_partial(:, c) = filter(b_l, a_l, active_partial(:, c) );
    end
end

%% MODULATION
mod_out = zeros(nsamples(1), channels);
f1 = start_freq;
logf2 = log10(f1) + log_gap_freq;
f2 = (10^logf2);

if get(handles.frequency_spread, 'Value') % -- spread around A --
    n_spread = str2num(get(handles.n_freqspread, 'String'));
    spread = zeros(nsamples(1), channels);
    log_gap_fspread = log_gap_freq/(2*n_spread); % OR /(2*n_spread*4) comprime lo spread around less frequencies.
end % -- end spread A --

for c = 1:channels
    sinefreq = (f1 + f2)/2;
    s(:,1) = sin(2 * pi * sinefreq * t);
    mod_out(:, c) = s(:,1) .* active_partial(:, c);
    
    if get(handles.frequency_spread, 'Value') % -- spread around B --
        fsin_log = log10(sinefreq);
        for i=(-n_spread):(-1) %towards -ve frequencies
            spread_amp = (n_spread+i)/n_spread;
            fspread_log = fsin_log + i*log_gap_fspread;
            fspread = 10^fspread_log;
            spread_carrier(:,1) = rot90( spread_amp*sin(2 * pi * fspread * t) );
            spread(:,c) = spread(:,c) + (active_partial(:,c) .* spread_carrier(:,1));
        end
        for i=1:n_spread %towards +ve frequencies
            spread_amp = (n_spread-i)/n_spread;
            fspread_log = fsin_log + i*log_gap_fspread;
            fspread = 10^fspread_log;
            spread_carrier(:,1) = rot90( spread_amp*sin(2 * pi * fspread * t) );
            spread(:,c) = spread(:,c) + (active_partial(:,c) .* spread_carrier(:,1));
        end
    end     % -- end spread B --
    
    f1 = f2;
    logf2 = logf2 + log_gap_freq;
    f2 = (10^logf2);
end

%% OUT 
% add freq spread
if get(handles.frequency_spread, 'Value')
    for c = 1:channels
        mod_out(:, c) = (2/n_spread)*(mod_out(:, c) + spread(:, c));
    end
end
% sum them up
out = zeros(nsamples(1), 1);
for c = 1:channels
    out(:,1) = out(:,1) + mod_out(:,c);
end
% output gain
if get(handles.out_gain, 'Value')
    out = out * str2double(get(handles.out_gainvalue, 'String'));
end
% plot
plot(handles.output_plot, t, out);
title(handles.output_plot, 'CIS Output Signal'); xlabel(handles.output_plot, 'time'); ylabel(handles.output_plot, 'amplitude'); zoom on;


%%%%%%%%%%%%
% gui cose %
%%%%%%%%%%%%
function inputfile_name_Callback(hObject, eventdata, handles)
% hObject    handle to inputfile_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of inputfile_name as text
%        str2double(get(hObject,'String')) returns contents of inputfile_name as a double
clearvars -except hObject eventdata handles
update(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function inputfile_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inputfile_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function outputfile_name_Callback(hObject, eventdata, handles)
% hObject    handle to outputfile_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of outputfile_name as text
%        str2double(get(hObject,'String')) returns contents of outputfile_name as a double

% --- Executes during object creation, after setting all properties.


function outputfile_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to outputfile_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save_out.
function save_out_Callback(hObject, eventdata, handles)
% hObject    handle to save_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global out fs;
out_filename = get(handles.outputfile_name, 'String');
audiowrite(out_filename,out,fs);


function channels_Callback(hObject, eventdata, handles)
% hObject    handle to channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of channels as text
%        str2double(get(hObject,'String')) returns contents of channels as a double
input = str2num(get(hObject, 'String'));
if isempty(input)
    set (hObject, 'String','8');
end
update(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function channels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function cutoff1_Callback(hObject, eventdata, handles)
% hObject    handle to cutoff1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cutoff1 as text
%        str2double(get(hObject,'String')) returns contents of cutoff1 as a double
input = str2num(get(hObject, 'String'));
if isempty(input)
    set (hObject, 'String', '300');
end
update(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function cutoff1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cutoff1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function cutoff2_Callback(hObject, eventdata, handles)
% hObject    handle to cutoff2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cutoff2 as text
%        str2double(get(hObject,'String')) returns contents of cutoff2 as a double
input = str2num(get(hObject, 'String'));
if isempty(input) %or smaller than the 1st cutoff -- !! --
    set (hObject, 'String', '3500');
end
update(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function cutoff2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cutoff2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function amplitudes_Callback(hObject, eventdata, handles)
% hObject    handle to amplitudes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of amplitudes as text
%        str2double(get(hObject,'String')) returns contents of amplitudes as a double
input = str2num(get(hObject, 'String'));
if isempty(input)
    set (hObject, 'String', '16');
end
update(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function amplitudes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to amplitudes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function min_amp_Callback(hObject, eventdata, handles)
% hObject    handle to min_amp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of min_amp as text
%        str2double(get(hObject,'String')) returns contents of min_amp as a double
input = str2num(get(hObject, 'String'));
if isempty(input) %or smaller than the 1st cutoff -- !! --
    set (hObject, 'String', '0.0001');
end
update(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function min_amp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to min_amp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function max_amp_Callback(hObject, eventdata, handles)
% hObject    handle to max_amp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of max_amp as text
%        str2double(get(hObject,'String')) returns contents of max_amp as a double
input = str2num(get(hObject, 'String'));
if isempty(input) %or smaller than the 1st cutoff -- !! --
    set (hObject, 'String', '1');
end
update(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function max_amp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to max_amp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pps_channel_Callback(hObject, eventdata, handles)
% hObject    handle to pps_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pps_channel as text
%        str2double(get(hObject,'String')) returns contents of pps_channel as a double
input = str2num(get(hObject, 'String'));
if isempty(input) %or smaller than the 1st cutoff -- !! --
    set (hObject, 'String', '1800');
end
update(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function pps_channel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pps_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in time_spread.
function time_spread_Callback(hObject, eventdata, handles)
% hObject    handle to time_spread (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of time_spread
update(hObject, eventdata, handles);


function time_cutoff_Callback(hObject, eventdata, handles)
% hObject    handle to time_cutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of time_cutoff as text
%        str2double(get(hObject,'String')) returns contents of time_cutoff as a double
input = str2num(get(hObject, 'String'));
if isempty(input) %or smaller than the 1st cutoff -- !! --
    set (hObject, 'String', '5000');
end
if get(handles.time_spread, 'Value')
    update(hObject, eventdata, handles);
end

% --- Executes during object creation, after setting all properties.
function time_cutoff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_cutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in frequency_spread.
function frequency_spread_Callback(hObject, eventdata, handles)
% hObject    handle to frequency_spread (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of frequency_spread
update(hObject, eventdata, handles);


function n_freqspread_Callback(hObject, eventdata, handles)
% hObject    handle to n_freqspread (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of n_freqspread as text
%        str2double(get(hObject,'String')) returns contents of n_freqspread as a double
input = str2num(get(hObject, 'String'));
if isempty(input) %or smaller than the 1st cutoff -- !! --
    set (hObject, 'String', '4');
end
if get(handles.frequency_spread, 'Value')
    update(hObject, eventdata, handles);
end


% --- Executes during object creation, after setting all properties.
function n_freqspread_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n_freqspread (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in out_gain.
function out_gain_Callback(hObject, eventdata, handles)
% hObject    handle to out_gain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of out_gain
update(hObject, eventdata, handles);


function out_gainvalue_Callback(hObject, eventdata, handles)
% hObject    handle to out_gainvalue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of out_gainvalue as text
%        str2double(get(hObject,'String')) returns contents of out_gainvalue as a double
input = str2double(get(hObject, 'String'));
if isempty(input) %or smaller than the 1st cutoff -- !! --
    set (hObject, 'String', '1');
end
if get(handles.out_gain, 'Value')
    update(hObject, eventdata, handles);
end


% --- Executes during object creation, after setting all properties.
function out_gainvalue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to out_gainvalue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in play_input.
function play_input_Callback(hObject, eventdata, handles)
% hObject    handle to play_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global y fs;
if exist('y')
    sound(y, fs);
end


% --- Executes on button press in play_output.
function play_output_Callback(hObject, eventdata, handles)
% hObject    handle to play_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global out fs;
if exist('out')
    sound(out, fs);
end


% --- Executes on button press in plot1.
function plot1_Callback(hObject, eventdata, handles)
% hObject    handle to plot1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global t y band channels;
if exist('y') && exist('band')
    figure;
    subplot(channels+1, 1, 1);
    plot(t, y, '-r');
    title ('original signal + band pass filtered');
    for c = 1:channels
        subplot(channels+1, 1, c+1);
        plot(t, band(:,c));
    end
end


% --- Executes on button press in plot2.
function plot2_Callback(hObject, eventdata, handles)
% hObject    handle to plot2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global t rect channels;
figure;
for c = 1:channels
    subplot(channels, 1, c);
    plot(t, rect(:,c), '-b');
end


% --- Executes on button press in plot3.
function plot3_Callback(hObject, eventdata, handles)
% hObject    handle to plot3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global t lpf_out channels
figure;
for c = 1:channels
    subplot(channels, 1, c);
    plot(t, lpf_out(:,c), '-r');
end


% --- Executes on button press in plot4.
function plot4_Callback(hObject, eventdata, handles)
% hObject    handle to plot4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global t map_out channels
figure;
for c = 1:channels
    subplot(channels, 1, c);
    plot(t, map_out(:,c), '-k');
end


% --- Executes on button press in plot5.
function plot5_Callback(hObject, eventdata, handles)
% hObject    handle to plot5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global t active_partial channels
figure;
for c = 1:channels
    subplot(channels, 1, c);
    plot(t, active_partial(:,c), '-b');
end


% --- Executes on button press in plot6.
function plot6_Callback(hObject, eventdata, handles)
% hObject    handle to plot6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global t mod_out channels
figure;
for c = 1:channels
    subplot(channels, 1, c);
    plot(t, mod_out(:,c), '-m');
end


% --- Executes on button press in listen1.
function listen1_Callback(hObject, eventdata, handles)
% hObject    handle to listen1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global band channels nsamples fs;
if exist('band')
    band_sound = zeros(nsamples(1), 1);
    for c=1:channels
        band_sound(:,1) = band_sound(:,1) + band(:,c);
    end
    sound(band_sound, fs);
end


% --- Executes on button press in listen2.
function listen2_Callback(hObject, eventdata, handles)
% hObject    handle to listen2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global rect channels nsamples fs;
if exist('rect')
    rect_sound = zeros(nsamples(1), 1);
    for c=1:channels
        rect_sound(:,1) = rect_sound(:,1) + rect(:,c);
    end
    sound(rect_sound, fs);
end


% --- Executes on button press in listen3.
function listen3_Callback(hObject, eventdata, handles)
% hObject    handle to listen3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global lpf_out channels nsamples fs;
if exist('lpf_out')
    lpf_sound = zeros(nsamples(1), 1);
    for c=1:channels
        lpf_sound(:,1) = lpf_sound(:,1) + lpf_out(:,c);
    end
    sound(lpf_sound, fs);
end


% --- Executes on button press in listen4.
function listen4_Callback(hObject, eventdata, handles)
% hObject    handle to listen4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global map_out channels nsamples fs;
if exist('map_out')
    map_sound = zeros(nsamples(1), 1);
    for c=1:channels
        map_sound(:,1) = map_sound(:,1) + map_out(:,c);
    end
    sound(map_sound, fs);
end


% --- Executes on button press in listen5.
function listen5_Callback(hObject, eventdata, handles)
% hObject    handle to listen5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global active_partial channels nsamples fs;
if exist('active_partial')
    active_sound = zeros(nsamples(1), 1);
    for c=1:channels
        active_sound(:,1) = active_sound(:,1) + active_partial(:,c);
    end
    sound(active_sound, fs);
end


% --- Executes on button press in listen6.
function listen6_Callback(hObject, eventdata, handles)
% hObject    handle to listen6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global mod_out channels nsamples fs;
if exist('mod_out')
    mod_sound = zeros(nsamples(1), 1);
    for c=1:channels
        mod_sound(:,1) = mod_sound(:,1) + mod_out(:,c);
    end
    sound(mod_sound, fs);
end


% --- Executes on button press in plot_filters.
function plot_filters_Callback(hObject, eventdata, handles)
% hObject    handle to plot_filters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global a b;
if exist('a') && exist('b')
    fvtool(b(1,:),a(1,:),b(2,:),a(2,:),b(3,:),a(3,:),b(4,:),a(4,:),b(5,:),a(5,:),b(6,:),a(6,:),b(7,:),a(7,:),b(8,:),a(8,:))
end


% --- Executes on button press in plotspec.
function plotspec_Callback(hObject, eventdata, handles)
% hObject    handle to plotspec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global y map_out out fs;
figure;
subplot(3, 1, 1);
spectrogram(y,256,250,256,fs,'yaxis'); %parameters scelti un po' a caso aggiusta %vorrei solo 300-5000??
title ('original signal');
test = compressed_analogue(map_out);  
subplot(3, 1, 2);
spectrogram(test,256,250,256,fs,'yaxis');
title ('output signal - all channels active');
subplot(3, 1, 3);
spectrogram(out,256,250,256,fs,'yaxis');
title ('output signal - one active channel at a time');


% --- Executes on button press in plotfft.
function plotfft_Callback(hObject, eventdata, handles)
% hObject    handle to plotfft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global y out map_out nsamples channels;

Fy = fft(y);
ampin = sqrt (Fy .* conj(Fy) );
figure
subplot(3,1,1)
plot(ampin)
title('Input'); grid on; xlabel('Frequency'); ylabel('Amplitude');

test = compressed_analogue(map_out);
Ftest = fft(test);
amptest = sqrt (Ftest .* conj(Ftest) );
subplot(3,1,2)
plot(amptest)
title('Output - all channels active'); grid on; xlabel('Frequency'); ylabel('Amplitude');

Fout = fft(out);
ampout = sqrt (Fout .* conj(Fout) );
subplot(3,1,3)
plot(ampout)
title('Output - one active channel at a time'); grid on; xlabel('Frequency'); ylabel('Amplitude');
%X AXIS non dovrebbe ripertersi dopo 22k?? 

function [ca_out] = compressed_analogue(map_out) % generate compressed analogue signal
global channels nsamples fs t start_freq log_gap_freq pulses_per_sec_per_channel;
ca_out(:,1) = zeros(nsamples(1), 1);
f1 = start_freq;
logf2 = log10(f1) + log_gap_freq;
f2 = (10^logf2);

pulses_per_sec = pulses_per_sec_per_channel * channels;
t_pulse = 1 / pulses_per_sec;
ts = 1/fs;
n = round( t_pulse / ts );
n_period = n*channels;

for c = 1:channels
    sinefreq = (f1 + f2)/2;
    sr = sin(2 * pi * sinefreq * t);
    s = rot90(sr);
    ca_out(:,1) = ca_out(:,1) + (s(:,1) .* map_out(:, c));
    
    f1 = f2;
    logf2 = logf2 + log_gap_freq;
    f2 = (10^logf2);
end


% --- Executes on button press in listen_ca.
function listen_ca_Callback(hObject, eventdata, handles)
% hObject    handle to listen_ca (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global map_out fs;
ca = compressed_analogue(map_out);
sound(ca, fs);