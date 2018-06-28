function varargout = sochro(varargin)
% Social Chronnectome Toolbox
% Developed by Chase Sherwell (2017)
%      The Social Chronnectome (sochro_v101) Toolbox was developed to analyse
%      EDA data using graph metrics. Analysis is optimised for data
%      collected using Empatica Wristbands.


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @sochro_v101_OpeningFcn, ...
    'gui_OutputFcn',  @sochro_v101_OutputFcn, ...
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


% --- Executes just before sochro_v101 is made visible.
function sochro_v101_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

% warning if below minimum screen size
handles.screen_size = get(0,'ScreenSize');
if handles.screen_size(3) < 1280 || handles.screen_size(4) < 800
    uiwait(errordlg('Minimum recommended screen resolution when using this program is 1280 x 800 pixels..','ERROR','modal'));
end

% POSSIBLE PLOTLINE COLOURS -----------------------------------------------
handles.plotline_colours=[];
for colours_list_repeat = 1:100
    handles.plotline_colours= [handles.plotline_colours;
        0 0  1; 0  1 0;  1 0 0; .5 .5 0; .5  0 .5; 0 .5 .5; .5  .5  .5; ...
        0 0 .5; 0 .5 0; .5 0 0;  1 .5 0;  1  0 .5; 0  1 .5;  1  .5  .5; ...
        .5  1 0; .5  0  1; 0 .5  1; .5   1  .5; ...
        1  1 0;  1  0  1; 0  1  1; .5  .5   1; ...
        1  .75 .25; ...
        .75 .25  1; ...
        .25  1  .75;];
end

%define subj selection window parameters
handles.uipanel_left_y_offset = get(handles.vis_uipanel,'Position');
handles.uipanel_left_y_offset = handles.uipanel_left_y_offset(2);
select_subj_dropdown_position = get(handles.vissubj_toggle,'Position');
handles.select_subj_dropdown_y = select_subj_dropdown_position(2)+handles.uipanel_left_y_offset;
handles.hidden_subjs_listbox_position = get(handles.vissubj_list,'Position');
handles.height_per_listbox_line = 13.6;

% %group uipanel selection window parameters
handles.uipanel_group_y_offset = get(handles.group_uipanel,'Position');
handles.uipanel_group_y_offset = handles.uipanel_group_y_offset(2);
select_group_dropdown_position = get(handles.groupsubj_toggle,'Position');
handles.select_group_dropdown_y = select_group_dropdown_position(2)+handles.uipanel_group_y_offset;
handles.hidden_group_listbox_position = get(handles.groupsubj_list,'Position');

% group selection dropdown window
select_groups_dropdown_position = get(handles.visgroups_toggle,'Position');
handles.select_groups_dropdown_y = select_groups_dropdown_position(2)+handles.uipanel_left_y_offset;
handles.hidden_groups_listbox_position = get(handles.visgroups_list,'Position');
handles.height_per_listbox_line = 13.6;

%% CREATE EMPTY VARIABLES
handles.sdir=0;
handles.data=[];
handles.parameters=[];
handles.time=[];
handles.display=[];
handles.EDA=[];
handles.ACC=[];
handles.groups=[];
handles.group_subj = [];
handles.conditions = [];
handles.surrogate = [];
handles.corr = [];
handles.pos = get(handles.edaaxes,'Position');

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = sochro_v101_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in loadraw_push.
function loadraw_push_Callback(hObject, eventdata, handles)
if handles.sdir~=0,
    buttonresponse = questdlg('Do you wish to clear current workspace to load new data?',...
        'Load raw data',...
        'Yes','No','Yes');
    switch buttonresponse,
        case 'Yes', handles.sdir=0; handles = clearvariables(handles);
        case 'No', return
    end
end

handles.sdir = uigetdir(cd,'Select directory containing biometric data...');
if isequal(handles.sdir,0), return; end %if nothing selected, stop function
clc
if ispc, fs='\'; else fs='/'; end %switch file separator for different OS
%% Get directory and ID names
subdir = dir(handles.sdir); %get file in directory
subdir = {subdir.name};
for ii = 1:length(subdir),
    ff = strncmp('.',subdir{ii},1);
    if ff==1, fn(ii)=0; else fn(ii)=ii; end
end
subj_id={subdir{fn>0}};
clear subdir fn

%% Get individual data
fn = []; %counter for removing subjects
for ii = 1:length(subj_id), %for each subject
    set(handles.status_edit,'String',['Loading data for subject ' subj_id{ii}]); pause(0.01);
    subdir = [handles.sdir fs subj_id{ii}];
    if isdir(subdir), %if listing is a directory
        fl = dir(subdir); %get file list of everything in that directory
        handles.data(ii).subj = subj_id{ii};
        EDAcheck = 0; ACCcheck = 0; %check variables for data
        for kk = 1:length(fl), %go through every file in folder fl
            if strcmp('EDA.csv',fl(kk).name), %if EDA file
                disp(['Loading subject: ' subj_id{ii}])
                handles.data(ii).urEDA=importdata([subdir fs 'EDA.csv']);
                if length(handles.data(ii).urEDA)>2, %if no data
                    EDAcheck = 1;
                end
            elseif strcmp('ACC.csv',fl(kk).name), %if ACC file
                handles.data(ii).urACC=importdata([subdir fs 'ACC.csv']);
                ACCcheck = 1;
            end
        end %end of file loop
        if EDAcheck==0,
            disp(['NO DATA FOUND FOR ' subj_id{ii} ' - removing from participant list']);
            fn(end+1)=ii;
        elseif ACCcheck==0,
            disp(['WARNING - NO ACCELEROMETER DATA FOUND FOR  ' subj_id{ii}]);
        end
    end
end
handles.data(fn) = []; subj_id(fn)=[]; %Remove missing datasets

%% Construct timelines
set(handles.status_edit,'String','Temporal alignment of datasets...'); pause(0.01);
for ii = 1:length(handles.data),
    %EDA
    handles.data(ii).EDAsamplerate=handles.data(ii).urEDA(2,1); %get samplerate
    handles.data(ii).EDAon=handles.data(ii).urEDA(1,1); %EDA time on
    handles.data(ii).EDAtimestamps=handles.data(ii).EDAon:1/handles.data(ii).EDAsamplerate: handles.data(ii).EDAon + (length(handles.data(ii).urEDA(3:end,1))/handles.data(ii).EDAsamplerate)-(1/handles.data(ii).EDAsamplerate);
    handles.data(ii).EDAsecs=handles.data(ii).EDAtimestamps-handles.data(ii).EDAon;
    handles.data(ii).EDAmins=handles.data(ii).EDAsecs./60;
    handles.data(ii).EDAoff=handles.data(ii).EDAtimestamps(end);
    
    %ACC
    if ~isempty(handles.data(ii).urACC), %if there's ACC data
        handles.data(ii).ACCsamplerate=handles.data(ii).urACC(2,1); %get samplerate
        handles.data(ii).ACCon=handles.data(ii).urACC(1,1); %EDA time on
        if handles.data(ii).EDAon~=handles.data(ii).ACCon,
            errordlg(['Warning - EDA and accelerometer data for subject ' handles.data(ii).subj ' have non-matching timestamps. Errors in computation may result. Please check that data files are correct']);
        end
        handles.data(ii).ACCtimestamps=handles.data(ii).ACCon:1/handles.data(ii).ACCsamplerate: handles.data(ii).ACCon + (length(handles.data(ii).urACC(3:end,1))/handles.data(ii).ACCsamplerate)-(1/handles.data(ii).ACCsamplerate);
        handles.data(ii).ACCsecs=handles.data(ii).ACCtimestamps-handles.data(ii).ACCon;
        handles.data(ii).ACCmins=handles.data(ii).ACCsecs./60;
        handles.data(ii).ACCoff=handles.data(ii).ACCtimestamps(end);
    end
    
    %GROUP DEFINITION (empty)
    handles.data(ii).group = [];
end
%Check samplerate consistency
if std([handles.data.EDAsamplerate])~=0,
    errordlg('Warning - Inconsistent sample rates between EDA files. Errors in computation will result. Please check data files');
end
if std([handles.data.ACCsamplerate])~=0,
    errordlg('Warning - Inconsistent sample rates between ACC files. Errors in computation will result. Please check data files');
end
%Copy samplerate to handles
handles.parameters.EDAsamplerate=handles.data(1).EDAsamplerate;
handles.parameters.ACCsamplerate=mean([handles.data.ACCsamplerate]);
%Construct master timeline
EDAtime=(min([handles.data.EDAon]):(1/handles.parameters.EDAsamplerate):max([handles.data.EDAoff]))';
ACCtime=(min([handles.data.ACCon]):(1/handles.parameters.ACCsamplerate):max([handles.data.ACCoff]))';
%Temporal alignment / pad with NaN
for ii = 1:length(handles.data),
    [inx, ~] = find(EDAtime==handles.data(ii).EDAon);
    handles.data(ii).EDA(1:length(EDAtime),1)=NaN;
    handles.data(ii).EDA(inx:inx+length(handles.data(ii).urEDA(3:end-1)),1)=handles.data(ii).urEDA(3:end);
    if ~isempty(handles.data(ii).urACC),
        [inx, ~] = find(ACCtime==handles.data(ii).ACCon);
        handles.data(ii).ACC(1:length(ACCtime),1:3)=NaN;
        handles.data(ii).ACC(inx:inx+length(handles.data(ii).urACC(3:end-1,1)),:)=handles.data(ii).urACC(3:end,:);
    end
end

handles.time.EDAtimestamps=EDAtime;
handles.time.EDAsecs=EDAtime-EDAtime(1);
handles.time.EDAmins=handles.time.EDAsecs./60;
handles.time.ACCtimestamps=ACCtime;
handles.time.ACCsecs=ACCtime-ACCtime(1);
handles.time.ACCmins=handles.time.ACCsecs./60;

%Compute summed accelerometer signals
for ii=1:length(handles.data),
    if ~isempty(handles.data(ii).urACC),
        handles.data(ii).ACC=sum((handles.data(ii).ACC).^2,2);
    end
end

%Copy original data to current data handles
handles.EDA=[handles.data.EDA];
handles.ACC=[handles.data.ACC];

%Set default times
set(handles.vistimemin_edit,'String',num2str(handles.time.EDAmins(1)));
set(handles.vistimemax_edit,'String',num2str(handles.time.EDAmins(end)));

%Populate Subj Lists
set(handles.vissubj_list,'String',{handles.data.subj});
set(handles.vissubj_list,'Value',1);
handles.group_subj = {handles.data.subj};
set(handles.groupsubj_list,'String',handles.group_subj); %update available subjects
set(handles.groupsubj_list,'Value',1); %update list selection
set(handles.status_edit,'String',['Dataset - ' handles.sdir ' ready for preprocessing']); pause(0.01);
set(handles.vis_uipanel,'Visible','on');
set(handles.preprocess_uipanel,'Visible','on');
handles.parameters.corr=0;
update_EDAdisplay(handles);
guidata(hObject, handles);

function [handles] = update_EDAdisplay(handles)
if handles.parameters.corr==0,
    axes(handles.edaaxes);
    %% DISPLAY EDA/ACC
    % Visibility of average & remove subject button
    if ~isempty(get(handles.vissubj_list,'Value'))
        if length(get(handles.vissubj_list,'Value')) > 1
            set(handles.visavg_check,'Visible','on');
            set(handles.visremovesubj_push,'Visible','off');
        else
            set(handles.visavg_check,'Visible','off');
            set(handles.visavg_check,'Value',0);
            set(handles.visremovesubj_push,'Visible','on');
        end
    else
        uiwait(errordlg('Please select subjects','ERROR','modal'));
        set(handles.edaaxes,'NextPlot','replace'); plot(handles.edaaxes,[0 0],[0 0]); return;
    end
    % Move appropriate data into display handle
    if get(handles.visEDA_radio,'Value')==1, %display EDA
        inx = get(handles.vissubj_list,'Value');
        handles.display.subj = {handles.data(inx).subj};
        if get(handles.visavg_check,'Value')==1,
            handles.display.data=nanmean(handles.EDA(:,inx),2);
        else handles.display.data=handles.EDA(:,inx);
        end
        handles.display.time=handles.time.EDAmins;
    elseif get(handles.visACC_radio,'Value')==1,
        inx = get(handles.vissubj_list,'Value');
        handles.display.subj = {handles.data(inx).subj};
        if get(handles.visavg_check,'Value')==1,
            handles.display.data=nanmean(handles.ACC(:,inx),2);
        else handles.display.data=handles.ACC(:,inx);
        end
        handles.display.time=handles.time.ACCmins;
    end
    
    % Set report list with selected subject information
    if ~isempty(get(handles.vissubj_list,'Value'))
        clear report
        if length(get(handles.vissubj_list,'Value')) > 1
            report{1}= [num2str(length(get(handles.vissubj_list,'Value'))) ' Subjects selected'];
            inx = get(handles.vissubj_list,'Value');
            firstdate=min([handles.data(inx).EDAon]);
            lastdate=max([handles.data(inx).EDAoff]);
            report{2}= ['Start: ' num2str(datestr(firstdate/(60*60*24)+datenum('01-Jan-1970')))];
            report{3}= ['End: ' num2str(datestr(lastdate/(60*60*24)+datenum('01-Jan-1970')))];
        else
            inx = get(handles.vissubj_list,'Value');
            report{1}= ['Subject ' handles.data(inx).subj ' selected'];
             ontime= find(~isnan(handles.display.data(:,1)),1); %get index of on time
             offtime= find(~isnan(handles.display.data(:,1)),1,'last'); %get index of off time
% %             ontime
% %             offtime
% %             length(handles.display.time)
  %          report{2}= ['Recording times: ' num2str(handles.display.time(1)) ' - ' num2str(handles.display.time(end)) ];
            report{2}= ['Recording times: ' num2str(handles.display.time(ontime)) ' - ' num2str(handles.display.time(offtime)) ];
            report{4}= ['Start: ' num2str(datestr(handles.data(inx).EDAtimestamps(1)/(60*60*24)+datenum('01-Jan-1970')))];
            report{5}= ['End: ' num2str(datestr(handles.data(inx).EDAtimestamps(end)/(60*60*24)+datenum('01-Jan-1970')))];
            
            datalength = handles.display.data;
            datalength = length(datalength(~isnan(datalength)));
            timelength = length(handles.display.time);
%             datalength
%             timelength
            report{3}= ['Data: ' num2str((datalength/timelength)*100) '% of session'];
        end
    else report = 'Status...';
    end
    set(handles.visreport_list,'String',report);
    % Check max amplitude
    if get(handles.visamp_check,'Value')==0,
        if str2double(get(handles.visampmin_edit,'String'))~=0 || str2double(get(handles.visampmax_edit,'String'))~=0,
            if str2double(get(handles.visampmin_edit,'String')) >= str2double(get(handles.visampmax_edit,'String')),
                uiwait(errordlg('Maximum amplitude must be larger than minimum amplitude','ERROR','modal'));
                set(handles.edaaxes,'NextPlot','replace'); plot(handles.edaaxes,[0 0],[0 0]); return;
            else handles.display.minamp = str2double(get(handles.visampmin_edit,'String'));
                handles.display.maxamp = str2double(get(handles.visampmax_edit,'String'));
            end
        else handles.display.minamp = min(min(handles.display.data));
            handles.display.maxamp = max(max(handles.display.data));
        end
    else handles.display.minamp = min(min([handles.display.data]));
        handles.display.maxamp = max(max([handles.display.data]));
    end
    % Check time range
    if get(handles.vistime_check,'Value')==0,
        if str2double(get(handles.vistimemin_edit,'String'))~=0 || str2double(get(handles.vistimemax_edit,'String'))~=0,
            if str2double(get(handles.vistimemin_edit,'String')) >= str2double(get(handles.vistimemax_edit,'String')),
                uiwait(errordlg('Maximum time range must be larger than minimum time range','ERROR','modal'));
                set(handles.edaaxes,'NextPlot','replace'); plot(handles.edaaxes,[0 0],[0 0]); return;
            elseif str2double(get(handles.vistimemin_edit,'String')) < min(handles.display.time),
                uiwait(errordlg(['Minimum display time must be greater than ' num2str(handles.display.time(1))],'ERROR','modal'));
                set(handles.edaaxes,'NextPlot','replace'); plot(handles.edaaxes,[0 0],[0 0]); return;
            elseif str2double(get(handles.vistimemax_edit,'String')) > max(handles.display.time)+1,
                uiwait(errordlg(['Maximum display time must be less than ' num2str(handles.display.time(end))],'ERROR','modal'));
                set(handles.edaaxes,'NextPlot','replace'); plot(handles.edaaxes,[0 0],[0 0]); return;
            else
                handles.display.mintime = str2double(get(handles.vistimemin_edit,'String'));
                handles.display.maxtime = str2double(get(handles.vistimemax_edit,'String'));
            end
        else handles.display.mintime = min(min(handles.display.time));
            handles.display.maxtime = max(max(handles.display.time));
        end
    else handles.display.mintime = min(min(handles.display.time));
        handles.display.maxtime = max(max(handles.display.time));
    end
    %Get data for correct times
    differences1 = handles.display.time - handles.display.mintime;
    smallest_difference = min(abs(differences1));
    lower_datapoint = find(abs(differences1) == smallest_difference);
    differences2 = handles.display.time - handles.display.maxtime;
    smallest_difference = min(abs(differences2));
    upper_datapoint = find(abs(differences2) == smallest_difference);
    handles.display.time = handles.display.time(lower_datapoint:upper_datapoint);
    handles.display.data = handles.display.data(lower_datapoint:upper_datapoint,:);
    %Check enough data is being plotted
    if handles.display.time(1) >= handles.display.time(end)
        uiwait(errordlg('Please select a wider data range to view','ERROR','modal'));
        set(handles.edaaxes,'NextPlot','replace'); plot(handles.edaaxes,[0 0],[0 0]); return;
    end
    %Plot lines
    set(handles.edaaxes,'NextPlot','replace');
    plot(handles.edaaxes,handles.display.time,handles.display.data(:,1), ...
        'LineStyle','-', ...
        'Color',handles.plotline_colours(1,:));
    set(handles.edaaxes,'NextPlot','add');
    for ii = 2:length(handles.display.data(1,:)),
        plot(handles.edaaxes,handles.display.time,handles.display.data(:,ii), ...
            'LineStyle','-', ...
            'Color',handles.plotline_colours(ii,:));
    end
    %Plot limits
    set(handles.edaaxes,'XLim',[handles.display.mintime handles.display.maxtime]);
    set(handles.edaaxes,'YLim',[handles.display.minamp handles.display.maxamp]);
    Ylimit = [handles.display.minamp handles.display.maxamp];
    %Axis labels
    set(get(handles.edaaxes,'XLabel'),'String','Time (minutes)');
    set(get(handles.edaaxes,'YLabel'),'String','Amplitude');
    % On/Off times
    if get(handles.visonoff_toggle,'Value')==1,
        for ii = 1:length(handles.display.data(1,:)), %for each participant dataset
            ontime= find(~isnan(handles.display.data(:,ii)),1); %get index of on time
            offtime= find(~isnan(handles.display.data(:,ii)),1,'last'); %get index of off time
            line( [handles.display.time(ontime) handles.display.time(ontime)], Ylimit, 'Color', [0 1 0]);
            line( [handles.display.time(offtime) handles.display.time(offtime)], Ylimit, 'Color', [1 0 0]);
        end
    end
    % Condition windows
    if get(handles.conddisp_toggle,'Value')==1
        if ~isempty(handles.conditions),
            selectedcondition = get(handles.cond_list,'Value');
            if ~isempty(selectedcondition),
                for ii = 1:length(handles.conditions(selectedcondition).windows),
                    t1 = handles.conditions(selectedcondition).windows{ii}{1}(1);
                    t2 = handles.conditions(selectedcondition).windows{ii}{1}(2);
                    
                    patch([t1 t2 t2 t1], [handles.display.minamp handles.display.minamp handles.display.maxamp handles.display.maxamp], 'red', 'FaceAlpha',.4, 'EdgeAlpha',.6);
                end
            end
        end
    end
else
    %% DISPLAY GRAPH METRICS
    %Get data type
    clear title
    title{2}=' ';
    kk = handles.parameters.condition-1; %remove 1 for title of poplist
    if kk>length(handles.conditions),
        title{3} = ' ';
    else title{3} = handles.conditions(kk).name;
    end
    if handles.parameters.connectivity==1, %STATIC connectivity
        title{4} = '[static] ';
        switch handles.parameters.level, %find specified level
            case 1, %Covariance
                plottype = 'covar';
                title{1} = 'Covariance ';
                switch handles.parameters.thresh,
                    case 1, handles.display.data = handles.corr.static(kk).R;
                        title{2} = ' ';
                    case 2, handles.display.data = handles.corr.static(kk).Rt;
                        title{2} = '(session threshold) ';
                end
                
            case 2, %Subject
                plottype = 'bar';
                switch handles.parameters.metric,
                    case 2, %Mean connectivity
                        title{1} = 'Mean connectivity ';
                        if kk > length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = vertcat(handles.corr.static.FC);
                                    title{2} = ' ';
                                case 2, handles.display.data = vertcat(handles.corr.static.FCt);
                                    title{2} = '(session threshold) ';
                            end
                        else
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = vertcat(handles.corr.static(kk).FC);
                                    title{2} = ' ';
                                case 2, handles.display.data = vertcat(handles.corr.static(kk).FCt);
                                    title{2} = '(session threshold) ';
                            end
                        end
                    case 3, %connectivity strength
                        title{1} = 'Connectivity strength';
                        if kk > length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = vertcat(handles.corr.static.S);
                                    title{2} = ' ';
                                case 2, handles.display.data = vertcat(handles.corr.static.St);
                                    title{2} = '(session threshold) ';
                            end
                        else
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = vertcat(handles.corr.static(kk).S);
                                    title{2} = ' ';
                                case 2, handles.display.data = vertcat(handles.corr.static(kk).St);
                                    title{2} = '(session threshold) ';
                            end
                        end
                    case 4, %Node degree
                        title{1} = 'Degree ';
                        if kk > length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            switch handles.parameters.thresh,
                                case 2, handles.display.data = vertcat(handles.corr.static.NDt);
                                    title{2} = '(session threshold) ';
                            end
                        else
                            switch handles.parameters.thresh,
                                case 2, handles.display.data = vertcat(handles.corr.static(kk).NDt);
                                    title{2} = '(session threshold) ';
                            end
                        end
                    case 5, %Node density
                        title{1} = 'Density ';
                        if kk > length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            switch handles.parameters.thresh,
                                case 2, handles.display.data = vertcat(handles.corr.static.NPt);
                                    title{2} = '(session threshold) ';
                            end
                        else
                            switch handles.parameters.thresh,
                                case 2, handles.display.data = vertcat(handles.corr.static(kk).NPt);
                                    title{2} = '(session threshold) ';
                            end
                        end
                end
                
            case 3, %Group dynamics
                plottype = 'bar';
                switch handles.parameters.metric,
                    case 6, %Within group connectivity
                        title{1} = 'Within group connectivity ';
                        if kk > length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = vertcat(handles.corr.static.FC_within);
                                    title{2} = ' ';
                                case 2, handles.display.data = vertcat(handles.corr.static.FCt_within);
                                    title{2} = '(session threshold) ';
                            end
                        else
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = vertcat(handles.corr.static(kk).FC_within);
                                    title{2} = ' ';
                                case 2, handles.display.data = vertcat(handles.corr.static(kk).FCt_within);
                                    title{2} = '(session threshold) ';
                            end
                        end
                    case 7, %Between group connectivity
                        title{1} = 'Between group connectivity ';
                        if kk > length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = vertcat(handles.corr.static.FC_between);
                                    title{2} = ' ';
                                case 2, handles.display.data = vertcat(handles.corr.static.FCt_between);
                                    title{2} = '(session threshold) ';
                            end
                        else
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = vertcat(handles.corr.static(kk).FC_between);
                                    title{2} = ' ';
                                case 2, handles.display.data = vertcat(handles.corr.static(kk).FCt_between);
                                    title{2} = '(session threshold) ';
                            end
                        end
                    case 8, %Within group strength
                        title{1} = 'Within group connectivity strength ';
                        if kk > length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = vertcat(handles.corr.static.S_within);
                                    title{2} = ' ';
                                case 2, handles.display.data = vertcat(handles.corr.static.St_within);
                                    title{2} = '(session threshold) ';
                            end
                        else
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = vertcat(handles.corr.static(kk).S_within);
                                    title{2} = ' ';
                                case 2, handles.display.data = vertcat(handles.corr.static(kk).St_within);
                                    title{2} = '(session threshold) ';
                            end
                        end
                    case 9, %Between group strength
                        title{1} = 'Between group connectivity strength ';
                        if kk > length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = vertcat(handles.corr.static.S_between);
                                    title{2} = ' ';
                                case 2, handles.display.data = vertcat(handles.corr.static.St_between);
                                    title{2} = '(session threshold) ';
                            end
                        else
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = vertcat(handles.corr.static(kk).S_between);
                                    title{2} = ' ';
                                case 2, handles.display.data = vertcat(handles.corr.static(kk).St_between);
                                    title{2} = '(session threshold) ';
                            end
                        end
                    case 10, %Within group degree
                        title{1} = 'Within group degree ';
                        if kk > length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            switch handles.parameters.thresh,
                                case 2, handles.display.data = vertcat(handles.corr.static.NDt_within);
                                    title{2} = '(session threshold) ';
                            end
                        else
                            switch handles.parameters.thresh,
                                case 2, handles.display.data = vertcat(handles.corr.static(kk).NDt_within);
                                    title{2} = '(session threshold) ';
                            end
                        end
                    case 11, %Between group degree
                        title{1} = 'Between group degree ';
                        if kk > length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            switch handles.parameters.thresh,
                                case 2, handles.display.data = vertcat(handles.corr.static.NDt_between);
                                    title{2} = '(session threshold) ';
                            end
                        else
                            switch handles.parameters.thresh,
                                case 2, handles.display.data = vertcat(handles.corr.static(kk).NDt_between);
                                    title{2} = '(session threshold) ';
                            end
                        end
                    case 12, %Within group density
                        title{1} = 'Within group density ';
                        if kk > length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            switch handles.parameters.thresh,
                                case 2, handles.display.data = vertcat(handles.corr.static.NPt_within);
                                    title{2} = '(session threshold) ';
                            end
                        else
                            switch handles.parameters.thresh,
                                case 2, handles.display.data = vertcat(handles.corr.static(kk).NPt_within);
                                    title{2} = '(session threshold) ';
                            end
                        end
                    case 13, %Between group density
                        title{1} = 'Between group density ';
                        if kk > length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            switch handles.parameters.thresh,
                                case 2, handles.display.data = vertcat(handles.corr.static.NPt_between);
                                    title{2} = '(session threshold) ';
                            end
                        else
                            switch handles.parameters.thresh,
                                case 2, handles.display.data = vertcat(handles.corr.static(kk).NPt_between);
                                    title{2} = '(session threshold) ';
                            end
                        end
                end
                
            case 4, %Group comparison
                plottype = 'bar';
                switch handles.parameters.metric,
                    case 6, %Within group connectivity
                        title{1} = 'Within group connectivity ';
                        if kk > length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = vertcat(handles.corr.static.FC_gwithin);
                                    title{2} = ' ';
                                case 2, handles.display.data = vertcat(handles.corr.static.FCt_gwithin);
                                    title{2} = '(session threshold) ';
                            end
                        else
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = vertcat(handles.corr.static(kk).FC_gwithin);
                                    title{2} = '- ';
                                case 2, handles.display.data = vertcat(handles.corr.static(kk).FCt_gwithin);
                            end
                        end
                    case 7, %Between group connectivity
                        title{1} = 'Between group connectivity ';
                        if kk > length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = vertcat(handles.corr.static.FC_gbetween);
                                    title{2} = ' ';
                                case 2, handles.display.data = vertcat(handles.corr.static.FCt_gbetween);
                                    title{2} = '(session threshold) ';
                            end
                        else
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = vertcat(handles.corr.static(kk).FC_gbetween);
                                    title{2} = ' ';
                                case 2, handles.display.data = vertcat(handles.corr.static(kk).FCt_gbetween);
                                    title{2} = '(session threshold) ';
                            end
                        end
                    case 8, %Within group connectivity strength
                        title{1} = 'Within group connectivity strength ';
                        if kk > length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = vertcat(handles.corr.static.S_gwithin);
                                    title{2} = ' ';
                                case 2, handles.display.data = vertcat(handles.corr.static.St_gwithin);
                                    title{2} = '(session threshold) ';
                            end
                        else
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = vertcat(handles.corr.static(kk).S_gwithin);
                                    title{2} = '- ';
                                case 2, handles.display.data = vertcat(handles.corr.static(kk).St_gwithin);
                            end
                        end
                    case 9, %Between group connectivity
                        title{1} = 'Between group connectivity strength ';
                        if kk > length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = vertcat(handles.corr.static.S_gbetween);
                                    title{2} = ' ';
                                case 2, handles.display.data = vertcat(handles.corr.static.St_gbetween);
                                    title{2} = '(session threshold) ';
                            end
                        else
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = vertcat(handles.corr.static(kk).S_gbetween);
                                    title{2} = ' ';
                                case 2, handles.display.data = vertcat(handles.corr.static(kk).St_gbetween);
                                    title{2} = '(session threshold) ';
                            end
                        end
                    case 10, %Within group degree
                        title{1} = 'Within group degree ';
                        if kk > length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            switch handles.parameters.thresh,
                                case 2, handles.display.data = vertcat(handles.corr.static.NDt_gwithin);
                                    title{2} = '(session threshold) ';
                            end
                        else
                            switch handles.parameters.thresh,
                                case 2, handles.display.data = vertcat(handles.corr.static(kk).NDt_gwithin);
                                    title{2} = '(session threshold) ';
                            end
                        end
                    case 11, %Between group degree
                        title{1} = 'Between group degree ';
                        if kk > length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            switch handles.parameters.thresh,
                                case 2, handles.display.data = vertcat(handles.corr.static.NDt_gbetween);
                                    title{2} = '(session threshold) ';
                            end
                        else
                            switch handles.parameters.thresh,
                                case 2, handles.display.data = vertcat(handles.corr.static(kk).NDt_gbetween);
                                    title{2} = '(session threshold) ';
                            end
                        end
                    case 12, %Within group density
                        title{1} = 'Within group density ';
                        if kk > length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            switch handles.parameters.thresh,
                                case 2, handles.display.data = vertcat(handles.corr.static.NPt_gwithin);
                                    title{2} = '(session threshold) ';
                            end
                        else
                            switch handles.parameters.thresh,
                                case 2, handles.display.data = vertcat(handles.corr.static(kk).NPt_gwithin);
                                    title{2} = '(session threshold)- ';
                            end
                        end
                    case 13, %Between group density
                        title{1} = 'Between group density ';
                        if kk > length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            switch handles.parameters.thresh,
                                case 2, handles.display.data = vertcat(handles.corr.static.NPt_gbetween);
                                    title{2} = '(session threshold) ';
                            end
                        else
                            switch handles.parameters.thresh,
                                case 2, handles.display.data = vertcat(handles.corr.static(kk).NPt_gbetween);
                                    title{2} = '(session threshold) ';
                            end
                        end
                end
                
            case 5, %Graph
                plottype = 'bar';
                switch handles.parameters.metric,
                    case 2, %Mean connectivity
                        title{1} = 'Mean connectivity ';
                        if kk > length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = vertcat(handles.corr.static.GC);
                                    title{2} = ' ';
                                case 2, handles.display.data = vertcat(handles.corr.static.GCt);
                                    title{2} = '(session threshold) ';
                            end
                        else
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = vertcat(handles.corr.static(kk).GC);
                                    title{2} = ' ';
                                case 2, handles.display.data = vertcat(handles.corr.static(kk).GCt);
                                    title{2} = '(session threshold) ';
                            end
                        end
                    case 3, %connectivity strength
                        title{1} = 'Connectivity strength ';
                        if kk > length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = vertcat(handles.corr.static.GS);
                                    title{2} = ' ';
                                case 2, handles.display.data = vertcat(handles.corr.static.GSt);
                                    title{2} = '(session threshold) ';
                            end
                        else
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = vertcat(handles.corr.static(kk).GS);
                                    title{2} = ' ';
                                case 2, handles.display.data = vertcat(handles.corr.static(kk).GSt);
                                    title{2} = '(session threshold) ';
                            end
                        end
                    case 4, %Node degree
                        title{1} = 'Degree ';
                        if kk > length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            switch handles.parameters.thresh,
                                case 2, handles.display.data = vertcat(handles.corr.static.GDt);
                                    title{2} = '(session threshold) ';
                            end
                        else
                            switch handles.parameters.thresh,
                                case 2, handles.display.data = vertcat(handles.corr.static(kk).GDt);
                                    title{2} = '(session threshold) ';
                            end
                        end
                    case 5, %Node density
                        title{1} = 'Density ';
                        if kk > length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            switch handles.parameters.thresh,
                                case 2, handles.display.data = vertcat(handles.corr.static.GPt);
                                    title{2} = '(session threshold) ';
                            end
                        else
                            switch handles.parameters.thresh,
                                case 2, handles.display.data = vertcat(handles.corr.static(kk).GPt);
                                    title{2} = '(session threshold) ';
                            end
                        end
                    case 6, %Within group connectivity
                        title{1} = 'Within group connectivity ';
                        if kk > length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = vertcat(handles.corr.static.GC_within);
                                    title{2} = ' ';
                                case 2, handles.display.data = vertcat(handles.corr.static.GCt_within);
                                    title{2} = '(session threshold) ';
                            end
                        else
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = vertcat(handles.corr.static(kk).GC_within);
                                    title{2} = ' ';
                                case 2, handles.display.data = vertcat(handles.corr.static(kk).GCt_within);
                                    title{2} = '(session threshold) ';
                            end
                        end
                    case 7, %Between group connectivity
                        title{1} = 'Between group connectivity ';
                        if kk > length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = vertcat(handles.corr.static.GC_between);
                                    title{2} = ' ';
                                case 2, handles.display.data = vertcat(handles.corr.static.GCt_between);
                                    title{2} = '(session threshold) ';
                            end
                        else
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = vertcat(handles.corr.static(kk).GC_between);
                                    title{2} = ' ';
                                case 2, handles.display.data = vertcat(handles.corr.static(kk).GCt_between);
                                    title{2} = '(session threshold) ';
                            end
                        end
                    case 8, %Within group strength
                        title{1} = 'Within group connectivity strength ';
                        if kk > length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = vertcat(handles.corr.static.GS_within);
                                    title{2} = ' ';
                                case 2, handles.display.data = vertcat(handles.corr.static.GSt_within);
                                    title{2} = '(session threshold) ';
                            end
                        else
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = vertcat(handles.corr.static(kk).GS_within);
                                    title{2} = ' ';
                                case 2, handles.display.data = vertcat(handles.corr.static(kk).GSt_within);
                                    title{2} = '(session threshold) ';
                            end
                        end
                    case 9, %Between group strength
                        title{1} = 'Between group connectivity strength ';
                        if kk > length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = vertcat(handles.corr.static.GS_between);
                                    title{2} = ' ';
                                case 2, handles.display.data = vertcat(handles.corr.static.GSt_between);
                                    title{2} = '(session threshold) ';
                            end
                        else
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = vertcat(handles.corr.static(kk).GS_between);
                                    title{2} = ' ';
                                case 2, handles.display.data = vertcat(handles.corr.static(kk).GSt_between);
                                    title{2} = '(session threshold) ';
                            end
                        end
                    case 10, %Within group degree
                        title{1} = 'Within group degree ';
                        if kk > length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            switch handles.parameters.thresh,
                                case 2, handles.display.data = vertcat(handles.corr.static.GDt_within);
                                    title{2} = '(session threshold) ';
                            end
                        else
                            switch handles.parameters.thresh,
                                case 2, handles.display.data = vertcat(handles.corr.static(kk).GDt_within);
                                    title{2} = '(session threshold) ';
                            end
                        end
                    case 11, %Between group degree
                        title{1} = 'Between group degree ';
                        if kk > length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            switch handles.parameters.thresh,
                                case 2, handles.display.data = vertcat(handles.corr.static.GDt_between);
                                    title{2} = '(session threshold) ';
                            end
                        else
                            switch handles.parameters.thresh,
                                case 2, handles.display.data = vertcat(handles.corr.static(kk).GDt_between);
                                    title{2} = '(session threshold) ';
                            end
                        end
                    case 12, %Within group density
                        title{1} = 'Within group density ';
                        if kk > length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            switch handles.parameters.thresh,
                                case 2, handles.display.data = vertcat(handles.corr.static.GPt_within);
                                    title{2} = '(session threshold) ';
                            end
                        else
                            switch handles.parameters.thresh,
                                case 2, handles.display.data = vertcat(handles.corr.static(kk).GPt_within);
                                    title{2} = '(session threshold) ';
                            end
                        end
                    case 13, %Between group density
                        title{1} = 'Between group density ';
                        if kk > length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            switch handles.parameters.thresh,
                                case 2, handles.display.data = vertcat(handles.corr.static.GPt_between);
                                    title{2} = '(session threshold) ';
                            end
                        else
                            switch handles.parameters.thresh,
                                case 2, handles.display.data = vertcat(handles.corr.static(kk).GPt_between);
                                    title{2} = '(session threshold) ';
                            end
                        end
                end
        end
    elseif handles.parameters.connectivity==2, %DYNAMIC connectivity
        title{4} = '[dynamic] ';
        switch handles.parameters.level, %find specified level
            case 1, %Covariance
                plottype = 'covar';
                switch handles.parameters.metric,
                    case 1, %mean connectivity
                        title{1} = 'Mean connectivity ';
                        switch handles.parameters.thresh,
                            case 1, handles.display.data = handles.corr.dyn(kk).R;
                                title{2} = ' ';
                            case 3, handles.display.data = handles.corr.dyn(kk).Rt;
                                title{2} = '(condition threshold) ';
                        end
                    case 2, %mean connectivity strength
                        title{1} = 'Mean connectivity strength ';
                        switch handles.parameters.thresh,
                            case 1, handles.display.data = handles.corr.dyn(kk).S;
                                title{2} = ' ';
                            case 3, handles.display.data = handles.corr.dyn(kk).St;
                                title{2} = '(condition threshold) ';
                        end
                    case 3, %Temporal degree
                        handles.display.data = handles.corr.dyn(kk).Tt;
                        title{1} = 'Temporal degree ';
                    case 4, %Temporal density
                        handles.display.data = handles.corr.dyn(kk).Tp;
                        title{1} = 'Temporal density ';
                end
            case 2, %Subject
                plottype = 'bar';
                switch handles.parameters.metric,
                    case 1, %mean connectivity
                        title{1} = 'Mean connectivity ';
                        if kk >length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = vertcat(handles.corr.dyn.FC);
                                    title{2} = ' ';
                                case 3, handles.display.data = vertcat(handles.corr.dyn.FCt);
                                    title{2} = '(condition threshold) ';
                            end
                        else
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = handles.corr.dyn(kk).FC;
                                    title{2} = ' ';
                                case 3, handles.display.data = handles.corr.dyn(kk).FCt;
                                    title{2} = '(condition threshold) ';
                            end
                        end
                    case 2, %Mean connectivity strength
                        title{1} = 'Mean connectivity strength ';
                        if kk >length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = vertcat(handles.corr.dyn.FS);
                                    title{2} = ' ';
                                case 3, handles.display.data = vertcat(handles.corr.dyn.FSt);
                                    title{2} = '(condition threshold) ';
                            end
                        else
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = handles.corr.dyn(kk).FS;
                                    title{2} = ' ';
                                case 3, handles.display.data = handles.corr.dyn(kk).FSt;
                                    title{2} = '(condition threshold) ';
                            end
                        end
                    case 3, %Mean temporal degree
                        title{1} = 'Mean temporal degree ';
                        title{2} = '(condition threshold) ';
                        if kk >length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            handles.display.data = vertcat(handles.corr.dyn.FTt);
                        else
                            handles.display.data = handles.corr.dyn(kk).FTt;
                        end
                    case 4, %Mean temporal density
                        title{1} = 'Mean temporal density ';
                        title{2} = '(condition threshold) ';
                        if kk >length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            handles.display.data = vertcat(handles.corr.dyn.FTp);
                        else
                            handles.display.data = handles.corr.dyn(kk).FTp;
                        end
                    case 5, %Mean degree
                        title{1} = 'Mean degree ';
                        title{2} = '(condition threshold) ';
                        if kk >length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            handles.display.data = vertcat(handles.corr.dyn.ND);
                        else
                            handles.display.data = handles.corr.dyn(kk).ND;
                        end
                    case 6, %Sum degree
                        title{1} = 'Sum degree ';
                        title{2} = '(condition threshold) ';
                        if kk >length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            handles.display.data = vertcat(handles.corr.dyn.NDs);
                        else
                            handles.display.data = handles.corr.dyn(kk).NDs;
                        end
                    case 7, %Mean density
                        title{1} = 'Mean density ';
                        title{2} = '(condition threshold) ';
                        if kk >length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            handles.display.data = vertcat(handles.corr.dyn.NDp);
                        else
                            handles.display.data = handles.corr.dyn(kk).NDp;
                        end
                end
            case 3, %Group Dynamics
                plottype = 'bar';
                switch handles.parameters.metric,
                    case 8, %Within group connectivity
                        title{1} = 'Within group connectivity ';
                        if kk >length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = vertcat(handles.corr.dyn.FC_within);
                                    title{2} = ' ';
                                case 3, handles.display.data = vertcat(handles.corr.dyn.FCt_within);
                                    title{2} = '(condition threshold) ';
                            end
                        else
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = handles.corr.dyn(kk).FC_within;
                                    title{2} = ' ';
                                case 3, handles.display.data = handles.corr.dyn(kk).FCt_within;
                                    title{2} = '(condition threshold) ';
                            end
                        end
                    case 9, %Between group connectivity
                        title{1} = 'Between group connectivity ';
                        if kk >length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = vertcat(handles.corr.dyn.FC_between);
                                    title{2} = ' ';
                                case 3, handles.display.data = vertcat(handles.corr.dyn.FCt_between);
                                    title{2} = '(condition threshold) ';
                            end
                        else
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = handles.corr.dyn(kk).FC_between;
                                    title{2} = ' ';
                                case 3, handles.display.data = handles.corr.dyn(kk).FCt_between;
                                    title{2} = '(condition threshold) ';
                            end
                        end
                    case 10, %Within group connectivity strength
                        title{1} = 'Within group connectivity strength ';
                        if kk >length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = vertcat(handles.corr.dyn.FS_within);
                                    title{2} = ' ';
                                case 3, handles.display.data = vertcat(handles.corr.dyn.FSt_within);
                                    title{2} = '(condition threshold) ';
                            end
                        else
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = handles.corr.dyn(kk).FS_within;
                                    title{2} = ' ';
                                case 3, handles.display.data = handles.corr.dyn(kk).FSt_within;
                                    title{2} = '(condition threshold) ';
                            end
                        end
                    case 11, %Between group connectivity strength
                        title{1} = 'Between group connectivity strength ';
                        if kk >length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = vertcat(handles.corr.dyn.FS_between);
                                    title{2} = ' ';
                                case 3, handles.display.data = vertcat(handles.corr.dyn.FSt_between);
                                    title{2} = '(condition threshold) ';
                            end
                        else
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = handles.corr.dyn(kk).FC_between;
                                    title{2} = ' ';
                                case 3, handles.display.data = handles.corr.dyn(kk).FCt_between;
                                    title{2} = '(condition threshold) ';
                            end
                        end
                    case 12, %Within group degree
                        title{1} = 'Within group degree ';
                        title{2} = '(condition threshold) ';
                        if kk >length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            handles.display.data = vertcat(handles.corr.dyn.NDt_within);
                        else
                            handles.display.data = handles.corr.dyn(kk).NDt_within;
                        end
                    case 13, %Between group degree
                        title{1} = 'Between group degree ';
                        title{2} = '(condition threshold) ';
                        if kk >length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            handles.display.data = vertcat(handles.corr.dyn.NDt_between);
                        else
                            handles.display.data = handles.corr.dyn(kk).NDt_between;
                        end
                    case 14, %Within group density
                        title{1} = 'Within group density ';
                        title{2} = '(condition threshold) ';
                        if kk >length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            handles.display.data = vertcat(handles.corr.dyn.NPt_within);
                        else
                            handles.display.data = handles.corr.dyn(kk).NPt_within;
                        end
                    case 15, %Between group density
                        title{1} = 'Between group density ';
                        title{2} = '(condition threshold) ';
                        if kk >length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            handles.display.data = vertcat(handles.corr.dyn.NPt_between);
                        else
                            handles.display.data = handles.corr.dyn(kk).NPt_between;
                        end
                end
            case 5, %Graph Measures
                plottype = 'bar';
                switch handles.parameters.metric,
                    
                    case 1, %mean connectivity
                        title{1} = 'Mean connectivity ';
                        if kk >length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = vertcat(handles.corr.dyn.GC);
                                    title{2} = ' ';
                                case 3, handles.display.data = vertcat(handles.corr.dyn.GCt);
                                    title{2} = '(condition threshold) ';
                            end
                        else
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = handles.corr.dyn(kk).GC;
                                    title{2} = ' ';
                                case 3, handles.display.data = handles.corr.dyn(kk).GCt;
                                    title{2} = '(condition threshold) ';
                            end
                        end
                    case 2, %Mean connectivity strength
                        title{1} = 'Mean connectivity strength ';
                        if kk >length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = vertcat(handles.corr.dyn.GS);
                                    title{2} = ' ';
                                case 3, handles.display.data = vertcat(handles.corr.dyn.GSt);
                                    title{2} = '(condition threshold) ';
                            end
                        else
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = handles.corr.dyn(kk).GS;
                                    title{2} = ' ';
                                case 3, handles.display.data = handles.corr.dyn(kk).GSt;
                                    title{2} = '(condition threshold) ';
                            end
                        end
                    case 3, %Mean temporal degree
                        title{1} = 'Mean temporal degree ';
                        title{2} = '(condition threshold) ';
                        if kk >length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            handles.display.data = vertcat(handles.corr.dyn.GTt);
                        else
                            handles.display.data = handles.corr.dyn(kk).GTt;
                        end
                    case 4, %Mean temporal density
                        title{1} = 'Mean temporal density ';
                        title{2} = '(condition threshold) ';
                        if kk >length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            handles.display.data = vertcat(handles.corr.dyn.GTp);
                        else
                            handles.display.data = handles.corr.dyn(kk).GTp;
                        end
                    case 5, %Mean degree
                        title{1} = 'Mean degree ';
                        title{2} = '(condition threshold) ';
                        if kk >length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            handles.display.data = vertcat(handles.corr.dyn.GD);
                        else
                            handles.display.data = handles.corr.dyn(kk).GD;
                        end
                    case 6, %Sum degree
                        title{1} = 'Sum degree ';
                        title{2} = '(condition threshold) ';
                        if kk >length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            handles.display.data = vertcat(handles.corr.dyn.GDs);
                        else
                            handles.display.data = handles.corr.dyn(kk).GDs;
                        end
                    case 7, %Mean density
                        title{1} = 'Mean density ';
                        title{2} = '(condition threshold) ';
                        if kk >length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            handles.display.data = vertcat(handles.corr.dyn.GDp);
                        else
                            handles.display.data = handles.corr.dyn(kk).GDp;
                        end
                    case 8, %Within group connectivity
                        title{1} = 'Within group connectivity ';
                        if kk >length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = vertcat(handles.corr.dyn.GC_within);
                                    title{2} = ' ';
                                case 3, handles.display.data = vertcat(handles.corr.dyn.GCt_within);
                                    title{2} = '(condition threshold) ';
                            end
                        else
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = handles.corr.dyn(kk).GC_within;
                                    title{2} = ' ';
                                case 3, handles.display.data = handles.corr.dyn(kk).GCt_within;
                                    title{2} = '(condition threshold) ';
                            end
                        end
                    case 9, %Between group connectivity
                        title{1} = 'Between group connectivity ';
                        if kk >length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = vertcat(handles.corr.dyn.GC_between);
                                    title{2} = ' ';
                                case 3, handles.display.data = vertcat(handles.corr.dyn.GCt_between);
                                    title{2} = '(condition threshold) ';
                            end
                        else
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = handles.corr.dyn(kk).GC_between;
                                    title{2} = ' ';
                                case 3, handles.display.data = handles.corr.dyn(kk).GCt_between;
                                    title{2} = '(condition threshold) ';
                            end
                        end
                    case 10, %Within group connectivity strength
                        title{1} = 'Within group connectivity strength ';
                        if kk >length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = vertcat(handles.corr.dyn.GS_within);
                                    title{2} = ' ';
                                case 3, handles.display.data = vertcat(handles.corr.dyn.GSt_within);
                                    title{2} = '(condition threshold) ';
                            end
                        else
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = handles.corr.dyn(kk).GS_within;
                                    title{2} = ' ';
                                case 3, handles.display.data = handles.corr.dyn(kk).GSt_within;
                                    title{2} = '(condition threshold) ';
                            end
                        end
                    case 11, %Between group connectivity strength
                        title{1} = 'Between group connectivity strength ';
                        if kk >length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = vertcat(handles.corr.dyn.GS_between);
                                    title{2} = ' ';
                                case 3, handles.display.data = vertcat(handles.corr.dyn.GSt_between);
                                    title{2} = '(condition threshold) ';
                            end
                        else
                            switch handles.parameters.thresh,
                                case 1, handles.display.data = handles.corr.dyn(kk).GC_between;
                                    title{2} = ' ';
                                case 3, handles.display.data = handles.corr.dyn(kk).GCt_between;
                                    title{2} = '(condition threshold) ';
                            end
                        end
                    case 12, %Within group degree
                        title{1} = 'Within group degree ';
                        title{2} = '(condition threshold) ';
                        if kk >length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            handles.display.data = vertcat(handles.corr.dyn.GDt_within);
                        else
                            handles.display.data = handles.corr.dyn(kk).GDt_within;
                        end
                    case 13, %Between group degree
                        title{1} = 'Between group degree ';
                        title{2} = '(condition threshold) ';
                        if kk >length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            handles.display.data = vertcat(handles.corr.dyn.GDt_between);
                        else
                            handles.display.data = handles.corr.dyn(kk).GDt_between;
                        end
                    case 14, %Within group density
                        title{1} = 'Within group density ';
                        title{2} = '(condition threshold) ';
                        if kk >length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            handles.display.data = vertcat(handles.corr.dyn.GPt_within);
                        else
                            handles.display.data = handles.corr.dyn(kk).GPt_within;
                        end
                    case 15, %Between group density
                        title{1} = 'Between group density ';
                        title{2} = '(condition threshold) ';
                        if kk >length(handles.conditions), %if compare conditions
                            plottype = 'compare';
                            handles.display.data = vertcat(handles.corr.dyn.GPt_between);
                        else
                            handles.display.data = handles.corr.dyn(kk).GPt_between;
                        end
                end
            case 6, %Whole session inspection
                plottype = 'covartime';
                switch handles.parameters.metric,
                    case 1, %Connectivity
                        title{1} = 'Connectivity ';
                        switch handles.parameters.thresh,
                            case 1, handles.display.data = handles.corr.dynwhole.FC;
                                title{2} = ' ';
                            case 2, handles.display.data = handles.corr.dynwhole.FCt;
                                title{2} = '(whole session threshold) ';
                        end
                    case 2, %Connectivity strength
                        title{1} = 'Connectivity strength ';
                        switch handles.parameters.thresh,
                            case 1, handles.display.data = handles.corr.dynwhole.S;
                                title{2} = ' ';
                            case 2, handles.display.data = handles.corr.dynwhole.St;
                                title{2} = '(whole session threshold) ';
                        end
                    case 3, %Degree
                        title{1} = 'Degree ';
                        switch handles.parameters.thresh,
                            case 1, handles.display.data = handles.corr.dynwhole.degree;
                                title{2} = '(whole session threshold) ';
                        end
                    case 4, %Density
                        title{1} = 'Density ';
                        switch handles.parameters.thresh,
                            case 1, handles.display.data = handles.corr.dynwhole.degreep;
                                title{2} = '(whole session threshold) ';
                        end
                    case 5, %Within group connectivity
                        title{1} = 'Within group connectivity ';
                        switch handles.parameters.thresh,
                            case 1, handles.display.data = handles.corr.dynwhole.FC_within;
                                title{2} = ' ';
                            case 2, handles.display.data = handles.corr.dynwhole.FCt_within;
                                title{2} = '(whole session threshold) ';
                        end
                    case 6, %Between group connectivity
                        title{1} = 'Between group connectivity ';
                        switch handles.parameters.thresh,
                            case 1, handles.display.data = handles.corr.dynwhole.FC_between;
                                title{2} = ' ';
                            case 2, handles.display.data = handles.corr.dynwhole.FCt_between;
                                title{2} = '(whole session threshold) ';
                        end
                    case 7, %Within group degree
                        title{1} = 'Within group degree ';
                        switch handles.parameters.thresh,
                            case 1, handles.display.data = handles.corr.dynwhole.NDt_within;
                                title{2} = '(whole session threshold) ';
                        end
                    case 8, %Between group degree
                        title{1} = 'Between group degree ';
                        switch handles.parameters.thresh,
                            case 1, handles.display.data = handles.corr.dynwhole.NDt_between;
                                title{2} = '(whole session threshold) ';
                        end
                    case 9, %Within group density
                        title{1} = 'Within group density ';
                        switch handles.parameters.thresh,
                            case 1, handles.display.data = handles.corr.dynwhole.NPt_within;
                                title{2} = '(whole session threshold) ';
                        end
                    case 10, %Between group density
                        title{1} = 'Between group density ';
                        switch handles.parameters.thresh,
                            case 1, handles.display.data = handles.corr.dynwhole.NPt_between;
                                title{2} = '(whole session threshold) ';
                        end
                end
        end
    end
    
    %% PLOT GRAPH METRICS
    % Select correct subplot to plot in
    if isfield(handles,'edaaxes'),
        set(handles.edaaxes,'Visible','off');
        axes(handles.edaaxes); cla;
    end
    if get(handles.corleft_radio,'Value')==1,
        axes(handles.GM1axes);
        set(handles.GM1axes,'Visible','on');
        set(handles.gmlefttext,'Visible','on');
        set(handles.gmlefttext,'String',cell2mat(title));
    elseif get(handles.corright_radio,'Value')==1,
        axes(handles.GM2axes);
        set(handles.GM2axes,'Visible','on');
        set(handles.gmrighttext,'Visible','on');
        set(handles.gmrighttext,'String',cell2mat(title));
    end
    %Select plot type
    switch plottype
        case 'bar',
            bar(gca,handles.display.data);
            if handles.parameters.level==5, %Graph level
                set(gca,'XTickLabel',{handles.conditions.name});
            end
        case 'compare',
            bar(gca,handles.display.data');
            if handles.parameters.level==5, %Graph level
                set(gca,'XTickLabel',{handles.conditions.name});
            else
                legend(gca,{handles.conditions.name});
            end
        case 'covar',
            imagesc(handles.display.data);
            colorbar; colormap('jet');
        case 'covartime',
            imagesc(handles.display.data);
            colorbar; colormap('jet');
            set(gca,'YTick',1:10:length(handles.time.DYNmins));
            inx = get(gca,'YTick');
            Ylabel = (round(handles.time.DYNmins(inx)*10))./10;
            set(gca,'YTickLabel',Ylabel);
            lim = caxis;
            if lim(2)<1,
            caxis([-1 1]);
            end
            
    end
    
    %% END PLOT
end

function loadwork_push_Callback(hObject, eventdata, handles)
filename = 0;
[filename, pathname] = uigetfile({'*.mat'},'Select file containing saved data...');
if isequal(filename,0), return; end %if nothing selected, stop function
handles = clearvariables(handles);
load([pathname filename]);
handles.sdir = OUTPUT.sdir;
handles.data = OUTPUT.data;
handles.parameters = OUTPUT.parameters;
handles.time = OUTPUT.time;
handles.EDA = OUTPUT.EDA;
handles.ACC = OUTPUT.ACC;
handles.groups = OUTPUT.groups;
handles.group_subj = OUTPUT.group_subj;
handles.conditions = OUTPUT.conditions;
handles.surrogate = OUTPUT.surrogate;
handles.corr = OUTPUT.corr;
clear OUTPUT

%fill subject list
set(handles.vissubj_list,'String',{handles.data.subj});
set(handles.vissubj_list,'Value',1);
%fill group lists (if groups exist)
if ~isempty(handles.groups),
    grouplist = [];
    for gg = 1:length(handles.groups), %for each group
        templist = [];
        for ii = 1:length(handles.groups{gg}), %for each groupsubj index selected
            if ii<length(handles.groups{gg}),
                templist = [templist handles.groups{gg}{ii} ',']; %create templist for group string
            else  templist = [templist handles.groups{gg}{ii}];
            end
        end
        templist = cellstr(templist);
        grouplist = [grouplist; templist];
    end
    set(handles.groupdefined_list,'String',grouplist);
    set(handles.visgroups_list,'String',grouplist);
    set(handles.groupsubj_list,'String',handles.group_subj); %update available subjects
set(handles.groupsubj_list,'Value',1); %update list selection
end

%fill condition list
if ~isempty(handles.conditions),
    set(handles.cond_list,'String',{handles.conditions.string});
end
%Visibility of panels depending on step of preprocessing
set(handles.vis_uipanel,'Visible','on');
if isempty(handles.conditions) && isempty(handles.groups) && ~isfield(handles.parameters,accpeak),
    set(handles.preprocess_uipanel,'Visible','on');
else
    set(handles.align_uipanel,'Visible','on');
    set(handles.group_uipanel,'Visible','on');
    set(handles.condition_uipanel,'Visible','on');
    set(handles.preprocess_uipanel,'Visible','on');
    set(handles.GM_uipanel,'Visible','on');
end
update_EDAdisplay(handles);
guidata(hObject, handles);

function savework_push_Callback(hObject, eventdata, handles)
OUTPUT.sdir = handles.sdir;
OUTPUT.data = handles.data;
OUTPUT.parameters = handles.parameters;
OUTPUT.time = handles.time;
OUTPUT.EDA = handles.EDA;
OUTPUT.ACC = handles.ACC;
OUTPUT.groups = handles.groups;
OUTPUT.group_subj = handles.group_subj;
OUTPUT.conditions = handles.conditions;
OUTPUT.surrogate = handles.surrogate;
OUTPUT.corr = handles.corr;
uisave('OUTPUT','preprocessed_data.mat');
clear OUTPUT
guidata(hObject, handles);

function push2mat_push_Callback(hObject, eventdata, handles)
OUTPUT.sdir = handles.sdir;
OUTPUT.data = handles.data;
OUTPUT.parameters = handles.parameters;
OUTPUT.time = handles.time;
OUTPUT.EDA = handles.EDA;
OUTPUT.ACC = handles.ACC;
OUTPUT.groups = handles.groups;
OUTPUT.group_subj = handles.group_subj;
OUTPUT.conditions = handles.conditions;
OUTPUT.surrogate = handles.surrogate;
OUTPUT.corr = handles.corr;
OUTPUT.display = handles.display;
assignin('base','OUTPUT',OUTPUT);
clear OUTPUT
guidata(hObject, handles);

function vissubj_toggle_Callback(hObject, eventdata, handles)
if get(handles.vissubj_toggle,'Value')==1,
    ideal_height = length([handles.data])*handles.height_per_listbox_line+1;
    if ideal_height < handles.select_subj_dropdown_y
        listbox_height = ideal_height;
    else
        listbox_height = handles.select_subj_dropdown_y - 1;
    end
    listbox_y = handles.select_subj_dropdown_y - listbox_height - handles.uipanel_left_y_offset;
    set(handles.vissubj_list,'Position',[handles.hidden_subjs_listbox_position(1) listbox_y handles.hidden_subjs_listbox_position(3) listbox_height]);
    set(handles.vissubj_list,'Visible','on');
else
    set(handles.vissubj_list,'Visible','off');
    set(handles.vissubj_list,'Position',handles.hidden_subjs_listbox_position);
end
guidata(hObject, handles);

function visgroups_toggle_Callback(hObject, eventdata, handles)
if get(handles.visgroups_toggle,'Value')==1,
    ideal_height = length(handles.groups)*handles.height_per_listbox_line+1;
    if ideal_height < handles.select_groups_dropdown_y
        listbox_height = ideal_height;
    else
        listbox_height = handles.select_groups_dropdown_y - 1;
    end
    listbox_y = handles.select_groups_dropdown_y - listbox_height - handles.uipanel_left_y_offset;
    set(handles.visgroups_list,'Position',[handles.hidden_groups_listbox_position(1) listbox_y handles.hidden_groups_listbox_position(3) listbox_height]);
    set(handles.visgroups_list,'Visible','on');
else
    set(handles.visgroups_list,'Visible','off');
    set(handles.visgroups_list,'Position',handles.hidden_groups_listbox_position);
end
guidata(hObject, handles);

function visreset_push_Callback(hObject, eventdata, handles)
set(handles.GM1axes,'Visible','off');
set(handles.GM2axes,'Visible','off');
axes(handles.GM1axes); cla
axes(handles.GM2axes); cla
set(handles.gmlefttext,'Visible','off');
set(handles.gmrighttext,'Visible','off');

handles.edaaxes = axes('Position',handles.pos);
handles.parameters.corr=0;
set(handles.visamp_check,'Value',1);
set(handles.vistime_check,'Value',1);
set(handles.vissubj_list,'Value',1);
set(handles.visEDA_radio,'Value',1);
handles = update_EDAdisplay(handles);
set(handles.vistimemin_edit,'String',num2str(handles.display.mintime));
set(handles.vistimemax_edit,'String',num2str(handles.display.maxtime));
set(handles.visampmin_edit,'String',num2str(handles.display.minamp));
set(handles.visampmax_edit,'String',num2str(handles.display.maxamp));
guidata(hObject, handles);

function visonoff_toggle_Callback(hObject, eventdata, handles)
handles = update_EDAdisplay(handles);
guidata(hObject, handles);

function visamp_check_Callback(hObject, eventdata, handles)
handles = update_EDAdisplay(handles);
guidata(hObject, handles);

function visampmax_edit_Callback(hObject, eventdata, handles)
handles = update_EDAdisplay(handles);
guidata(hObject, handles);

function visampmax_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function visampmin_edit_Callback(hObject, eventdata, handles)
handles = update_EDAdisplay(handles);
guidata(hObject, handles);

function visampmin_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function vistimemin_edit_Callback(hObject, eventdata, handles)
handles = update_EDAdisplay(handles);
guidata(hObject, handles);

function vistimemin_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function vistimemax_edit_Callback(hObject, eventdata, handles)
handles = update_EDAdisplay(handles);
guidata(hObject, handles);

function vistimemax_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function vistime_check_Callback(hObject, eventdata, handles)
handles = update_EDAdisplay(handles);
guidata(hObject, handles);

function vissubj_list_Callback(hObject, eventdata, handles)
handles = update_EDAdisplay(handles);
guidata(hObject, handles);


function vissubj_list_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function visgroups_list_Callback(hObject, eventdata, handles)
selectedgroup = get(handles.visgroups_list,'Value');
subjs = handles.groups{selectedgroup};
subj = [];
for ii = 1:length(subjs),
    subj = [subj find(strcmp(subjs{ii},{handles.data.subj}))];
end
set(handles.vissubj_list,'Value',subj);
handles = update_EDAdisplay(handles);
guidata(hObject, handles);

function visgroups_list_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function status_edit_Callback(hObject, eventdata, handles)

function status_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function visremovesubj_push_Callback(hObject, eventdata, handles)
if length(get(handles.vissubj_list,'Value'))==1,
    inx = get(handles.vissubj_list,'Value');
    deletesubj = handles.data(inx).subj;
    choice = questdlg(['Do you wish to remove subject ' deletesubj '? This will reset defined groups and preprocessing.'], 'Yes', 'No');
    switch choice
        case 'Yes'
             handles.data(inx) = [];
             handles.EDA = [];
             handles.ACC = [];
             handles.time = [];
             handles.display = [];
             
             for ii = 1:length(handles.data),
                 %EDA
                 handles.data(ii).EDAsamplerate=handles.data(ii).urEDA(2,1); %get samplerate
                 handles.data(ii).EDAon=handles.data(ii).urEDA(1,1); %EDA time on
                 handles.data(ii).EDAtimestamps=handles.data(ii).EDAon:1/handles.data(ii).EDAsamplerate: handles.data(ii).EDAon + (length(handles.data(ii).urEDA(3:end,1))/handles.data(ii).EDAsamplerate)-(1/handles.data(ii).EDAsamplerate);
                 handles.data(ii).EDAsecs=handles.data(ii).EDAtimestamps-handles.data(ii).EDAon;
                 handles.data(ii).EDAmins=handles.data(ii).EDAsecs./60;
                 handles.data(ii).EDAoff=handles.data(ii).EDAtimestamps(end);
                 
                 %ACC
                 if ~isempty(handles.data(ii).urACC), %if there's ACC data
                     handles.data(ii).ACCsamplerate=handles.data(ii).urACC(2,1); %get samplerate
                     handles.data(ii).ACCon=handles.data(ii).urACC(1,1); %EDA time on
                     if handles.data(ii).EDAon~=handles.data(ii).ACCon,
                         errordlg(['Warning - EDA and accelerometer data for subject ' handles.data(ii).subj ' have non-matching timestamps. Errors in computation may result. Please check that data files are correct']);
                     end
                     handles.data(ii).ACCtimestamps=handles.data(ii).ACCon:1/handles.data(ii).ACCsamplerate: handles.data(ii).ACCon + (length(handles.data(ii).urACC(3:end,1))/handles.data(ii).ACCsamplerate)-(1/handles.data(ii).ACCsamplerate);
                     handles.data(ii).ACCsecs=handles.data(ii).ACCtimestamps-handles.data(ii).ACCon;
                     handles.data(ii).ACCmins=handles.data(ii).ACCsecs./60;
                     handles.data(ii).ACCoff=handles.data(ii).ACCtimestamps(end);
                 end
                 
                 %GROUP DEFINITION (empty)
                 handles.data(ii).group = [];
                 handles.data(ii).EDA = [];
                 handles.data(ii).ACC = [];
             end
            
            %% RELOAD DATA VARIABLES
            %% Construct timelines
            set(handles.status_edit,'String','Temporal alignment of datasets...'); pause(0.01);
            %Check samplerate consistency
            if std([handles.data.EDAsamplerate])~=0,
                errordlg('Warning - Inconsistent sample rates between EDA files. Errors in computation will result. Please check data files');
            end
            if std([handles.data.ACCsamplerate])~=0,
                errordlg('Warning - Inconsistent sample rates between EDA files. Errors in computation will result. Please check data files');
            end
            %Copy samplerate to handles
            handles.parameters.EDAsamplerate=handles.data(1).EDAsamplerate;
            handles.parameters.ACCsamplerate=mean([handles.data.ACCsamplerate]);
            %Construct master timeline
            EDAtime=(min([handles.data.EDAon]):(1/handles.parameters.EDAsamplerate):max([handles.data.EDAoff]))';
            ACCtime=(min([handles.data.ACCon]):(1/handles.parameters.ACCsamplerate):max([handles.data.ACCoff]))';
            %Temporal alignment / pad with NaN
            for ii = 1:length(handles.data),
                
                [inx, ~] = find(EDAtime==handles.data(ii).EDAon);
                handles.data(ii).EDA(1:length(EDAtime),1)=NaN;
                handles.data(ii).EDA(inx:inx+length(handles.data(ii).urEDA(3:end-1)),1)=handles.data(ii).urEDA(3:end);
                if ~isempty(handles.data(ii).urACC),
                    
                    [inx, ~] = find(ACCtime==handles.data(ii).ACCon);
                    handles.data(ii).ACC(1:length(ACCtime),1:3)=NaN;
                    handles.data(ii).ACC(inx:inx+length(handles.data(ii).urACC(3:end-1,1)),:)=handles.data(ii).urACC(3:end,:);
                end
            end
            
            handles.time.EDAtimestamps=EDAtime;
            handles.time.EDAsecs=EDAtime-EDAtime(1);
            handles.time.EDAmins=handles.time.EDAsecs./60;
            handles.time.ACCtimestamps=ACCtime;
            handles.time.ACCsecs=ACCtime-ACCtime(1);
            handles.time.ACCmins=handles.time.ACCsecs./60;
            
            %Compute summed accelerometer signals
            for ii=1:length(handles.data),
                if ~isempty(handles.data(ii).urACC),
                    handles.data(ii).ACC=sum(abs(handles.data(ii).ACC),2);
                end
            end
            
            %Copy original data to current data handles
            handles.EDA=[handles.data.EDA];
            handles.ACC=[handles.data.ACC];
            
            %Set default times
            set(handles.vistimemin_edit,'String',num2str(handles.time.EDAmins(1)));
            set(handles.vistimemax_edit,'String',num2str(handles.time.EDAmins(end)));
            
            %Populate Subj Lists
            set(handles.vissubj_list,'String',{handles.data.subj});
            set(handles.vissubj_list,'Value',1);
            handles.group_subj = {handles.data.subj};
            set(handles.groupsubj_list,'String',handles.group_subj); %update available subjects
            set(handles.groupsubj_list,'Value',1); %update list selection
            set(handles.visgroups_list,'String','No groups defined');
            set(handles.status_edit,'String',['Dataset - ' handles.sdir ' ready for preprocessing']); pause(0.01);
            guidata(hObject, handles);
            update_EDAdisplay(handles);
            
        case 'No'
            return
    end
end
guidata(hObject, handles);

function visavg_check_Callback(hObject, eventdata, handles)
handles = update_EDAdisplay(handles);
guidata(hObject, handles);

function visreport_list_Callback(hObject, eventdata, handles)

function visreport_list_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function vistime_uibutton_SelectionChangedFcn(hObject, eventdata, handles)
handles = update_EDAdisplay(handles);
guidata(hObject, handles);

function visamp_uibutton_SelectionChangedFcn(hObject, eventdata, handles)
handles = update_EDAdisplay(handles);
guidata(hObject, handles);

function uibuttongroup6_SelectionChangedFcn(hObject, eventdata, handles)
handles = update_EDAdisplay(handles);
guidata(hObject, handles);

function groupsubj_toggle_Callback(hObject, eventdata, handles)
%Adjust list length
if get(handles.groupsubj_toggle,'Value')==1,
    listbox_height = .5;
    listbox_y = handles.select_group_dropdown_y - listbox_height - handles.uipanel_group_y_offset;
    set(handles.groupsubj_list,'Position',[handles.hidden_group_listbox_position(1) listbox_y handles.hidden_group_listbox_position(3) listbox_height]);
    set(handles.groupsubj_list,'Visible','on');
else
    set(handles.groupsubj_list,'Visible','off');
end
guidata(hObject, handles);

function groupsubj_list_Callback(hObject, eventdata, handles)

function groupsubj_list_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function groupdefined_list_Callback(hObject, eventdata, handles)

function groupdefined_list_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function groupcreate_push_Callback(hObject, eventdata, handles)
tempgroup = get(handles.groupsubj_list,'Value'); %indices of selected subj
if ~isempty(get(handles.groupsubj_list,'Value')),
    %Make new list entry for group
    templist = [];
    for ii = 1:length(tempgroup), %for each groupsubj index selected
        if ii<length(tempgroup),
            templist = [templist handles.group_subj{tempgroup(ii)} ',']; %create templist for group string
        else  templist = [templist handles.group_subj{tempgroup(ii)}];
        end
    end
    templist = cellstr(templist);
    if isempty(handles.groups),
        handles.groups{1} = {handles.group_subj{tempgroup}};
        set(handles.groupdefined_list,'String', templist);
        set(handles.visgroups_list,'String',templist);
    else
        handles.groups{end+1} = {handles.group_subj{tempgroup}}; %update groups
        initiallist = cellstr(get(handles.groupdefined_list,'String'));
        newlist = [initiallist; templist];
        set(handles.groupdefined_list,'String',newlist);
        set(handles.visgroups_list,'String',newlist);
    end
    %Update subject list at the end
    handles.group_subj(tempgroup)=[]; %remove subjects from list
    set(handles.groupsubj_list,'String',handles.group_subj); %update available subjects
    set(handles.groupsubj_list,'Value',1); %update list selection
end
guidata(hObject, handles);

function groupremove_push_Callback(hObject, eventdata, handles)
if ~isempty(handles.groups),
    deletegroup = get(handles.groupdefined_list,'Value'); %index of group to delete
    oldlist = get(handles.groupdefined_list,'String'); %list of groups
    oldlist(deletegroup)=[]; %removes group from list
    if ~isempty(oldlist), %if list has remaining groups
        set(handles.groupdefined_list,'String',oldlist); %update list
        set(handles.visgroups_list,'String',oldlist);
    else set(handles.groupdefined_list,'String','No groups defined'); %default line
        set(handles.visgroups_list,'String','No groups defined');
    end
    handles.group_subj = [handles.group_subj handles.groups{deletegroup}]; %put subjects back into available list
    handles.groups(deletegroup) = []; %delete group from master list
    set(handles.groupsubj_list,'String',handles.group_subj);
    set(handles.groupdefined_list,'Value',1);
end
guidata(hObject, handles);

function zscoreon_toggle_Callback(hObject, eventdata, handles)
if get(handles.zscoreoff_toggle,'Value')==1,
    set(handles.zscoreoff_toggle,'Value',0);
    handles=preprocess(handles);
    handles = update_EDAdisplay(handles);
else set(handles.zscoreon_toggle,'Value',1); %do nothing
end
guidata(hObject, handles);

function zscoreoff_toggle_Callback(hObject, eventdata, handles)
if get(handles.zscoreon_toggle,'Value')==1,
    set(handles.zscoreon_toggle,'Value',0);
    handles=preprocess(handles);
    handles = update_EDAdisplay(handles);
else set(handles.zscoreoff_toggle,'Value',1); %do nothing
end
guidata(hObject, handles);

function handles = preprocess(handles)
%% RESET TIME
%Construct master timeline
EDAtime=(min([handles.data.EDAon]):(1/handles.parameters.EDAsamplerate):max([handles.data.EDAoff]))';
ACCtime=(min([handles.data.ACCon]):(1/handles.parameters.ACCsamplerate):max([handles.data.ACCoff]))';
%Temporal alignment / pad with NaN
for ii = 1:length(handles.data),
    [inx, ~] = find(EDAtime==handles.data(ii).EDAon);
    handles.data(ii).EDA(1:length(EDAtime),1)=NaN;
    handles.data(ii).EDA(inx:inx+length(handles.data(ii).urEDA(3:end-1)),1)=handles.data(ii).urEDA(3:end);
    if ~isempty(handles.data(ii).urACC),
        [inx, ~] = find(ACCtime==handles.data(ii).ACCon);
        handles.data(ii).ACC(1:length(ACCtime),1:3)=NaN;
        handles.data(ii).ACC(inx:inx+length(handles.data(ii).urACC(3:end-1,1)),:)=handles.data(ii).urACC(3:end,:);
    end
end

handles.time.EDAtimestamps=EDAtime;
handles.time.EDAsecs=EDAtime-EDAtime(1);
handles.time.EDAmins=handles.time.EDAsecs./60;
handles.time.ACCtimestamps=ACCtime;
handles.time.ACCsecs=ACCtime-ACCtime(1);
handles.time.ACCmins=handles.time.ACCsecs./60;

%Compute summed accelerometer signals
for ii=1:length(handles.data),
    if ~isempty(handles.data(ii).urACC),
        handles.data(ii).ACC=sum(abs(handles.data(ii).ACC),2);
    end
end

%Copy original data to current data handles
handles.EDA=[handles.data.EDA];
handles.ACC=[handles.data.ACC];

%% Check signal processing steps
if get(handles.zero_checkbox,'Value')==1
    removezero=1;
else removezero=0;
end
if get(handles.z_checkbox,'Value')==1
    zscore=1;
else zscore=0;
end
if get(handles.detrend_checkbox,'Value')==1
    lindetrend=1;
else lindetrend=0;
end
if get(handles.tempsmooth_checkbox,'Value')==1
    smoothing=str2double(get(handles.tempsmooth_edit,'String'));
else smoothing=0;
end
if get(handles.filter_checkbox,'Value')==1
    filters(1)=str2double(get(handles.filtorder_edit,'String'));
    filters(2)=str2double(get(handles.filthigh_edit,'String'));
    filters(3)=str2double(get(handles.filtlow_edit,'String'));
else filters=[0 0 0];
end
if get(handles.croptime_checkbox,'Value')==1
    croptimestart=str2double(get(handles.crop1_edit,'String'));
    croptimeend=str2double(get(handles.crop2_edit,'String'));
else croptimestart=0; croptimeend=0;
end
[handles.EDA,handles.ACC,handles.time] = sochro_preprocess(handles.EDA,handles.ACC,handles.time,handles.parameters.EDAsamplerate,zscore,filters,smoothing,lindetrend,removezero,croptimestart,croptimeend);

function filtorder_edit_Callback(hObject, eventdata, handles)
    handles = preprocess(handles);
    handles = update_EDAdisplay(handles);
    guidata(hObject, handles);

function filtorder_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function filthigh_edit_Callback(hObject, eventdata, handles)
    handles = preprocess(handles);
    handles = update_EDAdisplay(handles);
    guidata(hObject, handles);

function filthigh_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function filtlow_edit_Callback(hObject, eventdata, handles)
    handles = preprocess(handles);
    handles = update_EDAdisplay(handles);
    guidata(hObject, handles);

function filtlow_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function accpeaks_push_Callback(hObject, eventdata, handles)
handles = accpeaks(handles);
guidata(hObject,handles);

function handles = accpeaks(handles)
%Sum and normalise
handles.display.data = nansum([handles.ACC],2);
handles.display.data = ((handles.display.data-min(handles.display.data))/(max(handles.display.data)-min(handles.display.data))).^2;
handles.display.time = handles.time.ACCmins;
%Find spikes
rule=1;
ACCthresh = str2double(get(handles.accthresh_edit,'String'));
while rule==1, %continue until rule changes
    ts = find(handles.display.data>ACCthresh);
    if isempty(ts), %if no peaks found
        disp('No thresholded accelerometer activity peaks found - threshold lowered');
        ACCthresh=ACCthresh-.1;
        set(handles.accthresh_edit,'String',num2str(ACCthresh));
    else rule=0;
    end
end
%Define peaks
inx=[]; TS=[];
if length(ts)>1,
    for ii = 1:length(ts),
        if ii==1, inx(1)=ts(ii); %if first ts
        else
            if ts(ii)-ts(ii-1)==1, inx(end+1)=ts(ii); %add ts to list if continued peak
            else
                [~,ix]=max(handles.display.data(inx)); %find peak within list of ts
                TS(end+1)=inx(ix);
                inx=ts(ii);
            end
        end
    end
    if ~isempty(inx),
        [~,ix]=max(handles.display.data(inx));
        TS(end+1)=inx(ix);
    end
else TS=ts;
end
%Plot data
set(handles.edaaxes,'NextPlot','replace');
plot(handles.edaaxes,handles.display.time,handles.display.data, ...
    'LineStyle','-', ...
    'Color',handles.plotline_colours(1,:));
set(handles.edaaxes,'NextPlot','add');

%[handles.display.time(1) handles.display.time(end)]
set(handles.edaaxes,'XLim',[handles.display.time(1) handles.display.time(end) ]);
set(handles.edaaxes,'YLim',[0 1]);
set(get(handles.edaaxes,'XLabel'),'String','Time (minutes)');
set(get(handles.edaaxes,'YLabel'),'String','Summed Accelerometer Activity');
set(handles.accpeak_text,'String',['of ' num2str(length(TS)) ' peaks']);
%% Insert peak lines (TS is in datapoints)
for ii = 1:length(TS),
    line( [handles.display.time(TS(ii)) handles.display.time(TS(ii))], [0 1], 'Color', [1 0 0]);
end
currentpeak = str2double(get(handles.accpeak_edit,'String'));
if currentpeak>length(TS),
    set(handles.accpeak_edit,'String','1');
    TSx = TS(1);
else
    TSx = TS(currentpeak);
end
line( [handles.display.time(TSx) handles.display.time(TSx)], [0 1], 'Color', [0 1 0]);
handles.parameters.accpeak = handles.display.time(TSx);
handles.parameters.peakn = length(TS);
set(handles.accpeakamp_text,'String',['Amplitude: ' num2str(handles.display.data(TSx))]);
set(handles.accpeakmins_text,'String',['Peak at ' num2str(handles.display.time(TSx)) ' minutes']);

function accsync_push_Callback(hObject, eventdata, handles)
delay = str2double(get(handles.acctime_edit,'String')) - handles.parameters.accpeak;
% Create shifted times
handles.time.EDAsecs = handles.time.EDAsecs + (delay*60);
handles.time.EDAmins = handles.time.EDAmins + delay;
handles.time.ACCsecs = handles.time.ACCsecs + (delay*60);
handles.time.ACCmins = handles.time.ACCmins + delay;
handles = accpeaks(handles);
guidata(hObject, handles);

function accreset_push_Callback(hObject, eventdata, handles)
EDAtime=(min([handles.data.EDAon]):(1/handles.parameters.EDAsamplerate):max([handles.data.EDAoff]))';
ACCtime=(min([handles.data.ACCon]):(1/handles.parameters.ACCsamplerate):max([handles.data.ACCoff]))';
handles.time.EDAtimestamps=EDAtime;
handles.time.EDAsecs=EDAtime-EDAtime(1);
handles.time.EDAmins=handles.time.EDAsecs./60;
handles.time.ACCtimestamps=ACCtime;
handles.time.ACCsecs=ACCtime-ACCtime(1);
handles.time.ACCmins=handles.time.ACCsecs./60;
handles = preprocess(handles);
handles = accpeaks(handles);
guidata(hObject,handles);

function accpeak_edit_Callback(hObject, eventdata, handles)
current = round(str2double(get(handles.accpeak_edit,'String')));
if current<1 || current>handles.parameters.peakn,
    set(handles.accpeak_edit,'String','1');
end
handles = accpeaks(handles);
guidata(hObject,handles);

function accpeak_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function accleft_push_Callback(hObject, eventdata, handles)
current = str2double(get(handles.accpeak_edit,'String'));
if current==1,
    return
else set(handles.accpeak_edit,'String',num2str(current-1));
end
handles = accpeaks(handles);
guidata(hObject,handles);

function accright_push_Callback(hObject, eventdata, handles)
current = str2double(get(handles.accpeak_edit,'String'));
if current==handles.parameters.peakn,
    return
else set(handles.accpeak_edit,'String',num2str(current+1));
end
handles = accpeaks(handles);
guidata(hObject,handles);

function accthresh_edit_Callback(hObject, eventdata, handles)

function accthresh_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function acctime_edit_Callback(hObject, eventdata, handles)

function acctime_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function cond_list_Callback(hObject, eventdata, handles)
handles = update_EDAdisplay(handles);
guidata(hObject,handles);

function cond_list_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function condname_edit_Callback(hObject, eventdata, handles)
if isempty(get(handles.condname_edit,'String'));
    set(handles.condname_edit,'String','Enter condition name...');
else
    check=0;
    for ii=1:length(handles.conditions),
        if strcmpi(get(handles.condname_edit,'String'),handles.conditions(ii).name);
            set(handles.condcreate_push,'String','Add to condition');
            check=1;
        end
    end
    if check==0,
        set(handles.condcreate_push,'String','Create condition');
    end
end
guidata(hObject,handles);

function condname_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function condstart_edit_Callback(hObject, eventdata, handles)
handles = update_EDAdisplay(handles);
guidata(hObject,handles);

function condstart_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function condend_edit_Callback(hObject, eventdata, handles)
handles = update_EDAdisplay(handles);
guidata(hObject,handles);

function condend_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function condcreate_push_Callback(hObject, eventdata, handles)
t1 = str2double(get(handles.condstart_edit,'String'));
t2 = str2double(get(handles.condend_edit,'String'));
%Check to see if t1 and t2 are in range
if t1<handles.time.EDAmins(1) || t2>handles.time.EDAmins(end),
    errordlg('Times entered are outside of data range');
    return
elseif t1>t2,
    errordlg('Beginning of time window must come before end of time window');
    return
end
[c1, inx1] = min(abs(handles.time.EDAmins-t1));
[c2, inx2] = min(abs(handles.time.EDAmins-t2));
c1 = handles.time.EDAmins(inx1); %convert c to time
c2 = handles.time.EDAmins(inx2); %convert c to time
if (inx2-inx1)<2,
    errordlg('Time window is too narrow');
    return
end
condname=get(handles.condname_edit,'String');
if strcmp(condname,'Enter condition name...');
    errordlg('Please enter a name for the current condition');
    return
else %Check if condition already exists
    check=0;
    for ii = 1:length(handles.conditions),
        if strcmpi(condname,handles.conditions(ii).name)==1, %if preexisting condition
            handles.conditions(ii).windows{end+1}={[c1 c2]}; %assign times
            handles.conditions(ii).windowsinx{end+1}={[inx1 inx2]}; %assign indices (based on handles.time.EDAmins)
            handles.conditions(ii).string = [condname ': ' num2str(length(handles.conditions(ii).windows)) ' time windows'];
            check=1;
        end
    end
    if check==0, %if new condition
        handles.conditions(end+1).name = condname; %only add one to end in first line
        handles.conditions(end).windows{1}={[c1 c2]}; %assign times
        handles.conditions(end).windowsinx{1}={[inx1 inx2]}; %assign indices (based on handles.data.times)
        handles.conditions(end).string = [condname ': 1 time window'];
    end
end
set(handles.condstart_edit,'String',num2str(0));
set(handles.condend_edit,'String',num2str(0));
set(handles.cond_list,'String',{handles.conditions.string});
set(handles.condcreate_push,'String','Add to condition');
handles = update_EDAdisplay(handles);
guidata(hObject,handles);

function condremove_push_Callback(hObject, eventdata, handles)
selectedcondition = get(handles.cond_list,'Value');
handles.conditions(selectedcondition)=[];
set(handles.condstart_edit,'String',num2str(0));
set(handles.condend_edit,'String',num2str(0));
set(handles.cond_list,'String',{handles.conditions.string});
set(handles.cond_list,'Value',1);
%if no conditions left
if isempty(handles.conditions),
    set(handles.cond_list,'String','No conditions defined');
    set(handles.conddisp_toggle,'Value',0);
end
handles = update_EDAdisplay(handles);
guidata(hObject,handles);

function conddisp_toggle_Callback(hObject, eventdata, handles)
handles = update_EDAdisplay(handles);
guidata(hObject,handles);


% --- Executes on button press in connectivity_push.
function connectivity_push_Callback(hObject, eventdata, handles)
%SORT BY SUBJECT
if ~isfield(handles.parameters,'grouporder'), %If not done, sort subjects by group specification
    %Assign subjects to groups if they have no specification
    gg = length(handles.groups)+1;
    for ii = 1:length(handles.group_subj),
        inx = find(strcmp(handles.group_subj{ii},{handles.data.subj}));
        handles.data(inx).group=gg; %assign group to handles.data
        handles.groups{gg}{1}=handles.group_subj{ii};
        gg=gg+1;
    end
    for gg = 1:length(handles.groups), %for each group
        for ii = 1:length(handles.groups{gg}), %for each subject in group gg
            inx = find(strcmp(handles.groups{gg}{ii},{handles.data.subj}));
            handles.data(inx).group=gg; %assign group to handles.data
        end
    end
    %sort data by group
    for ii = 1:length(handles.data), %copy processed EDA to handles.data
        handles.data(ii).EDA = handles.EDA(:,ii);
        handles.data(ii).ACC = handles.ACC(:,ii);
    end
    fields=fieldnames(handles.data);
    cell=struct2cell(handles.data);
    sz =size(cell);
    cell=reshape(cell,sz(1),[]);
    cell=cell';
    inx = find(strcmp('group',fields));
    cell=sortrows(cell,inx);
    cell=reshape(cell',sz);
    handles.data=cell2struct(cell,fields,1);
    handles.EDA = [handles.data.EDA];
    handles.ACC = [handles.data.ACC];
    handles.parameters.grouporder = [handles.data.group];
end

if isempty(handles.conditions),
    errordlg('Conditions must be defined before running connectivity analysis');
    return
end
%% CREATE SURROGATES
if isempty(handles.surrogate), %if no surrogate data exists
    handles = corsurrogate(handles);
end
%% CALCULATE CONDITION THRESHOLD
handles.surrogate.statcond = []; %clear thresholds
Y = [];
for pp=1:1000, %for each permutation
    set(handles.status_edit,'String',['Calculating condition significance threshold: Correlating surrogate dataset ' num2str(pp) '/1000']);
    pause(0.001);
    for kk = 1:length(handles.conditions),
        win = [];
        for ww = 1:length(handles.conditions(kk).windowsinx), %for each time window
            win = [win handles.conditions(kk).windowsinx{ww}{1}(1):handles.conditions(kk).windowsinx{ww}{1}(2)];
        end
        X = tril(corr(handles.surrogate.data(win,:,pp)),-1); %concatonate and correlated
        Y = [Y X(X~=0)']; %for std calculation
        handles.surrogate.statcond.R(pp,kk)=mean(X(X~=0));
    end
end
for kk = 1:length(handles.conditions),
    handles.surrogate.statcond.mean(kk) = mean(handles.surrogate.statcond.R(:,kk));
    handles.surrogate.statcond.meanvar(kk) = std(handles.surrogate.statcond.R(:,kk));
    handles.surrogate.statcond.std(kk) = std(Y);
    handles.surrogate.statcond.thresh(kk) = (1.96*handles.surrogate.statcond.std(kk))+handles.surrogate.statcond.mean(kk);
end

%% CALCULATE TIME-LAG CORRECTED CORRELATION
handles.parameters.lag = 0;
if get(handles.corlag_toggle,'Value')==1, %if time-lag correlation
    handles.parameters.lag = 1;
    lagsamples=handles.parameters.EDAsamplerate*str2double(get(handles.corlag_edit,'String'));
    for kk = 1:length(handles.conditions), %for each condition
        set(handles.status_edit,'String',['Calculating time-lagged correlations: Condition "' handles.conditions(kk).name '"'  ]); pause(0.001);
        for ii = 1:length(handles.data), %for each participant
            for jj = 1:length(handles.data), %for each other participant
                if ii==jj, %if autocorrelate
                    handles.corr.static(kk).R(ii,jj)=NaN;
                else
                    wincorr = [];
                    wincorrlag = [];
                    X = []; Y = [];
                    for ww = 1:length(handles.conditions(kk).windowsinx), %for each time window
                        win1 = handles.conditions(kk).windowsinx{ww}{1}(1):handles.conditions(kk).windowsinx{ww}{1}(2);
                        [acor,lag]=xcorr(handles.EDA(win1,ii),handles.EDA(win1,jj),lagsamples,'coeff');
                        [~,I]=max(acor);
                        win2 = win1 - lag(I);
                        X = [X handles.EDA(win1,ii)];
                        Y = [Y handles.EDA(win2,jj)];
                    end
                    handles.corr.static(kk).R(ii,jj)=corr(handles.EDA(X,Y));
                end
            end
        end
    end
else
    %% CALCULATE CORRELATION (NO LAG)
    for kk = 1:length(handles.conditions), %for each condition
        set(handles.status_edit,'String',['Calculating correlations: Condition "' handles.conditions(kk).name '"'  ]); pause(0.001);
        win = [];
        for ww = 1:length(handles.conditions(kk).windowsinx), %for each time window
            win = [win handles.conditions(kk).windowsinx{ww}{1}(1):handles.conditions(kk).windowsinx{ww}{1}(2)];
        end
        handles.corr.static(kk).R=corr(handles.EDA(win,:));
    end
end

%% Calculate Z threshold
X = tril(corr(handles.EDA));
X=X(X~=0); X=X(X~=1);
handles.surrogate.statwhole.mean = mean(X);
handles.surrogate.statwhole.std = std(X);
handles.surrogate.statwhole.thresh = (1.96*handles.surrogate.statwhole.std)+handles.surrogate.statwhole.mean;

%% Update status
readout = {};
readout{end+1} = 'STATIC CONNECTIVITY NULL:';
for kk = 1:length(handles.conditions),
    readout{end+1} = ['Condition: ' handles.conditions(kk).name];
    readout{end+1} = ['Null mean: ' num2str(handles.surrogate.statcond.mean(kk))];
    readout{end+1} = ['Null SD: ' num2str(handles.surrogate.statcond.std(kk))];
    readout{end+1} = ['Surrogate variance: ' num2str(handles.surrogate.statcond.meanvar(kk))];
    readout{end+1} = ['Threshold R value: ' num2str(handles.surrogate.statcond.thresh(kk))];
    readout{end+1} = ['Data mean value: ' num2str(handles.surrogate.statwhole.mean)];
    readout{end+1} = ['Data SD: ' num2str(handles.surrogate.statwhole.std)];
    readout{end+1} = ['Z threshold: ' num2str(handles.surrogate.statwhole.thresh)];
end
set(handles.corstatus_list,'String',readout);


%% CALCULATE GRAPH METRICS PER CONDITION
for kk = 1:length(handles.corr.static), %for each condition
    set(handles.status_edit,'String',['Calculating graph metrics: Condition "' handles.conditions(kk).name '"'  ]); pause(0.001);
    %covariance
    handles.corr.static(kk).Rt = handles.corr.static(kk).R; %copy raw covariance
    handles.corr.static(kk).Rt(handles.corr.static(kk).R<handles.surrogate.statcond.thresh(kk))=0; %condition specific threshold
    handles.corr.static(kk).Rz = handles.corr.static(kk).R; %copy raw covariance
    handles.corr.static(kk).Rz(handles.corr.static(kk).R<handles.surrogate.statwhole.thresh)=0; %condition specific threshold
    %Subject level
    %mean connectivity, connectivity strength, and degree measures
    X=handles.corr.static(kk).R; %copy covariance matrix
    X(X==0)=NaN; X(X==1)=NaN; %make zeros and ones into NaN
    handles.corr.static(kk).FC=nanmean(X); %no threshold
    handles.corr.static(kk).FS = nansum(X); %no threshold
    X=handles.corr.static(kk).Rt;
    X(X==0)=NaN; X(X==1)=NaN;
    handles.corr.static(kk).FCt=nanmean(X); %threshold
    handles.corr.static(kk).FSt = nansum(X); %threshold
    handles.corr.static(kk).NDt=nansum(ceil(X)); %threshold
    handles.corr.static(kk).NPt=nansum(ceil(X))./(length(X(:,1))-1);
    X=handles.corr.static(kk).Rz;
    X(X==0)=NaN; X(X==1)=NaN;
    handles.corr.static(kk).FCz=nanmean(X); %threshold
    handles.corr.static(kk).FSz = nansum(X); %threshold
    handles.corr.static(kk).NDz=nansum(ceil(X)); %threshold
    handles.corr.static(kk).NPz=nansum(ceil(X))./(length(X(:,1))-1);
    %Group Dynamics level
    %Mean connectivity / Degree / Proportion connectivity per participant
    for ii = 1:length(handles.corr.static(kk).R), %for each participant
        inx = find([handles.data.group]==handles.data(ii).group); %find ingroup members
        iinx = setdiff(1:length(handles.corr.static(kk).R),inx); %remove ingroup members to make indices of outgroup members
        inx(inx==ii)=[]; %remove current subject (autocorrelate)
        handles.corr.static(kk).FC_within(ii) = nanmean(handles.corr.static(kk).R(inx,ii)); %within group connectivity
        handles.corr.static(kk).FC_between(ii) = nanmean(handles.corr.static(kk).R(iinx,ii)); %between group connectivity
        handles.corr.static(kk).FS_within(ii) = nansum(handles.corr.static(kk).R(inx,ii)); %within group connectivity
        handles.corr.static(kk).FS_between(ii) = nansum(handles.corr.static(kk).R(iinx,ii)); %between group connectivity
        
        handles.corr.static(kk).FCt_within(ii) = nanmean(handles.corr.static(kk).Rt(inx,ii)); %within group connectivity
        handles.corr.static(kk).NDt_within(ii) = nansum(ceil(handles.corr.static(kk).Rt(inx,ii))); %within group degree
        handles.corr.static(kk).NPt_within(ii) = nansum(ceil(handles.corr.static(kk).Rt(inx,ii)))/length(inx); %within group proportion
        handles.corr.static(kk).FCt_between(ii) = nanmean(handles.corr.static(kk).Rt(iinx,ii)); %between group connectivity
        handles.corr.static(kk).NDt_between(ii) = nansum(ceil(handles.corr.static(kk).Rt(iinx,ii))); %between group degree
        handles.corr.static(kk).NPt_between(ii) = nansum(ceil(handles.corr.static(kk).Rt(iinx,ii)))/length(iinx); %between group proportion
        handles.corr.static(kk).FSt_within(ii) = nansum(handles.corr.static(kk).Rt(inx,ii)); %within group connectivity
        handles.corr.static(kk).FSt_between(ii) = nansum(handles.corr.static(kk).Rt(iinx,ii)); %between group connectivity
        
        handles.corr.static(kk).FCz_within(ii) = nanmean(handles.corr.static(kk).Rz(inx,ii)); %within group connectivity
        handles.corr.static(kk).NDz_within(ii) = nansum(ceil(handles.corr.static(kk).Rz(inx,ii))); %within group degree
        handles.corr.static(kk).NPz_within(ii) = nansum(ceil(handles.corr.static(kk).Rz(inx,ii)))/length(inx); %within group proportion
        handles.corr.static(kk).FCz_between(ii) = nanmean(handles.corr.static(kk).Rz(iinx,ii)); %between group connectivity
        handles.corr.static(kk).NDz_between(ii) = nansum(ceil(handles.corr.static(kk).Rz(iinx,ii))); %between group degree
        handles.corr.static(kk).NPz_between(ii) = nansum(ceil(handles.corr.static(kk).Rz(iinx,ii)))/length(iinx); %between group proportion
        handles.corr.static(kk).FSz_within(ii) = nansum(handles.corr.static(kk).Rz(inx,ii)); %within group connectivity
        handles.corr.static(kk).FSz_between(ii) = nansum(handles.corr.static(kk).Rz(iinx,ii)); %between group connectivity
    end
    %Group Comparison
    for gg = 1:length(handles.groups),
        inx = find([handles.data.group]==gg); %find group members
        iinx = setdiff(1:length(handles.corr.static(kk).R),inx);
        X = handles.corr.static(kk).R; %RAW R
        X(X==1)=NaN; %remove autocorrelations
        Y = X(inx,inx);
        Z = X(iinx,inx);
        handles.corr.static(kk).FC_gwithin(gg) = nanmean(nanmean(Y)); %within group mean connectivity
        handles.corr.static(kk).FC_gbetween(gg) = nanmean(nanmean(Z)); %between group mean connectivity
        handles.corr.static(kk).FS_gwithin(gg) = nansum(nansum(Y)); %within group mean connectivity
        handles.corr.static(kk).FS_gbetween(gg) = nansum(nansum(Z)); %between group mean connectivity
        
        X = handles.corr.static(kk).Rt; %THRESHOLD R
        X(X==1)=NaN; %remove autocorrelations
        Y = X(inx,inx);
        Z = X(iinx,inx);
        handles.corr.static(kk).FCt_gwithin(gg) = nanmean(nanmean(Y)); %within group mean connectivity 
        handles.corr.static(kk).NDt_gwithin(gg) = mean(nansum(ceil(Y))); %within group mean degree 
        handles.corr.static(kk).NPt_gwithin(gg) = nansum(nansum(ceil(Y))) / (length(inx)*(length(inx)-1)); %within group proportion of possible degree 
        handles.corr.static(kk).FCt_gbetween(gg) = nanmean(nanmean(Z)); %between group mean connectivity 
        handles.corr.static(kk).NDt_gbetween(gg) = mean(nansum(ceil(Z))); %between group mean degree 
        handles.corr.static(kk).NPt_gbetween(gg) = nansum(nansum(ceil(Z))) / (length(iinx)*(length(inx))); %between group proportion of possible degree 
        handles.corr.static(kk).FSt_gwithin(gg) = nansum(nansum(Y)); %within group mean connectivity
        handles.corr.static(kk).FSt_gbetween(gg) = nansum(nansum(Z)); %between group mean connectivity
        
        X = handles.corr.static(kk).Rz; %THRESHOLD R
        X(X==1)=NaN; %remove autocorrelations
        Y = X(inx,inx);
        Z = X(iinx,inx);
        handles.corr.static(kk).FCz_gwithin(gg) = nanmean(nanmean(Y)); %within group mean connectivity (whole thresh)
        handles.corr.static(kk).NDz_gwithin(gg) = mean(nansum(ceil(Y))); %within group mean degree (whole thresh)
        handles.corr.static(kk).NPz_gwithin(gg) = nansum(nansum(ceil(Y))) / (length(inx)*(length(inx)-1)); %within group proportion of possible degree (whole thresh)
        handles.corr.static(kk).FCz_gbetween(gg) = nanmean(nanmean(Z)); %between group mean connectivity (whole thresh)
        handles.corr.static(kk).NDz_gbetween(gg) = mean(nansum(ceil(Z))); %between group mean degree (whole thresh)
        handles.corr.static(kk).NPz_gbetween(gg) = nansum(nansum(ceil(Z))) / (length(iinx)*(length(inx))); %between group proportion of possible degree (whole thresh)
        handles.corr.static(kk).FSz_gwithin(gg) = nansum(nansum(Y)); %within group mean connectivity
        handles.corr.static(kk).FSz_gbetween(gg) = nansum(nansum(Z)); %between group mean connectivity
    end
    %Graph level
    handles.corr.static(kk).GC= nanmean(handles.corr.static(kk).FC); %mean connectivity
    handles.corr.static(kk).GCt = nanmean(handles.corr.static(kk).FCt); %mean connectivity
    handles.corr.static(kk).GCz = nanmean(handles.corr.static(kk).FCz); %mean connectivity
    handles.corr.static(kk).GS= nanmean(handles.corr.static(kk).FS); %mean connectivity strength
    handles.corr.static(kk).GSt = nanmean(handles.corr.static(kk).FSt); %mean connectivity strength
    handles.corr.static(kk).GDt = nanmean(handles.corr.static(kk).NDt); %mean degree
    handles.corr.static(kk).GDt_within = nanmean(handles.corr.static(kk).NDt_within);
    handles.corr.static(kk).GDt_between = nanmean(handles.corr.static(kk).NDt_between);
    handles.corr.static(kk).GPt = nanmean(handles.corr.static(kk).NPt);
    handles.corr.static(kk).GPt_within = nanmean(handles.corr.static(kk).NPt_within);
    handles.corr.static(kk).GPt_between = nanmean(handles.corr.static(kk).NPt_between);
    
    handles.corr.static(kk).GSz = nanmean(handles.corr.static(kk).FSz); %mean connectivity strength
    handles.corr.static(kk).GDz = nanmean(handles.corr.static(kk).NDz); %mean degree
    handles.corr.static(kk).GDz_within = nanmean(handles.corr.static(kk).NDz_within);
    handles.corr.static(kk).GDz_between = nanmean(handles.corr.static(kk).NDz_between);
    handles.corr.static(kk).GPz = nanmean(handles.corr.static(kk).NPz);
    handles.corr.static(kk).GPz_within = nanmean(handles.corr.static(kk).NPz_within);
    handles.corr.static(kk).GPz_between = nanmean(handles.corr.static(kk).NPz_between);
    
    handles.corr.static(kk).GC_within = nanmean(handles.corr.static(kk).FC_within); %within group connectivity
    handles.corr.static(kk).GC_between = nanmean(handles.corr.static(kk).FC_between); %between group connectivity
    handles.corr.static(kk).GCt_within = nanmean(handles.corr.static(kk).FCt_within); %within group connectivity
    handles.corr.static(kk).GCt_between = nanmean(handles.corr.static(kk).FCt_between); %between group connectivity
    handles.corr.static(kk).GS_within = nanmean(handles.corr.static(kk).FS_within); %within group connectivity strength
    handles.corr.static(kk).GS_between = nanmean(handles.corr.static(kk).FS_between); %between group connectivity strength
    handles.corr.static(kk).GSt_within = nanmean(handles.corr.static(kk).FSt_within); %within group connectivity strength
    handles.corr.static(kk).GSt_between = nanmean(handles.corr.static(kk).FSt_between); %between group connectivity strength
    
    handles.corr.static(kk).GCz_within = nanmean(handles.corr.static(kk).FCz_within); %within group connectivity
    handles.corr.static(kk).GCz_between = nanmean(handles.corr.static(kk).FCz_between); %between group connectivity
    handles.corr.static(kk).GSz_within = nanmean(handles.corr.static(kk).FSz_within); %within group connectivity strength
    handles.corr.static(kk).GSz_between = nanmean(handles.corr.static(kk).FSz_between); %between group connectivity strength
    
end

%% Finalise and setup strings/update GUI
set(handles.status_edit,'String','Static connectivity analysis complete'); pause(0.001);
set(handles.corstat_text,'String','Complete'); %change text to complete
string(1).t = 'Static';
if isfield(handles.corr,'dynamic'),
    string(end+1).t = 'Dynamic';
end
set(handles.corconn_pop,'String',{'Connectivity: ' string.t});
handles.parameters.corr=1;
if isfield(handles,'edaaxes'),
    set(handles.edaaxes,'Visible','off');
    axes(handles.edaaxes); cla;
end
guidata(hObject,handles)

function cordyn_push_Callback(hObject, eventdata, handles)
%SORT BY SUBJECT
if ~isfield(handles.parameters,'grouporder'), %If not done, sort subjects by group specification
    set(handles.status_edit,'String','Sorting subjects by group specification'); pause(0.001);
    %Assign subjects to groups if they have no specification
    gg = length(handles.groups)+1;
    for ii = 1:length(handles.group_subj),
        inx = find(strcmp(handles.group_subj{ii},{handles.data.subj}));
        handles.data(inx).group=gg; %assign group to handles.data
        handles.groups{gg}{1}=handles.group_subj{ii};
        gg=gg+1;
    end
    for gg = 1:length(handles.groups), %for each group
        for ii = 1:length(handles.groups{gg}), %for each subject in group gg
            inx = find(strcmp(handles.groups{gg}{ii},{handles.data.subj}));
            handles.data(inx).group=gg; %assign group to handles.data
        end
    end
    %sort data by group
    for ii = 1:length(handles.data),
        handles.data(ii).EDA = handles.EDA(:,ii);
        handles.data(ii).ACC = handles.ACC(:,ii);
    end
    fields=fieldnames(handles.data);
    cell=struct2cell(handles.data);
    sz =size(cell);
    cell=reshape(cell,sz(1),[]);
    cell=cell';
    inx = find(strcmp('group',fields));
    cell=sortrows(cell,inx);
    cell=reshape(cell',sz);
    handles.data=cell2struct(cell,fields,1);
    handles.EDA = [handles.data.EDA];
    handles.ACC = [handles.data.ACC];
    handles.parameters.grouporder = [handles.data.group];
end

if isempty(handles.conditions),
    errordlg('Conditions must be defined before running connectivity analysis');
    return
end
%% CREATE SURROGATES
if isempty(handles.surrogate), %if no surrogate data exists
    handles = corsurrogate(handles);
end

%Sliding window parameters and timelimes
winlength = handles.parameters.EDAsamplerate*str2double(get(handles.corwin_edit,'String')); 
winslide = handles.parameters.EDAsamplerate*str2double(get(handles.corslide_edit,'String'));
handles.time.DYNmins = handles.time.EDAmins(1):((winslide/handles.parameters.EDAsamplerate)/60):handles.time.EDAmins(length(handles.time.EDAmins)-winlength); %whole session window inx

% Create indices for dynamic timeseries per condition
for kk = 1:length(handles.conditions), %for each condition kk
    for ww = 1:length(handles.conditions(kk).windows) %for each window for condition kk
        handles.conditions(kk).DYNinx{ww}=handles.conditions(kk).windowsinx{ww}{1}(1):winslide:handles.conditions(kk).windowsinx{ww}{1}(2); %condition window inx
        handles.conditions(kk).DYNmins{ww}=handles.conditions(kk).windowsinx{ww}{1}(1):((winslide/handles.parameters.EDAsamplerate)/60):handles.conditions(kk).windowsinx{ww}{1}(2); %condition window inx
    end
end

%% SURROGATE DYNAMIC CONNECTIVITY
if ~isfield(handles.surrogate,'dynwhole'), %if correlations not run
    handles.surrogate.dyndata = zeros(length(handles.data),length(handles.data),length(1:winslide:length(handles.surrogate.data(:,1,1))-(winlength)),1000);
    for pp = 1:1000,
        set(handles.status_edit,'String',['Calculating significance threshold: Correlating surrogate dataset ' num2str(pp) '/1000']); pause(0.001);
        cc=1;
        for ww = 1:winslide:length(handles.surrogate.data(:,1,1))-(winlength);
            handles.surrogate.dyndata(:,:,cc,pp)=tril(corr(handles.surrogate.data(ww:ww+winlength,:,pp)),-1);
            cc=cc+1;
        end
    end
    %%  CALCULATE WHOLE SESSION THRESHOLD
    Y = reshape(handles.surrogate.dyndata,[1,length(handles.data)*length(handles.data)*(cc-1)*pp]);
    Y(Y==0)=NaN; Y(Y==1)=NaN;
    handles.surrogate.dynwhole.mean = nanmean(Y);
    handles.surrogate.dynwhole.std = nanstd(Y);
    handles.surrogate.dynwhole.thresh = (1.96*handles.surrogate.dynwhole.std)+handles.surrogate.dynwhole.mean;
    
    %% CALCULATE CONDITION THRESHOLD
    % Get condition specific correlations
    for kk = 1:length(handles.conditions),
        wininx = [handles.conditions(kk).DYNinx{:}];
        X=zeros(length(handles.data),length(handles.data),length(wininx),1000);
        for pp = 1:1000,
            cc=1;
            for ww = length(wininx);
                set(handles.status_edit,'String',['Calculating condition specific significance threshold: surrogate set ' num2str(pp) '/1000 for condition ' num2str(kk)]); pause(0.001);
                if wininx(ww)+winlength<length(handles.surrogate.data(:,1,pp))
                X(:,:,cc,pp) = tril(corr(handles.surrogate.data(wininx(ww):wininx(ww)+winlength,:,pp)));
                cc=cc+1;
                end
            end
        end
        Y = reshape(X,[1,length(handles.data)*length(handles.data)*length(wininx)*1000]);
        Y(Y==0)=NaN; Y(Y==1)=NaN;
        handles.surrogate.dyncond.mean(kk)=nanmean(Y);
        handles.surrogate.dyncond.std(kk)=nanstd(Y);
        handles.surrogate.dyncond.thresh(kk) = (1.96*handles.surrogate.dyncond.std(kk))+handles.surrogate.dyncond.mean(kk);
    end
    %% CALCULATE DYNAMIC CONNECTIVITY (CORRELATION PER TIME WINDOW)
    %Whole session
    wininx = 1:winslide:length(handles.surrogate.data(:,1,1))-(winlength);
    for ww = 1:length(wininx)
        handles.corr.dynamicwhole(:,:,ww)= corr(handles.EDA(wininx(ww):wininx(ww)+winlength,:));
    end
    
    for kk = 1:length(handles.conditions),
        wininx = [handles.conditions(kk).DYNinx{:}];
        for ww = 1:length(wininx);
            set(handles.status_edit,'String',['Calculating dynamic correlations: ' num2str(ww) ' of ' num2str(length(wininx)) ' for condition: ' num2str(kk)]); pause(0.001);
            if wininx(ww)+winlength<length(handles.surrogate.data(:,1,pp))
            handles.corr.dynamic(kk).R(:,:,ww) = corr(handles.EDA(wininx(ww):wininx(ww)+winlength,:));
            end
        end
    end
    
    %Within/without group
    for kk = 1:length(handles.conditions),
        wininx = [handles.conditions(kk).DYNinx{:}];
        for ii = 1:length(handles.data)
            inx = find([handles.data.group]==handles.data(ii).group); %find ingroup members
            iinx = setdiff(1:length(handles.data),inx); %remove ingroup members to make indices of outgroup members
            inx(inx==ii)=[]; %remove current subject (autocorrelate)
            for ww = 1:length(handles.corr.dynamic(kk).R(1,1,:))
               handles.corr.dynamic(kk).Rwithin(ww,ii)= nanmean(handles.corr.dynamic(kk).R(inx,ii,ww));
               handles.corr.dynamic(kk).Rbetween(ww,ii)= nanmean(handles.corr.dynamic(kk).R(iinx,ii,ww)); 
            end
        end
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE GRAPH METRICS FOR WHOLE SESSION VISUALISATION (GM OVER TIME)
set(handles.status_edit,'String','Calculating whole session graph metrics for visualisation'); pause(0.001);
X = handles.corr.dynamicwhole; %Raw dynamic correlation matrices
X(isnan(X))=0; %Reduce NaN to zero
Xt = X;
Xt(Xt<handles.surrogate.dynwhole.thresh)=NaN; %threshold by whole session threshold
% Node degree over time
for cc = 1:length(Xt(1,1,:)), %for each time window (cc)
    handles.corr.dynwhole.degree(cc,:)=nansum(~isnan(Xt(:,:,cc)))-1; %dynamicND(subject,timepoint,condition)
end
handles.corr.dynwhole.degreep= handles.corr.dynwhole.degree/(length(handles.data)-1); %As percentage of possible connections
handles.corr.dynwhole.degree(handles.corr.dynwhole.degree<0)=0;
handles.corr.dynwhole.degreep(handles.corr.dynwhole.degreep<0)=0;
% CONNECTIVITY BASED MEASURES
%Matrix measures - connectivity (R), summed connectivity (S), time of connection (C)
X(X==1) = NaN; %make autocorrelates NaNs
Xt(Xt==1) = NaN; %make autocorrelates NaNs
for cc = 1:length(Xt(1,1,:)), %for each time window (cc)
    for ii = 1:length(handles.data),
        %Subject measures
        handles.corr.dynwhole.FC(cc,ii) = nanmean(X(ii,:,cc)); %mean connectivity over time window
        handles.corr.dynwhole.FCt(cc,ii) = nanmean(Xt(ii,:,cc)); %mean connectivity over time window
        handles.corr.dynwhole.S(cc,ii) = nansum(X(ii,:,cc)); %summed connectivity over time window
        handles.corr.dynwhole.St(cc,ii) = nansum(Xt(ii,:,cc)); %summed connectivity over time window
        %Group measures
        inx = find([handles.data.group]==handles.data(ii).group); %find ingroup members
        iinx = setdiff(1:length(handles.data),inx); %remove ingroup members to make indices of outgroup members
        inx(inx==ii)=[]; %remove current subject (autocorrelate)
        handles.corr.dynwhole.FC_within(cc,ii) = mean(X(ii,inx,cc)); %within group connectivity
        handles.corr.dynwhole.FCt_within(cc,ii) = mean(Xt(inx,ii,cc)); %within group connectivity
        handles.corr.dynwhole.NDt_within(cc,ii) = sum(ceil(Xt(inx,ii,cc))); %within group degree
        handles.corr.dynwhole.NPt_within(cc,ii) = sum(ceil(Xt(inx,ii,cc)))./length(inx); %within group proportion
        handles.corr.dynwhole.FC_between(cc,ii) = mean(X(iinx,ii,cc)); %between group connectivity
        handles.corr.dynwhole.FCt_between(cc,ii) = mean(Xt(iinx,ii,cc)); %within group connectivity
        handles.corr.dynwhole.NDt_between(cc,ii) = sum(ceil(Xt(iinx,ii,cc))); %within group degree
        handles.corr.dynwhole.NPt_between(cc,ii) = sum(ceil(Xt(iinx,ii,cc)))./length(iinx); %within group proportion
    end
end
%Graph measures
handles.corr.dynwhole.GC = nanmean(handles.corr.dynwhole.FC,2);
handles.corr.dynwhole.GCt = nanmean(handles.corr.dynwhole.FCt,2);
handles.corr.dynwhole.GS = nanmean(handles.corr.dynwhole.S,2);
handles.corr.dynwhole.GSt = nanmean(handles.corr.dynwhole.St,2);
handles.corr.dynwhole.ND = nanmean(handles.corr.dynwhole.degree,2); %mean degree
handles.corr.dynwhole.NDs = nansum(handles.corr.dynwhole.degree,2); %summed degree
handles.corr.dynwhole.NDp = nanmean(handles.corr.dynwhole.degreep,2); %mean density
handles.corr.dynwhole.GC_within = nanmean(handles.corr.dynwhole.FC_within,2);
handles.corr.dynwhole.GCt_within = nanmean(handles.corr.dynwhole.FCt_within,2);
handles.corr.dynwhole.GC_between = nanmean(handles.corr.dynwhole.FC_between,2);
handles.corr.dynwhole.GCt_between = nanmean(handles.corr.dynwhole.FCt_between,2);
handles.corr.dynwhole.NDt_within = nanmean(handles.corr.dynwhole.NDt_within,2); %mean degree
handles.corr.dynwhole.NDt_between = nansum(handles.corr.dynwhole.NDt_between,2); %summed degree
handles.corr.dynwhole.NPt_within = nanmean(handles.corr.dynwhole.NPt_within,2); %mean density
handles.corr.dynwhole.NPt_between = nanmean(handles.corr.dynwhole.NPt_between,2); %mean density

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE GRAPH METRICS PER CONDITION
%Get condition data
for kk = 1:length(handles.conditions),
    set(handles.status_edit,'String',['Calculating graph metrics for condition: ' handles.conditions(kk).name]); pause(0.001);
%     inx = [];
%     for ii = 1:length(handles.conditions(kk).DYNinx), %Get indices for data for condition kk
%         inx = [inx handles.conditions(kk).DYNinx{ii}{1}(1):handles.conditions(kk).DYNinx{ii}{1}(2)];
%     end
%     X = handles.corr.dynamic(:,:,inx); %Get data for condition kk in X
    X = handles.corr.dynamic(kk).R; %get conditions timewindow correlations
    X(isnan(X))=0; %remove NaNs
    Xt = X;
    Xz = X;
    Xt(Xt<handles.surrogate.dyncond.thresh(kk))=NaN; %Remove subthreshold values
    Xz(Xz<handles.surrogate.dynwhole.thresh)=NaN; %Remove subthreshold values
    % Node degree per time sample
    for cc = 1:length(Xt(1,1,:)),
        handles.corr.dyn(kk).degree(cc,:)=nansum(Xt(:,:,cc)~=0)-1; %dynamicND(subject,timepoint,condition)
        handles.corr.dyn(kk).degreez(cc,:)=nansum(Xz(:,:,cc)~=0)-1; %dynamicND(subject,timepoint,condition)
    end
    handles.corr.dyn(kk).degreep(:,:)= handles.corr.dyn(kk).degree./(length(handles.data)-1); %As percentage of possible connections
    handles.corr.dyn(kk).degreepz(:,:)= handles.corr.dyn(kk).degreez./(length(handles.data)-1);
    % CONNECTIVITY BASED MEASURES
    %Covariance measures
    for ii = 1:length(handles.data),
        for jj = 1:length(handles.data),
            if ii==jj,
                handles.corr.dyn(kk).R(ii,jj) = NaN;
                handles.corr.dyn(kk).Rt(ii,jj) = NaN;
                 handles.corr.dyn(kk).Rz(ii,jj) = NaN;
                handles.corr.dyn(kk).S(ii,jj) = NaN;
                handles.corr.dyn(kk).St(ii,jj) = NaN;
                handles.corr.dyn(kk).Tt(ii,jj) = NaN;
                handles.corr.dyn(kk).Tp(ii,jj) = NaN;
                handles.corr.dyn(kk).Sz(ii,jj) = NaN;
                handles.corr.dyn(kk).Tz(ii,jj) = NaN;
                handles.corr.dyn(kk).Tpz(ii,jj) = NaN;
            else
                handles.corr.dyn(kk).R(ii,jj) = nanmean(X(ii,jj,:)); %mean connectivity over time window
                handles.corr.dyn(kk).Rt(ii,jj) = nanmean(Xt(ii,jj,:)); %mean connectivity over time window
                handles.corr.dyn(kk).Rz(ii,jj) = nanmean(Xz(ii,jj,:)); %mean connectivity over time window
                handles.corr.dyn(kk).S(ii,jj) = nansum(X(ii,jj,:)); %summed connectivity over time window
                handles.corr.dyn(kk).St(ii,jj) = nansum(Xt(ii,jj,:)); %summed connectivity over time window
                handles.corr.dyn(kk).Sz(ii,jj) = nansum(Xz(ii,jj,:)); %summed connectivity over time window
                Y = squeeze(Xt(ii,jj,:));
                handles.corr.dyn(kk).Tt(ii,jj) = length(Y(Y>0)); %Time of connectivity (temporal degree)
                handles.corr.dyn(kk).Tp(ii,jj) = length(Y(Y>0))/length(Y); %Time of connectivity (% of activity) temporal density
                Y = squeeze(Xz(ii,jj,:));
                handles.corr.dyn(kk).Tz(ii,jj) = length(Y(Y>0)); %Time of connectivity (temporal degree)
                handles.corr.dyn(kk).Tpz(ii,jj) = length(Y(Y>0))/length(Y); %Time of connectivity (% of activity) temporal density
            end
        end
    end
    %Subject measures
    handles.corr.dyn(kk).FC = nanmean(handles.corr.dyn(kk).R); %mean connectivity
    handles.corr.dyn(kk).FCt = nanmean(handles.corr.dyn(kk).Rt);
    handles.corr.dyn(kk).FCz = nanmean(handles.corr.dyn(kk).Rz);
    handles.corr.dyn(kk).FS = nanmean(handles.corr.dyn(kk).S); %mean connectivity strength
    handles.corr.dyn(kk).FSt = nanmean(handles.corr.dyn(kk).St);
    handles.corr.dyn(kk).FSz = nanmean(handles.corr.dyn(kk).Sz);
    handles.corr.dyn(kk).FTt = nanmean(handles.corr.dyn(kk).Tt); %mean temporal degree
    handles.corr.dyn(kk).FTp = nanmean(handles.corr.dyn(kk).Tp); %mean temporal density
    handles.corr.dyn(kk).ND = nanmean(handles.corr.dyn(kk).degree); %mean degree per participant
    handles.corr.dyn(kk).NDs = nansum(handles.corr.dyn(kk).degree); %summed degree per participant
    handles.corr.dyn(kk).NDp = nanmean(handles.corr.dyn(kk).degreep); %mean density per participant
    handles.corr.dyn(kk).FTz = nanmean(handles.corr.dyn(kk).Tz); %mean temporal degree
    handles.corr.dyn(kk).FTpz = nanmean(handles.corr.dyn(kk).Tpz); %mean temporal density
    handles.corr.dyn(kk).NDz = nanmean(handles.corr.dyn(kk).degreez); %mean degree per participant
    handles.corr.dyn(kk).NDsz = nansum(handles.corr.dyn(kk).degreez); %summed degree per participant
    handles.corr.dyn(kk).NDpz = nanmean(handles.corr.dyn(kk).degreepz); %mean density per participant
    %Group measures
    for ii = 1:length(handles.data), %for each participant
        inx = find([handles.data.group]==handles.data(ii).group); %find ingroup members
        iinx = setdiff(1:length(handles.data),inx); %remove ingroup members to make indices of outgroup members
        inx(inx==ii)=[]; %remove current subject (autocorrelate)
        handles.corr.dyn(kk).FC_within(ii) = mean(handles.corr.dyn(kk).R(inx,ii)); %within group connectivity
        handles.corr.dyn(kk).FC_between(ii) = mean(handles.corr.dyn(kk).R(iinx,ii)); %between group connectivity
        handles.corr.dyn(kk).FCt_within(ii) = mean(handles.corr.dyn(kk).Rt(inx,ii)); %within group connectivity
        handles.corr.dyn(kk).FCt_between(ii) = mean(handles.corr.dyn(kk).Rt(iinx,ii)); %between group connectivity
        handles.corr.dyn(kk).FS_within(ii) = mean(handles.corr.dyn(kk).S(inx,ii)); %mean within group connectivity strength
        handles.corr.dyn(kk).FS_between(ii) = mean(handles.corr.dyn(kk).S(iinx,ii)); %mean between group connectivity strength
        handles.corr.dyn(kk).FSt_within(ii) = mean(handles.corr.dyn(kk).St(inx,ii)); %mean within group connectivity strength
        handles.corr.dyn(kk).FSt_between(ii) = mean(handles.corr.dyn(kk).St(iinx,ii)); %mean between group connectivity strength
        handles.corr.dyn(kk).FTt_within(ii) = mean(handles.corr.dyn(kk).Tt(inx,ii)); %mean within group temporal degree
        handles.corr.dyn(kk).FTt_between(ii) = mean(handles.corr.dyn(kk).Tt(iinx,ii)); %mean between group temporal degree
        handles.corr.dyn(kk).FTp_within(ii) = mean(handles.corr.dyn(kk).Tp(inx,ii)); %mean within group temporal density
        handles.corr.dyn(kk).FTp_between(ii) = mean(handles.corr.dyn(kk).Tp(iinx,ii)); %mean between group temporal density
        handles.corr.dyn(kk).ND_within(ii) = mean(sum(ceil(Xt(inx,ii,:)))); %within group degree
        handles.corr.dyn(kk).ND_between(ii) = mean(sum(ceil(Xt(iinx,ii,:)))); %between group degree
        handles.corr.dyn(kk).NDs_within(ii) = nansum(sum(ceil(Xt(inx,ii,:)))); %within group summed degree
        handles.corr.dyn(kk).NDs_between(ii) = nansum(sum(ceil(Xt(iinx,ii,:)))); %between group summed degree
        handles.corr.dyn(kk).NDp_within(ii) = mean(sum(ceil(Xt(inx,ii,:)))./length(inx)); %within group proportion
        handles.corr.dyn(kk).NDp_between(ii) = mean(sum(ceil(Xt(iinx,ii,:)))./length(iinx)); %within group proportion
    end
    
    
    %Graph measures
    handles.corr.dyn(kk).GC = nanmean(handles.corr.dyn(kk).FC);
    handles.corr.dyn(kk).GCt = nanmean(handles.corr.dyn(kk).FCt);
    handles.corr.dyn(kk).GS = nanmean(handles.corr.dyn(kk).FS);
    handles.corr.dyn(kk).GSt = nanmean(handles.corr.dyn(kk).FSt);
    handles.corr.dyn(kk).GTt = nanmean(handles.corr.dyn(kk).FTt);
    handles.corr.dyn(kk).GTp = nanmean(handles.corr.dyn(kk).FTp);
    handles.corr.dyn(kk).GD = nanmean(handles.corr.dyn(kk).ND);
    handles.corr.dyn(kk).GDs = nanmean(handles.corr.dyn(kk).NDs);
    handles.corr.dyn(kk).GDp = nanmean(handles.corr.dyn(kk).NDp);
    handles.corr.dyn(kk).GC_within = nanmean(handles.corr.dyn(kk).FC_within);
    handles.corr.dyn(kk).GC_between = nanmean(handles.corr.dyn(kk).FC_between);
    handles.corr.dyn(kk).GCt_within = nanmean(handles.corr.dyn(kk).FCt_within);
    handles.corr.dyn(kk).GCt_between = nanmean(handles.corr.dyn(kk).FCt_between);
    handles.corr.dyn(kk).GS_within = nanmean(handles.corr.dyn(kk).FS_within);
    handles.corr.dyn(kk).GS_between = nanmean(handles.corr.dyn(kk).FS_between);
    handles.corr.dyn(kk).GSt_within = nanmean(handles.corr.dyn(kk).FSt_within);
    handles.corr.dyn(kk).GSt_between = nanmean(handles.corr.dyn(kk).FSt_between);
    handles.corr.dyn(kk).GTt_within = nanmean(handles.corr.dyn(kk).FTt_within);
    handles.corr.dyn(kk).GTt_between = nanmean(handles.corr.dyn(kk).FTt_between);
    handles.corr.dyn(kk).GTp_within = nanmean(handles.corr.dyn(kk).FTp_within);
    handles.corr.dyn(kk).GTp_between = nanmean(handles.corr.dyn(kk).FTp_between);
    handles.corr.dyn(kk).GD_within = nanmean(handles.corr.dyn(kk).ND_within);
    handles.corr.dyn(kk).GD_between = nanmean(handles.corr.dyn(kk).ND_between);
    handles.corr.dyn(kk).GDs_within = nanmean(handles.corr.dyn(kk).NDs_within);
    handles.corr.dyn(kk).GDs_between = nanmean(handles.corr.dyn(kk).NDs_between);
    handles.corr.dyn(kk).GDp_within = nanmean(handles.corr.dyn(kk).NDp_within);
    handles.corr.dyn(kk).GDp_between = nanmean(handles.corr.dyn(kk).NDp_between);
end
%% Set up strings for connectivity menu
string(1).t= [];
if isfield(handles.corr,'static'),
    string(1).t='Static';
    if isfield(handles.corr,'dynamic'),
        string(end+1).t = 'Dynamic';
    end
else string(1).t='Dynamic';
end
set(handles.corconn_pop,'String',{'Connectivity: ' string.t});
handles.parameters.corr=1;
if isfield(handles,'edaaxes'),
    set(handles.edaaxes,'Visible','off');
    axes(handles.edaaxes); cla;
end
set(handles.status_edit,'String','Dynamic connectivity analysis complete'); pause(0.001);
set(handles.cordyn_text,'String','Complete'); %change text to complete
guidata(hObject,handles)

function handles = corsurrogate(handles)
data = handles.EDA;
data(isnan(data))=0;
[points,subjects] = size(data); %get size of data
if rem(points,2)==0
    points = points-1;
    data=data(1:points,:);
end
K = (points-1)/2;
M = 2:K+1;
N = K+2:points;
fft_data = fft(data); %perform FFT
handles.surrogate.data = zeros(points,subjects,1000);
for kk = 1:1000,
    set(handles.status_edit,'String',['Creating null distribution: Surrogate dataset ' num2str(kk) '/1000']);
    pause(0.001);
    rnddata = rand([K 1]);
    rndphase1 = repmat(exp (2*pi*1i*rnddata),1,subjects);
    rndphase2 = conj( flipud( rndphase1));
    fft_data_surr = fft_data;
    fft_data_surr(M,:) = fft_data(M,:).*rndphase1;
    fft_data_surr(N,:) = fft_data(N,:).*rndphase2;
    handles.surrogate.data(:,:,kk) = real(ifft(fft_data_surr));
end


function corlag_toggle_Callback(hObject, eventdata, handles)

function corlag_edit_Callback(hObject, eventdata, handles)

function corlag_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function corwin_edit_Callback(hObject, eventdata, handles)

function corwin_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function corslide_edit_Callback(hObject, eventdata, handles)

function corslide_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function croptime_toggle_Callback(hObject, eventdata, handles)
handles = preprocess(handles);
handles = update_EDAdisplay(handles);
guidata(hObject,handles)


% --- Executes on button press in clearwork_push.
function clearwork_push_Callback(hObject, eventdata, handles)
handles = clearvariables(handles);
guidata(hObject, handles);

function handles = clearvariables(handles)
%% CREATE EMPTY VARIABLES
handles.sdir=0;
handles.data=[];
handles.parameters=[];
handles.time=[];
handles.display=[];
handles.EDA=[];
handles.ACC=[];
handles.groups=[];
handles.group_subj = [];
handles.conditions = [];
handles.surrogate = [];
handles.corr = [];

set(handles.edaaxes,'Visible','off');
set(handles.GM1axes,'Visible','off');
set(handles.GM2axes,'Visible','off');
set(handles.vis_uipanel,'Visible','off');
set(handles.preprocess_uipanel,'Visible','off');
set(handles.align_uipanel,'Visible','off');
set(handles.group_uipanel,'Visible','off');
set(handles.condition_uipanel,'Visible','off');
set(handles.GM_uipanel,'Visible','off');
set(handles.corstat_text,'String','Ready...');
set(handles.cordyn_text,'String','Ready...');



function condload_push_Callback(hObject, eventdata, handles)
filename = 0;
[filename, pathname] = uigetfile({'*.xls'},'Select edit list file...');
if isequal(filename,0), return; end %if nothing selected, stop function
handles.conditions = [];
group_t = importdata([pathname filename]);
for kk = 2:length(group_t.textdata(:,1)), %skip header, for each instance
    on=group_t.textdata{kk,1};
    off = group_t.textdata{kk,2};
    t1= (str2double(on(1:2))*60)+(str2double(on(4:5)))+(str2double(on(7:8))/60)+(str2double(on(10:11))/600);
    t2= (str2double(off(1:2))*60)+(str2double(off(4:5)))+(str2double(off(7:8))/60)+(str2double(on(10:11))/600);
    %find index of closest matching time
    [c1, inx1] = min(abs(handles.time.EDAmins-t1));
    [c2, inx2] = min(abs(handles.time.EDAmins-t2));
    c1 = handles.time.EDAmins(inx1); %convert c to time
    c2 = handles.time.EDAmins(inx2); %convert c to time
    condname = group_t.textdata{kk,3};
    %Check if condition exists
    check=0;
    for ii = 1:length(handles.conditions),
        if strcmpi(condname,handles.conditions(ii).name)==1,
            handles.conditions(ii).windows{end+1}={[c1 c2]}; %assign times
            handles.conditions(ii).windowsinx{end+1}={[inx1 inx2]};
            handles.conditions(ii).string = [condname ': ' num2str(length(handles.conditions(ii).windows)) ' time windows'];
            check=1;
        end
    end
    if check==0,
        handles.conditions(end+1).name = condname;
        handles.conditions(end).windows{1}={[c1 c2]}; %assign times
        handles.conditions(end).windowsinx{1}={[inx1 inx2]};
        handles.conditions(end).string = [condname ': 1 time window'];
    end
end
set(handles.condstart_edit,'String',num2str(0));
set(handles.condend_edit,'String',num2str(0));
set(handles.cond_list,'String',{handles.conditions.string});
set(handles.condcreate_push,'String','Add to condition');
handles = update_EDAdisplay(handles);
guidata(hObject, handles);

function corleft_radio_Callback(hObject, eventdata, handles)

function corright_radio_Callback(hObject, eventdata, handles)

function corconn_pop_Callback(hObject, eventdata, handles)
inx = get(handles.corconn_pop,'Value');
string = get(handles.corconn_pop,'String');
if strcmp(string{inx},'Connectivity:'),
    set(handles.corlevel_pop,'Value',1);
    set(handles.corlevel_pop,'String','Level:');
    handles.parameters.connectivity = 0;
else
    switch string{inx},
        case 'Static', handles.parameters.connectivity = 1;
            set(handles.corlevel_pop,'String',{'Level:' 'Covariance' 'Subject' 'Group dynamics' 'Group comparison' 'Graph'});
        case 'Dynamic', handles.parameters.connectivity = 2;
            set(handles.corlevel_pop,'String',{'Level:' 'Covariance' 'Subject' 'Group dynamics' 'Group comparison' 'Graph' 'Inspect session'});
    end
    set(handles.corlevel_pop,'Value',1);
    handles.parameters.level=0;
    set(handles.cormetric_pop,'Value',1);
    set(handles.cormetric_pop,'String','Metric:');
    handles.parameters.metric=0;
    set(handles.corcond_pop,'Value',1);
    set(handles.corcond_pop,'String','Condition:');
    handles.parameters.condition=0;
    set(handles.corthresh_pop,'Value',1);
    set(handles.corthresh_pop,'String','Edge threshold:');
    handles.parameters.condition=0;
end
guidata(hObject,handles)

function corconn_pop_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function corlevel_pop_Callback(hObject, eventdata, handles)
inx = get(handles.corlevel_pop,'Value');
string = get(handles.corlevel_pop,'String');
if strcmp(string{inx},'Level:'),
    set(handles.cormetric_pop,'Value',1);
    set(handles.cormetric_pop,'String','Metric:');
    handles.parameters.level = 0;
elseif handles.parameters.connectivity==1, %Static connectivity
    switch string{inx},
        case 'Covariance', handles.parameters.level = 1;
            set(handles.cormetric_pop,'Value',1);
            set(handles.cormetric_pop,'String',{'Metric:' 'Correlation'});
        case 'Subject', handles.parameters.level = 2;
            set(handles.cormetric_pop,'Value',1);
            set(handles.cormetric_pop,'String',{'Metric:' 'Mean connectivity' ...
                'Connectivity strength' 'Degree' 'Density'});
        case 'Group dynamics', handles.parameters.level = 3;
            set(handles.cormetric_pop,'Value',1);
            set(handles.cormetric_pop,'String',{'Metric:' 'Within group connectivity' ...
                'Between group connectivity' 'Within group strength' 'Between group strength' ...
                'Within group degree' 'Between group degree' 'Within group density' ...
                'Between group density'});
        case 'Group comparison', handles.parameters.level = 4;
            set(handles.cormetric_pop,'Value',1);
            set(handles.cormetric_pop,'String',{'Metric:' 'Within group connectivity' ...
                'Between group connectivity' 'Within group strength' 'Between group strength' ...
                'Within group degree' 'Between group degree' 'Within group density' ...
                'Between group density'});
        case 'Graph', handles.parameters.level = 5;
            set(handles.cormetric_pop,'Value',1);
            set(handles.cormetric_pop,'String',{'Metric:' 'Mean connectivity' ...
                'Connectivity strength' 'Degree' 'Density' 'Within group connectivity' ...
                'Between group connectivity' 'Within group strength' 'Between group strength' ...
                'Within group degree' 'Between group degree' 'Within group density' ...
                'Between group density'});
    end
elseif handles.parameters.connectivity==2, %Dynamic connectivity
    switch string{inx},
        case 'Covariance', handles.parameters.level = 1;
            set(handles.cormetric_pop,'String',{'Metric:' 'Mean connectivity' ...
                'Mean connectivity strength' 'Temporal degree' 'Temporal density'});
        case 'Subject', handles.parameters.level = 2;
            set(handles.cormetric_pop,'String',{'Metric:' 'Mean connectivity' ...
                'Mean connectivity strength' 'Mean temporal degree' ...
                'Mean temporal density' 'Mean degree' 'Sum degree' ...
                'Mean density'});
        case 'Group dynamics', handles.parameters.level = 3;
            set(handles.cormetric_pop,'String',{'Metric:' 'Within group connectivity' ...
                'Between group connectivity' 'Within group strength' ...
                'Between group strength' 'Within group temporal degree' ...
                'Between group temporal degree ' 'Within group temporal density' ...
                'Between group temporal density' 'Within group mean degree' ...
                'Between group mean degree' 'Within group summed degree' ...
                'Within group density' 'Between group density'});
        case 'Group comparison', handles.parameters.level = 4;
            set(handles.cormetric_pop,'String',{'Metric:' 'Within group connectivity' ...
                'Between group connectivity' 'Within group degree' 'Between group degree' ...
                'Within group density' 'Between group density'});
        case 'Graph', handles.parameters.level = 5;
            set(handles.cormetric_pop,'String',{'Metric:' 'Mean connectivity' ...
                'Mean connectivity strength' 'Mean temporal degree' ...
                'Mean temporal density' 'Mean degree' 'Sum degree' ...
                'Mean density' 'Within group connectivity' ...
                'Between group connectivity' 'Within group strength' ...
                'Between group strength' 'Within group temporal degree' ...
                'Between group temporal degree ' 'Within group temporal density' ...
                'Between group temporal density' 'Within group mean degree' ...
                'Between group mean degree' 'Within group summed degree' ...
                'Within group density' 'Between group density'});
        case 'Inspect session', handles.parameters.level = 6;
            set(handles.cormetric_pop,'String',{'Metrics:' 'Connectivity' ...
                'Connectivity strength' 'Degree' 'Density' ...
                'Within group connectivity' 'Between group connectivity' ...
                'Within group degree' 'Between group degree' 'Within group density'...
                'Between group density' });
    end
end
set(handles.cormetric_pop,'Value',1);
handles.parameters.metric=0;
set(handles.corcond_pop,'Value',1);
set(handles.corcond_pop,'String','Condition:');
handles.parameters.condition=0;
set(handles.corthresh_pop,'Value',1);
set(handles.corthresh_pop,'String','Edge threshold:');
handles.parameters.thresh=0;
guidata(hObject,handles)

function corlevel_pop_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function cormetric_pop_Callback(hObject, eventdata, handles)
inx = get(handles.cormetric_pop,'Value');
string = get(handles.cormetric_pop,'String');
if strcmp(string{inx},'Metric:'),
    set(handles.corcond_pop,'Value',1);
    set(handles.corcond_pop,'String','Condition:');
    handles.parameters.metric = 0;
elseif handles.parameters.connectivity==1, %Static metrics
    switch string{inx},
        case 'Correlation', handles.parameters.metric=1; handles.parameters.thresht = 3;
        case 'Mean connectivity', handles.parameters.metric=2; handles.parameters.thresht = 3;
        case 'Connectivity strength', handles.parameters.metric=3; handles.parameters.thresht = 3;
        case 'Degree', handles.parameters.metric=4; handles.parameters.thresht = 2;
        case 'Density', handles.parameters.metric=5; handles.parameters.thresht = 2;
        case 'Within group connectivity', handles.parameters.metric=6; handles.parameters.thresht = 3;
        case 'Between group connectivity', handles.parameters.metric=7; handles.parameters.thresht = 3;
        case 'Within group strength', handles.parameters.metric=8; handles.parameters.thresht = 3;
        case 'Between group strength', handles.parameters.metric=9; handles.parameters.thresht = 3;
        case 'Within group degree', handles.parameters.metric=10; handles.parameters.thresht = 2;
        case 'Between group degree', handles.parameters.metric=11; handles.parameters.thresht = 2;
        case 'Within group density', handles.parameters.metric=12; handles.parameters.thresht = 2;
        case 'Between group density', handles.parameters.metric=13; handles.parameters.thresht = 2;
    end
    if handles.parameters.level==1,
        set(handles.corcond_pop,'String',{'Conditions: ' handles.conditions.name}); %set up string for condition dropdown
    else set(handles.corcond_pop,'String',{'Conditions: ' handles.conditions.name 'Compare'});
    end
elseif handles.parameters.connectivity==2, %Dynamic metrics
    if handles.parameters.level == 6, %if whole session inspection
        switch string{inx},
            case 'Connectivity', handles.parameters.metric=1; handles.parameters.thresht = 6;
            case 'Connectivity strength', handles.parameters.metric=2; handles.parameters.thresht = 6;
            case 'Degree', handles.parameters.metric=3; handles.parameters.thresht = 7;
            case 'Density', handles.parameters.metric=4; handles.parameters.thresht = 7;
            case 'Within group connectivity', handles.parameters.metric=5; handles.parameters.thresht = 6;
            case 'Between group connectivity', handles.parameters.metric=6; handles.parameters.thresht = 6;
            case 'Within group degree', handles.parameters.metric=7; handles.parameters.thresht = 7;
            case 'Between group degree', handles.parameters.metric=8; handles.parameters.thresht = 7;
            case 'Within group density', handles.parameters.metric=9; handles.parameters.thresht = 7;
            case 'Between group density', handles.parameters.metric=10; handles.parameters.thresht = 7;
        end
    else
        switch string{inx},
            case 'Mean connectivity', handles.parameters.metric=1; handles.parameters.thresht = 4;
            case 'Mean connectivity strength', handles.parameters.metric=2; handles.parameters.thresht = 4;
            case 'Mean temporal degree', handles.parameters.metric=3; handles.parameters.thresht = 5;
            case 'Mean temporal density', handles.parameters.metric=4; handles.parameters.thresht = 5;
            case 'Mean degree', handles.parameters.metric=5; handles.parameters.thresht = 5;
            case 'Sum degree', handles.parameters.metric=6; handles.parameters.thresht = 5;
            case 'Mean density', handles.parameters.metric=7; handles.parameters.thresht = 5;
            case 'Within group connectivity', handles.parameters.metric=8; handles.parameters.thresht = 4;
            case 'Between group connectivity', handles.parameters.metric=9; handles.parameters.thresht = 4;
            case 'Within group strength', handles.parameters.metric=10; handles.parameters.thresht = 4;
            case 'Between group strength', handles.parameters.metric=11; handles.parameters.thresht = 4;
            case 'Within group degree', handles.parameters.metric=12; handles.parameters.thresht = 5;
            case 'Between group degree', handles.parameters.metric=13; handles.parameters.thresht = 5;
            case 'Within group density', handles.parameters.metric=14; handles.parameters.thresht = 5;
            case 'Between group density', handles.parameters.metric=15; handles.parameters.thresht = 5;
            case 'Within group temporal density', handles.parameters.metric=16; handles.parameters.thresht = 5;
            case 'Between group temporal density', handles.parameters.metric=17; handles.parameters.thresht = 5;
        end
    end
    if handles.parameters.level==1,
        set(handles.corcond_pop,'String',{'Conditions: ' handles.conditions.name}); %set up string for condition dropdown
    elseif handles.parameters.level==6,
        set(handles.corcond_pop,'String',{'Conditions: ' 'Whole session'});
    else set(handles.corcond_pop,'String',{'Conditions: ' handles.conditions.name 'Compare'});
    end
end
set(handles.corcond_pop,'Value',1);
handles.parameters.condition=0;
set(handles.corthresh_pop,'Value',1);
set(handles.corthresh_pop,'String','Edge threshold:');
handles.parameters.thresh=0;
guidata(hObject,handles)

function cormetric_pop_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function corcond_pop_Callback(hObject, eventdata, handles)
inx = get(handles.corcond_pop,'Value');
string = get(handles.corcond_pop,'String');
if strcmp(string{inx},'Condition:'),
    set(handles.corthresh_pop,'Value',1);
    set(handles.corthesh_pop,'String','Significance threshold:');
    handles.parameters.condition = 0;
elseif strcmp(string{inx},'Whole session'),
    handles.parameters.condition = length(handles.conditions)+2;
else
    handles.parameters.condition = get(handles.corcond_pop,'Value');
end
switch handles.parameters.thresht,
    case 2, set(handles.corthresh_pop,'String',{'Edge threshold:' ...
            'Whole session null distribution' 'Condition specific null distribution'});
    case 3, set(handles.corthresh_pop,'String',{'Edge threshold:' ...
            'No threshold - correlation strength' 'Whole session null distribution' ...
            'Condition specific null distribution'});
    case 4, set(handles.corthresh_pop,'String',{'Edge threshold:' ...
            'No threshold - correlation strength' 'Condition specific null distribution'});
    case 5, set(handles.corthresh_pop,'String',{'Edge threshold:' ...
            'Condition specific null distribution'});
    case 6, set(handles.corthresh_pop,'String',{'Edge threshold:' ...
            'No threshold - correlation strength' 'Whole session null distribution'});
    case 7, set(handles.corthresh_pop,'String',{'Edge threshold:' ...
            'Whole session null distribution'});
end
set(handles.corthresh_pop,'Value',1);
handles.parameters.thresh=0;
guidata(hObject,handles)

function corcond_pop_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function corthresh_pop_Callback(hObject, eventdata, handles)
inx = get(handles.corthresh_pop,'Value');
string = get(handles.corthresh_pop,'String');
if strcmp(string{inx},'Edge threshold:'),
    handles.parameters.thresh = 0;
else
    switch string{inx},
        case 'Whole session null distribution', handles.parameters.thresh = 2;
        case 'Condition specific null distribution', handles.parameters.thresh = 3;
        case 'No threshold - correlation strength', handles.parameters.thresh = 1;
    end
end
handles.parameters.corr=1;
handles = update_EDAdisplay(handles);
guidata(hObject,handles)

function corthresh_pop_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function predone_button_Callback(hObject, eventdata, handles)
set(handles.preprocess_uipanel,'Visible','off');
set(handles.align_uipanel,'Visible','on');
set(handles.group_uipanel,'Visible','on');
set(handles.condition_uipanel,'Visible','on');
set(handles.GM_uipanel,'Visible','on');

function crop1_edit_Callback(hObject, eventdata, handles)
    handles = preprocess(handles);
    handles = update_EDAdisplay(handles);
    guidata(hObject, handles);

function crop1_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function crop2_edit_Callback(hObject, eventdata, handles)
    handles = preprocess(handles);
    handles = update_EDAdisplay(handles);
    guidata(hObject, handles);

function crop2_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function corstatus_list_Callback(hObject, eventdata, handles)

function corstatus_list_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in static_checkbox.
function static_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to static_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of static_checkbox


% --- Executes on button press in dynamic_checkbox.
function dynamic_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to dynamic_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dynamic_checkbox


% --- Executes on button press in mtd_checkbox.
function mtd_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to mtd_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of mtd_checkbox


% --- Executes on button press in surrogate_checkbox.
function surrogate_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to surrogate_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of surrogate_checkbox


% --- Executes on button press in croptime_checkbox.
function croptime_checkbox_Callback(hObject, eventdata, handles)
    handles = preprocess(handles);
    handles = update_EDAdisplay(handles);
    guidata(hObject, handles);
    
function filter_checkbox_Callback(hObject, eventdata, handles)
    handles = preprocess(handles);
    handles = update_EDAdisplay(handles);
    guidata(hObject, handles);
    
function tempsmooth_edit_Callback(hObject, eventdata, handles)
    handles = preprocess(handles);
    handles = update_EDAdisplay(handles);
    guidata(hObject, handles);

function tempsmooth_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function tempsmooth_checkbox_Callback(hObject, eventdata, handles)
    handles = preprocess(handles);
    handles = update_EDAdisplay(handles);
    guidata(hObject, handles);

function detrend_checkbox_Callback(hObject, eventdata, handles)
    handles = preprocess(handles);
    handles = update_EDAdisplay(handles);
    guidata(hObject, handles);

function z_checkbox_Callback(hObject, eventdata, handles)
    handles = preprocess(handles);
    handles = update_EDAdisplay(handles);
    guidata(hObject, handles);

function zero_checkbox_Callback(hObject, eventdata, handles)
    handles = preprocess(handles);
    handles = update_EDAdisplay(handles);
    guidata(hObject, handles);


% --- Executes on button press in batch_push.
function batch_push_Callback(hObject, eventdata, handles)
% hObject    handle to batch_push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
