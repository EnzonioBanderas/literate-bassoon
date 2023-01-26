function func_plotPhysioByTrial(filePathPhysio,filePathEvents,filePathEyetrack,varargin)
%% 
%% Load physio, eyetracking and event data for one participant

fn_events = convertStringsToChars(filePathEvents);
fn_eyetrack = convertStringsToChars(filePathEyetrack);
fn_physio = convertStringsToChars(filePathPhysio);
% fn_diamond = '/media/diskEvaluation/Evaluation/sfb1280a05study7/misc/diamond.png';
% fn_diamond = fullfile('sfb1280a05study7', 'misc', 'diamond.png');
% fn_square  = '/media/diskEvaluation/Evaluation/sfb1280a05study7/misc/square.png';
% img_d = imread(fn_diamond);
% img_s = imread(fn_square);

bf = essbids_parseLabel(fn_physio);
figTitle = sprintf('sub-%s  ses-%s  run-%s  trial-',bf.sub,bf.ses,bf.run);

if isempty(varargin)
    figure_setting = 1;
elseif length(varargin) == 1
    figure_setting = varargin{1};
end

curr_idx = 1;
%% Loading in Physio data
t_physio = essbids_readTsv(fn_physio);
t_physio.time = t_physio.Properties.CustomProperties.Time;

%% Loading in Eyetrack data
[t_eyetrack, time_eyetrack] = essbids_readTsv(fn_eyetrack);
t_eyetrack.time = t_eyetrack.Properties.CustomProperties.Time;
if figure_setting == 2
    t_eyetrack_pro = cleanEyeTrackData(t_eyetrack,time_eyetrack,0,.5);
    t_eyetrack_pro7 = cleanEyeTrackData(t_eyetrack,time_eyetrack,0.6,.5);
    t_eyetrack_pro8 = cleanEyeTrackData(t_eyetrack,time_eyetrack,0.8,.5);
end

%% Loading in Events data & process it to be plotted
t_events = essbids_readTsv(fn_events);
t_eventmarkers = get_Eventmarkers(t_events);
log_type = ismember(t_events.trial_type,{'US'});
t_events(log_type,:) = [];


hFig = figure('Units','normalized','OuterPosition',[0.1 0 1 1],'Name',bf.sub);
%% Plotting Screen
% if figure_setting == 3
%     t_eyetrack.EyeA_X_Gaze = t_eyetrack.EyeA_X_Gaze .* 1980;
%     t_eyetrack.EyeB_X_Gaze = t_eyetrack.EyeB_X_Gaze .* 1980;
%     
%     t_eyetrack.EyeA_Y_Gaze = t_eyetrack.EyeA_Y_Gaze .* 1080;
%     t_eyetrack.EyeB_Y_Gaze = t_eyetrack.EyeB_Y_Gaze .* 1080;
%     
%     ax_screen = subplot(4,2,1:4);
%     temp1 = t_events.stim_file{1};
%     log_trial = t_eyetrack.time > t_events.onset(curr_idx) - 4 & t_eyetrack.time < t_events.onset(curr_idx) + 6;
%     if contains(temp1,'diamond')
%         imshow(img_d);
%         axis on;
%     elseif contains(temp1,'square')
%         imshow(img_s);
%         axis on;
%     end
%     hold on
%     scatter(t_eyetrack.EyeA_X_Gaze(log_trial),t_eyetrack.EyeA_Y_Gaze(log_trial),'green','filled','o','MarkerFaceAlpha',0.1);
%     scatter(t_eyetrack.EyeB_X_Gaze(log_trial),t_eyetrack.EyeB_Y_Gaze(log_trial),'blue','filled','o','MarkerFaceAlpha',0.1);
%     legend('Eye A','Eye B')
%     title(string(figTitle) + string(t_events.trial_index(1)));
%     ax_screen.XLim = [-200 2180];
%     ax_screen.YLim = [-200 1280];
% end    
%% Plotting the Eventmarkers
if figure_setting < 3
    ax_events = subplot(5,2,9:10);
elseif figure_setting == 3
    ax_events = subplot(4,2,7:8);
end
plot(t_eventmarkers.timestamps(t_eventmarkers.marker_value == 1 | t_eventmarkers.marker_value == 0),t_eventmarkers.marker_value(t_eventmarkers.marker_value == 1 | t_eventmarkers.marker_value == 0));
hold on
plot(t_eventmarkers.timestamps(t_eventmarkers.marker_value == 1.1 | t_eventmarkers.marker_value == 0),t_eventmarkers.marker_value(t_eventmarkers.marker_value == 1.1 | t_eventmarkers.marker_value == 0));

if nnz(t_eventmarkers.marker_value == 1.3) > 1
    plot(t_eventmarkers.timestamps(t_eventmarkers.marker_value == 1.3 | t_eventmarkers.marker_value == 0), t_eventmarkers.marker_value(t_eventmarkers.marker_value == 1.3 | t_eventmarkers.marker_value == 0));
    legend('CSminus','CSplus','US');
else
    legend('CSminus','CSplus');
end
hold off
%Changing lims to represent the first trial.
ax_events.XLim = [t_events.onset(curr_idx)-4,t_events.onset(curr_idx)+t_events.duration(curr_idx)+6];
ax_events.YLim = [0,1.5];
ax_events.YTick = [];
xlabel("Time [s]");
%% Plot Eyetrack Area
if figure_setting == 1
    ax_eyetrack = subplot(5,2,7:8);
elseif figure_setting == 3
    ax_eyetrack = subplot(4,2,5:6);
end
plot(t_eyetrack.time,t_eyetrack.EyeA_Area,'Color','green');
hold on 
plot(t_eyetrack.time,t_eyetrack.EyeB_Area,'Color','blue');
legend('EyeA Area','EyeB Area');
hold off
%Changing lims to represent the first trial
ax_eyetrack.XLim = [t_events.onset(curr_idx)-4,t_events.onset(curr_idx)+t_events.duration(curr_idx)+6];
ax_eyetrack.XTick = 0;
ylabel("[a.u]");
if figure_setting == 1
    %% Plot EDA
    ax_eda = subplot(5,2,5:6);
    plot(t_physio.time,t_physio.skinconductance,'Color','red');
    legend('EDA');
    %Changing lims to represent the first trial
    ax_eda.XLim = [t_events.onset(curr_idx)-4,t_events.onset(curr_idx)+t_events.duration(curr_idx)+6];
    ax_eda.XTick = 0;
    ylabel('[µS]');
    
    %% Plot Pulseoxy
    ax_pulse = subplot(5,2,3:4);
    plot(t_physio.time,t_physio.pulseoximeter,'Color','blue');
    legend('Pulse');
    %Changing lims to represent the first trial
    ax_pulse.XLim = [t_events.onset(curr_idx)-4,t_events.onset(curr_idx)+t_events.duration(curr_idx)+6];
    ax_pulse.XTick = 0;
    ylabel('[V]');
    
    %% Plot Respiratory
    ax_resp = subplot(5,2,1:2);
    plot(t_physio.time,t_physio.respiratory,'Color','green');
    legend('Respiratory');
    %Changing lims to represent the first trial
    ax_resp.XLim = [t_events.onset(1)-4,t_events.onset(1)+t_events.duration(1)+6];
    ax_resp.XTick = 0;
    ylabel("[cmH2O]");
    title(string(figTitle) + string(t_events.trial_index(1)));
elseif figure_setting == 2
    %% Plot Eyetrack Raw Area
    ax_eyetrack_raw = subplot(5,2,1:2);
    plot(t_eyetrack.time,t_eyetrack.EyeA_Area,'Color','green');
    hold on 
    plot(t_eyetrack.time,t_eyetrack.EyeB_Area,'Color','blue');
    legend('EyeA Area','EyeB Area');
    hold off
    %Changing lims to represent the first trial
    ax_eyetrack_raw.XLim = [t_events.onset(curr_idx)-4,t_events.onset(curr_idx)+t_events.duration(curr_idx)+6];
    ax_eyetrack_raw.XTick = 0;
    title(string(figTitle) + string(t_events.trial_index(1)));
    ylabel("[a.u]");
     %% Plot Eyetrack Processed Area (Calculated)
    ax_eyetrack_pro1 = subplot(5,2,3:4);
    plot(t_eyetrack_pro.time,t_eyetrack_pro.EyeA_Area,'Color','green');
    hold on 
    plot(t_eyetrack_pro.time,t_eyetrack_pro.EyeB_Area,'Color','blue');
    legend('EyeA Area','EyeB Area');
    hold off
    %Changing lims to represent the first trial
    ax_eyetrack_pro1.XLim = [t_events.onset(curr_idx)-4,t_events.onset(curr_idx)+t_events.duration(curr_idx)+6];
    ax_eyetrack_pro1.Title.String = sprintf("After Enzos Algorithm (AR=%s, Blink-time: %ds)",t_eyetrack_pro.Properties.CustomProperties.AR,t_eyetrack_pro.Properties.CustomProperties.MaxBlinkTime);
    ax_eyetrack_pro1.XTick = 0;
    ylabel("[a.u]");
     %% Plot Eyetrack Processed Area (AR: 0.7)
    ax_eyetrack_pro2 = subplot(5,2,5:6);
    plot(t_eyetrack_pro7.time,t_eyetrack_pro7.EyeA_Area,'Color','green');
    hold on 
    plot(t_eyetrack_pro7.time,t_eyetrack_pro7.EyeB_Area,'Color','blue');
    legend('EyeA Area','EyeB Area');
    hold off
    %Changing lims to represent the first trial
    ax_eyetrack_pro2.XLim = [t_events.onset(curr_idx)-4,t_events.onset(curr_idx)+t_events.duration(curr_idx)+6];
    ax_eyetrack_pro2.Title.String = sprintf("After Enzos Algorithm (AR=%s, Blink-time: %ds)",t_eyetrack_pro7.Properties.CustomProperties.AR,t_eyetrack_pro7.Properties.CustomProperties.MaxBlinkTime);
    ax_eyetrack_pro2.XTick = 0;
    ylabel("[a.u]");
     %% Plot Eyetrack Processed Area (AR: 0.8)
    ax_eyetrack_pro3 = subplot(5,2,7:8);
    plot(t_eyetrack_pro8.time,t_eyetrack_pro8.EyeA_Area,'Color','green');
    hold on 
    plot(t_eyetrack_pro8.time,t_eyetrack_pro8.EyeB_Area,'Color','blue');
    legend('EyeA Area','EyeB Area');
    hold off
    %Changing lims to represent the first trial
    ax_eyetrack_pro3.XLim = [t_events.onset(curr_idx)-4,t_events.onset(curr_idx)+t_events.duration(curr_idx)+6];
    ax_eyetrack_pro3.Title.String = sprintf("After Enzos Algorithm (AR=%s, Blink-time: %ds)",t_eyetrack_pro8.Properties.CustomProperties.AR,t_eyetrack_pro7.Properties.CustomProperties.MaxBlinkTime);
    ax_eyetrack_pro3.XTick = 0;
    ylabel("[a.u]");
end
%% Button: Create buttons and set callback functions
button_load = uicontrol("Style","pushbutton","String","Load New","Position",[0 0 70 20],"Callback",@loadPushButton);
button_prev = uicontrol("Style","pushbutton","String","Previous","Position",[70 0 60 20],"Callback",@prevPushButton);
button_continue = uicontrol("Style","pushbutton","String","Next","Position",[130 0 60 20],"Callback",@nxtPushButton);
button_print = uicontrol("Style","pushbutton","String","Print","Position",[190 0 60 20],"Callback",@printPushButton);
max_idx = numel(t_events.trial_index(:,1));
%% Define Callback function for "Load New"-Button
function loadPushButton(~,~                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               )
    [select_file,select_path] = uigetfile('*_physio.tsv.gz','Select a physio file','/media/diskEvaluation/Evaluation/sfb1280a05study7/rawdata/');
    if select_file == 0
        warndlg('You have no selected a file.','Warning','modal')
        return
    end
    f_fn_physio = fullfile(select_path,select_file);
    f_fn_events = strrep(f_fn_physio,"_physio.tsv.gz","_events.tsv");
    f_fn_eyetrack = strrep(f_fn_physio,"_physio.tsv.gz","_eyetrack.tsv.gz");
    if ~exist(f_fn_events,"file") || ~exist(f_fn_eyetrack,"file")
        warndlg('There we no corresponding eyetrack and event datafiles found.','Warning','modal')
        return
    else
        fu_bf = essbids_parseLabel(f_fn_physio);
        msgbox(sprintf("Loading sub-%s ses-%s run-%s.",fu_bf.sub,fu_bf.ses,fu_bf.run));
        close(hFig);
        func_plotPhysioByTrial(f_fn_physio,f_fn_events,f_fn_eyetrack,figure_setting);
    end
end
%% Define Callback function for "Previous"-Button
function prevPushButton(~,~)

        if (curr_idx-1) < 1
            warndlg("The current trial selected is the last trial of this session.",'Warning','modal')
            return
        end
        curr_idx = curr_idx-1;
        imgTrial = t_events.stim_file{curr_idx};
        if figure_setting == 1
            ax_resp.Title.String = sprintf('sub-%s  ses-%s  run-%s  trial-',bf.sub,bf.ses,bf.run)+string(t_events.trial_index(curr_idx));
            ax_events.XLim = [t_events.onset(curr_idx)-4,t_events.onset(curr_idx)+t_events.duration(curr_idx)+6];
            ax_eyetrack.XLim = [t_events.onset(curr_idx)-4,t_events.onset(curr_idx)+t_events.duration(curr_idx)+6];
            ax_eda.XLim = [t_events.onset(curr_idx)-4,t_events.onset(curr_idx)+t_events.duration(curr_idx)+6];
            ax_pulse.XLim = [t_events.onset(curr_idx)-4,t_events.onset(curr_idx)+t_events.duration(curr_idx)+6];
            ax_resp.XLim = [t_events.onset(curr_idx)-4,t_events.onset(curr_idx)+t_events.duration(curr_idx)+6];
        elseif figure_setting == 2
            ax_eyetrack_raw.Title.String = sprintf('sub-%s  ses-%s  run-%s  trial-',bf.sub,bf.ses,bf.run)+string(t_events.trial_index(curr_idx));
            ax_eyetrack_raw.XLim = [t_events.onset(curr_idx)-4,t_events.onset(curr_idx)+t_events.duration(curr_idx)+6];
            ax_eyetrack_pro1.XLim = [t_events.onset(curr_idx)-4,t_events.onset(curr_idx)+t_events.duration(curr_idx)+6];
            ax_eyetrack_pro2.XLim = [t_events.onset(curr_idx)-4,t_events.onset(curr_idx)+t_events.duration(curr_idx)+6];
            ax_eyetrack_pro3.XLim = [t_events.onset(curr_idx)-4,t_events.onset(curr_idx)+t_events.duration(curr_idx)+6];
            ax_events.XLim = [t_events.onset(curr_idx)-4,t_events.onset(curr_idx)+t_events.duration(curr_idx)+6];
        elseif figure_setting == 3
            updateScreenPlot(imgTrial);
            ax_events.XLim = [t_events.onset(curr_idx)-4,t_events.onset(curr_idx)+t_events.duration(curr_idx)+6];
            ax_eyetrack.XLim = [t_events.onset(curr_idx)-4,t_events.onset(curr_idx)+t_events.duration(curr_idx)+6];
        end

end
%% Define Callback function for "Next"-Button
function nxtPushButton(~,~)
        wHand = []; 
        if curr_idx == max_idx+1
            warndlg('The current trial selected is the last trial of this session.','Warning','modal');
            return
        end
        imgTrial = t_events.stim_file{curr_idx};
        if figure_setting == 1
            ax_resp.Title.String = sprintf('sub-%s  ses-%s  run-%s  trial-',bf.sub,bf.ses,bf.run)+string(t_events.trial_index(curr_idx));
            ax_events.XLim = [t_events.onset(curr_idx)-4,t_events.onset(curr_idx)+t_events.duration(curr_idx)+6];
            ax_eyetrack.XLim = [t_events.onset(curr_idx)-4,t_events.onset(curr_idx)+t_events.duration(curr_idx)+6];
            ax_eda.XLim = [t_events.onset(curr_idx)-4,t_events.onset(curr_idx)+t_events.duration(curr_idx)+6];
            ax_pulse.XLim = [t_events.onset(curr_idx)-4,t_events.onset(curr_idx)+t_events.duration(curr_idx)+6];
            ax_resp.XLim = [t_events.onset(curr_idx)-4,t_events.onset(curr_idx)+t_events.duration(curr_idx)+6];
        elseif figure_setting == 2
            ax_eyetrack_raw.Title.String = sprintf('sub-%s  ses-%s  run-%s  trial-',bf.sub,bf.ses,bf.run)+string(t_events.trial_index(curr_idx));
            ax_eyetrack_raw.XLim = [t_events.onset(curr_idx)-4,t_events.onset(curr_idx)+t_events.duration(curr_idx)+6];
            ax_eyetrack_pro1.XLim = [t_events.onset(curr_idx)-4,t_events.onset(curr_idx)+t_events.duration(curr_idx)+6];
            ax_eyetrack_pro2.XLim = [t_events.onset(curr_idx)-4,t_events.onset(curr_idx)+t_events.duration(curr_idx)+6];
            ax_eyetrack_pro3.XLim = [t_events.onset(curr_idx)-4,t_events.onset(curr_idx)+t_events.duration(curr_idx)+6];
            ax_events.XLim = [t_events.onset(curr_idx)-4,t_events.onset(curr_idx)+t_events.duration(curr_idx)+6];
        elseif figure_setting == 3
            updateScreenPlot(imgTrial);
            ax_events.XLim = [t_events.onset(curr_idx)-4,t_events.onset(curr_idx)+t_events.duration(curr_idx)+6];
            ax_eyetrack.XLim = [t_events.onset(curr_idx)-4,t_events.onset(curr_idx)+t_events.duration(curr_idx)+6];
        end
        curr_idx = curr_idx +1;
end
%% Define Callback function for "Print"-Button
function printPushButton(~,~)
    dir_save = uigetdir('sfb1280a05study7','Please select a folder to save to.');
    if dir_save == 0
        warndlg('You have not selected a folder.','No folder selected','modal')
        return
    else
        try 
            set(hFig,'Units','centimeters','Position',[-400 0 19 27]);
            set(button_load,'Visible','off');
            set(button_print,'Visible','off');
            set(button_continue,'Visible','off');
            set(button_prev,'Visible','off');
            saveas(gcf,sprintf('%s/sub-%s_ses-%s_run-%s_trial-%s.pdf',dir_save,bf.sub,bf.ses,bf.run,string(t_events.trial_index(curr_idx))));
            set(hFig,'Units','normalized','OuterPosition',[0.1 0 1 1]);
            set(button_load,'Visible','on');
            set(button_print,'Visible','on');
            set(button_continue,'Visible','on');
            set(button_prev,'Visible','on');

        catch
            errordlg('The file could not be saved.','Error: Save file.')
            return
        end
        msgbox('File was saved. Path: ' + string(dir_save));
    end
end
%% Use Enzos Algorithm to clean data
function tableProcessedEyetrack = cleanEyeTrackData(tableRawEyetrack,ti_eyetrack,val_AR,val_MaxBlinkTime)
% Process eyetracking data to check for blinks
% Copied from the script_ETanalyze by Enzo
tableProcessedEyetrack = tableRawEyetrack;
et_sampling_freq = 1000/tableRawEyetrack.Properties.CustomProperties.JsonSidecar.SamplingFrequency;
tableRawEyetrack.EyeA_AspectRatio = tableRawEyetrack.EyeA_PupilHeight ./ tableRawEyetrack.EyeA_PupilWidth;
tableRawEyetrack.EyeB_AspectRatio = tableRawEyetrack.EyeB_PupilHeight ./ tableRawEyetrack.EyeB_PupilWidth;

mean_AR_EA = mean(tableRawEyetrack.EyeA_AspectRatio,'omitnan');
std_AR_EA = std(tableRawEyetrack.EyeA_AspectRatio,'omitnan');

mean_AR_EB = mean(tableRawEyetrack.EyeB_AspectRatio,'omitnan');
std_AR_EB = std(tableRawEyetrack.EyeB_AspectRatio,'omitnan');

calc_AR_EA = mean_AR_EA - (2*std_AR_EA);
calc_AR_EB = mean_AR_EB - (2*std_AR_EB);

if val_AR > 0
    calc_AR_EA = val_AR;
    calc_AR_EB = val_AR;
end

tableRawEyetrack.EyeA_Blink = tableRawEyetrack.EyeA_AspectRatio < calc_AR_EA;
tableRawEyetrack.EyeB_Blink = tableRawEyetrack.EyeB_AspectRatio < calc_AR_EB;


tableRawEyetrack.EyeA_exclu_log = conv(tableRawEyetrack.EyeA_AspectRatio < calc_AR_EA, true(round(.5*et_sampling_freq), 1), 'same') > 0;
tableRawEyetrack.EyeB_exclu_log = conv(tableRawEyetrack.EyeB_AspectRatio < calc_AR_EB, true(round(.5*et_sampling_freq), 1), 'same') > 0;

t_exclusion_EyeA = func_signal2table(tableRawEyetrack.EyeA_exclu_log, et_sampling_freq, ti_eyetrack(1) * 1000);
t_exclusion_EyeB = func_signal2table(tableRawEyetrack.EyeB_exclu_log, et_sampling_freq, ti_eyetrack(1) * 1000);

t_blink_EyeA = t_exclusion_EyeA(t_exclusion_EyeA.duration < val_MaxBlinkTime * 1000,:);
t_blink_EyeB = t_exclusion_EyeB(t_exclusion_EyeB.duration < val_MaxBlinkTime * 1000,:);

log_EyeA_blink = func_table2signal(t_blink_EyeA,et_sampling_freq, ti_eyetrack*1000);
log_EyeB_blink = func_table2signal(t_blink_EyeB,et_sampling_freq, ti_eyetrack*1000);
log_Eyes_blink = log_EyeA_blink & log_EyeB_blink;

log_EyeA_OOF = tableRawEyetrack.EyeA_exclu_log(:);
log_EyeA_OOF(log_Eyes_blink) = false;
log_EyeB_OOF = tableRawEyetrack.EyeB_exclu_log(:);
log_EyeB_OOF(log_Eyes_blink) = false;

tableProcessedEyetrack.EyeA_Area(log_EyeA_OOF) = NaN;
tableProcessedEyetrack.EyeB_Area(log_EyeB_OOF) = NaN;
tableProcessedEyetrack.EyeA_Area(log_Eyes_blink) = NaN;
tableProcessedEyetrack.EyeB_Area(log_Eyes_blink) = NaN;

tableProcessedEyetrack = addprop(tableProcessedEyetrack,{'AR','MaxBlinkTime'},{'table','table'});
tableProcessedEyetrack.Properties.CustomProperties.AR = string(calc_AR_EA);
tableProcessedEyetrack.Properties.CustomProperties.MaxBlinkTime = val_MaxBlinkTime;
end
%% Create Eventmarkers out of the Events listed in the input
function tableEventmarkers = get_Eventmarkers(tableEvents)
    preal_a = zeros(numel(tableEvents(:,1))*4,3);
    for idx = 1:numel(tableEvents(:,1))
        idx_bf = (idx*4)-3;
        idx_on = (idx*4)-2;
        idx_off = (idx*4)-1;
        idx_nx = idx*4;
        c_trial_type = tableEvents.trial_type(idx);
         if isequal(c_trial_type{1},'CSminus')
            preal_a(idx_bf,:) = [(tableEvents.onset(idx)-0.001),0,tableEvents.trial_index(idx)];
            preal_a(idx_on,:) = [tableEvents.onset(idx),1,tableEvents.trial_index(idx)];
            preal_a(idx_off,:) = [tableEvents.onset(idx)+tableEvents.duration(idx),1,tableEvents.trial_index(idx)];
            preal_a(idx_nx,:) = [tableEvents.onset(idx)+tableEvents.duration(idx)+0.001,0,tableEvents.trial_index(idx)];
         elseif isequal(c_trial_type{1},'CSplus')
            preal_a(idx_bf,:) = [(tableEvents.onset(idx)-0.001),0,tableEvents.trial_index(idx)];
            preal_a(idx_on,:) = [tableEvents.onset(idx),1.1,tableEvents.trial_index(idx)];
            preal_a(idx_off,:) = [tableEvents.onset(idx)+tableEvents.duration(idx),1.1,tableEvents.trial_index(idx)];
            preal_a(idx_nx,:) = [tableEvents.onset(idx)+tableEvents.duration(idx)+0.001,0,tableEvents.trial_index(idx)];
         elseif isequal(c_trial_type{1},'US')
            preal_a(idx_bf,:) = [(tableEvents.onset(idx)-0.001),0,tableEvents.trial_index(idx)];
            preal_a(idx_on,:) = [tableEvents.onset(idx),1.3,tableEvents.trial_index(idx)];
            preal_a(idx_off,:) = [tableEvents.onset(idx)+tableEvents.duration(idx),1.3,tableEvents.trial_index(idx)];
            preal_a(idx_nx,:) = [tableEvents.onset(idx)+tableEvents.duration(idx)+0.001,0,tableEvents.trial_index(idx)];
         end
    end

    tableEventmarkers.timestamps = preal_a(:,1);
    tableEventmarkers.marker_value = preal_a(:,2);
    tableEventmarkers.trial_index = preal_a(:,3);
end
%% Update Screenview
function updateScreenPlot(trial_desc)
    axes(ax_screen);
    cla(ax_screen);
    log_trial = t_eyetrack.time > t_events.onset(curr_idx) - 4 & t_eyetrack.time < t_events.onset(curr_idx) + 6;
    if contains(trial_desc,'diamond')
        imshow(img_d);
        axis on;
    elseif contains(trial_desc,'square')
        imshow(img_s);
        axis on;
    end
    hold on
    scatter(t_eyetrack.EyeA_X_Gaze(log_trial),t_eyetrack.EyeA_Y_Gaze(log_trial),'green','filled','o','MarkerFaceAlpha',0.1);
    scatter(t_eyetrack.EyeB_X_Gaze(log_trial),t_eyetrack.EyeB_Y_Gaze(log_trial),'blue','filled','o','MarkerFaceAlpha',0.1);
    legend('Eye A','Eye B')
    ax_screen.Title.String = sprintf('sub-%s  ses-%s  run-%s  trial-',bf.sub,bf.ses,bf.run)+string(t_events.trial_index(curr_idx));
    ax_screen.XLim = [-200 2180];
    ax_screen.YLim = [-200 1280];
    hold off
end
end