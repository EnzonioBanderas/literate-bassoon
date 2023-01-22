%%
script_init_study7;

%% Get all physio files
fl_physio = dir(fullfile(fp_d,"**","*_physio.tsv.gz"));


%%Figure Options
a_options = ["Default","Eyetracking Comparason","Screen Preview"];
%%Preallocate space for the information gathered in the for loop
prealloc_a = strings(numel(fl_physio),5);
t_options = table();
for idx = 1:numel(fl_physio)
    %Locate equally named files with different suffixes, check if they
    %exist. If not move on to the next found physio file.
    fn_physio = fullfile(fl_physio(idx).folder,fl_physio(idx).name);
    bf = essbids_parseLabel(fn_physio);
    fn_events = strrep(fn_physio,"_physio.tsv.gz","_events.tsv");
    fn_eyetrack = strrep(fn_physio,"_physio.tsv.gz","_eyetrack.tsv.gz");
    if ~exist(fn_events,"file") || ~exist(fn_eyetrack,"file")
        continue
    end
    %% Parse information on session and run into descriptions of the trial.
    if bf.ses == '1' && bf.run == '1'
        trial_desc = "Habituation";
    elseif bf.ses == '1' && bf.run == '2'
        trial_desc = "Acquisition";
    elseif bf.ses == '2' && bf.run == '1'
        trial_desc = "Extinction";
    elseif bf.ses == '3' && bf.run == '1'
        trial_desc = "Recall";
    elseif bf.ses == '3' && bf.run == '2'
        trial_desc = "Volatile";
    end
    subject_num = erase(string(bf.sub),"Z7T7");

    %% Add gathered data into preallocated string array
    prealloc_a(idx,1) = subject_num;
    prealloc_a(idx,2) = trial_desc;
    prealloc_a(idx,3) = fn_physio;
    prealloc_a(idx,4) = fn_events;
    prealloc_a(idx,5) = fn_eyetrack;
end
%% Create a table that contains the data gathered in string array.
t_options.subject_id = str2double(prealloc_a(:,1));
t_options.trial_type = prealloc_a(:,2);
t_options.fn_physio = prealloc_a(:,3);
t_options.fn_events = prealloc_a(:,4);
t_options.fn_eyetrack = prealloc_a(:,5);
% Remove empty values
t_options(any(ismissing(t_options),2),:) = [];

%% Start of Menu
dir_choice = uigetdir('/media/diskEvaluation/Evaluation/sfb1280a05study7/rawdata/','Please select Participant folder');
if isequal(dir_choice,0)
    input_subject = inputNumber('Please select a subject: Z7T7');
else
    input_subject = string(extractAfter(dir_choice,'sub-Z7T7'));
    if strlength(input_subject) > 3
        input_subject = extractBefore(input_subject,filesep);
    end
    input_subject = double(input_subject);
end


input_choice = selectRunMenu(t_options(t_options.subject_id == input_subject,:),input_subject);
t_choice = t_options(t_options.subject_id == input_subject,:);
if input_choice ~= 0
    input_figure = figureMenu(a_options);
else
    clear;
    script_plotPhysioByTrial;
end
t_choice = t_choice(input_choice,:);
for idx=1:30
    fprintf('\n');
end

fprintf('Loading plots for sub-Z7T7%d. Selected run: %s\n',input_subject,t_choice.trial_type(1));
fprintf('Selected figure type: %s\n',a_options(input_figure));
func_plotPhysioByTrial(t_choice.fn_physio(1),t_choice.fn_events(1),t_choice.fn_eyetrack(1),input_figure);

clear;
function num = inputNumber(prompt)
    while true
        num = str2double(input(prompt,"s"));
        if ~isnan(num)
            break
        end
    end
end

function choice = selectRunMenu(options,subject_num)
    if isempty(options)
        fprintf('There were no sessions found for that subject.\n');
        choice = 0;
        return
    end
    for idx = 1:30
        fprintf('\n')
    end
    fprintf('The following runs were found for sub-Z7T7%d\n',subject_num)
    for idx=1:height(options)
    fprintf('(%d) %s\n',idx,options.trial_type(idx));
    end
    fprintf('Please select a run by entering the corresponding number.\n')
    choice = 0;
    while ~any(choice == 1:height(options))
        choice = inputNumber('>> ');
    end
end

function figure_setting = figureMenu(options)
    for idx = 1:10
        fprintf('\n')
    end
    if isempty(options)
        fprintf('There seems to be a problem with the figure type selection. Loading: Default.\n');
        figure_setting = 1;
        return
    end
    fprintf('The following figure options are available.\n')
    for idx=1:length(options)
    fprintf('(%d) %s\n',idx,options(idx));
    end
    fprintf('Please select a figure type by entering the corresponding number.\n')
    figure_setting = 0;
    while ~any(figure_setting == 1:length(options))
        figure_setting = inputNumber('>> ');
    end
end

