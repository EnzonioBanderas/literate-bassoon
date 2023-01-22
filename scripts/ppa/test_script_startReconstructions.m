%% Initialize Package Paths
fp_s = '/media/diskEvaluation/Evaluation/sfb1280a05study7/sourcedata/';
addpath /media/diskEvaluation/Evaluation/sfb1280a05study7/scripts/packages/spm12
addpath /media/diskEvaluation/Evaluation/sfb1280a05study7/scripts/packages/mapVBVD_20150918withFatNavs
addpath(genpath('/media/diskEvaluation/Evaluation/sfb1280a05study7/scripts/packages/retroMoCoBox/'));

%% Create Reconstruction Logs directories for the current session
fp_curr = cd;
fp_er_logs = fullfile(cd,'reconstruction_logs');
if ~exist(fp_er_logs,"dir")
    mkdir(fp_er_logs)
end
mkdir(fullfile(fp_er_logs,string("logs_" + string(datetime('now','Format','d-MM-y_HH-mm')))))
fp_curr_logs = fullfile(fp_er_logs,string("logs_" + string(datetime('now','Format','d-MM-y_HH-mm'))));

cd(fp_curr_logs);
diary command_line_logs;
diary on;
cd(fp_curr);
%% Get file lists
%  Script looks for all .dat files that contain the prefix sub-
%  Create a table (tab_recons) that will later hold information on the MID
%  path and whether this reconstruction should be run.
fl_raw = dir(fullfile(fp_s,'**','sub-*.dat')); %creates file list
tab_recons = table; %creates the table for later

%% Aquiring information about files contained in the file list
%  There are two loops: First this script loops through the previously 
%  definied file list (fl_raw) and the second will gather information on
%  each element that is in the same folder as the currently indexed file.
%  Goal: Find the values that are needed for reconstruction
for idx = 1:numel(fl_raw)
    fl_subfold = dir(fl_raw(idx).folder); %creates file list of the subfolder of the currently open file
    fol_name = "MPRAGErecon_" + string(extractBetween(fl_raw(idx).name,"meas_","_FID"));
    %fol_name is a string that mimics the name of the folder that will be
    %creatd during reconstruction. This allows for easy comparason if that
    %foldername does already exist or not and thus reconstruction has been
    %run
    
    fol_exists = zeros(1,numel(fl_subfold)); %preallocating space before foreloop
    for i = 1:numel(fl_subfold) %looping over each file/folder where the file is locating
        if strcmp(fol_name,fl_subfold(i).name) %comparing if the current element is the reconstruction folder
            fol_exists(i) = true; %assigning true into the preallocated array if thats the case
        else
            fol_exists(i) = false; %assing false if thats not the case
        end
    end

    warnstate = warning('query','MATLAB:table:RowsAddedExistingVars');
    warning('off','MATLAB:table:RowsAddedExistingVars'); 
    %assigning information gathered in the previous step to the table
    tab_recons.MID(idx) = fol_name; 
    tab_recons.PATH(idx) = {[fl_raw(idx).folder filesep fl_raw(idx).name]}; %fullpath!

    %checking whether reconstruction should take place or not by using nnz
    %which will count the numbers that are not zero (false) in the array
    %fol_exists. If there already was a reconstruction folder the value
    %will be 1 or higher or thus i assign true to show that the
    %reconstruction has already taken place
    if nnz(fol_exists) > 0 || contains(fl_raw(idx).folder,'sourcedata/mri_rawdata') %optional or condition for not sorted files
        tab_recons.RECON(idx) = true;
    else
        tab_recons.RECON(idx) = false;
    end
    warning(warnstate.state,'MATLAB:table:RowsAddedExistingVars');
end


%% Starting the reconstruction for the files
%  Now looping over the table that holds the values relevant for
%  reconstruction: Path and whether I want reconstruction to take place
for idx = size(tab_recons,1):-1:1
    if tab_recons.RECON(idx) == true
        %if the file has been reconstructed OR there is a reconstruction
        %run by this script just go to the next file (continue)
        disp("Reconstruction already found or file has not been sorted for: " + tab_recons.MID(idx))
        continue
    %elseif isfile(fullfile(fp_er_logs,"recon_in_progress.txt"))
       % disp("There is already a reconstruction in process.")
       % continue
    end
    starttime = string(datetime("now"));
    fid = fopen(fullfile(fp_er_logs,"recon_in_progress.txt"),"w");
    fprintf(fid,"Reconstruction started for: " + tab_recons.MID(idx) + "\n"); %add a datestring function now (isostandard)
    fprintf(fid,"Reconstruction started on the: " + starttime + "\n");
    fprintf(fid,"Filepath: " + tab_recons.PATH(idx));
    fclose(fid);
    %open a file to indiciate other users that somebody is running a
    %reconstruction
    fp_recon = convertStringsToChars(tab_recons.PATH(idx)); %use brace indexing here to then not use it in line 77
    disp("Reconstruction started for: " + tab_recons.MID(idx))
    try
        reconstructSiemensMP2RAGEwithFatNavs(fp_recon{1})
    catch ME
        fer_id = fopen(fullfile(fp_curr_logs,"ERROR_LOG_" + tab_recons.MID(idx)),"w");
        fprintf(fer_id,"Could not complete reconstruction for " + tab_recons.MID(idx) + "\n");
        fprintf(fer_id,"Filepath: " + tab_recons.PATH(idx) + "\n");
        fprintf(fer_id,"Starttime: " + starttime + "\n");
        fprintf(fer_id,"Error occured: " + string(datetime("now")) + "\n");
        fprintf(fer_id, "Error message: " + ME.message);
        fclose(fer_id);
        delete(fullfile(fp_curr_logs,"ERROR_LOG_" + tab_recons.MID(idx)));
        continue
    end
    fsc_id = fopen(fullfile(fp_curr_logs,"SUCESS_LOG_" + tab_recons.MID(idx)),"w");
    fprintf(fsc_id,"Sucessfully completed reconstruction for " + tab_recons.MID(idx) + "\n");
    fprintf(fsc_id,"Starttime: " + starttime + "\n");
    fprintf(fsc_id,"Completed at: " + string(datetime("now")) + "\n");
    fprintf(fsc_id,"Filepath: " + tab_recons.PATH(idx));
    delete(fullfile(fp_curr_logs,"ERROR_LOG_" + tab_recons.MID(idx)));
    %deleting the file to allow for a new reconstruction
 
end

cd(fp_curr_logs);
diary off;
cd(fp_curr);