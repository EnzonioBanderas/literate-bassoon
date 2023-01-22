script_init_study7

path_to_file = mfilename("fullpath") % path to this codefile 

fp_scr = fileparts(mfilename("fullpath")) %path to dir where this codefile is

fp0 = fileparts(fp_scr)% path to one dir out '/media/diskEvaluation/Evaluation/sfb1280a05study7'

path_t= mfilename("class")

fp_d = fullfile(fp0, 'rawdata') % path to the rawdata dir in the fpO dir '/media/diskEvaluation/Evaluation/sfb1280a05study7/rawdata'

supertab = essbids_readTsv([fp_d filesep 'participants.tsv'])

fp_m = fullfile(fp0, 'misc') % path to the misc dir in the fp0 dir '/media/diskEvaluation/Evaluation/sfb1280a05study7/misc'

path_to_tab = fullfile(fp_m,"SFB-Z7T7-List.xlsx") % path to tab in fp_m dir /media/diskEvaluation/Evaluation/sfb1280a05study7/misc/SFB-Z7T7-List.xlsx

read_this_tab = readtable(path_to_tab)

two_columns_tab =read_this_tab(:,["Z7T7_Code" "SFB_Code"])

three_columns_tab =read_this_tab(:,["Z7T7_Code" "Comment" "SFB_Code"])

three_raws_three_columns_tab =read_this_tab([7 18 23], ["Z7T7_Code" "Comment" "SFB_Code"]) 

three_raws_four_columns_tab =read_this_tab([7 18 23], ["Nr_" "Z7T7_Code" "Comment" "SFB_Code"]) % ??why the number column has name Nr_ instead of Nr.

renamed_twocolumns_tab =renamevars(two_columns_tab,["Z7T7_Code" "SFB_Code"],["incode" "outcode"])

fp_s = fullfile(fp0, 'sourcedata') %path to the sourcedata dir in the fp0 dir '/media/diskEvaluation/Evaluation/sfb1280a05study7/sourcedata'



%DONT UNDERSTAND HOW THIS WORKS FOR NOW
liststruct = dir(fullfile(fp_s,'**','*_setme.csv'))

listtable = struct2table(dir(fullfile(fp_s,'**','*_setstim.txt')))

liststuct2 = dir(fullfile(fp_s,'**','*_quest.xlxs'))

ID = strings(size(list,1),1)
CS = strings(size(list,1),1)
Stim = nan(size(list,1),1)
%DONT UNDERSTAND HOW THIS WORKS FOR NOW

for idx = 1:size(list,1)  
    
    fcur = fileread([liststruct(idx).folder filesep list(idx).name])
   
    f_participant_id = regexp (fcur,'ID\s*: (?<val>[a-zA-Z0-9\-]+)','names')
    
    array_ParticipantID(idx,1) = string(f_participant_id.val)
   
    f_cs_plus = regexp (fcur,'CSplus\s*: (?<val>[a-z]+)','names') 
    
    array_CSplus(idx,1) = string(f_cs_plus.val)

    f_stim_intensity = regexp (fcur,'US_mA\s*: (?<val>[0-9\.]+)','names')
    
    array_StimIntensity(idx,1) = double(string(f_stim_intensity.val))
end












