clearvars
addpath /opt/MATLAB/toolbox_ter/
%fp0 = '/media/diskEvaluation/Evaluation/sfb1280a05study6/derivatives/EDAevaluation';
fp0 = '/media/diskEvaluation/Evaluation/sfb1280a05study7/dumpHereForSorting/testEDAoutput/EDAevaluation';
flist = ter_listFiles(fullfile(fp0,''),'*_EDA.mat');

%remove already evaluated files
flist(contains(flist,"_EDA_EDA.mat"))=[];


 
for i=1:numel(flist)

edaFN = flist{i};


eventtable = loadEventTable(edaFN);
load (edaFN, 'data')
%we only need samplingrate from trialdefinition file, essentially this should be skipped

TDefFile = strrep(edaFN,"_EDA.mat","_trialDefinition.m");
if isfile(TDefFile) 
eval(fileread(fullfile(TDefFile))); 
else 
    disp('Trial Definition file not found for: ' + edaFN )
end

sPos = [];
trialStartMarkers = [];
trialStopMarkers = [];
EIRlabels = [];
EIRStartMarkers = [];
TIRlabels = [];
TIRStartMarkers = [];


for i=1:height(eventtable)
label = eventtable.event_label{i};
labelparts = split(label,"_");
suffix = "-NoUS";
      if contains(label, "Paired") 
              suffix = "-US";
         end
label = strcat(labelparts{1},'-',labelparts{2},suffix);

switch labelparts {1}
    case 'EIR'
        EIRlabels = [EIRlabels label];
        EIRStartMarkers = [EIRStartMarkers eventtable.event_onsets(i)];
        trialStartMarkers = [trialStartMarkers eventtable.event_onsets(i)];
       
    case 'TIR'
        TIRlabels = [TIRlabels label];
        TIRStartMarkers = [TIRStartMarkers eventtable.event_onsets(i)];
        trialStopMarkers = [trialStopMarkers eventtable.evalBlock_stop(i)];

end




end


%generating trial definition values


%marking important events
latencyEIR   = 1;  
latencyUS    = 0.5;

durationUS   = 4;
durationEIR  = 6 - latencyUS;

numTrials = numel (trialStartMarkers);
recordingSystem = 'GSR100C'; %needs to be checked
stimOnset = trialStartMarkers;
trialSep = [stimOnset trialStopMarkers(end)];
trialStart = [];
trialStop = [];


 


uniquelabels = horzcat(unique (EIRlabels), unique(TIRlabels));


interval_boundaries_lower = EIRStartMarkers + latencyEIR;

sPos{1}.color = 'c';
sPos{1}.name = 'intervalBoundary_l_o_w';
sPos{1}.time = interval_boundaries_lower;

  interval_boundaries_upper = EIRStartMarkers   +  latencyEIR + durationEIR + durationUS;

sPos{2}.color = 'g';
sPos{2}.name = 'intervalBoundary_u_p';
sPos{2}.time = interval_boundaries_upper;


%comment this loop out for blind evaluation, it marks specific CS/US window
%intervals

for j=1:numel(uniquelabels)
    currentLabel = uniquelabels(j);
    currentWindow = split(currentLabel,'-');
    switch currentWindow{1}
        case 'EIR'
            L = EIRStartMarkers;
            L(EIRlabels ~= currentLabel) = -1;
        case 'TIR'
            L = TIRStartMarkers;
            L(TIRlabels~=currentLabel) = -1;
    end
    sPos{j+2}.time = L;
    sPos{j+2}.color = 'b';
    sPos{j+2}.name = char(currentLabel);
end





tDEFfName = strrep(edaFN,"_EDA.mat","_trialDef.mat");
save(tDEFfName,'numTrials','recordingSystem','samplRate','sPos','stimOnset','trialSep','trialStart','trialStop');
end

function eventtable = loadEventTable (edaFilename)
disp ("Generating new TrialDefinition file for :" + edaFilename);
  if isfile(strrep(edaFilename,"_EDA.mat","_eventTable.mat"))
        FN = strrep(edaFilename,"_EDA.mat","_eventTable.mat");
        
    else
        [FN,PN] = uigetfile('*_eventTable.mat','Select event file');
    end
    
     try
        %load eventtable and remove all non-FIR events
        
        load(FN,'eventtable');
 %       [I, ~] = find(cellfun(@(s) ~contains(s, 'FIR'), eventtable{:,4}));
  %      eventtable(I, :) = [];
        
    catch ME
       disp ('ERROR: eventtable variable not found')
     end    
    

end

