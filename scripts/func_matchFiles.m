function [fp1, fp2] = func_matchFiles(fl1, fl2, numMatch, outputFL)
% FUNC_MATCHFILES: given two BIDS file lists, find best file match
% TD:
%   1. Change positional input parameters to inputParser
%   2. Add fmap bool, which tells the function whether it should
%      concatenate task and acq in the non-fmap fl
%   3. Generalize function so that the fields over which it matches files
%      is an input in the form of a cell array of character vectors 
if ~iscell(fl1)
    fl1 = {fl1};
end
if ~iscell(fl2)
    fl2 = {fl2};
end
if ~exist('numMatch', 'var')
    numMatch = 14;
end
if ~exist('outputFL', 'var')
    outputFL = false;
end
% in case of an fmap, the fmap fl has an acq field which is a concatenation
% of task and acq fields, which means that the non-fmap field's acq field
% should be changed so that it is also a concatenation of the two

matchMat_sub = false(length(fl1), length(fl2));
matchMat_run = false(length(fl1), length(fl2));
matchMat_ses = false(length(fl1), length(fl2));
matchMat_acq = false(length(fl1), length(fl2));
for iFL1 = 1:length(fl1)
    fl1_parsed = essbids_parseLabel(fl1{iFL1});
    for iFL2 = 1:length(fl2)
        fl2_parsed = essbids_parseLabel(fl2{iFL2});
            
            if isfield(fl1_parsed, 'sub') && isfield(fl2_parsed, 'sub')
                if strcmp(fl1_parsed.sub, fl2_parsed.sub)
                    matchMat_sub(iFL1, iFL2) = true;
                end
            end
            if isfield(fl1_parsed, 'ses') && isfield(fl2_parsed, 'ses')
                if strcmp(fl1_parsed.ses, fl2_parsed.ses)
                    matchMat_ses(iFL1, iFL2) = true;
                end
            end
            if isfield(fl1_parsed, 'run') && isfield(fl2_parsed, 'run')
                if strcmp(fl1_parsed.run, fl2_parsed.run)
                    matchMat_run(iFL1, iFL2) = true;
                end
            end
            if isfield(fl1_parsed, 'acq') && isfield(fl2_parsed, 'acq')
                if strcmp(fl1_parsed.acq, fl2_parsed.acq)
                    matchMat_acq(iFL1, iFL2) = true;
                end
            end
    end
end

% If nothing matches, skip all other steps and return empty output
if ~any(matchMat_sub(:)) && ~any(matchMat_run(:)) && ...
   ~any(matchMat_ses(:)) && ~any(matchMat_acq(:))
    fp1=[]; fp2=[];
else

% If nothing matches for a field, essentially skip this field by setting all to true
% if ~any(matchMat_sub(:))
%     matchMat_sub(:) = true;
% end
% if ~any(matchMat_run(:))
%     matchMat_run(:) = true;
% end
% if ~any(matchMat_ses(:))
%     matchMat_ses(:) = true;
% end
if ~any(matchMat_acq(:))
    matchMat_acq(:) = true;
end
matchMat = matchMat_sub;
matchMat = matchMat & matchMat_ses;
matchMat = matchMat & matchMat_run;
matchMat = matchMat & matchMat_acq;

% find matches, if there is only one match this is the output!
[iFL1, iFL2] = find(matchMat);

% if there is no match, try to find the best match
if isempty(iFL1)
    matchMat_sum = matchMat_sub*8 + matchMat_ses*4 + matchMat_run*2 + matchMat_acq;
    maxMatch = max(matchMat_sum);
    if maxMatch >= numMatch % if there is a match
        [iFL1, iFL2] = find(matchMat_sum == maxMatch);
        if isempty(iFL1)
            error('no matches found, maxMatch if statement')
        end
    else
        fp1=[]; fp2=[];
        return
    end
end
    
%     [iFL1, iFL2] = find(matchMat);
%     if isempty(iFL1)
%         error('no match found')
%     else
%         iFL1 = iFL1(1);
%         iFL2 = iFL2(1);
%     end

% if there are still multiple matches, give warning and arbitrarily choose first match
if length(iFL1)>1
    if outputFL
        fp1 = fl1(iFL1);
        fp2 = fl2(iFL2);
        return
    else
        warning('multiple best matches found, arbitrarily choosing first best match')
        iFL1 = iFL1(1);
        iFL2 = iFL2(1);
    end
%     [iFL1_acq, iFL2_acq] = find(matchMat);
%     iFL_acq_logical = ismember([iFL1_acq, iFL2_acq], [iFL1, iFL2]);
%     if sum(iFL_acq_logical)==0
%         iFL1 = iFL1(1);
%         iFL2 = iFL2(1);
%     else
%         iFL1 = iFL1_acq(find(iFL_acq_logical, 1));
%         iFL2 = iFL2_acq(find(iFL_acq_logical, 1));
%     end
end

% this code is bad, define matches at start and then loop through different
% match levels, if no matches remain after a refinement go ahead and keep
% old matches, go through all levels, if multiple matches remain at the end
% just take the first match and give a warning, if at any point only one
% match remains, quit refining and take that match to be true.

% also code does not work for fmaps where task is added to acq

% 05-10-2022 update: Code might still be bad, but warning added in case of
% multiple best matches. Could be better to generalize and add a
% match_fields parameter which you loop over e.g. {'sub', 'ses', 'run', ...}

% Should this function be able to return multiple filenames instead?

fp1 = fl1{iFL1};
fp2 = fl2{iFL2};




% % find matching fmap
% for iFL1 = 1:length(fl1)
% for iFMAP = 1:length(fl_fmap)
%     fp_fmap = fullfile(fl_fmap(iFMAP).folder, fl_fmap(iFMAP).name);
%     fp_fmap_parsed = essbids_parseLabel(fp_fl);
%     try
%         if strcmp(fp_fmap_parsed.run, fp_func_preprocessed_parsed.run) && ...
%            strcmp(fp_fmap_parsed.acq, [fp_func_preprocessed_parsed.acq, ...
%                 fp_func_preprocessed_parsed.task])
%             disp('break');
%             break
%         end
%     catch
%     end
% end
end

end

