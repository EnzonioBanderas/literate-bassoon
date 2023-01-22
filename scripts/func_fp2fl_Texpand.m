function fl = func_fp2fl_Texpand(fp)
%FUNC_FP2FL_TEXPAND Expand image file path to image file path list of time
% series
%   Takes a character vector image file path input and reas the numbber of
%   images in the time series, then expands the file path to a file list
%   with an element for each slice in time of the image
fp_info = niftiinfo(fp);
if length(fp_info.ImageSize)<4
    nT = 1;
else
    nT = fp_info.ImageSize(4);
end
fl = repmat({fp}, [nT, 1]);
for iFL=1:length(fl)
    fl{iFL} = [fl{iFL}, ',', num2str(iFL)];
end

end

