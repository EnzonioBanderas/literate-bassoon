function fp_info = func_niftiinfo(fp)
% same as niftiinfo, but gunzips file beforehand if niftiinfo fails
    try
        fp_info = niftiinfo(fp);
    catch
        warning('cannot unzip via niftiinfo, unzipping beforehand')
        fl = gunzip(fp);
        fp_info = niftiinfo(fl{1});
        S = recycle;
        recycle('off')
        delete(fl{1})
        recycle(S)
    end
end