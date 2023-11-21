function dHighres_click_batch(fullFiles,fullLabels,p,...
    tfFullFile, encounterTimes)

N = length(fullFiles);
previousFs = 0; % make sure we build filters on first pass
for idx1 =  1:N % for each data file
    fprintf('beginning file %d of %d \n',idx1,N)
    %(has to be inside loop for parfor, ie, filters are rebuilt every time,
    % can be outside for regular for)
    
    recFile = fullFiles{idx1};
    labelFile = fullLabels{idx1};
    
    % read file header
    hdr = dInput_HR_files(recFile,p);
    
    if isempty(hdr)
        continue % skip if you couldn't read a header
    elseif hdr.fs ~= previousFs
        % otherwise, if this is the first time through, build your filters,
        % only need to do this once though, so if you already have this
        % info, this step is skipped
        [previousFs,p] = dBuild_filters(p,hdr.fs);
        
        p = interp_tf(p,tfFullFile);
        if ~isfield(p,'countThresh') || isempty(p.countThresh)
            p.countThresh = 10^((p.dBpp - median(p.xfrOffset))/20)*(1/2);
        end
    end
    
    if exist(labelFile,'file')
        % Read in the .c file produced by the short term detector.
        [starts,stops] = ioReadLabelFile(labelFile);
    else
        continue
    end
    % Open xwav file
    fid = fopen(recFile, 'r');
    
    % Look for clicks, hand back parameters of retained clicks
    [cParams,f] = dProcess_HR_starts(fid,starts,stops,...
        p,hdr,recFile,labelFile);
    
    % Done with that file
    fclose(fid);
    fclose all;
    fprintf('done with %s\n', recFile);
    
    % Run post processing to remove rogue loner clicks, prior to writing
    % the remaining output files.
    clickTimes = sortrows(cParams.clickTimes);
    
    keepFlag = clickInlinePProc(labelFile,clickTimes,p,hdr,encounterTimes);
    keepIdx = find(keepFlag==1);
    
    % save a mat file now, rather than recalculating later
    cParams.clickTimes = clickTimes(keepIdx,:);
    cParams.ppSignalVec = cParams.ppSignalVec(keepIdx,:);
    cParams.durClickVec = cParams.durClickVec(keepIdx,:);
    cParams.bw3dbVec = cParams.bw3dbVec(keepIdx,:);
    
    cParams.specClickTfVec = cParams.specClickTfVec(keepIdx,:);
    if p.saveNoise
        cParams.specNoiseTfVec = cParams.specNoiseTfVec(keepIdx,:);
    end
    cParams.peakFrVec = cParams.peakFrVec(keepIdx,:);
    cParams.deltaEnvVec = cParams.deltaEnvVec(keepIdx,:);
    cParams.nDurVec = cParams.nDurVec(keepIdx,:);
    
    if ~isempty(keepIdx)
        cParams.yFiltVec = cParams.yFiltVec(keepIdx);
        cParams.yFiltBuffVec = cParams.yFiltBuffVec(keepIdx);
        if p.saveNoise
            cParams.yNFiltVec = cParams.yNFiltVec(keepIdx);
        end
    else
        cParams.yFiltVec = {};
        cParams.yFiltBuffVec = {};
        cParams.yNFiltVec = {};
    end
    
    save_dets2mat(strrep(labelFile,'.c','.mat'),cParams,f,hdr,p);
end