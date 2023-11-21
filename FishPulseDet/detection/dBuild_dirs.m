function [metaDir,storeDir] = dBuild_dirs(baseDir, outDir)
% build output directories

if ~isempty(outDir) % use outDir if specified
    metaDir = ([outDir,'metadata']);
    storeDir = outDir;
else  % otherwise use baseDir
    metaDir = ([baseDir,'metadata']);
    storeDir = baseDir;
end

if ~isdir(metaDir)
    mkdir(metaDir)
end