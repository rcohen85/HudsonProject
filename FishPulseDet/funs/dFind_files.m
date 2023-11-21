function fileNames = dFind_files(baseDir,depl,metaDir,fileExt)
% also returns .wav files

% Find folders in baseDir
folders = dir(baseDir);
trueIdx = [];
% Remove those that don't belong to data
for fidx = 1:length(folders)
    true = strfind(folders(fidx).name, depl);
    decim = strfind(folders(fidx).name, 'd100');
    other = strfind(folders(fidx).name, 'other');
    % (you can use "other" in name to avoid scanning for xwavs)
    if isempty(true) || ~isempty(decim) || ~isempty(other)
        trueIdx(fidx) = 0;
    else
        trueIdx(fidx) = 1;
    end
end

keep = find(trueIdx==1);

% Build file structure or use guided detection spreadsheet to identify files of
% interest.
folderNames = {};
m = 1;
for fidx = 1:length(keep)
    if isdir(fullfile(baseDir,folders(keep(fidx)).name)) == 1
        folderNames{m,1} = char(folders(keep(fidx)).name);
        m = m+1;
    end
end

% Pull out files from all folders, combine full paths into one long list
fileNames = [];
for fidx = 1:size(folderNames,1)
    fileDir = fullfile(baseDir,folderNames{fidx,1});
    metaSubDir = fullfile(metaDir,folderNames{fidx,1});
    if ~isdir(metaSubDir)
        mkdir(metaSubDir)
    end
    % list of files
    d = dir(fullfile(fileDir,['*',fileExt])); % list of files
    files = char(d.name);      % file names in directory
    
    % make full path filenames
    fileList = [];
    for s = 1:size(files,1)
        fileList(s,:) = fullfile(baseDir,folderNames{fidx,1},files(s,:));
    end
    fileNames = [fileNames;char(fileList)];
end

