function [fullFiles,fullLabels] = get_fileset(baseDir,metaDir,detFiles,fileExt)
% Make list of what you're going to name your output files, for easy
% reference later.
fullFiles = []; % 
fullLabels = []; % .c files

for f2 = 1:size(detFiles,1)
    thisFile = detFiles(f2,:);
    fullFiles{f2}= thisFile;
    [pathStr, thisName, ext] = fileparts(thisFile);
    [~,subDir] = fileparts(pathStr);
    thisName2 = [thisName,ext];
    if strfind(thisName2,'.x.wav')
        thisLabel = strrep(thisName2,'.x.wav','.c');
    else
        thisLabel = strrep(thisName2,fileExt,'.c');
    end
    fullLabels{f2} = fullfile(metaDir,subDir,thisLabel);
end
