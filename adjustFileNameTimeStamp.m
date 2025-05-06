% Change file name timestamps by adding/subtracting hours; expects
% timestamps to be in 'yyyymmdd_HHMMSS' format
% Directory containing file to be renamed
inDir = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_WAVE\Stockport\ST300\SCh03';
% Where to save renamed files; if same as inDir, files will be overwritten
outDir = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_WAVE\Stockport\ST300\SCh03';
% Wildcard to identify files to rename
fileExt = '*.wav';
% Regular expression for time stamp in file names
timeFormat = '\d{8}_\d{6}';
tzAdjust = -4; % Hours to add to current time stamp

%%
fileList = dir(fullfile(inDir,fileExt));
fileNames = cellstr(vertcat(fileList.name));
timeStampsOld = regexp([fileNames{:}],timeFormat,'match')';
timeStampsNew = datetime(timeStampsOld,'InputFormat','uuuuMMdd_HHmmss') + hours(tzAdjust);
timeStr = cellstr(datestr(timeStampsNew,'yyyymmdd_HHMMSS'));

[startInd,endInd] = cellfun(@(x) regexp(x,timeFormat),fileNames);

newNames = cellfun(@(x,y,z) strrep(x,y,z),fileNames,timeStampsOld,timeStr,'UniformOutput',false);

for i=2:length(fileList)
    movefile(fullfile(inDir,fileNames{i}),fullfile(outDir,newNames{i}));
end