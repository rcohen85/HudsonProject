% Find files with no calls annotated and annotate as many noise clips as possible

audioDir = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_FLAC\AtlSturSpw\ST300'; % Path to audio files
fileExt = '.flac'; % file extension of audio files
selDir = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AtlStur\ThunderClassifier_Output\CleanedForTimeseries'; % Path to save selection tables
clip_len = 2.05; % duration of clips (s)
label = 'noise'; % desired noise label (e.g. 'noise', 'background', 'other')

%%

audioList = dir(fullfile(audioDir,['**\*',fileExt]));
audioPaths = {audioList.folder};
audioList = strrep(cellstr(vertcat(audioList.name)),fileExt,'');
selList = dir(fullfile(selDir,'*.txt'));
selList = strrep(cellstr(vertcat(selList.name)),'.txt','');

if ~isfolder(fullfile(selDir,'Noise'))
    mkdir(fullfile(selDir,'Noise'));
end

[noCallFiles, ia] = setdiff(audioList,selList);
noCallFiles = strcat(noCallFiles,fileExt);
audioPaths = audioPaths(ia);

for i=1:numel(noCallFiles)

    info = audioinfo(fullfile(audioPaths{i},noCallFiles{i}));
    endSec = repelem([clip_len:clip_len:(info.TotalSamples/info.SampleRate)]',info.NumChannels);
    startSec = endSec-clip_len;
    lowFreq = repelem(zeros(numel(startSec),1),numel(info.NumChannels));
    highFreq = repelem((info.SampleRate/2)*ones(numel(startSec),1),numel(info.NumChannels));
    fileOffset = startSec;
    beginFile = cellstr(repmat(noCallFiles{i},numel(startSec),1));
    tags = cellstr(repmat(label,numel(startSec),1));
    chan = repmat([1:info.NumChannels]',numel(startSec)/info.NumChannels,1);

    selTab = table(chan,startSec,endSec,lowFreq,highFreq,beginFile,fileOffset,tags,...
        'VariableNames',{'Channel','Begin Time (s)','End Time (s)','Low Freq (Hz)',...
        'High Freq (Hz)','Begin File','File Offset (s)','Tags'});

    saveName = fullfile(selDir,'Noise',strrep(noCallFiles{i},fileExt,'.txt'));
    writetable(selTab,saveName,'Delimiter','\t');

end