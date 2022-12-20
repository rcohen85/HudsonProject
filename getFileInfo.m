% bigStats = []; 
inDir = 'W:\projects\2022_NOAA-NERRS(UMich)_HudsonNY_144488\AcousticData\BlackCreek\Swift\WAVE\BC05';
% depName = '';
fileList = dir(fullfile(inDir,'*.wav'));
badFiles = [];
badFiles = find(vertcat(fileList.bytes)==0);
if ~isempty(badFiles)
sprintf('Warning: Empty files found. Skipping these\n')
goodFiles = setdiff(1:length(fileList),badFiles);
fileList = fileList(goodFiles);
end
% fileStats = struct('Dep',{},'FileName',{},'FileStart',[],'FileStop',[],'Fs',[],'NumChan',[],'Bits',[]);
fileStats = struct('FileName',{},'FileStart',[],'FileStop',[],'Fs',[],'NumChan',[],'Bits',[]);

for i = 1:length(fileList)
info = audioinfo(fullfile(inDir,fileList(i).name));
st = datetime(wavname2dnum(fileList(i).name),'ConvertFrom','datenum');

% fileStats(i).Dep = depName;
fileStats(i).FileName = fileList(i).name;
fileStats(i).FileStart = st;
fileStats(i).FileStop = st + seconds(info.Duration);
fileStats(i).Fs = info.SampleRate;
fileStats(i).NumChan = info.NumChannels;
fileStats(i).Bits = info.BitsPerSample;

end

% bigStats = [bigStats;struct2table(fileStats)];
figure(202)
plot(datenum(vertcat(fileStats.FileStart)),vertcat(fileStats.Fs),'o');
datetick
