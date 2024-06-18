% Script to get basic info about acoustic files; requires Triton to be on
% your path (for function wavname2dnum, to parse timestamps from file
% names)

% Directory containing sound files (wave or flac)
inDir = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_FLAC\Stockport\Swift\SCr02';
% Uncomment these two lines and give a deployment name if you want to
% collect information across multiple deployments in a single matrix (bigStats)
% bigStats = []; 
% depName = '';

% list files in the directory
fileList = dir(fullfile(inDir,'*.flac')); 
badFiles = [];
badFiles = find(vertcat(fileList.bytes)==0); % find any empty files
if ~isempty(badFiles)
sprintf('Warning: Empty files found. Skipping these\n')
goodFiles = setdiff(1:length(fileList),badFiles);
fileList = fileList(goodFiles); % remove empty files from list to be iterated over
end

% These are the things that will be extracted from/calculated for each file
% For a single deployment
fileStats = struct('FileName',{},'FileStart',[],'FileStop',[],'Fs',[],'NumChan',[],'Bits',[]); 
% when collecting info across multiple deployments
% fileStats = struct('Dep',{},'FileName',{},'FileStart',[],'FileStop',[],'Fs',[],'NumChan',[],'Bits',[]);

% Iterate over file list and get desired information about each file
for i = 1:length(fileList)
info = audioinfo(fullfile(inDir,fileList(i).name));
st = datetime(wavname2dnum(fileList(i).name),'ConvertFrom','datenum');

% fileStats(i).Dep = depName; % uncomment when collecting info across multiple deployments
fileStats(i).FileName = fileList(i).name;
fileStats(i).FileStart = st;
fileStats(i).FileStop = st + seconds(info.Duration);
fileStats(i).Fs = info.SampleRate;
fileStats(i).TotSamps = info.TotalSamples;
fileStats(i).NumChan = info.NumChannels;
fileStats(i).Bits = info.BitsPerSample;

end

% bigStats = [bigStats;struct2table(fileStats)]; % uncomment when collecting info across multiple deployments

figure(202) % Simple plot showing file start times; if files are expected to be evenly spaced 
% (e.g. one file every hour), this can be used to identify gaps
plot(datenum(vertcat(fileStats.FileStart)),vertcat(fileStats.Fs),'o');
datetick
