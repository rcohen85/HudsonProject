% Compare manually-annotated pulses to automatically detected ones

startEff = datetime(2021,11,19,12,10,0); % timestamp at which manual annotation effort began
endEff = datetime(2021,11,19,15,0,0); % timestamp at which manual annotation effort stopped
man = readtable('W:\projects\2022_NOAA-NERRS_HudsonNY_144488\ShortStur\144488Hudson_024K_SS02.ManualSelections.txt','Delimiter','tab');
det = readtable('W:\projects\2022_NOAA-NERRS_HudsonNY_144488\ShortStur\SelectionTables\SS02_SelectionTable_20211119_112835_20211119_231735.txt','Delimiter','tab');

% Remove duplicate entries (due to multiple sound views in Raven)
manWav = find(strcmp(table2array(man(:,'View')),'Waveform 1'));
detWav = find(strcmp(table2array(det(:,'View')),'Waveform 1'));
man = man(manWav,:);
det = det(detWav,:);

% Get timestamps of manually annotated pulses and automated detections within effort period
manInds = find(strcmp(table2array(man(:,'SoundType')),'Chirp')|strcmp(table2array(man(:,'SoundType')),'Squirt'));
manTimes = sort(datetime(datevec(table2array(man(manInds,'BeginDateTime')))));
manTimes = manTimes(manTimes>=startEff & manTimes<=endEff);

detInds = find(strcmp(table2array(det(:,'SoundType')),'Pulse'));
detDate_strs = regexp(table2array(det(detInds,'BeginFile')),'\d{8}[_]\d{6}','match');
detTimes = datetime(cell2mat(cellfun(@(x)datevec(x,'yyyymmdd_HHMMSS'),detDate_strs,'UniformOutput',false)));
detTimes = sort(detTimes+seconds(table2array(det(detInds,'FileOffset_s_'))));
detTimes = detTimes(detTimes>=startEff & detTimes<=endEff);

% Get rid of detections during boat passage
manBoat = (manTimes<=datetime(2021,11,19,13,59,0) | manTimes>=datetime(2021,11,19,14,7,0));
manTimes = manTimes(manBoat);

detBoat = (detTimes<=datetime(2021,11,19,13,59,0) | detTimes>=datetime(2021,11,19,14,7,0));
detTimes = detTimes(detBoat);

% Compare number of manual vs automated detections in short time bins
bins = startEff:minutes(1):endEff;
[NMan,bins,manBin] = histcounts(manTimes,bins);
[NDet,bins,detBin] = histcounts(detTimes,bins);
counts = [NMan',NDet'];

rho = corr(counts);
txt = sprintf('Rho = %0.3f',rho(1,2));

figure(9)
plot(bins(1:(end-1)),counts);
legend({'Manual Detections','Automated Detections'});
text(bins(floor(length(bins)*0.75)),max(max(counts))*0.95,txt);

figure(99)
bar(bins(1:(end-1)),counts);
legend({'Manual Detections','Automated Detections'});
text(bins(floor(length(bins)*0.75)),max(max(counts))*0.95,txt);

% Find detections which refer to the same pulse
pulseDiffs = TKTK;

