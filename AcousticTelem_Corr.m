telemPath = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\ShortStur\esopus telemetry subset.csv';
autDir = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\ShortStur\SelectionTables\SS02\It_21';
manPath = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\ShortStur\SelectionTables\144488Hudson_024K_SS02.ManualSelections.txt';
ShStRec = 550786;
spec = 'shortnose sturgeon';
startEff = datetime(2021,11,19,12,0,0); % timestamp at which to begin comparison
endEff = datetime(2021,12,02,18,0,0); % timestamp at which manual annotation effort stopped
fileSetStart = datetime(2021,11,19,11,28,17);

%% Load data
% Automated detections
autFiles = dir(fullfile(autDir,'*.txt'));
autDet = [];
for i = 1:length(autFiles)
    tb = readtable(fullfile(autDir,autFiles(i).name),'Delimiter','tab','VariableNamingRule','preserve');
    autDet = [autDet;tb];
end

% Manual detections
manDet = readtable(manPath,'Delimiter','tab','VariableNamingRule','preserve');

% Remove duplicate detection entries (due to multiple sound views in Raven)
pulseWav = find(strcmp(table2array(autDet(:,'View')),'Waveform 1'));
autDet = autDet(pulseWav,:);
manWav = find(strcmp(table2array(manDet(:,'View')),'Waveform 1'));
manDet = manDet(manWav,:);

% Load telmetry data
telemDat = readtable(telemPath,'Delimiter',',','VariableNamingRule','preserve');

%% Wrangle telemetry tag data
% get timestamps, transmitter #s, and species IDs for pings on receiver of
% interest
recInd = find(table2array(telemDat(:,'receiver_sn'))==ShStRec);
pingDates = table2array(telemDat(recInd,'Date and Time (UTC)'));
trans = table2array(telemDat(recInd,'Transmitter'));
species = table2array(telemDat(recInd,'Species'));

% filter out non-sturgeon pings
ShStPings = find(strmatch(spec,species));
pingDates = pingDates(ShStPings);
trans = trans(ShStPings);

% Bin detections by hour, only count each individual fish once per bin
timeBins = startEff:minutes(60):endEff;
[N,timeBins,bin] = histcounts(pingDates,timeBins);
func = @(x) length(unique(x));
[G,ID] = findgroups(bin);
counts = splitapply(func,trans,G);

telemCounts.TimeBin = timeBins(1:(end-1))';
telemCounts.Count = zeros(length(timeBins)-1,1);
telemCounts.Count(ID(2:end)) = counts(2:end);

%% Wrangle automated pulse detector output

% Find timestamps of automated pulse detections within effort period
% pulseInds = find(strcmp(table2array(autDet(:,'SoundType')),'Pulse'));
% pulseTimes = sort(datetime(datevec(table2array(autDet(pulseInds,'BeginDateTime')))));
timeStrs = regexp(table2array(autDet(:,'Begin File')),'\d{8}_\d{6}','match');
fileTimes = (datetime([timeStrs{:}],'InputFormat','yyyyMMdd_HHmmss'))';
pulseTimes = sort(fileTimes+seconds(table2array(autDet(:,'Begin Time (s)'))));
pulseTimes = pulseTimes(pulseTimes>=startEff & pulseTimes<=endEff);
onEffInd = find(pulseTimes>=startEff & pulseTimes<=endEff);

% % Get rid of detections during boat passages using manual boat annotations
% boatInds = find(strcmp(table2array(manDet(:,'SoundType')),'Boat') | strcmp(table2array(manDet(:,'SoundType')),'Mechanical Noise'));
% % boatStarts = datevec(table2array(manDet(boatInds,'BeginDateTime')));
% % boatEnds = datevec(cellstr(char(table2array(manDet(boatInds,'EndClockTime')))));
% % boatEnds(:,1:3) = boatStarts(:,1:3);
% boatStarts = fileSetStart+seconds(table2array(manDet(boatInds,'Begin Time (s)')));
% boatEnds = fileSetStart+seconds(table2array(manDet(boatInds,'End Time (s)')));
% boatTimes = sortrows([datetime(boatStarts),datetime(boatEnds)]);
% for i = 1:length(boatTimes)
%     badDets = find(pulseTimes>=boatTimes(i,1) & pulseTimes<=boatTimes(i,2));
%     pulseTimes(badDets) = [];
% end

% Find detections which are labeled and not junk
goodInds = find(~isnan(autDet.Label) & autDet.Label~=99);

% Filter detections
keepInd = intersect(onEffInd,goodInds);

% Bin automated pulse detections by 1hr bin
[Npulse,timeBins,pulseBin] = histcounts(pulseTimes(keepInd),timeBins);
pulseCounts.TimeBin = timeBins(1:(end-1))';
pulseCounts.Count = Npulse';

%% Wrangle rumble annotations
% Get timestamps of automated pulse detections within effort period
rumbleInds = find(strcmp(table2array(manDet(:,'SoundType')),'Rumble'));
% rumbleTimes = sort(datetime(datevec(table2array(manDet(rumbleInds,'BeginDateTime')))));
rumbleTimes = sort(fileSetStart+seconds(table2array(manDet(rumbleInds,'Begin Time (s)'))));
rumbleTimes = rumbleTimes(rumbleTimes>=startEff & rumbleTimes<=endEff);

% Bin rumbles by 1hr bin
[Nrumble,timeBins,rumbleBin] = histcounts(rumbleTimes,timeBins);
rumbleCounts.TimeBin = timeBins(1:(end-1))';
rumbleCounts.Count = Nrumble';

%% Plot
% plot # individuals fish per hour and # of each call type per hour
figure(1)
subplot(3,1,1)
plot(telemCounts.TimeBin,telemCounts.Count,'LineWidth',2);
xlim([startEff,endEff]);
% ylim([0,max(telemCounts.Count)])
title(sprintf('Shortnose Sturgeon Counts on Receiver %d',ShStRec));
subplot(3,1,2)
plot(pulseCounts.TimeBin,pulseCounts.Count,'LineWidth',2)
% patch([boatTimes';flipud(boatTimes')],[repmat(0,2,length(boatTimes));repmat(max(pulseCounts.Count),2,length(boatTimes))],[224,224,224]./255,'EdgeColor','none');
xlim([startEff,endEff]);
% ylim([0,max(pulseCounts.Count)])
title(sprintf('Automated Pulse Detections'));
subplot(3,1,3)
plot(rumbleCounts.TimeBin,rumbleCounts.Count,'LineWidth',2)
xlim([startEff,endEff]);
% ylim([0,max(pulseCounts.Count)])
title(sprintf('Rumbles'));


%% Boxplots of rumbles vs telem pings vs net catches vs sonar counts

days = [datetime(2021,11,19,0,0,0),datetime(2021,12,2,0,0,0)];

% count rumbles during survey hours of each day
nov19rumb = sum(rumbleTimes>=datetime(2021,11,19,8,0,0) & rumbleTimes<=datetime(2021,11,19,16,0,0));
dec2rumb = sum(rumbleTimes>=datetime(2021,12,2,8,0,0) & rumbleTimes<=datetime(2021,12,2,16,0,0));

pingInd = find(pingDates>=datetime(2021,11,19,8,0,0) & pingDates<=datetime(2021,11,19,16,0,0));
nov19telem = length(unique(trans(pingInd)));
pingInd = find(pingDates>=datetime(2021,12,2,8,0,0) & pingDates<=datetime(2021,12,2,16,0,0));
dec2telem = length(unique(trans(pingInd)));

nov19net = 21+32+24;
dec2net = 0;

nov19sonar = mean([30+4+33+6;20+1+25+1]);
dec2sonar = mean([167+54+196+49;63+21+103+12]);

figure(2)
plot(days,[nov19rumb,dec2rumb],'-o')
hold on
plot(days,[nov19telem,dec2telem],'-o')
plot(days,[nov19net,dec2net],'-o')
plot(days,[nov19sonar,dec2sonar],'-o')
hold off
xlim([datetime(2021,11,18),datetime(2021,12,3)])
legend({'Rumbles','Telemetry Tags','Net Catches','Sonar Counts'},'Location','northwest')


