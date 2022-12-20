% Test for correlation between fish counts and soundscape metrics

% Import fish count data
fishData20 = readtable('C:\Users\rec297\CCB\HudsonProject\Black Creek Herring Counts 2020.xlsx','Sheet','Daily_Counts.hr');
fishDates = [datetime(table2array(fishData20(:,1)),'InputFormat','d-MMM')];
fishCounts = table2array(fishData20(:,17));

% Import broadband noise level data
BBdata = readtable('W:\projects\2022_NOAA-NERRS(UMich)_HudsonNY_144488\SoundscapeMetrics\BlackCreek\BC05_BB_1h.csv');
BBdates = datetime(table2array(BBdata(:,1)),'InputFormat','yyyy-MM-dd''T''HH'); % dates get read in messed up from this one
BBlevels = str2double(extractAfter(table2array(BBdata(:,3)),'00.000Z,'));

% Import octave levels
OLdata = readtable('W:\projects\2022_NOAA-NERRS(UMich)_HudsonNY_144488\SoundscapeMetrics\BlackCreek\BC05_OL_1h.csv');
OLdates = datetime(table2array(OLdata(:,1)),'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z');
OLevels = table2array(OLdata(:,2:end));
OCenters = OLdata.Properties.VariableNames(2:end);
OCenters = str2double(strrep(cellfun(@(x) extractAfter(x,3),OCenters,'UniformOutput',0),'_','.'));

% Import third octave levels
TOLdata = readtable('W:\projects\2022_NOAA-NERRS(UMich)_HudsonNY_144488\SoundscapeMetrics\BlackCreek\BC05_TOL_1h.csv');
TOLdates = datetime(table2array(TOLdata(:,1)),'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z');
TOLevels = table2array(TOLdata(:,2:end));
TOCenters = TOLdata.Properties.VariableNames(2:end);
TOCenters = str2double(strrep(cellfun(@(x) extractAfter(x,4),TOCenters,'UniformOutput',0),'_','.'));

% Compute daily average BB, O, and TO levels
dayBins = dateshift(BBdates(1),'start','day'):1:dateshift(BBdates(end),'end','day');
[N, edges,BBbin] = histcounts(BBdates,dayBins);
dailyBBlevels = grpstats(table(BBlevels,BBbin),'BBbin','mean');

[N, edges,OLbin] = histcounts(OLdates,dayBins);
dailyOlevels = grpstats(table(OLevels,OLbin),'OLbin','mean');

[N, edges,TOLbin] = histcounts(TOLdates,dayBins);
dailyTOlevels = grpstats(table(TOLevels,TOLbin),'TOLbin','mean');



