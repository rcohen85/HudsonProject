% Test for correlation between fish counts, soundscape metrics, and
% environmental conditions

%% Import and organize fish count data
fishData20 = readtable('C:\Users\rec297\CCB\HudsonProject\BlackCreekHerringCountsPerHour.xlsx','Sheet',1);
fishData21 = readtable('C:\Users\rec297\CCB\HudsonProject\BlackCreekHerringCountsPerHour.xlsx','Sheet',2);

% pull counts and weir set starts/ends
countsPerHour20 = table2array(fishData20(:,11));
weirSets20 = table2array(fishData20(:,'daySet'))+days(table2array(fishData20(:,'timeSet')));
weirPulls20 = table2array(fishData20(:,'dayPull'))+days(table2array(fishData20(:,'timePull')));

countsPerHour21 = table2array(fishData21(:,11));
weirSets21 = table2array(fishData21(:,'daySet'))+days(table2array(fishData21(:,'timeSet')));
weirPulls21 = table2array(fishData21(:,'dayPull'))+days(table2array(fishData21(:,'timePull')));

% quick plot to check that data is continuous, no gaps (start of each set is end of previous set)
figure(1)
plot(weirSets20,ones(length(weirSets20),1),'o')
hold on
plot(weirPulls20,ones(length(weirPulls20),1),'*')
hold off

figure(2)
plot(weirSets21,ones(length(weirSets21),1),'o')
hold on
plot(weirPulls21,ones(length(weirPulls21),1),'*')
hold off

disp('Examine plots for data continuity, then press any key to continue')
pause;

% calculate average fish/hr/day
weirDurs20 = hours(split(between(weirSets20,weirPulls20),'Time'));
weirDurs20(:,2) = floor(weirDurs20(:,1));
weirDurs20(:,3) = rem(weirDurs20(:,1),1);
allHours20 = (dateshift(weirSets20(1),'start','hour'):hours(1):dateshift(weirPulls20(end),'end','hour'))';
hourlyCounts20 = zeros(length(allHours20),1);

weirDurs21 = hours(split(between(weirSets21,weirPulls21),'Time'));
weirDurs21(:,2) = floor(weirDurs21(:,1));
weirDurs21(:,3) = rem(weirDurs21(:,1),1);
allHours21 = (dateshift(weirSets21(1),'start','hour'):hours(1):dateshift(weirPulls21(end),'end','hour'))';
hourlyCounts21 = zeros(length(allHours21),1);

for i = 1:length(weirSets20)

    % find which hour bins overlap with this set
    thisChunk = find(allHours20>=dateshift(weirSets20(i),'start','hour') & allHours20<dateshift(weirPulls20(i),'end','hour'));

    if length(thisChunk)>=2
        % scale counts in first hour and last hours proportional to effort
        firstEff = hours(dateshift(weirSets20(i),'end','hour') - weirSets20(i));
        scaledFirstHrCount = countsPerHour20(i)*firstEff;
        lastEff = hours(weirPulls20(i) - dateshift(weirPulls20(i),'start','hour'));
        scaledLastHrCount = countsPerHour20(i)*lastEff;
        hourlyCounts20([thisChunk(1),thisChunk(end)],1) = [hourlyCounts20(thisChunk(1)) + scaledFirstHrCount,...
                                                            hourlyCounts20(thisChunk(end)) + scaledLastHrCount];

        if length(thisChunk)>=3
            % plug unscaled countsPerHour into all complete hours
            hourlyCounts20(thisChunk(2):thisChunk(end-1),1) = countsPerHour20(i);
        end

    elseif length(thisChunk)==1
        hourEff = hours(weirPulls20(i) - weirSets20(i));
        scaledCounts = countsPerHour20(i)*hourEff;
        hourlyCounts20(thisChunk) = hourlyCounts20(thisChunk) + scaledCounts;
    end
end

for i = 1:length(weirSets21)

    % find which hour bins overlap with this set
    thisChunk = find(allHours21>=dateshift(weirSets21(i),'start','hour') & allHours21<dateshift(weirPulls21(i),'end','hour'));

    if length(thisChunk)>=2
        % scale counts in first hour and last hours proportional to effort
        firstEff = hours(dateshift(weirSets21(i),'end','hour') - weirSets21(i));
        scaledFirstHrCount = countsPerHour21(i)*firstEff;
        lastEff = hours(weirPulls21(i) - dateshift(weirPulls21(i),'start','hour'));
        scaledLastHrCount = countsPerHour21(i)*lastEff;
        hourlyCounts21([thisChunk(1),thisChunk(end)],1) = [hourlyCounts21(thisChunk(1)) + scaledFirstHrCount,...
                                                            hourlyCounts21(thisChunk(end)) + scaledLastHrCount];

        if length(thisChunk)>=3
            % plug unscaled countsPerHour into all complete hours
            hourlyCounts21(thisChunk(2):thisChunk(end-1),1) = countsPerHour21(i);
        end

    elseif length(thisChunk)==1
        hourEff = hours(weirPulls21(i) - weirSets21(i));
        scaledCounts = countsPerHour21(i)*hourEff;
        hourlyCounts21(thisChunk) = hourlyCounts21(thisChunk) + scaledCounts;
    end
end

% Tally fish counts per day (only for days with 100% effort)
dayBins20 = dateshift(weirSets20(1),'end','day'):days(1):dateshift(weirPulls20(end),'start','day');
dayBins21 = dateshift(weirSets21(1),'end','day'):days(1):dateshift(weirPulls21(end),'start','day');
[N, dayBins,fishBin20] = histcounts(allHours20,dayBins20);
[N, dayBins,fishBin21] = histcounts(allHours21,dayBins21);
dailyFishCounts20 = grpstats(table(hourlyCounts20,fishBin20),'fishBin20','sum');
dailyFishCounts21 = grpstats(table(hourlyCounts21,fishBin21),'fishBin21','sum');

fishDates = [dayBins20(1:(end-1))';dayBins21(1:(end-1))'];
dailyFishCounts = [table2array(dailyFishCounts20(2:end,3));table2array(dailyFishCounts21(2:end,3))]; % hours in bin zero were from partial effort days, disregard

%% Import soundscape metrics

% Broadband (BB) noise level data
BBdata = [readtable('W:\projects\2022_NOAA-NERRS_HudsonNY_144488\SoundscapeMetrics\BlackCreek\BC02\BC02_BB_1h.csv');
    readtable('W:\projects\2022_NOAA-NERRS_HudsonNY_144488\SoundscapeMetrics\BlackCreek\BC04\BC04_BB_1h.csv');
    readtable('W:\projects\2022_NOAA-NERRS_HudsonNY_144488\SoundscapeMetrics\BlackCreek\BC05\BC05_BB_1h.csv');
    readtable('W:\projects\2022_NOAA-NERRS_HudsonNY_144488\SoundscapeMetrics\BlackCreek\BC06\BC06_BB_1h.csv')];
BBdates = datetime(table2array(BBdata(:,1)),'InputFormat','yyyy-MM-dd''T''HH'); % dates get read in messed up from this one
BBlevels = str2double(extractAfter(table2array(BBdata(:,3)),'00.000Z,'));

% Import octave (O) levels
OLdata = [readtable('W:\projects\2022_NOAA-NERRS_HudsonNY_144488\SoundscapeMetrics\BlackCreek\BC02\BC02_OL_1h.csv');
    readtable('W:\projects\2022_NOAA-NERRS_HudsonNY_144488\SoundscapeMetrics\BlackCreek\BC04\BC04_OL_1h.csv');
    readtable('W:\projects\2022_NOAA-NERRS_HudsonNY_144488\SoundscapeMetrics\BlackCreek\BC05\BC05_OL_1h.csv');
    readtable('W:\projects\2022_NOAA-NERRS_HudsonNY_144488\SoundscapeMetrics\BlackCreek\BC06\BC06_OL_1h.csv')];
OLdates = datetime(table2array(OLdata(:,1)),'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z');
OLevels = table2array(OLdata(:,2:end));
OCenters = OLdata.Properties.VariableNames(2:end);
% OCenters = str2double(strrep(cellfun(@(x) extractAfter(x,3),OCenters,'UniformOutput',0),'_','.'));
OCenters = strrep(OCenters,'L_','');

% Import third octave (TO) levels
TOLdata = [readtable('W:\projects\2022_NOAA-NERRS_HudsonNY_144488\SoundscapeMetrics\BlackCreek\BC02\BC02_TOL_1h.csv');
    readtable('W:\projects\2022_NOAA-NERRS_HudsonNY_144488\SoundscapeMetrics\BlackCreek\BC04\BC04_TOL_1h.csv');
    readtable('W:\projects\2022_NOAA-NERRS_HudsonNY_144488\SoundscapeMetrics\BlackCreek\BC05\BC05_TOL_1h.csv');
    readtable('W:\projects\2022_NOAA-NERRS_HudsonNY_144488\SoundscapeMetrics\BlackCreek\BC06\BC06_TOL_1h.csv')];
TOLdates = datetime(table2array(TOLdata(:,1)),'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z');
TOLevels = table2array(TOLdata(:,2:end));
TOCenters = TOLdata.Properties.VariableNames(2:end);
% TOCenters = str2double(strrep(cellfun(@(x) extractAfter(x,4),TOCenters,'UniformOutput',0),'_','.'));
TOCenters = strrep(TOCenters,'OL_','');

% Compute daily average BB, O, and TO levels
dayBins = dateshift(BBdates(1),'start','day'):1:dateshift(BBdates(end),'end','day');
[N, dayBins,BBbin] = histcounts(BBdates,dayBins);
dailyBBlevels = grpstats(table(BBlevels,BBbin),'BBbin','mean');

[N, dayBins,OLbin] = histcounts(OLdates,dayBins);
dailyOlevels = grpstats(table(OLevels,OLbin),'OLbin','mean');

[N, dayBins,TOLbin] = histcounts(TOLdates,dayBins);
dailyTOlevels = grpstats(table(TOLevels,TOLbin),'TOLbin','mean');

% only retain days with 100% effort
fullDays = find(dailyBBlevels.GroupCount==24);
noiseDates = dayBins(dailyBBlevels.BBbin(fullDays));
dailyNoiseData = [table2array(dailyBBlevels(fullDays,3)),table2array(dailyOlevels(fullDays,3)),table2array(dailyTOlevels(fullDays,3))];
noiseBands = ['BB',OCenters,TOCenters];


%% Calculate correlation w soundscape metrics
% Align dates and create matrix of counts, soundscape metrics, and environmental covariates
sharedDates = intersect(fishDates,noiseDates);
fishInd = ismember(fishDates,sharedDates);
noiseInd = ismember(noiseDates,sharedDates);

fishNoiseMat = [dailyFishCounts(fishInd),dailyNoiseData(noiseInd,:)];
fishNoiseTable = array2table(fishNoiseMat,'VariableNames',['FishCounts',noiseBands]);
% dataTable = array2table(dataMat);

% Compute correlation matrix and plot (2 approaches)
figure(1), clf % heatmap w correlation coefficient
hotCold = interp1([-1,0,1]',[204 0 0; 255 255 255; 0 0 204]./255,linspace(-1,1,51));
c = corr(fishNoiseMat);
isupper = logical(triu(ones(size(c)),1));
c(isupper) = NaN;
h = heatmap(c,'MissingDataColor','w');
colormap(hotCold)
clim([-1 1])
labels = ['FishCounts',noiseBands];
h.XDisplayLabels = labels;
h.YDisplayLabels = labels;
saveas(gcf,'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\BC_Herring\FishCounts_Noise_Corr1.fig')
saveas(gcf,'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\BC_Herring\FishCounts_Noise_Corr1.png')

figure(2), clf % scatter plots w regression line & correlation coefficient
[Cr,~,h] = corrplot_RC(fishNoiseTable);
islower = logical(tril(ones(size(Cr)),-1));
mirrors = find(islower == 1);
for i=1:size(mirrors)
    delete(subplot(31,31,mirrors(i)));    
end
saveas(gcf,'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\BC_Herring\FishCounts_Noise_Corr2.fig')
saveas(gcf,'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\BC_Herring\FishCounts_Noise_Corr2.png')

%% Import in environmental data
% meteorological data
metData = readtable('C:\Users\rec297\CCB\HudsonProject\NERRS_MeteorologicalData\HUDNPMET.csv');
metDates = datetime(table2array(metData(:,3)),'InputFormat','MM/dd/yyyy HH:mm');
metDatesInd = find(~isnat(metDates));
% hydrology & water chemistry data
waterData = readtable('C:\Users\rec297\CCB\HudsonProject\NERRS_WaterQualityData\HUDNPWQ.csv');
waterDates = datetime(table2array(waterData(:,3)),'InputFormat','MM/dd/yyyy HH:mm');
waterDatesInd = find(~isnat(waterDates));

% Compute daily averages of air temp (C), mean wind speed (m/s), wind direction 
% (degrees), wind dir std dev, and precip (mm)
metDayBins = dateshift(metDates(1),'start','day'):1:dateshift(metDates(metDatesInd(end)),'end','day');
[N, dayBins,metBin] = histcounts(metDates(metDatesInd),metDayBins);
dailyMetVals = grpstats([metData(metDatesInd,[8,14,19,21,25]),table(metBin)],'metBin','mean');
dailyMetVals.Properties.VariableNames(3:7) = {'AirTmp','WinSpd','WinDir','WDSD','Precip'};
metDays = metDayBins(1:(end-1));

% Compute daily averages of water temp (C), DO (mg/L), depth (m), pH, turbidity
watDayBins = dateshift(waterDates(1),'start','day'):1:dateshift(waterDates(waterDatesInd(end)),'end','day');
[N, dayBins,waterBin] = histcounts(waterDates(waterDatesInd),watDayBins);
dailyWaterVals = grpstats([waterData(waterDatesInd,[7,15,17,21,23]),table(waterBin)],'waterBin','mean');
dailyWaterVals.Properties.VariableNames(3:7) = {'WatTmp','DO','Depth','pH','Turb'};
watDays = watDayBins(1:(end-1));

%% Calculate correlation btwn fish counts & environmental data
% align dates
envDates = intersect(metDays,watDays);
dailyMetVals = dailyMetVals(ismember(metDays,envDates),:);
metDays = metDays(ismember(metDays,envDates));
dailyWaterVals = dailyWaterVals(ismember(watDays,envDates),:);
watDays = watDays(ismember(watDays,envDates));

sharedDates = intersect(fishDates,envDates);
fishInd = ismember(fishDates,sharedDates);
metInd = ismember(metDays,sharedDates);
watInd = ismember(watDays,sharedDates);

fishEnvMat = [dailyFishCounts(fishInd),table2array(dailyWaterVals(watInd,3:7)),table2array(dailyMetVals(metInd,3:7))];
fishEnvTable = array2table(fishEnvMat,'VariableNames',['FishCounts',dailyWaterVals.Properties.VariableNames(3:7),dailyMetVals.Properties.VariableNames(3:7)]);

% Compute correlation matrix and plot (2 approaches)
figure(3), clf % heatmap w correlation coefficient
hotCold = interp1([-1,0,1]',[204 0 0; 255 255 255; 0 0 204]./255,linspace(-1,1,51));
c = corr(fishEnvMat);
isupper = logical(triu(ones(size(c)),1));
c(isupper) = NaN;
h = heatmap(c,'MissingDataColor','w');
colormap(hotCold)
clim([-1 1])
labels = ['FishCounts',dailyWaterVals.Properties.VariableNames(3:7),dailyMetVals.Properties.VariableNames(3:7)];
h.XDisplayLabels = labels;
h.YDisplayLabels = labels;
saveas(gcf,'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\BC_Herring\FishCounts_Env_Corr1.fig')
saveas(gcf,'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\BC_Herring\FishCounts_Env_Corr1.png')

figure(4), clf % scatter plots w regression line & correlation coefficient
[Cr,~,h] = corrplot_RC(fishEnvTable);
islower = logical(tril(ones(size(Cr)),-1));
mirrors = find(islower == 1);
for i=1:size(mirrors)
    delete(subplot(11,11,mirrors(i)));    
end
saveas(gcf,'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\BC_Herring\FishCounts_Env_Corr2.fig')
saveas(gcf,'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\BC_Herring\FishCounts_Env_Corr2.png')

%% Calculate correlation btwn soundscape metrics & environmental data
sharedDates = intersect(noiseDates,watDays);
noiseInd = ismember(noiseDates,sharedDates);
watInd = ismember(watDays,sharedDates); % can use same indices for meteorological data, since the date vectors were aligned above

envNoiseMat = [dailyNoiseData(noiseInd,:),table2array(dailyWaterVals(watInd,3:7)),table2array(dailyMetVals(watInd,3:7))];
envNoiseTable = array2table(envNoiseMat,'VariableNames',[noiseBands,dailyWaterVals.Properties.VariableNames(3:7),dailyMetVals.Properties.VariableNames(3:7)]);

% Compute correlation matrix and plot (2 approaches)
figure(5), clf % heatmap w correlation coefficient
hotCold = interp1([-1,0,1]',[204 0 0; 255 255 255; 0 0 204]./255,linspace(-1,1,51));
c = corr(envNoiseMat);
isupper = logical(triu(ones(size(c)),1));
c(isupper) = NaN;
h = heatmap(c,'MissingDataColor','w','CellLabelFormat','%.2f');
colormap(hotCold)
clim([-1 1])
labels = [noiseBands,dailyWaterVals.Properties.VariableNames(3:7),dailyMetVals.Properties.VariableNames(3:7)];
h.XDisplayLabels = labels;
h.YDisplayLabels = labels;
saveas(gcf,'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\BC_Herring\Noise_Env_Corr1.fig')
saveas(gcf,'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\BC_Herring\Noise_Env_Corr1.png')

figure(6), clf % scatter plots w regression line & correlation coefficient
[Cr,~,h] = corrplot_RC(envNoiseTable);
islower = logical(tril(ones(size(Cr)),-1));
mirrors = find(islower == 1);
for i=1:size(mirrors)
    delete(subplot(40,40,mirrors(i)));    
end
saveas(gcf,'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\BC_Herring\Noise_Env_Corr2.fig')
saveas(gcf,'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\BC_Herring\Noise_Env_Corr2.png')

%% Import acoustic indices



%% Calcuate correlation btwn fish counts & acoustic indices

%% Calculate correltion btwn acoustic indices & environmental data
