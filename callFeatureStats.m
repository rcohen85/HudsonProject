
% selDir = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AtlStur';
% fileWildcard = '*cleanDetsForStats.txt';
% labCol = 'SoundType';
selDir = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AtlStur\ManualSelectionTables';
fileWildcard = 'BB24.ManualSelections.txt';
labCol = 'Tags';
labs = {'Thunder','Noise'};
saveDir = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AtlStur\Figures';

Fs = 24000;
cal = 176.3; % end-to-end calibration value (RL resulting in full scale signal)

fileList = dir(fullfile(selDir,fileWildcard));
cal_muPa = power (10, cal / 20);
Waves = {};
Waves_muPa = {};
Specs = [];
PSD = [];
Labels = {};
startTimes = [];
endTimes = [];
lowFreq = [];
highFreq = [];
% dur90 = [];
dur = [];
for i=1:length(fileList)

    tab = readtable(fullfile(selDir,fileList(i).name),'Delimiter',"\t",'VariableNamingRule',"preserve");
    tab = tab(strmatch('Spectrogram 1',tab.View),:);
% if ~contains(tab.Properties.VariableNames,'Delta Time (s)')
%     tab(:,'Delta Time (s)') = array2table(table2array(tab(:,'End Time (s)')) - table2array(tab(:,'Begin Time (s)')));
% end

    soundFiles = table2array(unique(tab(:,'Begin Path')));

    for j=1:length(soundFiles)
        % get info on this audio file
        info = audioinfo(soundFiles{j});
        % start datetime of file
        fileStart = datetime(regexp(soundFiles{j},'\d{8}_\d{6}','match'),'InputFormat','uuuuMMdd_HHmmss');

        % find selections in this sound file
        thisFileInd = strmatch(soundFiles{j},table2array(tab(:,'Begin Path')));
        Labels = [Labels;table2array(tab(thisFileInd,labCol))];
        startTimes = [startTimes; seconds(table2array(tab(thisFileInd,'File Offset (s)')))+fileStart];
        endTimes = [endTimes; seconds(table2array(tab(thisFileInd,'File Offset (s)'))+table2array(tab(thisFileInd,'Delta Time (s)')))+fileStart];
        lowFreq = [lowFreq;table2array(tab(thisFileInd,'Low Freq (Hz)'))];
        highFreq = [highFreq;table2array(tab(thisFileInd,'High Freq (Hz)'))];
%         dur90 = [dur90;table2array(tab(thisFileInd,'Dur 90% (s)'))];
        dur = [dur;table2array(tab(thisFileInd,'Delta Time (s)'))];
        
        startSamps = table2array(tab(thisFileInd,'File Offset (s)')).*info.SampleRate;
        endSamps = startSamps + table2array(tab(thisFileInd,'Delta Time (s)')).*info.SampleRate;

        for k = 1:length(thisFileInd)
            % get waveform clip
            clip = double(audioread(soundFiles{j},[round(startSamps(k))-10000,round(endSamps(k))]));
            Waves = [Waves;clip'];
            clip_muPa = clip.*cal_muPa;
            Waves_muPa = [Waves_muPa;clip_muPa'];

            % Calculate power spectrum
            nfft = Fs;
            noverlap = floor(nfft*0.75);
            if length(clip_muPa)>=nfft
                win = hann(nfft);
            else
                win = hann(length(clip_muPa));
            end
            [s,f,t,ps] = spectrogram(clip_muPa,win,noverlap,nfft,Fs,'psd');
            s = s/(sum(win)/nfft);    % accound for reduction of energy due to windowing
            s = s/nfft; % normalize by nfft to preserve physical amplitude
            P = (abs(s).^2)*2; % convert to power and double to conserve energy from double-sided spectrum
            Specs = [Specs;mean(P,2)'];
            %PSD = [PSD;mean(ps,2)'];

        end

    end

end
 Specs_dB = 10*log10(Specs); % 10log10 for power spectrum (would be 20log10 for magnitude spectrum)

% Calculate spectra
% window_length = max(cell2mat(cellfun(@(x) length(x(10000:end)),Waves,'UniformOutput',0)));
% nfft = 2^nextpow2(window_length);
% cal_muPa = power (10, cal/20);
% Waves_muPa = cellfun(@(x) x(10000:end).*cal_muPa,Waves,'UniformOutput',0);
% % Waves_muPa = Waves.*cal_muPa;
% Specs = (cell2mat(cellfun(@(x) abs(fft(x,nfft)),Waves_muPa,'UniformOutput',0))./floor(nfft))*2;
% % Specs = (abs(fft(Waves_muPa,nfft,2))./floor(nfft))*2;
% Specs = Specs(:,1:(floor(nfft/2))+1);
% PSD = Specs.^2;
% % Specs_dB = 20.*log10(abs(Specs));
% % Specs = cell2mat(cellfun(@(x) pwelch(x,8,0,nfft),Waves,'UniformOutput',0));
% f = info.SampleRate/nfft*(0:(nfft/2));

% Pad waveforms to a consistent length
window_length = max(cell2mat(cellfun(@(x) length(x),Waves_muPa,'UniformOutput',0)));
for i=1:length(Waves_muPa)
    zeroPad = window_length - length(Waves_muPa{i,1});
    Waves_muPa{i,1} = [Waves_muPa{i,1},ones(1,zeroPad)*mean(Waves_muPa{i,1})];
end
Waves_muPa = cat(1,Waves_muPa{:});

% Calculate features for each label
allData = struct('Label',[],'Wave',[],'PowSpec_lin',[],'FreqVec',[],'PeakFreqs',[],'Durations',[],'Bandwidths',[]);
for i=1:length(labs)
    thisLab = strmatch(labs{i},Labels);

    pkFreqs = [];
    BWs = [];
    for j=1:length(thisLab)

        smoothSpec = fastsmooth(Specs_dB(thisLab(j),:),10);
        [m,whichMinF] = min(abs(lowFreq(thisLab(j))-f)); % only look at spectrum within freqs described by bounding box
        [m,whichMaxF] = min(abs(highFreq(thisLab(j))-f));
        [valMx,posMx] = max(smoothSpec(1,whichMinF:whichMaxF)); % find amplitude and value of peak freq
%         [valMx,posMx] = max(Specs_dB(thisLab(j),whichMinF:whichMaxF)); % find amplitude and value of peak freq
        posMx = posMx+whichMinF-1; % correct index of maximum to refer to entire freq vector, not subset of freqs in bounding box
        pkFreqs = [pkFreqs;f(posMx)];

        % Calculate -10dB bandwidth
        low = valMx-10; % -10dB value relative to peak
        %from peak freq, walk along spectrogram until low is reached on either side
        slopeup=fliplr(smoothSpec(1:posMx));
        for e10dB=1:length(slopeup)
            if slopeup(e10dB)<low %stop at value < -10dB: point of lowest frequency
                break
            end
        end
        slopedown=smoothSpec(posMx:length(smoothSpec));
        for o10dB=1:length(slopedown)
            if slopedown(o10dB)<low %stop at value < -10dB: point of highest frequency
                break
            end
        end
        
        high10dB = f(posMx+o10dB-1); %-10dB highest frequency in Hz
        low10dB = f(posMx-e10dB+1); %-10dB lowest frequency in Hz
        BWs = [BWs,high10dB-low10dB];
        if high10dB-low10dB==0
            fprintf('Warning: 0Hz Bandwidth computed\n')
            pause
        end

    end
    
%     durs = endTimes(thisLab)-startTimes(thisLab);
%     durs = dur90(thisLab);
    durs = dur(thisLab);

    allData(i).Label = labs{i};
    allData(i).Wave = Waves_muPa(thisLab,:);
    allData(i).PowSpec_lin = Specs(thisLab,:);
%     allData(i).PSD = PSD(thisLab,:);
    allData(i).FreqVec = f;
    allData(i).PeakFreqs = pkFreqs;
    allData(i).Durations = durs;
    allData(i).Bandwidths = BWs;
    

end

% Save call stats
featureDir = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AtlStur\Thunders';
save(fullfile(featureDir,'FieldThunders.m'),'allData','-v7.3');

%% Create Figures
% Plot mean (filtered) waveform
% [b,a] = butter(5,0.0083,'low'); % lowpass filter
% filtWaves = filtfilt(b,a,allData(1).Wave');
% [b,a] = butter(5,0.00041667,'high'); % highpass filter
% filtWaves = (filtfilt(b,a,filtWaves))';
% meanWav = mean(filtWaves,1);
% meanWav = meanWav./1e6;
% % meanWav_norm = meanWav-min(meanWav);
% % meanWav_norm = meanWav_norm./max(meanWav_norm);
% figure(1),clf
% plot((1:length(meanWav))/Fs,meanWav,'LineWidth',1.5)
% xlim([0,3.5])
% % ylim([min(meanWav),max(meanWav)]);
% ylabel('Pressure (Pa)')
% xlabel('Time (s)')
% set(gca,'FontSize',14)
% saveas(gcf,fullfile(saveDir,strcat(strrep(labs{1,1},' ',''),'_MeanWave_Pa.png')));
% exportgraphics(gcf,fullfile(saveDir,strcat(strrep(labs{1,1},' ',''),'_MeanWave_Pa.pdf')),'ContentType','vector');


% Plot mean/median spectra
% stDev = std(allData(1).Spec_lin(:,3:end)); % Looks wonky!
% stDev = std(20.*log10(allData(1).Spec_lin(:,3:end)));
figure(4),clf
% plot(allData(1).FreqVec(3:end),20*log10(mean(allData(1).Spec_lin(:,3:end),1)),'LineWidth',1.5)
plot(allData(1).FreqVec(3:end),10*log10(quantile(allData(1).PowSpec_lin(:,3:end),0.5)),'LineWidth',1.5)
hold on
% plot(allData(1).FreqVec(3:end),10*log10(mean(allData(1).PowSpec_lin(:,3:end),1))+stDev,'--','Color','#9E9E9E','LineWidth',1.5)
% plot(allData(1).FreqVec(3:end),10*log10(mean(allData(1).PowSpec_lin(:,3:end),1))-stDev,'--','Color','#9E9E9E','LineWidth',1.5)
plot(allData(1).FreqVec(3:end),10*log10(quantile(allData(1).PowSpec_lin(:,3:end),0.25)),'--','Color','#9E9E9E','LineWidth',1.5)
plot(allData(1).FreqVec(3:end),10*log10(quantile(allData(1).PowSpec_lin(:,3:end),0.75)),'--','Color','#9E9E9E','LineWidth',1.5)
plot(allData(1).FreqVec(3:end),10*log10(quantile(allData(2).PowSpec_lin(:,3:end),0.5)),'r','LineWidth',1.5)
hold off
legend([labs{1,1},', N = ' num2str(length(allData(1).PeakFreqs))],'','',['Background Noise, N = ' num2str(length(allData(2).PeakFreqs))])
xlim([3 250])
ylim([65 120])
ylabel('RL (dB re \muPa^2)')
xlabel('Frequency (Hz)')
set(gca,'FontSize',14)
% saveas(gcf,fullfile(saveDir,strcat(strrep(labs{1,1},' ',''),'_MeanSpec.png')))
saveas(gcf,fullfile(saveDir,strcat('Captive',strrep(labs{1,1},' ',''),'_MedianSpec.png')))
% exportgraphics(gcf,fullfile(saveDir,strcat(strrep(labs{1,1},' ',''),'_MeanSpec.pdf')),'ContentType','vector');
exportgraphics(gcf,fullfile(saveDir,strcat('Captive',strrep(labs{1,1},' ',''),'_MedianSpec.pdf')),'ContentType','vector');

% Plot mean spectrogram
% nfft = 2^nextpow2(3000);
% win = hann(nfft);
% overlap = floor(nfft*0.75);
% [s,f,t] = spectrogram(meanWav,win,overlap,nfft,Fs);
% s = s/(sum(win)/nfft);    % accound for reduction of energy due to windowing
% s = s/nfft; % normalize by nfft to preserve physical amplitude
% P = (abs(s).^2)*2; % convert to power and double to conserve energy from double-sided spectrum
% s = 10*log10(P);
% figure(3)
% colormap jet
% imagesc(t,f,s)
% set(gca,'YDir','normal')
% xlabel('Time (s)')
% ylabel('Frequency (Hz)')
% ylim([0 350])
% xlim([0 4])
% clim([-130,45])
% set(gca,'FontSize',14)
% saveas(gcf,fullfile(saveDir,strcat(strrep(labs{1,1},' ',''),'_MeanSpectrogram.png')))

% % Plot PSD
% P = prctile(allData(1).PSD,[0.01,0.05,0.25,0.5,0.75,0.95,0.99]);
% figure(4),clf
% plot(allData(1).FreqVec,10*log10(P),'LineWidth',1.5)
% legend({'1%','5%','25%','50%','75%','95%','99%'})
% xlim([0 450])
% ylim([45,80])
% ylabel('Power Spectral Density')
% xlabel('Frequency (Hz)')
% set(gca,'FontSize',14)
% saveas(gcf,fullfile(saveDir,strcat(strrep(labs{1,1},' ',''),'_PSD.png')))

% Plot central tendency & variability of other features
figure(5),clf
subplot(1,3,1)
boxchart(allData(1).PeakFreqs)
title('Peak Frequency');
ylabel('Hz')
xticks([])
set(gca,'FontSize',14)
subplot(1,3,2)
boxchart(allData(1).Bandwidths)
title('-10 dB Bandwidth');
ylabel('Hz')
xticks([])
set(gca,'FontSize',14)
subplot(1,3,3)
boxchart(allData(1).Durations)
title('Duration');
ylabel('Seconds')
xticks([])
set(gca,'FontSize',14)
saveas(gcf,fullfile(saveDir,strcat(strrep(labs{1,1},' ',''),'_FeatureVariability.png')))


%% Plot representative spectrogram & waveform

% inFile = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_FLAC\AtlSturSpw\ST300\AS03\144488Hudson_024K_AS03_ST300-5683_20210610_234010Zm0500.flac';
% startS = 3209;
% dur = 40;
% inFile = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_FLAC\AtlSturSpw\ST300\AS03\144488Hudson_024K_AS03_ST300-5683_20210610_114025Zm0500.flac';
% startS = 42109;
% dur = 25;
% inFile = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_FLAC\AtlSturSpw\ST300\AS03\144488Hudson_024K_AS03_ST300-5683_20210617_233832Zm0500.flac';
% startS = 9866.75;
% dur = 7;
% Thunder w boat
% inFile = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_FLAC\AtlSturSpw\ST300\AS03\144488Hudson_024K_AS03_ST300-5683_20210610_114025Zm0500.flac';
% startS = 4110;
% dur = 90;
% Thunder w train
% inFile = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_FLAC\AtlSturSpw\ST300\AS03\144488Hudson_024K_AS03_ST300-5683_20210609_234043Zm0500.flac';
% startS = 40600;
% dur = 60;
% Thunder w drum chorus
% inFile = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_FLAC\AtlSturSpw\ST300\AS03\144488Hudson_024K_AS03_ST300-5683_20210610_114025Zm0500.flac';
% startS = 16133;
% dur = 45;
% Captive thunder
% inFile = "W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_FLAC\BearsBluff\BB24\144488Hudson_024K_BB24_ST300-6530_20240819_062958Zm0400.flac";
% startS = 498;
% dur = 35;
% inFile = "W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_FLAC\BearsBluff\BB24\144488Hudson_024K_BB24_ST300-6530_20240819_165958Zm0400.flac";
% startS = 370;
% dur = 30;
% Wild thunder representative example detail
% inFile = "W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_FLAC\AtlSturSpw\ST300\AS03\144488Hudson_024K_AS03_ST300-5683_20210609_234043Zm0500.flac";
% startS = 40175;
% dur = 3;
% Captive thunder representative example detail
% inFile = "W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_FLAC\BearsBluff\BB24\144488Hudson_024K_BB24_ST300-6530_20240819_115958Zm0400.flac";
% startS = 770;
% dur = 3;
% inFile = "W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_FLAC\BearsBluff\BB24\144488Hudson_024K_BB24_ST300-6530_20240819_092958Zm0400.flac";
% startS = 662;
% dur = 3;
inFile = "W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_FLAC\BearsBluff\BB24\144488Hudson_024K_BB24_ST300-6530_20240819_122958Zm0400.flac";
startS = 54.1;
dur = 3;
% inFile = "W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_FLAC\BearsBluff\BB24\144488Hudson_024K_BB24_ST300-6530_20240819_142958Zm0400.flac";
% startS = 130;
% dur = 3;
% inFile = "W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_FLAC\BearsBluff\BB24\144488Hudson_024K_BB24_ST300-6530_20240819_142958Zm0400.flac";
% startS = 381;
% dur = 3;
% inFile = "W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_FLAC\BearsBluff\BB24\144488Hudson_024K_BB24_ST300-6530_20240819_142958Zm0400.flac";
% startS = 594.8;
% dur = 3;

cal = 176.4;
cal_muPa = power (10, cal / 20);

info = audioinfo(inFile);
[y,Fs] = audioread(inFile,floor([startS*info.SampleRate,(startS*info.SampleRate)+(dur*info.SampleRate)]));
y = y.*cal_muPa;

% bandpass filter data (passband: [10,200]Hz)
[b,a] = butter(2,[8.3333e-04,0.0167]);
filtWave = filtfilt(b,a,y');
filtWave = filtWave./1e6;

figure(1),clf
plot((1:length(filtWave))/Fs,filtWave,'LineWidth',1.5)
xlim([0,3])
ylim([min(filtWave)+(min(filtWave)*0.1),max(filtWave)+(max(filtWave)*0.1)])
ylabel('Pressure (Pa)')
xlabel('Time (s)')
set(gca,'FontSize',14)
saveas(gcf,'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AtlStur\Figures\Revision\CaptiveThunderRepEx_Wave.png');
exportgraphics(gcf,'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AtlStur\Figures\Revision\CaptiveThunderRepEx_Wave.pdf','ContentType','vector');


nfft = 2^nextpow2(4000);
noverlap = floor(nfft*0.95);
win = hann(nfft);
[s,f,t] = spectrogram(y,win,noverlap,nfft,Fs);
s = s/(sum(win)/nfft);    % accound for reduction of energy due to windowing
s = s/nfft; % normalize by nfft to preserve physical amplitude
P = (abs(s).^2)*2; % convert to power and double to conserve energy from double-sided spectrum
s = 10*log10(P);

% Spectrogram
figure(2),clf
colormap jet
imagesc(t,f,s)
set(gca,'YDir','normal');
xlabel('Time (s)')
ylabel('Frequency (Hz)')
ylim([0 500])
yticks([0,500,1000]);
yticklabels({'0','500','1000'});
c = colorbar;
ylabel(c,'dB re: \muPa^2');
clim([50,125])
set(gca,'FontSize',14)
saveas(gcf,'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AtlStur\Figures\Revision\CaptiveThunderRepEx_nfft4096_overlap95.png');
exportgraphics(gcf,'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AtlStur\Figures\Revision\CaptiveThunderRepEx_nfft4096_overlap95.pdf','ContentType','vector');
