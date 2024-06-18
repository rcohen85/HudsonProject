
selDir = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AtlStur';
saveDir = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AtlStur\Figures';
fileWildcard = '*cleanDetsForStats.txt';
labs = {'Thunder','Noise'};
Fs = 24000;
cal = 176.3; % end-to-end calibration value (RL resulting in full scale signal)

fileList = dir(fullfile(selDir,fileWildcard));

Waves = {};
Waves_muPa = {};
Specs = [];
PSD = [];
Labels = {};
startTimes = [];
endTimes = [];
lowFreq = [];
highFreq = [];
for i=1:length(fileList)

    tab = readtable(fullfile(selDir,fileList(i).name),'Delimiter',"\t",'VariableNamingRule',"preserve");
    tab = tab(strmatch('Spectrogram 1',tab.View),:);

    soundFiles = table2array(unique(tab(:,'Begin Path')));

    for j=1:length(soundFiles)
        % get info on this audio file
        info = audioinfo(soundFiles{j});
        % start datetime of file
        fileStart = datetime(regexp(soundFiles{j},'\d{8}_\d{6}','match'),'InputFormat','uuuuMMdd_HHmmss');

        % find selections in this sound file
        thisFileInd = strmatch(soundFiles{j},table2array(tab(:,'Begin Path')));
        Labels = [Labels;table2array(tab(thisFileInd,'SoundType'))];
        startTimes = [startTimes; seconds(table2array(tab(thisFileInd,'File Offset (s)')))+fileStart];
        endTimes = [endTimes; seconds(table2array(tab(thisFileInd,'File Offset (s)'))+table2array(tab(thisFileInd,'Delta Time (s)')))+fileStart];
        lowFreq = [lowFreq;table2array(tab(thisFileInd,'Low Freq (Hz)'))];
        highFreq = [highFreq;table2array(tab(thisFileInd,'High Freq (Hz)'))];
        
        startSamps = table2array(tab(thisFileInd,'File Offset (s)')).*info.SampleRate;
        endSamps = startSamps + table2array(tab(thisFileInd,'Delta Time (s)')).*info.SampleRate;

        for k = 1:length(thisFileInd)
            % get waveform clip
            clip = double(audioread(soundFiles{j},[round(startSamps(k))-10000,round(endSamps(k))]));
            Waves = [Waves;clip'];
            clip_muPa = clip.*cal_muPa;
            Waves_muPa = [Waves_muPa;clip_muPa];

            % Calculate power spectrum
            nfft = Fs;
            noverlap = floor(nfft*0.75);
            win = hann(nfft);
            [s,f,t,ps] = spectrogram(clip_muPa,win,noverlap,nfft,Fs,'psd');
            Specs = [Specs;mean((abs(s)./nfft)*2,2)'];
            PSD = [PSD;mean(ps,2)'];

        end

    end

end
Specs_dB = 20*log10(Specs);

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
window_length = max(cell2mat(cellfun(@(x) length(x),Waves,'UniformOutput',0)));
for i=1:length(Waves)
    zeroPad = window_length - length(Waves{i,1});
    Waves{i,1} = [Waves{i,1},ones(1,zeroPad)*mean(Waves{i,1})];
end
Waves = cell2mat(Waves);

% Calculate features for each label
allData = struct('Label',[],'Wave',[],'Spectrum',[],'PSD',[],'FreqVec',[],'PeakFreqs',[],'Durations',[],'Bandwidths',[]);
for i=1:length(labs)
    thisLab = strmatch(labs{i},Labels);

    pkFreqs = [];
    BWs = [];
    for j=1:length(thisLab)

        [m,whichMinF] = min(abs(lowFreq(thisLab(j))-f)); % only look at spectrum within freqs described by bounding box
        [m,whichMaxF] = min(abs(highFreq(thisLab(j))-f));
        [valMx,posMx] = max(Specs_dB(thisLab(j),whichMinF:whichMaxF)); % find amplitude and value of peak freq
        posMx = posMx+whichMinF-1; % correct index of maximum to refer to entire freq vector, not subset of freqs in bounding box
        pkFreqs = [pkFreqs;f(posMx)];

        % Calculate -10dB bandwidth
        low = valMx-10; % -10dB value relative to peak
        smoothSpec = fastsmooth(Specs_dB(thisLab(j),:),10);
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

    end
    
    durs = endTimes(thisLab)-startTimes(thisLab);

    allData(i).Label = labs{i};
    allData(i).Wave = Waves(thisLab,:);
    allData(i).Spectrum = Specs(thisLab,:);
    allData(i).PSD = PSD(thisLab,:);
    allData(i).FreqVec = f;
    allData(i).PeakFreqs = pkFreqs;
    allData(i).Durations = durs;
    allData(i).Bandwidths = BWs;
    

end

% Plot mean (filtered) waveform
[b,a] = butter(5,0.0083,'low'); % lowpass filter
filtWaves = filtfilt(b,a,allData(1).Wave');
[b,a] = butter(5,0.00041667,'high'); % highpass filter
filtWaves = (filtfilt(b,a,filtWaves))';
meanWav = mean(filtWaves,1);
meanWav_norm = meanWav-min(meanWav);
meanWav_norm = meanWav_norm./max(meanWav_norm);
figure(1),clf
plot((1:length(meanWav_norm))/Fs,meanWav_norm,'LineWidth',1.5)
xlim([0,3.5])
ylabel('Normalized Amplitude')
xlabel('Time (s)')
set(gca,'FontSize',14)
saveas(gcf,fullfile(saveDir,'MeanWave.png'))

% Plot mean spectra
figure(2),clf
plot(allData(1).FreqVec(3:end),20*log10(mean(allData(1).Spectrum(:,3:end),1)),'LineWidth',1.5)
hold on
% plot(f,quantile(allData(1).Spectrum,0.25),'k:','LineWidth',1.5)
% plot(f,quantile(allData(1).Spectrum,0.75),'k:','LineWidth',1.5)
plot(allData(1).FreqVec(3:end),20*log10(mean(allData(2).Spectrum(:,3:end),1)),'LineWidth',1.5)
hold off
legend(['Sturgeon Thunder, N = ' num2str(length(allData(1).PeakFreqs))],['Background Noise, N = ' num2str(length(allData(2).PeakFreqs))])
xlim([0 450])
ylim([55 105])
ylabel('Amplitude (dB re \muPa)')
xlabel('Frequency (Hz)')
set(gca,'FontSize',14)
saveas(gcf,fullfile(saveDir,'MeanSpecs.png'))

% Plot mean spectrogram
nfft = 2^nextpow2(3000);
win = hann(nfft);
overlap = floor(nfft*0.75);
[s,f,t] = spectrogram(meanWav,win,overlap,nfft,Fs);
s = 20*log10((abs(s)*2).^2);
figure(3)
colormap jet
imagesc(t,f,s)
set(gca,'YDir','normal')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
ylim([0 350])
xlim([0 4])
clim([-130,15])
set(gca,'FontSize',14)
saveas(gcf,fullfile(saveDir,'MeanSpectrogram.png'))

% Plot PSD
P = prctile(allData(1).PSD,[0.01,0.05,0.25,0.5,0.75,0.95,0.99]);
figure(4),clf
plot(allData(1).FreqVec,10*log10(P))
legend({'1%'})

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
boxchart(seconds(allData(1).Durations))
title('Duration');
ylabel('Seconds')
xticks([])
set(gca,'FontSize',14)
saveas(gcf,fullfile(saveDir,'FeatureVariability.png'))


