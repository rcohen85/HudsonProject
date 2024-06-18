% directory containing sound files
inDir = 'P:\users\cohen_rebecca_rec297\CCB\MarisolTest Data\SoundFiles';
% directory to save selection tables to
outDIr = 'P:\users\cohen_rebecca_rec297\CCB\MarisolTest Data\EnergyDet_SelTabs';
% file extension of audio files
fileExt = '.wav';
% frequency band of interest
freqBounds = [500,24000];
% desired frequency resolution (Hz)
freqRes = 10;
% energy threshold
eThresh = 1e5;
% duration threshold (s)
durThresh = 0.2;

%%
audioList = dir(fullfile(inDir,['**\*',fileExt]));

for i=1:size(audioList,1)

    info = audioinfo(fullfile(inDir,audioList(i).name));
    % Calculate FFT length neede dto obtain desired freq resolution
    NFFT = info.SampleRate./freqRes;
    % create window function
    win = hann(NFFT);
    % overlap of sliding FFT windows (samples)
    noverlap = floor(NFFT*0.75);
    FFT_time = NFFT/info.SampleRate;

    % segment audio data into chunks to avoid reading in too much data at once
    audioSegs = 1:1.5e6:info.TotalSamples;

    % initialize variables to save selection table
    beginTime = [];
    endTime = [];

    for j=1:(size(audioSegs,2)-1)

        segStart = audioSegs(j);

        if j<(size(audioSegs,2)-1)
            [y,Fs] = audioread(fullfile(inDir,audioList(i).name),[audioSegs(j),audioSegs(j+1)-1],'native');
        else
            [y,Fs] = audioread(fullfile(inDir,audioList(i).name),[audioSegs(j),info.TotalSamples],'native');
        end

        [s,f,t] = spectrogram(double(y),win,noverlap,NFFT,info.SampleRate); % calculate spectrogram
        s = 20*log10(abs(s)); % convert FFT output to dB

        if ~isempty(freqBounds)
            [minimum,lowInd] = min(abs(freqBounds(1)-f));
            [minimum,highInd] = min(abs(freqBounds(2)-f));
        else
            lowInd = 1;
            highInd = length(f);
        end

        energySum = sum(s(lowInd:highInd,:),1);

        hotspots = find(energySum>eThresh);
        diffInds = diff(hotspots);
        gapInds = find(diffInds>1);
        gapInds = [1,gapInds,length(diffInds)];
        peakCounts = [];
        for k=1:(numel(gapInds)-1)
            peakCounts = [peakCounts;sum(diffInds(gapInds(k):gapInds(k+1)))];
        end


        figure(3)
        colormap jet
        subplot(2,1,1)
        imagesc(t,f,s) % plot spectrogram
        set(gca,'YDir','normal')
        xlabel('Time (s)')
        ylabel('Frequency (Hz)')
        clim([30,80])
        ylim(freqBounds)
        xlim([0,t(end)])
        subplot(2,1,2)
        plot(t,energySum,'-','LineWidth',1.5) % plot energy sum
        xlim([0,t(end)])
        xlabel('Time (s)')
        ylabel('Energy Sum')

    end


end