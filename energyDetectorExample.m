% Simple implementation of a basic energy detector

% Example sound file to read in
inFile = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AtlStur\SoundsOfInterest\AS03_20210614_145138_AtlStur_FreshDrum_LongerClip.wav';
% If desired, select frequency bounds to restrict energy summation to just a region of the spectrogram; otherwise, leave empty
freqBounds = [0,200];

%%

info = audioinfo(inFile); % get info about file
% if it's a long file, just read a segment where the sounds of interest are known to be present
% [y,Fs] = audioread(inFile,[(60*info.SampleRate)+1,(2*60*info.SampleRate)+1],'native');

% if it's a short file, read the whole thing
[y,Fs] = audioread(inFile,'native'); 

NFFT = 4500; % FFT length (samples)
win = hann(NFFT); % create window function
noverlap = floor(NFFT*0.75); % overlap of sliding FFT windows (samples)

[s,f,t] = spectrogram(double(y),win,noverlap,NFFT,info.SampleRate); % calculate spectrogram
s = 20*log10(abs(s)); % convert FFT output to dB

if ~isempty(freqBounds)
[minimum,lowInd] = min(abs(freqBounds(1)-f));
[minimum,highInd] = min(abs(freqBounds(2)-f));
else
    lowInd = 1;
    highInd = length(f);
end

energySum = sum(s(lowInd:highInd,:),1); % sum energy in each spectrogram slice in freq range of interest

figure(3)
colormap jet
subplot(2,1,1) 
imagesc(t,f,s) % plot spectrogram
yline(0,'-m','LineWidth',3) % add lines demarcating frequency range energy was summed within
yline(200,'-m','LineWidth',3)
set(gca,'YDir','normal')
xlabel('Time (s)')
ylabel('Frequency (Hz)'),
clim([45,100])
ylim([0,2000])
xlim([0,t(end)])
subplot(2,1,2)
plot(t,movmean(energySum,5),'-','LineWidth',1.5) % plot energy sum
yline(2500,'-r','LineWidth',2) % add line denoting threshold for "present"
xlim([0,t(end)])
xlabel('Time (s)')
ylabel('Energy Sum')

saveas(gcf,'C:\Users\rec297\CCB\Teaching\EnergySumDetector.png')
