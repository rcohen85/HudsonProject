% Sturgeon thunder spectrogram cross-correlation detector

% Create template
inFile = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_FLAC\AtlSturSpw\ST300\AS03\144488Hudson_024K_AS03_ST300-5683_20210609_234043Zm0500.flac';
info = audioinfo(inFile);
st = 41139*info.SampleRate;
ed = st + 2*info.SampleRate;
[y,Fs] = audioread(inFile,[st,ed]);

% Get call to use as specctrogram template
[b,a] = butter(7,0.013);
cleanCall = filter(b,a,y);

NFFT = 2048;
win = hann(NFFT);
overlap = floor(NFFT*0.75);

[s,f,t] = spectrogram(cleanCall,win,overlap,NFFT,info.SampleRate);
s = 20*log10(abs(s));

figure(1)
colormap jet
imagesc(t,f,s)
set(gca,'YDir','normal')
ylim([0,1000])
colorbar
clim([-100,0])


%% Test data to determine threshold to trigger detection
testFile = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_FLAC\AtlSturSpw\ST300\AS03\144488Hudson_024K_AS03_ST300-5683_20210609_234043Zm0500.flac';
info = audioinfo(testFile);
st = 43140*info.SampleRate;
ed = st + 180*info.SampleRate;
[y,Fs] = audioread(testFile,[st,ed]);

NFFT = 2048;
win = hann(NFFT);
overlap = floor(NFFT*0.75);

[s,f,t] = spectrogram(y,win,overlap,NFFT,info.SampleRate);
s = 20*log10(abs(s));

win = 801;
NFFT = 256;
noverlap = floor(win*0.75);

S = [];
T = [];
k = 1;
for i=1:(length(y)-win)
s = xspectrogram(cleanCall,y(i:(i+win-1)),win,noverlap,NFFT);
S = [S;max(s(2:end))];
end

figure(5),clf
subplot(2,1,1)
imagesc(t,f,s)
set(gca,'YDir','normal')
ylim([0,1000])
colorbar
clim([-100,0])
subplot(2,1,2)
plot(S,'-','LineWidth',1.5)
yline(0.2e-4,'-r','LineWidth',1.5)
xlabel('Lag')
ylabel('Cross-Correlation')
xlim([0 length(S)])
set(gca,'FontSize',14)



%% Loop through data and perform spectrogram cross-correlation

inDir = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_FLAC\AtlSturSpw\ST300\AS03'; % directory containing audio data (.wav or .flac)
outDir = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AtlStur\ThunderClassifier_Training\SpectroXCorr'; % directory to save 
readSegments = 3; % duration (mins) of audio data segments to successively read in
