inFile = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AcousticData_FLAC\ShrtSturWnt\ST300\SS03\144488Hudson_024K_SS03_ST300-5683_20220325_080438Zm0500.flac';
info = audioinfo(inFile);
startTime = datetime(2022,03,25,08,04,38);
timeOfInterest = datetime(2022,03,25,12,53,45);

timeDiff = seconds(timeOfInterest-startTime);

[y,Fs] = audioread(inFile,[info.SampleRate*timeDiff,(info.SampleRate*timeDiff)+(info.SampleRate*15)],'native');
y = double(y);
y = [y;y(1:5*info.SampleRate);;y(1:5*info.SampleRate)]; % add some noise at end of pulse train
y_norm = y-min(y);
y_norm = y_norm./max(y_norm);

NFFT = 512;
win = hann(NFFT);
noverlap = floor(NFFT*0.75);

[s,f,t] = spectrogram(double(y),win,noverlap,NFFT,info.SampleRate);
s = 20*log10(abs(s));

energySum = sum(s(87:129,:),1);
winLen = 100;
winAC = [];
tAC = [];
k=1;
for i=1:floor(length(t)/winLen)
[acf,lags] = autocorr(energySum(k:k+winLen-1),'NumLags',98);
winAC = [winAC; acf(2:end)];
tAC = [tAC;t(k)];
k=k+winLen;
end

figure(90),clf
subplot(4,1,1)
plot([1:length(y_norm)]/info.SampleRate,y_norm)
xlabel('Time (s)')
ylabel('Normalized Amplitude')
xlim([0,25])
subplot(4,1,2)
imagesc(t,f,s)
set(gca,'YDir','normal')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
clim([0,50])
xlim([0,25])
subplot(4,1,3)
plot(t,energySum)
xlim([0,25])
xlabel('Time (s)');
ylabel('Energy Sum')
subplot(4,1,4)
plot(tAC,max(abs(winAC),[],2),'-','LineWidth',1.5)
yline(0.6,'-r','LineWidth',2)
xlabel('Time (s)')
ylabel('Detection Function')
xlim([0,25])

saveas(gcf,'C:\Users\rec297\CCB\Teaching\SequenceDetector.png')

figure(5)
plot(lags(2:end),abs(winAC(21,:)),'-','LineWidth',1.5)
xlabel('Lag')
ylabel('Autocorrelation')
xlim([0,98])