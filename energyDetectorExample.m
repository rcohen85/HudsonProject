inFile = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\AtlStur\SoundsOfInterest\AS03_20210610_133825_AtlStur.wav';
info = audioinfo(inFile);
[y,Fs] = audioread(inFile,[(60*info.SampleRate)+1,(2*60*info.SampleRate)+1],'native');

NFFT = 4500;
win = hann(NFFT);
noverlap = floor(NFFT*0.75);

[s,f,t] = spectrogram(double(y),win,noverlap,NFFT,info.SampleRate);
s = 20*log10(abs(s));

energySum = sum(s(1:39,:),1);

figure(3)
colormap jet
subplot(2,1,1)
imagesc(t,f,s)
yline(0,'-m','LineWidth',3)
yline(200,'-m','LineWidth',3)
set(gca,'YDir','normal')
xlabel('Time (s)')
ylabel('Frequency (Hz)'),
clim([45,100])
ylim([0,2000])
subplot(2,1,2)
plot(t,movmean(energySum,5),'-','LineWidth',1.5)
yline(2600,'-r','LineWidth',2)
xlabel('Time (s)')
ylabel('Energy Sum')

saveas(gcf,'C:\Users\rec297\CCB\Teaching\EnergySumDetector.png')
