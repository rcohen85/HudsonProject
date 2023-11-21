%% Plot spectrogram of a nice tonal call
% inFile = 'C:\Users\rec297\CCB\Teaching\S1099Hokkaido01_250K_S01_RH403_20210708_081107Z.flac';
% info = audioinfo(inFile);
% [y,Fs] = audioread(inFile,[185.5*info.SampleRate,189*info.SampleRate]);

% inFile = 'C:\Users\rec297\CCB\Teaching\S1099Hokkaido01_250K_S01_RH403_20210708_084107Z.flac';
% info = audioinfo(inFile);
% % [y,Fs] = audioread(inFile,[145*info.SampleRate,149*info.SampleRate]);
% [y,Fs] = audioread(inFile,[153*info.SampleRate,157*info.SampleRate]);

inFile = 'C:\Users\rec297\CCB\Teaching\S1099Hokkaido01_250K_S01_RH403_20210708_091107Z.flac';
info = audioinfo(inFile);
[y,Fs] = audioread(inFile,[192.5*info.SampleRate,195.5*info.SampleRate]);


% filter to attenuate clicks
[b,a] = butter(5,0.0867);
y_filt = filter(b,a,y);

NFFT = 4960;
win = hann(NFFT);
overlap = floor(NFFT*0.75);

[s,f,t] = spectrogram(y_filt,win,overlap,NFFT,info.SampleRate);
s = 20*log10(abs(s));

figure(5),clf
colormap jet
imagesc(t,f,s)
set(gca,'YDir','normal')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
set(gca,'FontSize',14)
ylim([0 13000])
clim([-30 10])

saveas(gcf,'C:\Users\rec297\CCB\Teaching\TonalSpec.png')

%% Plot spectra of a few slices
ts = 126;
[M I] = max(s(5:end,ts));

figure(1),clf
plot(f,s(:,ts),'-','LineWidth',1.5)
hold on
scatter(f(I+4),M,150,'filled','Color','red')
hold off
xlim([0 13000])
xlabel('Frequency (Hz)')
ylabel('Amplitude')
set(gca,'FontSize',14)

saveas(gcf,'C:\Users\rec297\CCB\Teaching\TonalPeakFreq4.png')

%% Plot peak freqs on top of spectrogram

tInds = 110:8:360;
M = [];

for i=tInds
[m ind] = max(s(5:end,i));
M = [M;ind];
end

figure(5),clf
colormap jet
imagesc(t,f,s)
hold on
scatter(t(tInds),f(M+4),60,'white','filled')
hold off
set(gca,'YDir','normal')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
set(gca,'FontSize',14)
ylim([0 13000])
clim([-30 10])

saveas(gcf,'C:\Users\rec297\CCB\Teaching\TracedContour.png')
