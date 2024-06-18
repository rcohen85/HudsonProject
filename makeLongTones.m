[y,fs] = audioread('/Users/rec297/Downloads/Tonal_sounds_0_20.wav');
longTones = [];
step = 250;
samp = 132301;

for i = 1:32

    thisTone = y(samp:(samp+(fs*2)-1));
    longTones = [longTones;thisTone;thisTone];
    samp = samp+(fs*3);

end


NFFT = 3500;
win = hann(NFFT);
noverlap = floor(NFFT*0.75);

[s,f,t] = spectrogram(longTones,win,noverlap,NFFT,fs);
s = 20*log10(abs(s));

imagesc(t,f,s)
set(gca,'YDir','normal')

audiowrite('LongTones_250Hz_8kHz.wav',longTones,fs);