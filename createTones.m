outFile = '/Users/rec297/Documents/CCB/WRI_2024/ShortTones_ManAmpMod_30mins.wav';
Fs = 10000;
dur = 0.2; % duration of each tone (s)
freqs = [150,300,1000,2000,3000]; % desired frequencies
x = [120.311, 120.311, 111.076, 147.724, 135.036]; % observed output amplitude at each frequency
bitDepth = 16;
repeat = 900;

x_lin = power (10, x ./ 20);
x_lin_prop = x_lin./max(x_lin);
q_fac = 1./x_lin_prop;
q_fac(q_fac>1) = q_fac(q_fac>1)./2;
q = q_fac./max(q_fac);
q = [0.9,0.9,1,0.1,0.12];

y=zeros(Fs/2,1);

for i=1:length(freqs)
sine = dsp.SineWave('Amplitude',q(i),'Frequency',freqs(i),'SampleRate',Fs,'SamplesPerFrame',Fs*dur);
ySub = sine();
y=[y;ySub(1:end-1);0];
end

y = [y;zeros(Fs/2,1)];

y = repmat(y,repeat,1);

figure(3),plot(y,'*-')
audiowrite(outFile,y,Fs,'BitsPerSample',bitDepth);