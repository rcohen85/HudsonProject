function noise = dHR_get_noise(candidatesRel,stop,p,hdr)
% Get noise

% minClickSamples = ceil(hdr.fs / 1e6 * p.minClick_us);
maxClickSamples = ceil(hdr.fs  /1e6 * p.maxClick_us);
noise = [];
candidatesRelwEnds = [1,candidatesRel,stop];
dCR = diff(candidatesRelwEnds);
[mC,mI] = max(dCR);
if mC > 2*maxClickSamples
    noiseStart = candidatesRelwEnds(mI)+maxClickSamples;
    noise = [noiseStart, noiseStart+maxClickSamples];
end

if isempty(noise) % if it didn't find any noise, grab some at random.
    noise = [1,maxClickSamples/2];
    %disp('Warning: No noise sample available')
end