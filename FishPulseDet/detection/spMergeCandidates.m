function [c_stopsM ,c_startsM] = spMergeCandidates(mergeThr,fs,stops,starts)

% merge candidates that are too close together so they will be considered
% to be one larger detection
minGap_samples = ceil(mergeThr*fs/1e6);
c_startsM(1) = starts(1);
c_stopsM = [];
itr = 1;
startCtr = 2;
stopCtr = 1;
while itr <= length(starts)
    k = 0;
    merg = itr;
    while (k+itr)<(length(starts)) && starts(itr+k+1,1) - stops(itr+k,1)< minGap_samples    
        k = k+1;
        merg = [merg,itr+k];
    end
    if(k+itr)==(length(starts))
         c_stopsM(stopCtr,1) = max(stops(merg));
         break
    else
        c_startsM(startCtr,1)= starts(itr+k+1,1);
        c_stopsM(stopCtr,1) = stops(itr+k);
        itr = itr+1+k;
        startCtr = startCtr + 1;
        stopCtr = stopCtr + 1;
        k=0;
    end
end