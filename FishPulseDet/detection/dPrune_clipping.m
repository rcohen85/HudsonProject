function validClicks = dPrune_clipping(clicks,p,hdr,teagerBandData)
% Perform a series of tests to retain only promising clicks:

validClicks = ones(1, size(clicks,1));  % assume all are good to begin

for c = 1:size(clicks, 1);
    % Go through the segments of interest one by one, and make sure they
    % don't exceed clip threshold.
    if any(abs(teagerBandData(clicks(c,1):clicks(c,2))) > p.clipThreshold *(2^hdr.nBits)/2)
        validClicks(c) = 0;
    end
end