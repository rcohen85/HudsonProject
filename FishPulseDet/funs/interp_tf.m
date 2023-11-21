function p = interp_tf(p,tfFullFile)

% If a transfer function is provided, interpolate to desired frequency bins

% Determine the frequencies for which we need the transfer function
p.xfr_f = (p.specRange(1)-1)*p.binWidth_Hz:p.binWidth_Hz:...
    (p.specRange(end)-1)*p.binWidth_Hz;
if ~isempty(tfFullFile)
    [p.xfr_f, p.xfrOffset] = dtf_map(tfFullFile, p.xfr_f);    
else
    % if you didn't provide a tf function, then just create a
    % vector of zeros of the right size.
    p.countThresh = p.countThresh;
    p.ppThresh = 20*log10(2*sqrt(p.countThresh));
    p.xfrOffset = zeros(size(p.xfr_f));
end