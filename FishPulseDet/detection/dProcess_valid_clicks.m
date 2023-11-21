function [clkStart,clkEnd]= dProcess_valid_clicks(clicks,clickInd,startsK,hdr,...
    fidOut,p)
% Write click times to .ctg label file

clkStart = nan(length(clickInd),1);
clkEnd = nan(length(clickInd),1);

for c = 1:length(clickInd)
    cI = clickInd(c);
    currentClickStart = startsK + clicks(cI,1)/hdr.fs; % start time
    currentClickEnd = startsK + clicks(cI,2)/hdr.fs;
    if ~isempty(currentClickEnd)
        % If true write individual click annotations to .ctg file
        fprintf(fidOut, '%f %f \n', ...
            currentClickStart + length(p.fB)/2/hdr.fs, ...
            currentClickEnd + length(p.fB)/2/hdr.fs);
    end
    
    % Compute parameters of click frames
    clkStart(c,1) = currentClickStart;
    clkEnd(c,1) = currentClickEnd;
    
end
