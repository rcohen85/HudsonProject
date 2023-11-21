function [hdr] = dInput_HR_files(fullFile,p)
% figure out how to read header info


% Is it a wav or an xwav file?
[fileType, ~] = ioGetFileType(fullFile);

% Retrieve header information for this file
if fileType == 1
    hdr = ioReadWavHeader(fullFile, p.DateRE);
else
    hdr = ioReadXWAVHeader(fullFile, 'ftype', fileType);
end
