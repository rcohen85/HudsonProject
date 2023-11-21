function hdr = ioGetFlacInfo(Filename, DateRE)
% Adapted from ioReadWavHeader.m
%
% Read header of Microsoft RIFF wav header
% See http://www.sonicspot.com/guide/wavefiles.html for
% layout of Microsoft RIFF wav files
%
% CAVEATS:  Assumes a single DATA chunk.
% To modify to handle multiple data chunks, be sure to also
% consider ioReadWav which will need modifications as well.
%
% Attempts to infer the timestamp of the recording based
% upon the filename and the regular expression(s) DateRE
% which must conform to the standards in function dateregexp.
%
% Do not modify the following line, maintained by CVS
% $Id: ioReadWavHeader.m,v 1.6 2008/12/09 19:35:38 mroch Exp $

global PARAMS

error(nargchk(1,2,nargin));
if nargin < 2
    % Use global timestamp if available
    if exist('PARAMS', 'var') && isfield(PARAMS, 'fnameTimeRegExp')
        DateRE = PARAMS.fnameTimeRegExp;
    else
        DateRE = [];
    end
end

% Check file type
f_handle = fopen(Filename, 'rb', 'b'); % open big-endian binary file
[id, bytes] = fread(f_handle, 4,'char');

if f_handle == -1
    error('io:Unable to open file %s', Filename);
end

fclose(f_handle);
fType = deblank(char(id'));

if strmatch(fType,'fLaC')
    hdr.fType = ['.' lower(fType)];
else
    error('%s is not a .flac file', Filename)
end

fileInfo = audioinfo(Filename);

if fileInfo.TotalSamples==1
    error('Empty file')
end

hdr.fs = fileInfo.SampleRate;
hdr.nch = fileInfo.NumChannels;
hdr.nbits = fileInfo.BitsPerSample;
hdr.samp.byte = fileInfo.BitsPerSample/8;
hdr.xhd.ByteRate = [];
hdr.xhd.byte_length = [];
hdr.xhd.byte_loc = [];

% Add HARP data structures for uniform access
hdr.xgain = 1;          % gain (1 = no change)
[~,shortName,~] = fileparts(Filename);
% determine timestamp
[~,~,~,~,k] = regexp(shortName, DateRE);
catDate = cell2mat(k{1});
hdr.start.dvec = [str2double(catDate(1:4)),str2double(catDate(5:6)),...
    str2double(catDate(7:8)),str2double(catDate(9:10)),...
    str2double(catDate(11:12)),str2double(catDate(13:14))];
hdr.start.dnum = datenum(hdr.start.dvec);
hdr.xhd.year = hdr.start.dvec(1);          % Year
hdr.xhd.month = hdr.start.dvec(2);         % Month
hdr.xhd.day = hdr.start.dvec(3);           % Day
hdr.xhd.hour = hdr.start.dvec(4);          % Hour
hdr.xhd.minute = hdr.start.dvec(5);        % Minute
hdr.xhd.secs = hdr.start.dvec(6);          % Seconds

samplesN = fileInfo.TotalSamples/fileInfo.NumChannels;
hdr.end.dnum = hdr.start.dnum + datenum([0 0 0 0 0 samplesN/hdr.fs]); % THIS MIGHT NOT BE 100% RIGHT!

hdr.start.dvec = datevec(hdr.start.dnum);




