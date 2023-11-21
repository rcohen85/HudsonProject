function [xwavNames,matlabDates] = guidedDetection(detFiles,gDxls)
% Use to increase efficiency and accuracy by only running detector over 
% xwav files spanned by a previously defined "detection", requires .xls
% input file, with start/end times of encounters formatted as numbers.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read in .xls file containing 2 columns: interval start, interval times
% you must format your Excel dates/times as NUMBERS, not dates

% read the file into 3 matrices-- numeric, text, and raw cell array
[num, ~, ~] = xlsread(gDxls);

% error check 
[~,y] = size(num);
if y < 2; %start and end dates not formatted as numbers
    error('Dates in guided detection spreadsheet must be saved in NUMBER format');
end  

excelDates = num(:,1:2); % numeric array contains datenums

% convert excel datenums to matlab datenums (different pivot year)
matlabDates = ones(size(excelDates)).*datenum('30-Dec-1899') ...
    + excelDates; % x2mdate does this, but requires financial toolbox

% read xwav headers to determine start of each xwav file
startFile = ones(size(detFiles,1),1);
yr2000 = datenum([2000,00,00]);
fprintf('Reading xwav headers to identify files for guided detection. This may take awhile.\n')
for m = 1:size(detFiles,1)
    thisXwav = detFiles(m,:);
    fileHead = ioReadXWAVHeader(thisXwav);
    startFile(m,1) = fileHead.start.dnum + yr2000;
    endFile(m,1) = fileHead.end.dnum + yr2000;
end
fprintf('Done reading xwav headers.\n')

% in case files are not ordered by time for some reason:
[startFileSort, rIDX] = sortrows(startFile);
endFileSort = endFile(rIDX,:);
detFilesSort = detFiles(rIDX,:); % make a sorted list of file names 

detXwavNames = []; % holder for names of files to process

%take each detection, check which xwav files are associated with the detection
for i = 1:size(matlabDates,1)   %%%%%%%% problems with this logic.
    % find which xwav file(s) correspond(s) with manual detection start 
    thisBoutStart = matlabDates(i,1);
    thisBoutEnd = matlabDates(i,2);

    % find the xwav file that starts before this bout.
    startIdx = find(startFileSort < thisBoutStart,1,'last');
    % if there is none, find the earliest file that matches. If the bout
    % start before recording starts, this could happen.
    if isempty(startIdx)
        startIdx = find(startFileSort > thisBoutStart,1,'first');
    end
    
    endIdx = find(endFileSort > thisBoutEnd,1,'first');
    if isempty(startIdx)
        startIdx = find(endFileSort < thisBoutEnd,1,'last');
    end
    
    if isempty([startIdx, endIdx])
        fprintf('No Recordings during defined detection time #%d \n', i);
    else 
        fIdx = startIdx:endIdx; %combine all indices of files associate with this detection
        fIdx = unique(fIdx);        
        detXwavNames = [detXwavNames; detFilesSort(fIdx,:)];
    end
end

xwavNames = unique(detXwavNames,'rows');