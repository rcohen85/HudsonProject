function percLabeled(inDir,varargin)

% inDir: directory containing selection tables with numeric labels in a 
% column named "Labels"
% varargin: name-value pair 'junkLabels' and a 1xN cell array containing 
% the labels to be considered junk classes

for i = 1:2:length(varargin)
    if ischar(varargin{i}) && strmatch('junkLabels',varargin{i})
        sortJunk = 1;
        badLabs = varargin{i+1};
    else
        sortJunk = 0;
    end
end

fileList = dir(fullfile(inDir,'*.txt'));
labels = [];
count = 0;

for i = 1:size(fileList,1)

    thisTable = readtable(fullfile(inDir,fileList(i).name),'Delimiter','\t','VariableNamingRule','preserve');
    names = thisTable.Properties.VariableNames;
    if any(strcmp(names,'Label'))
        labels = [labels;thisTable.Label];
    else
        labels = [labels;cell(size(thisTable,1),1)];
    end
    count = count + size(thisTable,1);
    thisTable = [];
end

if isnumeric(labels)
    labels = strtrim(cellstr(num2str(labels)));
    labels = strrep(labels,' 0','');
    labels = strrep(labels,'NaN','');
end

noLabInds = cell2mat(cellfun(@(x) isempty(x),labels,'UniformOutput',0));
labInds = setdiff(1:length(labels),find(noLabInds));
propLabeled = size(labInds,2)/size(labels,1);
fprintf('%0.2f%% of detections labeled, %0.2f%% are unlabeled (%d and %d of %d)\n',round(propLabeled*100,2),round((1-propLabeled)*100,2),size(labInds,2),sum(noLabInds),size(labels,1));

if sortJunk
    junkInd = strcmp(labels,badLabs); % find indices where labels match bad labels specified by user
    notJunkLabs = setdiff(labInds,find(junkInd));
    propJunkLabels = sum(junkInd)/size(labInds,2);
    propGoodLabels = size(notJunkLabs,2)/size(labInds,2);
    fprintf('%0.2f%% of labels are good, %0.2f%% of labels are junk (%d and %d of %d)\n',round(propGoodLabels*100,2),round(propJunkLabels*100,2),size(notJunkLabs,2),sum(junkInd),size(labInds,2));

    % find indices where labels exist but do not match junk labels
    notJunk = (setdiff(1:size(labels,1),find(junkInd)))';
    notJunk_butLabeled = cell2mat(cellfun(@(x) ~isempty(x),labels(notJunk),'UniformOutput',0));

    propNotJunk_butLabeled = sum(notJunk_butLabeled)/size(notJunk,1);
    fprintf('%0.2f%% of non-junk detections labeled (%d of %d)\n',round(propNotJunk_butLabeled*100,2),sum(notJunk_butLabeled),size(notJunk,1));
    
end
end