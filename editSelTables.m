% Clear or change all or a subset of labels in selection tables

selTableDir = 'P:\users\cohen_rebecca_rec297\CCB\DCLDE2024\BirdNET\Testing\V2.4\20240501\Detections';
outDir = ''; % if empty, tables in selTableDir will be overwritten
labCol = "Species Code";
refresh = 0; % 1 = get rid of existing labels, 0 = keep existing labels; use "retain" below to specify a subset of labels to keep while refreshing the rest
retain = {}; % specify subset of labels to keep when refresh = 1; may be numeric or string
% renumber = {'3','Pm';'5','Zc'}; % change values of specified labels; may be numeric or string; if refresh=1, these labels must also be in "retain"
renumber = {'Upcall','NARW'};
reChan = []; % change channel numbers; first col should be existing channel #, second should be new channel #
%%
selList = dir(fullfile(selTableDir,'**\*.txt'));
if isnumeric(retain)
retain = cellfun(@(x) num2str(x),retain,'UniformOutput',0)
end
if isnumeric([renumber{:}])
renumber = cellfun(@(x) num2str(x),renumber,'UniformOutput',0);
end

for j = 1:length(selList) % for each selection table

    % open selection table
    thisTable = readtable(fullfile(selTableDir,selList(j).name),'Delimiter','\t','VariableNamingRule','preserve');
    if isnumeric(thisTable(:,labCol))
        thisTable.Label = strtrim(cellstr(num2str(thisTable(:,labCol))));
        thisTable.Label = strrep(thisTable(:,labCol),'NaN','');
    end

    if refresh & isempty(retain) % Get rid of all existing labels
        thisTable(:,labCol) = [];
    elseif refresh & ~isempty(retain) % Get rid of all BUT the specified labels, which will be renumbered
%         keepInd = strcmp(thisTable.Label,retain);
%         refreshInd = setdiff(1:size(thisTable,1),find(keepInd));
%         thisTable.Label(refreshInd) = '';
        labs = unique(thisTable(:,labCol));
        discard = setdiff(labs,retain);
        thisTable.Label = strrep(thisTable(:,labCol),discard,'');
    end

    if ~isempty(renumber) % Renumber specified labels
        for k = 1:size(renumber,1)
            renumInd = strcmp(table2array(thisTable(:,labCol)),renumber{k,1});
            thisTable(renumInd,labCol) = cellstr(renumber(k,2));
        end
    end

    if ~isempty(reChan) % Change channel numbers as specified
        for k=1:size(reChan,1)
            reChanInd = find(thisTable.Channel==reChan(k,1));
            thisTable.Channel(reChanInd) = reChan(k,2);
        end
    end

    % Save updated selection table
    if isempty(outDir)
        saveName = fullfile(fullfile(selTableDir,selList(j).name));
    else
        saveName = fullfile(fullfile(outDir,selList(j).name));
    end
    writetable(thisTable,saveName,'Delimiter','\t');

    thisTable = [];

end