% Clear or change all or a subset of labels in selection tables

selTableDir = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\ShortStur\SelectionTables\SS04\It_03';
outDir = ''; % if empty, tables in selTableDir will be overwritten
refresh = 0; % 1 = get rid of existing labels, 0 = keep existing labels; use "retain" below to specify a subset of labels to keep while refreshing the rest
retain = {}; % specify subset of labels to keep when refresh = 1; may be numeric or string
renumber = {3,99}; % change values of specified labels; may be numeric or string; if refresh=1, these labels must also be in "retain"

%%
selList = dir(fullfile(selTableDir,'*.txt'));
if isnumeric(retain)
retain = cellfun(@(x) num2str(x),retain,'UniformOutput',0)
end
if isnumeric([renumber{:}])
renumber = cellfun(@(x) num2str(x),renumber,'UniformOutput',0);
end

for j = 1:length(selList) % for each selection table

    % open selection table
    thisTable = readtable(fullfile(selTableDir,selList(j).name),'Delimiter','\t','VariableNamingRule','preserve');
    if isnumeric(thisTable.Label)
        thisTable.Label = strtrim(cellstr(num2str(thisTable.Label)));
        thisTable.Label = strrep(thisTable.Label,'NaN','');
    end

    if refresh & isempty(retain) % Get rid of all existing labels
        thisTable.Label = [];
    elseif refresh & ~isempty(retain) % Get rid of all BUT the specified labels, which will be renumbered
%         keepInd = strcmp(thisTable.Label,retain);
%         refreshInd = setdiff(1:size(thisTable,1),find(keepInd));
%         thisTable.Label(refreshInd) = '';
        labs = unique(thisTable.Label);
        discard = setdiff(labs,retain);
        thisTable.Label = strrep(thisTable.Label,discard,'');
    end

    if ~isempty(renumber) % Renumber specified labels
        for k = 1:size(renumber,1)
            renumInd = strcmp(thisTable.Label,renumber(k,1));
            thisTable.Label(renumInd) = cellstr(renumber(k,2));
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