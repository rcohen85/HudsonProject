
inDir = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\ShortStur\SelectionTables\NorriePoint_Trains\Hydro';
fileStr = '*Train_selections.txt';

fileList = dir(fullfile(inDir,fileStr));
Trains = struct();
for i = 1:length(fileList)
    
    % read selection table
    tb = readtable(fullfile(inDir,fileList(i).name),'Delimiter','tab','VariableNamingRule','preserve');

    % remove duplicate views
    specs = find(strcmp(table2array(tb(:,'View')),'Spectrogram 1'));
    tb = tb(specs,:);

    % get times
    timeStrs = regexp(table2array(tb(:,'Begin File')),'\d{8}_\d{6}','match');
    fileStarts = (datetime([timeStrs{:}],'InputFormat','yyyyMMdd_HHmmss'))';

    % detection start datetimes
    detSt = fileStarts + seconds(table2array(tb(:,'Begin Time (s)')));

    % get durations
    durs = table2array(tb(:,'End Time (s)')) - table2array(tb(:,'Begin Time (s)'));

    % detection end datetimes
    detEnd = detSt + seconds(durs);

    % detection start times
    stTime = datetime(0,0,0,hour(detSt),minute(detSt),second(detSt));

    % detection end times
    endTime = datetime(0,0,0,hour(detEnd),minute(detEnd),second(detEnd));

    if i==1
        Trains.StartDT = detSt;
        Trains.EndDT = detEnd;
        Trains.StartTime = stTime;
        Trains.EndTime = endTime;
        Trains.Dur = durs;
    else
        Trains.StartDT = [Trains.StartDT;detSt];
        Trains.EndDT = [Trains.EndDT;detEnd];
        Trains.StartTime = [Trains.StartTime;stTime];
        Trains.EndTime = [Trains.EndTime;endTime];
        Trains.Dur = [Trains.Dur;durs];
    end

end

figure(99)
scatter(Trains.StartTime,Trains.Dur,40,'MarkerFaceColor','b','MarkerFaceAlpha',0.3,'MarkerEdgeColor','b','MarkerEdgeAlpha',0.3)
xticklabels({'00:00','06:00','12:00','18:00','24:00'})
ylabel('Duration (s)')
title('Train Transit Durations at Norrie Point (Hydrophone)')
