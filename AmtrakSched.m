% Amtrak Train Schedule; times trains arrive at (northbound) or depart from
% (southbound) the Rhinecliff station

inDat = sortrows({'01:13','N','261 Empire';
'07:09','S','250 Empire';
'07:41','S','234 Empire';
'08:54','S','260 Empire';
'08:58','N','63 Maple';
'09:09','S','236 Empire'
'10:28','N','69 Adirondack';
'10:51','S','280 Empire';
'12:00','N','281 Empire';
'12:58','N','233 Empire';
'14:58','N','283 Empire';
'15:58','N','291 Ethan Allen';
'16:58','N','235 Empire';
'17:24','N','49 Lake Shore';
'18:15','N','237 Empire';
'18:56','N','253 Empire';
'19:36','N','239 Empire';
'20:55','N','241 Empire';
'22:36','N','243 Empire';
'22:55','N','259 Empire';
'11:54','S','240 Empire';
'12:54','S','238 Empire;';
'13:51','S','284 Empire';
'15:54','S','290 Ethan Allen';
'16:55','S','48 Lake Shore';
'17:14','S','244 Empire';
'18:20','S','246 Empire';
'19:20','S','64 Maple';
'20:25','S','68 Adirondack'});

Trains = struct('Times',[],'Dir',[],'Route',[]);
for i = 1:length(inDat)

    Trains(i).Times = datetime(inDat{i,1},'InputFormat','HH:mm');
    Trains(i).Dir = inDat{i,2};
    Trains(i).Route = inDat{i,3};

end

%% Plot daily train occurrence at Rhinecliff

figure(88)
scatter(vertcat(Trains.Times),categorical(cellstr(vertcat(Trains.Dir))),60,[255,0,127]/255,'filled','MarkerEdgeColor','b')
% scatter(vertcat(Trains.Times),ones(length(Trains),1),60,[255,0,127]/255,'filled','MarkerEdgeColor','b')
% scatter(vertcat(Trains.Times),ones(length(Trains),1),60,categorical(cellstr(vertcat(Trains.Dir))),'filled','MarkerEdgeColor','m')
% xticklabels({'00:00','06:00','12:00','18:00','24:00'})
xticks([dateshift(datetime(),'start','day')+hours([0:2:24])]);
xticklabels({'00:00','02:00','04:00','06:00','08:00','10:00','12:00','14:00','16:00','18:00','20:00','22:00','24:00'})
set(gca,'FontSize',16)
grid on
