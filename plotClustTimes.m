inDir = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\ShortStur\Clusters\compositeClusters_combined\IT02';
% inDir = 'W:\projects\2022_NOAA-NERRS_HudsonNY_144488\ShortStur\Clusters\SS03\CompositeClusters\It_02';
site = 'Shortnose Sturgeon Overwintering Site';
types = [1,2,3,4,7,8];
% Shortnose sturgeon effort
effort = [datenum({'2021/11/19 11:28:17';'2022/03/09 07:09:07';'2022/10/12 13:17:04'}),datenum({'2021/12/09 06:25:45';'2022/03/29 01:51:05';'2022/10/27 05:39:38'})];

%%
load(fullfile(inDir,'_types_all.mat'));

% N = size(Tfinal,1);
N = length(types);
% names = flipud(join([repmat('Cluster',N,1),(string(1:N))']));
names = flipud(join([repmat('Cluster',N,1),(string(types))']));
figure(99);clf
starts = [];
ends = [];
hold on
for i = 1:N
    starts = [starts;min(Tfinal{types(i),7})];
    ends = [ends;max(Tfinal{types(i),7})];
scatter(Tfinal{types(i),7},ones(length(Tfinal{types(i),7}),1)*(N-i+1),60,'filled','MarkerFaceAlpha',0.4,'MarkerEdgeAlpha',0.3)
end
for i = 1:(length(effort)-1)
    noEff = patch([effort(i,2),effort(i,2),effort(i+1,1),effort(i+1,1)],[N,0,0,N],[128,128,128]./255);
    noEff.FaceAlpha = 0.2;
    noEff.EdgeColor = 'none';
end
hold off
grid on
grid minor
yticks([1:N])
yticklabels(names)
xlim([min(starts),max(ends)])
datetick('x',12)
xlim([min(starts),max(ends)])
title(site)
