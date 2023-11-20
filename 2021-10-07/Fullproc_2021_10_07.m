
%% Process the data
fts(1).t = 'Min'; fts(1).c = 'Intensity_IntegratedIntensity_GFP'; fts(1).p = 0;
% fts(2).t = 'Max'; fts(2).c = 'EKAR'; fts(2).p = 1;
% fts(3).t = '25prctmin'; fts(3).c = 'YFP_Nuc'; fts(3).p = 300;
filterp.dlength = 5;
filterp.fts = fts;

%% Load the data
%dataloc = JB_DataHandler(); %reorganizedata

%filter the data
dataloc = JB_DataHandler_V2('reorganizedata', true, 'filterp',filterp);%'reorganizedata',true

%% Plot the data
plotme = {'Intensity_IntegratedIntensity_GFP', 'Intensity_MeanIntensity_GFP', 'Intensity_StdIntensity_GFP'};
plot_by_ND('treatment', dataloc,'plottype',{'meanslope', 'heatmap'},'channel',plotme,'looptime',4)

