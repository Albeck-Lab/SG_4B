addpath('\\albecklab.mcb.ucdavis.edu\data\Notebooks\Nick DeCuzzi\Papers\TF paper with Jessica\Code','\\albecklab.mcb.ucdavis.edu\data\Code\Image Analysis','\\albecklab.mcb.ucdavis.edu\data\Code\Cell Trace','\\albecklab.mcb.ucdavis.edu\data\Code\Nick')

% 3 minutes between each loop!

%% Load the data
% dataloc = JB_DataHandler_V3();

%% Extract the data
dataloc = JB_DataHandler_V3('extractdata',true);

%% Filter the data
fts(1).t = 'anymax'; fts(1).c = 'AreaShape_Area'; fts(1).p = 3;
filterp.dlength = 40;
filterp.fts = fts;
dataloc = JB_DataHandler_V3('dataloc',dataloc,'filterp',filterp);


%% Plot the data
%plotme = {'Intensity_IntegratedIntensity_GFP', 'Intensity_MeanIntensity_GFP', 'Intensity_StdIntensity_GFP'};
%plotme = {'granspercell','AreaShape_Area'};
%plot_by_ND('cell', dataloc,'plottype',{'meanslope','mean'},'channel',plotme,'looptime',3,'combinexys',true)

plotme = {'t_granspercell','t_count_cells','NumGrans'};
%plotme = {};
plottype = {'mean'}; % 
plot_by_ND_forJB('treatment', dataloc,'plottype',plottype,'channel',plotme,'looptime',3,'saveloc',dataloc.fold.fig,'subset','24')
%plot_by_ND_forJB('treatment', dataloc,'plottype','sum mean','channel','count_cells','looptime',3,'saveloc',dataloc.fold.fig,'subset','24')

%plot_by_ND_forJB('cell', dataloc,'plottype',plottype,'channel',plotme,'looptime',3,'saveloc',dataloc.fold.fig)


