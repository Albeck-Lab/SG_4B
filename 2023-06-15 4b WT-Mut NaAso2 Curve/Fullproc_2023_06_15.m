addpath('\\albecklab.mcb.ucdavis.edu\data\Notebooks\Nick DeCuzzi\Papers\TF paper with Jessica\Code','\\albecklab.mcb.ucdavis.edu\data\Code\Image Analysis','\\albecklab.mcb.ucdavis.edu\data\Code\Cell Trace','\\albecklab.mcb.ucdavis.edu\data\Code\Nick')

% 3 minutes between each loop!

%% Load the data
% dataloc = JB_DataHandler_V3();

%% Extract and Filter the data
fts(1).t = 'twin'; fts(1).c = 'AreaShape_Area'; fts(1).p = [3,38];
%fts(2).t = 'min'; fts(2).c = 'NumGrans'; fts(2).p = 2;
filterp.fts = fts;
dataloc = JB_DataHandler_V3('extractdata',true,'filterp',filterp);

%% Plot the data
plotme = {'GranMeanGFP','t_granspercell','t_count_cells','NumGrans','GranAreaShape_Area','t_granspercell','t_count_cells','NumGrans','GranAreaShape_Area','GranIntegratedGFP'};
plottype = {'mean'};
plot_by_ND_forJB('treatment', dataloc,'plottype',plottype,'channel',plotme,'looptime',3,'saveloc',dataloc.fold.fig)
plot_by_ND_forJB('cell', dataloc,'plottype',plottype,'channel',plotme,'looptime',3,'saveloc',dataloc.fold.fig)

% Fit the data
plotme = {'t_granspercell','NumGrans'};
plottype = {'albeck mean fit'}; % 
plot_by_ND_forJB('treatment', dataloc,'plottype',plottype,'channel',plotme,'looptime',3,'saveloc',dataloc.fold.fig)
plot_by_ND_forJB('cell', dataloc,'plottype',plottype,'channel',plotme,'looptime',3,'saveloc',dataloc.fold.fig)
