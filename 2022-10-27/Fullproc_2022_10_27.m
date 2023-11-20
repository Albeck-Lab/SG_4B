addpath('\\albecklab.mcb.ucdavis.edu\data\Notebooks\Nick DeCuzzi\Papers\TF paper with Jessica\Code','\\albecklab.mcb.ucdavis.edu\data\Code\Image Analysis','\\albecklab.mcb.ucdavis.edu\data\Code\Cell Trace','\\albecklab.mcb.ucdavis.edu\data\Code\Nick')

% 3 minutes between each loop!

%% Load the data
dataloc = JB_DataHandler_V2();

%% Extract the data
%dataloc = JB_DataHandler_V2('extractdata',true);

%% Filter the data
% fts(1).t = 'Min'; fts(1).c = 'Intensity_IntegratedIntensity_GFP'; fts(1).p = 0;
% fts(2).t = 'Max'; fts(2).c = 'EKAR'; fts(2).p = 1;
% fts(3).t = '25prctmin'; fts(3).c = 'YFP_Nuc'; fts(3).p = 300;
% filterp.dlength = 5;
% filterp.fts = fts;

%% Plot the data
%plotme = {'Intensity_IntegratedIntensity_GFP', 'Intensity_MeanIntensity_GFP', 'Intensity_StdIntensity_GFP'};
%plotme = {'granspercell','AreaShape_Area'};
%plot_by_ND('cell', dataloc,'plottype',{'meanslope','mean'},'channel',plotme,'looptime',3,'combinexys',true)

% plotme = {'granspercell'};
% 
% plottype = {'mean','albeck mean fit'}; % 'albeck mean fit'
% plot_by_ND('treatment', dataloc,'plottype',plottype,'channel',plotme,'looptime',3,'combinexys',true)
% plot_by_ND('cell', dataloc,'plottype',plottype,'channel',plotme,'looptime',3,'combinexys',true)
% 

