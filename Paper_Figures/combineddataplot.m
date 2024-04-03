addpath('\\albecklab.mcb.ucdavis.edu\data\Code\Nick','N:\Code\Cell Trace')
%run it all together
dataloc{1} = "\\albecklab.mcb.ucdavis.edu\data\Notebooks\Nick DeCuzzi\Papers\TF paper with Jessica\Experiments\2022-11-02\2022-11-02_Processed.mat";
dataloc{2} = "\\albecklab.mcb.ucdavis.edu\data\Notebooks\Nick DeCuzzi\Papers\TF paper with Jessica\Experiments\2022-10-27\2022-10-27_Processed.mat";
%dataloc{3} = "\\albecklab.mcb.ucdavis.edu\data\Notebooks\Nick DeCuzzi\Papers\TF paper with Jessica\Experiments\2022-11-03 4Brrmmut\2022-11-03 4Brrmmut_Processed.mat";

datalocDF = makeLiveCellDataframe(dataloc);

plotme = {'granspercell'}; %,'granspercell'
plottype = {'albeck mean fit fixed f'}; % 'albeck mean fit'

plot_by_ND_forJB('treatment', datalocDF,'plottype',plottype,'channel',plotme,'looptime',3,'font_size',8)
