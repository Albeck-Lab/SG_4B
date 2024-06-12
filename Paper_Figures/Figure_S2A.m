addpath('Z\data\Code\Nick')
%run it all together
dataloc{1} = "Z:\imageData\SG_4B\2022-11-02\2022-11-02_Processed.mat";
dataloc{2} = "Z:\imageData\SG_4B\2022-10-27\2022-10-27_Processed.mat";
dataloc{3} = "Z:\imageData\SG_4B\2022-11-03 4Brrmmut\2022-11-03 4Brrmmut_Processed.mat";

datalocDF = makeLiveCellDataframe(dataloc,'exclude','washoff');

plotme = {'granspercell'}; %,'granspercell'
plottype = {'albeck mean fit fixed f'}; % 'albeck mean fit'

plot_by_ND_forJB('treatment', datalocDF,'plottype',plottype,'channel',plotme,'looptime',3,'font_size',8,'tmaxaftertx',1,'saveloc','Z:\imageData\SG_4B\Paper_Figures\Output_Figures','closefigs',false,'exclude','washoff')
cFig = gcf;
set(cFig,'Units','inches','Position',[0.5,0.5,3,3])
cAx = gca;
xlabel('minutes')
xticks([2     6    10    14    18    22]) % tp 2 aka 6 min,etc
xticklabels(([2     6    10    14    18    22]-2)*3) % convert to minutes
lgy = legend;
lgy.String={'WT','M2','M1','M3'}; % was WT; 99,102,135; 139; 135,137,139

saveas(cFig,'Z:\imageData\SG_4B\Paper_Figures\Output_Figures\Figure_S2A.svg')
saveas(cFig,'Z:\imageData\SG_4B\Paper_Figures\Output_Figures\Figure_S2A.fig')

close all