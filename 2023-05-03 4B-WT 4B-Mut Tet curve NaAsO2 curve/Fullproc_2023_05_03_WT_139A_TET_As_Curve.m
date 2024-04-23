addpath('\\albecklab.mcb.ucdavis.edu\data\Code\Image Analysis','\\albecklab.mcb.ucdavis.edu\data\Code\Cell Trace','\\albecklab.mcb.ucdavis.edu\data\Code\Nick','\\albecklab.mcb.ucdavis.edu\data\imageData\SG_4B\Code\')

% 3 minutes between each loop!

%% Load the data
%dataloc = JB_DataHandler_V3();

%% Extract the data
fts(1).t = 'anymax'; fts(1).c = 'AreaShape_Area'; fts(1).p = 3;
filterp.dlength = 40;
filterp.fts = fts;
dataloc = SG_Datahandler('dataloc',dataloc,'extractif',true);

%% Filter the data
goodD = all([dataloc.ifd.Children_Grans_4B_Count >= 3, dataloc.ifd.Children_Grans_G3BP1_Count >= 3, contains(dataloc.ifd.treatment,'hour -24'), contains(dataloc.ifd.treatment,'and 62')],2);
subD = dataloc.ifd(goodD,:);
boxplot(subD.Children_Grans_G3BP1_Count,subD.cell,'Notch','on','Symbol','')
% xlabel('4B'); ylabel('G3BP1')
% xlim([0,30]); ylim([0,30])
grpstats(subD,"cell","mean","DataVars",["Children_Grans_4B_Count","Children_Grans_G3BP1_Count"])

[~,~,statsF] = anova1(subD.Children_Grans_G3BP1_Count,subD.cell,'off');
[resultsMaxG,~,~,gnamesF] = multcompare(statsF,"CriticalValueType","dunnett",'ControlGroup',find(matches(statsF.gnames,'HeLa eIF4BGFP')),'Display','off','Approximate',false); 
MaxGranulesFormed = array2table(resultsMaxG,"VariableNames", ["Group","Control Group","Lower Limit","Difference","Upper Limit","P-value"]);
MaxGranulesFormed.("Group") = gnamesF(MaxGranulesFormed.("Group"));
MaxGranulesFormed.("Control Group") = gnamesF(MaxGranulesFormed.("Control Group"))

%% Plot the data
%plotme = {'Intensity_IntegratedIntensity_GFP', 'Intensity_MeanIntensity_GFP', 'Intensity_StdIntensity_GFP'};
%plotme = {'granspercell','AreaShape_Area'};
%plot_by_ND('cell', dataloc,'plottype',{'meanslope','mean'},'channel',plotme,'looptime',3,'combinexys',true)

% plotme = {'t_granspercell','t_count_cells','NumGrans'};
% %plotme = {};
% plottype = {'mean'}; % 
% plot_by_ND_forJB('treatment', dataloc,'plottype',plottype,'channel',plotme,'looptime',3,'saveloc',dataloc.fold.fig,'subset','24')
%plot_by_ND_forJB('treatment', dataloc,'plottype','sum mean','channel','count_cells','looptime',3,'saveloc',dataloc.fold.fig,'subset','24')

%plot_by_ND_forJB('cell', dataloc,'plottype',plottype,'channel',plotme,'looptime',3,'saveloc',dataloc.fold.fig)


