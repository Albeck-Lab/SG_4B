%% Load the model fit data for the live cell
addpath('\\albecklab.mcb.ucdavis.edu\data\code\Nick')
%% 24h induced Wild-Type 4B Data
% 2022-11-03
%dataset1 = load('');
%dataset1 = dataset1.dataloc; % pull the loaded dataloc structure

% 2023-05-03
dataset2 = load('\\albecklab.mcb.ucdavis.edu\Data\imageData\SG_4B\2023-05-03 4B-WT 4B-Mut Tet curve NaAsO2 curve\2023-05-03 4B-WT 4B-Mut Tet curve NaAsO2 curve_Processed.mat');
dataset2 = dataset2.dataloc; % pull the loaded dataloc structure

% 2023-06-15
dataset3 = load('\\albecklab.mcb.ucdavis.edu\Data\imageData\SG_4B\2023-06-15 4b WT-Mut NaAso2 Curve\2023-06-15 4b WT-Mut NaAso2 Curve_Processed.mat');
dataset3 = dataset3.dataloc; % pull the loaded dataloc structure

% 2023-06-29
dataset4 = load('\\albecklab.mcb.ucdavis.edu\Data\imageData\SG_4B\2023-06-29 4B WT vs Mut TET curve NaAsO2 curve\2023-06-29 4B WT vs Mut TET curve NaAsO2 curve_Processed.mat');
dataset4 = dataset4.dataloc; % pull the loaded dataloc structure

% 2023-08-23
%dataset5 = load('');
%dataset5 = dataset5.dataloc; % pull the loaded dataloc structure


%%  24h induced 139A mutant 4B data
% 2023-05-03
dataset6 = load('\\albecklab.mcb.ucdavis.edu\Data\imageData\SG_4B\2023-05-03 4B-WT 4B-Mut Tet curve NaAsO2 curve\2023-05-03 4B-WT 4B-Mut Tet curve NaAsO2 curve_Processed.mat');
dataset6 = dataset6.dataloc; % pull the loaded dataloc structure

% 2023-06-15 (already loaded)
% 2023-06-29 (already loaded)
% 2023-08-23 (already loaded)
%% Make the fit models
datalocDF = makeLiveCellDataframe({dataset2,dataset4,dataset6},'subset','TET100n24t_NaAsO2125u2t');

plotme = {'NumGrans'}; %,'granspercell'
plottype = {'mean'}; % 'albeck mean fit'

plot_by_ND_forJB('treatment', datalocDF,'plottype',plottype,'channel',plotme,'looptime',3,'font_size',8,'tmaxaftertx',1,'saveloc','Z:\imageData\SG_4B\Paper Figures')
plot_by_ND_lcDF('treatment', datalocDF,'plottype',plottype,'channel',plotme,'looptime',3,'font_size',8,'tmaxaftertx',1,'saveloc','Z:\imageData\SG_4B\Paper Figures')
%% fit the model to the dataset
[fitData2,~] = convertDatalocToModelFit({dataset2,dataset4,dataset6}, 'NumGrans');

fitData = fitData2; % work with duplicated data (for safety)

gFitData = fitData((fitData.NumGrans_rsquared > 0.8),:); % look for an r squared greater than 0.8 ? 

% boxplot(gfitData.NumGrans_min_to_respond,gFitData.treatment)

%% Max number of grans FROM MODEL (f) versus treatment and cell line
minGrans = 3;
tetTime = '-24';
ylimz = [0, 50]; % axis limits
naAsO2 = '125';

subz = all([contains(gFitData.treatment,['TET at hour ', tetTime]),contains(gFitData.treatment,[naAsO2,'uM NaAsO2']),...
    (gFitData.NumGrans_f >= minGrans)],2); % filter for the parameters set above

% subset the data 
subData = gFitData(subz,:);

% plot all 4 as box and wisker plots

figgy = figure;

% plot the rate 
subplot(1,4,1);
boxplot(subData.NumGrans_rate_in_min,subData.cell,'Notch','on');
ylim([0,10])
xticklabels({'WT','Mut'})
xlabel('4B Cell line'); ylabel('Rate of 4B Granule Formation (minutes)')

% plot the max granules 
subplot(1,4,2);
boxplot(subData.NumGrans_f,subData.cell,'Notch','on');
ylim([0,52])
xticklabels({'WT','Mut'})
xlabel('4B Cell line'); ylabel('Max Number of 4B Granules (f)')

% plot the time 2 respond
subplot(1,4,3);
boxplot(subData.NumGrans_min_to_respond,subData.cell,'Notch','on');
xticklabels({'WT','Mut'})
ylim([0,50])
xlabel('4B Cell line'); ylabel('Time to respond (minutes)')

% plot the average granule area at max f
% subplot(1,4,1);
% boxplot(gFitData.NumGrans_f,gFitData.cell,'Notch','on');
% ylim(ylimz)
% xlabel('4B Cell line'); ylabel('f (Number of 4B Granules)')


% barHold = grpstats(subData,"cell",["mean","median","sem","std"],"DataVars","NumGrans_f")

% Test if Rate is significantly different btwn the 4B cell lines
[~,~,statsR] = anova1(subData.NumGrans_rate_in_min,subData.cell,'off');
[resultsRate,~,~,gnamesRate] = multcompare(statsR,"CriticalValueType","dunnett",'ControlGroup',find(matches(statsR.gnames,'HeLa_eIF4BGFP')),'Display','off'); 
RateOfGranuleFormation = array2table(resultsRate,"VariableNames", ["Group","Control Group","Lower Limit","Difference","Upper Limit","P-value"]);
RateOfGranuleFormation.("Group") = gnamesRate(RateOfGranuleFormation.("Group"));
RateOfGranuleFormation.("Control Group") = gnamesRate(RateOfGranuleFormation.("Control Group"))

% Test if number (f) is significantly different btwn the 4B cell lines
[~,~,statsF] = anova1(subData.NumGrans_f,subData.cell,'off');
[resultsMaxG,~,~,gnamesF] = multcompare(statsF,"CriticalValueType","dunnett",'ControlGroup',find(matches(statsF.gnames,'HeLa_eIF4BGFP')),'Display','off'); 
MaxGranulesFormed = array2table(resultsMaxG,"VariableNames", ["Group","Control Group","Lower Limit","Difference","Upper Limit","P-value"]);
MaxGranulesFormed.("Group") = gnamesF(MaxGranulesFormed.("Group"));
MaxGranulesFormed.("Control Group") = gnamesF(MaxGranulesFormed.("Control Group"))

% Test if time to respond is significantly different btwn the 4B cell lines
[~,~,statsT2R] = anova1(subData.NumGrans_min_to_respond,subData.cell,'off');
[resultsTime2Resp,~,~,gnamesT2R] = multcompare(statsT2R,"CriticalValueType","dunnett",'ControlGroup',find(matches(statsT2R.gnames,'HeLa_eIF4BGFP')),'Display','off'); 
Time2Respond = array2table(resultsTime2Resp,"VariableNames", ["Group","Control Group","Lower Limit","Difference","Upper Limit","P-value"]);
Time2Respond.("Group") = gnamesT2R(Time2Respond.("Group"));
Time2Respond.("Control Group") = gnamesT2R(Time2Respond.("Control Group"))

%% pull from combined data figures
% Test if average size is significantly different btwn the 4B cell lines
[~,~,stats] = anova1(subData.NumGrans_f,subData.cell,'off');
[resultsTime2Resp,~,~,gnamesT2R] = multcompare(stats,"CriticalValueType","dunnett",'ControlGroup',find(matches(stats.gnames,'HeLa_eIF4BGFP')),'Display','off'); 
Time2Respond = array2table(resultsTime2Resp,"VariableNames", ["Group","Control Group","Lower Limit","Difference","Upper Limit","P-value"]);
Time2Respond.("Group") = gnamesT2R(Time2Respond.("Group"));
Time2Respond.("Control Group") = gnamesT2R(Time2Respond.("Control Group"))


% title([tetTime,'h WT4b Tet induced and all NaAsO2 Doses'])
% 
% subz = all([contains(gFitData.treatment,['TET at hour ', tetTime]),contains(gFitData.cell,'139A'),...
%     (gFitData.NumGrans_f >= minGrans)],2);
% figure;
% boxplot(gFitData.NumGrans_f(subz),gFitData.treatment(subz));
% xlabel('treatment'); ylabel('f (Number of 4B Granules)')
% ylim(ylimz)
% title([tetTime,'h 139A mut Tet induced and all NaAsO2 Doses'])


%% Max number of grans FROM MODEL (f) versus treatment and cell line
minGrans = 3;
tetTime = '-24';
ylimz = [0, 20]; % Opp axis limits
% 
% figure
% hp1 = uipanel('position',[0,   0, 0.5, 1]);
% hp2 = uipanel('position',[0.5, 0, 0.5, 1]);

figure;
subz = all([contains(gFitData.treatment,['TET at hour ', tetTime]),contains(gFitData.cell,'4BGFP'),...
    (gFitData.NumGrans_f >= minGrans)],2);
boxplot(gFitData.NumGrans_f(subz),gFitData.treatment(subz));
ylim(ylimz)
xlabel('treatment'); ylabel('f (Number of 4B Granules)')
title([tetTime,'h WT4b Tet induced and all NaAsO2 Doses'])

subz = all([contains(gFitData.treatment,['TET at hour ', tetTime]),contains(gFitData.cell,'139A'),...
    (gFitData.NumGrans_f >= minGrans)],2);
figure;
boxplot(gFitData.NumGrans_f(subz),gFitData.treatment(subz));
xlabel('treatment'); ylabel('f (Number of 4B Granules)')
ylim(ylimz)
title([tetTime,'h 139A mut Tet induced and all NaAsO2 Doses'])















%% 24 hour tet wt 4b grans or 139a vs opp vs naaso2 dose/vehicle
minGrans = 0;

figure
hp1 = uipanel('position',[0,   0, 0.5, 1]);
hp2 = uipanel('position',[0.5, 0, 0.5, 1]);

subz = all([contains(dataloc.ifd.treatment,'TET at hour -24'),contains(dataloc.ifd.cell,'4BGFP'),...
    (dataloc.ifd.Children_Grans_4B_Count >= minGrans),],2);
scatterhist(dataloc.ifd.Children_Grans_4B_Count(subz),dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subz)*65000,...
    'Group',dataloc.ifd.treatment(subz),'Kernel','on','Direction','Out','Parent',hp1,'Location','SouthEast' )
xlim([0,15]); ylim([0,10000]); xlabel('4B Granules'); ylabel('Mean OPP Intensity')
title('24h WT4b Tet induced and all NaAsO2 Doses')


subz = all([contains(dataloc.ifd.treatment,'TET at hour -24'),contains(dataloc.ifd.cell,'139A'),...
    (dataloc.ifd.Children_Grans_4B_Count >= minGrans)],2);
scatterhist(dataloc.ifd.Children_Grans_4B_Count(subz),dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subz)*65000,...
    'Group',dataloc.ifd.treatment(subz),'Kernel','on','Direction','Out','Parent',hp2)
xlim([0,15]); ylim([0,10000]); xlabel('4B Granules'); ylabel('Mean OPP Intensity')
title('24h 137A Tet induced and all NaAsO2 Doses')

%% X hour tet wt 4b grans or 139a vs opp vs naaso2 dose/vehicle with box and wiskers?
minGrans = 0;
tetTime = '-24';
xlimz = [0,30]; % granule axis limits
ylimz = [0, 10000]; % Opp axis limits

figure
hp1 = uipanel('position',[0,   0, 0.5, 1]);
hp2 = uipanel('position',[0.5, 0, 0.5, 1]);

subz = all([contains(dataloc.ifd.treatment,['TET at hour ', tetTime]),contains(dataloc.ifd.cell,'4BGFP'),...
    (dataloc.ifd.Children_Grans_4B_Count >= minGrans),],2);
h = scatterhist(dataloc.ifd.Children_Grans_4B_Count(subz),dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subz)*65000,...
    'Group',dataloc.ifd.treatment(subz),'Parent',hp1,'Location','SouthEast');
xlabel('4B Granules'); ylabel('Mean OPP Intensity')
title([tetTime,'h WT4b Tet induced and all NaAsO2 Doses'])
hold on;
clr = get(h(1),'colororder');
boxplot(h(2),dataloc.ifd.Children_Grans_4B_Count(subz),dataloc.ifd.treatment(subz),'orientation','horizontal',...
     'label',{'','','',''},'color',clr);
boxplot(h(3),(dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subz)*65000),dataloc.ifd.treatment(subz),'orientation','horizontal',...
     'label', {'','','',''},'color',clr);
set(h(2:3),'XTickLabel','');
view(h(3),[270,90]);  % Rotate the Y plot
xlim(h(1),xlimz); ylim(h(1),ylimz);
xlim(h(2),xlimz); xlim(h(3),ylimz); % Sync axes
hold off;


subz = all([contains(dataloc.ifd.treatment,['TET at hour ', tetTime]),contains(dataloc.ifd.cell,'139A'),...
    (dataloc.ifd.Children_Grans_4B_Count >= minGrans)],2);
h=scatterhist(dataloc.ifd.Children_Grans_4B_Count(subz),dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subz)*65000,...
    'Group',dataloc.ifd.treatment(subz),'Parent',hp2);
xlabel('4B Granules'); ylabel('Mean OPP Intensity')
title([tetTime,'h 137A Tet induced and all NaAsO2 Doses'])
hold on;
clr = get(h(1),'colororder');
boxplot(h(2),dataloc.ifd.Children_Grans_4B_Count(subz),dataloc.ifd.treatment(subz),'orientation','horizontal',...
     'label',{'','','',''},'color',clr);
boxplot(h(3),dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subz)*65000,dataloc.ifd.treatment(subz),'orientation','horizontal',...
     'label', {'','','',''},'color',clr);
set(h(2:3),'XTickLabel','');
view(h(3),[270,90]);  % Rotate the Y plot
xlim(h(1),xlimz); ylim(h(1),ylimz);
xlim(h(2),xlimz); xlim(h(3),ylimz);  % Sync axes
hold off;

%% X hour tet wt 4b grans or 139a vs opp vs naaso2 dose/vehicle with box and wiskers > 0 uM
minGrans = 2;
tetTime = '-24';
xlimz = [0,35]; % granule axis limits
ylimz = [0, 5000]; % Opp axis limits

figure
hp1 = uipanel('position',[0,   0, 0.5, 1]);
hp2 = uipanel('position',[0.5, 0, 0.5, 1]);

subz = all([contains(dataloc.ifd.treatment,['TET at hour ', tetTime]),contains(dataloc.ifd.cell,'4BGFP'),...
    ~contains(dataloc.ifd.treatment,' 0uM'),(dataloc.ifd.Children_Grans_4B_Count >= minGrans),],2);
h = scatterhist((dataloc.ifd.Children_Grans_4B_Count(subz)./(dataloc.ifd.AreaShape_Area(subz)*0.1089/1000)),((dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subz)*65000)./(dataloc.ifd.AreaShape_Area(subz)*0.1089/1000)),...
    'Group',dataloc.ifd.treatment(subz),'Parent',hp1,'Location','SouthEast');
xlabel('4B Granules'); ylabel('Mean OPP Intensity')
title([tetTime,'h WT4b Tet induced and all NaAsO2 Doses'])
hold on;
clr = get(h(1),'colororder');
boxplot(h(2),(dataloc.ifd.Children_Grans_4B_Count(subz)./(dataloc.ifd.AreaShape_Area(subz)*0.1089/1000)),dataloc.ifd.treatment(subz),'orientation','horizontal',...
     'label',{'','',''},'color',clr,'Notch','on');
boxplot(h(3),((dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subz)*65000)./(dataloc.ifd.AreaShape_Area(subz)*0.1089/1000)),dataloc.ifd.treatment(subz),'orientation','horizontal',...
     'label', {'','',''},'color',clr,'Notch','on');
set(h(2:3),'XTickLabel','');
view(h(3),[270,90]);  % Rotate the Y plot
xlim(h(1),xlimz); ylim(h(1),ylimz);
xlim(h(2),xlimz); xlim(h(3),ylimz); % Sync axes
hold off;


subz = all([contains(dataloc.ifd.treatment,['TET at hour ', tetTime]),contains(dataloc.ifd.cell,'139A'),...
    ~contains(dataloc.ifd.treatment,' 0uM'),(dataloc.ifd.Children_Grans_4B_Count >= minGrans)],2);
h=scatterhist((dataloc.ifd.Children_Grans_4B_Count(subz)./(dataloc.ifd.AreaShape_Area(subz)*0.1089/1000)),((dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subz)*65000)./(dataloc.ifd.AreaShape_Area(subz)*0.1089/1000)),...
    'Group',dataloc.ifd.treatment(subz),'Parent',hp2);
xlabel('4B Granules'); ylabel('Mean OPP Intensity')
title([tetTime,'h 137A Tet induced and all NaAsO2 Doses'])
hold on;
clr = get(h(1),'colororder');
boxplot(h(2),(dataloc.ifd.Children_Grans_4B_Count(subz)./(dataloc.ifd.AreaShape_Area(subz)*0.1089/1000)),dataloc.ifd.treatment(subz),'orientation','horizontal',...
     'label',{'','',''},'color',clr,'Notch','on');
boxplot(h(3),((dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subz)*65000)./(dataloc.ifd.AreaShape_Area(subz)*0.1089/1000)),dataloc.ifd.treatment(subz),'orientation','horizontal',...
     'label', {'','',''},'color',clr,'Notch','on');
set(h(2:3),'XTickLabel','');
view(h(3),[270,90]);  % Rotate the Y plot
xlim(h(1),xlimz); ylim(h(1),ylimz);
xlim(h(2),xlimz); xlim(h(3),ylimz);  % Sync axes
hold off;
%% X hour tet wt 4b grans or 139a vs opp vs naaso2 dose > 0
minGrans = 1;
tetTime = '-24';
figure
hp1 = uipanel('position',[0,   0, 0.5, 1]);
hp2 = uipanel('position',[0.5, 0, 0.5, 1]);

subz = all([contains(dataloc.ifd.treatment,['TET at hour ', tetTime]),contains(dataloc.ifd.cell,'4BGFP'),~contains(dataloc.ifd.treatment,' 0uM'),...
    (dataloc.ifd.Children_Grans_4B_Count >= minGrans),],2);
scatterhist((dataloc.ifd.Children_Grans_4B_Count(subz)./(dataloc.ifd.AreaShape_Area(subz)/1000)),dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subz)*65000,...
    'Group',dataloc.ifd.treatment(subz),'Kernel','on','Direction','Out','Parent',hp1,'Location','SouthEast' )
xlim([0,5]); ylim([0,3500]); xlabel('4B Granules'); ylabel('Mean OPP Intensity')
title([tetTime,'h WT4b Tet induced and all NaAsO2 Doses'])


subz = all([contains(dataloc.ifd.treatment,['TET at hour ', tetTime]),contains(dataloc.ifd.cell,'139A'),~contains(dataloc.ifd.treatment,' 0uM'),...
    (dataloc.ifd.Children_Grans_4B_Count >= minGrans)],2);
scatterhist((dataloc.ifd.Children_Grans_4B_Count(subz)./(dataloc.ifd.AreaShape_Area(subz)/1000)),dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subz)*65000,...
    'Group',dataloc.ifd.treatment(subz),'Kernel','on','Direction','Out','Parent',hp2)
xlim([0,5]); ylim([0,3500]); xlabel('4B Granules'); ylabel('Mean OPP Intensity')
title([tetTime,'h 137A Tet induced and all NaAsO2 Doses'])

%% X hour tet wt 4b grans or 139a vs opp vs naaso2 dose > 0
minGrans = 0;
tetTime = '-24';
figure
hp1 = uipanel('position',[0,   0, 0.5, 1]);
hp2 = uipanel('position',[0.5, 0, 0.5, 1]);

subz = all([contains(dataloc.ifd.treatment,['TET at hour ', tetTime]),contains(dataloc.ifd.cell,'4BGFP'),~contains(dataloc.ifd.treatment,' 0uM'),...
    (dataloc.ifd.Children_Grans_4B_Count >= minGrans),],2);
scatterhist(dataloc.ifd.Children_Grans_4B_Count(subz),dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subz)*65000,...
    'Group',dataloc.ifd.treatment(subz),'Kernel','on','Direction','Out','Parent',hp1,'Location','SouthEast' )
xlim([0,20]); ylim([0,4000]); xlabel('4B Granules'); ylabel('Mean OPP Intensity')
title([tetTime,'h WT4b Tet induced and all NaAsO2 Doses'])


subz = all([contains(dataloc.ifd.treatment,['TET at hour ', tetTime]),contains(dataloc.ifd.cell,'139A'),~contains(dataloc.ifd.treatment,' 0uM'),...
    (dataloc.ifd.Children_Grans_4B_Count >= minGrans)],2);
scatterhist(dataloc.ifd.Children_Grans_4B_Count(subz),dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subz)*65000,...
    'Group',dataloc.ifd.treatment(subz),'Kernel','on','Direction','Out','Parent',hp2)
xlim([0,20]); ylim([0,4000]); xlabel('4B Granules'); ylabel('Mean OPP Intensity')
title([tetTime,'h 137A Tet induced and all NaAsO2 Doses'])

%% X hour tet wt 4b grans or 139a vs opp @ NO naaso2 dose
naDose = ' 0uM';
figure;
subz = all([contains(dataloc.ifd.cell,'4BGFP'),contains(dataloc.ifd.treatment,naDose)],2); % subset the data to be wt 4b
catz = categorical(dataloc.ifd.treatment(subz));
catz = reordercats(catz,unique((dataloc.ifd.treatment(subz)),"stable"));
boxplot(dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subz)*65000,catz,'Notch','on')
ylim([0, 22000])
xlabel('TET Induction'); ylabel('Mean OPP Intensity')
title('WT4b Tet induction versus OPP')

figure;
subz = all([contains(dataloc.ifd.cell,'139A'),contains(dataloc.ifd.treatment,naDose)],2); % subset the data to be 139A 4b
catz = categorical(dataloc.ifd.treatment(subz));
catz = reordercats(catz,unique((dataloc.ifd.treatment(subz)),"stable"));
boxplot(dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subz)*65000,catz,'Notch','on')
xlabel('TET Induction'); ylabel('Mean OPP Intensity')
ylim([0, 22000])
title('139A Tet induction versus OPP')
%% OPP vs G3BP1 Grans by TET induction
tetDose = '250uM';
minGrans = 3;
ylimz = [0, 4000];
xlimz = [0, 20];

figure
hp1 = uipanel('position',[0,   0, 0.5, 1]);
hp2 = uipanel('position',[0.5, 0, 0.5, 1]);

subz = all([contains(dataloc.ifd.treatment,tetDose),contains(dataloc.ifd.cell,'4BGFP'), ...
    (dataloc.ifd.Children_Grans_G3BP1_Count >= minGrans)],2);
h = scatterhist(dataloc.ifd.Children_Grans_G3BP1_Count(subz),dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subz)*65000,...
    'Group',dataloc.ifd.full(subz),'Parent',hp1,'Location','SouthEast'); %'Color','kcbgmr',
xlabel('G3BP1 Granules');ylabel('Mean OPP Intensity')
title([tetDose,' NaAsO2 Treated OPP vs Grans by WT 4B by TET induction'])

hold on;
clr = get(h(1),'colororder');
boxplot(h(2),dataloc.ifd.Children_Grans_G3BP1_Count(subz),dataloc.ifd.treatment(subz),'orientation','horizontal',...
     'label',{'','','','','',''},'color',clr);

boxplot(h(3),(dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subz)*65000),dataloc.ifd.treatment(subz),'orientation','horizontal',...
     'label', {'','','','','',''},'color',clr);
set(h(2:3),'XTickLabel','');
view(h(3),[270,90]);  % Rotate the Y plot
hold off;
xlim(h(1),xlimz); ylim(h(1),ylimz);
xlim(h(2),xlimz); xlim(h(3),ylimz);

subz = all([contains(dataloc.ifd.treatment,tetDose),contains(dataloc.ifd.cell,'139A'),...
    (dataloc.ifd.Children_Grans_G3BP1_Count >= minGrans)],2);
h = scatterhist(dataloc.ifd.Children_Grans_G3BP1_Count(subz),dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subz)*65000,...
    'Group',dataloc.ifd.full(subz),'Parent',hp2); %'Color','kcbgmr',
xlabel('G3BP1 Granules');ylabel('Mean OPP Intensity')
title([tetDose,' NaAsO2 Treated OPP vs Grans by 139A Mut 4B by TET induction'])
hold on;
clr = get(h(1),'colororder');
boxplot(h(2),dataloc.ifd.Children_Grans_G3BP1_Count(subz),dataloc.ifd.treatment(subz),'orientation','horizontal',...
     'label',{'','','','','',''},'color',clr);
boxplot(h(3),dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subz)*65000,dataloc.ifd.treatment(subz),'orientation','horizontal',...
     'label', {'','','','','',''},'color',clr);
set(h(2:3),'XTickLabel','');
view(h(3),[270,90]);  % Rotate the Y plot
xlim(h(1),xlimz); ylim(h(1),ylimz);
xlim(h(2),xlimz); xlim(h(3),ylimz);
hold off;

%% 4B wt or mut grans vs OPP grouped by TET induction without 0
figure;
subz = all([contains(dataloc.ifd.cell,'4BGFP'),~contains(dataloc.ifd.treatment,' 0uM')],2); % subset the data to be wt 4b
catz = categorical(dataloc.ifd.treatment(subz));
catz = reordercats(catz,unique((dataloc.ifd.treatment(subz)),"stable"));
boxplot(dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subz)*65000,catz,"ColorGroup",catz,'PlotStyle','compact')
ylim([0, 4000])
xlabel('TET Induction'); ylabel('Mean OPP Intensity')
title('WT4b Tet induction versus OPP')

figure;
subz = all([contains(dataloc.ifd.cell,'139A'),~contains(dataloc.ifd.treatment,' 0uM')],2); % subset the data to be 139A 4b
catz = categorical(dataloc.ifd.treatment(subz));
catz = reordercats(catz,unique((dataloc.ifd.treatment(subz)),"stable"));
boxplot(dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subz)*65000,catz,"ColorGroup",catz,'PlotStyle','compact')
xlabel('TET Induction'); ylabel('Mean OPP Intensity')
ylim([0, 4000])
title('139A Tet induction versus OPP')


%% 4B wt or mut grans vs OPP grouped by TET induction (only 24hr and 0 hr induced)
minGrans = 0;
tetTime = ("hour -24");%|"hour 0"
naAsO2 = (" 125uM"|" 0uM");
% hp1 = uipanel('position',[0,   0, 0.5, 1]);
% hp2 = uipanel('position',[0.5, 0, 0.5, 1]);
% dataloc.ifd.ptx = extractBefore(dataloc.ifd.treatment,' and ');
% dataloc.ifd.tx1 = extractAfter(dataloc.ifd.treatment,' and ');


subz = all([contains(dataloc.ifd.treatment,tetTime),~contains(dataloc.ifd.cell,'bad'),contains(dataloc.ifd.treatment,naAsO2),...
    (dataloc.ifd.Children_Grans_G3BP1_Count >= minGrans)],2);

catOrder = ["0ug/mL TET at hour 0 and 0uM NaAsO2","0.1ug/mL TET at hour -24 and 0uM NaAsO2",...
    "0.1ug/mL TET at hour -24 and 125uM NaAsO2"]; % "0ug/mL TET at hour 0 and 250uM NaAsO2",
statsOut = grpstats(dataloc.ifd(subz,:),["cell","treatment"],["mean","median","sem","std"],"DataVars",["Intensity_MeanIntensity_Masked_OPP","Children_Grans_4B_Count","Children_Grans_G3BP1_Count"]);

figure
boxchart(categorical(dataloc.ifd.cell(subz),unique(dataloc.ifd.cell(subz),"stable")),log2((dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subz)*65535)),'GroupByColor',categorical(dataloc.ifd.treatment(subz),catOrder),'Notch','on','MarkerStyle','.',"BoxWidth",1)
ylim([7,14.5])
legend('Location','best'); ylabel("Log_2 Mean OPP Intensity per Cell")


figure
catOrder2={'HeLa eIF4BGFP_0.1ug/mL TET at hour -24 and 0uM NaAsO2','HeLa eIF4BGFP_0.1ug/mL TET at hour -24 and 125uM NaAsO2','HeLa 4B139AGFP_0.1ug/mL TET at hour -24 and 125uM NaAsO2'};
bar(categorical(statsOut.Properties.RowNames,catOrder2),statsOut.mean_Intensity_MeanIntensity_Masked_OPP*65525)
hold on;
errorbar(categorical(statsOut.Properties.RowNames,catOrder2),(statsOut.mean_Intensity_MeanIntensity_Masked_OPP*65525),(statsOut.sem_Intensity_MeanIntensity_Masked_OPP*65525),"LineStyle","none","LineWidth",2,'Color','k')
ylabel("Mean OPP Intensity per Cell")

[~,~,stats] = anova1((dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subz)*65535),dataloc.ifd.full(subz,:),'off');
[results,~,~,gnames] = multcompare(stats,"CriticalValueType","dunnett",'ControlGroup',find(contains(stats.gnames,'HeLa eIF4BGFP20k TET100n24t NaAsO2125u2t')),'Display','off'); 
granstats = array2table(results,"VariableNames", ["Group","Control Group","Lower Limit","Difference","Upper Limit","P-value"]);
granstats.("Group") = gnames(granstats.("Group"));
granstats.("Control Group") = gnames(granstats.("Control Group"));

figure
catOrder2={'HeLa eIF4BGFP_0.1ug/mL TET at hour -24 and 0uM NaAsO2','HeLa eIF4BGFP_0.1ug/mL TET at hour -24 and 125uM NaAsO2','HeLa 4B139AGFP_0.1ug/mL TET at hour -24 and 125uM NaAsO2'};
bar(categorical(statsOut.Properties.RowNames,catOrder2),statsOut.mean_Intensity_MeanIntensity_Masked_OPP*65525)
hold on;
errorbar(categorical(statsOut.Properties.RowNames,catOrder2),(statsOut.mean_Intensity_MeanIntensity_Masked_OPP*65525),(statsOut.sem_Intensity_MeanIntensity_Masked_OPP*65525),"LineStyle","none","LineWidth",2,'Color','k')
ylabel("Mean OPP Intensity per Cell")

[~,~,stats] = anova1((dataloc.ifd.Children_Grans_4B_Count(subz)),dataloc.ifd.full(subz,:),'off');
[results,~,~,gnames] = multcompare(stats,"CriticalValueType","dunnett",'ControlGroup',find(contains(stats.gnames,'HeLa eIF4BGFP20k TET100n24t NaAsO2125u2t')),'Display','off'); 
granstats = array2table(results,"VariableNames", ["Group","Control Group","Lower Limit","Difference","Upper Limit","P-value"]);
granstats.("Group") = gnames(granstats.("Group"));
granstats.("Control Group") = gnames(granstats.("Control Group"));

%% 4B wt or mut cells G3BP1 grans treated with 125uM NaAsO2 (only 24hr induced)
minGrans = 1;
tetTime = ("hour -24");
naAsO2 = (" 125uM");

subz = all([contains(dataloc.ifd.treatment,tetTime),~contains(dataloc.ifd.cell,'bad'),contains(dataloc.ifd.treatment,naAsO2),...
    (dataloc.ifd.Children_Grans_4B_Count >= minGrans)],2);

fourBstats = grpstats(dataloc.ifd(subz,:),"cell",["mean","median","sem","std","meanci"],"DataVars","Children_Grans_4B_Count");

figure;
plt = bar(categorical(fourBstats.cell,["HeLa eIF4BGFP","HeLa 4B139AGFP"]),fourBstats.mean_Children_Grans_4B_Count,1,"LineWidth",1);
hold on;
errorbar(categorical(fourBstats.cell,["HeLa eIF4BGFP","HeLa 4B139AGFP"]),fourBstats.mean_Children_Grans_4B_Count,fourBstats.sem_Children_Grans_4B_Count,"LineStyle","none","Color","k","LineWidth",1)

plt.CData(2,:) = [1, 0, 0];

anova1(dataloc.ifd.Children_Grans_4B_Count(subz,:),dataloc.ifd.cell(subz,:))


ylim([7,14.5])
ylabel("Log_2 Mean OPP Intensity per Cell")



%% EXTRA
statsOut = grpstats(dataloc.ifd(subz,:),["cell","ptx","tx1"],["mean","median","sem","std"],"DataVars",["Intensity_MeanIntensity_Masked_OPP","Children_Grans_4B_Count","Children_Grans_G3BP1_Count"]);
scatterhist(categorical(dataloc.ifd.cell(subz),unique(dataloc.ifd.cell(subz),"stable")),(dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subz)*65535),"Group",dataloc.ifd.treatment(subz))
boxchart(categorical(dataloc.ifd.cell(subz),unique(dataloc.ifd.cell(subz),"stable")),log2((dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subz)*65535)),'GroupByColor',dataloc.ifd.treatment(subz),'Notch','on','MarkerStyle','.')
legend('Location','best'); 



bar(categorical(statsOut.Properties.RowNames,unique(statsOut.Properties.RowNames,"stable")),statsOut.mean_Intensity_MeanIntensity_Masked_OPP*65525)

scatterhist(dataloc.ifd.Children_Grans_4B_Count(subz),dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subz)*65535,...
    'Group',dataloc.ifd.treatment(subz),'Kernel','on','Direction','Out','Parent',hp1,'Location','SouthEast' )
xlim([0,20]); ylim([0,4000]); xlabel('4B Granules'); ylabel('Mean OPP Intensity')
title([tetTime,'h WT4b Tet induced and all NaAsO2 Doses'])


subz = all([contains(dataloc.ifd.treatment,['TET at hour ', tetTime]),contains(dataloc.ifd.cell,'139A'),~contains(dataloc.ifd.treatment,' 0uM'),...
    (dataloc.ifd.Children_Grans_4B_Count >= minGrans)],2);
scatterhist(dataloc.ifd.Children_Grans_4B_Count(subz),dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subz)*65535,...
    'Group',dataloc.ifd.treatment(subz),'Kernel','on','Direction','Out','Parent',hp2)
xlim([0,20]); ylim([0,4000]); xlabel('4B Granules'); ylabel('Mean OPP Intensity')
title([tetTime,'h 137A Tet induced and all NaAsO2 Doses'])
