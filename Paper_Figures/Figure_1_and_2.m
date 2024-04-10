%% Figure 1 and 2 
% Do the analysis and make both of these figures. 

%% Load the model fit data for the live cell
addpath('Z:\code\Nick')
%% 24h induced Wild-Type 4B Data

% 2023-05-03
dataset1 = load('Z:\imageData\SG_4B\2023-05-03 4B-WT 4B-Mut Tet curve NaAsO2 curve\2023-05-03 4B-WT 4B-Mut Tet curve NaAsO2 curve_Processed_Copy.mat');
dataset1 = dataset1.dataloc; % pull the loaded dataloc structure

% append the movie's metadata
dataset1.movieinfo.PixSizeX = 0.16; % um/px
dataset1.movieinfo.PixSizeY = 0.16; % um/px
dataset1.movieinfo.PixNumX = 3200; % pixels
dataset1.movieinfo.PixNumY = 3200; % pixels
dataset1.movieinfo.tsamp = 3; % minutes

% 2023-06-15
dataset2 = load('Z:\imageData\SG_4B\2023-06-15 4b WT-Mut NaAso2 Curve\2023-06-15 4b WT-Mut NaAso2 Curve_Processed_Copy.mat');
dataset2 = dataset2.dataloc; % pull the loaded dataloc structure

% append the movie's metadata
dataset2.movieinfo.PixSizeX = 0.33; % um/px
dataset2.movieinfo.PixSizeY = 0.33; % um/px
dataset2.movieinfo.PixNumX = 1280; % pixels
dataset2.movieinfo.PixNumY = 1080; % pixels
dataset2.movieinfo.tsamp = 3; % minutes

% 2023-06-29
dataset3 = load('Z:\imageData\SG_4B\2023-06-29 4B WT vs Mut TET curve NaAsO2 curve\2023-06-29 4B WT vs Mut TET curve NaAsO2 curve_Processed_Copy.mat');
dataset3 = dataset3.dataloc; % pull the loaded dataloc structure

% append the movie's metadata
dataset3.movieinfo.PixSizeX = 0.33; % um/px
dataset3.movieinfo.PixSizeY = 0.33; % um/px
dataset3.movieinfo.PixNumX = 1600; % pixels
dataset3.movieinfo.PixNumY = 1600; % pixels
dataset3.movieinfo.tsamp = 3; % minutes

% 2023-08-23
%dataset4 = load('');
%dataset4 = dataset5.dataloc; % pull the loaded dataloc structure


%%  24h induced 139A mutant 4B data

% 2023-05-03 (already loaded)
% 2023-06-15 (already loaded)
% 2023-06-29 (already loaded)

%% fit the model to the datasets
[fitData2,~] = convertDatalocToModelFit({dataset3,dataset3,dataset3}, 'NumGrans','pulsepars',{'f','td','ts','rate_in_min','min_to_respond','rsquared','livecell','granarea'});

%% Subset only the good data 
fitData = fitData2; % work with duplicated data (for safety)
gFitData = fitData((fitData.NumGrans_rsquared > 0.85),:); % look for an r squared greater than 0.8 ? 

gFitData.cell = strrep(gFitData.cell,'Hela_eIF4BGFP','HeLa_eIF4BGFP');
gFitData.cell = strrep(gFitData.cell,'Hela_eIF4BF139A','HeLa_4B139AGFP');

%% Now plot the data
% Make the fit models
% datalocDF = makeLiveCellDataframe({dataset1,dataset2,dataset3},'subset','TET100n24t_NaAsO2125u2t');
% 
% plotme = {'NumGrans'}; %,'granspercell'
% plottype = {'albeck mean fit fixed f'}; % 'albeck mean fit'
% 
% plot_by_ND_forJB('treatment', datalocDF,'plottype',plottype,'channel',plotme,'looptime',3,'font_size',8)

%% Focus on the wt 4b cells for now
% Get the Min number of grans FROM MODEL (f), how long the cells were tet induced, and the dose of NaAsO2 for the Wt 4B cells
minGrans = 3;
tetTime = '-24';
naAsO2 = ("62.5"|"125"|"250");

WTsubz = all([(any([contains(gFitData.treatment,['TET at hour ', tetTime]),~contains(gFitData.treatment,'TET')],2)),contains(gFitData.treatment,naAsO2),...
    (gFitData.NumGrans_f >= minGrans), contains(gFitData.cell,'HeLa_eIF4BGFP'),~contains(gFitData.treatment,'15.125uM')],2); % filter for the parameters set above

wtSubData = gFitData(WTsubz,:);
wtSubData.treatment = strrep(wtSubData.treatment,'0.1ug/mL TET at hour -24 and ','');
wtSubData.treatment = strrep(wtSubData.treatment,' at hour 0','');
wtSubData.cell = strrep(wtSubData.cell,'HeLa_eIF4BGFP','Wt');
% Run the statistics comparing the wt-4B treated cells accross various concentrations of NaAsO2 and control

% get averything for the reader if needed
grpstats(wtSubData,"treatment",["mean","median","sem","std"],"DataVars",["NumGrans_rate_in_min","NumGrans_f","NumGrans_min_to_respond"])

% Now just print the means
grpstats(wtSubData,"treatment","mean","DataVars",["NumGrans_rate_in_min","NumGrans_f","NumGrans_min_to_respond"])

%% Start with wt 4b cells treated with 125uM NaAsO2 as control vs all other NaAsO2 concentrations

% Test if number (f) is significantly different btwn the 4B cell lines
[~,~,statsF] = anova1(wtSubData.NumGrans_f,wtSubData.treatment,'off');
[resultsMaxG,~,~,gnamesF] = multcompare(statsF,"CriticalValueType","dunnett",'ControlGroup',find(matches(statsF.gnames,'125uM NaAsO2')),'Display','off','Approximate',false); 
MaxGranulesFormed = array2table(resultsMaxG,"VariableNames", ["Group","Control Group","Lower Limit","Difference","Upper Limit","P-value"]);
MaxGranulesFormed.("Group") = gnamesF(MaxGranulesFormed.("Group"));
MaxGranulesFormed.("Control Group") = gnamesF(MaxGranulesFormed.("Control Group"))

% Test if time to respond is significantly different btwn the 4B cell lines
[~,~,statsT2R] = anova1(wtSubData.NumGrans_min_to_respond,wtSubData.treatment,'off');
[resultsTime2Resp,~,~,gnamesT2R] = multcompare(statsT2R,"CriticalValueType","dunnett",'ControlGroup',find(matches(statsT2R.gnames,'125uM NaAsO2')),'Display','off','Approximate',false); 
Time2Respond = array2table(resultsTime2Resp,"VariableNames", ["Group","Control Group","Lower Limit","Difference","Upper Limit","P-value"]);
Time2Respond.("Group") = gnamesT2R(Time2Respond.("Group"));
Time2Respond.("Control Group") = gnamesT2R(Time2Respond.("Control Group"))

% Test if Rate is significantly different btwn the 125 and other 2 doses
[~,~,statsR] = anova1(wtSubData.NumGrans_rate_in_min,wtSubData.treatment,'off');
[resultsRate,~,~,gnamesRate] = multcompare(statsR,"CriticalValueType","dunnett",'ControlGroup',find(matches(statsR.gnames,'125uM NaAsO2')),'Display','off','Approximate',false); 
RateOfGranuleFormation = array2table(resultsRate,"VariableNames", ["Group","Control Group","Lower Limit","Difference","Upper Limit","P-value"]);
RateOfGranuleFormation.("Group") = gnamesRate(RateOfGranuleFormation.("Group"));
RateOfGranuleFormation.("Control Group") = gnamesRate(RateOfGranuleFormation.("Control Group"))

%% Make Figure 1 - eIF4B-GFP (WT)

% Figure 1D
% plot the max grans, Rate of granule formation, and time to respond to treatment

figure1 = figure;
set(figure1,'Units',"Inches",'Position',[0,0,8.5,11],'PaperPosition',[0,0,8.5,11]);

f1ax = []; clear f1ax;
f1ax = axes;
% plot the rate 
ax(4) = subplot(4,4,13);
boxplot(wtSubData.NumGrans_rate_in_min,wtSubData.treatment,'Symbol','','Notch','on','colorgroup',wtSubData.cell,'Colors',clrz); %
ylim([0,9])
ylabel('Rate of 4B Granule Formation (minutes)')

% plot the max granules 
ax(5) = subplot(4,4,14);
boxplot(wtSubData.NumGrans_f,wtSubData.treatment,'Notch','on','Symbol','','colorgroup',wtSubData.cell,'Colors',clrz);
ylim([0,55])
ylabel('Max Number of 4B Granules (f)')

% plot the time 2 respond
ax(6) = subplot(4,4,15);
boxplot(wtSubData.NumGrans_min_to_respond,wtSubData.treatment,'Notch','on','Symbol','','colorgroup',wtSubData.cell,'Colors',clrz);
ylim([0,55])
ylabel('Time to respond (minutes)')

% plot the average granule area at max f
ax(7) = subplot(4,4,16);
% boxplot(gFitData.NumGrans_f,gFitData.cell,'Notch','on','colorgroup',subData.cell,'Colors',clrz);
% ylim([,])
ylabel('Granule Area (um^2) at f')

% loop through the subplots that make up figure 2D, size and space them appropriately 
pWidth = 1.5; % plot width in inches
pHeight = 3; % plot height in inches
sWidth = (7.25-(pWidth*4))/3; % gap btwn plots in inches
for iSub = 4:7
    set(ax(iSub),'Units','Inches','Position',[0.75+(pWidth*(iSub-4))+(sWidth*(iSub-4)), 0.75, pWidth, pHeight]) 
end

fontsize(8,"points"); fontname("Arial");

saveas(figure1,'Z:\imageData\SG_4B\Paper_Figures\Output_Figures\Figure_1.fig')
saveas(figure1,'Z:\imageData\SG_4B\Paper_Figures\Output_Figures\Figure_1.svg')


%% Figure 2 - RRM of eIF4B suppresses stress granule formation
% Max number of grans FROM MODEL (f) versus treatment and cell line
minGrans = 3;
tetTime = '-24';
naAsO2 = ("62.5"|"125"|"250");

subz = all([(any([contains(gFitData.treatment,['TET at hour ', tetTime]),~contains(gFitData.treatment,'TET')],2)),contains(gFitData.treatment,naAsO2),...
    (gFitData.NumGrans_f >= minGrans),~contains(gFitData.treatment,'15.125uM')],2); % filter for the parameters set above

% subset the data 
subData = gFitData(subz,:);

% simplify the name since all the data is 24hr tet induced, NaAsO2 is at hour 0, and we know HeLa_eIF4BGFP is wt/HeLa_4B139AGFP is mut
subData.treatment = strrep(subData.treatment,'0.1ug/mL TET at hour -24 and ','');
subData.treatment = strrep(subData.treatment,' NaAsO2 at hour 0','');
subData.cell = strrep(subData.cell,'HeLa_4B139AGFP','Mut');
subData.cell = strrep(subData.cell,'HeLa_eIF4BGFP','Wt');

%% Run the statistics comparing the cell lines using wt treated with various concentrations of NaAsO2 as control

% get averything for the reader if needed
grpstats(subData,["treatment","cell"],["mean","median","sem","std"],"DataVars",["NumGrans_rate_in_min","NumGrans_f","NumGrans_min_to_respond","NumGrans_granarea"])

% Now just print the means
grpstats(subData,["treatment","cell"],"mean","DataVars",["NumGrans_rate_in_min","NumGrans_f","NumGrans_min_to_respond","NumGrans_granarea"])

[~,~,statsR] = anova1(subData.NumGrans_granarea,strcat(subData.treatment,{' '},subData.cell),'off');
[resultsRate,~,~,gnamesRate] = multcompare(statsR,"CriticalValueType","dunnett",'ControlGroup',find(matches(statsT2R.gnames,'250uM Wt')),'Display','off','Approximate',false); 
RateOfGranuleFormation = array2table(resultsRate,"VariableNames", ["Group","Control Group","Lower Limit","Difference","Upper Limit","P-value"]);
RateOfGranuleFormation.("Group") = gnamesRate(RateOfGranuleFormation.("Group"));
RateOfGranuleFormation.("Control Group") = gnamesRate(RateOfGranuleFormation.("Control Group"))


%% Start with wt 4b cells treated with 62.5uM NaAsO2 as control vs all other NaAsO2 concentrations and btwn cell lines
% Test if Rate is significantly different btwn the 4B cell lines
[~,~,statsR] = anova1(subData.NumGrans_rate_in_min,strcat(subData.treatment,{' '},subData.cell),'off');
[resultsRate,~,~,gnamesRate] = multcompare(statsR,"CriticalValueType","dunnett",'ControlGroup',find(matches(statsR.gnames,'62.5uM Wt')),'Display','off','Approximate',false); 
RateOfGranuleFormation = array2table(resultsRate,"VariableNames", ["Group","Control Group","Lower Limit","Difference","Upper Limit","P-value"]);
RateOfGranuleFormation.("Group") = gnamesRate(RateOfGranuleFormation.("Group"));
RateOfGranuleFormation.("Control Group") = gnamesRate(RateOfGranuleFormation.("Control Group"))

% Test if number (f) is significantly different btwn the 4B cell lines
[~,~,statsF] = anova1(subData.NumGrans_f,strcat(subData.treatment,{' '},subData.cell),'off');
[resultsMaxG,~,~,gnamesF] = multcompare(statsF,"CriticalValueType","dunnett",'ControlGroup',find(matches(statsF.gnames,'62.5uM Wt')),'Display','off','Approximate',false); 
MaxGranulesFormed = array2table(resultsMaxG,"VariableNames", ["Group","Control Group","Lower Limit","Difference","Upper Limit","P-value"]);
MaxGranulesFormed.("Group") = gnamesF(MaxGranulesFormed.("Group"));
MaxGranulesFormed.("Control Group") = gnamesF(MaxGranulesFormed.("Control Group"))

% Test if time to respond is significantly different btwn the 4B cell lines
[~,~,statsT2R] = anova1(subData.NumGrans_min_to_respond,strcat(subData.treatment,{' '},subData.cell),'off');
[resultsTime2Resp,~,~,gnamesT2R] = multcompare(statsT2R,"CriticalValueType","dunnett",'ControlGroup',find(matches(statsT2R.gnames,'62.5uM Wt')),'Display','off','Approximate',false); 
Time2Respond = array2table(resultsTime2Resp,"VariableNames", ["Group","Control Group","Lower Limit","Difference","Upper Limit","P-value"]);
Time2Respond.("Group") = gnamesT2R(Time2Respond.("Group"));
Time2Respond.("Control Group") = gnamesT2R(Time2Respond.("Control Group"))

%% Now use wt 4b cells treated with 125uM as control vs all other NaAsO2 concentrations and btwn cell lines
% Test if Rate is significantly different btwn the 4B cell lines
[~,~,statsR] = anova1(subData.NumGrans_rate_in_min,strcat(subData.treatment,{' '},subData.cell),'off');
[resultsRate,~,~,gnamesRate] = multcompare(statsR,"CriticalValueType","dunnett",'ControlGroup',find(matches(statsR.gnames,'125uM Wt')),'Display','off','Approximate',false); 
RateOfGranuleFormation = array2table(resultsRate,"VariableNames", ["Group","Control Group","Lower Limit","Difference","Upper Limit","P-value"]);
RateOfGranuleFormation.("Group") = gnamesRate(RateOfGranuleFormation.("Group"));
RateOfGranuleFormation.("Control Group") = gnamesRate(RateOfGranuleFormation.("Control Group"))

% Test if number (f) is significantly different btwn the 4B cell lines
[~,~,statsF] = anova1(subData.NumGrans_f,strcat(subData.treatment,{' '},subData.cell),'off');
[resultsMaxG,~,~,gnamesF] = multcompare(statsF,"CriticalValueType","dunnett",'ControlGroup',find(matches(statsF.gnames,'125uM Wt')),'Display','off','Approximate',false); 
MaxGranulesFormed = array2table(resultsMaxG,"VariableNames", ["Group","Control Group","Lower Limit","Difference","Upper Limit","P-value"]);
MaxGranulesFormed.("Group") = gnamesF(MaxGranulesFormed.("Group"));
MaxGranulesFormed.("Control Group") = gnamesF(MaxGranulesFormed.("Control Group"))

% Test if time to respond is significantly different btwn the 4B cell lines
[~,~,statsT2R] = anova1(subData.NumGrans_min_to_respond,strcat(subData.treatment,{' '},subData.cell),'off');
[resultsTime2Resp,~,~,gnamesT2R] = multcompare(statsT2R,"CriticalValueType","dunnett",'ControlGroup',find(matches(statsT2R.gnames,'125uM Wt')),'Display','off','Approximate',false); 
Time2Respond = array2table(resultsTime2Resp,"VariableNames", ["Group","Control Group","Lower Limit","Difference","Upper Limit","P-value"]);
Time2Respond.("Group") = gnamesT2R(Time2Respond.("Group"));
Time2Respond.("Control Group") = gnamesT2R(Time2Respond.("Control Group"))

%% Now use wt 4b cells treated with 250uM as control vs all other NaAsO2 concentrations and btwn cell lines
% Test if Rate is significantly different btwn the 4B cell lines
[~,~,statsR] = anova1(subData.NumGrans_rate_in_min,strcat(subData.treatment,{' '},subData.cell),'off');
[resultsRate,~,~,gnamesRate] = multcompare(statsR,"CriticalValueType","dunnett",'ControlGroup',find(matches(statsR.gnames,'250uM Wt')),'Display','off','Approximate',false); 
RateOfGranuleFormation = array2table(resultsRate,"VariableNames", ["Group","Control Group","Lower Limit","Difference","Upper Limit","P-value"]);
RateOfGranuleFormation.("Group") = gnamesRate(RateOfGranuleFormation.("Group"));
RateOfGranuleFormation.("Control Group") = gnamesRate(RateOfGranuleFormation.("Control Group"))

% Test if number (f) is significantly different btwn the 4B cell lines
[~,~,statsF] = anova1(subData.NumGrans_f,strcat(subData.treatment,{' '},subData.cell),'off');
[resultsMaxG,~,~,gnamesF] = multcompare(statsF,"CriticalValueType","dunnett",'ControlGroup',find(matches(statsF.gnames,'250uM Wt')),'Display','off','Approximate',false); 
MaxGranulesFormed = array2table(resultsMaxG,"VariableNames", ["Group","Control Group","Lower Limit","Difference","Upper Limit","P-value"]);
MaxGranulesFormed.("Group") = gnamesF(MaxGranulesFormed.("Group"));
MaxGranulesFormed.("Control Group") = gnamesF(MaxGranulesFormed.("Control Group"))

% Test if time to respond is significantly different btwn the 4B cell lines
[~,~,statsT2R] = anova1(subData.NumGrans_min_to_respond,strcat(subData.treatment,{' '},subData.cell),'off');
[resultsTime2Resp,~,~,gnamesT2R] = multcompare(statsT2R,"CriticalValueType","dunnett",'ControlGroup',find(matches(statsT2R.gnames,'250uM Wt')),'Display','off','Approximate',false); 
Time2Respond = array2table(resultsTime2Resp,"VariableNames", ["Group","Control Group","Lower Limit","Difference","Upper Limit","P-value"]);
Time2Respond.("Group") = gnamesT2R(Time2Respond.("Group"));
Time2Respond.("Control Group") = gnamesT2R(Time2Respond.("Control Group"))
%% Plot the data for Figure 2

% make a new figure
figure2 = figure;
set(figure2,'Units',"Inches",'Position',[0,0,8.5,11],'PaperPosition',[0,0,8.5,11]);

ax=[]; clear ax;

% Figure 2A (4B mutant diagram) - make an axis handle for the 4B RRM diagram and load the svg into it
ax(1) = subplot(4,4,[1,4]);
set(ax(1),'Units','Inches','Position',[0.5, 9.2, 7.5, 1]) 
F2A = imread('Z:\imageData\SG_4B\Paper_Figures\Images_For_Figures\Figure_2A_4B_RRM_Diagram.jpg');
imshow(F2A,"Parent",ax)

% Figure 2B - Images of 4B mutant images - [Nuclei, 4B mutant GFP, G3BP1, merged]
ax(2) = subplot(4,4,[5,8]);
set(ax(2),'Units','Inches','Position',[0.5, 7.1, 7.5, 2]) 
xticklabels({''}); yticklabels({''});

% Figure 2C - ???
ax(3) = subplot(4,4,[9,12]);
set(ax(3),'Units','Inches','Position',[0.5, 4, 7.5, 3]) 
xticklabels({''}); yticklabels({''});

% Figure 4D - plot all rate, f, time 2 respond, and max granule area as box and wisker plots

% set the wt vs mutant colors
clrz = [0, 0.5, 0; ...
        0.5, 0, 0];

% plot the rate 
ax(4) = subplot(4,4,13);
boxplot(subData.NumGrans_rate_in_min,{subData.treatment,subData.cell},'Symbol','','Notch','on','colorgroup',subData.cell,'Colors',clrz); %
ylim([0,9])
ylabel('Rate of 4B Granule Formation (minutes)')

% plot the max granules 
ax(5) = subplot(4,4,14);
boxplot(subData.NumGrans_f,{subData.treatment,subData.cell},'Notch','on','Symbol','','colorgroup',subData.cell,'Colors',clrz);
ylim([0,55])
ylabel('Max Number of 4B Granules (f)')

% plot the time 2 respond
ax(6) = subplot(4,4,15);
boxplot(subData.NumGrans_min_to_respond,{subData.treatment,subData.cell},'Notch','on','Symbol','','colorgroup',subData.cell,'Colors',clrz);
ylim([0,55])
ylabel('Time to respond (minutes)')

% plot the average granule area at max f
ax(7) = subplot(4,4,16);
boxplot(subData.NumGrans_granarea,{subData.treatment,subData.cell},'Symbol','','Notch','on','colorgroup',subData.cell,'Colors',clrz);
ylim([1.5,5.5])
ylabel('Granule Area (um^2) at f')

% loop through the subplots that make up figure 2D, size and space them appropriately 
pWidth = 1.5; % plot width in inches
pHeight = 3; % plot height in inches
sWidth = (7.25-(pWidth*4))/3; % gap btwn plots in inches
for iSub = 4:7
    set(ax(iSub),'Units','Inches','Position',[0.75+(pWidth*(iSub-4))+(sWidth*(iSub-4)), 0.75, pWidth, pHeight]) 
end

fontsize(8,"points"); fontname("Arial");
%saveas(figure2,'Z:\imageData\SG_4B\Paper_Figures\Output_Figures\Figure_2.fig')
%saveas(figure2,'Z:\imageData\SG_4B\Paper_Figures\Output_Figures\Figure_2.svg')


