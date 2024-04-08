%% Figure 1: 

%% Load the model fit data for the live cell
addpath('Z:\code\Nick')
%% 24h induced Wild-Type 4B Data

% 2023-05-03
dataset1 = load('Z:\imageData\SG_4B\2023-05-03 4B-WT 4B-Mut Tet curve NaAsO2 curve\2023-05-03 4B-WT 4B-Mut Tet curve NaAsO2 curve_Processed_Copy.mat');
dataset1 = dataset1.dataloc; % pull the loaded dataloc structure

% 2023-06-15
dataset2 = load('Z:\imageData\SG_4B\2023-06-15 4b WT-Mut NaAso2 Curve\2023-06-15 4b WT-Mut NaAso2 Curve_Processed_Copy.mat');
dataset2 = dataset2.dataloc; % pull the loaded dataloc structure

% 2023-06-29
dataset3 = load('Z:\imageData\SG_4B\2023-06-29 4B WT vs Mut TET curve NaAsO2 curve\2023-06-29 4B WT vs Mut TET curve NaAsO2 curve_Processed_Copy.mat');
dataset3 = dataset3.dataloc; % pull the loaded dataloc structure

% 2023-08-23
%dataset4 = load('');
%dataset4 = dataset5.dataloc; % pull the loaded dataloc structure

%% Make the fit models
% datalocDF = makeLiveCellDataframe({dataset1,dataset2,dataset3},'subset','TET100n24t_NaAsO2125u2t');
% 
% plotme = {'NumGrans'}; %,'granspercell'
% plottype = {'albeck mean fit fixed f'}; % 'albeck mean fit'
% 
% plot_by_ND_forJB('treatment', datalocDF,'plottype',plottype,'channel',plotme,'looptime',3,'font_size',8)

%% fit the model to the datasets
[fitData2,~] = convertDatalocToModelFit({dataset1,dataset2,dataset3}, 'NumGrans');

%% Subset only the WT-4B-GFP data with a model fit r^2 above 0.8

fitData = fitData2; % work with duplicated data (for safety)
gFitData = fitData((fitData.NumGrans_rsquared > 0.8),:); % look for an r squared greater than 0.8 ? 

% Max number of grans FROM MODEL (f) versus treatment and cell line
minGrans = 3;
tetTime = '-24';
naAsO2 = ("62.5"|"125"|"250");

subz = all([contains(gFitData.treatment,['TET at hour ', tetTime]),contains(gFitData.treatment,naAsO2),...
    (gFitData.NumGrans_f >= minGrans), contains(gFitData.cell,'HeLa_eIF4BGFP')],2); % filter for the parameters set above

% subset the data 
subData = gFitData(subz,:);

% simplify the name since all the data is 24hr tet induced, NaAsO2 is at hour 0, and we know HeLa_eIF4BGFP is wt
subData.treatment = strrep(subData.treatment,'0.1ug/mL TET at hour -24 and ','');
subData.treatment = strrep(subData.treatment,' NaAsO2 at hour 0','');
subData.cell = strrep(subData.cell,'HeLa_eIF4BGFP','Wt');

%% Run the statistics comparing the cell lines using wt treated with various concentrations of NaAsO2 as control

% get averything for the reader if needed
grpstats(subData,"treatment",["mean","median","sem","std"],"DataVars",["NumGrans_rate_in_min","NumGrans_f","NumGrans_min_to_respond"])

% Now just print the means
grpstats(subData,"treatment","mean","DataVars",["NumGrans_rate_in_min","NumGrans_f","NumGrans_min_to_respond"])

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
ylim([0,60])
ylabel('Max Number of 4B Granules (f)')

% plot the time 2 respond
ax(6) = subplot(4,4,15);
boxplot(subData.NumGrans_min_to_respond,{subData.treatment,subData.cell},'Notch','on','Symbol','','colorgroup',subData.cell,'Colors',clrz);
ylim([0,60])
ylabel('Time to respond (minutes)')

% plot the average granule area at max f
ax(7) = subplot(4,4,16);
% boxplot(gFitData.NumGrans_f,gFitData.cell,'Notch','on','colorgroup',subData.cell,'Colors',clrz);
% ylim([,])
ylabel('Max 4B Granule Size (pixels)')

% loop through the subplots that make up figure 2D, size and space them appropriately 
pWidth = 1.5; % plot width in inches
pHeight = 3; % plot height in inches
sWidth = (7.25-(pWidth*4))/3; % gap btwn plots in inches
for iSub = 4:7
    set(ax(iSub),'Units','Inches','Position',[0.75+(pWidth*(iSub-4))+(sWidth*(iSub-4)), 0.75, pWidth, pHeight]) 
end

fontsize(8,"points"); fontname("Arial");
saveas(figure2,'Z:\imageData\SG_4B\Paper_Figures\Output_Figures\Figure_2.fig')
saveas(figure2,'Z:\imageData\SG_4B\Paper_Figures\Output_Figures\Figure_2.svg')

%% pull from combined data figures
% Test if average size is significantly different btwn the 4B cell lines
% [~,~,stats] = anova1(subData.NumGrans_f,subData.cell,'off');
% [resultsTime2Resp,~,~,gnamesT2R] = multcompare(stats,"CriticalValueType","dunnett",'ControlGroup',find(matches(stats.gnames,'HeLa_eIF4BGFP')),'Display','off'); 
% Time2Respond = array2table(resultsTime2Resp,"VariableNames", ["Group","Control Group","Lower Limit","Difference","Upper Limit","P-value"]);
% Time2Respond.("Group") = gnamesT2R(Time2Respond.("Group"));
% Time2Respond.("Control Group") = gnamesT2R(Time2Respond.("Control Group"))


% title([tetTime,'h WT4b Tet induced and all NaAsO2 Doses'])
% 
% subz = all([contains(gFitData.treatment,['TET at hour ', tetTime]),contains(gFitData.cell,'139A'),...
%     (gFitData.NumGrans_f >= minGrans)],2);
% figure;
% boxplot(gFitData.NumGrans_f(subz),gFitData.treatment(subz));
% xlabel('treatment'); ylabel('f (Number of 4B Granules)')
% ylim(ylimz)
% title([tetTime,'h 139A mut Tet induced and all NaAsO2 Doses'])

