%% Figure S1 

% Load the model fit data for the live cell
addpath('Z:\code\Nick')

%% Plot the Pearson's for G3BP1 and 4B-GFP for both cell lines
coLoc = readtable('Z:\imageData\SG_4B\Paper_Figures\Output_Figures\G3BP1_vs_4B_Colocalization_Data.xlsx');

% pull the data for the 125uM NaAsO2 treated cells
coLoc = coLoc(coLoc.NaAsO2Dose == 125,:);

coLoc.Cell = categorical(coLoc.Cell,{'Wt-4B','Mut-4B'})

% now make a bar graph with the two cell lines and their pearsons with each
% dot being a replicate

coLocStats = grpstats(coLoc,"Cell",["mean","sem"],"DataVars",["Pearson_sCorrelation","Mander_sOverlap"])

supF1 = figure;
bar(coLocStats.Cell,coLocStats.mean_Mander_sOverlap)
hold on;
errorbar(coLocStats.mean_Mander_sOverlap,coLocStats.sem_Mander_sOverlap,'LineStyle','none','LineWidth',2,'CapSize',15,'Color','k')
swarmchart(coLoc,"Cell","Mander_sOverlap",'XJitter','density','YJitter','none','MarkerEdgeColor','k','MarkerFaceColor',[0.5,0.5,0.5],'XJitterWidth',0.25)
ylim([0,1])

saveas(supF1,'Z:\imageData\SG_4B\Paper_Figures\Output_Figures\Figure_S1.fig')
saveas(supF1,'Z:\imageData\SG_4B\Paper_Figures\Output_Figures\Figure_S1.svg')

