%% Load the appended and aligned Live_cell last frame and IF data

% 2023-06-29 Data
basePath = 'Z:\imageData\SG_4B\';
dataset2=load([basePath,'2023-06-29 4B WT vs Mut TET curve NaAsO2 curve\2023-06-29 4B WT vs Mut TET curve NaAsO2 curve_Processed_Copy.mat']);
dataset2 = dataset2.dataloc; % pull the loaded dataloc structure
dataloc = dataset2;
% 
% % Pull and combine the IF datasets
% pooledData = [dataset1.ifd,dataset2.ifd];


%% Figure 3 - RRM of 4B affects Translation and Translational haulting in stress conditions

% Figure 3A - % OPP intensity +/- stress for wt-4b and mut-4b vs wt-4B-GFP (as control) (24 hour induced only)

% 24h TET w/ no NaAsO2 and 24h TET with 125uM NaAsO2
% 4B wt or mut grans vs OPP %

tetTime = "hour -24"; % amount of tet induction
naAsO2 = (" 0uM"|" 125uM"); % NaAsO2 treatment

subF3A = all([contains(dataloc.ifd.treatment,tetTime),~contains(dataloc.ifd.cell,'bad'),contains(dataloc.ifd.treatment,naAsO2) ...
    ],2); % , (dataloc.ifd.Children_Grans_G3BP1_Count >= minGrans)

subF3Adata = dataloc.ifd(subF3A,:); % pull that subset of data

subF3Adata.cell = categorical(subF3Adata.cell,{'HeLa eIF4BGFP','HeLa 4B139AGFP'}); % make the cell line a categorical

% make the treatments a categorical
catOrder = ["0.1ug/mL TET at hour -24 and 0uM NaAsO2","0.1ug/mL TET at hour -24 and 125uM NaAsO2"]; % give the order I want for the plot
subF3Adata.treatment = categorical(subF3Adata.treatment,catOrder); % set the order
subF3Adata.OPP = subF3Adata.Intensity_MeanIntensity_Masked_OPP*65535; % get the real opp values

% Get the mean OPP intensity of everyone by treatment and cell line
generalOPPData = grpstats(subF3Adata,["cell","treatment"],["mean","median","sem","std"],"DataVars","OPP")

%% Do the % OPP intensity of everyone vs the mean OPP intensity of the 24 TET induced vehicle control wt-4b cells

% pull the mean OPP of vehicle treated wt-4b cells
meanOppWt4bNoAs = generalOPPData{"HeLa eIF4BGFP_0.1ug/mL TET at hour -24 and 0uM NaAsO2","mean_OPP"}; 

% divide all of the OPP data by this number * 100 for % mean contol opp
subF3Adata.PercOPP = (subF3Adata.OPP/meanOppWt4bNoAs)*100;

% get the general values of the % differences in OPP
generalPercOPPData = grpstats(subF3Adata,["cell","treatment"],["mean","median","sem","std"],"DataVars","PercOPP")

% Do the statistics for % OPP versus each cell line or treatment (WT-4B 0uM NaAsO2 as control)
[~,~,statsOPP] = anova1(subF3Adata.PercOPP,join(string([subF3Adata.cell,subF3Adata.treatment])),'off');
[resultsOPP,~,~,gnamesOPP] = multcompare(statsOPP,"CriticalValueType","dunnett",'ControlGroup',find(matches(statsOPP.gnames,'HeLa eIF4BGFP 0.1ug/mL TET at hour -24 and 0uM NaAsO2')),'Display','off','Approximate',false);
PercOPP = array2table(resultsOPP,"VariableNames", ["Group","Control Group","Lower Limit","Difference","Upper Limit","P-value"]);
PercOPP.("Group") = gnamesOPP(PercOPP.("Group"));
PercOPP.("Control Group") = gnamesOPP(PercOPP.("Control Group"))

% Do the statistics for % OPP versus each cell line or treatment (MUT-4B 125uM NaAsO2 as control)
[~,~,statsOPP] = anova1(subF3Adata.PercOPP,join(string([subF3Adata.cell,subF3Adata.treatment])),'off');
[resultsOPP,~,~,gnamesOPP] = multcompare(statsOPP,"CriticalValueType","dunnett",'ControlGroup',find(matches(statsOPP.gnames,'HeLa 4B139AGFP 0.1ug/mL TET at hour -24 and 125uM NaAsO2')),'Display','off','Approximate',false);
PercOPP = array2table(resultsOPP,"VariableNames", ["Group","Control Group","Lower Limit","Difference","Upper Limit","P-value"]);
PercOPP.("Group") = gnamesOPP(PercOPP.("Group"));
PercOPP.("Control Group") = gnamesOPP(PercOPP.("Control Group"))


%% Make the plots
% make a new figure
figure3 = figure;
figure3.Units = "Inches";
figure3.Position = [0.05,0.05,8.5,11];

topP =  uipanel('Units','inches','position',[0.5,  6.75,  7.6,     2.5]); % location for F3A
botLP = uipanel('Units','inches','position',[0.5,   3,    3.75,    3.5]); % bottom left panel (left half of F3b)
botRP = uipanel('Units','inches','position',[4.35,  3,    3.75,    3.5]); % bottom right panel (right half of F3b)

% make a new container for axes 
f3axes = axes(topP);

% Unique categories for legend
uniqueTx = unique(generalPercOPPData.treatment,'stable');

% Colors for each treatment conc
colorz = [0,       0.4470,  0.7410 ;... % 0uM Color
          0.8500,  0.3250,  0.0980]; % 125uM Color

colorzScatter = [0.6350, 0.0780,  0.1840; ...	% 250uM  Color
                 0.8500,  0.3250,  0.0980;...    % 125uM Color
                 0.9290, 0.6940,  0.1250;] ;  	% 62.5uM Color


% Bar locations
barLocations = [1, 2, 4, 5];

hold on;
% Plot F3A
for i = 1:height(generalPercOPPData)
    bar(barLocations(i), generalPercOPPData.mean_PercOPP(i), 'facecolor', colorz(matches(string(uniqueTx),string(generalPercOPPData.treatment(i))),:));
end
% Error bars
errorbar(barLocations, generalPercOPPData.mean_PercOPP, generalPercOPPData.std_PercOPP, 'k', 'linestyle', 'none');

% Customizing the plot
set(gca, 'XTick', [1.5, 4.5], 'XTickLabel', {'Wt', 'Mut'}); % Adjusting X-ticks to match cell lines
legend(strrep(string(uniqueTx),'0.1ug/mL TET at hour -24 and ',''), 'Location', 'Best'); % Only show legend for unique Group tx categories
ylabel("Percent OPP Intensity vs Wt-4b Control"); xlabel('4B Cell line') 
hold off;
ylim([0,110]); xlim([0,6]);
fontname('Arial'); fontsize(8,"points");


% Plot F3b left half (24 hr tet wt-4b-gfp)

% 24 hour tet wt 4b grans or 139a vs opp vs naaso2 dose/vehicle with box and wiskers?
minGrans = 3;
tetTime = '-24';
xlimz = [0,30]; % granule axis limits
ylimz = [7.5, 12.5]; % Opp axis limits


subzWT = all([contains(dataloc.ifd.treatment,['TET at hour ', tetTime]),contains(dataloc.ifd.cell,'4BGFP')...
    , (dataloc.ifd.Children_Grans_4B_Count >= minGrans), ~contains(dataloc.ifd.treatment,'and 0uM')],2); % 

% pull the wt data
wtData = dataloc.ifd(subzWT,:);

% make the scatter hist for the wt data (bottom left panel)
h=scatterhist(wtData.Children_Grans_4B_Count,log2(wtData.Intensity_MeanIntensity_Masked_OPP*65535),...
    'Group',wtData.treatment,'parent',botLP,'MarkerSize',3,'color',colorzScatter);
xlabel('WT-4B Granules'); ylabel('Log_2 OPP Intensity')
title('WT4b all NaAsO2 Doses')
hold on;
boxplot(h(2),wtData.Children_Grans_4B_Count,wtData.treatment,'orientation','horizontal',...
     'label',{'','',''},'color',colorzScatter);
boxplot(h(3),log2(wtData.Intensity_MeanIntensity_Masked_OPP*65535),wtData.treatment,'orientation','horizontal',...
     'label', {'','',''},'color',colorzScatter);
set(h(2:3),'XTickLabel','');
view(h(3),[270,90]);  % Rotate the Y plot
xlim(h(1),xlimz); ylim(h(1),ylimz);
xlim(h(2),xlimz); xlim(h(3),ylimz); % Sync axes
hold off;
fontname('Arial'); fontsize(8,"points");
lg = legend(strrep(string(unique(wtData.treatment,'stable')),'0.1ug/mL TET at hour -24 and ',''), 'Location', 'Best'); % Only show legend for unique Group tx categories
lg.Position = [0.1,    0.1327,    0.1,    0.1325]; % set the legend size

subzMut = all([contains(dataloc.ifd.treatment,['TET at hour ', tetTime]),contains(dataloc.ifd.cell,'139A')...
    ,(dataloc.ifd.Children_Grans_4B_Count >= minGrans), ~contains(dataloc.ifd.treatment,'and 0uM')],2); % Mutant Data

% pull the mutant data
mutData = dataloc.ifd(subzMut,:);

% make the scatter hist for the mutant data (bottom right panel)
h=scatterhist(mutData.Children_Grans_4B_Count,log2(mutData.Intensity_MeanIntensity_Masked_OPP*65535),...
    'Group',mutData.treatment,'parent',botRP,'MarkerSize',3,'color',colorzScatter);
xlabel('Mut 4B Granules'); ylabel('Log_2 OPP Intensity')
title('137A all NaAsO2 Doses')
hold on;

boxplot(h(2),mutData.Children_Grans_4B_Count,mutData.treatment,'orientation','horizontal',...
     'label',{'','',''},'color',colorzScatter);
boxplot(h(3),log2(mutData.Intensity_MeanIntensity_Masked_OPP*65535),mutData.treatment,'orientation','horizontal',...
     'label', {'','',''},'color',colorzScatter);
set(h(2:3),'XTickLabel','');
view(h(3),[270,90]);  % Rotate the Y plot
xlim(h(1),xlimz); ylim(h(1),ylimz);
xlim(h(2),xlimz); xlim(h(3),ylimz);  % Sync axes
hold off;
fontname('Arial'); fontsize(8,"points");
legend off;

figure3.Color = [1,1,1];
topP.BackgroundColor = [1,1,1];
botLP.BackgroundColor = [1,1,1];
botRP.BackgroundColor = [1,1,1];


saveas(figure3,'Z:\imageData\SG_4B\Paper_Figures\Output_Figures\Figure_3.fig')
saveas(figure3,'Z:\imageData\SG_4B\Paper_Figures\Output_Figures\Figure_3.svg')

%% Now print the saw OPP values per cell line per NaAsO2 treatment

subz = ~contains(dataloc.ifd.cell,'bad','ignorecase',true); % filter for the parameters set above

orderOfPlot = {...
    '0ug/mL TET at hour 0 and 0uM NaAsO2', '0.1ug/mL TET at hour -2 and 0uM NaAsO2', '0.1ug/mL TET at hour -4 and 0uM NaAsO2','0.1ug/mL TET at hour -8 and 0uM NaAsO2', '0.1ug/mL TET at hour -14 and 0uM NaAsO2', '0.1ug/mL TET at hour -24 and 0uM NaAsO2',...
    '0ug/mL TET at hour 0 and 62.5uM NaAsO2', '0.1ug/mL TET at hour -2 and 62.5uM NaAsO2', '0.1ug/mL TET at hour -4 and 62.5uM NaAsO2','0.1ug/mL TET at hour -8 and 62.5uM NaAsO2', '0.1ug/mL TET at hour -14 and 62.5uM NaAsO2', '0.1ug/mL TET at hour -24 and 62.5uM NaAsO2',...
    '0ug/mL TET at hour 0 and 125uM NaAsO2', '0.1ug/mL TET at hour -2 and 125uM NaAsO2', '0.1ug/mL TET at hour -4 and 125uM NaAsO2','0.1ug/mL TET at hour -8 and 125uM NaAsO2', '0.1ug/mL TET at hour -14 and 125uM NaAsO2', '0.1ug/mL TET at hour -24 and 125uM NaAsO2',...
    '0ug/mL TET at hour 0 and 250uM NaAsO2', '0.1ug/mL TET at hour -2 and 250uM NaAsO2', '0.1ug/mL TET at hour -4 and 250uM NaAsO2','0.1ug/mL TET at hour -8 and 250uM NaAsO2', '0.1ug/mL TET at hour -14 and 250uM NaAsO2', '0.1ug/mL TET at hour -24 and 250uM NaAsO2',...
    };

allOPPData = dataloc.ifd(subz,:);
allOPPData.txcat = categorical(allOPPData.treatment,orderOfPlot);
allOPPData.OPP = allOPPData.Intensity_MeanIntensity_Masked_OPP*65535;


allOPPStats = grpstats(allOPPData,["txcat","cell"],["mean","median","sem","std"],"DataVars","OPP")

% save it as a csv 
writetable(allOPPStats,'Z:\imageData\SG_4B\Paper_Figures\Output_Figures\allOppStats.csv')

%% Run the statistics comparing the wt-4B treated cells accross various concentrations of NaAsO2 and control

% get averything for the reader if needed
grpstats(wtSubData,"treatment",["mean","median","sem","std"],"DataVars",["NumGrans_rate_in_min","NumGrans_f","NumGrans_min_to_respond"])

% Now just print the means
grpstats(wtSubData,"treatment","mean","DataVars",["NumGrans_rate_in_min","NumGrans_f","NumGrans_min_to_respond"])


%% X hour tet wt 4b or 139a 4b vs opp @ each naaso2 dose


dataloc.ifd.txcat = categorical(dataloc.ifd.treatment,orderOfPlot);

subzwt = contains(dataloc.ifd.cell,'4BGFP'); % subset the data to be wt 4b
subzmut = contains(dataloc.ifd.cell,'139A'); % subset the data to be 139A 4b

barColz = [ 0.3010, 0.7450, 0.9330;...
            0.6350, 0.0780, 0.1840;...
            0.9290, 0.6940, 0.1250;...
            0, 0.4470, 0.7410;...
            0.4660, 0.6740, 0.1880;...
            0.4940, 0.1840, 0.5560]; % colors for bars

figure;
boxplot(log2(dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subzwt)*65535),dataloc.ifd.txcat(subzwt),'Notch','on','Symbol','.','Colors',barColz)
ylim([7, 15])
xlabel('TET Induction'); ylabel('Log_2 OPP Intensity')
title('WT4b Tet induction versus OPP')
legend

figure;
boxplot(log2(dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subzmut)*65535),dataloc.ifd.txcat(subzmut),'Notch','on','Symbol','.','Colors',barColz)
xlabel('TET Induction'); ylabel('Log_2 OPP Intensity')
ylim([7,15])
title('139A Tet induction versus OPP')
%% % OPP by dose
% figure;
% minGrans = 0; % min number of granules
% tetTime = "hour -24"; % amount of tet induction
% % all NaAsO2 treatments
% 
% subF3A = all([contains(dataloc.ifd.treatment,tetTime),~contains(dataloc.ifd.cell,'bad'),...
%     (dataloc.ifd.Children_Grans_G3BP1_Count >= minGrans)],2);
% 
% subF3Adata = dataloc.ifd(subF3A,:); % pull that subset of data
% 
% subF3Adata.cell = categorical(subF3Adata.cell,{'HeLa eIF4BGFP','HeLa 4B139AGFP'}); % make the cell line a categorical
% 
% % make the treatments a categorical
% catOrder = ["0.1ug/mL TET at hour -24 and 0uM NaAsO2","0.1ug/mL TET at hour -24 and 125uM NaAsO2"]; % give the order I want for the plot
% subF3Adata.treatment = categorical(subF3Adata.treatment,catOrder); % set the order
% subF3Adata.OPP = subF3Adata.Intensity_MeanIntensity_Masked_OPP*65535; % get the real opp values
% 
% % Get the mean OPP intensity of everyone by treatment and cell line
% generalOPPData = grpstats(subF3Adata,["cell","treatment"],["mean","median","sem","std"],"DataVars","OPP")