%% Process the data
fts(1).t = 'Min'; fts(1).c = 'Granularity_1_EnhancedGreen'; fts(1).p = 0;
filterp.dlength = 2;
filterp.fts = fts;

dataloc = JB_DataHandler('reorganizedata',true,'filterp',filterp); 

%% Load the data
%dataloc = JB_DataHandler(); 

%% Plot the data 
plotme = {'Granularity_10_AlignedGreen', 'Granularity_11_AlignedGreen', 'Granularity_12_AlignedGreen', 'Granularity_13_AlignedGreen', 'Granularity_14_AlignedGreen', 'Granularity_15_AlignedGreen', 'Granularity_16_AlignedGreen', 'Granularity_1_AlignedGreen', 'Granularity_2_AlignedGreen', 'Granularity_3_AlignedGreen', 'Granularity_4_AlignedGreen', 'Granularity_5_AlignedGreen', 'Granularity_6_AlignedGreen', 'Granularity_7_AlignedGreen', 'Granularity_8_AlignedGreen', 'Granularity_9_AlignedGreen','Granularity_10_EnhancedGreen', 'Granularity_9_EnhancedGreen', 'Granularity_8_EnhancedGreen', 'Granularity_7_EnhancedGreen', 'Granularity_6_EnhancedGreen', 'Granularity_5_EnhancedGreen', 'Granularity_4_EnhancedGreen', 'Granularity_3_EnhancedGreen', 'Granularity_2_EnhancedGreen', 'Granularity_1_EnhancedGreen', 'Granularity_16_EnhancedGreen', 'Granularity_15_EnhancedGreen', 'Granularity_14_EnhancedGreen', 'Granularity_13_EnhancedGreen', 'Granularity_12_EnhancedGreen', 'Granularity_11_EnhancedGreen'};
plot_by_ND('treatment', dataloc,'plottype',{'meanslope', 'heatmap'},'channel',plotme,'looptime',4)