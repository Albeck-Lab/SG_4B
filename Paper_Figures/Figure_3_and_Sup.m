%% Load the appended and aligned Live_cell last frame and IF data

% 2023-06-29 Data
basePath = 'Z:\imageData\SG_4B\';
dataset1=load([basePath,'2023-06-29 4B WT vs Mut TET curve NaAsO2 curve\2023-06-29 4B WT vs Mut TET curve NaAsO2 curve_Processed.mat']); % _Copy
dataset1 = dataset1.dataloc; % pull the loaded dataloc structure

% append the movie's metadata
dataset1.movieinfo.PixSizeX = 0.33; % um/px
dataset1.movieinfo.PixSizeY = 0.33; % um/px
dataset1.movieinfo.PixNumX = 1600; % pixels
dataset1.movieinfo.PixNumY = 1600; % pixels
dataset1.movieinfo.tsamp = 3; % minutes


%% Figure 3 - RRM of 4B affects Translation and Translational haulting in stress conditions

% Figure 3A - % OPP intensity +/- stress for wt-4b and mut-4b vs wt-4B-GFP (as control) (24 hour induced only)

% 24h TET w/ no NaAsO2 and 24h TET with 125uM NaAsO2
% 4B wt or mut grans vs OPP %

tetTime = "hour -24"; % amount of tet induction
naAsO2 = (" 0uM"|" 125uM"); % NaAsO2 treatment

subF3A = all([contains(dataset1.ifd.treatment,tetTime),~contains(dataset1.ifd.cell,'bad'),contains(dataset1.ifd.treatment,naAsO2) ...
    ],2); % , (dataloc.ifd.Children_Grans_G3BP1_Count >= minGrans)

subF3Adata = dataset1.ifd(subF3A,:); % pull that subset of data

% make the treatments a categorical and simplify them
subF3Adata.treatment = strrep(subF3Adata.treatment,'0.1ug/mL TET at hour -24 and ','');
catOrder = ["0uM NaAsO2","125uM NaAsO2"]; % give the order I want for the plot
subF3Adata.treatment = categorical(subF3Adata.treatment,catOrder); % set the order

% make the cell lines a categorical and simplify them
subF3Adata.cell = strrep(subF3Adata.cell,'HeLa eIF4BGFP','Wt');
subF3Adata.cell = strrep(subF3Adata.cell,'HeLa 4B139AGFP','Mt');
subF3Adata.cell = categorical(subF3Adata.cell,{'Wt','Mt'}); % make the cell line a categorical

% get the real opp values
subF3Adata.OPP = subF3Adata.Intensity_MeanIntensity_Masked_OPP*65535; 

% and the log2 OPP
subF3Adata.l2OPP = log2(subF3Adata.OPP); 

% Get the mean OPP intensity of everyone by treatment and cell line
generalOPPData = grpstats(subF3Adata,["cell","treatment"],["mean","median","sem","std"],"DataVars",["OPP","l2OPP"])

%% Do the % OPP intensity of everyone vs the mean OPP intensity of the 24 TET induced vehicle control wt-4b cells

% pull the mean OPP of vehicle treated wt-4b cells
meanOppWt4bNoAs = generalOPPData{"Wt_0uM NaAsO2","mean_OPP"}; 

% divide all of the OPP data by this number * 100 for % mean contol opp
subF3Adata.PercOPP = (subF3Adata.OPP/meanOppWt4bNoAs)*100;

% get the general values of the % differences in OPP
generalPercOPPData = grpstats(subF3Adata,["cell","treatment"],["mean","median","sem","std"],"DataVars","PercOPP")

% Do the statistics for % OPP versus each cell line or treatment (WT-4B 0uM NaAsO2 as control)
[~,~,statsOPP] = anova1(subF3Adata.PercOPP,join(string([subF3Adata.cell,subF3Adata.treatment])),'off');
[resultsOPP,~,~,gnamesOPP] = multcompare(statsOPP,"CriticalValueType","dunnett",'ControlGroup',find(matches(statsOPP.gnames,'Wt 0uM NaAsO2')),'Display','off','Approximate',false);
PercOPP = array2table(resultsOPP,"VariableNames", ["Group","Control Group","Lower Limit","Difference","Upper Limit","P-value"]);
PercOPP.("Group") = gnamesOPP(PercOPP.("Group"));
PercOPP.("Control Group") = gnamesOPP(PercOPP.("Control Group"))

% Do the statistics for % OPP versus each cell line or treatment (MUT-4B 125uM NaAsO2 as control)
[~,~,statsOPP] = anova1(subF3Adata.PercOPP,join(string([subF3Adata.cell,subF3Adata.treatment])),'off');
[resultsOPP,~,~,gnamesOPP] = multcompare(statsOPP,"CriticalValueType","dunnett",'ControlGroup',find(matches(statsOPP.gnames,'Mt 125uM NaAsO2')),'Display','off','Approximate',false);
PercOPP = array2table(resultsOPP,"VariableNames", ["Group","Control Group","Lower Limit","Difference","Upper Limit","P-value"]);
PercOPP.("Group") = gnamesOPP(PercOPP.("Group"));
PercOPP.("Control Group") = gnamesOPP(PercOPP.("Control Group"))

%% Figure 3D - PLSR of Granule kinetics and OPP outcome 

[fitData2,~] = convertDatalocToModelFit({dataset1}, 'NumGrans','pulsepars',{'f','td','ts','rate_in_min','min_to_respond','rsquared','granarea'},'afterf',1);

%% Subset only the good data 
fitData = fitData2; % work with duplicated data (for safety)
gFitData2 = fitData((fitData.NumGrans_rsquared > 0.8),:); % look for an r squared greater than 0.8 ? 

gFitData2.cell = strrep(gFitData2.cell,'HeLa_eIF4BGFP','Wt');
gFitData2.cell = strrep(gFitData2.cell,'HeLa_4B139AGFP','Mut');

gFitData2.xy = gFitData2.NumGrans_xy;
gFitData2.cellid = gFitData2.NumGrans_cellid;

ifD = dataset1.ifd;
ifD.xy = ifD.ImageNumber;
missinD = isnan(ifD.cellid);
ifD2 = ifD(~missinD,:);

[~, idxG] = ismember(gFitData2(:,["cellid","xy"]), ifD2(:,["cellid","xy"]), 'rows');
[~, idxI] = ismember(ifD2(:,["cellid","xy"]),gFitData2(:,["cellid","xy"]), 'rows');

gFitData2 = gFitData2(idxG > 0, :);
ifD2 = ifD2(idxI > 0, :);

fAndIF = join(gFitData2,ifD2,'LeftKeys',{'xy','cellid'},'RightKeys',{'xy','cellid'},'KeepOneCopy',{'treatment','cell','full'});

%% Collect the data for use in PLSR

minGrans = 3;
tetTime = "hour -24"; % amount of tet induction
naAsO2 = (" 62.5uM"|" 125uM"|" 250uM"); % NaAsO2 treatment 

subF3D = all([contains(fAndIF.treatment,tetTime),~contains(fAndIF.cell,'bad'),contains(fAndIF.treatment,naAsO2), ...
    (fAndIF.NumGrans_f >= minGrans)],2); %,contains(fAndIF.cell,'HeLa_eIF4BGFP')

subF3Ddata = fAndIF(subF3D,:); % pull that subset of data
subF3Ddata.cell = categorical(subF3Ddata.cell,{'Wt','Mut'}); % make the cell line a categorical
subF3Ddata.treatment = strrep(subF3Ddata.treatment,'0.1ug/mL TET at hour -24 and ','');
subF3Ddata.treatment = strrep(subF3Ddata.treatment,' NaAsO2 at hour 0','');
subF3Ddata.treatment = categorical(subF3Ddata.treatment,{'62.5uM','125uM','250uM'}); % set the order '62.5uM'

% get the real opp values and the log2 opp values
subF3Ddata.OPP = subF3Ddata.Intensity_MeanIntensity_Masked_OPP*65535; % get the real opp values
subF3Ddata.l2OPP = log2(subF3Ddata.OPP);

% make the treatment a real number
subF3Ddata.dose = str2double(strrep(string(subF3Ddata.treatment),'uM',''));

% make Wt vs Mut a number for plsr
subF3Ddata.celln = zeros([height(subF3Ddata),1]); % wt is 0
subF3Ddata{(subF3Ddata.cell=='Mut'),"celln"} = ones([sum((subF3Ddata.cell=='Mut')),1]); % mut is 1

plsOut = []; 
plsOut{1} = pls([subF3Ddata.dose,subF3Ddata.NumGrans_f,subF3Ddata.NumGrans_rate_in_min,subF3Ddata.NumGrans_min_to_respond,subF3Ddata.celln],...
    subF3Ddata.l2OPP, 'params', {'dose','F','rate','t2r','cellLine'},'ploton',false);

% remove only dose
plsOut{2} = pls([subF3Ddata.NumGrans_f,subF3Ddata.NumGrans_rate_in_min,subF3Ddata.NumGrans_min_to_respond,subF3Ddata.celln],...
    subF3Ddata.l2OPP, 'params', {'F','rate','t2r','cellLine'},'ploton',false);

% just dose
plsOut{3} = pls(subF3Ddata.dose,...
    subF3Ddata.l2OPP, 'params', {'dose'},'ploton',false);

% just F
plsOut{4} = pls(subF3Ddata.NumGrans_f,...
    subF3Ddata.l2OPP, 'params', {'F'},'ploton',false);

% just rate
plsOut{5} = pls(subF3Ddata.NumGrans_rate_in_min,...
    subF3Ddata.l2OPP, 'params', {'rate'},'ploton',false);

% just t2r
plsOut{6} = pls(subF3Ddata.NumGrans_min_to_respond,...
    subF3Ddata.l2OPP, 'params', {'t2r'},'ploton',false);

% just cell line
plsOut{7} = pls(subF3Ddata.celln,...
    subF3Ddata.l2OPP, 'params', {'cellLine'},'ploton',false);


% make a list for the legend
legList = {'dose+F+rate+t2r+cell','F+rate+t2r+cell','dose','F','rate','t2r','cell'};

% Now print out the variance explained as a table
varExpTable = legList;

for iPLSR = 1:numel(plsOut)
    varExpTable{2,iPLSR} = max(cumsum(plsOut{iPLSR}.PCTVAR(2,:))*100,[],2,"omitnan");
end
varExpTable = cell2table(varExpTable);
writetable(varExpTable,'Z:\imageData\SG_4B\Paper_Figures\Output_Figures\PlsrPercVarExplained.csv')

%% Prepare wt and mut data for Figure 3B and 3C plots
minGrans = 3;
tetTime = '-24';

% wt data
subzWT = all([contains(dataset1.ifd.treatment,['TET at hour ', tetTime]),contains(dataset1.ifd.cell,'4BGFP')...
    , (dataset1.ifd.Children_Grans_4B_Count >= minGrans), ~contains(dataset1.ifd.treatment,'and 0uM')],2); % 

% pull the wt data, get the real opp and l2opp
wtData = dataset1.ifd(subzWT,:);
wtData.OPP = wtData.Intensity_MeanIntensity_Masked_OPP*65535;
wtData.l2OPP = log2(wtData.OPP);

%simplify treatment name
wtData.treatment = strrep(wtData.treatment,'0.1ug/mL TET at hour -24 and ','');
wtData.treatment = strrep(wtData.treatment,' NaAsO2','');

% get the general values of the wt l2OPP
WTl2OppData = grpstats(wtData,"treatment",["mean","median","sem","std"],"DataVars","l2OPP")

% Do the statistics for l2 OPP versus treatment (WT-4B 125uM NaAsO2 as control)
[~,~,statsOPP] = anova1(wtData.l2OPP,wtData.treatment,'off');
[resultsOPP,~,~,gnamesOPP] = multcompare(statsOPP,"CriticalValueType","dunnett",'ControlGroup',find(matches(statsOPP.gnames,'125uM')),'Display','off','Approximate',false);
wtL2OPP = array2table(resultsOPP,"VariableNames", ["Group","Control Group","Lower Limit","Difference","Upper Limit","P-value"]);
wtL2OPP.("Group") = gnamesOPP(wtL2OPP.("Group"));
wtL2OPP.("Control Group") = gnamesOPP(wtL2OPP.("Control Group"))

% mutant data
subzMut = all([contains(dataset1.ifd.treatment,['TET at hour ', tetTime]),contains(dataset1.ifd.cell,'139A')...
    ,(dataset1.ifd.Children_Grans_4B_Count >= minGrans), ~contains(dataset1.ifd.treatment,'and 0uM')],2); % Mutant Data

% pull the mutant data, get the real opp and l2opp
mutData = dataset1.ifd(subzMut,:);
mutData.OPP = mutData.Intensity_MeanIntensity_Masked_OPP*65535;
mutData.l2OPP = log2(mutData.OPP);

%simplify treatment name
mutData.treatment = strrep(mutData.treatment,'0.1ug/mL TET at hour -24 and ','');
mutData.treatment = strrep(mutData.treatment,' NaAsO2','');

% get the general values of the mut l2OPP
MUTl2OppData = grpstats(mutData,"treatment",["mean","median","sem","std"],"DataVars","l2OPP")

% Do the statistics for l2 OPP versus treatment (WT-4B 125uM NaAsO2 as control)
[~,~,statsOPP] = anova1(mutData.l2OPP,mutData.treatment,'off');
[resultsOPP,~,~,gnamesOPP] = multcompare(statsOPP,"CriticalValueType","dunnett",'ControlGroup',find(matches(statsOPP.gnames,'125uM')),'Display','off','Approximate',false);
mutL2OPP = array2table(resultsOPP,"VariableNames", ["Group","Control Group","Lower Limit","Difference","Upper Limit","P-value"]);
mutL2OPP.("Group") = gnamesOPP(mutL2OPP.("Group"));
mutL2OPP.("Control Group") = gnamesOPP(mutL2OPP.("Control Group"))

%% Check all l2 opp data against eachother

minGrans = 3;
tetTime = '-24';

% wt and mut data
subzBoth = all([contains(dataset1.ifd.treatment,['TET at hour ', tetTime]), ~contains(dataset1.ifd.cell,'bad')...
    , (dataset1.ifd.Children_Grans_4B_Count >= minGrans), ~contains(dataset1.ifd.treatment,'and 0uM')],2); % 


% pull the mutant data, get the real opp and l2opp
allData = dataset1.ifd(subzBoth,:);
allData.OPP = allData.Intensity_MeanIntensity_Masked_OPP*65535;
allData.l2OPP = log2(allData.OPP);

%simplify treatment name
allData.treatment = strrep(allData.treatment,'0.1ug/mL TET at hour -24 and ','');
allData.treatment = strrep(allData.treatment,' NaAsO2','');

%simplify cell names
allData.cell = strrep(allData.cell,"HeLa 4B139AGFP",'mut');
allData.cell = strrep(allData.cell,"HeLa eIF4BGFP",'wt');

% get the general values of the mut l2OPP
MUTl2OppData = grpstats(allData,["treatment","cell"],["mean","median","sem","std"],"DataVars","l2OPP")

% Do the statistics for l2 OPP versus treatment (WT-4B 125uM NaAsO2 as control)
[~,~,statsOPP] = anova1(allData.l2OPP,allData.treatment,'off');
[resultsOPP,~,~,gnamesOPP] = multcompare(statsOPP,"CriticalValueType","dunnett",'ControlGroup',find(matches(statsOPP.gnames,'125uM')),'Display','off','Approximate',false);
mutL2OPP = array2table(resultsOPP,"VariableNames", ["Group","Control Group","Lower Limit","Difference","Upper Limit","P-value"]);
mutL2OPP.("Group") = gnamesOPP(mutL2OPP.("Group"));
mutL2OPP.("Control Group") = gnamesOPP(mutL2OPP.("Control Group"))


%% Make the plots
% make a new figure
figure3 = figure;
figure3.Units = "Inches";
figure3.Position = [0.05,0.05,8.5,11];

topP =  uipanel('Units','inches','position',[0.5,  7.75,  7.6,     2.5]); % location for F3A
botLP = uipanel('Units','inches','position',[0.5,   4.125,    3.75,    3.5]); % bottom left panel (left half of F3b)
botRP = uipanel('Units','inches','position',[4.35,  4.125,    3.75,    3.5]); % bottom right panel (right half of F3b)
botBot = uipanel('Units','inches','position',[0.5,  1.5,    7.6,    2.5]); % PLSR data

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
errorbar(barLocations, generalPercOPPData.mean_PercOPP, generalPercOPPData.sem_PercOPP, 'k', 'linestyle', 'none','CapSize',15,'LineWidth',1);

% Customizing the plot
set(gca, 'XTick', [1.5, 4.5], 'XTickLabel', {'Wt', 'Mut'}); % Adjusting X-ticks to match cell lines
legend(string(uniqueTx), 'Location', 'Best'); % Only show legend for unique Group tx categories
ylabel("Percent OPP Intensity vs Wt-4b Control"); xlabel('4B Cell line') 
hold off;
ylim([0,110]); xlim([0,6]);
fontname('Arial'); fontsize(8,"points");


% Plot F3b left half (24 hr tet wt-4b-gfp)

xlimz = [0,30]; % granule axis limits
ylimz = [7.5, 12]; % Opp axis limits

% Make the scatter hist for the wt data (bottom left panel)
h=scatterhist(wtData.Children_Grans_4B_Count,wtData.l2OPP,...
    'Group',wtData.treatment,'parent',botLP,'MarkerSize',3,'color',colorzScatter);
xlabel('WT-4B Granules'); ylabel('Log_2 OPP Intensity')
title('WT4b all NaAsO2 Doses')
hold on;
boxplot(h(2),wtData.Children_Grans_4B_Count,wtData.treatment,'orientation','horizontal',...
     'label',{'','',''},'color',colorzScatter,'Symbol','');
boxplot(h(3),wtData.l2OPP,wtData.treatment,'orientation','horizontal',...
     'label', {'','',''},'color',colorzScatter,'Symbol','');
set(h(2:3),'XTickLabel','');
view(h(3),[270,90]);  % Rotate the Y plot
xlim(h(1),xlimz); ylim(h(1),ylimz);
xlim(h(2),xlimz); xlim(h(3),ylimz); % Sync axes
hold off;
fontname('Arial'); fontsize(8,"points");
lg = legend(unique(wtData.treatment), 'Location', 'Best'); % Only show legend for unique Group tx categories
lg.Position = [0.1,    0.1327,    0.1,    0.1325]; % set the legend size

% make the scatter hist for the mutant data (bottom right panel)
h=scatterhist(mutData.Children_Grans_4B_Count,mutData.l2OPP,...
    'Group',mutData.treatment,'parent',botRP,'MarkerSize',3,'color',colorzScatter);
xlabel('Mut 4B Granules'); ylabel('Log_2 OPP Intensity')
title('137A all NaAsO2 Doses')
hold on;

boxplot(h(2),mutData.Children_Grans_4B_Count,mutData.treatment,'orientation','horizontal',...
     'label',{'','',''},'color',colorzScatter,'Symbol','');
boxplot(h(3),mutData.l2OPP,mutData.treatment,'orientation','horizontal',...
     'label', {'','',''},'color',colorzScatter,'Symbol',''); 
set(h(2:3),'XTickLabel','');
view(h(3),[270,90]);  % Rotate the Y plot
xlim(h(1),xlimz); ylim(h(1),ylimz);
xlim(h(2),xlimz); xlim(h(3),ylimz);  % Sync axes
hold off;
fontname('Arial'); fontsize(8,"points");
legend off;


% Plot PLS model vs. data
%Plot: Scatter, Variance explained, Component weights, Parameter weights
%   Scatter plot
% 
% This is code directly taken from plsplot written by M. Pargett from the Albeck Lab

p.nout = 1;
p.alpha = 0.1;

%   Get modeled values and axis limits for both model and data
mapped_vals = [ones(size(plsOut{1}.X,1),1), plsOut{1}.X]*plsOut{1}.BETA;
axmin = min([plsOut{1}.Y(:,p.nout); mapped_vals(:,p.nout)]);
axmax = max([plsOut{1}.Y(:,p.nout); mapped_vals(:,p.nout)]);

% loop through the PLS plots

for iPLS = 1:numel(plsOut)
    %   Cumulative Variance explained per component (2nd plot)
    vh = subplot(1,4,1,'Parent',botBot); 
    plot(vh, cumsum(plsOut{iPLS}.PCTVAR(2,:))*100,'-o');
    title(vh, '% Variance explained'); xlabel(vh, 'Component'); ylabel(vh, '%OPP Variance Explained'); 
    axis(vh, 'square'); 
    
    
    %   Parameter weights
    if iPLS < 2
        xlim(vh, [1,max(2,size(plsOut{iPLS}.PCTVAR,2))]);
        set(vh,'XTick',1:size(plsOut{iPLS}.PCTVAR,2)); ylim(vh, [0,100]);

        ph = subplot(1,4,2,'Parent',botBot); bar(ph, plsOut{iPLS}.BETA(2:end,p.nout));
        set(ph,'XTick',[1:size(plsOut{iPLS}.param,2)],'XTickLabel',plsOut{iPLS}.param,'XTickLabelRotation',45);
        xlim(ph, [0,size(plsOut{iPLS}.BETA,1)]);
        %   Set an overall title
        title(['PLS: ',plsOut{iPLS}.input(:)',', ' ,plsOut{iPLS}.output(:)']);

        hold(vh, 'on');  hold(ph, 'on');
    
        %Add significance thresholds from bootstrapping, if available
        if isfield(plsOut{iPLS},'Boot');    cc = [0.4,0.4,0.4];
            %Plot Variance Explained random region
            rry  = cumsum(plsOut{iPLS}.Boot.tpv,1)*100; ncmp = size(rry,1); %Y-Coords
            if ncmp == 1; rry = [rry;rry]; ncmp = 2; end
            rry = [rry(1:ncmp,1)', rry(ncmp:-1:1,2)'];           %  Wrapped
            rrx = [1:ncmp, ncmp:-1:1];                           %X-Coords
            %   Draw patch for random region
            patch(vh, rrx, rry, 'k', 'FaceAlpha', 0.2, 'LineStyle', 'none');
            
        %     ylp = ylim(p3h); %Store Parameter Y Limits            
            %Get average thresholds for Parameters
            thr{1} = mean(plsOut{iPLS}.Boot.tbeta([false, plsOut{iPLS}.Boot.cat],p.nout,:),1);
            thr{2} = mean(plsOut{iPLS}.Boot.tbeta([false, plsOut{iPLS}.Boot.dis],p.nout,:),1);
            thr{3} = mean(plsOut{iPLS}.Boot.tbeta([false, ~plsOut{iPLS}.Boot.cat & ~plsOut{iPLS}.Boot.dis],p.nout,:),1);
            %   Plot Parameter thresholds
            np = numel(plsOut{iPLS}.BETA)-1;
            
            %Set up patch definition for random region
            rry = nan(2,np);
            for s = 1:np
                %   Determine threshold for this parameter
                if plsOut{iPLS}.Boot.cat(s);       thrt = thr{1};
                elseif plsOut{iPLS}.Boot.dis(s);   thrt = thr{2};
                else;                    thrt = thr{3};
                end
                %   Assemble segment for this parameter
                rry(:,s) = thrt;
            end
            jmp = find(diff(rry(1,:))); %Get sites of changes
            jmpx = [jmp;jmp]; jmpx = jmpx(:)';   %Duplicate for X-Axis
            jmpy = [jmp;jmp+1]; jmpy = jmpy(:)'; %Include step for Y-Axis
            rry = rry(:,[1,jmpy,end]);   rry = [rry(1,:), rry(2,end:-1:1)];
            rrx = [-0.5,jmpx,np+0.5] + 0.5;     rrx = [rrx, rrx(end:-1:1)];
            %   Draw patch for random region
            patch(ph, rrx, rry, 'k', 'FaceAlpha', 0.2, 'LineStyle', 'none');
        end
    end
end % pls loop
lego = legend(legList);

% since PLC1 and PC2 explain most variance, show thier weights
sh = []; clear sh;
for s = 1:2
    sh(s) = subplot(1, 4, s+2,'Parent',botBot); 
    bar(sh(s), plsOut{1}.XLinv(:,s));  title(sh(s), ['PC', num2str(s),' weights']);
    set(sh(s),'XTick',[1:size(plsOut{1}.param,2)],'XTickLabel', plsOut{1}.param, 'XTickLabelRotation',45);
    xlim(sh(s), [0,size(plsOut{1}.XLinv,1)+1]);
end

vh.Units="inches"; ph.Units="inches"; sh(1).Units="inches"; sh(2).Units="inches";
fontname('Arial'); fontsize(8,"points");

vh.InnerPosition = [0.25,0.2,1.5,2.1];
ph.InnerPosition = [2.15,0.475,1.5,1.525];
sh(1).InnerPosition = [4.1,0.475,1.5,1.525];
sh(2).InnerPosition = [6,0.475,1.5,1.525];


% change the background colors
figure3.Color = [1,1,1];
topP.BackgroundColor = [1,1,1];
botLP.BackgroundColor = [1,1,1];
botRP.BackgroundColor = [1,1,1];
botBot.BackgroundColor = [1,1,1];

lego.FontSize = 3;

saveas(figure3,'Z:\imageData\SG_4B\Paper_Figures\Output_Figures\Figure_3.fig')
saveas(figure3,'Z:\imageData\SG_4B\Paper_Figures\Output_Figures\Figure_3.svg')


%% Now print the raw OPP values per cell line per NaAsO2 treatment

subz = ~contains(dataset1.ifd.cell,'bad','ignorecase',true); % filter for the parameters set above

orderOfPlot = {...
    '0ug/mL TET at hour 0 and 0uM NaAsO2', '0.1ug/mL TET at hour -2 and 0uM NaAsO2', '0.1ug/mL TET at hour -4 and 0uM NaAsO2','0.1ug/mL TET at hour -8 and 0uM NaAsO2', '0.1ug/mL TET at hour -14 and 0uM NaAsO2', '0.1ug/mL TET at hour -24 and 0uM NaAsO2',...
    '0ug/mL TET at hour 0 and 62.5uM NaAsO2', '0.1ug/mL TET at hour -2 and 62.5uM NaAsO2', '0.1ug/mL TET at hour -4 and 62.5uM NaAsO2','0.1ug/mL TET at hour -8 and 62.5uM NaAsO2', '0.1ug/mL TET at hour -14 and 62.5uM NaAsO2', '0.1ug/mL TET at hour -24 and 62.5uM NaAsO2',...
    '0ug/mL TET at hour 0 and 125uM NaAsO2', '0.1ug/mL TET at hour -2 and 125uM NaAsO2', '0.1ug/mL TET at hour -4 and 125uM NaAsO2','0.1ug/mL TET at hour -8 and 125uM NaAsO2', '0.1ug/mL TET at hour -14 and 125uM NaAsO2', '0.1ug/mL TET at hour -24 and 125uM NaAsO2',...
    '0ug/mL TET at hour 0 and 250uM NaAsO2', '0.1ug/mL TET at hour -2 and 250uM NaAsO2', '0.1ug/mL TET at hour -4 and 250uM NaAsO2','0.1ug/mL TET at hour -8 and 250uM NaAsO2', '0.1ug/mL TET at hour -14 and 250uM NaAsO2', '0.1ug/mL TET at hour -24 and 250uM NaAsO2',...
    };

allOPPData = dataset1.ifd(subz,:);
allOPPData.txcat = categorical(allOPPData.treatment,orderOfPlot);
allOPPData.OPP = allOPPData.Intensity_MeanIntensity_Masked_OPP*65535;


allOPPStats = grpstats(allOPPData,["txcat","cell"],["mean","median","sem","std"],"DataVars","OPP")

% save it as a csv 
writetable(allOPPStats,'Z:\imageData\SG_4B\Paper_Figures\Output_Figures\allOppStats2.csv')

%% Run the statistics comparing the wt-4B treated cells accross various concentrations of NaAsO2 and control

% get averything for the reader if needed
grpstats(wtSubData,"treatment",["mean","median","sem","std"],"DataVars",["NumGrans_rate_in_min","NumGrans_f","NumGrans_min_to_respond"])

% Now just print the means
grpstats(wtSubData,"treatment","mean","DataVars",["NumGrans_rate_in_min","NumGrans_f","NumGrans_min_to_respond"])


%% X hour tet wt 4b or 139a 4b vs opp @ each naaso2 dose


dataset1.ifd.txcat = categorical(dataset1.ifd.treatment,orderOfPlot);

subzwt = contains(dataset1.ifd.cell,'4BGFP'); % subset the data to be wt 4b
subzmut = contains(dataset1.ifd.cell,'139A'); % subset the data to be 139A 4b

barColz = [ 0.3010, 0.7450, 0.9330;...
            0.6350, 0.0780, 0.1840;...
            0.9290, 0.6940, 0.1250;...
            0, 0.4470, 0.7410;...
            0.4660, 0.6740, 0.1880;...
            0.4940, 0.1840, 0.5560]; % colors for bars

supF3 = figure;
subplot(2,1,1)
boxplot(log2(dataset1.ifd.Intensity_MeanIntensity_Masked_OPP(subzwt)*65535),dataset1.ifd.txcat(subzwt),'Notch','on','Symbol','.','Colors',barColz)
ylim([7, 15])
xlabel('TET Induction'); ylabel('Log_2 OPP Intensity')
title('WT4b Tet induction versus OPP')
legend

subplot(2,1,2)
boxplot(log2(dataset1.ifd.Intensity_MeanIntensity_Masked_OPP(subzmut)*65535),dataset1.ifd.txcat(subzmut),'Notch','on','Symbol','.','Colors',barColz)
xlabel('TET Induction'); ylabel('Log_2 OPP Intensity')
ylim([7,15])
title('139A Tet induction versus OPP')


saveas(supF3,'Z:\imageData\SG_4B\Paper_Figures\Output_Figures\Figure_S3.fig')
saveas(supF3,'Z:\imageData\SG_4B\Paper_Figures\Output_Figures\Figure_S3.svg')

%% Get the number of G3BP1 Granules @ 125uM NaAsO2

minGrans = 3;
tetTime = '-24';

% wt and mut data
subzBoth = all([contains(dataset1.ifd.treatment,['TET at hour ', tetTime]), ~contains(dataset1.ifd.cell,'bad')...
    , (dataset1.ifd.Children_Grans_G3BP1_Count >= minGrans), contains(dataset1.ifd.treatment,'and 125uM')],2); % 


% pull the mutant data, get the real opp and l2opp
allData = dataset1.ifd(subzBoth,:);

%simplify treatment name
allData.treatment = strrep(allData.treatment,'0.1ug/mL TET at hour -24 and ','');
allData.treatment = strrep(allData.treatment,' NaAsO2','');

%simplify cell names
allData.cell = strrep(allData.cell,"HeLa 4B139AGFP",'mut');
allData.cell = strrep(allData.cell,"HeLa eIF4BGFP",'wt');

% get the general values of the mut l2OPP
G3BP1data = grpstats(allData,["treatment","cell"],["mean","median","sem","std"],"DataVars","Children_Grans_G3BP1_Count")


% Do the statistics for l2 OPP versus treatment (WT-4B 125uM NaAsO2 as control)
[~,~,statsG3BP1] = anova1(allData.Children_Grans_G3BP1_Count,allData.cell,'off');
[resultsG3bp1,~,~,gnamesG3bp1] = multcompare(statsG3BP1,"CriticalValueType","dunnett",'ControlGroup',find(matches(statsG3BP1.gnames,'wt')),'Display','off','Approximate',false);
g3bp1 = array2table(resultsG3bp1,"VariableNames", ["Group","Control Group","Lower Limit","Difference","Upper Limit","P-value"]);
g3bp1.("Group") = gnamesG3bp1(g3bp1.("Group"));
g3bp1.("Control Group") = gnamesG3bp1(g3bp1.("Control Group"))


%
% Now plot it
figG3bp1 = figure;
boxplot(allData.Children_Grans_G3BP1_Count,allData.cell,'Notch','on','Symbol','')
xlabel('Cell'); ylabel('G3BP1 Grans')
title('G3BP1 Grans At 125 uM')
ylim([0,17])

saveas(figG3bp1,'Z:\imageData\SG_4B\Paper_Figures\Output_Figures\Figure_125uM_G3BP1.fig')
saveas(figG3bp1,'Z:\imageData\SG_4B\Paper_Figures\Output_Figures\Figure_125uM_G3BP1.svg')

