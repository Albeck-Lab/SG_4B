dataloc.ifd % and 125uM

%% Do 4B granules and g3bp1 granules occur at a similar frequency within the cell?
figure
subplot(2,2,1)
subz = all([contains(dataloc.ifd.treatment,'TET at hour -24 and 250uM'),contains(dataloc.ifd.cell,'4BGFP'), ...
    (dataloc.ifd.Children_Grans_4B_Count > 1),(dataloc.ifd.Children_Grans_G3BP1_Count > 1), ...
    (dataloc.ifd.Children_Grans_4B_Count < 35),(dataloc.ifd.Children_Grans_G3BP1_Count < 35), ...
    ],2); %~contains(dataloc.ifd.treatment,' 0uM')
mdl = fitlm(dataloc.ifd.Children_Grans_G3BP1_Count(subz),dataloc.ifd.Children_Grans_4B_Count(subz),'Intercept',false,'VarNames',{'G3BP1 Granules','4B Granules'}); %
plot(mdl)
title('24 hour TET induced wt-4B-GFP vs G3BP1 Granules in the cell')
xlim([0,40]); ylim([0,40])
subplot(2,2,2)
histogram2(dataloc.ifd.Children_Grans_G3BP1_Count(subz),dataloc.ifd.Children_Grans_4B_Count(subz),'DisplayStyle','tile','ShowEmptyBins','on')

subplot(2,2,3)
subz = all([contains(dataloc.ifd.treatment,'TET at hour -24 and 250uM'),contains(dataloc.ifd.cell,'139A'), ...
    (dataloc.ifd.Children_Grans_4B_Count > 1),(dataloc.ifd.Children_Grans_G3BP1_Count > 1), ...
    (dataloc.ifd.Children_Grans_4B_Count < 35),(dataloc.ifd.Children_Grans_G3BP1_Count < 35), ...
    ],2); %~contains(dataloc.ifd.treatment,' 0uM')
mdl = fitlm(dataloc.ifd.Children_Grans_G3BP1_Count(subz),dataloc.ifd.Children_Grans_4B_Count(subz),'Intercept',false,'VarNames',{'G3BP1 Granules','4B Granules'}); %
plot(mdl)
title('24 hour TET induced 139A-4B-GFP vs G3BP1 Granules in the cell')
xlim([0,40]); ylim([0,40])
subplot(2,2,4)
histogram2(dataloc.ifd.Children_Grans_G3BP1_Count(subz),dataloc.ifd.Children_Grans_4B_Count(subz),'DisplayStyle','tile','ShowEmptyBins','on')


%% 24 hour tet grans by OPP
figure
subplot(1,2,1)
subz = all([contains(dataloc.ifd.treatment,'TET at hour -24 and 125uM'), (dataloc.ifd.Children_Grans_4B_Count > 3),(dataloc.ifd.Children_Grans_G3BP1_Count > 1)],2);
gscatter(dataloc.ifd.Children_Grans_4B_Count(subz),dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subz)*65000,dataloc.ifd.cell(subz))
lsline
xlabel('4B Granules')
ylabel('Mean OPP Intensity')
title('24h Tet induced and treated with 125uM NaAsO2')

subplot(1,2,2)
subz = all(contains(dataloc.ifd.treatment,'TET at hour -24'),2);
gscatter(dataloc.ifd.Children_Grans_4B_Count(subz),dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subz),dataloc.ifd.cell(subz))
lsline
xlabel('4B Granules')
ylabel('Mean OPP Intensity')
title('24h Tet induced and all NaAsO2 Doses')




%% 24 hour tet (wt 4b or 139a) g3bp1 grans vs opp vs naaso2 dose
minGrans = 2;

figure
subz = all([contains(dataloc.ifd.treatment,'TET at hour -24'),contains(dataloc.ifd.cell,'4BGFP'),~contains(dataloc.ifd.treatment,' 0uM'),...
    (dataloc.ifd.Children_Grans_G3BP1_Count > minGrans),],2);
scatterhist(dataloc.ifd.Children_Grans_G3BP1_Count(subz),dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subz)*65000,...
    'Group',dataloc.ifd.treatment(subz),'Kernel','on','Location','SouthEast')
xlim([0,40]);
ylim([0,0.07*65000])
xlabel('G3BP1 Granules')
ylabel('Mean OPP Intensity')
title('24h WT4b Tet induced and all NaAsO2 Doses')

figure
subz = all([contains(dataloc.ifd.treatment,'TET at hour -24'),contains(dataloc.ifd.cell,'139A'),~contains(dataloc.ifd.treatment,' 0uM'),...
    (dataloc.ifd.Children_Grans_G3BP1_Count > minGrans)],2);
scatterhist(dataloc.ifd.Children_Grans_G3BP1_Count(subz),dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subz)*65000,...
    'Group',dataloc.ifd.treatment(subz),'Kernel','on')
xlim([0,40])
ylim([0,0.07*65000])
xlabel('G3BP1 Granules')
ylabel('Mean OPP Intensity')
title('24h 137A Tet induced and all NaAsO2 Doses')


%% 24 hour tet g3bp1 grans by OPP WT
figure
subplot(1,2,1)
subz = all([contains(dataloc.ifd.treatment,'TET at hour -24 and 62.5uM'),contains(dataloc.ifd.cell,'4BGFP'), (dataloc.ifd.Children_Grans_G3BP1_Count > 1)],2);
gscatter(dataloc.ifd.Children_Grans_G3BP1_Count(subz),dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subz),dataloc.ifd.cell(subz))
lsline
xlabel('G3BP1 Granules')
ylabel('Mean OPP Intensity')
title('24h Tet induced and treated with 62.5uM NaAsO2')

subplot(1,2,2)
subz = all([contains(dataloc.ifd.treatment,'TET at hour -24'),contains(dataloc.ifd.cell,'4BGFP'), (dataloc.ifd.Children_Grans_G3BP1_Count > 1)],2);
gscatter(dataloc.ifd.Children_Grans_G3BP1_Count(subz),dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subz),dataloc.ifd.cell(subz))
lsline
xlabel('G3BP1 Granules')
ylabel('Mean OPP Intensity')
title('24h Tet induced and all NaAsO2 Doses')

%% 24 hour tet g3bp1 grans by OPP
figure
subplot(1,2,1)
subz = all([contains(dataloc.ifd.treatment,'TET at hour -24 and 62.5uM'), (dataloc.ifd.Children_Grans_G3BP1_Count > 1)],2);
gscatter(dataloc.ifd.Children_Grans_G3BP1_Count(subz),dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subz),dataloc.ifd.cell(subz))
lsline
xlabel('G3BP1 Granules')
ylabel('Mean OPP Intensity')
title('24h Tet induced and treated with 62.5uM NaAsO2')

subplot(1,2,2)
subz = all([contains(dataloc.ifd.treatment,'TET at hour -24'), (dataloc.ifd.Children_Grans_G3BP1_Count > 1)],2);
gscatter(dataloc.ifd.Children_Grans_G3BP1_Count(subz),dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subz),dataloc.ifd.cell(subz))
lsline
xlabel('G3BP1 Granules')
ylabel('Mean OPP Intensity')
title('24h Tet induced and all NaAsO2 Doses')


%% OPP vs G3BP1 Grans by TET induction
figure
subz = all([contains(dataloc.ifd.treatment,'125uM'), (dataloc.ifd.Children_Grans_G3BP1_Count > 1)],2);
gscatter(dataloc.ifd.Children_Grans_G3BP1_Count(subz),dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subz)*65000,...
    {dataloc.ifd.treatment(subz),dataloc.ifd.cell(subz)}, 'kcbgmr','.*',6)
lsline
xlabel('G3BP1 Granules')
ylabel('Mean OPP Intensity')
title('OPP vs Grans by TET induction')

%% OPP vs G3BP1 Grans by TET induction
tetDose = ' 0uM';
minGrans = 0;
maxG = 10;
ylimz = 25000;

figure
subz = all([contains(dataloc.ifd.treatment,tetDose),contains(dataloc.ifd.cell,'4BGFP'), (dataloc.ifd.Children_Grans_G3BP1_Count >= minGrans)],2);
scatterhist(dataloc.ifd.Children_Grans_G3BP1_Count(subz),dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subz)*65000,...
    'Group',dataloc.ifd.full(subz),'Color','kcbgmr','Kernel','on','Direction','out','Location','SouthEast')
xlabel('G3BP1 Granules');ylabel('Mean OPP Intensity')
xlim([minGrans,maxG]);ylim([0,ylimz]);
title([tetDose,' NaAsO2 Treated OPP vs Grans by WT 4B by TET induction'])

figure
subz = all([contains(dataloc.ifd.treatment,tetDose),contains(dataloc.ifd.cell,'139A'), (dataloc.ifd.Children_Grans_G3BP1_Count >= minGrans)],2);
scatterhist(dataloc.ifd.Children_Grans_G3BP1_Count(subz),dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subz)*65000,...
    'Group',dataloc.ifd.full(subz),'Color','kcbgmr','Kernel','on','Direction','out')
xlabel('G3BP1 Granules');ylabel('Mean OPP Intensity')
xlim([minGrans,maxG]);ylim([0,ylimz]);
title([tetDose,' NaAsO2 Treated OPP vs Grans by 139A Mut 4B by TET induction'])


%% OPP vs G3BP1 Grans by TET induction boxplots
tetDose = ' 0uM';
minGrans = 0;
gLimz = [0, 5];
ylimz = [0, 20000];

figure
hp1 = uipanel('position',[0,   0, 0.5, 1]);
hp2 = uipanel('position',[0.5, 0, 0.5, 1]);

subz = all([contains(dataloc.ifd.treatment,tetDose),contains(dataloc.ifd.cell,'4BGFP'), (dataloc.ifd.Children_Grans_G3BP1_Count >= minGrans)],2);

h = scatterhist(dataloc.ifd.Children_Grans_G3BP1_Count(subz),dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subz)*65000,...
    'Group',dataloc.ifd.full(subz),'Color','kcbgmr','Location','SouthEast','Parent',hp1);
xlabel('G3BP1 Granules'); ylabel('Mean OPP Intensity')
title([tetDose,' NaAsO2 Treated OPP vs Grans by WT 4B by TET induction'])
hold on;
clr = get(h(1),'colororder');
boxplot(h(2),dataloc.ifd.Children_Grans_4B_Count(subz),dataloc.ifd.treatment(subz),'orientation','horizontal',...
     'label',{'','','','','',''},'color',clr);
boxplot(h(3),(dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subz)*65000),dataloc.ifd.treatment(subz),'orientation','horizontal',...
     'label', {'','','','','',''},'color',clr);
set(h(2:3),'XTickLabel','');
view(h(3),[270,90]);  % Rotate the Y plot
xlim(h(1),gLimz); ylim(h(1),ylimz);
xlim(h(2),gLimz); xlim(h(3),ylimz); % Sync axes
hold off;


subz = all([contains(dataloc.ifd.treatment,tetDose),contains(dataloc.ifd.cell,'139A'), (dataloc.ifd.Children_Grans_G3BP1_Count >= minGrans)],2);

h = scatterhist(dataloc.ifd.Children_Grans_G3BP1_Count(subz),dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subz)*65000,...
    'Group',dataloc.ifd.full(subz),'Color','kcbgmr','Location','SouthWest','Parent',hp2);
xlabel('G3BP1 Granules'); ylabel('Mean OPP Intensity')
title([tetDose,' NaAsO2 Treated OPP vs Grans by Mutant-4B by TET induction'])
hold on;
clr = get(h(1),'colororder');
boxplot(h(2),dataloc.ifd.Children_Grans_4B_Count(subz),dataloc.ifd.treatment(subz),'orientation','horizontal',...
     'label',{'','','','','',''},'color',clr);
boxplot(h(3),(dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subz)*65000),dataloc.ifd.treatment(subz),'orientation','horizontal',...
     'label', {'','','','','',''},'color',clr);
set(h(2:3),'XTickLabel','');
view(h(3),[270,90]);  % Rotate the Y plot
xlim(h(1),gLimz); ylim(h(1),ylimz);
xlim(h(2),gLimz); xlim(h(3),ylimz); % Sync axes
hold off;

%% Does OPP 
figure
subz = all([contains(dataloc.ifd.treatment,'24 hour'),contains(dataloc.ifd.cell,'139A'), (dataloc.ifd.Children_Grans_G3BP1_Count > 1)],2);
scatterhist(dataloc.ifd.Children_Grans_G3BP1_Count(subz),dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subz)*65000,...
    'Group',dataloc.ifd.full(subz),'Color','kcbgmr')
xlabel('G3BP1 Granules')
ylabel('Mean OPP Intensity')
title('OPP vs Grans by 139A Mut 4B TET induction')

%%


subplot(1,2,2)
subz = contains(dataloc.ifd.treatment,'TET at hour -24');
gscatter(dataloc.ifd.Children_Grans_4B_Count(subz),dataloc.ifd.Intensity_MeanIntensity_Masked_OPP(subz),dataloc.ifd.cell(subz))
lsline
xlabel('4B Granules')
ylabel('Mean OPP Intensity')
title('24h Tet induced and all NaAsO2 Doses')
%%

scatter(dataloc.ifd.Children_Grans_4B_Count(subz),dataloc.ifd.Children_Grans_G3BP1_Count(subz),'filled')
lsline
xlabel('4B Granules')
ylabel('G3BP1 Granules')

scatter3(dataloc.ifd(subz,:),{'Children_Grans_4B_Count','Children_Grans_G3BP1_Count'},'cell');
legend
xlabel('4B Granules')
ylabel('G3BP1 Granules')