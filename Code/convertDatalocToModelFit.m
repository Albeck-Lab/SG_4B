%% Convert Dataloc to a Dataframe by which every cell's data is fit using the following function
% [modelDataframe] = convertDatalocToModelFit(dataloc, channel, varargin)
% Note: All output data will be in terms of minutes (such as channel per min)
%
% dataloc - should be dataloc object(s) or a list of paths to dataloc files
% channel - cell array of channel(s) to be parameterized and put into the table (aka "dataframe")
%
% varargin options:
% combinereplicates - combine technical replicates from experiments before adding to the table?  (default = false)
%                     (ex: 3 xys from one experiment become one xy in the final table)?
% 
% aftertx - after which treatment should data be considered? (default = 1)
%           Note: (ptx DOES NOT COUNT FOR THIS)
% 
% fillnans - fill nans with appropriate filler values? default = false
%
% tstartaftertx - how many hours after treatment should data start to be considered? (default = 0 - Immediately after tx chosen in aftertx)
%           Note: Give this value in hours!
%
% p.tmaxaftertx - how many hours after treatment should be considered? (default = [] - consider all time after tx chosen in aftertx)
%           Note: Give this value in hours!
%
%
% Made by Nick DeCuzzi
% Albeck Lab - 2023

function [pulseDF, txDF] = convertDatalocToModelFit(datalocs, channel, varargin)
% Inputs for truncating data
p.combinereplicates = false; % combine technical replicates from experiments before adding to the table (ex: 3 xys from one experiment become one xy in the final table)?
p.aftertx = 1; % after which treatment should data be considered? (ptx DOES NOT COUNT FOR THIS
p.tstartaftertx = 0; % how many hours after treatment should data start to be considered? (default = immediately after (0))
p.tmaxaftertx = []; % how many hours after treatment should be considered? (default = [] (until end of movie))
p.mintracklength = 20; % how many tps the track needs to exist in btwn tStart and tEnd to be used at all (default = 20tps)
p.normalize = false; % zscore normalize each column?
p.tmaxback = 1; % how many hours to go backward from treatment (used for "pre-treatment" baselines)
% Inputs for selecting what data to pull
p.subset = [];
p.exclude = [];
p.fillnans = false; % replace nans with values?
p.responderdelta = 0.4; % how large should the delta be to be considered a responder to the treatment?
p.respondermaxtx = 1; % how many hours after treatment should be considered for responder delta?
p.looptime = 3; % what looptime should be assumed if you dont have one?
p.presetf = true; % when solving the equation, should the f be fixed at the highest grans the cell gets and truncated to that time?
p.meanexperiments = false; % mean all the experiment's data before fitting?
p.afterf = 0; % how many tps after f to give for model fitment?
% I'm missing something here, idk what yet...

% How the table is constructed
p.expdata = {'exp',    'full',  'cell', 'treatment',  'xy', 'cellid'}; % variable names for the table
p.exptypes = {'string','string','string', 'string',    'string','double'}; % variable types for the varNames in the table

% Fix - add underscores in names here (and remove the names from the ('') below and add the _
p.pulsepars =   {'f','td','ts','rate_in_min','min_to_respond','rsquared'}; % pulse metrics to pull from each channel

p.expar = {'Cell', 'pTx', 'Tx1', 'Tx2','Tx3'}; % expar for ct_maptx

nin = length(varargin); % Map input-value pairs
if rem(nin,2) ~= 0; warning('Additional inputs must be provided as option, value pairs');  end %Check for even number of add'l inputs
for s = 1:2:nin;   p.(lower(varargin{s})) = varargin{s+1};end %Splits pairs to a structure

p.pulsetypes = repelem({'double'},numel(p.pulsepars));
if any(matches(p.pulsepars,'livecell')); p.pulsetypes{matches(p.pulsepars,'livecell')} = 'cell'; end

% Assert things are in cells
if ~iscell(datalocs); datalocs = {datalocs}; end
if ~iscell(channel); channel = {channel}; end
if ~isempty(p.subset) && ~iscell(p.subset); p.subset = {p.subset}; end 
if ~isempty(p.exclude) && ~iscell(p.exclude); p.exclude = {p.exclude}; end 

addpath('\\albecklab.mcb.ucdavis.edu\data\Code\Cell Trace','\\albecklab.mcb.ucdavis.edu\data\Code\Nick')

% Make a PA Dataframe for each dataloc given and combine them
pulseDF = table;
txDF = table;

for iData = 1:numel(datalocs) 
    dHold = datalocs{iData};
    % Check if dataloc is a dataloc struct or is a path to a dataloc file
    if isstring(dHold) || ischar(dHold) 
        dHold = load(dHold);
        dHold = dHold.dataloc;
    end
    if isstruct(dHold)
        [dHold, txDFhold]= makeModelDF(dHold, channel, p); % make the pulse analysis dataframe
        
        % append the dHold dataframe to datalocDF
        pulseDF = [pulseDF; dHold];

        % append the treament mapping 
        if isempty(txDF); txDF = [txDF; txDFhold];
        else
            missingTxs = ~ismember(txDFhold.Properties.VariableNames,txDF.Properties.VariableNames);
            if any(missingTxs)
                for iTx = txDFhold.Properties.VariableNames(missingTxs)
                    txDF.(iTx{:}) = zeros(size(txDF,1),1);
                end
            end
            clear missingTxs;

            missingTxs = ~ismember(txDF.Properties.VariableNames,txDFhold.Properties.VariableNames);
            if any(missingTxs)
                for iTx = txDF.Properties.VariableNames(missingTxs)
                    txDFhold.(iTx{:}) = zeros(size(txDFhold,1),1);
                end
            end
            clear missingTxs;

            txDF = [txDF; txDFhold];

        end 
    end
    clear dHold txDFhold;
end 

% add a column giving every cell a unique identifier for use to reference in plots
pulseDF.cellID = (1:size(pulseDF,1))';
pulseDF = movevars(pulseDF,{'cellID'},'After','tx_num');

[~,I] = sortrows(txDF);
txDF = txDF(I,:);
pulseDF = pulseDF(I,:);


end % [pulseDataframe] = convertPulseToDataframe(dataloc, channel, varargin)

%% Actual function for making each table or "dataframe"
function [pdDF, txxDF] = makeModelDF(dataloc,channel,p)
pdDF = table; txxDF = table;
if isfield(dataloc,'movieinfo'); p.timePerLoop = unique(dataloc.movieinfo.tsamp); % get the mins btwn loops
else; p.timePerLoop = p.looptime;
end
p.tpPerHr = 60/p.timePerLoop; %set up how many tps per hour

% FIX, maybe make this floor? or ceiling? for weird times per loop (CUA, IM LOOKNIG AT YOU , MR. 7 minutes per loop!!!)
p.tstartaftertx = p.tstartaftertx * p.tpPerHr; % convert tstartafter tx from hours into tps
if ~isempty(p.tmaxaftertx); p.tmaxaftertx = p.tmaxaftertx * p.tpPerHr; end % convert respond time from hours to tps

%% Map all the data in this experiment to their treatments 
mappedTxs = ct_maptx(dataloc.platemapd.idx,'expar',p.expar,'reduce',false); % map the platemap
txNames = fieldnames(mappedTxs);% get the unique datasets to either combine or not

%% Subset and exlude things

% get the tx that DO have these things...
subI = []; % set up the list of things to keep
if ~isempty(p.subset)
    for iSubset = p.subset
        subI2 = contains(txNames,iSubset{:},'IgnoreCase',true);
        subI =[subI, subI2];
    end
    keepMe = any(subI,2);
    if any(keepMe)
        txNames = txNames(keepMe);
    end
end

% get the txs that DONT have these things...
dexI = []; % set up the list of things to not exclude
if ~isempty(p.exclude)
    for iExclude = p.exclude
        dexI2 = ~contains(txNames,iExclude{:},'IgnoreCase',true);
        dexI = [dexI, dexI2];
    end
    keepMe = all(dexI,2);
    if any(keepMe)
        txNames = txNames(keepMe);
    end
end

if isempty(txNames); return; end % if there is nothing left, return
%% Set up info for the dataframe
 %('Size',[0,numel(p.varNames)],'VariableNames',p.varNames,'VariableTypes', p.varTypes); % Set up the dataframe


%% Figure out which xys have data we (used later) 
xyList = 1:size(dataloc.d,2); % set up xyList
xyList2 = xyList; % make a copy for empty xys
hasD = arrayfun(@(x) ~isempty(dataloc.d{x}), xyList); % get the xys for this dataset and check that the d is not empty
xyList = xyList(hasD); % pull non empty xys
hasChans = arrayfun(@(xx)all([cellfun(@(x)isfield(dataloc.d{xx}.data,x),channel)],2),xyList)'; % data needs all channels youre trying to pull
if any(~hasChans); fprintf('These xys dont have data and/or the channels you asked for: %s', num2str(xyList(~hasChans))); end % tell you about empty xys
xyList = xyList(hasChans); % pull the data with data and chans we want

% Make sure you have live cell data too!
hasChans = arrayfun(@(xx)all([cellfun(@(x)isfield(dataloc.d{xx}.data,x),channel)],2),xyList)'; %check for live cell data to go with the pulse data
xyList = xyList(hasChans); % pull the data with data and chans we want

emptyXYs = ~ismember(xyList2',xyList'); %find xys that were empty
emptyXYs = xyList2(emptyXYs)'; % get those xys

% give warning about empty xys if there are any.
if ~isempty(emptyXYs)
    emptyXYs = num2str(emptyXYs); % make it a string
    % if emptyXYs > 1; emptyXYs = strcat(emptyXYs, {', '}); end %combine
    % the empty xys FIX
    warning('The following xys from %s are empty: %s.',dataloc.file.base, emptyXYs);
end

if isempty(xyList); return; end % if no xys have your data, get out of here.

% Establish the movie length
movieLength = size(dataloc.d{xyList(1)}.data.(channel{1}),2); % Establish entire movie length in tps

%% Make the treatment dataframe to be filled as each treatment is added
allTxs = cell2mat(cellfun(@(x){mappedTxs.(x).tx},txNames)); % get all the txs that exist in the treatments
notCells = ~contains({allTxs.dunit},'cells')'; % find the rows that dont have cell line info (FIX - for now?)
allTxs = allTxs(notCells); % pull the treatment data
allTxs = struct2table(allTxs);
allTxs = allTxs(~cellfun(@isempty, allTxs.name),:);
if ~isa(allTxs.time,'double')
    allTxs.time = cell2mat(allTxs{:,'time'});
end
allTxs = unique(allTxs,'rows'); % pull the unqiue possible treatments
allTxs = sortrows(allTxs,["name","dose"]); % FIX? sort them by treatment and dose? idk

uniqueTxsWithoutDose = unique(allTxs(:,[1,3:5]),'rows','stable'); % find unique treatments considering timing, but ignoring dose
uniqueTxsWithoutDose = uniqueTxsWithoutDose(~ismissing(uniqueTxsWithoutDose.tunit),:); % remove any empty fields
% txTableFieldNames = {}; % set up the table names for txxDF
% for iTx = 1:size(uniqueTxsWithoutDose,1)
%     if uniqueTxsWithoutDose{iTx,"tunit"} == "h"; txTableFieldNames(iTx) = append(uniqueTxsWithoutDose{iTx,"dunit"}, ' ', uniqueTxsWithoutDose{iTx,"name"},' at hour ', num2str(uniqueTxsWithoutDose{iTx,"time"}));
%     elseif uniqueTxsWithoutDose{iTx,"tunit"} == "tp"; txTableFieldNames(iTx) = append(uniqueTxsWithoutDose{iTx,"dunit"}, ' ', uniqueTxsWithoutDose{iTx,"name"},' at tp ', num2str(uniqueTxsWithoutDose{iTx,"time"}));
%     end
% end
% txTableFieldNames = strtrim(txTableFieldNames);  % clean them up
% txxDF = table('Size',[1,numel(txTableFieldNames)],'VariableNames',txTableFieldNames,'VariableTypes',repelem({'double'},numel(txTableFieldNames)));

% create the txxDF with one row that needs to be removed at the end!
clear uniqueTxsWithoutDose allTxs notCells;

%% loop through every uniquely mapped cell/condition combo (keeping them together in the table makes it easier to check if done properly)
for iMap = 1:size(txNames,1)

    % Make a temp dataframe for this condition (each channel's data will be appended to this, then this df will be appended to the final dataframe)
    txTempDF = table;

    txXys = intersect(mappedTxs.(txNames{iMap}).xy,xyList'); % get the xys for this mapped info and make sure they have data

    if isempty(txXys); continue; end % skip this tx if there is no data

    hasTxs = contains({mappedTxs.(txNames{iMap}).tx(:).tunit},'tp')';
    theseTxTps = unique(cell2mat({mappedTxs.(txNames{iMap}).tx(hasTxs).time}')); % get the treatment tps (that happen during the movie)

    % figure out what the first tp is for this treatment set or make it tp 1
    if ~isempty(p.aftertx)
        if numel(theseTxTps) < p.aftertx; tStart=theseTxTps(end); else; tStart=theseTxTps(p.aftertx); end % get treatment number or last treatment
        txTP = tStart;
    
    else; tStart = 1; txTP = 1;
    end

    % establish tStart and tEnd from aftertx, tstartaftertx, and tmaxaftertx
    tStart = tStart + p.tstartaftertx;
    if tStart > movieLength; tStart = movieLength-1; elseif tStart < 1; tStart = 1; end % check tStart is possible

    if isempty(p.tmaxaftertx); tEnd = movieLength; else; tEnd = p.tmaxaftertx+tStart; end % check tEnd is possible
    if tEnd > movieLength; tEnd = movieLength; end
    
    %% now go through each channel and pull everything from combined datasets
    for oniChan = channel
        % make a table for the pulse data
        tempLC = cell2mat(arrayfun(@(xx){dataloc.d{xx}.data.(oniChan{:})},txXys)); %pull the live cell data
        tempCellIDs = cell2mat(arrayfun(@(xx){dataloc.d{xx}.cellindex},txXys)); %pull the cell ids
        tempXYs = [];
        for iXY = txXys'
            tempXYs = [tempXYs; repelem(iXY,size(dataloc.d{iXY}.cellindex,1))']; %pull the xys the data came from
        end

        %% Truncate the Live cell data
        tempLC = tempLC(:,(tStart:tEnd)); % Pulls the live cell data from that time window   

        % truncate the data so that only the data thats long enough gets
        % through.
        tklngth = sum(~isnan(tempLC),2);
        goodData = (tklngth >= p.mintracklength); % good to set this min to 1 hr
        tempLC = tempLC(goodData,:);
        tempCellIDs = tempCellIDs(goodData);
        tempXYs= tempXYs(goodData);

        if ismember('granarea',p.pulsepars)
            tempGranArea = cell2mat(arrayfun(@(xx){dataloc.d{xx}.data.('GranAreaShape_Area')},txXys)); %pull the live cell data
            tempGranArea = tempGranArea(goodData,:);
            
            % convert the area from pixels to area
            tempGranArea = tempGranArea * (dataloc.movieinfo.PixSizeY * dataloc.movieinfo.PixSizeX);
        end

        %% Fit the live cell data for each cell

        avgData = nanmean(tempLC,1); % use averages for approximations
        maxY = max(avgData); % get max of current data
        minY = min(avgData); % get min of current data
        [~,Td0] = find(diff((tempLC >= maxY),1,2),1); % approximate Td0 by finding first time y crosses the half max (going either up or down)
        if isempty(Td0); Td0 = size(tempLC,2)/2; end % give a value if nan
        [~,crossNearMax] = find(diff((avgData > maxY*0.90),1,2),1); % approximate crossing values
        if isempty(crossNearMax); crossNearMax = 1; end % just in case 
        [~,crossNearMin] = find(diff((avgData > (minY + (maxY-minY)*0.1)),1,2),1);
        if isempty(crossNearMin); crossNearMin = 1; end % just in case 
        Ts0 = abs(crossNearMax - crossNearMin) ; % approximate Ts0 by finding when y crosses 10% of the max (from the min) and when y crosses 90% of the max and getting that time span
        if isempty(Ts0); Ts0 = Td0/2; end % guess
        
        Td0 = Td0 * p.timePerLoop; % convert the value to minutes
        Ts0 = Ts0 * p.timePerLoop; % convert the value to minutes
        f0 = maxY; % give an approximate f

        xs = (tStart-tStart):(tEnd-tStart); % get what data you could consider
        xs = xs*p.timePerLoop; % set the data to minutes

        % fits an equation based on Albeck et al (2008) Modeling a Snap-Action, Variable-Delay Switch Controlling Extrinsic Cell Death (PLoS)
        % f is max signal, t is time, Td is time to half max signal, Ts is time of activation (time from signal starts to go up until it reaches max), 
        % original equation: c(t) = f-(f/(1+e^((t-Td)/(Ts/4))))
        funTime = 'f-f/(1+exp(((xdata-Td)/(Ts/4))))';
        optz = fitoptions('Method','NonlinearLeastSquares', 'MaxIter',50000,'MaxFunEvals',50000);

        % if you want to force f to be the max value of grans the cell has
        if p.presetf
            z0 = [Td0,Ts0]; % give the initial guesses to start with  
            optz.StartPoint = z0;
            optz.Upper = [max(xs)*1.2,max(xs)*1.2]; % set the upper bounds to a real number
            optz.Lower = [1,1]; % set the lower bounds to a real number
            coeffs = {'Td','Ts'};

            % set up the equation, but allow f to be fixed
            F = fittype(funTime,'dependent',{'y'},'independent',{'xdata'},'coefficients',coeffs,'options',optz,'problem','f'); 

        else % otherwise determine it via fitting
            z0 = [f0,Td0,Ts0]; % give the initial guesses to start with                                        
            optz.StartPoint = z0;
            optz.Upper = [(max(tempLC,[],'all')*1.2),max(xs)*1.2,max(xs)*1.2]; % set the upper bounds to a real number
            optz.Lower = [0,1,1]; % set the lower bounds to a real number
            coeffs = {'f','Td','Ts'}; % set up the equation coeffs

            % set up the equation, but allow f to be solved for
            F = fittype(funTime,'dependent',{'y'},'independent',{'xdata'},'coefficients',coeffs,'options',optz); 

        end
        
        %% assign the data we have to the dataframe

        % mean that experiment's data before fitting it?
        if p.meanexperiments
            tempLC = mean(tempLC,1,"omitnan");
            tempXYs = iMap; % the data is being meaned
            tempCellIDs = iMap;% the data is being meaned
        end
        tempDF = table('Size',[size(tempLC,1),numel(p.pulsepars)],'VariableNames',p.pulsepars,'VariableTypes', p.pulsetypes); % Generates empty table of parameters
        tempDF.xy = tempXYs; % assign the xys
        tempDF.cellid = tempCellIDs; % assign the cell IDs

        %% Now do a per cell fitting

        for iCell = 1:size(tempLC,1)
            ys = tempLC(iCell,:);
            goodTps = ~isnan(ys); % pull only real data
            tempX = xs(goodTps); tempY = ys(goodTps);
            if p.presetf % if f is preset, truncate the data then fit giving the F as max grans
                maxF = find(tempY == max(tempY),1); % find where the F reaches max, but exclude the first tp
                if isempty(maxF) || tempY(maxF) == 0|| maxF == 1; maxF = size(tempY,2); % if no max give max amount of tps allowed
                end
                if (maxF+p.afterf) <= (size(tempY,2)-p.afterf); maxF = maxF+p.afterf; end
                [modelOut,gof] = fit(tempX(1:maxF)',tempY(1:maxF)',F,'problem',tempY(maxF)); % fit the model
                clear maxF;
            else; [modelOut,gof] = fit(tempX',tempY',F); % fit the model
            end

            tempDF.f(iCell) = modelOut.f;
            tempDF.td(iCell) = modelOut.Td;
            tempDF.ts(iCell) = modelOut.Ts;
            tempDF.td(iCell) = modelOut.Td;
            tempDF.rate_in_min(iCell) = (modelOut.f/modelOut.Ts); % solve for rate (grans per min)
            tempDF.min_to_respond(iCell) = modelOut.Td-(modelOut.Ts/2); % solve for time to respond (in minutes)
            tempDF.rsquared(iCell) = gof.adjrsquare;

            if ismember('livecell',p.pulsepars) % pull the live cell data if requested and put it in a cell
               tempDF.livecell(iCell) = {ys};
            end

            if ismember('granarea',p.pulsepars)
                maxF = ceil(((modelOut.Td+(modelOut.Ts/2))/p.timePerLoop)); % tp of the max granule
                if maxF > size(tempY,2); maxF = size(tempY,2); end % if f isn't reached before the end of time then get the last tp
                tempDF.granarea(iCell) = tempGranArea(iCell,maxF);
                clear maxF;
            end

            clear modelOut gof ys goodTps tempX tempY
            
        end
        
%         whereMax = find(ys == max(ys),1); % find the point where the number of granules reaches max
%         xs = xs(1:whereMax); % truncate data to span until max
%         ys = ys(1:whereMax)
        
        tempDF = movevars(tempDF,{'cellid','xy'},'Before',1); % move the data information to the first columns of the dataframe

        % Append the channel name to every parameter in the table
        newVarNames = strcat(oniChan,'_',tempDF.Properties.VariableNames);  % append channel name to variable names
        tempDF.Properties.VariableNames = newVarNames;  % assign the new variable names to the table

        % Append this tempDF to the condition's datatable
        txTempDF = [txTempDF, tempDF];

        clear tempDF;

    end % Channel loop

    % add the data source (the file the data came from), and the treatments the cells in this condition got (FIX: maybe add xys? idk)
    numData = size(txTempDF,1); % how many rows are in this dataset

    txTempDF.full = repelem(string(txNames{iMap}),numData)'; % add the full sorted name
    holdThis = '';

    tTxList = mappedTxs.(txNames{iMap}).tx;
    notCells = ~contains({tTxList.dunit},'cells')'; % find the rows that dont have cell line info (FIX - for now?)
    if any(~notCells)
        txTempDF.cell = repelem(string(mappedTxs.(txNames{iMap}).tx(~notCells).name),numData)';
    else; txTempDF.cell = repelem(" ",numData)'; % if no cell line fill it with nothing
    end
    
    tTxList = tTxList(notCells); % pull the treatment data
    for iTx = 1:size(tTxList,1)
        if tTxList(iTx).tunit == "h"; holdThis = append(holdThis, append(num2str(tTxList(iTx).dose), tTxList(iTx).dunit, ' ', tTxList(iTx).name,' at hour ', num2str(tTxList(iTx).time)),' and ');
        elseif tTxList(iTx).tunit == "tp" 
            if tTxList(iTx).time < tEnd
                holdThis = append(holdThis, append(num2str(tTxList(iTx).dose), tTxList(iTx).dunit, ' ', tTxList(iTx).name,' at hour ', num2str(round(((tTxList(iTx).time-tStart)/p.tpPerHr),0)),' and ')); % round to the nearest hour (.5 >= rounds up to 1, <0.5 rounds down to 0)
            end
        end
    end
    holdThis = char(holdThis);
    holdThis = string(holdThis(1:end-5));
    holdThis = strrep(holdThis,'_',' ');
    txTempDF.treatment = repelem(holdThis,numData)';
    txTempDF.tx_num = repelem(iMap,numData)';
    txTempDF = movevars(txTempDF,{'full','cell','treatment','tx_num'},'Before',1); % move the data information to the first columns of the dataframe
    pdDF = [pdDF; txTempDF]; % append the temp dataframe to the main one
    
    % make the teatment dose temp dataframe
    temptxxDF = table('Size',[numData,numel(txxDF.Properties.VariableNames)],'VariableNames',txxDF.Properties.VariableNames,'VariableTypes',repelem({'double'},numel(txxDF.Properties.VariableNames)));
    %loop through the txs and add their dose to the proper line
    for iTx = 1:size(tTxList,1)
        if tTxList(iTx).tunit == "h"; txfieldHold = append(tTxList(iTx).dunit, ' ', tTxList(iTx).name,' at hour ', num2str(tTxList(iTx).time));
        elseif tTxList(iTx).tunit == "tp"
            if tTxList(iTx).time > tEnd
            elseif tTxList(iTx).time < tStart; txfieldHold = append(tTxList(iTx).dunit, ' ', tTxList(iTx).name);
            else; txfieldHold = append(tTxList(iTx).dunit, ' ', tTxList(iTx).name,' at tp ', num2str(tTxList(iTx).time));
            end
        end
        if exist("txfieldHold",'var')
            txfieldHold = strtrim(txfieldHold);
            temptxxDF.(txfieldHold) = repelem(tTxList(iTx).dose,numData)'; % put the dose into its proper column
            clear txfieldHold;
        end
    end
    % append the temp treatment table to the real treatment table
    if isempty(txxDF); txxDF = [txxDF; temptxxDF];
    else; missingTxs = ~ismember(temptxxDF.Properties.VariableNames,txxDF.Properties.VariableNames);
        if any(missingTxs)
            for iTx = temptxxDF.Properties.VariableNames(missingTxs)
                txxDF.(iTx{:}) = zeros(size(txxDF,1),1);
            end
        end
        txxDF = [txxDF; temptxxDF];
    end 


    clear txTempDF tTxList holdThis numData temptxxDF;
end % txNames loop

% remove the first data point from the txxDF as it is fake FIX
%txxDF = txxDF(2:end,:);

% remove any columns where none all the cells got the same thing
badCols = false([1,size(txxDF,2)]);
for iTxCol = 1:size(txxDF,2)
    uniqueConcs = unique(txxDF{:,iTxCol});
    if uniqueConcs < 2
        badCols(iTxCol) = true;
    end
end

if any(badCols)
    txxDF = txxDF(:,~badCols);
end
%Mean experiment Combine experimental replicates 

%% z score normalize the dataset within itself if requested
if p.normalize
    % Normalize all columns except the full,cell,treatment,txnum, and any live cell data
    colNames = pdDF.Properties.VariableNames; % find the columns that aren't the ones above
    colNames = colNames(~contains(colNames,["full","cell","treatment","tx_num", "LiveCell"]));

    for i = 1:size(colNames,2)
    
         % Get the data for this column
        data = pdDF{:,colNames{i}};

        % get the indicies of non-nan data
        notTrash = ~isnan(data);

        % pull the not trash data
        goodData = data(notTrash);

        % Normalize the data to have zero mean and unit variance
        normalizedData = zscore(goodData);

        % return the normalized data to the places it came from?
        data(notTrash) = normalizedData;
    
        % Write the normalized data back into the table
        pdDF{:,colNames{i}} = data;
    
        clear data normalizedData notTrash goodData;
    end
    clear colNames;
end


end %end of makePADF(dataloc, channel, p) subfunction                  
                      
                            
                        