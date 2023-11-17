% dataloc = SG_Datahandler(varagin)
% SG_Datahandler was designed to add accountabilty to where data comes from
% and how it is handled. 
%
% You dont need to put anything in varargin if you just want the data in
% the folder you're in. If you want a different folder's data make sure to
% set difffolder.
% 
% The output of all the data SG_Datahandler puts together is saved as one
% variable in the processed file's folder with the name ending with _Processed.mat
%
% INPUTS:
% inputs are given as pairs into SG_Datahandler('like', like, 'this', this)
% The first part of the pair in the ''s matches one of the inputs below,
% The second part of the pair is the information you want to pair with it.
%
% Setting up data / image processing:
% saveit - if you want to save the output of SG_Datahandler (default = true)
% difffolder - a string containing a different folder to get the dataloc data from.
% savename - if you want to save the processed file under a different name
% (but in the same folder)
%
% Raw data processing:
% filterp - filter parameters to use to filter the data 
% filerp.dlength, .ft.t ft.c and ft.p are expected
% EX: filterp.fts(1).t = 'minmax';  filterp.fts(1).c = 'Med_RFP_Nuc';  filterp.fts(1).p = 800;
% filterp.dlength = 261;
% filters avaliable - 'Max', 'Min', 'dMax', 'dMin', 'tWin', 'mindyn', 'minmax', 'anymax', 'custom', 'areafilt','25prctlmin'
%
% If you pass parameters in for filterp then filter will switch to true
% and it will filer the output from the cellprofiler pipeline.
%
% certainxys - for processing only certain xys that are given (default = all xys)
%
%
% Data normalization:
% normalizedata - cell array of live movie data channels to make normalized versions of (for Live Cell)
% ex: 'normalizedata', 'AMPKAR', will produce a normalized version of
% the AMPKAR called 'nAMPKAR' using the normalization method selected below
%
% normalizetype - type of normalization (minmax (default), minsub, znorm, robust, none) 
% normto - specific xys to use for data normalization (defaut = use all xys)
%
% Plotting functions:
% updatemap - update the platemap (default = false)
% Note: This is required for all plot_by_ND functions to work, if you save
% changes to your platemap this will update automatically.
% Note2: if not done before (or the platemap is new), this will set true automatically
%
%  Copyright.
%% N. DeCuzzi, Albeck Lab, UC Davis. 2019-2023

function dataloc = SG_Datahandler(varargin)
%set up dataproc parameters
inp.sensors = ''; inp.filterp = []; inp.leavegaps = false; inp.gapmax = 2;

inp.fixfields = false; % make sure all of the data fields match
inp.saveit = true; %if you want to save the output
%set up (re)processing parameters
inp.extractdata = false; inp.updatemap = false; inp.filter = false; 
inp.certainxys = []; %do this only to certainxys
inp.remakegranspercell = false; % remake grans per cell
inp.extractif = false; % extract IF data and make IF dataframe
inp.normalizedata = {}; inp.normalizetype = 'minsub'; inp.normto = []; inp.normmax = 98; inp.normmin = 2;
inp.ifnorm = struct;
inp.datafields = {'Intensity_IntegratedIntensity_GFP', 'Intensity_MeanIntensity_GFP', 'AreaShape_Area', 'Location_Center_X', 'Location_Center_Y','Mean_Grans_AreaShape_Area','Mean_Grans_Intensity_MaxIntensityEdge_GFP','Mean_Grans_Intensity_MaxIntensity_GFP','Mean_Grans_Intensity_MeanIntensityEdge_GFP','Mean_Grans_Intensity_MeanIntensity_GFP','Mean_Grans_Intensity_IntegratedIntensity_GFP'};
inp.ifdatafields = {'Intensity_IntegratedIntensity_GFP', 'Intensity_MeanIntensity_GFP', 'AreaShape_Area', 'Location_Center_X', 'Location_Center_Y','Mean_Grans_AreaShape_Area','Mean_Grans_Intensity_MaxIntensityEdge_GFP','Mean_Grans_Intensity_MaxIntensity_GFP','Mean_Grans_Intensity_MeanIntensityEdge_GFP','Mean_Grans_Intensity_MeanIntensity_GFP','Mean_Grans_Intensity_IntegratedIntensity_GFP'};
inp.subset = [];
inp.looptime = 3; % how many min per loop

% set up dataloc
inp.dataloc = [];
dataloc = struct('fold',  struct('data','','proc','','fig','','if',''),... %Initialize the folders
    'file', struct('data','', 'base','', 'proc','', 'platemap',''),... %Initialize the files
    'procdate','', 'filterp',[], 'normalizetype', [], 'version','1', 'd', [], 'platemapd',struct('pmd',[],'idx',[],'date','')); %Initialize everything else
dataloc.sensorinfo = {}; dataloc.xys = {}; dataloc.appendedchans = {}; 

inp.difffolder = ''; %set this to get data fom a different folder
inp.savename = ''; %set this to save dataloc as a different file name

%% Check input varargin parameters
nin = length(varargin);     %Check for even number of add'l inputs
if rem(nin,2) ~= 0; warning(['Additional inputs must be provided as ',...
        'option, value pairs']);  end%Splits pairs to a structure
for s = 1:2:nin;   inp.(lower(varargin{s})) = varargin{s+1};   end; clear nin s varargin;

%% Determine what called SG_Datahandler and what files exist already
[dataloc, inp] = FolderFinderMaker(dataloc, inp);

saveIt = false; saveString = '';

%% See if XYs have been converted into time series and saved
if inp.extractdata || isempty(dataloc.d)
    [dataloc, inp] = ExtractData(dataloc, inp);
    saveString = [saveString 'reorganized, '];
    save(fullfile(dataloc.fold.data,dataloc.file.proc), 'dataloc')
end
% See if data proc perameters exist or have changed
if ~isempty(inp.filterp) 
   dataloc.filterp = inp.filterp;
   inp.filter = true;
end


%% See if the data has been processed
if inp.filter
    fprintf('Data processing has been requested. \n');       
     
    %Proc the files
    fprintf('Running data through ct_dataproc. \n');
    dataloc = ProcMe(dataloc,inp); 
    inp.remakegranspercell = true;

    % Change procdate and save it bro
    dataloc.procdate = datetime;
    saveString = [saveString 'filtered, '];
    saveIt = true; 
end
% filter the data from cellprofiler


%% Normalize and save data
if ~isempty(inp.normalizedata)
    dataloc = Normalize_Data(dataloc, inp);
    saveString = [saveString 'normalized, '];
    saveIt = true; 
end %normalize data

%% Pull IF and save data
if inp.extractif
    if isfield(dataloc.fold,'if') && ~isempty(dataloc.fold.if) && exist([dataloc.fold.if,'\','Cells.csv'],"file")
        dataloc = makeIFDF(dataloc, inp);
        saveIt = true;
        saveString = [saveString, 'IF,' ];
    end
end

%% Remake granspercell and save data
if inp.remakegranspercell
    dataloc = remakeGransPerCell_Data(dataloc, inp);
    saveString = [saveString 'remade granspercell, '];
    saveIt = true; 
end %remake granspercell data

%% Check for a save call
if inp.saveit && saveIt
    fprintf(['Saving the updated ', saveString, 'data. \n']);
    save(fullfile(dataloc.fold.data, dataloc.file.proc), 'dataloc')
end
end %JB_DatalocHandler_V2 End

%% All the functions DatalocHanlder will call
%% Proc function 
function dataloc = ProcMe(dataloc, inp)

if ~isempty(inp.certainxys)
    goodXYs = ~arrayfun(@(x)isempty(dataloc.d{x}),inp.certainxys);
    XYnums = inp.certainxys(goodXYs);
    NumXYs = numel(XYnums);
else
    %get the xy names from the data
    XYnums = find(~cellfun(@(x)isempty(x),dataloc.d));
    NumXYs = numel(XYnums);
end

PercentBar = floor(NumXYs/4); PercentBar = [1:PercentBar:PercentBar*4];
Quarters = 0;

for ii = XYnums
    if any(ii == PercentBar)
       fprintf('Data Processed: %d/4th done \n', Quarters); 
       Quarters = Quarters +1;
    end
    
    if ~isempty(dataloc.d{ii}) && isfield(dataloc.d{ii},'data')
        if isfield(dataloc.d{ii}.data,'t_granspercell')
            dataloc.d{ii}.data = rmfield(dataloc.d{ii}.data,{'t_granspercell','t_count_cells'});
        end
        dataloc.d(ii) = ct_datafilt_ND(dataloc.d(ii), dataloc.filterp.fts,'LeaveGaps',inp.leavegaps);
    else; dataloc.d{ii} = [];
    end %fix for empty xys
end %numxys 

fprintf('SG_Datahandler has processed %s. \n', dataloc.file.proc)
end

%% Normalize Data Function
function dataloc = Normalize_Data(dataloc, inp)
% dataloc = Normalize_Data(dataloc, inp)
% Takes in dataloc structure and returns the data normalized to either all of the data or using specific control xys.
% 'normto' = []; % what xys to normalize the data to (default == all xys)
% 'certainxys' = []; % normalize only certain xys (default == all xys)

%process either certain xys or all of them
if ~isempty(inp.certainxys)
    goodXYs = ~arrayfun(@(x)isempty(dataloc.d{x}),inp.certainxys);
    XYnums = inp.certainxys(goodXYs);
else
    %get the xy names from the data
    XYnums = extractAfter(dataloc.xys,'XY '); XYnums = str2double(XYnums)';
    NumXYs = numel(XYnums);
end

% How many Xys to use for the normalization
if ~isempty(inp.normto)
    NumXYs2Norm2 = numel(inp.normto);
else
    XYnums = extractAfter(dataloc.xys,'XY '); XYnums = str2double(XYnums)';
    NumXYs = numel(XYnums);
end

if ischar(inp.normalizedata)
    inp.normalizedata = {inp.normalizedata};
end

ChanNames = inp.normalizedata;
NumChans = numel(ChanNames);

%loop over all (or some) of the channels
for iChan = 1:NumChans
    CurChan = ChanNames{iChan}; % get the channel name
    fprintf('The %s channel will be normalized using %s method. \n', CurChan, inp.normalizetype)
    
    NormChan = ['n', CurChan]; % get the normalization channel name
    dataloc.normalizetype.(NormChan) = inp.normalizetype;
    MegaDataHolder = []; %This will hold all the data for a second
    
    %loop over all the xys and combine them into a mega data
    for iXY = 1:NumXYs2Norm2
        tXY = inp.normto(iXY);
        
        if ~isempty(dataloc.d{tXY}) %check the XY isn't empty
        if isfield(dataloc.d{tXY}, 'data') %check the XY has data
        if isfield(dataloc.d{tXY}.data, CurChan) %check if that XY has the channel in it
            MegaDataHolder = [MegaDataHolder; dataloc.d{tXY}.data.(CurChan)]; %append all the data!
        end %if the chan exists
        end %if data exists
        end %if dataloc.d is not empty
        
    end %end normdata collection loop
    
    %now calculate the method of normalization requested
    switch inp.normalizetype
        case 'znorm' %Z-Score normalize (subtract mean of all data, divide by std dev)
            MegaDataHolder = MegaDataHolder - min(MegaDataHolder,[],2); %subtract the min first
            Subtract = nanmean(MegaDataHolder, 'all'); %get the mean
            DivideBy = std(MegaDataHolder, 0, 'all', 'omitnan'); %get std dev
           
        case 'minmax' %Min/Max (percentile) normalize (subtract min, divide by diff btwn max and min
            MegaDataHolder = MegaDataHolder - min(MegaDataHolder,[],2); %subtract the min first
            ChanMax = prctile(prctile(MegaDataHolder, inp.normmax,2), inp.normmax); % get the max for normalizing data
            Subtract = prctile(prctile(MegaDataHolder,inp.normmin,2), inp.normmin); % get the min for normalizing data
            DivideBy = ChanMax - Subtract; % get the diff for normalizing data
       
        case 'minmaxself' %Min/Max (percentile) normalize (subtract min, divide by diff btwn max and min
        
        case 'minsub' %If you have a known control set of wells for Min/Max (percentile) normalize (subtract min, divide by diff btwn max and min
            Subtract = 0; % get the min for normalizing data
            DivideBy = 1; % divide by 1
            
        case 'robust' %Robust Scalar (scales to median and quantiles)
            Subtract = nanmedian(MegaDataHolder, 'all'); %get the median
            Quantilz = quantile(MegaDataHolder, (0:0.25:1),'all');
            Q25 = Quantilz(2); Q75 = Quantilz(4); %get the 25th and 75th quantiles
            DivideBy = Q75 - Q25; %get the interquantile range
            
        case 'controlrobust' %Robust Scalar (scales to median and quantiles) with known control wells
            MegaDataHolder = MegaDataHolder - min(MegaDataHolder,[],2); %subtract the min first
            Subtract = nanmedian(MegaDataHolder, 'all'); %get the median
            Quantilz = quantile(MegaDataHolder, (0:0.25:1),'all');
            Q25 = Quantilz(2); Q75 = Quantilz(4); %get the 25th and 75th quantiles
            DivideBy = Q75 - Q25; %get the interquantile range
            
        otherwise %use the default of zscore, but warn them about it
            fprintf('Normalization method %s given as input not recognized, using Z-Score normalization instead. \n', p.normalize)
            inp.normalize = 'znorm';
            Subtract = nanmean(MegaDataHolder, 'all'); %get the mean
            DivideBy = std(MegaDataHolder, 0, 'all', 'omitnan'); %get std dev
    end

    %loop over all the xys and again to normalize them
    for iXY = 1:NumXYs
        tXY = XYnums(iXY);
        if ~isempty(dataloc.d{tXY}) %check the XY isn't empty
        if isfield(dataloc.d{tXY}, 'data') %check the XY has data
        if isfield(dataloc.d{tXY}.data, CurChan) %check if that XY has the channel in it
            
            switch inp.normalizetype
                case 'minsub'
                    minSub = min(dataloc.d{tXY}.data.(CurChan),[],2); %get every cell's min
                    dataloc.d{tXY}.data.(NormChan) = dataloc.d{tXY}.data.(CurChan) - minSub; %subtract every cell's min from itself
                case 'minmaxself'
                    Subtract = min(dataloc.d{tXY}.data.(CurChan),[],2);
                    DivideBy = max(dataloc.d{tXY}.data.(CurChan),[],2) - Subtract;
                    dataloc.d{tXY}.data.(NormChan) = ((dataloc.d{tXY}.data.(CurChan)) - Subtract) ./ DivideBy; %make normalized data
                otherwise
                    dataloc.d{tXY}.data.(NormChan) = ((dataloc.d{tXY}.data.(CurChan)) - Subtract) ./ DivideBy; %make normalized data

            end %end switch
            
        end %if the chan exists
        end %if data exists
        end %if dataloc.d is not empty
    end %end second iXY loop
end %iChan loop

end %for normalize data function

%% Remake grans per cell Data Function
function dataloc = remakeGransPerCell_Data(dataloc, inp)
% dataloc = remakeGransPerCell_Data(dataloc, inp)
% Takes in dataloc structure and remakes the granspercell field for each xy
% by taking the total number of granules per tp (in that xy) and divides it by the total number of cells in that tp

%process either certain xys or all of them
if ~isempty(inp.certainxys)
    XYnums = inp.certainxys;
else
    %get the xy names from the data
    XYnums = extractAfter(dataloc.xys,'XY '); XYnums = str2double(XYnums)';
end
goodXYs = ~arrayfun(@(x)isempty(dataloc.d{x}),XYnums); % check there is data
XYnums = XYnums(goodXYs); % and only use those xys

% check the xys have NumGrans field
goodXYs = arrayfun(@(x)isfield(dataloc.d{x}.data,'NumGrans'),XYnums);
% only do xys with the right data
XYnums = XYnums(goodXYs);

for iXY = XYnums
    dataloc.d{iXY}.data.t_count_cells = sum(~isnan(dataloc.d{iXY}.data.NumGrans),1); % how many cells
    sumGrans = sum(dataloc.d{iXY}.data.NumGrans,1,'omitnan'); % how many grans
    dataloc.d{iXY}.data.t_granspercell = sumGrans./dataloc.d{iXY}.data.t_count_cells; % how many granulespercells
end %iXY loop

end

%% Folder Finder/Maker
function [dataloc, inp]= FolderFinderMaker(dataloc, inp)
saveFlag = false; loadFlag = false;

%add these paths 
addpath('\\albecklab.mcb.ucdavis.edu\Data\imageData\SG_4B','\\albecklab.mcb.ucdavis.edu\data\Code\Image Analysis','\\albecklab.mcb.ucdavis.edu\data\Code\Cell Trace','\\albecklab.mcb.ucdavis.edu\data\Code\Nick')
if isempty(inp.dataloc)

%check for different input folder, else figure it out
if ~isempty(inp.difffolder); icalledyou = inp.difffolder; 
    if icalledyou(end) == '\'; icalledyou = icalledyou(1:end-1); end  %make sure you didn't add a \ at the end./
else; whocalledme = dbstack('-completenames'); whocalledme = whocalledme(end).file; [icalledyou, ~, ~] = fileparts(whocalledme);
end 

%% Get the Base file name
[~,b] = system('net use');
tScan = textscan(b,'%s');
tScan = tScan{1};
hasLetters = cellfun(@(x)regexpi(x,'\S:'),tScan,'UniformOutput',false);
nothingHere = cellfun(@isempty,hasLetters);
hasLetters(nothingHere)=deal({0});
hasLetters = logical(cell2mat(hasLetters));
theseLetters = tScan(hasLetters);
theseLetters(:,2) = tScan(find(hasLetters)+1);
itHasThis = cell2mat(cellfun(@(x)contains(icalledyou,x),theseLetters(:,1),'UniformOutput',false));
if any(itHasThis)
    replaceThisLetter = find(itHasThis);
    icalledyou = strrep(icalledyou,theseLetters(replaceThisLetter,1),theseLetters{replaceThisLetter,2});
end
if iscell(icalledyou); icalledyou = icalledyou{1}; end
filebase = extractAfter(icalledyou,find(icalledyou == '\', 1,'last'));
if isempty(dataloc.file.base); dataloc.file.base = filebase; end

thisFolder = icalledyou;

%% Make sure we are not accidentally processing a random folder or the code folder
if ~contains(thisFolder,{'imageData','Processed Data'})
    error('DatalocHandler cannot load image data from %s.\nPlease check where you are calling the fullproc file from.\n', thisFolder)
end

%% Check for, collect, and/or make folder and file names 
% Get the Base file name
if isempty(dataloc.file.base); dataloc.file.base = extractAfter(thisFolder,find(thisFolder == '\', 1,'last')); end

% Get the data folder
if isempty(dataloc.fold.data); dataloc.fold.data = thisFolder; end

% Make or Get the processed data file
if isempty(dataloc.file.proc) && isempty(inp.savename)
    a = dir([dataloc.fold.data,'\',dataloc.file.base,'_Processed.mat']);
    if ~isempty(a); dataloc.file.proc = a.name; loadFlag = true;
    else; dataloc.file.proc = [dataloc.file.base,'_Processed.mat']; saveFlag = true;
    end 
    clear a;
elseif ~isempty(inp.savename)
    a = dir([dataloc.fold.data,'\',inp.savename,'_Processed.mat']);
    if ~isempty(a); dataloc.file.proc = a.name; loadFlag = true;
    else; dataloc.file.proc = [inp.savename,'_Processed.mat']; saveFlag = true; 
    end 
    clear a;
end

if loadFlag; a = load(fullfile(dataloc.fold.data, dataloc.file.proc),'dataloc'); a = a.dataloc;
    b = VersionCheck(dataloc, a, inp); if ~isempty(b); dataloc = b; end; clear a b %FIX ELSE BREAK
end

if isempty(dataloc.procdate); dataloc.procdate = '1-May-2000 01:01:01'; end

% Get the XY data folders from the processed data folder
if isempty(dataloc.fold.proc); dataloc.fold.proc = strrep(dataloc.fold.data,'imageData','Processed Data'); end
if ~isfolder(dataloc.fold.proc); mkdir(dataloc.fold.proc); end % if the folder doesn't exist already, make it
  
% look for XYs data in the Processed Data folder
a = dir(dataloc.fold.proc);
if ~isempty(a)
    a = a([a.isdir]); a = {a(3:end).name}'; 
    a = a(contains(a,'XY ')); a = string(a);
    if ~isempty(dataloc.xys); diffXYs = setdiff(a,dataloc.xys); 
        if ~isempty(diffXYs); inp.extractdata = true; saveFlag = true; end
    end
end
    
if inp.extractdata || isempty(dataloc.xys)
    if ~isempty(a); dataloc.xys = a; end
    clear a; saveFlag = true; 
end

% Check for IF data, if it exists append the path
if exist([dataloc.fold.proc,'\IF'],"dir") && exist([dataloc.fold.proc,'\IF\Cells.csv'],"file") %FIX
    dataloc.fold.if = [dataloc.fold.proc,'\IF'];
end

% Get / make the figure folder
if isempty(dataloc.fold.fig)
    dataloc.fold.fig = append(dataloc.fold.data,'\Figures\');
end

if ~isfolder(dataloc.fold.fig); mkdir(dataloc.fold.fig); end

%FIX: add check for path, if not add it, addpath(dataloc.fold.data); addpath(dataloc.fold.fig);

% Get the platemap info, update if needed
a = dir([dataloc.fold.data,'\','*','_platemap.xlsx']); 
if ~isempty(a); a = a(~contains({a.name},'~')); end

if isempty(dataloc.platemapd.date); dataloc.platemapd.date = '1-May-2000 01:01:01'; end % if nothing was found, give it a date
if ~isempty(dataloc.platemapd.date) && ~isdatetime(dataloc.platemapd.date); dataloc.platemapd.date = datetime(dataloc.platemapd.date); end
if dataloc.platemapd.date < datetime(a.date) || inp.updatemap; dataloc.platemapd.pmd = []; fprintf('Platemap old or is requested to update. \n'); end %replace the old platemap if requested or it has been updated 

if isempty(dataloc.file.platemap) %add the platemap's directory if it doesn't have one
    if isempty(a); fprintf('Error: Could not locate a platemap for this experiment.\n'); %FIX FOR BREAK
    else; dataloc.file.platemap = a.name; dataloc.platemapd.date = a.date; end 
end; clear a

% Map the platemap if it hasn't been done yet or has been updated
if ~isempty(dataloc.file.platemap) && isempty(dataloc.platemapd.pmd)
        a = fullfile(dataloc.fold.data, dataloc.file.platemap);
        [dataloc.platemapd.pmd, dataloc.platemapd.idx] = iman_readdatasheet(a); dataloc.platemapd.date = datetime; 
        fprintf('Updating Platemap. \n'); saveFlag = true; clear a
end

if saveFlag; save(fullfile(dataloc.fold.data, dataloc.file.proc), 'dataloc'); end
elseif ~isempty(inp.dataloc); dataloc = inp.dataloc; inp = rmfield(inp,'dataloc');
end % loaded data bypass

end %folderfindermaker function

%% extractdata
function [dataloc, inp] = ExtractData(dataloc, inp)

XYf = dataloc.xys; 
XYn = extractAfter(XYf, 'XY '); 
XYn = str2double(XYn)'; %get the xy folder and the xy number

for iXY = 1:numel(XYn)

    tXY = XYn(iXY); % xy number
    fXY = XYf(iXY); % xy file

    dataloc.d{tXY} = {}; %clear that xy
    if exist(fullfile(dataloc.fold.proc,fXY,'Cells.csv'),"file")
        cellData = readtable(fullfile(dataloc.fold.proc,fXY,'Cells.csv'),'VariableNamingRule','preserve','Delimiter',',');

        if size(cellData,1) < 1; continue; end % if there is no data skip this XY
        
        cFields = fieldnames(cellData); % get the fields
        cTrackField = cFields{startsWith(cFields,'TrackObjects_Label')}; % find the track objects field

        cellIdxs = unique(cellData.(cTrackField)); % get the cell idxs  (this indexes starting at 1)
        dataloc.d{tXY}.cellindex = cellIdxs; % assign them as cell idxs

        maxT = max(cellData.Metadata_SizeT); % how long is the data (this indexes starting at 1)
        maxCells = max(cellIdxs); % get the max number of tracked grans (this indexes starting at 1)
        nanBoxCells = nan(maxCells, maxT); % create a NaN box for granules

        % map the cells and time to a matrix
        theseTps = (cellData.Metadata_T)+1; % get the tps the values for column assignment (and add 1 to make tp 0 into 1, yay python vs matlab!)
        theseCells = cellData.(cTrackField); % get the tracked gran ids for row assignment (this indexes starting at 1)
        dataHome = sub2ind(size(nanBoxCells),theseCells,theseTps); % convert the row vectors into matrix locations
        goodData = ~isnan(dataHome); %find the nan data that won't be kept
        dataHome = dataHome(goodData); % get rid of the bad data

        % get only the fields we want
        isFieldies = ismember(cFields,inp.datafields);
        cFields = cFields(isFieldies); 

        for iData = 1:numel(cFields)
            tField = cFields{iData}; % field name
            tField = strrep(tField,'Intensity_','');
            tField = strrep(tField,'Mean_Grans_','Gran');
            tField = strrep(tField,'_GFP','');
            dataloc.d{tXY}.data.(tField) = nanBoxCells; % create a nan box to put the data into
            dataloc.d{tXY}.data.(tField)(dataHome) = cellData.(cFields{iData})(goodData); % pull that cells data ignoring bad data
        end
        
        %% change the x and y coord names 
        dataloc.d{tXY}.data.XCoord = dataloc.d{tXY}.data.Location_Center_X; dataloc.d{tXY}.data = rmfield(dataloc.d{tXY}.data,'Location_Center_X');
        dataloc.d{tXY}.data.YCoord = dataloc.d{tXY}.data.Location_Center_Y; dataloc.d{tXY}.data = rmfield(dataloc.d{tXY}.data,'Location_Center_Y');

        %% get the single cell grans per cell (if avaliable)
        if any(ismember(fieldnames(cellData),'Children_Grans_Count'))
            dataloc.d{tXY}.data.NumGrans = nanBoxCells; % make a nan box for the cells data
            dataloc.d{tXY}.data.NumGrans(dataHome) = cellData.Children_Grans_Count(goodData); % now assign the data to the d structure
        end
        % dataloc.d{tXY}.data.CytoIntegratedIntensity = cellData Intensity_IntegratedIntensity_GFP
    end % pulling the cell.csv data if it exists
    
    if exist(fullfile(dataloc.fold.proc,fXY,'Image.csv'),"file")
        imageData = readtable(fullfile(dataloc.fold.data,fXY,'Image.csv'),'VariableNamingRule','preserve','Delimiter',',');
        if any(ismember(fieldnames(imageData),'Count_Grans'))
            dataloc.d{tXY}.data.t_granspercell = (imageData.Count_Grans./imageData.Count_Cells)'; % grans per cell (mean)
        end
        if any(ismember(fieldnames(imageData),'Count_Grans'))
            dataloc.d{tXY}.data.t_count_cells = (imageData.Count_Cells)'; % grans per cell (mean)
        end
    end % image data file

end %for each XY
end % extractdata

%% Check for dataloc, its version, and allocate the info accordingly
function dout = VersionCheck(dout,din,inp)
    if isfield(din,'dataloc'); din = din.dataloc; end
    if isfield(din, 'version')
        switch din.version 
            case '1'
            otherwise %Start alllll over 
                dout = [];
        end
        tFields = fieldnames(dout.fold);
        for iFold = 1:numel(tFields)
            cFold = char(dout.fold.(tFields{iFold}));
            if ~isempty(cFold) && cFold(end) == '\'; dout.fold.(tFields{iFold}) = dout.fold.(tFields{iFold})(1:end-1); end  %make sure you didn't add a \ at the end./
        end
    end
end

%% Make IF Dataframe
% [datalocTable] = makeIFDF(dataloc, varargin)
%
% Optional Inputs:
% IFpath = required input if IF data is not in folder called "IF" in the
% experiment folder.
% 
% Copyright. Nicholaus DeCuzzi. 2023
function dataloc = makeIFDF(dataloc, p)
p.varNames = { 'exp',   'full',  'cell', 'treatment',  'xy',   'data', 'cellidx', 'txinfo'}; % variable names for the table
p.varTypes = {'string','string','string', 'string',    'string', 'cell',  'cell',  'cell'}; % variable types for the varNames in the table
p.expar = {'Cell', 'pTx', 'Tx1', 'Tx2','Tx3','Tx4'}; % expar for ct_maptx
p.reduce = false;
p.tpPerHr = 60/p.looptime; %set up how many tps per hour

%% Set up info for dataframe
% Set up the dataframe
mappedTxs = ct_maptx(dataloc.platemapd.idx,'expar',p.expar,'reduce',false); % map the platemap
txNames = fieldnames(mappedTxs);% get the unique datasets to either combine or not

emptyXYs = []; % set up empty XY list
xyList = 1:size(dataloc.d,2); % set up xyList
hasD = arrayfun(@(x) ~isempty(dataloc.d{x}), xyList); % get the xys for this dataset and check that the d is not empty
emptyXYs = [emptyXYs, xyList(~hasD)]; % get the xys that have no data
xyList = xyList(hasD); % get the non empty ds
hasCellIdx = arrayfun(@(x) ~isempty(dataloc.d{x}.cellindex), xyList); % get the xys that have cellindicies
emptyXYs = [emptyXYs, xyList(~hasCellIdx)]; % get the xys that have no cell data
xyList = xyList(hasCellIdx); % get the xys that have no data in them (perhaps over filtered data).

ifData = readtable(fullfile(dataloc.fold.if,'Cells.csv'),'VariableNamingRule','preserve','Delimiter',',');
nRow = size(ifData,1);
ifData.treatment = repmat("", nRow, 1); 
ifData.cell = repmat("", nRow, 1);
ifData.full = repmat("", nRow, 1); 


% loop through every uniquely mapped cell/condition combo (keeping them together in the table makes it easier to check if done properly)
for iMap = 1:size(txNames,1)

    txXys = intersect(mappedTxs.(txNames{iMap}).xy,xyList'); % get the xys for this mapped info and make sure they have data

    matchesXY = [];
    for iXY = txXys'
        matchesXY = [matchesXY, ifData.ImageNumber == iXY];
    end
    matchesXY = any(matchesXY,2);
    numData = sum(matchesXY);

    if numData > 0
        ifData.full(matchesXY) = repelem(string(strrep(txNames{iMap},'_',' ')),numData)'; % add the full sorted name
    
        holdThis = '';
        tTxList = mappedTxs.(txNames{iMap}).tx;
        notCells = ~contains({tTxList.dunit},'cells')'; % find the rows that dont have cell line info (FIX - for now?)
        if any(~notCells)
            ifData.cell(matchesXY) = repelem(strrep(string(mappedTxs.(txNames{iMap}).tx(~notCells).name),'_',' '),numData);
        else; ifData.cell(matchesXY) = repelem(" ",numData)'; % if no cell line fill it with nothing
        end
        
        tTxList = tTxList(notCells); % pull the treatment data
        for iTx = 1:size(tTxList,1)
            if tTxList(iTx).tunit == "h"; holdThis = append(holdThis, append(num2str(tTxList(iTx).dose), tTxList(iTx).dunit, ' ', tTxList(iTx).name,' at hour ', num2str(tTxList(iTx).time)),' and ');
            elseif tTxList(iTx).tunit == "tp" 
                    holdThis = append(holdThis, append(num2str(tTxList(iTx).dose), tTxList(iTx).dunit, ' ', tTxList(iTx).name,' and ')); % round to the nearest hour (.5 >= rounds up to 1, <0.5 rounds down to 0)
            end
        end
        holdThis = char(holdThis);
        holdThis = string(holdThis(1:end-5));
        holdThis = strrep(holdThis,'_',' ');
        ifData.treatment(matchesXY) = repelem(holdThis,numData)';
    end
end
ifData = movevars(ifData,{'full','cell','treatment'},'Before',1); % move the data information to the first columns of the dataframe
dataloc.ifd = ifData;

end %end of makeIFDF(dataloc, varargin) function