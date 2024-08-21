%% T-DISTRIBUTED STOCHASTIC NEIGHBOR EMBEDDING OF BEHAVIORAL DATA TPI:
%
% Dimensionality reduction (2D) of the individual TPI curves (within the 
% bottom arm) of humans and all fly lines via t-distributed Stochastic 
% Neighbor Embedding (tSNE). With cosine distance metric. 
%
% INPUT FILE: 
%   flies_TPIs.mat
%
% FIGURES:
%   T-SNE embedding of behavioral TPI's:
%       Fig. 4G: flies_tSNE_TPI[DOMAIN].fig.
%
% Copyright (c) Lior Lebovich, 2024
% lebovich.lior@gmail.com


    
%% Use relevant directory and add path to external fucntions and datafiles:

cdName = 'SpatiotempDM';
currentFolder = cd;
if ~endsWith(currentFolder,cdName)
    cd(cdName);
end
addpath( 'analyses', 'datafiles', 'analyses/customfuns', '-frozen');



%% Load indiv TPIs, define datatypes, col's and spatial/temporal bins:

% Load TPI curves of individual decision-makers (behavioral data):
load( 'flies_TPIs.mat', 'savedTPI', 'binEdgesTPI', 'savedTPIinf' );

% Define datasets and corresp. colors:
datatypes = {'Turndata_ShortYy', 'Turndata_LongYy', ...
    'Turndata_ShortNorpAYy', 'Turndata_ShortDumbYy', ...
    'Turndata_ShortFoxPYy', 'Turndata_ShortNompCYy', ...
    'Turndata_ShortABM_NoWF', ...
    'Turndata_ShortABM_WF_00100', 'Turndata_ShortABM_WF_00200', ...
    'Turndata_ShortABM_WF_00500', 'Turndata_ShortABM_WF_00700', ...
    'Turndata_ShortABM_WF_00900', ...
    'Turndata_ShortABM_WallFollowing', 'Turndata_ShortABM_brown5', ...
    'Turndata_HumanYy'};
datatypes_TSNE = datatypes([1,2,5,4,3,6,end]);
dataCols_TSNE = [0,0,0; 0,0,1; 0,1,0; 1,0,1; 1,0,0; .47,.67,.19; .9,.7,.1]; 
dataNames = {'WT','Long','FoxP','Dumb','NorpA','NompC','Human'};

% Define spatial bins of intereset (within bottom arm):
yEdges = -1.3:.1:1.3;
yCenters = yEdges(1:end-1) + .5*diff(yEdges(1:2));
relLocsTpiY = find(yCenters>=-1 & yCenters<=1);

% Define temporal bins of intereset (within bottom arm):
tEdges = -1.275:.075:1.275;
tCenters = tEdges(1:end-1) + .5*diff(tEdges(1:2));
relLocsTpiT = find(tCenters>=-1 & tCenters<=1);



%% TSNE:

% Define distance metric:
distMetric = 'cosine';

% Read datasets' TPIs and names for all decision-makers (without missing 
% values):
tpiY_TSNE = [];
tpiY_TSNE_names = cell(1,0);
tpiT_TSNE = [];
tpiT_TSNE_names = cell(1,0);
for d = 1:numel(datatypes_TSNE)
    datatype = datatypes_TSNE{d};
    % Spatial domain:
    tpiY_dataset = savedTPI.(datatype(10:end-2)).y(:,relLocsTpiY);
    notNansSubsY = ( sum( isnan( tpiY_dataset ), 2) == 0 );
    tpiY_TSNE = [tpiY_TSNE; tpiY_dataset(notNansSubsY,:)];
    tpiY_TSNE_names = [tpiY_TSNE_names; ...
        repmat( dataNames(d), [sum(notNansSubsY),1] )];
    % Temporal domain:
    tpiT_dataset = savedTPI.(datatype(10:end-2)).time(:,relLocsTpiT);
    notNansSubsT = ( sum( isnan( tpiT_dataset ), 2) == 0 );
    tpiT_TSNE = [tpiT_TSNE; tpiT_dataset(notNansSubsT,:)];
    tpiT_TSNE_names = [tpiT_TSNE_names; ...
        repmat( dataNames(d), [sum(notNansSubsT),1] )];
end


% Plot TPI (within bottom arm) T-SNE in the spatial domain:
currFig = figure;
rng default % for reproducibility/ fair comparison
[Y, loss] = tsne( tpiY_TSNE, 'Algorithm' ,'exact' ,'Distance', ...
    distMetric, 'Standardize', true );
clear ppp;
for dat = 1:numel(dataNames)
    dataName = dataNames{dat};
    dataCol = dataCols_TSNE(dat,:);
    dataLocs = strcmp(tpiY_TSNE_names,dataName);
    ppp(dat) = scatter( Y(dataLocs,1), Y(dataLocs,2), 20, dataCol, ...
        'filled', 'MarkerFaceAlpha', .5 );
    hold on;
end
title([distMetric ', ' num2str(loss)]);
legend(ppp,dataNames)
axis square;
figStartName = 'analyses/figures/flies_tSNE_TPIy';
saveas( currFig, [figStartName '.fig']);


% Plot TPI (within bottom arm) T-SNE in the temporal domain:
currFig = figure;
rng default % for reproducibility/ fair comparison
[Y, loss] = tsne( tpiT_TSNE, 'Algorithm' ,'exact' ,'Distance', ...
    distMetric , 'Standardize', true);
clear ppp;
for dat = 1:numel(dataNames)
    dataName = dataNames{dat};
    dataCol = dataCols_TSNE(dat,:);
    dataLocs = strcmp(tpiT_TSNE_names,dataName);
    ppp(dat) = scatter( Y(dataLocs,1), Y(dataLocs,2), 20, dataCol, ...
        'filled', 'MarkerFaceAlpha', .5 );
    hold on;
end
title([distMetric ', ' num2str(loss)]);
legend(ppp,dataNames)
axis square;
figStartName = 'analyses/figures/flies_tSNE_TPIt';
saveas( currFig, [figStartName '.fig']);


