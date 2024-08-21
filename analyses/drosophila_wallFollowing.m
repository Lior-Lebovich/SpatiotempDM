%% QUANTIFYING GLOBAL LATERAL TENDENCIES WITH MAD SCORES: 
%
% Computes flies' Median Absolute Deviation (MAD) around the horizontal (x)
% center during upward and/or downward motion in bottom arm. Plots 
% MAD-based comparison of TPI(y) and MADs datasets comparison.
%
% INPUT FILES: 
%   output_drosophila_main_[DATASET_NAME(10:END)].mat
%   flies_TPIs.mat
%
% FUNCTION:
%   TPI.m for computing the TPI.
%
% OUTPUT FILES: 
%   output_drosophila_MADs.mat: MAD scores of individual decision-makers 
%       within the bottom arm (excluding cul-de-sac). 
%       MADs.[DATASET].xMadArmMid.(DIR)(k) is the MAD around
%       the horizontal middle of decision-maker k in a given DATASET for 
%       walking dir.=DIR.
%   output_drosophila_MADs_order.mat: ordered MAD scores.
%       madOrdStrct.[DATASET].xMadArmMid.(DIR)(k) is the order index.
%
% FIGURES:
%   Within dataset MAD-based TPI(y) (Avg. +- SEM) comparison:
%       Figs. 2E (S2A), S3E, S11C: 
%       flies_MAD[DATASET]_BasedTPI_[DATATYPE(10:end-2)]_upArmWalk.fig
%   Between dataset MAD comparison - fly data in short maze:
%       Fig. 4B: flies_MADcenter_allShortsCompShortBox_upArmWalk
%   Between dataset MAD comparison - ABM simulation data:
%       Fig. S13B: flies_MADcenter_allShortsCompModelBox_upArmWalk
%
% More information about MAD score and corresp. analyses in the sections 
% below.
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



%% Define datasets and arm directions:

dataTypes_stable = {'Turndata_ShortYy', 'Turndata_LongYy', ...
    'Turndata_ShortNorpAYy', 'Turndata_ShortDumbYy', ...
    'Turndata_ShortFoxPYy', 'Turndata_ShortNompCYy', ...
    'Turndata_ShortABM_NoWF', ...
    'Turndata_ShortABM_WF_00100', 'Turndata_ShortABM_WF_00200', ...
    'Turndata_ShortABM_WF_00500', 'Turndata_ShortABM_WF_00700', ...
    'Turndata_ShortABM_WF_00900', ...
    'Turndata_ShortABM_WallFollowing', 'Turndata_ShortABM_brown5'};
dataCols_stable = {'k', 'b', 'r', 'm', 'g', '#77AC30', ... % fly data
    "#A2142F", ... % no WF 
    "#D95319", "#EDB120", "#77AC30", "#4DBEEE", "#0072BD", ... % WF parm.
    "#7E2F8E", "#7E2F8E" }; % WF, brownian

dataTypes = dataTypes_stable; 
dataCols = dataCols_stable;

% Define upper edge of cul-de-sac. Same for all short fly mazes
for d = 1:length(dataTypes)
    dataType = dataTypes{d};
    d_stable = d;
    if strcmp(dataType,'Turndata_LongYy')
        culDeSacEdge.(dataType) = 0.68/3.36;
    elseif strcmp(dataType,'Turndata_ShortYy') || ...
            strcmp(dataType,'Turndata_ShortNorpAYy') || ...
            strcmp(dataType,'Turndata_ShortDumbYy') || ...
            strcmp(dataType,'Turndata_ShortFoxPYy') || ...
            strcmp(dataType,'Turndata_ShortNompCYy') 
        culDeSacEdge.(dataType) = 0.69/2.03;
    elseif strcmp(dataType,'Turndata_ShortABM_NoWF') || ...
            strcmp(dataType,'Turndata_ShortABM_WallFollowing') || ...
            strcmp(dataType,'Turndata_ShortABM_WF_00100') || ...
            strcmp(dataType,'Turndata_ShortABM_WF_00200') || ...
            strcmp(dataType,'Turndata_ShortABM_WF_00500') || ...
            strcmp(dataType,'Turndata_ShortABM_WF_00700') || ...
            strcmp(dataType,'Turndata_ShortABM_WF_00900') || ...
            strcmp(dataType,'Turndata_ShortABM_brown5')
        culDeSacEdge.(dataType) = 7.5/30;
    end
end

clearvars -except dataTypes culDeSacEdge dataTypes_stable dataCols_stable;
close all; clc

% Define arm walking directions:
armDirNames = {'down', 'up', 'both'}; 
armDirLocs = [-2,-1; 1, 2; -2, 2]; 

% MAD types (can possibly replace 'MAD around center' with 'MAD around 
% median'):
namesMAD = {'MAD around center'}; 



%% Load data, compute MAD(x) for each fly when walking up/down/both the arm:

% MAD(x) is computed around the center [to compute MAD around the median 
    % change above: namesMAD = {'MAD around center'} to: 
    % {'MAD around median'}].
% MAD(x) is computed for trajectories within the bottom arm - between the 
    % cul-de-sac and the intersection - separately for 3 walking 
    % directions: (1) 'up': when the fly walks up the arm (toward 
    % intersection), 'down': when the fly walks down the arm (toward 
    % cul-de-sac) or 'both': both up and down the arm.
% flyTrlsXDir - all of the fly's x-locations between the cul-de-sac and the 
    % intersection, for the relevant walking dir..
% flyTrlsXDir_Nrmlzd - noramlized flyTrlsXDir (between 0 and 1).
% xMadArmMid - mean absolude deivation from center (1/2), based on 
    % flyTrlsXDir_Nrmlzd.
% xMadArmMed - mean absolude deivation around the median, based on 
    % flyTrlsXDir_Nrmlzd.

% subsampling frames (1= no subsamling):
kSubSamp = 1; 


% Compute MADs for each dataset, fly and arm walking dir.:

for d = 1:length(dataTypes)

    d_stable = d;
    dataType = dataTypes{d};

    % Load dataset:
    matFileName = ['output_drosophila_main_' dataType(10:end) '.mat'];
    load(matFileName);
    
    % Keep dataType and Cols despite val's in matFileName.m:
    d = d_stable;
    dataTypes = dataTypes_stable; 
    dataCols = dataCols_stable;
    dataType = dataTypes{d};
   
    % Load arm-edge(y) def. (MAD's will be based on coordinates within 
    % armEdges):
    armTopEnd = 1; 
    armEdges = [culDeSacEdge.(dataType), armTopEnd]; 
    
    % Create structures for storing flys' MADs in each walk dir.:
    for dirr = 1:numel(armDirNames)
        armDirName = armDirNames{dirr};
        MADs.([dataType(10:end-2)]).xMadArmMed.(armDirName) = ...
            nan( numel(selectedFlies), 1 );
        MADs.([dataType(10:end-2)]).xMadArmMid.(armDirName) = ...
            nan( numel(selectedFlies), 1 );
    end
    fliesPTrialsConsidered = nan( numel(selectedFlies), 1 ); 
    

    % Run over flies, compute MADs in each dir and plot MAD-based
    % separation of TPI(y)'s:

    for nF = 1:length(selectedFlies) 
    
        % Read fly's data:
        f = selectedFlies(nF);
        fly_xytyd = meanX_vects_ALL.(['f' num2str(f)]).xytyd;
        fly_locsBotTrialT = meanX_vects_ALL.(['f' num2str(f)]).locsBotTrialT;
        fly_Dec = meanX_vects_ALL.(['f' num2str(f)]).Dec;
        
        % Read fly's botttom trials:
        fly_xytyd_good = fly_xytyd( fly_locsBotTrialT );

        % Create structure for storing x coordinates for each walk dir.:
        for dir = 1:numel(armDirNames)
            armDirName = armDirNames{dir};
            flyTrlsXDir.(armDirName) = [];
        end
        
        flyTrialsConsidered = zeros( numel(fly_xytyd_good), 1); %%%

        for t = 1:numel(fly_xytyd_good)
    
            % read fly's x,y in trial (downsampling):
            fly_xytyd_good_trial = fly_xytyd_good{t};
            fly_xGood = fly_xytyd_good_trial(1:kSubSamp:end,1);
            fly_yGood = fly_xytyd_good_trial(1:kSubSamp:end,2);
    
            % Only use trials with full data:
            if min(fly_yGood)<=-armTopEnd && max(fly_yGood)>=armTopEnd
                
                flyTrialsConsidered(t) = 1; %%%

                % Run over arm walking dir's (up, down, both) and store 
                % corresp. x traj. data (constraints over y):
                for dir = 1:numel(armDirNames)
                    armDirName = armDirNames{dir};
                    armDirLoc = abs( armDirLocs(dir,:) );
                    armDirPol = sign( armDirLocs(dir,:) );

                    % Bottom arm x coordinates within y arm dir. def. 
                    % (3rd constraint relevant for "both" cond.):
                    xInArmYDir = fly_xGood( ...
                        (fly_yGood >= (armEdges(armDirLoc(1))*armDirPol(1))) & ...
                        (fly_yGood <= (armEdges(armDirLoc(2))*armDirPol(2))) & ...
                        (abs(fly_yGood) >= culDeSacEdge.(dataType)) );

                    % Store with data from additional trials:
                    flyTrlsXDir.(armDirName) = [flyTrlsXDir.(armDirName); ...
                        xInArmYDir];
                end
    
            end
    
        end

        fliesPTrialsConsidered(nF) = mean( flyTrialsConsidered );

        % Run over walking dir's and compute the fly's MAD scores:
        for dir = 1:numel(armDirNames)
            armDirName = armDirNames{dir};

            % Center x such that 0<=x''<=1:
            flyTrlsXDir_0correct = flyTrlsXDir.(armDirName) - ...
                min(flyTrlsXDir.(armDirName)); % --> min(x') = 0
            flyTrlsXDir_Nrmlzd = flyTrlsXDir_0correct / ...
                max(flyTrlsXDir_0correct); % --> x'' = [0,1]

            % Save Median Absolute Deviation around the MEDIAN: 
            MADs.([dataType(10:end-2)]).xMadArmMed.(armDirName)(nF) = ...
                mad(flyTrlsXDir_Nrmlzd,1); % = MEDIAN( ABS(X-MEDIAN(X)) )
            
            % Save Median Absolute Deviation around the MIDDLE: 
            MADs.([dataType(10:end-2)]).xMadArmMid.(armDirName)(nF) = ...
                median(abs(flyTrlsXDir_Nrmlzd-.5)); % = MEDIAN( ABS(X-.5) )
        end

        clear flyTrlsXDir;

    end

    clear meanX_vects_ALL;
    
end


% Save MADs:
save( 'datafiles/output_drosophila_MADs.mat', 'MADs', 'namesMAD' );



%% Compute short WT corr. for MAD(x|UP the arm) vs MAD(x|DOWN the arm):

[RHO,PVAL] = corr( MADs.Short.xMadArmMid.up, MADs.Short.xMadArmMid.down );  



%% Plot MAD(x)-based comparison of TPI(y) - only fly datasets:

% Grouping flies based on their MAD(x|UP the arm) and comparing the
% avg.+-SEM TPI's of the resulting groups.
% To plot MAD(x|DOWN the arm) instead change the value of armDirName below
% to 'down'.

armDirName = 'up';

dataTypes = {'Turndata_ShortYy', 'Turndata_LongYy', ...
    'Turndata_ShortNorpAYy', 'Turndata_ShortDumbYy', ...
    'Turndata_ShortFoxPYy', 'Turndata_ShortNompCYy'}; 

% Determine #flies in group:
nInGroupS = [11, 10, 50, 18, 49, 11]; % for MAD-based TPI comparison. 
                                  % For equal-size groups.

% Groups for which avg. TPI(y) will be plotted (use 'all' for all groups, 
% use 'extreme' for the 2 extreme groups):
groupTpiToPlot = 'extreme'; 

% Group colors:
colsGroups.Short = {'r', 'm', 'c', 'b', 'k' };
colsGroups.Long = {'r', 'm', 'c', 'b', 'k' };
colsGroups.ShortNorpA = {'r', 'm', 'c', 'b', 'k' };
colsGroups.ShortDumb = {'r', 'm', 'c', 'b', 'k' };
colsGroups.ShortFoxP = {'r', 'm', 'c', 'b', 'k' };
colsGroups.ShortNompC = {'r', 'm', 'c', 'b', 'k' };

% Load TPI(y)'s and corresp. bin-edges for flies in both datasets:
load('flies_TPIs.mat');

% Run over datasets and MAD types, MAD-base group flies, plot group TPI:

for d = 1:length(dataTypes)
    dataType = dataTypes{d};

    % Load TPI(y)'s and corresp. bin-edges for flies in dataset:
    TIP_y_rFly_cBin = savedTPI.([dataType(10:end-2)]).y;
    yEdges = binEdgesTPI.([dataType(10:end-2)]).y;
    % Compute bin-centers:
    yCenters = .5*diff(yEdges(1:2)) + yEdges(1:end-1);

    % Load groups' color for current dataset:
    colGroups = colsGroups.(dataType(10:end-2));

    % Load group size for current dataset and derive #groups:
    nInGroup = ceil( size(TIP_y_rFly_cBin,1) / 5 );
    nGroups = ceil( size(TIP_y_rFly_cBin,1)/nInGroup );
    
    % Determine the groups for which TPI's will be plotted:
    if strcmp(groupTpiToPlot,'extreme')
        numGrTpiToPlot = [1, nGroups];
    elseif strcmp(groupTpiToPlot,'all')
        numGrTpiToPlot = 1:nGroups;
    end


    % Load MAD values, order flies, group flies, plot groups' TPI(y):

    for m = 1:numel(namesMAD)
        nameMAD = namesMAD{m};
        clear plt;

        % Load MAD variable name (MAD type):
        if strcmp( nameMAD, 'MAD around median' )
            MADVarName = 'xMadArmMed';
        elseif strcmp( nameMAD, 'MAD around center' )
            MADVarName = 'xMadArmMid';
        end

        % Load MAD values for a given dataset and MAD type:
        thisMAD = MADs.([dataType(10:end-2)]).(MADVarName).(armDirName);
        fliesCorrespGroups = nan( size(thisMAD) );

        % Order flies in dataset based on their MAD scrores:
        [~,madOrd] = sort( thisMAD );

        % For storing within-group extreme MAD values (group edges in inset):
        minMaxMadGroup = nan(nGroups,2);

        % Create figure:
        figMADBasedTpi.([dataType(10:end-2)]).(MADVarName) = figure;


        % Plot (main): MAD(x)-based fly-groups' TPI(y)'s (avg.+-SEM):

        for gr = 1:nGroups

            % Load flies' locations in current group:
            if gr == nGroups
                locFliesInGroup = madOrd( (1+(gr-1)*nInGroup):end );
            else
                locFliesInGroup = madOrd( (1+(gr-1)*nInGroup):(gr*nInGroup) );
            end
            % Load flies' MADs in current group:
            xMadInGroup = thisMAD(locFliesInGroup);
            % Store current group's extreme MADs (for inset):
            minMaxMadGroup(gr,:) = [min(xMadInGroup), max(xMadInGroup)];
            % Load TPI(y) of current group:
            TPI_group = TIP_y_rFly_cBin( locFliesInGroup, : );
            
            % Plot TPI(y) avg. +- SEM of flies in group:
            if sum( gr == numGrTpiToPlot ) == 1
                figure( figMADBasedTpi.([dataType(10:end-2)]).(MADVarName) );
                plt(gr) = errorbar( yCenters, mean( TPI_group, 'omitnan' ), ...
                    std( TPI_group, 'omitnan' ) ./ sqrt(sum(~isnan(TPI_group))), ...
                    colGroups{gr}, 'lineWidth', 1 ); 
                hold on;
            end

            fliesCorrespGroups(locFliesInGroup) = gr;
                
        end
        
        % axes (main):
        xlabel('y'); 
        xlim( yEdges([1,end]) );
        ylim([-1,1]);
        ylabel('group average TPI(y)');
        grid on;
        title([dataType(10:end-2) ' - TPI(y) by ' nameMAD]);


        % Plot (inset): histogram of MAD values of each group:

        figure( figMADBasedTpi.([dataType(10:end-2)]).(MADVarName) );
        axes('Position',[.16 .8 .4 .1]);
        for gr = 1:nGroups
            if gr == nGroups
                locFliesInGroup = madOrd( (1+(gr-1)*nInGroup):end );
            else
                locFliesInGroup = madOrd( (1+(gr-1)*nInGroup):(gr*nInGroup) );
            end
            xMadInGroup = thisMAD(locFliesInGroup);
            histogram( xMadInGroup, 'BinEdges', minMaxMadGroup(gr,:), ...
                'FaceColor', colGroups{gr} ); hold on;
        end

        % Derive xticks (for inset):
        histEdges = [minMaxMadGroup(1,1), ...
            .5 * (minMaxMadGroup(1:end-1,2) + minMaxMadGroup(2:end,1) )', ...
            minMaxMadGroup(end,2)];
        theXticks = [floor(histEdges(1)*100)/100, ...
            round(histEdges(2:end-1)*100)/100, ...
            ceil(histEdges(end)*100)/100];

        % axes (inset):
        xlabel('MAD(x)');
        xticks( theXticks ); 
        xtickangle(45); 
        xlim( [theXticks(1), theXticks(end)] );
        ylabel('#');
        yticks([0,nInGroup]);
        ylim([0,nInGroup]);
        ggg = gca; 
        ggg.YAxisLocation = 'right';
        ggg.FontSize = 11;

        
        % Save figure:
        figStartName = ['analyses/figures/flies_MAD' nameMAD(12:end) ...
            '_BasedTPI_' dataType(10:end-2) '_' armDirName 'ArmWalk'];
        saveas( figMADBasedTpi.([dataType(10:end-2)]).(MADVarName), ...
            [figStartName '.fig']);

        % Store ordered location of MADs (for use in external script):
        
        % madOrdStrct.(...) is the ordered selectedFlies (e.g.,
            % thisMAD(madOrd(1)) and thisMAD(madOrd(end)) are min. and max. 
            % MAD, resp..
        madOrdStrct.([dataType(10:end-2)]).(MADVarName).(armDirName) = ...
            madOrd;
        % fliesMadGroups.(...) is the group associated with each fly, s.t.,
            % fliesMadGroups(i) is the group assoc. with fly i in 
            % thisMAD(i) and in selectedFly(i).
        fliesMADGroupStrct.([dataType(10:end-2)]).(MADVarName).(armDirName) ...
            = fliesCorrespGroups;

    end
end

% Save MADs orders:
save( 'datafiles/output_drosophila_MADs_order.mat', 'madOrdStrct', ...
    'fliesMADGroupStrct' );



%% Compare: MAD(x) of short-maze fly datasets:

load( 'output_drosophila_MADs.mat' );

armDirName = 'up';

dataTypesCols = {'k', 'g', 'm', 'r', '#77AC30'};
dataTypes = {'Turndata_ShortYy', 'Turndata_ShortFoxPYy', ...
    'Turndata_ShortDumbYy', 'Turndata_ShortNorpAYy', ...
    'Turndata_ShortNompCYy' };

clear forLeg dataTypePrintNames;

for m = 1:numel(namesMAD)
    nameMAD = namesMAD{m};
    dataTypePrintNames = cell(size(dataTypes));
    if strcmp( nameMAD, 'MAD around median' )
        MADVarName = 'xMadArmMed';
    elseif strcmp( nameMAD, 'MAD around center' )
        MADVarName = 'xMadArmMid';
    end
    maxForSig = -Inf;
    
    % Create figure:
    figMADCompBox.(MADVarName) = figure;

    for d = 1:length(dataTypes)
        dataType = dataTypes{d};
        MAD_up = MADs.([dataType(10:end-2)]).(MADVarName).(armDirName);
        maxForSig = max( [maxForSig; MAD_up] );
        dataTypePrintName = dataType(10:end-2);
        if strcmp( dataTypePrintName, 'Short' )
            dataTypePrintNames{d} = 'WT';
        else
            dataTypePrintNames{d} = dataTypePrintName(6:end);
        end
    
        % Plot boxchart of MAD for each dataset:
        figure( figMADCompBox.(MADVarName) );
        boxchart( d*ones(size(MAD_up)), MAD_up, 'BoxFaceColor', ...
            dataTypesCols{d}, 'MarkerColor', dataTypesCols{d} );
        hold on;

    end

    % Axes - boxchart:
    figure( figMADCompBox.(MADVarName) );
    xticks(1:length(dataTypes));
    xticklabels( dataTypePrintNames );
    xtickangle(-45);
    ylabel([namesMAD{m} ' - walking up the arm']);

    % Add significance between WT-mutant comparison:
    MAD_up_WT = MADs.Short.(MADVarName).(armDirName);
    maxForSig = 1.01 * maxForSig;
    for d = 2:length(dataTypes)
        dataType = dataTypes{d};
        MAD_up_mutant = MADs.([dataType(10:end-2)]).(MADVarName).(...
            armDirName);
        p = ranksum( MAD_up_WT, MAD_up_mutant);
        % Plot if sign. (after mult. comp. correction):
        if p<=0.05/(length(dataTypes)-1)
            maxForSig = maxForSig + .02;
            plot( [1, d], maxForSig * [1,1], 'k-', ...
                (1+d)/2, .005 + maxForSig, '*k' );
            hold on;
        end
    end
    ylim([0,inf]);

    % Save figure - boxchart:
    xtickangle(0);
    figStartName = ['analyses/figures/flies_MAD' nameMAD(12:end) ...
        '_allShortsCompShortBox_' armDirName 'ArmWalk'];
    saveas( figMADCompBox.(MADVarName), [figStartName '.fig']);

end



%% Compare: MAD(x) of ABM versions:

load( 'output_drosophila_MADs.mat' );

armDirName = 'up';

dataTypes = {'Turndata_ShortABM_brown5', 'Turndata_ShortABM_NoWF', ...
    'Turndata_ShortABM_WF_00100', 'Turndata_ShortABM_WF_00200', ...
    'Turndata_ShortABM_WF_00500', 'Turndata_ShortABM_WF_00700', ...
    'Turndata_ShortABM_WF_00900', 'Turndata_ShortABM_WallFollowing'};    
modelPrintNames = {'Brownian', 'No WF', 'WF 001', 'WF 002', 'WF 005', ...
    'WF 007', 'WF 009', 'WF'};
dataTypesCols = {"#7E2F8E", "#A2142F", "#D95319", "#EDB120", "#77AC30", ...
    "#4DBEEE", "#0072BD", "#7E2F8E"};

clear forLeg dataTypePrintNames;

for m = 1:numel(namesMAD)
    nameMAD = namesMAD{m};
    dataTypePrintNames = modelPrintNames;
    if strcmp( nameMAD, 'MAD around median' )
        MADVarName = 'xMadArmMed';
    elseif strcmp( nameMAD, 'MAD around center' )
        MADVarName = 'xMadArmMid';
    end
    maxForSig = -Inf;
    
    % Create figure:
    figMADCompBox.(MADVarName) = figure;

    for d = 1:length(dataTypes)
        dataType = dataTypes{d};
        MAD_up = MADs.([dataType(10:end-2)]).(MADVarName).(armDirName);
        maxForSig = max( [maxForSig; MAD_up] );
        dataTypePrintName = dataType(10:end-2);
        if strcmp( dataType, 'Turndata_ShortABM_brown5' )
            dataCol = [.5, .3, .1];
        else
            dataCol = dataTypesCols{d};
        end
    
        % Plot boxchart of MAD for each dataset:
        figure( figMADCompBox.(MADVarName) );
        boxchart( d*ones(size(MAD_up)), MAD_up, 'BoxFaceColor', ...
            dataCol, 'MarkerColor', dataCol );
        hold on;

    end

    % Axes - boxchart:
    figure( figMADCompBox.(MADVarName) );
    xticks(1:length(dataTypes));
    xticklabels( dataTypePrintNames );
    xtickangle(0);
    ylabel([namesMAD{m} ' - walking up the arm']);
    ylim([0,inf]);

    % Save figure - boxchart:
    xtickangle(0);
    figStartName = ['analyses/figures/flies_MAD' nameMAD(12:end) ...
        '_allShortsCompModelBox_' armDirName 'ArmWalk'];
    saveas( figMADCompBox.(MADVarName), [figStartName '.fig']);

end


