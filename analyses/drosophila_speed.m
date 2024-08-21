%% KINEMATICS BY SPATIAL AND TEMPORAL BINS IN THE SHORT WT DATASET:
%
% Plots avg. speed, distance and rel. time by y-bins and relative time-bins 
% for all flies in the short WT dataset.
%
% INPUT: 
%   output_drosophila_main_[DATASET_NAME(10:END)].mat
%
% FUNCTION:
%   ySpeed.m for computing kinematics over bins.
%
% OUTPUT: 
%   flies_Speed.mat: individual kinematics parameters (savedSpeed struct.) 
%       per bin edges (binEdgesSpeed). 
%       savedSpeed.DATATYPE(10:end-1).PARAMNAME.binType is 
%       numAnimals X numBinsBinType matrix.
%
% Figures:
%   Avg. +-% SEM kinematics param by spatial and temporal bins, overlaid 
%       on individual curves:
%       Fig. S1D: fliesShort_avg[PARAM_NAME]Y.fig
%       Fig. S1E: fliesShort_avg[PARAM_NAME]T.fig
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



%% Define relevant datasets and parameters of interest and load WT short data:

paramNames = {'speed', 'distance', 'relTime'};
paramTitles = {'speed', 'distance', 'relative time'};

datatypes = {'Turndata_ShortYy', 'Turndata_LongYy', ...
    'Turndata_ShortNorpAYy', 'Turndata_ShortDumbYy', ...
    'Turndata_ShortFoxPYy', 'Turndata_ShortNompCYy', ...
    'Turndata_humanMediumYy', 'Turndata_humanShortYy', ...
    'Turndata_ShortABM_NoWF', ...
    'Turndata_ShortABM_WF_00100', 'Turndata_ShortABM_WF_00200', ...
    'Turndata_ShortABM_WF_00500', 'Turndata_ShortABM_WF_00700', ...
    'Turndata_ShortABM_WF_00900', ...
    'Turndata_ShortABM_WallFollowing', 'Turndata_ShortABM_brown5'};
dataCols = {'k', 'b', 'r', 'm', 'g', '#77AC30', ... % fly data
    'c', 'y', ... % human data
    "#A2142F", ... % no WF 
    "#D95319", "#EDB120", "#77AC30", "#4DBEEE", "#0072BD", ... % WF parm.
    "#7E2F8E", "#7E2F8E" }; % WF, brownian

% Load processed fly short WT data:
matFileName = ['output_drosophila_main_' datatypes{1}(10:end) '.mat'];
load(matFileName);



%% Plot Speed / distance / relativeTime(bins) for the Short WT dataset:

datatype = datatypes{1};
dataCol = dataCols{1};

binTypes = {'y', 'time'};
binPrintNames = {'y', 'T'};

fps = 30;
s = RandStream('mlfg6331_64'); 

for p = 1:numel(paramNames)
    paramName = paramNames{p};
    paramTitle = paramTitles{p};

    if strcmp(datatype,'Turndata_ShortABM_brown5')
        dataCol = [.5, .3, .1];
    end
    
    % Define bottom arm cul-de-sac and maze ROIs:
    if strcmp(datatype,'Turndata_ShortYy') || ...
            strcmp(datatype,'Turndata_ShortNorpAYy') || ...
            strcmp(datatype,'Turndata_ShortDumbYy') || ...
            strcmp(datatype,'Turndata_ShortFoxPYy') || ...
            strcmp(datatype,'Turndata_ShortNompCYy')
        culDeSacY = 0.69/2.03;
        % Define maze ROIs:
        regionsBreak_Strct.y = [-1.3, -1, -culDeSacY, culDeSacY, 1, 1.3];
        regionPrintNames_Strct.y = {'Pre-arm', 'Down-the-arm', ...
            'Cul-de-sac', 'Up-the-arm', 'Post-arm'};
        regionsBreak_Strct.time = [-1.3, -1, 0, 1, 1.3];
        regionPrintNames_Strct.time = {'Pre-arm', ...
            'Down-the-arm & c.d.s<0', 'C.d.s>0 & up-the-arm', 'Post-arm' };
        regionBreakTicks.y = -1.:.5:1;
        regionBreakTicks.time = -1.:.5:1;
        infVal = 1.5;
        % Define (new) y-bin and time-bin edges:
        yEdges = -1.3:.1:1.3;
        tEdges = -1.275:.075:1.275;
    elseif strcmp(datatype,'Turndata_LongYy')
        culDeSacY = 0.68/3.36;
        % Define maze ROIs:
        regionsBreak_Strct.y = [-1.3, -1, -culDeSacY, culDeSacY, 1, 1.3];
        regionPrintNames_Strct.y = {'Pre-arm', 'Down-the-arm', ...
            'Cul-de-sac', 'Up-the-arm', 'Post-arm'};
        regionsBreak_Strct.time = [-1.3, -1, 0, 1, 1.3];
        regionPrintNames_Strct.time = {'Pre-arm', ...
            'Down-the-arm & c.d.s<0', 'C.d.s>0 & up-the-arm', 'Post-arm' };
        regionBreakTicks.y = -1.:.5:1;
        regionBreakTicks.time = -1.:.5:1;
        infVal = 1.5;
        % Define (new) y-bin and time-bin edges:
        yEdges = -1.3:.1:1.3;
        tEdges = -1.275:.075:1.275;
    elseif strcmp(datatype,'Turndata_ShortABM_NoWF') || ...
            strcmp(datatype,'Turndata_ShortABM_WallFollowing') || ...
            strcmp(datatype,'Turndata_ShortABM_Brownian') || ...
            strcmp(datatype,'Turndata_ShortABM_WF_00100') || ...
            strcmp(datatype,'Turndata_ShortABM_WF_00200') || ...
            strcmp(datatype,'Turndata_ShortABM_WF_00500') || ...
            strcmp(datatype,'Turndata_ShortABM_WF_00700') || ...
            strcmp(datatype,'Turndata_ShortABM_WF_00900') || ...
            strcmp(datatype,'Turndata_ShortABM_brown5')
        culDeSacY = 7.5/30;
        % Define maze ROIs:
        regionsBreak_Strct.y = [-culDeSacY, culDeSacY, 1, 1.3];
        regionPrintNames_Strct.y = {'Cul-de-sac', 'Up-the-arm', 'Post-arm'};
        regionsBreak_Strct.time = [0, 1, 1.3];
        regionPrintNames_Strct.time = {'C.d.s>0 & up-the-arm', 'Post-arm'};
        regionBreakTicks.y = regionsBreak_Strct.y(1:end-1);
        regionBreakTicks.time = regionsBreak_Strct.time(1:end-1);
        regionBreakTicks.y = 0:.5:1;
        regionBreakTicks.time = 0:.5:1;
        infVal = .2 + culDeSacY; 
        % Define (new) y-bin and time-bin edges:
        yEdges = .1*floor(-culDeSacY*10):.1:1.3; % start from -cul-de-sac
        tEdges = 0:.075:1.275; % start from t=y=0 (min. loc).
    else % human data
        regionsBreak_Strct.time = [-1.3, -1, 0, 1, 1.3];
        regionPrintNames_Strct.time = {'Pre-arm', ...
            'Down-the-arm & c.d.s<0', 'C.d.s>0 & up-the-arm', 'Post-arm'};
        regionsBreak_Strct.y = regionsBreak_Strct.time;
        regionPrintNames_Strct.y = regionPrintNames_Strct.time;
        regionBreakTicks.y = -1.:.5:1;
        regionBreakTicks.time = -1.:.5:1;
        infVal = 1.5;
        % Define (new) y-bin and time-bin edges:
        yEdges = -1.3:.1:1.3;
        tEdges = -1.275:.075:1.275;
    end
    

    % Create speed avg. figures (all datasets): 
    for bbb = 1:numel(binTypes)
        binType = binTypes{bbb};
        figSpeed.(datatype(10:end-2)).(binType) = figure;
        ax1.(binType) = axes();
        box(ax1.(binType));
        % Plot region-breaks:
        regionsBreak = regionsBreak_Strct.(binType);
    end

    binEdgesSpeed.(datatype(10:end-2)).y = yEdges; 
    binEdgesSpeed.(datatype(10:end-2)).time = tEdges;
    binEdgesSpeed.(datatype(10:end-2)).percentArm = -1.3:.1:1.3;

    % Read selected fly names:
    flyNames = fieldnames(meanX_vects_ALL);
    
    % Create structure for storing TPI's:
    savedSpeed.(datatype(10:end-2)).(paramName).y = nan( length(selectedFlies), ...
        numel(yEdges)-1 );
    savedSpeed.(datatype(10:end-2)).(paramName).time = nan( length(selectedFlies), ...
        numel(tEdges)-1 );
    savedTPIinf.(datatype(10:end-2)).(paramName).y = nan( length(selectedFlies), 2 );
    savedTPIinf.(datatype(10:end-2)).(paramName).time = nan( length(selectedFlies), 2 );

    % Create matrix for storing indiv. curve color:
    turboCols = turbo(numel(flyNames));
    turboOrds = datasample(s,1:numel(flyNames), numel(flyNames), ...
        'Replace', false );
    ppColors = turboCols(turboOrds,:);

    % Run over individual animals in dataset and compute avg. param over y 
    % and rel. time bins:

    for nF = 1:length(flyNames) % length(selectedFlies)
        locFly_InSelected = nF; 
        flyName = flyNames{nF};

        % Read individual animal's processed data for bottom trials: 
        
        % Location of bottom trials:
        locsBotTrialT = meanX_vects_ALL.(flyName).locsBotTrialT;
        % Corrected trajectories in bottom trials: 
        xytyd = meanX_vects_ALL.(flyName).xytyd;
        xytyd_bottom = xytyd( locsBotTrialT );
        

        % Run over param dimensions, w={y,time}, compute and plot param(w):

        for b = 1:numel(binTypes)
            binType = binTypes{b};

            wBinEdges = binEdgesSpeed.(datatype(10:end-2)).(binType);
            wBinCenters = .5*diff(wBinEdges(1:2)) + wBinEdges(1:end-1);
            
            % Compute avg. Param(w):
            [wAvgSpeed_fly, ~, wAvgDist_fly, wAvgRelTime_fly] = ySpeed( ...
                xytyd_bottom, wBinEdges, binType, fps );
            
            if strcmp( paramName, 'speed' )
                wAvgParam_fly = wAvgSpeed_fly;
            elseif strcmp( paramName, 'distance' )
                wAvgParam_fly = wAvgDist_fly;
            elseif strcmp( paramName, 'relTime' )
                wAvgParam_fly = wAvgRelTime_fly;
            end
    
            % Plot avg. param(w):
            figure( figSpeed.(datatype(10:end-2)).(binType) );
            patchline(  wBinCenters, wAvgParam_fly, 'linestyle', '-', ...
                'edgecolor', ppColors(nF,:), 'linewidth', 0.5, ...
                'edgealpha', 0.5 );
            hold on;

            % Save avg param(w):
            savedSpeed.(datatype(10:end-2)).(paramName).(binType)(nF,:) = wAvgParam_fly;
            
            clear wBinEdges wBinCenters;
        end

    end


    % Plot AVG. param. +- SEM (overlaid on individual param's) :

    for bb = 1:numel(binTypes)
        binType = binTypes{bb};
        binPrintName = binPrintNames{bb};
        wBinEdges = binEdgesSpeed.(datatype(10:end-2)).(binType);
        wBinCenters = .5*diff(wBinEdges(1:2)) + wBinEdges(1:end-1);

        figure( figSpeed.(datatype(10:end-2)).(binType) );
        errorbar( wBinCenters, ...
            mean( savedSpeed.(datatype(10:end-2)).(paramName).(binType), 'omitnan' ), ...
            std( savedSpeed.(datatype(10:end-2)).(paramName).(binType), 'omitnan' ) ./ ...
            sqrt( sum( ~isnan(savedSpeed.(datatype(10:end-2)).(paramName).(binType)) ) ), ...
            'Color', 'w', 'capSize', 0, 'lineWidth', 4 );
        hold on;
        errorbar( wBinCenters, ...
            mean( savedSpeed.(datatype(10:end-2)).(paramName).(binType), 'omitnan' ), ...
            std( savedSpeed.(datatype(10:end-2)).(paramName).(binType), 'omitnan' ) ./ ...
            sqrt( sum( ~isnan(savedSpeed.(datatype(10:end-2)).(paramName).(binType)) ) ), ...
            'Color', dataCol, 'capSize', 0, 'lineWidth', 2 );
        hold on;

        % axes: 
        xticks(regionBreakTicks.(binType));
        xlim( regionsBreak_Strct.(binType)([1,end]) ); 
        ggg = gca;
        ggg.XGrid = 'off';
        ggg.YGrid = 'off'; 
        xlabel( binPrintName ); 
        ylabel( ['avg. ' paramTitle ' (' binPrintName ') [a.u]'] ); 

        % Plot region breaks:
        ylims = ggg.YLim;
        regionsBreak = regionsBreak_Strct.(binType);
        regionPrintNames = regionPrintNames_Strct.(binType);
        box(ax1.(binType));
        for rr = 1:numel(regionsBreak)
            plot(regionsBreak(rr) * [1,1], ylims, 'k:' ); 
            hold on;
        end
        
        % Add region names:
        ax2 = axes('Position', get(ax1.(binType),'Position'), ...
            'XAxisLocation', 'top', 'Color', 'none', 'XColor', 'k' );
        ax2.YAxis.Visible = 'off';
        ax2.XLim = ax1.(binType).XLim;
        ax2.XTick = regionsBreak(1:end-1) + .5*diff(regionsBreak);
        ax2.XTickLabel = regionPrintNames;
        ax2.FontSize = 9;

        % Save figure;
        figStartName = ['analyses/figures/flies' datatype(10:end-2) '_avg' ...
            upper(paramTitle(1)) paramTitle(2:end) upper(binType(1))];
        saveas( figSpeed.(datatype(10:end-2)).(binType), ...
            [figStartName '.fig']);

        clear wBinEdges wBinCenters;
    end


end

% Save output parameters:
save( 'datafiles/flies_speed.mat', 'savedSpeed', 'binEdgesSpeed' );


