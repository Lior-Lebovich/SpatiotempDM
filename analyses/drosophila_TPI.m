%% SPATIAL AND TEMPORAL TURN PREDICTIVENESS INDEX (TPI): 
%
% Computes and plots TPI(y) and TPI(rel. trial time) of for all 
% decision-makers.
%
% INPUT FILES: 
%   output_drosophila_main_[DATASET_NAME(10:END)].mat
%
% FUNCTION:
%   TPI.m for computing the TPI.
%   TPI1.m for computing the TPI in human dataset (x-centering not required).
%
% EXTERNAL FUNCTION:
%   This code uses an external function, patchline.m, Reference: 
%   Brett Shoelson (2023). Patchline 
%   (https://www.mathworks.com/matlabcentral/fileexchange/36953-patchline), 
%   MATLAB Central File Exchange. Retrieved September 18, 2023.
%
% OUTPUT FILE: 
%   flies_TPIs.mat: structures of individual TPI's (over both y and time)
%       computed for selected flies in each dataset (savedTPI) and 
%       corresp. binEdges (binEdgesTPI). 
%       savedTPI.DATATYPE(10:end-1).binType is numAnimals X numBinsBinType 
%       matrix.
%
% FIGURES:
%   Spatial [or temporal] Turn-Predictiveness-Index, TPI(Y) [or TPI(T)], 
%       avg.+-SEM, overlaid on individual TPI's in dataset:
%       Figs. 1F, 1G, S11A: flies[DATATYPE(10:END-2)]_TPI_[DOMAIN].fig
%   Dataset comparison of spatial [or temporal] avg.+-SEM TPI(Y) 
%       [ TPI(T)]:
%       Figs. 4A, S3C, S4A, S8A, S8A, S13A, S15B: 
%       fliesAll_[COMPNAME]Comp_TPI[DOMAIN].fig
%   Dataset comparison of TPI(Y) cul-de-sac increase estimates:
%       Figs. S6A, S15A: fliesAll_[COMPNAME]Comp_slopeTPI[DOMAIN].fig
%   Region-aligned short/long comparison of TPI(Y):
%       Figs. 3B: fliesAll_LongAlignComp_TPIy.fig
%   Absolute size short/long comparison of TPI(Y):
%       Fig. S3B: fliesAll_LongAbsComp_TPIy.fig
%
% Note that other than in the region-aligned or absolute size short/long 
% comparisons, TPI(Y) of the long dataset appears with |y|=1 corrsp. to
% upper edge of the bottom arm (before the intersection).



%% Use relevant directory and add path to external fucntions and datafiles:

cdName = 'SpatiotempDM';
currentFolder = cd;
if ~endsWith(currentFolder,cdName)
    cd(cdName);
end
addpath( 'analyses', 'datafiles', 'analyses/customfuns', '-frozen');



%% Define datasets, col's, TPI bin types:

% Define datatypes and corresp. colors:
datatypes = {'Turndata_ShortYy', 'Turndata_LongYy', ...
    'Turndata_ShortNorpAYy', 'Turndata_ShortDumbYy', ...
    'Turndata_ShortFoxPYy', 'Turndata_ShortNompCYy', ...
    'Turndata_ShortABM_NoWF', ...
    'Turndata_ShortABM_WF_00100', 'Turndata_ShortABM_WF_00200', ...
    'Turndata_ShortABM_WF_00500', 'Turndata_ShortABM_WF_00700', ...
    'Turndata_ShortABM_WF_00900', ...
    'Turndata_ShortABM_WallFollowing', 'Turndata_ShortABM_brown5', ...
    'Turndata_HumanYy'};
dataCols = {'k', 'b', 'r', 'm', 'g', '#77AC30', ... % fly data
    "#A2142F", ... % no WF 
    "#D95319", "#EDB120", "#77AC30", "#4DBEEE", "#0072BD", ... % WF param's
    "#7E2F8E", "#7E2F8E", ... % WF, brownian
    'k'}; % human data

% Define TPI bin types
binTypes = {'y', 'time'};
binPrintNames = {'y', 'T'};

% Determine location of yInF or tInf:
infVal_global= 1.5;

% Random stream (for individual decision-makers' curves):
s = RandStream('mlfg6331_64'); 



%% Compute, store and plot TPI(bins) for all datasets:


for d = 1:numel(datatypes) 
    datatype = datatypes{d};
    dataCol = dataCols{d};
    if strcmp(datatype,'Turndata_ShortABM_brown5')
        dataCol = [.5, .3, .1];
    end
    
    % Load processed data:
    matFileName = ['output_drosophila_main_' datatype(10:end) '.mat'];
    load(matFileName);

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
    elseif strcmp(datatype,'Turndata_HumanYy') % human data WITH cul-de-sac
        culDeSacY = 0.3;
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
    end
    

    % Create TPI avg. figures (all datasets): 
    for bbb = 1:numel(binTypes)
        binType = binTypes{bbb};
        figTPI.(datatype(10:end-2)).(binType) = figure;
        ax1.(binType) = axes();
        box(ax1.(binType));
        % Plot region-breaks:
        regionsBreak = regionsBreak_Strct.(binType);
        for rr = 1:numel(regionsBreak)
            plot(regionsBreak(rr) * [1,1], [-1,1], 'k:' ); 
            hold on;
        end
    end

    binEdgesTPI.(datatype(10:end-2)).y = yEdges; 
    binEdgesTPI.(datatype(10:end-2)).time = tEdges;
    binEdgesTPI.(datatype(10:end-2)).percentArm = -1.3:.1:1.3;

    % Read selected fly names:
    flyNames = fieldnames(meanX_vects_ALL);
    
    % Create structure for storing TPI's:
    savedTPI.(datatype(10:end-2)).y = nan( length(selectedFlies), ...
        numel(yEdges)-1 );
    savedTPI.(datatype(10:end-2)).time = nan( length(selectedFlies), ...
        numel(tEdges)-1 );
    savedTPIinf.(datatype(10:end-2)).y = nan( length(selectedFlies), 2 );
    savedTPIinf.(datatype(10:end-2)).time = nan( length(selectedFlies), 2 );
    savedTPI.(datatype(10:end-2)).percentArm = nan( ...
        length(selectedFlies), ...
        numel(binEdgesTPI.(datatype(10:end-2)).percentArm)-1 );
    
    % Create vectors for "bottom" arm yEdges w/o intersection (for 
    % binTypes='percentArm'):
    strcBotArmLenY.(datatype(10:end-2)).wInsctn = ...
        nan( length(selectedFlies), 1 );
    strcBotArmLenY.(datatype(10:end-2)).woInsctn = ...
        nan( length(selectedFlies), 1 );

    % Create matrix for storing indiv. curve color:
    turboCols = turbo(numel(flyNames));
    turboOrds = datasample(s,1:numel(flyNames), numel(flyNames), ...
        'Replace', false );
    ppColors = turboCols(turboOrds,:);

    % Run over individual animals in dataset and compute TPI of y and time:

    for nF = 1:length(flyNames) % length(selectedFlies)
        locFly_InSelected = nF; 
        % locFly_InAllFlies = selectedFlies(locFly_InSelected); 
        % flyName = ['f' num2str(locFly_InAllFlies)];
        flyName = flyNames{nF};

        % Read individual animal's processed data for bottom trials: 
        
        % Location of bottom trials:
        locsBotTrialT = meanX_vects_ALL.(flyName).locsBotTrialT;
        % Corrected trajectories in bottom trials: 
        xytyd = meanX_vects_ALL.(flyName).xytyd;
        xytyd_bottom = xytyd( locsBotTrialT );
        % Decision vector for bototm trials:
        Dec = meanX_vects_ALL.(flyName).Dec;
        Dec_bottom = Dec( locsBotTrialT );

        % Compute by-subject "bottom" arm yEdges w/o the intersection (for 
        % binTypes='percentArm'):
        xytyd_bottom_mat = cell2mat( xytyd_bottom );
        selectedSubX = xytyd_bottom_mat(:,1);
        selectedSubY = xytyd_bottom_mat(:,2); 
        

        % Run over TPI dimensions, w={y,time}, compute and plot TPI(w):

        for b = 1:numel(binTypes)
            binType = binTypes{b};

            wBinEdges = binEdgesTPI.(datatype(10:end-2)).(binType);
            wBinCenters = .5*diff(wBinEdges(1:2)) + wBinEdges(1:end-1);
            
            % Compute TPI(w):
            if strcmp(datatype,'Turndata_HumanYy')
                [TPI_w_fly,~, ~, ~, TPI_w_fly_mpInf, ~] = TPI1( ...
                    xytyd_bottom, wBinEdges, Dec_bottom, binType );
            else
                [TPI_w_fly,~, ~, ~, TPI_w_fly_mpInf, ~] = TPI( ...
                    xytyd_bottom, wBinEdges, Dec_bottom, binType );
            end
    
            % Plot TPI(w):
            figure( figTPI.(datatype(10:end-2)).(binType) );
            patchline(  wBinCenters, TPI_w_fly, 'linestyle', '-', ...
                'edgecolor', ppColors(nF,:), 'linewidth', 0.5, ...
                'edgealpha', 0.5 );
            hold on;
            % Add TPI(+-wInf):
            patchline( [-infVal,infVal_global], TPI_w_fly_mpInf, ...
                'linestyle', 'none', 'marker', 'o', ...
                'edgecolor', ppColors(nF,:), 'linewidth', 0.5, ...
                'edgealpha', 0.5 ); 
            hold on;

            % Save TPI(w):
            savedTPI.(datatype(10:end-2)).(binType)(nF,:) = TPI_w_fly;
            savedTPIinf.(datatype(10:end-2)).(binType)(nF,:) = ...
                TPI_w_fly_mpInf;
            
            clear wBinEdges wBinCenters;
        end

    end


    % Plot avg. TPI over dataset +- SEM (overlaid on individual TPI's) :

    for bb = 1:numel(binTypes)
        binType = binTypes{bb};
        binPrintName = binPrintNames{bb};
        wBinEdges = binEdgesTPI.(datatype(10:end-2)).(binType);
        wBinCenters = .5*diff(wBinEdges(1:2)) + wBinEdges(1:end-1);

        figure( figTPI.(datatype(10:end-2)).(binType) );
        errorbar( wBinCenters, ...
            mean( savedTPI.(datatype(10:end-2)).(binType), 'omitnan' ), ...
            std( savedTPI.(datatype(10:end-2)).(binType), 'omitnan' ) ./ ...
            sqrt( sum( ~isnan(savedTPI.(datatype(10:end-2)).(binType)) ) ), ...
            'Color', 'w', 'capSize', 0, 'lineWidth', 4 );
        hold on;
        errorbar( wBinCenters, ...
            mean( savedTPI.(datatype(10:end-2)).(binType), 'omitnan' ), ...
            std( savedTPI.(datatype(10:end-2)).(binType), 'omitnan' ) ./ ...
            sqrt( sum( ~isnan(savedTPI.(datatype(10:end-2)).(binType)) ) ), ...
            'Color', dataCol, 'capSize', 0, 'lineWidth', 2 );
        hold on;
        
        % Add +-Inf:
        errorbar(  [-infVal,infVal_global], ...
            mean( savedTPIinf.(datatype(10:end-2)).(binType), 'omitnan' ), ...
            std( savedTPIinf.(datatype(10:end-2)).(binType), 'omitnan' ) ./ ...
            sqrt( sum( ~isnan(savedTPIinf.(datatype(10:end-2)).(binType)) ) ), ...
            'Color', dataCol, 'capSize', 0, 'lineWidth', 2, ...
            'LineStyle', 'none', 'Marker', 'o' );

        % axes:
        xticks([-infVal, regionBreakTicks.(binType), infVal_global]);
        midLabels = arrayfun(@num2str, regionBreakTicks.(binType), ...
            'UniformOutput', false);
        firstLabel = ['-' binPrintName '_{\infty}'];
        lastLabel = [binPrintName '_{\infty}'];
        xticklabels([firstLabel, midLabels,lastLabel]);
        xlim( [-infVal, infVal_global] ); 
        ylim([-1,1]);
        yticks(-1:.2:1);
        ggg = gca;
        ggg.XGrid = 'off';
        ggg.YGrid = 'on'; 
        xlabel( binPrintName ); 
        ylabel( ['TPI(' binPrintName ')'] ); 
        
        % Add region names:
        regionsBreak = regionsBreak_Strct.(binType);
        regionPrintNames = regionPrintNames_Strct.(binType);
        ax2 = axes('Position', get(ax1.(binType),'Position'), ...
            'XAxisLocation', 'top', 'Color', 'none', 'XColor', 'k' );
        ax2.YAxis.Visible = 'off';
        ax2.XLim = ax1.(binType).XLim;
        ax2.XTick = regionsBreak(1:end-1) + .5*diff(regionsBreak);
        ax2.XTickLabel = regionPrintNames;
        ax2.FontSize = 9;

        % Save figure;
        figStartName = ['analyses/figures/flies' datatype(10:end-2) '_TPI' ...
            upper(binType(1))];
        saveas( figTPI.(datatype(10:end-2)).(binType), ...
            [figStartName '.fig']);

        clear wBinEdges wBinCenters;
    end
    clearvars -except d cdName datatypes dataCols figTPI binEdgesTPI ...
        savedTPI savedTPIinf binTypes binPrintNames infVal ...
        regionPrintNames_Strct regionsBreak_Strct infVal_global s;
end


% Save output TPI's:
save( 'datafiles/flies_TPIs.mat', 'savedTPI', 'binEdgesTPI', 'savedTPIinf' );



%% Compare datasets' AVG.+-SEM TPI(y/time) AND cul-de-sac increase est's:

% Note that short/long comp. plots long in norm. rather than 
% region-aligned. For region-alignment - see next section. For comp. in 
% abs. size the section after.

% Define dataset comparisons:
compNames = {'Short', 'LongNorm', 'Models'};
compLoc.Short = [1, 5, 4, 3, 6]; % short fly comparison
compLoc.Models = [14, 7:13]; % model comparison
compLoc.LongNorm = 1:2; % short/long Norm. comparison

infVal_global = 1.5;
infVal = infVal_global;

for c = 1:numel(compNames)
    compName = compNames{c};
    compLocs = compLoc.(compName);
    datatypes_c = datatypes(compLocs);
    dataCols_c = dataCols(compLocs);
    
    for b = 1:numel(binTypes)
        binType = binTypes{b};
        binPrintName = binPrintNames{b};
    
        % Create avg. TPI figure:
        figTPI.all.(binType) = figure;

        % Also create TPI(y) slope-comp. fig. for all but human datasets:
        if strcmp(binType,'y') && ~strcmp(compName,'Human')
            figSlopeTPI.all.(binType) = figure;
            nComps = numel(datatypes_c)-1; % #comparisons for sig.
            pValsSlopeComps = nan(3,nComps+1);
        end

        
        % Initialize legends:
        clear legs legNames;
        
        % Plot avg. TPI for each datatype:
        for d = 1:numel(datatypes_c)
            datatype = datatypes_c{d};
            dataCol = dataCols_c{d};
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
                regionPrintNames_Strct.y = {'Cul-de-sac', 'Up-the-arm', ...
                    'Post-arm'};
                regionsBreak_Strct.time = [0, 1, 1.3];
                regionPrintNames_Strct.time = {'C.d.s>0 & up-the-arm', ...
                    'Post-arm'};
                infVal = .2 + culDeSacY; 
                regionBreakTicks.y = 0:.5:1;
                regionBreakTicks.time = 0:.5:1;
            else % human data
                regionsBreak_Strct.time = [-1.3, -1, 0, 1, 1.3];
                regionPrintNames_Strct.time = {'Pre-arm', ...
                    'Down-the-arm & c.d.s<0', 'C.d.s>0 & up-the-arm', ...
                    'Post-arm'};
                regionsBreak_Strct.y = regionsBreak_Strct.time;
                regionPrintNames_Strct.y = regionPrintNames_Strct.time;
                regionBreakTicks.y = -1.:.5:1;
                regionBreakTicks.time = -1.:.5:1;
                infVal = 1.5;
            end

            if strcmp(datatype,'Turndata_ShortYy')
                legNames{d} = 'WT';
            elseif strcmp(datatype,'Turndata_LongYy')
                legNames{d} = 'WT long';
            elseif strcmp(datatype(end-1:end),'Yy') || ...
                    strcmp(datatype(10:14),'human')
                legNames{d} = datatype(15:end-2);
            elseif strcmp(datatype(15:21),'ABM_WF_')
                legNames{d} = ['WF' datatype(22:24)];
            elseif strcmp(datatype,'Turndata_ShortABM_NoWF')
                legNames{d} = 'No WF';
            elseif strcmp(datatype,'Turndata_ShortABM_WallFollowing')
                legNames{d} = 'WF';
            elseif strcmp(datatype,'Turndata_ShortABM_brown5')
                legNames{d} = 'Brownian';
            end
    

            % TPI FIGURE:

            figure(figTPI.all.(binType));

            % Plot region breaks:
            regionsBreak = regionsBreak_Strct.(binType);
            regionPrintNames = regionPrintNames_Strct.(binType);
            if d==1
                ax1 = axes();
                box(ax1);
                
                for rr = 1:numel(regionsBreak)
                    plot(regionsBreak(rr) * [1,1], [-1,1], 'k:' ); 
                    hold on;
                end
            end
            % Also plot Long region breaks if Long/short comparsion:
            if strcmp(datatype,'Turndata_LongYy') && ...
                    strcmp(compName,'LongNorm')
                for rr = 1:numel(regionsBreak)
                    plot(regionsBreak(rr) * [1,1], [-1,1], 'b:' ); 
                    hold on;
                end
            end
    
            % Define edges/centers:
            wBinEdges = binEdgesTPI.(datatype(10:end-2)).(binType);
            wBinCenters = .5*diff(wBinEdges(1:2)) + wBinEdges(1:end-1);
            errorbar( wBinCenters, ...
                mean( savedTPI.(datatype(10:end-2)).(binType), 'omitnan' ), ...
                std( savedTPI.(datatype(10:end-2)).(binType), 'omitnan' ) ./ ...
                sqrt( sum( ~isnan(savedTPI.(datatype(10:end-2)).(binType)) ) ), ...
                'Color', 'w', 'capSize', 0, 'lineWidth', 4 );
            hold on;
            legs(d) = errorbar( wBinCenters, ...
                mean( savedTPI.(datatype(10:end-2)).(binType), 'omitnan' ), ...
                std( savedTPI.(datatype(10:end-2)).(binType), 'omitnan' ) ./ ...
                sqrt( sum( ~isnan(savedTPI.(datatype(10:end-2)).(binType)) ) ), ...
                'Color', dataCol, 'capSize', 0, 'lineWidth', 2 );
            hold on;
            % Plot avg. TPI(+-Inf):
            errorbar(  [-infVal,infVal_global], ...
                mean( savedTPIinf.(datatype(10:end-2)).(binType), 'omitnan' ), ...
                std( savedTPIinf.(datatype(10:end-2)).(binType), 'omitnan' ) ./ ...
                sqrt( sum( ~isnan(savedTPIinf.(datatype(10:end-2)).(binType)) ) ), ...
                'Color', dataCol, 'capSize', 0, 'lineWidth', 2, ...
                'LineStyle', 'none', 'Marker', 'o' );


            
            % TPI SLOPE FIGURE:

            % Also compute slopes in cul-de-sac:

            if strcmp(binType,'y') && ~strcmp(compName,'Human')

                figure(figSlopeTPI.all.(binType));

                % Read first/last TPI(y) values within cul-de-sac and -Inf:
                lastLocInCul = find( wBinCenters <= culDeSacY, 1, 'last' );
                firstLocInCul = find( wBinCenters >= -culDeSacY, 1, ...
                    'first' );
                TPI_lastInCul = savedTPI.(datatype(10:end-2)).( ...
                    binType)( :,lastLocInCul);
                TPI_firstInCul = savedTPI.(datatype(10:end-2)).( ...
                    binType)( :,firstLocInCul);
                TPI_minusInf = savedTPIinf.(datatype(10:end-2)).( ...
                    binType)( :,1);

                % Plot box chart (current dataset) of slope estimate:
                subplot(3,1,1);
                boxchart( d*ones(size(TPI_lastInCul)), TPI_lastInCul, ...
                    'BoxFaceColor', dataCol, 'MarkerColor', dataCol );
                hold on;
                subplot(3,1,2);
                boxchart( d*ones(size(TPI_lastInCul)), ...
                    TPI_lastInCul - abs( TPI_firstInCul ), ...
                    'BoxFaceColor', dataCol, 'MarkerColor', dataCol );
                hold on;
                subplot(3,1,3);
                boxchart( d*ones(size(TPI_lastInCul)), ...
                    TPI_lastInCul - abs( TPI_minusInf ), ...
                    'BoxFaceColor', dataCol, 'MarkerColor', dataCol );
                hold on;

                % Add sig. of Wilcoxon rank sum test (fly data) with WT:
                if d==1 && (strcmp(compName, 'Short') || ...
                        strcmp(compName, 'LongNorm'))
                    TPI_lastInCul_WT = TPI_lastInCul;
                    TPI_firstInCul_WT = TPI_firstInCul;
                    TPI_minusInf_WT = TPI_minusInf; 
                elseif d>1 && (strcmp(compName, 'Short') || ...
                        strcmp(compName, 'LongNorm'))
                    pVal1 = ranksum(TPI_lastInCul,TPI_lastInCul_WT);
                    pVal2 = ranksum(TPI_lastInCul - abs( TPI_firstInCul ),...
                        TPI_lastInCul_WT - abs( TPI_firstInCul_WT ));
                    pVal3 = ranksum(TPI_lastInCul - abs( TPI_minusInf ),...
                        TPI_lastInCul_WT - abs( TPI_minusInf_WT ));
                    pValsSlopeComps(1,d) = pVal1;
                    pValsSlopeComps(2,d) = pVal2;
                    pValsSlopeComps(3,d) = pVal3;
                end
            end
            
        end
    
        % TPI figure axes: 
        figure(figTPI.all.(binType));
        xticks([-infVal, regionBreakTicks.(binType), infVal_global]);
        midLabels = arrayfun(@num2str, regionBreakTicks.(binType), ...
            'UniformOutput', false);
        firstLabel = ['-' binPrintName '_{\infty}'];
        lastLabel = [binPrintName '_{\infty}'];
        xticklabels([firstLabel, midLabels,lastLabel]);
        ylim([-1,1]);
        xlim( [-infVal, infVal_global] ); 
        yticks(-1:.2:1);
        ggg = gca;
        ggg.XGrid = 'off';
        ggg.YGrid = 'on'; 
        xlabel( binPrintName ); 
        ylabel( ['TPI(' binPrintName ')'] ); 
        legend( legs, legNames, 'Location', 'NorthWest' );
        % Add region names: y-regions are only relevant for SHORT fly 
        % datasets; Time-regions are relevant for all datasets:
        if d == numel(datatypes_c)
            ax2 = axes('Position', get(ax1,'Position'), 'XAxisLocation', ...
                'top', 'Color', 'none', 'XColor', 'k' );
            ax2.YAxis.Visible = 'off';
            ax2.XLim = ax1.XLim;
            ax2.XTick = regionsBreak(1:end-1) + .5*diff(regionsBreak);
            ax2.XTickLabel = regionPrintNames;
            ax2.FontSize = 9;
        end
        % Save figure: 
        figStartName = ['analyses/figures/fliesAll_' compName 'Comp_TPI' ...
            upper(binType(1))];
        saveas( figTPI.all.(binType), [figStartName '.fig']);


        % TPI SLOPE figure axes: 
        if strcmp(binType,'y') && ~strcmp(compName,'Human')
            figure(figSlopeTPI.all.(binType));
            for sp = 1:3
                subplot(3,1,sp);
                plot([0,numel(datatypes_c)+1],[0,0],'k:');
                hold on;
                % Add significance of Wilcoxon rank sum test:
                if strcmp(compName,'Short') || strcmp(compName,'LongNorm')
                    maxYForSig = 1.01;
                    pVals = pValsSlopeComps(sp,:);
                    for dd = 2:numel(datatypes_c)
                        if pVals(dd)<=0.05/nComps
                            maxYForSig = maxYForSig + .2;
                            plot( [1, dd], maxYForSig * [1,1], 'k-', ...
                                (1+dd)/2, .005 + maxYForSig, '*k' );
                            hold on;
                        end
                    end
                end
                xticks(1:numel(datatypes_c));
                yticks(-1:.5:1);
                ylim([-1,.1 + 1e-1*ceil(1e1*maxYForSig)]); 
                xticklabels( legNames );
                xtickangle(0);
            end
            subplot(3,1,1); 
            ylabel('TPI(y|c.d.s end)');
            subplot(3,1,2); 
            ylabel('TPI(y|c.d.s end) - |TPI(y|c.d.s start)|');
            subplot(3,1,3); 
            ylabel('TPI(y|c.d.s end) - |TPI(-y_{\infty})|');
            % Save figure:
            figStartName = ['analyses/figures/fliesAll_' compName ...
                'Comp_slopeTPI' upper(binType(1))];
            saveas( figSlopeTPI.all.(binType), [figStartName '.fig']);
        end
        
    end
    
end



%% Compare short/long TPI(Y) - REGIONS-ALIGNMENT (of long):

% [Long is aligned to short]

clear figTPI;
figTPI = figure;

ax1 = axes();
box(ax1);

% Define regions (corresp. short):
culDeSacY = 0.69/2.03;
regionsBreak = [-1.3, -1, -culDeSacY, culDeSacY, 1, 1.3];
for rr = 1:numel(regionsBreak)
    plot(regionsBreak(rr) * [1,1], [-1,1], 'k:' ); 
    hold on;
end

infVal_global = 1.5;
infVal = infVal_global;

datatypes = {'Turndata_ShortYy', 'Turndata_LongYy'};
dataCols = {'k','b'};
regionNames = {'beforeBottom', 'downArm', 'culdesac', 'upArm', ...
    'afterBottom'};
regionPrintNames = {'Pre-arm', 'Down-the-arm', 'Cul-de-sac', ...
    'Up-the-arm', 'Post-arm'};
edgesy.Short.beforeBottom = -1.3:.1:-1;
edgesy.Short.downArm = linspace(-1,-culDeSacY,7);
edgesy.Short.culdesac = linspace(-culDeSacY,culDeSacY,5);
edgesy.Short.upArm = linspace(culDeSacY,1,7);
edgesy.Short.afterBottom = 1:.1:1.3;

ratioForLong = load( 'output_drosophila_main_LongYy.mat', 'ratioSize' );
ratioLong = ratioForLong.ratioSize;

edgesy.Long.beforeBottom = -(1+(0.3*ratioLong)):(.1*ratioLong):-1;
edgesy.Long.downArm = linspace(-1,-culDeSacY*ratioLong,11);
edgesy.Long.culdesac = ratioLong * linspace(-culDeSacY,culDeSacY,5);
edgesy.Long.upArm = linspace(culDeSacY*ratioLong,1,11);
edgesy.Long.afterBottom = 1:(.1*ratioLong):(1+(0.3*ratioLong));

longPostCrrctn.beforeBottom = 1 / ratioLong;
longPostCrrctn.downArm = 1;
longPostCrrctn.culdesac = 1 / ratioLong;
longPostCrrctn.upArm = 1;
longPostCrrctn.afterBottom = 1 / ratioLong;

for d = 1:numel(datatypes) 
    datatype = datatypes{d};
    datatype0 = datatype(10:end-2);
    dataCol = dataCols{d};
    % Load data:
    matFileName = ['output_drosophila_main_' datatype(10:end) '.mat'];
    load(matFileName);
    % Read selected fly names:
    flyNames = fieldnames(meanX_vects_ALL);
    for r = 1:numel(regionNames)
        regionName = regionNames{r};
        yEdges = edgesy.(datatype0).(regionName);
        if strcmp(datatype0,'Long')
            yEdgesCrrtn = longPostCrrctn.(regionName);
        else
            yEdgesCrrtn = 1;
        end
        if strcmp(datatype0,'Long') && ( ...
                strcmp(regionName,'beforeBottom') || ...
                strcmp(regionName,'afterBottom') )
            theShort = edgesy.Short.(regionName);
            yCentersFinal = theShort(1:end-1) + .5 * diff(theShort(1:2));
        elseif strcmp(datatype0,'Long') && ( ...
                strcmp(regionName,'downArm') || ...
                strcmp(regionName,'upArm') )
            theShort = edgesy.Short.(regionName);
            theShort2 = linspace( theShort(1), theShort(end), numel(yEdges) );
            yCentersFinal = theShort2(1:end-1) + .5 * diff(theShort2(1:2));
        else
            yCentersFinal = yEdgesCrrtn * ( yEdges(1:end-1) + .5* diff(yEdges(1:2)) );
        end
        
        % Run over flies in data * region --> compute TPI:
        TPI_region = nan( numel(selectedFlies), numel(yEdges)-1 );
        TPI_regionInf = nan( numel(selectedFlies), 1 );
        

        for nF = 1:length(flyNames) % length(selectedFlies)
            locFly_InSelected = nF; 
            flyName = flyNames{nF};
    
            % Read individual animal's processed data for bottom trials: 
            
            % Location of bottom trials:
            locsBotTrialT = meanX_vects_ALL.(flyName).locsBotTrialT;
            % Corrected trajectories in bottom trials: 
            xytyd = meanX_vects_ALL.(flyName).xytyd;
            xytyd_bottom = xytyd( locsBotTrialT );
            % Decision vector for bototm trials:
            Dec = meanX_vects_ALL.(flyName).Dec;
            Dec_bottom = Dec( locsBotTrialT );
    
            % Compute by-subject "bottom" arm yEdges w/o the intersection (for 
            % binTypes='percentArm'):
            xytyd_bottom_mat = cell2mat( xytyd_bottom );
            selectedSubX = xytyd_bottom_mat(:,1);
            selectedSubY = xytyd_bottom_mat(:,2); 
            [TPI_w_fly,~, ~, ~, TPI_w_fly_mpInf, ~] = TPI( ...
                    xytyd_bottom, yEdges, Dec_bottom, 'y' );
            
            % Store fly TPI:
            TPI_region(nF,:) = TPI_w_fly;
            if strcmp(regionName,'beforeBottom')
                TPI_regionInf(nF) = TPI_w_fly_mpInf(1);
            elseif strcmp(regionName,'afterBottom')
                TPI_regionInf(nF) = TPI_w_fly_mpInf(2);
            end
            
        end

        % Plot TPI:
        errorbar( yCentersFinal, ...
            mean( TPI_region, 'omitnan' ), ...
            std( TPI_region, 'omitnan' ) ./ ...
            sqrt( sum( ~isnan(TPI_region) ) ), ...
            'Color', dataCol, 'capSize', 0, 'lineWidth', 2 );
        hold on;
    
        % Plot TPI +-yinf:
        if strcmp(regionName,'beforeBottom') || ...
                strcmp(regionName,'afterBottom')
            errorbar( sign(yEdges(1)) * 1.5, ...
                mean( TPI_regionInf, 'omitnan' ), ...
                std( TPI_regionInf, 'omitnan' ) ./ ...
                sqrt( sum( ~isnan(TPI_regionInf) ) ), ...
                'Color', dataCol, 'capSize', 0, 'lineWidth', 2, ...
                'marker', 'o' );
            hold on;
        end

    end

end

% Axis:
ylim([-1,1]);
xlabel('y'); 
ylabel('TPI(y)');
xticks([-infVal, -1:.5:1, infVal]);
xticklabels({['-' binPrintName '_{\infty}'], '-1', '-0.5', '0', ...
    '0.5', '1', [binPrintName '_{\infty}']});
yticks(-1:.2:1);
ggg = gca;
ggg.XGrid = 'off';
ggg.YGrid = 'on';  

% Add another x-axis --> add region names:
ax2 = axes('Position', get(ax1,'Position'), ...
    'XAxisLocation','top', ...
    'Color','none', ...
    'XColor','k');
ax2.YAxis.Visible = 'off';
ax2.XLim = ax1.XLim;
ax2.XTick = regionsBreak(1:end-1) + .5*diff(regionsBreak);
ax2.XTickLabel = regionPrintNames;
ax2.FontSize = 9;

% Save figure:
figStartName = 'analyses/figures/fliesAll_LongAlignComp_TPIy';
saveas( figTPI, [figStartName '.fig']);



%% Also Plot TPI(y) of short and Long in ABSOLUTE SIZE


% Compute TPI(y) for Long in ABSOLUTE SIZE:

% Load TPI(y) computed earlier (in drosophila_TPI.m):
load('flies_TPIs.mat');

% Note that the above TPI values correspond to |y|=1 <--> bottom arm edge.

% Compute TPI(y) for LONG in ABSOLUTE SIZE:

datatype = 'Turndata_LongYy';
dataSaveName = 'LongAbs';
binType = 'y';
binPrintName = 'y';
matFileName = ['output_drosophila_main_' datatype(10:end) '.mat'];
load(matFileName);

% Define abs. size yEdges:
rel_yEdges = yEdges;
delta_rel_yEdges = diff(yEdges(1:2));
abs_yEdges_pos = delta_rel_yEdges:delta_rel_yEdges:1.9;
abs_yEdges = [-1*flip(abs_yEdges_pos), 0, abs_yEdges_pos];
yEdges_final_long = abs_yEdges;
binEdgesTPI.(dataSaveName).(binType) = abs_yEdges * ratioSize;

% Read selected fly names:
flyNames = fieldnames(meanX_vects_ALL);
    
% Create structure for storing TPI's:
savedTPI.(dataSaveName).(binType) = nan( length(selectedFlies), ...
    numel(binEdgesTPI.(dataSaveName).(binType))-1 );
savedTPIinf.(dataSaveName).(binType) = nan( length(selectedFlies), 2 );

for nF = 1:length(flyNames) % length(selectedFlies)
    locFly_InSelected = nF; 
    flyName = flyNames{nF};
    
    % Read individual animal's processed data for bottom trials: 
    
    % Location of bottom trials:
    locsBotTrialT = meanX_vects_ALL.(flyName).locsBotTrialT;
    % Corrected trajectories in bottom trials: 
    xytyd = meanX_vects_ALL.(flyName).xytyd;
    xytyd_bottom = xytyd( locsBotTrialT );
    % Decision vector for bototm trials:
    Dec = meanX_vects_ALL.(flyName).Dec;
    Dec_bottom = Dec( locsBotTrialT );

    wBinEdges = binEdgesTPI.(dataSaveName).(binType);
    wBinCenters = .5*diff(wBinEdges(1:2)) + wBinEdges(1:end-1);

    % Compute TPI(yAbs):
    [TPI_w_fly,~, ~, ~, TPI_w_fly_mpInf, ~] = TPI( xytyd_bottom, ...
        wBinEdges, Dec_bottom, 'y' );

    % Save TPI(w):
    savedTPI.(dataSaveName).(binType)(nF,:) = TPI_w_fly;
    savedTPIinf.(dataSaveName).(binType)(nF,:) = TPI_w_fly_mpInf;

end
binEdgesTPI.(dataSaveName).(binType) = yEdges_final_long;


% Plot short vs. long in ABSOLUTE SIZE:

datatypes = {'Turndata_ShortYy', 'Turndata_LongAbsYy'};
dataCols ={'k', 'b'};
dataRatios = [1,ratioSize];
% Use short cul-de-sac:
culDeSacY = 0.69/2.03; 
% Define short region-breaks
regionsBreak = [-1.3, -1, -culDeSacY, culDeSacY, 1, 1.3];
infVal = 2.3;
% Initialize legends and figure:
clear legs legNames figTPI;
figTPI = figure;
ax1 = axes();
box(ax1);

% Plot region breaks:
for d = 1:numel(datatypes)    
    for rr = 1:numel(regionsBreak)
        plot( ( regionsBreak(rr) / dataRatios(d) ) * [1,1], [-1,1], ...
            [dataCols{d} ':'] ); 
        hold on;
    end
end

% Plot avg. TPI(y) in Abs. values: 
for d = 1:numel(datatypes)    
    datatype = datatypes{d};
    dataCol = dataCols{d};
    % Define edges/centers:
    wBinEdges = binEdgesTPI.(datatype(10:end-2)).(binType);
    wBinCenters = .5*diff(wBinEdges(1:2)) + wBinEdges(1:end-1);
    errorbar( wBinCenters, ...
        mean( savedTPI.(datatype(10:end-2)).(binType), 'omitnan' ), ...
        std( savedTPI.(datatype(10:end-2)).(binType), 'omitnan' ) ./ ...
        sqrt( sum( ~isnan(savedTPI.(datatype(10:end-2)).(binType)) ) ), ...
        'Color', 'w', 'capSize', 0, 'lineWidth', 4 );
    hold on;
    legs(d) = errorbar( wBinCenters, ...
        mean( savedTPI.(datatype(10:end-2)).(binType), 'omitnan' ), ...
        std( savedTPI.(datatype(10:end-2)).(binType), 'omitnan' ) ./ ...
        sqrt( sum( ~isnan(savedTPI.(datatype(10:end-2)).(binType)) ) ), ...
        'Color', dataCol, 'capSize', 0, 'lineWidth', 2 );
    hold on;
    legNames{d} = datatype(10:end-2);
    % Plot avg. TPI(+-Inf):
    errorbar(  infVal * [-1,1], ...
        mean( savedTPIinf.(datatype(10:end-2)).(binType), 'omitnan' ), ...
        std( savedTPIinf.(datatype(10:end-2)).(binType), 'omitnan' ) ./ ...
        sqrt( sum( ~isnan(savedTPIinf.(datatype(10:end-2)).(binType)) ) ), ...
        'Color', dataCol, 'capSize', 0, 'lineWidth', 2, ...
        'LineStyle', 'none', 'Marker', 'o' );
end

% axes: 
ylim([-1,1]);
xticks([-infVal, -1.5:.5:1.5, infVal]);
xticklabels({['-' binPrintName '_{\infty}'], '-1.5', '-1', '-0.5', ...
    '0', '0.5', '1', '1.5', [binPrintName '_{\infty}']});
xlim( infVal * [-1,1] ); 
yticks(-1:.2:1);
ggg = gca;
ggg.XGrid = 'off';
ggg.YGrid = 'on'; 
xlabel( binPrintName ); 
ylabel( ['TPI(' binPrintName ')'] ); 
legend( legs, legNames, 'Location', 'NorthWest' );

% Save figure:
figStartName = 'analyses/figures/fliesAll_LongAbsComp_TPIy';
saveas( figTPI, [figStartName '.fig']);


