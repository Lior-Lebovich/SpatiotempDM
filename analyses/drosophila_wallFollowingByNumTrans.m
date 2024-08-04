%% CROSSINGS-DEPENDENT GLOBAL LATERAL TENDENCIES AND SPATIAL TPI IN SHORT MAZE:
%
% Computes and plots MAD around the center as a function of 
% #midline-crossings post-cul-de-sac. Computes and plots spatial TPI and 
% estimates for TPI increase for trials with non-zero post-cul-de-sac 
% crossings.
%
% INPUT FILES: 
%   output_drosophila_main_[DATASET_NAME(10:END)].mat
%   flies_pEvenTrans.mat
%
% FUNCTION:
%   TPI.m for computing the TPI.
%
% OUTPUT FILE: 
%   flies_TPIs_noZero.mat: TPI and bin structures, after exclusion of
%       zero post cul-de-sac midline-crossings trials.
%
% Figures:
%   Datasets Avg.+-SEM MAD as a function of #midline-crossings post 
%       cul-de-sac:
%       Fig. 4B: flies_MAD_byTranNum_shortComp_upArmWalk.fig
%   Datasets Avg.+-SEM TPI(y) for non-zero post cul-de-sac #crossings
%       trials:
%       Fig. S9C: fliesAll_ShortComp_TPIY_noZeroTrans.fig
%   Dataset comparison of TPI(Y) cul-de-sac increase estimates for 
%       non-zero post cul-de-sac #crossings trials:
%       Fig. S9D: fliesAll_ShortComp_slopeTPIY_noZeroTrans.fig
%
% Note in the current midline-crossings are termed "transitions".



%% Use relevant directory and add path to external fucntions and datafiles:

cdName = 'SpatiotempDM';
currentFolder = cd;
if ~endsWith(currentFolder,cdName)
    cd(cdName);
end
addpath( 'analyses', 'datafiles', 'analyses/customfuns', '-frozen');



%% Define datasets, arm directions, cul-de-sac loc., load transition data:
% exist. Otherwise run codes producing them:

dataTypes_stable = {'Turndata_ShortYy', 'Turndata_ShortFoxPYy', ...
    'Turndata_ShortDumbYy', 'Turndata_ShortNorpAYy', ...
    'Turndata_ShortNompCYy'};
dataCols_stable = {'k', 'g', 'm', 'r', '#77AC30' }; 

dataTypes = dataTypes_stable; 
dataCols = dataCols_stable;

% Define upper edge of cul-de-sac (same for all short fly mazes):
for d = 1:length(dataTypes)
    dataType = dataTypes{d};
    culDeSacEdge.(dataType) = 0.69/2.03;
end

% Define walking direction:
armDirNames = {'up'}; % names for walk directions
armDirLocs = [1, 2]; % edge def's for walk directions

namesMAD = {'MAD around center'}; % MAD types 

% Define transitions numbers for which the MAD will be computed: 
transNumsVect = 0:10;

% Load transition data for each trial, of each fly, in each dataset:
load( 'flies_pEvenTrans.mat', 'numTrans' );



%% Computes and plots MAD|#trans.=k, computes TPI|#trans.~=0:

% Computes and plots MAD around the center as a function of #Transitions.
% Computes TPI for non-zero transition trial.

% The MAD computation below is similiar to the one in 
% drosophila_wallFollowing.m, only that:
% (1) the MAD is computed separately for each #transitions=k.
% (2) MAD(x) is computed here ONLY as median abs. dev. around the center.
% (3) MAD(x) is computed for trajectories within the bottom arm - between 
    % the cul-de-sac and the intersection - ONLY for 'up' direction: when 
    % the fly walks up the arm (toward intersection).

% Note that for the x-centering, we use the entire data, but the actual MAD 
% is computed separately for each #trans=k.
% As well, flyTrlsXDir here is a cell rather than a vector.

% For more details see: drosophila_wallFollowing.m.

kSubSamp = 1; % subsampling frames (1= no subsamling)


% Create figures for dataset and averages:
figDataset = figure;
figAvgs = figure;


% Compute MADs for each dataset, fly and arm walking dir.:

for d = 1:length(dataTypes)

    d_stable = d;
    dataType = dataTypes{d};
    matFileName = ['output_drosophila_main_' dataType(10:end) '.mat'];

    % Load dataset data:
    load(matFileName);
    % Keep dataType and Cols:
    d = d_stable;
    dataTypes = dataTypes_stable; 
    dataCols = dataCols_stable;
    dataType = dataTypes{d};
   
    % Create structure for storing flies' TPI without zero transitions:
    yEdges = -1.3:.1:1.3;
    binEdgesTPI.(dataType(10:end-2)).y = yEdges; 
    savedTPI.(dataType(10:end-2)).y = nan( length(selectedFlies), ...
        numel(yEdges)-1 );
    savedTPIinf.(dataType(10:end-2)).y = nan( length(selectedFlies), 2 );

    % Load arm-edge(y) def. (MAD's will be based on coordinates within 
    % armEdges):
    armTopEnd = 1; 
    armEdges = [culDeSacEdge.(dataType), armTopEnd]; 
    
    % Create structures for storing flys' MADs in each walk dir.:
    for dirr = 1:numel(armDirNames)
        armDirName = armDirNames{dirr};
        MADs.([dataType(10:end-2)]).xMadArmMid.(armDirName) = ...
            nan( numel(selectedFlies), numel(transNumsVect) );
    end
    % Also create structures for #trials for each MAD|trans in each fly and
    % dataset:
    numTrialsConsidered.(dataType(10:end-2)) = nan( ...
        numel(selectedFlies), numel(transNumsVect) ); 
    

    % Run over flies, compute MADs for walking dir and plot MAD-based
    % separation of TPI(y)'s:

    for nF = 1:length(selectedFlies) 
    
        % Read fly's data:
        f = selectedFlies(nF);
        fly_xytyd = meanX_vects_ALL.(['f' num2str(f)]).xytyd;
        fly_locsBotTrialT = meanX_vects_ALL.(['f' num2str(f)]).locsBotTrialT;
        fly_Dec = meanX_vects_ALL.(['f' num2str(f)]).Dec;
        
        % Read fly's botttom trials:
        fly_xytyd_good = fly_xytyd( fly_locsBotTrialT );
        fly_Dec_good = fly_Dec( fly_locsBotTrialT );

        % Create structure for storing x coordinates for 'up' walk dir.:
        dir = 1;
        armDirName = armDirNames{dir};
        flyTrlsXDir.(armDirName) = cell(numel(fly_xytyd_good),1);
        
        flyTrialsConsidered = zeros( numel(fly_xytyd_good), 1); 

        for t = 1:numel(fly_xytyd_good)
    
            % read fly's x,y in trial (downsampling):
            fly_xytyd_good_trial = fly_xytyd_good{t};
            fly_xGood = fly_xytyd_good_trial(1:kSubSamp:end,1);
            fly_yGood = fly_xytyd_good_trial(1:kSubSamp:end,2);
    
            % Only use trials with full data:
            if min(fly_yGood)<=-armTopEnd && max(fly_yGood)>=armTopEnd 
                
                flyTrialsConsidered(t) = 1; 

                % Run over arm walking dir's (up, down, both) and store 
                % corresp. x traj. data (constraints over y):
                for dir = 1:numel(armDirNames)
                    armDirName = armDirNames{dir};
                    armDirLoc = abs( armDirLocs(dir,:) );
                    armDirPol = sign( armDirLocs(dir,:) );

                    % Bottom arm x coordinates within y arm dir. def. 
                    % (3rd constraint relevant for "both" cond.):
                    xInArmYDir = fly_xGood( ...
                        (fly_yGood >= ...
                        (armEdges(armDirLoc(1))*armDirPol(1))) & ...
                        (fly_yGood <= ...
                        (armEdges(armDirLoc(2))*armDirPol(2))) & ...
                        (abs(fly_yGood) >= culDeSacEdge.(dataType)) );

                    % Store x values (of walking dir.) in fly's cell array:
                    flyTrlsXDir.(armDirName){t} = xInArmYDir;
                end
    
            end
    
        end

        
        % Find normalization for x such that 0<=x''<=1:
        flyTrlsXDir_allFlyTrials = cell2mat( flyTrlsXDir.(armDirName) );
        % Find min(x|walk dir.) in all of the fly's trials:
        minX_flyAllTrial = min(flyTrlsXDir_allFlyTrials);
        % Find max( x- min(x|walk dir.) ) in all of the fly's trials:
        max_XMinusMinX_flyAllTrial = max( flyTrlsXDir_allFlyTrials - ...
            minX_flyAllTrial );

        % Load fly's number of transitions in each trial:
        flyNumTrans = numTrans.(dataType(10:end-2)).yAfterCul{nF};

        % Run over #trnasitions and compute the fly's MAD score for 'up' 
        % walking dir, conditioned on #trans=k:
        for tr = 1:numel(transNumsVect)
            transNumK = transNumsVect(tr);
            % Load relevant trials:
            locRelTrials = ( (flyNumTrans==transNumK) & ...
                (flyTrialsConsidered==1) );
            flyTrlsXDir_flyKTrans = cell2mat( flyTrlsXDir.(armDirName)( ...
                locRelTrials ) );
            % Noramlize x values s.t. x'' = [0,1]:
            flyTrlsXDir_Nrmlzd = ( flyTrlsXDir_flyKTrans - ...
                minX_flyAllTrial ) / max_XMinusMinX_flyAllTrial; 
            % Save Median Absolute Deviation around the MIDDLE: 
            MADs.([dataType(10:end-2)]).xMadArmMid.(armDirName)(nF,tr) = ...
                median(abs(flyTrlsXDir_Nrmlzd-.5)); % = MEDIAN( ABS(X-.5) )
            % Also store #trials considered for each #trans=k:
            numTrialsConsidered.(dataType(10:end-2))(nF,tr) = ...
                sum( locRelTrials ); 
        end

        % Plot fly's MAD|#transitions=k:
        figure(figDataset);
        subplot(1,numel(dataTypes_stable),d);
        plot( transNumsVect, ...
            MADs.([dataType(10:end-2)]).xMadArmMid.(armDirName)(nF,:) );
        hold on;


        % Also compute the TPI in each fly without zero transitions:
        locNon0Trials = ( (flyNumTrans~=0) & (flyTrialsConsidered==1) );
        [TPI_w_fly,~, ~, ~, TPI_w_fly_mpInf, ~] = TPI( ...
            fly_xytyd_good(locNon0Trials), yEdges, ...
            fly_Dec_good(locNon0Trials), 'y' );
        % Store fly TPI(y):
        savedTPI.(dataType(10:end-2)).y(nF,:) = TPI_w_fly;
        savedTPIinf.(dataType(10:end-2)).y(nF,:) = ...
            TPI_w_fly_mpInf;

        clear flyTrlsXDir;

    end

    clear meanX_vects_ALL;

    % Plot average MAD|#transitions=k over flies in dataset:
    MAD_mat_dataset = MADs.([dataType(10:end-2)]).xMadArmMid.(armDirName);
    figure(figDataset);
    errorbar( transNumsVect, mean(MAD_mat_dataset,'omitnan'), ...
        std(MAD_mat_dataset,'omitnan') ./ ...
        sqrt( sum(~isnan(MAD_mat_dataset)) ), 'Color', dataCols{d}, ...
        'lineWidth', 2, 'capSize', 0 );

    % Also plot on avg figure;
    figure(figAvgs);
    meanVect = mean(MAD_mat_dataset,'omitnan');
    semVect = std(MAD_mat_dataset,'omitnan') ./ ...
        sqrt( sum(~isnan(MAD_mat_dataset)) );
    fill( [transNumsVect, flip(transNumsVect)], ...
        [meanVect+semVect, flip(meanVect-semVect)], ...
        'k', 'faceColor', dataCols{d}, 'FaceAlpha', .5, ...
        'EdgeColor', 'none' );
    hold on;
    plot( transNumsVect, meanVect, '-', 'Color', dataCols{d}, ...
        'LineWidth', 1 );
    hold on;

end


% Add axes to figures:
figure(figDataset);
for d = 1:length(dataTypes)
    subplot(1,numel(dataTypes_stable),d);
    dataType = dataTypes{d};
    title( dataType(10:end-2) );
    xlabel('k');
    ylabel('MAD | #Cross. = k');
    ylim([0,Inf]);
end
figure(figAvgs);
xlabel('k');
ylabel('MAD | #Crossings = k');   
xticks(0:transNumsVect(end));

% Add significance for avg. figure - Wilcoxon Rank-sum test (corrected for 
% multiple comparisons):
pValMAt = nan( numel(dataTypes), numel(transNumsVect) );
MAD_mat_WT = MADs.Short.xMadArmMid.(armDirName);
alphaCorrected = 0.05 / ( numel(transNumsVect) * (numel(dataTypes)-1) );
for d = length(dataTypes):-1:2
    dataType = dataTypes{d};
    MAD_mat_mutant = MADs.([dataType(10:end-2)]).xMadArmMid.(armDirName);
    for tr = 1:numel(transNumsVect)
        pValMAt(d,tr) = ranksum( MAD_mat_WT(:,tr), MAD_mat_mutant(:,tr) );
        if pValMAt(d,tr) <= alphaCorrected
            meanOfSig = mean( MAD_mat_mutant(:,tr), 'omitnan' );
            plot( transNumsVect(tr), meanOfSig, 'k*', 'LineWidth', 2 );
            hold on;
            plot( transNumsVect(tr), meanOfSig, '*', 'Color', ...
                dataCols{d}, 'LineWidth', 1 );
            hold on;
        end
    end
end
axis tight;  
xlim([-.5,.5+transNumsVect(end)]);

% Save avg. figure:
xtickangle(0);
figStartName = ['analyses/figures/flies_MAD_byTranNum_shortComp_' armDirName ...
    'ArmWalk'];
saveas( figAvgs, [figStartName '.fig']);

% Save output TPI's:
save( 'datafiles/flies_TPIs_noZero.mat', 'savedTPI', 'binEdgesTPI', ...
    'savedTPIinf' );



%% NON-ZERO TRANSITIONS: Plots AVG.+-SEM TPI(y) & cul-de-sac increase est's:

% Compare shot-maze fly datasets' AVG.+-SEM TPI(y/time) AND cul-de-sac 
% increase est's for trials with non-zero post-cul-de-sac transitions.

% Load TPI(y) of WT and mutants in short maze EXCLUDING ZERO TRANSITIONS:
load( 'flies_TPIs_noZero.mat', 'savedTPI', 'binEdgesTPI', 'savedTPIinf' );

datatypes = {'Turndata_ShortYy', 'Turndata_LongYy', ...
    'Turndata_ShortNorpAYy', 'Turndata_ShortDumbYy', ...
    'Turndata_ShortFoxPYy', 'Turndata_ShortNompCYy'};
dataCols = {'k', 'b', 'r', 'm', 'g', '#77AC30'}; 


% Define dataset comparisons:
compNames = {'Short'};
compLoc.Short = [1, 5, 4, 3, 6]; % short fly comparison
infVal_global = 1.5;

% Define domain (spatial):
binTypes = {'y'};
binPrintNames = {'y'};

% Define bottom arm cul-de-sac and maze ROIs (similar for all short-maze 
% fly datyasets):
culDeSacY = 0.69/2.03;
% Define maze ROIs:
regionsBreak_Strct.y = [-1.3, -1, -culDeSacY, ...
    culDeSacY, 1, 1.3];
regionPrintNames_Strct.y = {'Pre-arm', 'Down-the-arm', ...
    'Cul-de-sac', 'Up-the-arm', 'Post-arm'};
regionsBreak_Strct.time = [-1.3, -1, 0, 1, 1.3];
regionPrintNames_Strct.time = {'Pre-arm', ...
    'Down-the-arm & c.d.s<0', 'C.d.s>0 & up-the-arm', ...
    'Post-arm' };
regionBreakTicks.y = -1.:.5:1;
regionBreakTicks.time = -1.:.5:1;
infVal = 1.5;


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

        % Also create TPI(y) slope-comp. fig.:
        if strcmp(binType,'y')
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

            if strcmp(datatype,'Turndata_ShortYy')
                legNames{d} = 'WT';
            elseif strcmp(datatype,'Turndata_LongYy')
                legNames{d} = 'WT long';
            elseif strcmp(datatype(end-1:end),'Yy')
                legNames{d} = datatype(15:end-2);
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
    
            % Define edges/centers:
            wBinEdges = binEdgesTPI.(datatype(10:end-2)).(binType);
            wBinCenters = .5*diff(wBinEdges(1:2)) + wBinEdges(1:end-1);
            errorbar( wBinCenters, ...
                mean( savedTPI.(datatype(10:end-2)).(binType), ...
                'omitnan' ), ...
                std( savedTPI.(datatype(10:end-2)).(binType), ...
                'omitnan' ) ./ ...
                sqrt( sum( ~isnan(savedTPI.(datatype(10:end-2)).( ...
                binType)) ) ), ...
                'Color', 'w', 'capSize', 0, 'lineWidth', 4 );
            hold on;
            legs(d) = errorbar( wBinCenters, ...
                mean( savedTPI.(datatype(10:end-2)).(binType), ...
                'omitnan' ), ...
                std( savedTPI.(datatype(10:end-2)).(binType), ...
                'omitnan' ) ./ ...
                sqrt( sum( ~isnan(savedTPI.(datatype(10:end-2)).( ...
                binType)) ) ), ...
                'Color', dataCol, 'capSize', 0, 'lineWidth', 2 );
            hold on;
            % Plot avg. TPI(+-Inf):
            errorbar(  [-infVal,infVal_global], ...
                mean( savedTPIinf.(datatype(10:end-2)).(binType), ...
                'omitnan' ), ...
                std( savedTPIinf.(datatype(10:end-2)).(binType), ...
                'omitnan' ) ./ ...
                sqrt( sum( ~isnan(savedTPIinf.(datatype(10:end-2)).( ...
                binType)) ) ), ...
                'Color', dataCol, 'capSize', 0, 'lineWidth', 2, ...
                'LineStyle', 'none', 'Marker', 'o' );


            
            % TPI SLOPE FIGURE:

            % Also compute slopes in cul-de-sac:

            if strcmp(binType,'y')

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
                if d==1 && strcmp(compName, 'Short')
                    TPI_lastInCul_WT = TPI_lastInCul;
                    TPI_firstInCul_WT = TPI_firstInCul;
                    TPI_minusInf_WT = TPI_minusInf; 
                elseif d>1 && strcmp(compName, 'Short')
                    pVal1 = ranksum(TPI_lastInCul,TPI_lastInCul_WT);
                    pVal2 = ranksum(TPI_lastInCul - abs( TPI_firstInCul ...
                        ), TPI_lastInCul_WT - abs( TPI_firstInCul_WT ));
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
            ax2 = axes('Position', get(ax1,'Position'), ...
                'XAxisLocation', 'top', 'Color', 'none', 'XColor', 'k' );
            ax2.YAxis.Visible = 'off';
            ax2.XLim = ax1.XLim;
            ax2.XTick = regionsBreak(1:end-1) + .5*diff(regionsBreak);
            ax2.XTickLabel = regionPrintNames;
            ax2.FontSize = 9;
        end
        % Save figure: 
        figStartName = ['analyses/figures/fliesAll_' compName ...
            'Comp_TPI' upper(binType(1)) '_noZeroTrans'];
        saveas( figTPI.all.(binType), [figStartName '.fig']);


        % TPI SLOPE figure axes: 
        if strcmp(binType,'y')
            figure(figSlopeTPI.all.(binType));
            for sp = 1:3
                subplot(3,1,sp);
                plot([0,numel(datatypes_c)+1],[0,0],'k:');
                hold on;
                % Add significance of Wilcoxon rank sum test:
                if strcmp(compName,'Short')
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
                'Comp_slopeTPI' upper(binType(1)) '_noZeroTrans'];
            saveas( figSlopeTPI.all.(binType), [figStartName '.fig']);
        end
        
    end
    
end


