%% QUANTIFYING LOCAL LATERAL TENDENCIES VIA MIDLINE_CROSSINGS & PEVEN: 
%
% Computes trial #midline-crossing in each maze region, typicaly post 
% cul-de-sac, and the corresponding fractions of even #crossings, pEven. 
% Plots #Crossings and pEven histogram and boxplots, comparison with 
% Poisson predictions, pEven-based TPIs, and probability density kernal 
% estimates of time spent in each maze region.
%
% INPUT FILES: 
%   output_drosophila_main_[DATASET_NAME(10:END)].mat
%   flies_TPIs.mat
%
% FUNCTION:
%   TPI.m for computing the TPI.
%
% EXTERNAL FUNCTION:
%   This code uses an external function, myBinomTest.m, Reference:
%   Matthew Nelson (2015). 
%   https://www.mathworks.com/matlabcentral/fileexchange/24813-mybinomtest-s-n-p-sided
%   MATLAB Central File Exchange. Retrieved February 9, 2016.
%
% OUTPUT FILE: 
%   flies_pEvenTrans.mat
%
% FIGURES:
%   Within dataset region-dep. pEven and #Crossings distributions:
%       Figs. 2C, 3C, 3D: flies_pEvenTrans_hists_[DATASET].fig
%   Within dataset pEven-based TPI(y) (Avg. +- SEM) comparison:
%       Figs. 2D, S3D, S11B: flies_pEvenTrans_extrmTPI_[DATASET].fig
%   Between dataset pEven histogram comparison:
%       Figs. S6D (S9A), S13D: flies_pEvenTrans_comp[COMPNAME]PEven.fig
%   Between dataset pEven boxplot comparison:
%       Figs. 4C (S9B), S13C: flies_pEvenTrans_comp[COMPNAME]PEven2.fig
%   Between dataset #crossings histogram comparison:
%       Figs. S6C, S13E: 
%       flies_pEvenTrans_comp[COMPNAME]Trans_[REGIONNAME].fig
%   Datasets' observed avg. #Crossings and pEven vs pEven expected under  
%       a Possion distribution:
%       Fig. S15C: flies_pEvenTransPoisson.fig
%   Dist. of avg. #Crossings within bottom arm for short-maze WT flies:
%       Fig. S1B: flies_histAvgTransBottomArm_Short.fig
%   Time in maze-region probability density estimate:
%       Figs. S3A, S6B: flies_pEvenTrans_TimeInRegion.fig
%
% Note in the current midline-crossings are termed "x transitions".
% Note that the manuscript primarily focuses on #Crossings and pEven for
% (upward) motion away from the cul-de-sac.
% More information about #crossings (#transitions) and pEven and corresp. 
% analyses - in the sections below.
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



%% Define datasets, transition type, maze-regions and turn-directions:

% Define datasets:
dataTypes = {'Turndata_ShortYy', 'Turndata_LongYy', ...
    'Turndata_ShortNorpAYy', 'Turndata_ShortDumbYy', ...
    'Turndata_ShortFoxPYy','Turndata_ShortNompCYy', ...
    'Turndata_ShortABM_NoWF', ...
    'Turndata_ShortABM_WF_00100', 'Turndata_ShortABM_WF_00200', ...
    'Turndata_ShortABM_WF_00500', 'Turndata_ShortABM_WF_00700', ...
    'Turndata_ShortABM_WF_00900', ...
    'Turndata_ShortABM_WallFollowing', 'Turndata_ShortABM_brown5'};
dataCols = {'k', 'b', 'r', 'm', 'g', '#77AC30', ...
    "#A2142F", ... % no WF 
    "#D95319", "#EDB120", "#77AC30", "#4DBEEE", "#0072BD", ... % WF parm.
    "#7E2F8E", "#7E2F8E" }; % WF, brownian

% Define transition types:
transEvenTypes = {'even', 'odd', 'zero'};
transEvenEvenDivs = [0,1,999];

% Define maze-regions (yRanges) for which transitions will be computed:
yRangesNames = {'yBeforeCul', 'yDownArm', 'culDeSac', 'yUpArm', ...
    'yAfterCul'};
yRangesPrintNames = {'y<-Cul', '-1<y<-Cul', '-Cul<y<Cul', 'Cul<y<1', ...
    'y>Cul'};
yRangesCols = {'g', 'r', 'k', 'b', 'c'};
yRangesShapes = {'h', 's', 'd', '*', 'o'};

% Define turn-directions for which transitions will be computed [change 
% decVals to 1 or -1 to compute for right or left turns, resp.]:
decVals = 10; % BOTH Left and Right turns
decCols = {'k'};
decCols2 = {'k'};
turnNames = {'BOTH'};

infVal_global= 1.5;



%% Computes & plots #transitions, pEven by region x motion-dir, pEven-based TPI:

% Computes #transitions per type and maze-region for each dataset.
% Computes pEven (fraction of even #transitions) per fly, region,
    % motion-dir. 
% Plots region x walking dir. pEven (including/excluding zero-transition 
    % trials) and #transitions histograms for each dataset.
% Plots within dataset pEven-based TPI comparison.
% Also computes test statistics (e.g., sig. of individual pEven, Wilcoxon 
    % ranksum test for the med. of pEven-0.5 in each dataset).
% Also computes parameters relevant for long-maze-based simulations (used 
    % in drosophila_longMazeTransAndSims.m).


% Load TPI data:
load('flies_TPIs.mat');

for d = 1:numel(dataTypes) 
    dataType = dataTypes{d};
    dataCol = dataCols{d};
    if strcmp(dataType, 'Turndata_ShortABM_brown5')
        dataCol = [.5, .3, .1];
    end

    % Load dataset:
    matFileName = ['output_drosophila_main_' dataType(10:end) '.mat'];
    load(matFileName);


    % Create cell storing inter-transition-distances for all short WT flies:
    if strcmp(dataType, 'Turndata_LongYy')
        deltaY_IntTransInt_flies = cell(length(selectedFlies),1);
        deltaY_armOnly_IntTransInt_flies = cell(length(selectedFlies),1);
        deltaY_IntTransInt_flies_cc = cell(length(selectedFlies),1);
        deltaY_armOnly_IntTransInt_flies_cc = cell(length(selectedFlies),1);
        wDeltaY_IntTransInt_flies = cell(length(selectedFlies),1);
        wDeltaY_armOnly_IntTransInt_flies = cell(length(selectedFlies),1);
        yVals_trans_flies = cell(length(selectedFlies),1);
        % Also define upper edge (end of intersection)
        if strcmp(dataType, 'Turndata_LongYy')
            interUpperY = 1.158;
        elseif strcmp(dataType, 'Turndata_ShortYy')
            interUpperY = 1.27;
        end
    end

    % Define cul-de-sac y and y regions break (for TPIs of extreme pEven 
    % groups):
    if strcmp(dataType,'Turndata_LongYy')
        culdesacY = 0.68/3.36;
        infVal = infVal_global;
        regionsBreak = [-1.3, -1, -culdesacY, culdesacY, 1, 1.3];
        regionPrintNames = {'Pre-arm', 'Down-the-arm', 'Cul-de-sac', ...
            'Up-the-arm', 'Post-arm'};
        regionsTicks = -1.:.5:1;
    elseif strcmp(dataType,'Turndata_ShortYy') || ...
            strcmp(dataType,'Turndata_ShortNorpAYy') || ...
            strcmp(dataType,'Turndata_ShortDumbYy') || ...
            strcmp(dataType,'Turndata_ShortFoxPYy') || ...
            strcmp(dataType,'Turndata_ShortNompCYy')
        culdesacY = 0.69/2.03;
        infVal = infVal_global;
        regionsBreak = [-1.3, -1, -culdesacY, culdesacY, 1, 1.3];
        regionPrintNames = {'Pre-arm', 'Down-the-arm', 'Cul-de-sac', ...
            'Up-the-arm', 'Post-arm'};
        regionsTicks = -1.:.5:1;
    elseif strcmp(dataType,'Turndata_ShortABM_NoWF') || ...
            strcmp(dataType,'Turndata_ShortABM_WallFollowing') || ...
            strcmp(dataType,'Turndata_ShortABM_WF_00100') || ...
            strcmp(dataType,'Turndata_ShortABM_WF_00200') || ...
            strcmp(dataType,'Turndata_ShortABM_WF_00500') || ...
            strcmp(dataType,'Turndata_ShortABM_WF_00700') || ...
            strcmp(dataType,'Turndata_ShortABM_WF_00900') || ...
            strcmp(dataType,'Turndata_ShortABM_brown5')
        culdesacY = 7.5/30;
        infVal = culdesacY+0.2;
        regionsBreak = [-culdesacY, culdesacY, 1, 1.3];
        regionPrintNames = {'Cul-de-sac', 'Up-the-arm', 'Post-arm'};
        regionsTicks = 0:.5:1;
    end


    % Create figure for [...] COMPLETE INFO.!
    histsFig = figure( 'Units', 'Normalized', 'Position', [0 .3 1 .7] );

    % Define y-ranges for which #transitions will be computed:
    yRanges4Trans = [-Inf, -culdesacY; ...
        -1, -culdesacY; ...
        -culdesacY, culdesacY; ...
        culdesacY, 1; ...
        culdesacY, Inf];

    % Run over y-ranges and compute #transitions for each:
    
    for yr = 1:size(yRanges4Trans,1)
        
        yRanges4Tran = yRanges4Trans(yr,:);
        yRangeL = yRanges4Tran(1);
        yRangeH = yRanges4Tran(2);
        yRangesName = yRangesNames{yr};
        yRangesCol = yRangesCols{yr};
        yRangesShape = yRangesShapes{yr};

        % Create cells for storing #transitions:
        numTrans.(dataType(10:end-2)).(yRangesName) = ...
            cell(size(selectedFlies));

        % Create vectors for storing avg.(#transitions):
        avgNumTrans.(dataType(10:end-2)).(yRangesName) = ...
            nan(numel(selectedFlies),1);
    
        % Create matrices for storing fly #odd/even/zero trans. by trial type:
        numTransType.(dataType(10:end-2)).(yRangesName) = ...
            nan(numel(selectedFlies), numel(transEvenTypes));


        for nF = 1:length(selectedFlies)
            f = selectedFlies(nF);
            flyName = ['f' num2str(f)];
    
            % Read individual animal's processed data for bottom trials: 
            
            % Location of bottom trials:
            locsBotTrialT = meanX_vects_ALL.(flyName).locsBotTrialT;
            % Corrected trajectories in bottom trials: 
            xytyd = meanX_vects_ALL.(flyName).xytyd;
            xytyd_bottom = xytyd( locsBotTrialT );
            % Decision vector for bottom trials:
            Dec = meanX_vects_ALL.(flyName).Dec;
    
            % Find fly's med x loc.:
            traj_cell_mat = cell2mat( xytyd_bottom );
            xy_bottom_mat = traj_cell_mat(:,1:2);
            uniqueXY = unique( xy_bottom_mat, 'rows' );
            xGoodMedian_allTrajs = mean( ...
                [min( uniqueXY( abs(uniqueXY(:,2)) < 1, 1 ) ), ...
                max( uniqueXY( abs(uniqueXY(:,2)) < 1, 1 ) )] );
    
            % Create cell storing inter-transition-distances for current 
            % short WT fly:
            if strcmp(yRangesName, 'yAfterCul') && ...
                    strcmp(dataType, 'Turndata_LongYy')
                deltaY_IntTransInt_fly = cell(size(xytyd_bottom,1),1);
                wDeltaY_IntTransInt_fly = cell(size(xytyd_bottom,1),1);
                yVals_trans_fly = cell(size(xytyd_bottom,1),1);
            elseif strcmp(yRangesName, 'yUpArm') && ...
                    strcmp(dataType, 'Turndata_LongYy')
                deltaY_armOnly_IntTransInt_fly = cell(size(xytyd_bottom,1),1);
                wDeltaY_armOnly_IntTransInt_fly = cell(size(xytyd_bottom,1),1);
            end

            for s = 1:length(decVals)
                % Read turn-dir. value and colors:
                decVal = decVals(s);
                decCol = decCols{s};
                decCol2 = decCols2{s};
                turnName = turnNames{s};
        
                % Use only fly's bottom trials that ended in specific turn-dir.:
                if decVal==10
                    xytydBotDec = xytyd( locsBotTrialT );   
                else
                    xytydBotDec = xytyd( locsBotTrialT & Dec == decVal);   
                end
                
                % TRANSITIONS:
                
                % Compute number of transitions in each trial WITHIN THE 
                % CULDESAC:
                numTransitions = nan( numel(xytydBotDec), 1 );
                isRepeatDec = nan( numel(xytydBotDec), 1 );
                for t = 1:numel(xytydBotDec)
                    xy_raw = xytydBotDec{t}( :, 1:2 );
                    xy_InReg = xy_raw( (xy_raw(:,2) < yRangeH) & ...
                        (xy_raw(:,2) > yRangeL) , : ); 
                    xInReg = xy_InReg( : , 1 ); 
                    yInReg = xy_InReg( 2:end , 2 ); 
                    xInRegCentered = xInReg - xGoodMedian_allTrajs;
                    xRightOrLeft = ( xInRegCentered >= 0 ) - ...
                        ( xInRegCentered < 0 );
                    isTransitionM1 = xRightOrLeft(1:end-1) .* ...
                        xRightOrLeft(2:end);
                    % Store fly's #transitions in current trial:
                    numTransitions(t) = sum(isTransitionM1 == -1);
                    % Compute inter-transition-distances (y dist., from cul 
                    % to intersection), only for WT short:
                    if strcmp(yRangesName, 'yAfterCul') && ...
                            strcmp(dataType, 'Turndata_LongYy')
                        yOfTransitions = yInReg(isTransitionM1 == -1);
                        deltaY_IntTransInt = diff( ...
                            [culdesacY; yOfTransitions; interUpperY] );
                        deltaY_IntTransInt_fly{t} = deltaY_IntTransInt;
                        wDeltaY_IntTransInt_fly{t} = ...
                            ones( size(deltaY_IntTransInt) ) / ...
                            length(deltaY_IntTransInt);
                        yVals_trans_fly{t} = yOfTransitions;
                    elseif strcmp(yRangesName, 'yUpArm') && ...
                            strcmp(dataType, 'Turndata_LongYy')
                        yOfTransitions = yInReg(isTransitionM1 == -1);
                        deltaY_IntTransInt = diff( ...
                            [culdesacY; yOfTransitions; 1] );
                        deltaY_armOnly_IntTransInt_fly{t} = deltaY_IntTransInt;
                        wDeltaY_armOnly_IntTransInt_fly{t} = ...
                            ones( size(deltaY_IntTransInt) ) / ...
                            length(deltaY_IntTransInt);
                    end
                end
    
                % Store fly's #transitions by trialType:
                numTrans.(dataType(10:end-2)).(yRangesName){nF} = ...
                    numTransitions;

                % Store fly's avg.(#transitions) by trialType:
                avgNumTrans.(dataType(10:end-2)).(yRangesName)(nF) = ...
                    mean(numTransitions);
    
                % Compute and store fly's fly #odd/even/zero trans. AND 
                % pRepeat|#odd/even/zero trans. by trial type:
                for tra = 1:length(transEvenTypes)
                    transEvenType = transEvenTypes{tra};
                    transEvenEvenDiv = transEvenEvenDivs(tra);
                    if strcmp( transEvenType, 'zero' )
                        numEvType = sum( numTransitions == 0);
                    else 
                        numTransitionsWith0 = numTransitions;
                        numTransitionsGood = numTransitions( ...
                            numTransitions ~= 0 );
                        numEvType = sum( rem(numTransitionsGood,2) == ...
                            transEvenEvenDiv );
                    end
                    numTransType.(dataType(10:end-2)).(yRangesName)(nF,tra) ...
                        = numEvType;
                end
    
            end

            % Store long dataset arrays (for 
                % drosophila_longMazeTransAndSims.m): 
            if strcmp(yRangesName, 'yAfterCul') && ...
                    strcmp(dataType, 'Turndata_LongYy')
                deltaY_IntTransInt_flies{nF} = cell2mat( ...
                    deltaY_IntTransInt_fly );
                deltaY_IntTransInt_flies_cc{nF} = deltaY_IntTransInt_fly;
                wDeltaY_IntTransInt_flies{nF} = cell2mat( ...
                    wDeltaY_IntTransInt_fly );
                yVals_trans_flies{nF} = cell2mat( ...
                    yVals_trans_fly );
            elseif strcmp(yRangesName, 'yUpArm') && ...
                    strcmp(dataType, 'Turndata_LongYy')
                deltaY_armOnly_IntTransInt_flies{nF} = cell2mat( ...
                    deltaY_armOnly_IntTransInt_fly );
                deltaY_armOnly_IntTransInt_flies_cc{nF} = ...
                    deltaY_armOnly_IntTransInt_fly;
                wDeltaY_armOnly_IntTransInt_flies{nF} = cell2mat( ...
                    wDeltaY_armOnly_IntTransInt_fly );
            end
    
        end
        

        % Plot histograms of:
            % Columns: by y-Range name (before/within/after cul-de-sac).
            % Row 1/2: counts of flies' pEven including/excluding zero 
                % Transitions (numBins = #flies).
            % Row 3: Prob. #Transition over all trials of all flies.

        figure(histsFig);
        
        nTrans_rFly_cEOZ = numTransType.(dataType(10:end-2)).(yRangesName);
        nTransFly_even = nTrans_rFly_cEOZ(:,1);
        nTransFly_odd = nTrans_rFly_cEOZ(:,2);
        nTransFly_zero = nTrans_rFly_cEOZ(:,3);
        nTransFly_evenAndZero = nTransFly_even + nTransFly_zero;
        pEven_with0 = nTransFly_evenAndZero ./ ...
            (nTransFly_evenAndZero + nTransFly_odd );
        pEven_without0 = nTransFly_even ./ (nTransFly_even + nTransFly_odd );
        % Store pEven w/o for later datasets comparison:
        pEven_with0_strct.(yRangesName).(dataType(10:end-2)) = pEven_with0;
        pEven_without0_strct.(yRangesName).(dataType(10:end-2)) = ...
            pEven_without0;
        
        % Plot pEven counts (INcluding zero transitions):
        subplot(3,size(yRanges4Trans,1),yr);
        histogram( pEven_with0, 'binEdges', -.025:.05:1.025, ...
            'FaceColor', dataCol );
        xlabel('pEven'); xticks(0:.1:1);
        ylabel('#Flies');
        pVal = signtest(pEven_with0-.5);
        title([dataType(10:end-2) ', ' yRangesName ...
            ': with0 signTestP=' num2str(pVal)]);
        
        % Plot pEven counts (Excluding zero transitions):
        subplot(3,size(yRanges4Trans,1),yr+size(yRanges4Trans,1));
        histogram( pEven_without0, 'binEdges', -.025:.05:1.025, ...
            'FaceColor', dataCol );
        xlabel('pEven'); 
        ylabel('#Flies');
        xticks(0:.1:1);
        pVal = signtest(pEven_without0-.5);
        title(['without0 signTestP=' num2str(pVal)]);
        
        % Plot Prob.(#Transitions=k) over ALL trials of all flies:
        numTransAllTrialsFlies = cell2mat( ...
            numTrans.(dataType(10:end-2)).(yRangesName) );
        subplot(3,size(yRanges4Trans,1),yr+2*size(yRanges4Trans,1)); 
        if strcmp(dataType,'Turndata_ShortABM_brown5')
            h = histogram( numTransAllTrialsFlies, 'normalization', ...
                'probability', 'FaceColor', dataCol ); 
        else
            h = histogram( numTransAllTrialsFlies, 'normalization', ...
                'probability', 'BinEdges', [-.5:29.5, Inf], 'FaceColor', ...
                dataCol ); 
        end
        hold on;
        ggg = gca;
        maxYLim = ggg.YLim(2);
        if ~strcmp(dataType,'Turndata_ShortABM_brown5')
            plot( repmat(-.5:29.5,[2,1]), repmat([0;maxYLim],[1,31]), 'k');
        end
        xlabel('#Crossings'); 
        ylabel('Prob.');
        numTransAllTrialsFlies_without0 = numTransAllTrialsFlies(...
            numTransAllTrialsFlies~=0 );
        title(['all trials,flies: %Even=' ...
            num2str( mean(rem(numTransAllTrialsFlies,2)==0) )  ' [' ...
            num2str( mean(rem(numTransAllTrialsFlies_without0,2)==0) ) ...
            ']' ]);
        % Store:
        if ~strcmp(dataType, 'Turndata_LongYy')
            numTransProb.(dataType(10:end-2)).(yRangesName) = h.Values;
        end
        clear h;
    

        % Compute statistics of #transitions:
            
        numTransCell = numTrans.([dataType(10:end-2)]).(yRangesName);

        % Compute statistics for #transitions [trials of all flies]:
        allFliesTransStats.([dataType(10:end-2)]).(yRangesName).N = ...
            numel( cell2mat( numTransCell ) );
        allFliesTransStats.([dataType(10:end-2)]).(yRangesName).N_even = ...
            sum( rem( cell2mat( numTransCell ), 2 ) == 0 );
        allFliesTransStats.([dataType(10:end-2)]).(yRangesName).N_zero = ...
            sum( cell2mat( numTransCell ) == 0 );

        % Compute statistics for #transitions [avg. fly]:
        avgTrans = nan(numel(numTransCell),1);
        for f = 1:numel(numTransCell)
            numTrans_fly = numTransCell{f};
            avgTrans(f) = mean( numTrans_fly );
        end
        avgFlyTransStats.([dataType(10:end-2)]).(yRangesName).mean = ...
            mean(avgTrans);
        avgFlyTransStats.([dataType(10:end-2)]).(yRangesName).std = ...
            std(avgTrans);
        avgFlyTransStats.([dataType(10:end-2)]).(yRangesName).median = ...
            median(avgTrans);
        avgFlyTransStats.([dataType(10:end-2)]).(yRangesName).IQR = ...
            quantile(avgTrans,[.25,.75]);


        % Sig. test of pEven and pEven-based TPI comparison of 2 extreme 
        % groups --> ONLY FOR Up-the-arm MOTION: 

        if strcmp( yRangesName, 'yAfterCul' )


            % Test significace of p(Even) > 0.5 (for each fly):

            pVal_pEven_with0_half = nan(numel(selectedFlies),1);
            for f = 1:numel(selectedFlies)
                nTrials_fly = numel(numTransCell{f});
                pVal_pEven_with0_half(f) = myBinomTest( nTrials_fly * ...
                    pEven_with0(f) , nTrials_fly, 0.5, 'two' );
            end
            isSig_pEven_with0_half = (pVal_pEven_with0_half <= 0.05);
            pEvenStats.([dataType(10:end-2)]).(yRangesName).nSig = ...
                sum( isSig_pEven_with0_half );
            pEvenStats.([dataType(10:end-2)]).(yRangesName).perSig = ...
                mean( isSig_pEven_with0_half );
            pEvenStats.([dataType(10:end-2)]).(yRangesName).maxPVal = ...
                max( pVal_pEven_with0_half(isSig_pEven_with0_half) );


            % Compare 2 extreme groups in their TPIs:
            
            figExtrmEvenTPIs = figure;
            % Add regions break:
            ax1 = axes();
            box(ax1);
            for rr = 1:numel(regionsBreak)
                plot(regionsBreak(rr) * [1,1], [-1,1], 'k:' ); 
                hold on;
            end

            % Compute avg. TPI(y) for each extreme pEven group (20% flies each):
            n20Per = round( 0.2 * numel(selectedFlies) );
            [B,I] = sort(pEven_with0);
            TPI_y = savedTPI.(dataType(10:end-2)).y;
            TPI_y_edges = binEdgesTPI.(dataType(10:end-2)).y;
            TPI_y_centers = TPI_y_edges(1:end-1) + ...
                .5 * diff(TPI_y_edges(1:2));
            TPI_y_pmInf = savedTPIinf.(dataType(10:end-2)).y;
            TPI_y_pmInf_centers = [-infVal,infVal_global];
            TPI_y_smallEv = TPI_y(I(1:n20Per),:);
            TPI_y_bigEv = TPI_y(I(end-n20Per+1:end),:);
            TPI_y_smallEv_pmInf = TPI_y_pmInf(I(1:n20Per),:);
            TPI_y_bigEv_pmInf = TPI_y_pmInf(I(end-n20Per+1:end),:);
            
            % Plot TPI's:
            figure(figExtrmEvenTPIs);
            errorbar( TPI_y_centers, mean(TPI_y_smallEv,'omitnan'), ...
                std(TPI_y_smallEv,'omitnan') ./ ...
                sqrt( sum(~isnan(TPI_y_smallEv)) ), 'r', 'lineWidth', 1 );
            hold on;
            errorbar( TPI_y_pmInf_centers, ...
                mean(TPI_y_smallEv_pmInf,'omitnan'), ...
                std(TPI_y_smallEv_pmInf,'omitnan') ./ ...
                sqrt( sum(~isnan(TPI_y_smallEv_pmInf)) ), 'r', ...
                'lineWidth', 1, 'Marker', 'o', 'lineStyle', 'none' );
            hold on;
            errorbar( TPI_y_centers, mean(TPI_y_bigEv,'omitnan'), ...
                std(TPI_y_bigEv,'omitnan') ./ ...
                sqrt( sum(~isnan(TPI_y_bigEv)) ), 'k', 'lineWidth', 1 );
            hold on;
            errorbar( TPI_y_pmInf_centers, ...
                mean(TPI_y_bigEv_pmInf,'omitnan'), ...
                std(TPI_y_bigEv_pmInf,'omitnan') ./ ...
                sqrt( sum(~isnan(TPI_y_bigEv_pmInf)) ), 'k', ...
                'lineWidth', 1, 'Marker', 'o', 'lineStyle', 'none' );
            hold on;

            ylim([-1,1]);
            xticks([-infVal, regionsTicks, infVal_global]);
            midLabels = arrayfun(@num2str, regionsTicks, ...
                'UniformOutput', false);
            firstLabel = '-y_{\infty}';
            lastLabel = 'y_{\infty}';
            xticklabels([firstLabel, midLabels,lastLabel]);
            xlim( [-infVal, infVal_global] ); 
            yticks(-1:.2:1);
            ggg = gca;
            ggg.XGrid = 'off';
            ggg.YGrid = 'on'; 
            xlabel('y'); 
            ylabel('TPI(y)'); 
           
            % Add region names: 
            ax2 = axes('Position', get(ax1,'Position'), 'XAxisLocation', ...
                'top', 'Color', 'none', 'XColor', 'k' );
            ax2.YAxis.Visible = 'off';
            ax2.XLim = ax1.XLim;
            ax2.XTick = regionsBreak(1:end-1) + .5*diff(regionsBreak);
            ax2.XTickLabel = regionPrintNames;
            ax2.FontSize = 9;

            % Inset: add #transitions histograms for each extreme group 
            % (over all trials of all flies in group):
            figure(figExtrmEvenTPIs);
            numTrans_bigEv = cell2mat(numTransCell( I(end-n20Per+1:end) ));
            numTrans_smallEv = cell2mat( numTransCell( I(1:n20Per) ) );
            maxTransExts = max([numTrans_bigEv; numTrans_smallEv]);
            ax3 = axes('Position',[.485 .338 .4 .1]);
            if strcmp(dataType,'Turndata_ShortABM_brown5')
                histogram( numTrans_bigEv, 'FaceColor', 'k', ...
                    'Normalization', 'probability');
            else
                histogram( numTrans_bigEv, 'binEdges', [-.5:29.5, Inf], ...
                    'FaceColor', 'k', 'Normalization', 'probability');
            end
            hold on;
            maxLim3 = ax3.YLim(2);
            if ~strcmp(dataType,'Turndata_ShortABM_brown5')
                plot( repmat(-.5:29.5,[2,1]), repmat([0;1],[1,31]), 'k');
                hold on;
            end
            ylabel('Prob.'); 
            xlabel(' ');
            xticks(0:5:30);
            ax3.FontSize = 9;
            ax4 = axes('Position',[.485 .2 .4 .1]);
            if ~strcmp(dataType,'Turndata_ShortABM_brown5')
                histogram( numTrans_smallEv, 'binEdges', [-.5:29.5, Inf], ...
                    'FaceColor', 'r', 'Normalization', 'probability' );
            else
                histogram( numTrans_smallEv, 'FaceColor', 'r', ...
                    'Normalization', 'probability' );
            end
            hold on;
            maxLim4 = ax4.YLim(2);
            ylabel('Prob.'); 
            xlabel('#Crossings');
            xticks(0:5:30);
            ax4.FontSize = 9;
            maxLim = max([maxLim3,maxLim4]);
            if ~strcmp(dataType,'Turndata_ShortABM_brown5')
                plot( repmat(-.5:29.5,[2,1]), repmat([0;1],[1,31]), 'k');
            end
            ylim([0, maxLim]);
            ax3.YLim = [0, maxLim];
            

        end

    end

    clear meanX_vects_ALL;

    % Save figures:
    figure(histsFig);
    figStartName = ['analyses/figures/flies_pEvenTrans_hists_' ...
        dataType(10:end-2)];
    saveas( histsFig, [figStartName '.fig']);
    figure(figExtrmEvenTPIs);
    figStartName = ['analyses/figures/flies_pEvenTrans_extrmTPI_' ...
        dataType(10:end-2)];
    saveas( figExtrmEvenTPIs, [figStartName '.fig']);


end


% Save data:
close all;
save( 'datafiles/flies_pEvenTrans.mat' );



%% Plots datasets comparison: pEven (post cul-de-sac) hist's and boxplot:

% Plots pEven distribution and boxplots across flies in each dataset.
% pEven is computed for upward motion away from the cul-de-sac.
% To produce these dist's and boxplots for non-zero #transitions (Fig. 
    % S9A-B) change the line below to: "TransZeroName = '_noZeroTrans';".

TransZeroName = '';

% Read relevant pEven structure (including/excluding zero-transition
% trials):
pEven_strct = pEven_with0_strct;
if strcmp( TransZeroName, '_noZeroTrans')
    pEven_strct = pEven_without0_strct;
else
    pEven_strct = pEven_with0_strct;
end

compNames = {'Short','Model'};
compVects = {[1,5,4,3,6], [14, 7:12, 13]};

modelPrintNames = {'Brownian', 'No WF', 'WF 001', 'WF 002', 'WF 005', ...
    'WF 007', 'WF 009', 'WF'};

for c = 1:numel(compNames)
    compName = compNames{c};
    compVect = compVects{c};

    figPEvenComp = figure;
    figPEvenComp2 = figure;
    dataTypes2 = dataTypes(compVect);
    dataCols2 = dataCols(compVect);
    data2TypePrintNames = cell(size(dataTypes2));
    for d = 1:numel(dataTypes2)
        dataType2 = dataTypes2{d};
        dataType3 = dataType2(10:end-2);
        dataCol2 = dataCols2{d};
        if strcmp(dataType2,'Turndata_ShortYy')
            data2TypePrintNames{d} = 'WT';
        elseif strcmp(compName,'Short')
            data2TypePrintNames{d} = dataType3(6:end);
        elseif strcmp(compName,'Model')
            data2TypePrintNames{d} = modelPrintNames{d};
        end
        if strcmp(dataType2, 'Turndata_ShortABM_brown5')
            dataCol2 = [.5, .3, .1];
        end
        % Plot histogram:
        figure(figPEvenComp);
        if strcmp(compName,'Short')
            subplot(numel(dataTypes2),1,d);
        elseif strcmp(compName,'Model')
            subplot(numel(dataTypes2),1,d);
        end
        histogram( pEven_strct.yAfterCul.(dataType3), ...
            'Normalization', 'count', 'BinEdges',   -.025:.05:1.025, ...
            'FaceColor', dataCol2 );
        ylabel('#Flies');
        title( data2TypePrintNames{d} );
        xlabel('pEven (Cross. post cul-de-sac)');
        xticks(0:.25:1);
        % Plot in box plot:
        figure(figPEvenComp2);
        boxchart( d*ones(size(pEven_strct.yAfterCul.(dataType3))), ...
            pEven_strct.yAfterCul.(dataType3), 'BoxFaceColor', ...
            dataCol2, 'MarkerColor', dataCol2 );
        hold on;
    end
    
    % Save pEven histograms comparison figure;
    figure(figPEvenComp);
    figStartName = ['analyses/figures/flies_pEvenTrans_comp' compName ...
        'PEven' TransZeroName];
    saveas( figPEvenComp, [figStartName '.fig']);
    
    % Save pEven boxplot comparison figure;
    figure(figPEvenComp2);
    xticks(1:length(dataTypes2));
    xticklabels( data2TypePrintNames );
    xtickangle(0);
    ylim([0,1])
    hold on; 
    ggg = gca;
    plot( ggg.XLim, [.5,.5], 'k:', 'lineWidth', 1 );
    ylabel('pEven (Cross. post cul-de-sac)');
    figStartName = ['analyses/figures/flies_pEvenTrans_comp' compName ...
        'PEven2' TransZeroName];
    saveas( figPEvenComp2, [figStartName '.fig']);

end



%% Plot datasets comparison: #Transitions (post/within cul-de-sac) hist's:

% Plots #Transitions distribution across all trials in a dataset,
    % separately for motion post cul-de-sac and motion within the
    % cul-de-sac.

yRangesNames2 = {'culDeSac', 'yAfterCul'};
yRangesPrintNames2 = {'in cul-de-sac', 'post cul-de-sac'};

compNames = {'Short','Model'};
compVects = {[1,5,4,3,6], [14, 7:12, 13]};
compEdges.Short = [-.5:29.5, Inf];
compEdges.Model = [-.5:14.5, Inf];
yUpLims2.Short = [.15, .6];
yUpLims2.Model = [1, 1];

modelPrintNames = {'Brownian', 'No WF', 'WF 001', 'WF 002', 'WF 005', ...
    'WF 007', 'WF 009', 'WF'};

for c = 1:numel(compNames)
    compName = compNames{c};
    compVect = compVects{c};

    for r = 1:numel(yRangesNames2)
        yRangesName = yRangesNames2{r};
        yRangesPrintName = yRangesPrintNames2{r};
    
        figTRansComp = figure;
        dataTypes2 = dataTypes(compVect);
        dataCols2 = dataCols(compVect);
        data2TypePrintNames = cell(size(dataTypes2));
        for d = 1:numel(dataTypes2)

            dataType2 = dataTypes2{d};
            dataType3 = dataType2(10:end-2);
            numTransAllTrialsFlies = cell2mat( numTrans.(dataType3).(...
                yRangesName) );
            dataCol2 = dataCols2{d};
            if strcmp(compName,'Short')
                subplot(numel(dataTypes2),1,d);
            elseif strcmp(compName,'Model')
                subplot(numel(dataTypes2),1,d);
                %subplot(1,numel(dataTypes2),d);
            end
            if strcmp(dataType2, 'Turndata_ShortABM_brown5')
                dataCol2 = [.5, .3, .1];
            end
            if ~strcmp(dataType2,'Turndata_ShortABM_brown5')
                histogram( numTransAllTrialsFlies, 'normalization', ...
                    'probability', 'BinEdges', compEdges.(compName), ...
                    'FaceColor', dataCol2 ); 
                hold on;
                plot( repmat(compEdges.(compName)(1:end-1),[2,1]), ...
                    repmat([0;yUpLims2.(compName)(r)], ...
                    [1,length(compEdges.(compName))-1]), 'k');
            else
                histogram( numTransAllTrialsFlies, 'normalization', ...
                    'probability', 'FaceColor', dataCol2 ); 
            end
            ylabel('Prob.');
            xlabel(['#Crossings ' yRangesPrintName]);
            if strcmp(dataType2,'Turndata_ShortYy')
                title('WT');
            elseif strcmp(compName,'Short')
                title(dataType3(6:end));
            elseif strcmp(compName,'Model')
                title(modelPrintNames{d});
            end

            % Also compute pEven over all trials of all flies in dataset:
            pEvenAll = mean( rem(numTransAllTrialsFlies,2) == 0 );
            nAll = length(numTransAllTrialsFlies);
            pEvenAllTrialsInData.(yRangesName).(dataType3) = pEvenAll;
            pEvenAllTrialsInData_SEM.(yRangesName).(dataType3) = ...
                sqrt( pEvenAll * (1-pEvenAll) / nAll );
            % Repeat for p0:
            p0All = mean( numTransAllTrialsFlies == 0 );
            p0AllTrialsInData.(yRangesName).(dataType3) = p0All;
            p0AllTrialsInData_SEM.(yRangesName).(dataType3) = ...
                sqrt( p0All * (1-p0All) / nAll );

        end
        % Save pEven comparison figure;
        figStartName = ['analyses/figures/flies_pEvenTrans_comp' ...
            compName 'Trans_' yRangesName];
        saveas( figTRansComp, [figStartName '.fig']);
    end

end



%% Compare observed pEven vs expected under a Poisson distribution:

% Define datasets' names and col's:
dataTypes = {'Turndata_ShortYy', 'Turndata_ShortFoxPYy', ...
    'Turndata_ShortDumbYy', 'Turndata_ShortNorpAYy', ...
    'Turndata_ShortNompCYy', 'Turndata_LongYy', ...
    'Turndata_ShortABM_NoWF', ...
    'Turndata_ShortABM_WF_00100', 'Turndata_ShortABM_WF_00200', ...
    'Turndata_ShortABM_WF_00500', 'Turndata_ShortABM_WF_00700', ...
    'Turndata_ShortABM_WF_00900', 'Turndata_ShortABM_WallFollowing'};
dataNames = {'Short', 'FoxP','Dumb', 'NorpA', 'NompC', 'Long', ...
    'No WF', 'WF 001', 'WF 002', 'WF 005', 'WF 007', 'WF 009', 'WF'};
dataCols = {'k', 'g', 'm', 'r', '#77AC30', 'b', ...
    "#A2142F", ... % no WF 
    "#D95319", "#EDB120", "#77AC30", "#4DBEEE", "#0072BD", ... % WF parm.
    "#7E2F8E"}; % WF, brownian


% Create function: pEven as a fun. of lambda (mean #transitions):
pEvenPoisson = @(lambda)(0.5*(1+exp(-2*lambda)));

% Plot the function: 
figPoisson = figure;
lambda = 0:.01:5;
plot( lambda, pEvenPoisson(lambda), '-', 'Color', .5*[1,1,1], ...
    'lineWidth', 2 );
hold on;
plot( lambda([1,end]), [.5,.5], ':', 'Color', .5*[1,1,1], 'lineWidth', 1 );
hold on;

% Consider all trials in each dataset:
clear plt;
for d = 1:numel(dataTypes)
    dataType = dataTypes{d};
    dataCol = dataCols{d};
    % Compute the avg. number of transitions over the entire dataset:
    meanTransNumData = mean( cell2mat( ...
        numTrans.(dataType(10:end-2)).yAfterCul ) );
    % Compute obs. pEven, considering all trials in the dataset:
    pEvenData = mean( rem( cell2mat( ...
        numTrans.(dataType(10:end-2)).yAfterCul ), 2 ) == 0 );
    % Plot a line for the expected (Poisson) pEven given obs. <#trnas.>:
    %plot( meanTransNumData * [1,1], [0,1], '-', 'Color', dataCol );
    hold on;
    % Plot the observed pEven:
    if d<=6
        marShape = 'o';
    else
        marShape = 'd';
    end
    plt(d) = plot( meanTransNumData, pEvenData, marShape, ...
        'MarkerFaceColor', dataCol, 'MarkerEdgeColor', 'w' );
    hold on;
end
xlabel('Avg. #Crossings (or \lambda)');
ylabel('p_{Even}');
ylim([0,1]);
legend(plt([7:end,1:6]),dataNames([7:end,1:6]))

% save the figure:
figStartName = 'analyses/figures/flies_pEvenTransPoisson';
saveas( figPoisson, [figStartName '.fig']);



%% Plot hist. of avg.(#Transitions | -1<y<1) across WT flies in short maze:

dataType = 'Turndata_ShortYy';
dataCol = 'k';

flyAvgTransInBotArm = nan(numel(numTrans.(dataType(10:end-2)).yBeforeCul),1);

for nF = 1:numel(flyAvgTransInBotArm)
    flyAvgTransInBotArm(nF) =  mean( numTrans.(dataType(10:end-2)).yDownArm{nF} + ...
        numTrans.(dataType(10:end-2)).culDeSac{nF} + ...
        numTrans.(dataType(10:end-2)).yUpArm{nF} );
end

fig_HistAvgTrans = figure;
histogram( flyAvgTransInBotArm, 'Normalization', 'count', 'FaceColor', dataCol );
xlabel('Avg. #Crossings (bottom arm)');
ylabel('#Flies');
% Save figure:
figStartName = ['analyses/figures/flies_histAvgTransBottomArm_' ...
    dataType(10:end-2)];
saveas( fig_HistAvgTrans, [figStartName '.fig']);



%% Plots time in maze-region probability density estimate w/ Kernel smoothing:

% Plots the estimated probability density functions (PDFs) of the time 
% spent walking down the arm, within the cul-de-sac, walking up the arm
% and the fraction of time spent within the cul-de-sac (out of the total 
% time spent between -1<y<1) for WT flies in the short maze, WT flies in 
% the long maze and FoxP flies in the short maze.

fps = 30;

dataTypes = {'Turndata_ShortYy', 'Turndata_LongYy', 'Turndata_ShortFoxPYy'};
dataCols = {'k', 'b', 'g'};

figTimeInRegion = figure;

for d = 1:numel(dataTypes)
    dataType = dataTypes{d};
    dataCol = dataCols{d};
    % Load data (if exists) or run code to create data (otherwise):
    matFileName = ['output_drosophila_main_' dataType(10:end) '.mat'];
    load(matFileName);
    % Define cul-de-sac location:
    if strcmp(dataType,'Turndata_LongYy')
        culdesacY = 0.68/3.36;
    elseif strcmp(dataType,'Turndata_ShortYy') || ...
            strcmp(dataType,'Turndata_ShortNorpAYy') || ...
            strcmp(dataType,'Turndata_ShortDumbYy') || ...
            strcmp(dataType,'Turndata_ShortFoxPYy') 
        culdesacY = 0.69/2.03;
    elseif strcmp(dataType,'Turndata_ShortABM_NoWF') || ...
            strcmp(dataType,'Turndata_ShortABM_WallFollowing') || ...
            strcmp(dataType,'Turndata_ShortABM_WF_00100') || ...
            strcmp(dataType,'Turndata_ShortABM_WF_00200') || ...
            strcmp(dataType,'Turndata_ShortABM_WF_00500') || ...
            strcmp(dataType,'Turndata_ShortABM_WF_00700') || ...
            strcmp(dataType,'Turndata_ShortABM_WF_00900') || ...
            strcmp(dataType,'Turndata_ShortABM_brown5')
        culdesacY = 7.5/30;
    end
    nFramesInRegion.(dataType(10:end-2)).culAll = cell( ...
        numel(selectedFlies), 1 );
    nFramesInRegion.(dataType(10:end-2)).armUp = cell( ...
        numel(selectedFlies), 1 );
    nFramesInRegion.(dataType(10:end-2)).armDown = cell( ...
        numel(selectedFlies), 1 );
    nFramesInRegion.(dataType(10:end-2)).culAll_relative = cell( ...
        numel(selectedFlies), 1 );
    % Collect #frame in data: 
    for nF = 1:length(selectedFlies)
        f = selectedFlies(nF);
        flyName = ['f' num2str(f)];
        locsBotTrialT = meanX_vects_ALL.(flyName).locsBotTrialT;
        % Corrected trajectories in bottom trials: 
        xytyd = meanX_vects_ALL.(flyName).xytyd;
        xytyd_bottom = xytyd( locsBotTrialT );
        nFrames_flyTrial_culAll = nan(numel(xytyd_bottom),1);
        nFrames_flyTrial_armUp = nan(numel(xytyd_bottom),1);
        nFrames_flyTrial_armDown = nan(numel(xytyd_bottom),1);
        for t = 1:numel(xytyd_bottom)
            y_raw = xytyd_bottom{t}( :, 2 );
            nFrames_flyTrial_culAll(t) = sum( abs(y_raw)<culdesacY );
            nFrames_flyTrial_armUp(t) = sum( (y_raw>culdesacY) & ...
                (y_raw<1) );
            nFrames_flyTrial_armDown(t) = sum( (y_raw>-1) & ...
                (y_raw<-culdesacY) );
        end
        nFramesInRegion.(dataType(10:end-2)).culAll{nF} = ...
            nFrames_flyTrial_culAll;
        nFramesInRegion.(dataType(10:end-2)).armUp{nF} = ...
            nFrames_flyTrial_armUp;
        nFramesInRegion.(dataType(10:end-2)).armDown{nF} = ...
            nFrames_flyTrial_armDown;
        nFramesInRegion.(dataType(10:end-2)).culAll_relative{nF} = ...
            nFrames_flyTrial_culAll ./ (nFrames_flyTrial_culAll + ...
            nFrames_flyTrial_armUp + nFrames_flyTrial_armDown);
    end

    clear meanX_vects_ALL;

    culAll_datAll = cell2mat( nFramesInRegion.(dataType(10:end-2)).culAll );
    armUp_datAll = cell2mat( nFramesInRegion.(dataType(10:end-2)).armUp );
    armDown_datAll = cell2mat( nFramesInRegion.(dataType(10:end-2)).armDown );
    culAllRelative_datAll = cell2mat( nFramesInRegion.(dataType(10:end-2)...
        ).culAll_relative );

    subplot(1,4,1);
    theTs = armDown_datAll / fps;
    theTs = log10( theTs(theTs>0) );
    [fp,xfp] = ksdensity( theTs, 'function', 'pdf' );
    plot( xfp, fp, [dataCol '-'], 'lineWidth', 1 );
    xlim([-1,2]);
    xticks( -1:2 );
    xticklabels( 10 .^ (-1:2) );
    xlabel( 'time in region [s]' ); ylabel('pdf (kde)'); hold on;
    title('Down the arm');

    subplot(1,4,2);
    theTs = culAll_datAll / fps;
    theTs = log10( theTs(theTs>0) );
    [fp,xfp] = ksdensity( theTs, 'function', 'pdf' );
    plot( xfp, fp, [dataCol '-'], 'lineWidth', 1 );
    xlim([-1,2]);
    xticks( -1:2 );
    xticklabels( 10 .^ (-1:2) );
    xlabel( 'time in region [s]' ); ylabel('pdf (kde)'); hold on;
    title('Cul-de-sac');

    subplot(1,4,3);
    theTs = armUp_datAll / fps;
    theTs = log10( theTs(theTs>0) );
    [fp,xfp] = ksdensity( theTs, 'function', 'pdf' );
    plot( xfp, fp, [dataCol '-'], 'lineWidth', 1 );
    xlim([-1,2]);
    xticks( -1:2 );
    xticklabels( 10 .^ (-1:2) );
    xlabel( 'time in region [s]' ); ylabel('pdf (kde)'); hold on;    
    title('Up the arm');

    subplot(1,4,4);
    theTs = culAllRelative_datAll;
    theTs = theTs(theTs>0 & theTs<20);
    [fp,xfp] = ksdensity( theTs, 'Support', 'positive', 'function', 'pdf' );
    plot( xfp, fp, [dataCol '-'], 'lineWidth', 1 );
    xlabel( '% Time in cul-de-sac' ); ylabel('pdf (kde)'); hold on;    
    title('Cul-de-sac');

end


% Compute stats and tests:

armUp_datAll_short = cell2mat( nFramesInRegion.Short.armUp ) / fps;
armUp_datAll_long = cell2mat( nFramesInRegion.Long.armUp ) / fps;
[median(armUp_datAll_short), median(armUp_datAll_long)]
[p, ~, stats] = ranksum(armUp_datAll_short,armUp_datAll_long)

culAll_datAll_short = cell2mat( nFramesInRegion.Short.culAll ) / fps;
culAll_datAll_long = cell2mat( nFramesInRegion.Long.culAll ) / fps;
[median(culAll_datAll_short), median(culAll_datAll_long)]
[p, ~, stats] = ranksum(culAll_datAll_short,culAll_datAll_long)

armDown_datAll_short = cell2mat( nFramesInRegion.Short.armDown ) / fps;
armDown_datAll_long = cell2mat( nFramesInRegion.Long.armDown ) / fps;

perCul_short = culAll_datAll_short ./ (culAll_datAll_short + ...
    armUp_datAll_short + armDown_datAll_short);
perCul_long = culAll_datAll_long ./ (culAll_datAll_long + ...
    armUp_datAll_long + armDown_datAll_long);
[median(perCul_short), median(perCul_long)]
[p, ~, stats] = ranksum(perCul_short,perCul_long)

% Save the figure:
figStartName = 'analyses/figures/flies_pEvenTrans_TimeInRegion';
saveas( figTimeInRegion, [figStartName '.fig']);


