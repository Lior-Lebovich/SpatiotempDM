%% DATASET AVERAGE QUIVER PLOTS AND HEATMAPS:
%
% Computes and plots average quiver plots across flies in each fly dataset
% within the cul-de-sac, separately for left and right turns and motion
% towards or away from the cul-de-sac. 
% Also computes and plot average quivers exclusively for zero post 
% cul-de-sac midline-crossing trials for each fly dataset. 
% Computes and plots average heatmaps across decision-makers, separately 
% for left and right turns and the corresponding differnce. 
%
% INPUT FILES: 
%   output_drosophila_main_[DATASET_NAME(10:END)].mat
%   flies_boundPolProc_[DATA_NAME].mat
%   flies_pEvenTrans.mat
%
% FUNCTION:
%   quiverGivenTurns.m: computing a density-based quiver of an individual 
%   decision-maker, for a given maze-region of interest.
%
% EXTERNAL FUNCTION:
%   This code uses an external function, patchline.m, Reference: 
%   Brett Shoelson (2023). Patchline 
%   (https://www.mathworks.com/matlabcentral/fileexchange/36953-patchline), 
%   MATLAB Central File Exchange. Retrieved September 18, 2023.
%
% OUTPUT FILE: 
%   output_drosophila_Quivers.mat: with structures binsAngRad_dat and 
%       binsFreq_dat storing avg. directions and magnitudes for 
%       .[DATANAME].[TURN_DIRECTION].[WALKING_DIRECTION]
%
% FIGURES:
%   Average quiver plots for each fly datatype (averaged over flies in
%       datatype), separately for left/right turn and up/down walking 
%       directions:
%       Figs. 4D, S3F, S7: flies[DATA_NAME]_avgQuiver.fig
%   Heatmaps depicting bivariate probability(x,y) averaged over flies or 
%       agents in each fly or ABM dataset, separately for left/right turns, 
%       together with corresp. difference:
%       Figs. 2B, S13F: flies[DATA_NAME]_avgHeat.fig
%   Average quiver plots for each fly datatype, for zero post cul-de-sac 
%       midline-crossing trials:
%       Fig. S10: flies[DATA_NAME]_avgQuiverZeroTrans.fig
%
% For detailed information about how the quiver of an individual fly is 
% computed, visit: quiverGivenTurns.m (function) and
% drosophila_quiverIndivs.m (script).
% Note in the current midline-crossings are termed "transitions".



%% Use relevant directory and add path to external fucntions:

cdName = 'SpatiotempDM';
currentFolder = cd;
if ~endsWith(currentFolder,cdName)
    cd(cdName);
end
addpath( 'analyses', 'datafiles', 'analyses/customfuns', '-frozen');



%% Compute and plot avg. quiver over flies in each fly dataset:
 
% Average quivers are plotted separately for upward and downward motion
% (positive and negative y values) and separately for left and right turns.
% The current custom definition below, plot the avg. quiver with the
% cul-de-sac and the RoI (can be manually changed to 'intersection'),
% without any bckground image (can be changed: "overlaidOnHeat = 'yes'")
% and with an additional panel with right and left turn quivers overlaid
% (can be changed to all turns: "bothInterp = 'allData';").

% CUSTOM DEFINITIONS FOR PLOTS:

% To plot quivers for the intersection rather than cul-de-sac, change the 
    % value of regionName below from 'culdesac' to 'intersection':
regionName = 'culdesac';

% To plot overlaid on heat maps (redundant) change the value of 
    % overlaidOnHeat below from 'no' to 'yes':
overlaidOnHeat = 'no';

% To plot quivers for data of left, right and all turns (istead of left, 
    % right and left&right overlaid) change the value of bothInterp below 
    % from leftVsRight to 'allData'.
bothInterp = 'leftVsRight';


% Define turn-dir. names and colors:
decVals = [-1,1,10]; % Left, right, left and right turns
decCols = {'k','k','k'};
decCols2 = {'b','r','k'};
turnNames = {'LEFT', 'RIGHT', 'BOTH'};
if strcmp( overlaidOnHeat, 'no')
    boundCol = [0,0,0];
    decCols = decCols2;
elseif strcmp( overlaidOnHeat, 'yes')
    boundCol = [1,1,1];
end
walkDirs = {'down','up'};
walkDirPols = [-1,1];

% Define datasets:
dataTypes = {'Turndata_ShortYy', 'Turndata_LongYy', ...
    'Turndata_ShortNorpAYy', 'Turndata_ShortDumbYy', ...
    'Turndata_ShortFoxPYy', 'Turndata_ShortNompCYy'};

% Define subsampling rate:
subsamp = 1; % no subsampling


for d = 1:numel(dataTypes)
    dataType = dataTypes{d};
    
    % Load processed x,y turn decisions in dataset (output of 
    % drosophila_main.m):
    matFileName = ['output_drosophila_main_' dataType(10:end) '.mat'];
    load(matFileName);

    % Load avg. bounding polygon of all animals in dataset:
    if strcmp(dataType,'Turndata_ShortYy') || ...
            strcmp(dataType,'Turndata_LongYy') || ...
            strcmp(dataType,'Turndata_ShortNorpAYy') || ...
            strcmp(dataType,'Turndata_ShortDumbYy') || ...
            strcmp(dataType,'Turndata_ShortFoxPYy') || ...
            strcmp(dataType,'Turndata_ShortNompCYy')
        matFileName2 = ['flies_boundPolProc_' dataType(10:end-2) '.mat'];
        load( ['flies_boundPolProc_' dataType(10:end-2) '.mat'], ...
            'polyNormAvg' );
        polyDataset = polyNormAvg;
    else
        % Load the simulation's single bounding polygon:
        pol = polyAllFlies1; 
        polX = pol(:,1);
        polY = pol(:,2);
        minYPolyAllSims = min(polY);
        % Compute fly's bottom arm length WITH intersection:
        lenBotFlyAllSims = 34.7; 
        widCenFlyAllSims = .5 * max(polX);
        % Plot corrected bounding polygon of all sim. flies:
        polyDataset_x = (polX - widCenFlyAllSims ) / ...
            ( lenBotFlyAllSims * ratioLenWoW );
        polyDataset_y = (polY - minYPolyAllSims ) / ...
            ( lenBotFlyAllSims * ratioLenWoW );
        polyDataset = [polyDataset_x, polyDataset_y]; 
    end
    

    % Determine 2D grid for bins, based on maze size:

    % Read cul-de-sac ratios:
    if strcmp(dataType,'Turndata_LongYy')
        culDeSacEdge = 0.68/3.36;
        culPrcnt = .67;
        binMult = 1;
        yLimCul = .1;
    elseif strcmp(dataType,'Turndata_ShortYy') || ...
            strcmp(dataType,'Turndata_ShortNorpAYy') || ...
            strcmp(dataType,'Turndata_ShortDumbYy') || ...
            strcmp(dataType,'Turndata_ShortFoxPYy') || ...
            strcmp(dataType,'Turndata_ShortNompCYy')
        culDeSacEdge = 0.69/2.03;
        culPrcnt = .67;
        binMult = 1;
        yLimCul = .175;
    elseif strcmp(dataType,'Turndata_ShortABM_NoWF') || ...
            strcmp(dataType,'Turndata_ShortABM_WallFollowing') || ...
            strcmp(dataType,'Turndata_ShortABM_Brownian') || ...
            strcmp(dataType,'Turndata_ShortABM_WF_00100') || ...
            strcmp(dataType,'Turndata_ShortABM_WF_00200') || ...
            strcmp(dataType,'Turndata_ShortABM_WF_00500') || ...
            strcmp(dataType,'Turndata_ShortABM_WF_00700') || ...
            strcmp(dataType,'Turndata_ShortABM_WF_00900') || ...
            strcmp(dataType,'Turndata_ShortABM_brown5')
        culDeSacEdge = 7.5/30;
        culPrcnt = 1;
        binMult = 1.5;
        yLimCul = .15;
    end

    % Define X-edges: based on maze size:
    binWidthShort = .0162;
    binWidth = binMult * ratioSize * binWidthShort;
    xEdgesOneSide = 0:binWidth:culPrcnt*culDeSacEdge; 
    xEdges = [-flip(xEdgesOneSide(2:end)), xEdgesOneSide];

    % Define Y-edges: based on maze size:
    if strcmp(regionName,'culdesac')
        yEdgesRaw = 0:binWidth:1.05*culDeSacEdge; 
    elseif strcmp(regionName,'intersection')
        yEdgesRaw = linspace( 1, 1.05/ratioLenWoW, round(16/ratioSize) );
    end


    % Run over walking directions (down/up the arm) to compute the
    % left/right/both quivers separately for each: 

    dataTypeQuiverFig = figure( 'Units', 'Normalized', 'Position', ...
        [0 .3 1 .7] );

    for wd = 1:numel(walkDirs)

        walkDirPol = walkDirPols(wd);
        
        % Use y-edges according to walking direction (down/up the arm):
        walkDir = walkDirs{wd};
        if strcmp( walkDir, 'down' )
            yEdges = -1  * flip(yEdgesRaw);
        else
            yEdges = yEdgesRaw;
        end

        % For edges correction:
        dxEdges = diff( xEdges(1:2) );
        dyEdges = binWidth;
        
        % Create grid for 2D (x,y) bins:
        [X, Y] = meshgrid( xEdges, yEdges );


        % Create arrays for storing flys' magnitude-based avg. dir. and 
        % freq. (normalized counts) for each fly and each turn dir.:
        for ss = 1:length(decVals)
            binsDirRad.([dataType(10:end-2)]).(turnNames{ss}).(walkDir...
                ) = nan(size(X,1),size(X,2),numel(selectedFlies));
            binsFreq.([dataType(10:end-2)]).(turnNames{ss}).(walkDir) = ...
                nan(size(X,1),size(X,2),numel(selectedFlies));
        end
        xForCentering = nan( numel(selectedFlies), 1 );

        
        % Run over flies in dataset and compute the quiver of each fly:
    
        for nF = 1:length(selectedFlies)
            f = selectedFlies(nF);
            flyName = ['f' num2str(f)];
    
            % Read individual animal's processed data for bottom trials: 
            
            % Location of bottom trials:
            locsBotTrialT = meanX_vects_ALL.(flyName).locsBotTrialT;
            % Corrected trajectories in bottom trials: 
            xytyd = meanX_vects_ALL.(flyName).xytyd;
            xytyd_bottom = xytyd( locsBotTrialT );
            % Decision vector for bototm trials:
            Dec = meanX_vects_ALL.(flyName).Dec;
    
            % Compute #turn to each side from bottom arm:
            nTrialsLeftFly = sum( locsBotTrialT & Dec == -1);
            nTrialsRightFly = sum( locsBotTrialT & Dec == 1);
        
            % Compute quiver:
            if 1==1 
    
                % Find fly's med x loc.:
                traj_cell_mat = cell2mat( xytyd_bottom );
                xy_bottom_mat = traj_cell_mat(:,1:2);
                uniqueXY = unique( xy_bottom_mat, 'rows' );
                xGoodMedian_allTrajs = mean( ...
                    [min( uniqueXY( abs(uniqueXY(:,2)) < 1, 1 ) ), ...
                    max( uniqueXY( abs(uniqueXY(:,2)) < 1, 1 ) )] );
                % Store xMed required for x-centering:
                xForCentering(nF) = xGoodMedian_allTrajs;
    
                for s = 1:length(decVals)
                    % Read turn-dir. value and colors:
                    decVal = decVals(s);
                    turnName = turnNames{s};
            
                    % Use only fly's bottom trials that ended in specific 
                    % turn-dir.:
                    if decVal==10
                        xytydBotDec = xytyd( locsBotTrialT );   
                        TurnBotDec = Dec( locsBotTrialT );
                    else
                        xytydBotDec = xytyd( locsBotTrialT & Dec == ...
                            decVal);   
                        TurnBotDec = Dec( locsBotTrialT & Dec == decVal);   
                    end
    
                    % Compute counts and magnitude-based avg. dir. in 
                    % Radians for each 2D bin:
                    [counts, avgDirThetaRad] = quiverGivenTurns( X, Y, ...
                        xEdges, yEdges, xytydBotDec, ...
                        xGoodMedian_allTrajs, subsamp );
    
                    % Store fly's magnitude-based avg. dir. and freq.
                    % (normalized counts):
                    binsDirRad.([dataType(10:end-2)]).(turnName).(...
                        walkDir)(:,:,nF) = avgDirThetaRad;
                    binsFreq.([dataType(10:end-2)]).(turnName).(...
                        walkDir)(:,:,nF) = counts / sum(counts(:));
    
                end
    
            end
    
        end


        % Rerun over turn-directions --> compute & plot average quiver for 
        % each [note that this is conditioned on walking direction]
    
        for ssd = 1:length(decVals) 
    
            turnName = turnNames{ssd};
            decCol = decCols{ssd};
    
    
            % COMPUTE AVERAGE QUIVER IN DATATYPE * TURN-DIRECTION: 
    
            % Read the avg. dir's and freq's of all flies in all bins:
            binFreq_all = binsFreq.([dataType(10:end-2)]).(turnName).(...
                walkDir);
            binsDirRad_all = binsDirRad.([dataType(10:end-2)]).(...
                turnName).(walkDir);
            
            % Use complex representation for angles in bins:
            binsComplexAng_all =  exp(1i * binsDirRad_all); 
        
            % Compute avg. direction. For each bin, this is the flies' 
            % frequency-based weighted avg. of the angles (in complex 
            % representation):
            binsWeightedMeanAngleComplex = sum( binsComplexAng_all .* ...
                binFreq_all, 3, 'omitnan' ) ./ ...
                sum( binFreq_all, 3, 'omitnan' );
        
            % Convert avg. dir. to angle in radians:
            binsWeightedMeanAngleRad = angle(binsWeightedMeanAngleComplex);
        
            % Compute AVG. magnitude. For each bin, this is the relative 
            % frequency of visits of that bin, averaged over all flies.
            binFreq_allAvg = mean( binFreq_all, 3, 'omitnan' );
            binFreq_allNormAvg = binFreq_allAvg / sum( binFreq_allAvg(:) );
        
        
            % PLOT AVERAGE QUIVER:
        
            currSubPlot = subplot(1,2*length(decVals),2*(ssd-1)+wd);
            
            % Plot heat-map (bivariate PDF):
            if strcmp( overlaidOnHeat, 'yes')
                imagesc( (.5*dyEdges) + yEdges, (.5*dxEdges) + xEdges, ...
                    binFreq_allAvg' / sum(binFreq_allAvg(:)) );
                colorbar;
                set(gca, 'YDir', 'Normal');
                hold on;
            end
    
            % Also plot dataset avg. bounding polygons (if behavioral 
            % data):
            if ~contains( dataType, 'ABM')
                patchline( walkDirPol * polyDataset(:,2), ...
                    polyDataset(:,1), 'linestyle', '-', 'edgecolor', ...
                    boundCol, 'linewidth', 1, 'edgealpha', 0.5 ); 
                hold on;
            end
            
            % Plot the avg. vector field of all flies' trials w/ specific 
            % turn dir. [BOTH: plots either the entire data or left/right
            % overlaid]:
            if strcmp(turnName,'BOTH') && strcmp(bothInterp, 'leftVsRight')
                copyobj( quiveObj.LEFT.(walkDir), currSubPlot );
                copyobj( quiveObj.RIGHT.(walkDir), currSubPlot );
                hold on;
            else
                quiveObj.(turnName).(walkDir) = quiver( ...
                    (.5*dyEdges) + Y, (.5*dxEdges) + X, ...
                    ( binFreq_allAvg / sum(binFreq_allAvg(:)) ) .* ...
                    sin( binsWeightedMeanAngleRad ), ...
                    ( binFreq_allAvg / sum(binFreq_allAvg(:)) ) .* ...
                    cos( binsWeightedMeanAngleRad ), ...
                    decCol, 'LineWidth', 1 ); 
                hold on;
                % Store dataset's avg. Freq and angle:
                binsAngRad_dat.(dataType(10:end-2)).(turnName).(walkDir...
                    ) = binsWeightedMeanAngleRad;
                binsFreq_dat.(dataType(10:end-2)).(turnName).(walkDir...
                    ) = ( binFreq_allAvg / sum(binFreq_allAvg(:)) );
            end
    
            % Define subplot's axes:
            axis image;
            if strcmp(regionName,'culdesac')
                ylim( yLimCul * [-1,1]);
            elseif strcmp(regionName,'intersection')
                ylim(xEdges([1,end]));
            end
            xlim(yEdges([1,end]));
            xlabel('y');
            ggg = gca;
            ggg.YDir = 'reverse';
            ggg.GridAlpha = 0.2;
            grid on;
            if strcmp(walkDir,'down')
                ylabel('x');
            elseif strcmp(walkDir,'up')
                yticklabels([]);
            end
            if strcmp(turnName,'BOTH') && strcmp(bothInterp, 'leftVsRight')
                title([dataType(10:end-2) ' - LEFT vs RIGHT']);
            else
                title([dataType(10:end-2) ' - ' turnName;]);
            end

        end


    end
    

    % Save the figure:
    figure(dataTypeQuiverFig);
    figStartName = ['analyses/figures/flies' dataType(10:end-2) ...
        '_avgQuiver'];
    saveas( dataTypeQuiverFig, [figStartName '.fig']);
   
    clear meanX_vects_ALL quiveObj;
    
end

% Store datasets' avg. angles and magnitudes (freq.):
save('datafiles/output_drosophila_Quivers.mat', 'binsAngRad_dat', ...
    'binsFreq_dat' );



%% Avg. heat-maps over flies in each fly or ABM dataset:

dataTypes = {'Turndata_ShortYy', 'Turndata_LongYy', ...
    'Turndata_ShortNorpAYy', 'Turndata_ShortDumbYy', ...
    'Turndata_ShortFoxPYy', 'Turndata_ShortNompCYy', ...
    'Turndata_ShortABM_NoWF', 'Turndata_ShortABM_WallFollowing', ...
    'Turndata_ShortABM_WF_00100', ...
    'Turndata_ShortABM_WF_00200', 'Turndata_ShortABM_WF_00500', ...
    'Turndata_ShortABM_WF_00700','Turndata_ShortABM_WF_00900',  ...
    'Turndata_ShortABM_brown5'}; 

decVals = [10,-1,1]; % Left and right turns
decCols = {'k','k','k'};
turnNames = {'BOTH', 'LEFT', 'RIGHT'};
% Also plot delta(prob(x,y)):
turnNames2 = {'BOTH', 'LEFT', 'RIGHT', 'RIGHT-LEFT'};

% Define subsampling rate:
subsamp = 1; % no subsampling

for d = 1:numel(dataTypes) 
    dataType = dataTypes{d};
    
    % Load data (if exists) or run code to create data (otherwise):
    matFileName = ['output_drosophila_main_' dataType(10:end) '.mat'];
    load(matFileName);
    
    % Read cul-de-sac ratios:
    if strcmp(dataType,'Turndata_LongYy')
        culDeSacEdge = 0.68/3.36;
    elseif strcmp(dataType,'Turndata_ShortYy') || ...
            strcmp(dataType,'Turndata_ShortNorpAYy') || ...
            strcmp(dataType,'Turndata_ShortDumbYy') || ...
            strcmp(dataType,'Turndata_ShortFoxPYy') || ...
            strcmp(dataType,'Turndata_ShortNompCYy')
        culDeSacEdge = 0.69/2.03;
    elseif strcmp(dataType,'Turndata_ShortABM_NoWF') || ...
            strcmp(dataType,'Turndata_ShortABM_WallFollowing') || ...
            strcmp(dataType,'Turndata_ShortABM_Brownian') || ...
            strcmp(dataType,'Turndata_ShortABM_WF_00100') || ...
            strcmp(dataType,'Turndata_ShortABM_WF_00200') || ...
            strcmp(dataType,'Turndata_ShortABM_WF_00500') || ...
            strcmp(dataType,'Turndata_ShortABM_WF_00700') || ...
            strcmp(dataType,'Turndata_ShortABM_WF_00900') || ...
            strcmp(dataType,'Turndata_ShortABM_brown5')
        culDeSacEdge = 7.5/30;
    end

    % Determine 2D grid for bins, based on maze size:

    % Define X-edges and Y-edges: based on maze size:
    binWidthShort = 2 * .0162;
    binWidth = ratioSize * binWidthShort;
    %{
    xEdges = -1.2:binWidth:1.2; 
    yEdges = -1.5:binWidth:1.5; 
    %}
    xEdges = [-flip(binWidth:binWidth:.2), 0:binWidth:.2];  % -.2:binWidth:.2; 
    yEdges = 0:binWidth:1.5; %-.3:binWidth:.5; 

    % For edges correction:
    dxEdges = diff( xEdges(1:2) );
    dyEdges = diff( yEdges(1:2) );
    
    % Create grid for 2D (x,y) bins:
    [X, Y] = meshgrid( xEdges, yEdges );


    % Create arrays for storing flys' magnitude-based avg. dir. and freq.
    % (normalized counts) for each fly and each turn dir.:
    for ss = 1:length(decVals)
        binsProb.([dataType(10:end-2)]).(turnNames{ss}) = ...
            nan(size(X,1)-1,size(X,2)-1,numel(selectedFlies));
    end
    xForCentering = nan( numel(selectedFlies), 1 );

    
    % Run over flies in dataset and compute the quiver of each fly:

    for nF = 1:length(selectedFlies)
        f = selectedFlies(nF);
        flyName = ['f' num2str(f)];

        % Read individual animal's processed data for bottom trials: 
        
        % Location of bottom trials:
        locsBotTrialT = meanX_vects_ALL.(flyName).locsBotTrialT;
        % Corrected trajectories in bottom trials: 
        xytyd = meanX_vects_ALL.(flyName).xytyd;
        xytyd_bottom = xytyd( locsBotTrialT );
        % Decision vector for bototm trials:
        Dec = meanX_vects_ALL.(flyName).Dec;

        % Compute #turn to each side from bottom arm:
        nTrialsLeftFly = sum( locsBotTrialT & Dec == -1);
        nTrialsRightFly = sum( locsBotTrialT & Dec == 1);
    
        % Compute qiuver IFF, for bottom arm, #leftTurn>=30 and 
        % #rightTurn>=30:
        if 1==1 %%% nTrialsLeftFly>=30 && nTrialsRightFly>=30

            % Find fly's med x loc.:
            traj_cell_mat = cell2mat( xytyd_bottom );
            xy_bottom_mat = traj_cell_mat(:,1:2);
            uniqueXY = unique( xy_bottom_mat, 'rows' );
            xGoodMedian_allTrajs = mean( ...
                [min( uniqueXY( abs(uniqueXY(:,2)) < 1 & ...
                abs(uniqueXY(:,2)) > culDeSacEdge, 1 ) ), ...
                max( uniqueXY( abs(uniqueXY(:,2)) < 1 & ...
                abs(uniqueXY(:,2)) > culDeSacEdge, 1 ) )] );
            % Store xMed required for x-centering:
            xForCentering(nF) = xGoodMedian_allTrajs;

            for s = 1:length(decVals)
                % Read turn-dir. value and colors:
                decVal = decVals(s);
                turnName = turnNames{s};
        
                % Use only fly's bottom trials that ended in specific turn-dir.:
                if decVal==10
                    xytydBotDec = xytyd( locsBotTrialT );   
                    TurnBotDec = Dec( locsBotTrialT );
                else
                    xytydBotDec = xytyd( locsBotTrialT & Dec == decVal);   
                    TurnBotDec = Dec( locsBotTrialT & Dec == decVal);   
                end

                % Read x,y for specific turn-dir. and center x:
                xyElseNotCentered = cell2mat( xytydBotDec );
                xyCentered = [...
                    xyElseNotCentered(:,1)-xGoodMedian_allTrajs, ...
                    xyElseNotCentered(:,2)];
                % Compute bivariate Prob.:
                [N_dec,~,~,~,~] = histcounts2( xyCentered(:,1), ...
                    xyCentered(:,2), xEdges, yEdges, 'Normalization', ...
                    'probability' );

                % Store fly's bivariate Prob.::
                binsProb.([dataType(10:end-2)]).(turnName)(:,:,nF) = N_dec';

            end

        end

    end


    % Rerun over turn-directions --> compute & plot average heat-maps:
    
    dataTypeHeatFig = figure( 'Units', 'Normalized', 'Position', ...
            [0 .3 1 .7] );

    for ssd = 1:numel(turnNames2) 

        turnName2 = turnNames2{ssd};


        % COMPUTE AVERAGE BIVARIATE PROB. IN DATATYPE * TURN-DIRECTION: 

        % Read the avg. dir's and freq's of all flies in all bins:
        if strcmp(turnName2,'RIGHT-LEFT')
            binFreq_R = binsProb.([dataType(10:end-2)]).RIGHT;
            binFreq_L = binsProb.([dataType(10:end-2)]).LEFT;
            binFreq_allAvg = mean( binFreq_R, 3, 'omitnan' ) - ...
                mean( binFreq_L, 3, 'omitnan' );
        else
            binFreq_all = binsProb.([dataType(10:end-2)]).(turnName2);
            binFreq_allAvg = mean( binFreq_all, 3, 'omitnan' );
        end
    
        % Compute AVG. bivariate Prob. for specific dec:
        
        
    
        % PLOT AVG. bivariate Prob. for specific dec:
    
        subplot(2,2,ssd);
        
        % Plot bivariate PDF:
        normPDF = binFreq_allAvg' / abs(sum(binFreq_allAvg(:)));
        if strcmp(turnName2,'RIGHT-LEFT')
            maxAbsBin = max( abs( [min(normPDF(:)), max(normPDF(:))] ) );
            %{
            imagesc( (.5*dyEdges) + yEdges, (.5*dxEdges) + xEdges, ...
                normPDF, maxAbsBin * [-1,1] );
            %}
            g = imagesc( yEdges, xEdges, normPDF, maxAbsBin * [-1,1] );
            %g.Interpolation = 'bilinear';
            cc = colorbar;
            %maxCBLim = max( abs([cc.Limits(1), cc.Limits(2)]) );
            %cc.Limits = maxCBLim * [-1,1];
        else
            %{
            imagesc( (.5*dyEdges) + yEdges, (.5*dxEdges) + xEdges, ...
                normPDF );
            %}
            g = imagesc( yEdges, xEdges, normPDF );
            %g.Interpolation = 'bilinear';
            cc = colorbar;
        end
        hold on;

        % Define subplot's axes:
        axis image;
        ylim(xEdges([1,end]));
        xlim(yEdges([1,end]));
        xlabel('y');
        ylabel('x');
        title([dataType(10:end-2) ' - ' turnName2]);
        ggg = gca;
        ggg.YDir = 'reverse';
        ggg.GridAlpha = 0.2;
        grid on;
    end

    % Save the figure:
    figure(dataTypeHeatFig);
    figStartName = ['analyses/figures/flies' dataType(10:end-2) '_avgHeat'];
    saveas( dataTypeHeatFig, [figStartName '.fig']);
    clear meanX_vects_ALL;
    
end



%% Computes and plots avg. quiver | #Trans.=0 for each fly dataset:
 
% Average quivers (as the first quiver section above), only for trials with 
% zero post cul-de-sac transitions. 
% As in the previous quiver section, the RoI, background and third panel 
% may be manually modified below.

% CUSTOM DEFINITIONS FOR PLOTS:

% To plot quivers for the intersection rather than cul-de-sac, change the 
    % value of regionName below from 'culdesac' to 'intersection':
regionName = 'culdesac';

% To plot overlaid on heat maps (redundant) change the value of 
    % overlaidOnHeat below from 'no' to 'yes':
overlaidOnHeat = 'no';

% To plot quivers for data of left, right and all turns (istead of left, 
    % right and left&right overlaid) change the value of bothInterp below 
    % from leftVsRight to 'allData'.
bothInterp = 'leftVsRight';


% Define turn-dir. names and colors:
decVals = [-1,1,10]; % Left, right, left and right turns
decCols = {'k','k','k'};
decCols2 = {'b','r','k'};
turnNames = {'LEFT', 'RIGHT', 'BOTH'};
if strcmp( overlaidOnHeat, 'no')
    boundCol = [0,0,0];
    decCols = decCols2;
elseif strcmp( overlaidOnHeat, 'yes')
    boundCol = [1,1,1];
end
walkDirs = {'down','up'};
walkDirPols = [-1,1];

% Define datasets:
dataTypes = {'Turndata_ShortYy', 'Turndata_LongYy', ...
    'Turndata_ShortNorpAYy', 'Turndata_ShortDumbYy', ...
    'Turndata_ShortFoxPYy', 'Turndata_ShortNompCYy'};

% Define subsampling rate:
subsamp = 1; % no subsampling

% Load #Transitions made by each trial and each fly:
load( 'flies_pEvenTrans.mat', 'numTrans' );

for d = 1:numel(dataTypes)
    dataType = dataTypes{d};
    
    % Load processed x,y turn decisions in dataset (output of 
    % drosophila_main.m):
    matFileName = ['output_drosophila_main_' dataType(10:end) '.mat'];
    load(matFileName);

    % Load avg. bounding polygon of all animals in dataset:
    if strcmp(dataType,'Turndata_ShortYy') || ...
            strcmp(dataType,'Turndata_LongYy') || ...
            strcmp(dataType,'Turndata_ShortNorpAYy') || ...
            strcmp(dataType,'Turndata_ShortDumbYy') || ...
            strcmp(dataType,'Turndata_ShortFoxPYy') || ...
            strcmp(dataType,'Turndata_ShortNompCYy')
        matFileName2 = ['flies_boundPolProc_' dataType(10:end-2) '.mat'];
        load( ['flies_boundPolProc_' dataType(10:end-2) '.mat'], ...
            'polyNormAvg' );
        polyDataset = polyNormAvg;
    else
        % Read the simulation's single bounding polygon:
        pol = polyAllFlies1; 
        polX = pol(:,1);
        polY = pol(:,2);
        minYPolyAllSims = min(polY);
        % Compute fly's bottom arm length WITH intersection:
        lenBotFlyAllSims = 34.7; 
        widCenFlyAllSims = .5 * max(polX);
        % Plot corrected bounding polygon of all sim. flies:
        polyDataset_x = (polX - widCenFlyAllSims ) / ...
            ( lenBotFlyAllSims * ratioLenWoW );
        polyDataset_y = (polY - minYPolyAllSims ) / ...
            ( lenBotFlyAllSims * ratioLenWoW );
        polyDataset = [polyDataset_x, polyDataset_y]; 
    end
    

    % Determine 2D grid for bins, based on maze size:

    % Read cul-de-sac ratios:
    if strcmp(dataType,'Turndata_LongYy')
        culDeSacEdge = 0.68/3.36;
        culPrcnt = .67;
        binMult = 1;
        yLimCul = .1;
    elseif strcmp(dataType,'Turndata_ShortYy') || ...
            strcmp(dataType,'Turndata_ShortNorpAYy') || ...
            strcmp(dataType,'Turndata_ShortDumbYy') || ...
            strcmp(dataType,'Turndata_ShortFoxPYy') || ...
            strcmp(dataType,'Turndata_ShortNompCYy')
        culDeSacEdge = 0.69/2.03;
        culPrcnt = .67;
        binMult = 1;
        yLimCul = .175;
    elseif strcmp(dataType,'Turndata_ShortABM_NoWF') || ...
            strcmp(dataType,'Turndata_ShortABM_WallFollowing') || ...
            strcmp(dataType,'Turndata_ShortABM_Brownian') || ...
            strcmp(dataType,'Turndata_ShortABM_WF_00100') || ...
            strcmp(dataType,'Turndata_ShortABM_WF_00200') || ...
            strcmp(dataType,'Turndata_ShortABM_WF_00500') || ...
            strcmp(dataType,'Turndata_ShortABM_WF_00700') || ...
            strcmp(dataType,'Turndata_ShortABM_WF_00900') || ...
            strcmp(dataType,'Turndata_ShortABM_brown5')
        culDeSacEdge = 7.5/30;
        culPrcnt = 1;
        binMult = 1.5;
        yLimCul = .15;
    end

    % Define X-edges: based on maze size:
    binWidthShort = .0162;
    binWidth = binMult * ratioSize * binWidthShort;
    xEdgesOneSide = 0:binWidth:culPrcnt*culDeSacEdge; 
    xEdges = [-flip(xEdgesOneSide(2:end)), xEdgesOneSide];

    % Define Y-edges: based on maze size:
    if strcmp(regionName,'culdesac')
        yEdgesRaw = 0:binWidth:1.05*culDeSacEdge; 
    elseif strcmp(regionName,'intersection')
        yEdgesRaw = linspace( 1, 1.05/ratioLenWoW, round(16/ratioSize) );
    end


    % Run over walking directions (down/up the arm) to compute the
    % left/right/both quivers separately for each: 

    dataTypeQuiverFig = figure( 'Units', 'Normalized', 'Position', ...
        [0 .3 1 .7] );

    for wd = 1:numel(walkDirs)

        walkDirPol = walkDirPols(wd);
        
        % Use y-edges according to walking direction (down/up the arm):
        walkDir = walkDirs{wd};
        if strcmp( walkDir, 'down' )
            yEdges = -1  * flip(yEdgesRaw);
        else
            yEdges = yEdgesRaw;
        end

        % For edges correction:
        dxEdges = diff( xEdges(1:2) );
        dyEdges = binWidth;
        
        % Create grid for 2D (x,y) bins:
        [X, Y] = meshgrid( xEdges, yEdges );


        % Create arrays for storing flys' magnitude-based avg. dir. and 
        % freq. (normalized counts) for each fly and each turn dir.:
        for ss = 1:length(decVals)
            binsDirRad.([dataType(10:end-2)]).(turnNames{ss}).(walkDir) ...
                = nan(size(X,1),size(X,2),numel(selectedFlies));
            binsFreq.([dataType(10:end-2)]).(turnNames{ss}).(walkDir) ...
                = nan(size(X,1),size(X,2),numel(selectedFlies));
        end
        xForCentering = nan( numel(selectedFlies), 1 );

        
        % Run over flies in dataset and compute the quiver of each fly:
    
        for nF = 1:length(selectedFlies)
            f = selectedFlies(nF);
            flyName = ['f' num2str(f)];
    
            % Read individual animal's processed data for bottom trials: 
            
            % Location of bottom trials:
            locsBotTrialT = meanX_vects_ALL.(flyName).locsBotTrialT;
            % Corrected trajectories in bottom trials: 
            xytyd = meanX_vects_ALL.(flyName).xytyd;
            xytyd_bottom = xytyd( locsBotTrialT );
            % Decision vector for bototm trials:
            Dec = meanX_vects_ALL.(flyName).Dec;
            Dec_bottom = Dec(locsBotTrialT);

    
            % Compute #turn to each side from bottom arm:
            nTrialsLeftFly = sum( locsBotTrialT & Dec == -1);
            nTrialsRightFly = sum( locsBotTrialT & Dec == 1);

            % Also read fly's number of trnasitions after the cul-de-sac:
            nTransFly = numTrans.(dataType(10:end-2)).yAfterCul{nF};
        
            % Compute quiver:
            if 1==1 
    
                % Find fly's med x loc.:
                traj_cell_mat = cell2mat( xytyd_bottom );
                xy_bottom_mat = traj_cell_mat(:,1:2);
                uniqueXY = unique( xy_bottom_mat, 'rows' );
                xGoodMedian_allTrajs = mean( ...
                    [min( uniqueXY( abs(uniqueXY(:,2)) < 1, 1 ) ), ...
                    max( uniqueXY( abs(uniqueXY(:,2)) < 1, 1 ) )] );
                % Store xMed required for x-centering:
                xForCentering(nF) = xGoodMedian_allTrajs;
    
                for s = 1:length(decVals)
                    % Read turn-dir. value and colors:
                    decVal = decVals(s);
                    turnName = turnNames{s};
            
                    % Use only fly's bottom trials that ended in specific 
                    % turn-dir and only zero transition trials:
                    if decVal==10
                        goodLocs = (nTransFly==0);
                        xytydBotDecTrans0 = xytyd_bottom( goodLocs );   
                        TurnBotDec0Trans = Dec_bottom( goodLocs );
                    else
                        goodLocs = (nTransFly==0) & (Dec_bottom == decVal);
                        xytydBotDecTrans0 = xytyd_bottom( goodLocs );   
                        TurnBotDec0Trans = Dec_bottom( goodLocs );   
                    end

                    % Compute counts and magnitude-based avg. dir. in Radians 
                        % for each 2D bin:
                    [counts, avgDirThetaRad] = quiverGivenTurns( X, Y, ...
                        xEdges, yEdges, xytydBotDecTrans0, ...
                        xGoodMedian_allTrajs, subsamp );
    
                    % Store fly's magnitude-based avg. dir. and freq.
                    % (normalized counts):
                    binsDirRad.([dataType(10:end-2)]).(turnName).(...
                        walkDir)(:,:,nF) = avgDirThetaRad;
                    binsFreq.([dataType(10:end-2)]).(turnName).(...
                        walkDir)(:,:,nF) = counts / sum(counts(:));
    
                end
    
            end
    
        end


        % Rerun over turn-directions --> compute & plot average quiver for 
        % each [note that this is conditioned on walking direction]
    
        for ssd = 1:length(decVals) 
    
            turnName = turnNames{ssd};
            decCol = decCols{ssd};
    
    
            % COMPUTE AVERAGE QUIVER IN DATATYPE * TURN-DIRECTION: 
    
            % Read the avg. dir's and freq's of all flies in all bins:
            binFreq_all = binsFreq.([dataType(10:end-2)]).(turnName).( ...
                walkDir);
            binsDirRad_all = binsDirRad.([dataType(10:end-2)]).( ...
                turnName).(walkDir);
            
            % Use complex representation for angles in bins:
            binsComplexAng_all =  exp(1i * binsDirRad_all); 
        
            % Compute avg. direction. For each bin, this is the flies' 
            % frequency-based weighted avg. of the angles (in complex 
            % representation):
            binsWeightedMeanAngleComplex = sum( binsComplexAng_all .* ...
                binFreq_all, 3, 'omitnan' ) ./ ...
                sum( binFreq_all, 3, 'omitnan' );
        
            % Convert avg. dir. to angle in radians:
            binsWeightedMeanAngleRad = angle(binsWeightedMeanAngleComplex);
        
            % Compute AVG. magnitude. For each bin, this is the relative 
            % frequency of visits of that bin, averaged over all flies.
            binFreq_allAvg = mean( binFreq_all, 3, 'omitnan' );
            binFreq_allNormAvg = binFreq_allAvg / sum( binFreq_allAvg(:) );
        
        
            % PLOT AVERAGE QUIVER:
        
            currSubPlot = subplot(1,2*length(decVals),2*(ssd-1)+wd);
            
            % Plot heat-map (bivariate PDF):
            if strcmp( overlaidOnHeat, 'yes')
                imagesc( (.5*dyEdges) + yEdges, (.5*dxEdges) + xEdges, ...
                    binFreq_allAvg' / sum(binFreq_allAvg(:)) );
                colorbar;
                set(gca, 'YDir', 'Normal');
                hold on;
            end
    
            % Also plot dataset avg. bounding polygons (if behavioral 
            % data):
            if ~contains( dataType, 'ABM')
                patchline( walkDirPol * polyDataset(:,2), ...
                    polyDataset(:,1), 'linestyle', '-', 'edgecolor', ...
                    boundCol, 'linewidth', 1, 'edgealpha', 0.5 ); 
                hold on;
            end
            
            % Plot the avg. vector field of all flies' trials w/ specific 
            % turn dir. [BOTH: plots either the entire data or left/right
            % overlaid]:
            if strcmp(turnName,'BOTH') && strcmp(bothInterp, 'leftVsRight')
                copyobj( quiveObj.LEFT.(walkDir), currSubPlot );
                copyobj( quiveObj.RIGHT.(walkDir), currSubPlot );
                hold on;
            else
                quiveObj.(turnName).(walkDir) = quiver( ...
                    (.5*dyEdges) + Y, (.5*dxEdges) + X, ...
                    ( binFreq_allAvg / sum(binFreq_allAvg(:)) ) .* ...
                    sin( binsWeightedMeanAngleRad ), ...
                    ( binFreq_allAvg / sum(binFreq_allAvg(:)) ) .* ...
                    cos( binsWeightedMeanAngleRad ), ...
                    decCol, 'LineWidth', 1 ); 
                hold on;
            end
    
            % Define subplot's axes:
            axis image;
            if strcmp(regionName,'culdesac')
                ylim( yLimCul * [-1,1]);
            elseif strcmp(regionName,'intersection')
                ylim(xEdges([1,end]));
            end
            xlim(yEdges([1,end]));
            xlabel('y');
            ggg = gca;
            ggg.YDir = 'reverse';
            ggg.GridAlpha = 0.2;
            grid on;
            if strcmp(walkDir,'down')
                ylabel('x');
            elseif strcmp(walkDir,'up')
                yticklabels([]);
            end
            if strcmp(turnName,'BOTH') && strcmp(bothInterp, 'leftVsRight')
                title([dataType(10:end-2) ' - LEFT vs RIGHT']);
            else
                title([dataType(10:end-2) ' - ' turnName;]);
            end

        end


    end
    

    % Save the figure:
    figure(dataTypeQuiverFig);
    figStartName = ['analyses/figures/flies' dataType(10:end-2) ...
        '_avgQuiverZeroTrans'];
    saveas( dataTypeQuiverFig, [figStartName '.fig']);
    clear meanX_vects_ALL quiveObj;
    
end

