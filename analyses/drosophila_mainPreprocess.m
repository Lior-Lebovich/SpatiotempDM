%% MAIN PREPROCESSING STAGE OF FLY DATA AND ABM SIMULATION DATA:
%
% Normalizes trajectories, computes relative time in trial, determines turn 
% decisions, excludes invalid trials.
%
% INPUT FILES [all outputs of drosophila_selectFlies.m]:
%   Turndata_[DATASET_NAME]_visInsp.mat
%   Turndata_[DATASET_NAME]_edgesDef.mat
%   poly360_Turndata_[DATASET_NAME].mat
%   selectedFlies360DBTurndata_[DATASET_NAME].mat
%   saved_yMinByFly360_[DATASET_NAME].mat
%
% EXTERNAL FUNCTION:
%   This code uses an external function, p_poly_dist.m, Reference:
%   Michael Yoshpe (2008). 
%   https://www.mathworks.com/matlabcentral/fileexchange/19398-distance-from-a-point-to-polygon
%   MATLAB Central File Exchange. Retrieved November 8, 2022.
%
% OUTPUT FILE [required for all analyses in the manuscript]:
%   output_drosophila_main_[DATASET_NAME(10:END)].mat
%
% Copyright (c) Lior Lebovich, 2024
% lebovich.lior@gmail.com



%% Define dataset, use relevant dir. and add path to external fucntions/data:

% Select one of the following to manually run over (section by section):
% 'Turndata_ShortYy' 
% 'Turndata_LongYy';
% 'Turndata_ShortNorpAYy'; 
% 'Turndata_ShortDumbYy'; 
% 'Turndata_ShortFoxPYy'; 
% 'Turndata_ShortNompCYy'; 
% 'Turndata_ShortABM_NoWF';
% 'Turndata_ShortABM_WallFollowing';
% 'Turndata_ShortABM_WF_00100';
% 'Turndata_ShortABM_WF_00200';
% 'Turndata_ShortABM_WF_00500';
% 'Turndata_ShortABM_WF_00700';
% 'Turndata_ShortABM_WF_00900';
% 'Turndata_ShortABM_brown5';

dataType = 'Turndata_ShortABM_brown5';

if strcmp(dataType,'Turndata_ShortYy') || ...
    strcmp(dataType,'Turndata_LongYy') || ...
    strcmp(dataType,'Turndata_ShortNorpAYy') || ...
    strcmp(dataType,'Turndata_ShortDumbYy') || ...
    strcmp(dataType,'Turndata_ShortFoxPYy') || ...
    strcmp(dataType,'Turndata_ShortNompCYy')
    real1Sim0 = 1;
elseif strcmp(dataType,'Turndata_ShortABM_NoWF') || ...
    strcmp(dataType,'Turndata_ShortABM_WallFollowing') || ...
    strcmp(dataType,'Turndata_ShortABM_WF_00100') || ...
    strcmp(dataType,'Turndata_ShortABM_WF_00200') || ...
    strcmp(dataType,'Turndata_ShortABM_WF_00500') || ...
    strcmp(dataType,'Turndata_ShortABM_WF_00700') || ...
    strcmp(dataType,'Turndata_ShortABM_WF_00900') || ...
    strcmp(dataType,'Turndata_ShortABM_brown5')
    real1Sim0 = 0;
end

fps = 30;

% Note that the below code does not omit "change of mind" trials. That is, 
% if the fly leaves the bottom arm, then enters the right arm, but then
% switches to the left arm beofre reaching the middle of the right arm -
% then the turn from the bottom arm will be considered as a LEFT TURN.
% (note, that these "changes of mind" seem to be neglected in the human 
% data).

% Define current directory:
cdName = 'SpatiotempDM';
currentFolder = cd;
if ~endsWith(currentFolder,cdName)
    cd(cdName);
end

% Add path for functions and datafiles:
addpath( 'analyses', 'datafiles', 'analyses/customfuns', '-frozen');



%% Load data, apply dataset definitions, plot polygons of included flies:


% Load relevant datasets:
load([dataType '_edgesDef.mat']);
load([dataType '_visInsp']); 
load(['poly360_' dataType '.mat']);
load(['selectedFlies360DB' dataType '.mat']); 
load(['saved_yMinByFly360_' dataType(10:end) '.mat']);
maxNTrial = size(Turncell,2);

% Define ratios of bottom arm length (without/ with intersection):
if strcmp(dataType,'Turndata_ShortYy') || ...
        strcmp(dataType,'Turndata_ShortNorpAYy') || ...
        strcmp(dataType,'Turndata_ShortDumbYy') || ...
        strcmp(dataType,'Turndata_ShortFoxPYy') || ...
        strcmp(dataType,'Turndata_ShortNompCYy')
    ratioLenWoW = 0.8024;
    ratioSize = 1;
    culDeSacY = 0.69/2.03;
elseif strcmp(dataType,'Turndata_LongYy')
    ratioLenWoW = 0.8705;
    ratioSize = 2.53/3.86;
    culDeSacY = 0.68/3.36;
elseif strcmp(dataType,'Turndata_ShortABM_NoWF') || ...
        strcmp(dataType,'Turndata_ShortABM_WallFollowing') || ...
        strcmp(dataType,'Turndata_ShortABM_WF_00100') || ...
        strcmp(dataType,'Turndata_ShortABM_WF_00200') || ...
        strcmp(dataType,'Turndata_ShortABM_WF_00500') || ...
        strcmp(dataType,'Turndata_ShortABM_WF_00700') || ...
        strcmp(dataType,'Turndata_ShortABM_WF_00900') || ...
        strcmp(dataType,'Turndata_ShortABM_brown5')
    ratioLenWoW = 30/34.7;
    ratioSize = 2.53/3.86;
    culDeSacY = 7.5/30;
end

% Redefine xBottomEdge:
if strcmp(dataType,'Turndata_ShortYy')
    xBottomEdge = [21, 43]; 
elseif strcmp(dataType,'Turndata_ShortNorpAYy') || ...
        strcmp(dataType,'Turndata_ShortDumbYy') || ...
        strcmp(dataType,'Turndata_ShortFoxPYy') || ...
    strcmp(dataType,'Turndata_ShortNompCYy')
    xBottomEdge = [23, 40];
elseif strcmp(dataType,'Turndata_LongYy')
    xBottomEdge = [45, 67]; 
elseif strcmp(dataType,'Turndata_ShortABM_NoWF') || ...
        strcmp(dataType,'Turndata_ShortABM_WallFollowing') || ...
        strcmp(dataType,'Turndata_ShortABM_WF_00100') || ...
        strcmp(dataType,'Turndata_ShortABM_WF_00200') || ...
        strcmp(dataType,'Turndata_ShortABM_WF_00500') || ...
        strcmp(dataType,'Turndata_ShortABM_WF_00700') || ...
        strcmp(dataType,'Turndata_ShortABM_WF_00900') || ...
        strcmp(dataType,'Turndata_ShortABM_brown5')
    xBottomEdge = [24, 32.5]; 
end

% Create vectors for storing nTrials, pRight etc.:
yMinByFly = nan( length(selectedFlies), 1 );
pRepeatGood = nan( length(selectedFlies), 1 );
pAll.right = nan( length(selectedFlies), 1 );
pPrevAll.right = nan( length(selectedFlies), 1 );

% Define positions of Y, relative time:
dy = .1;
yEdgeCut = 1.5; 
yEdges = -yEdgeCut:dy:yEdgeCut;
yCenters = .5*diff(yEdges(1:2)) + yEdges(1:end-1);
dt = .075; 
tEdges = -1.5:dt:1.5;
tCenters = ((.5*(tEdges(2)-tEdges(1)))+tEdges(1)):(tEdges(2)-tEdges(1)...
    ):tEdges(end);

% Corrected x,y bins for possible later map:
[yMesh,tMesh] = meshgrid( yEdges, tEdges );
yEdges_new = linspace(yEdges(1),yEdges(end),50);
yCenters_new = (yEdges_new(1)+0.5*(yEdges_new(2)-yEdges_new(1))):...
    (yEdges_new(2)-yEdges_new(1)):yEdges_new(end);
tEdges_new = linspace(tEdges(1),tEdges(end),50);
tCenters_new = (tEdges_new(1)+0.5*(tEdges_new(2)-tEdges_new(1))):...
    (tEdges_new(2)-tEdges_new(1)):tEdges_new(end);
numBinsX = length(yCenters_new); 
xEdges = linspace(-.5, .5, numBinsX+1); 
xCenters = (xEdges(1)+0.5*(xEdges(2)-xEdges(1))):...
    (xEdges(2)-xEdges(1)):xEdges(end);


% BEHAVIORAL DATA: For later y-normalization (bottom arm UP UNTIL 
% intesection--> 1) and x-centering [to learn more, see: 
% yCorrectionExample.m]:
if real1Sim0 == 1
    lenBotFly = nan(size(selectedFlies));
    widCenFly = nan(size(selectedFlies));
    minYPolyFly = nan(size(selectedFlies));
    figure; 
    leg = cell(length(selectedFlies),1);
    for nF = 1:length(selectedFlies)
        f = selectedFlies(nF);
        % Compute fly's bottom arm length WITH intersection:
        pol = polyAllFlies{nF}; % The bounding polygon of the fly' traj's.
        polX = pol(:,1);
        polY = pol(:,2);
        polyYRight = polY;
        polyYRight( polX<=xBottomEdge(2) ) = 0;
        polyYLeft = polY;
        polyYLeft( polX>=xBottomEdge(1) ) = 0;
        [~,locYMaxRight] = max( polyYRight );
        [~,locYMaxLeft] = max( polyYLeft );
        lenPol = length(polyYRight);
        if locYMaxRight<locYMaxLeft
            polyLocs = locYMaxRight:1:locYMaxLeft;
        elseif locYMaxRight>locYMaxLeft
            polyLocs = [locYMaxRight:lenPol, 1:locYMaxLeft];
        end
        minYPolyFly(nF) = min(polY);
        % Compute fly's bottom arm length WITH intersection:
        lenBotFly(nF) = min( polY(polyLocs) ) - minYPolyFly(nF);
        widCenFly(nF) = .5 * max(polX);
        % Eventually, the corrected value would be:
            % y' --> ( y - min(y all trials) ) / (lenBotFly(nF) * ratioLenWoW) 
            % x' --> ( x - widCenFly(nF) ) / (lenBotFly(nF) * ratioLenWoW) 
        % To plot the short and long mazes with real ratio, the correction
        % would be:
            % y'' --> y' * ratioLenWoW / ratioSize
            % x'' --> x' * ratioLenWoW / ratioSize
        leg{nF} = num2str(f); 
        plot( (polX - widCenFly(nF) ) / (lenBotFly(nF) * ratioLenWoW), ...
            ( (polY-minYPolyFly(nF)) / ( lenBotFly(nF) * ratioLenWoW)) ); 
        hold on;
    end
    title(dataType(10:end-2)); 
    axis image; 
    polyAllFlies0 = polyAllFlies;
end

% SIMULATED DATA: y-normalization: Use only 1 polygon:
if real1Sim0 == 0
    % Compute ALL SIM. FLIES bottom arm length WITH intersection:
    figure; 
    pol = polyAllFlies1; % The bounding polygon of the fly' traj's.
    polX = pol(:,1);
    polY = pol(:,2);
    minYPolyAllSims = min(polY);
    % Compute fly's bottom arm length WITH intersection:
    lenBotFlyAllSims = 34.7; % min( polY(polyLocs) ) - minYPolyAllSims;
    widCenFlyAllSims = .5 * max(polX);
    % Plot corrected bounding polygon of all sim. flies:
    plot( (polX - widCenFlyAllSims ) / (lenBotFlyAllSims * ratioLenWoW), ...
        ( (polY-minYPolyAllSims) / ( lenBotFlyAllSims * ratioLenWoW)) ); 
    title(dataType(10:end-2));  
    axis image;
    % Store (identical) values per simulated fly:
    lenBotFly = lenBotFlyAllSims * ones(size(selectedFlies));
    widCenFly = widCenFlyAllSims * ones(size(selectedFlies));
    minYPolyFly = minYPolyAllSims * ones(size(selectedFlies));
    % Also store (identical) bounding polygon for all flies:
    polyAllFlies = repmat( {polyAllFlies1}, length(selectedFlies), 1 );
    polyAllFlies0 = repmat( {[0, 0; 60, 0; 60, 60 ;0, 60]}, ...
        length(selectedFlies), 1 );
end

xGoodMedian_allTrajs = nan(size(selectedFlies));
xGoodMean_allTrajs = nan(size(selectedFlies));



%% Run main preprocessing:

close all;

for nF = 1:length(selectedFlies)
    nF
    clear meanX_vects; 
    f = selectedFlies(nF);
    % Load fly's bounding polygon (irrelevant for SIMULATED flies):
    polyFlySav = polyAllFlies0{nF};
    % Define vectors and cells for individual fly:
    AllTurnsX = cell(maxNTrial,1);
    AllTurnsY = cell(maxNTrial,1);
    meanX_vects.Dec = nan(maxNTrial,1);
    meanX_vects.PrevDec = nan(maxNTrial,1);
    meanX_vects.xytyd = cell(maxNTrial,1);
    meanX_vects.forMedians.AllInRange = nan(maxNTrial,1);
    yMinFlyVect = nan(maxNTrial,1);

    for t = 2:maxNTrial
        x = Turncell{f,t,1}; 
        y = Turncell{f,t,2}; 
        xPrevForCheck = Turncell{f,t-1,1};
        if ~isempty(x) && sum(isnan(x)) == 0 && ~isempty(xPrevForCheck) 
            xCurEnd = x(end);
            % we save All turns to compute the x-center of the bottom leg.
            AllTurnsX{t} = x;
            AllTurnsY{t} = y;
            xPrev = Turncell{f,t-1,1}; 
            yPrev = Turncell{f,t-1,2};
            % Test whether trial lies within bounding polygon:
            trialInPoly = inpolygon( [double(xPrev); double(x)], ...
                [double(yPrev); double(y)], polyFlySav(:,1), polyFlySav(:,2) );
            % Further test whether trial lies within bounding polygon:
            xPrevX = double( [xPrev; x] );
            yPrevY = double( [yPrev; y] );
            trialInPoly2 = inpolygon( ...
                .5 * ( xPrevX(1:end-1)+xPrevX(2:end) ), ...
                .5 * ( yPrevY(1:end-1)+yPrevY(2:end) ), ...
                polyFlySav(:,1), polyFlySav(:,2) );
            
            % Use trial data IFF current turn is from bottom arm AND only
            % if trial is within bounding polygon:
            if (y(1) < yIniBottomEdges(2)) && ...
                    (min(double(y)) <= yIniBottomEdges(1)) && ...
                    (sum(trialInPoly==0)==0) && (sum(trialInPoly2==0)==0)
                
                xPrevStart = xPrev(1);


                % Add OTHER METHOD for brown5: 
                if strcmp(dataType,'Turndata_ShortABM_brown5')
                    % For brown5, we consider almost only the current trial. 
                    % That is, we only take the first location of xPrev,yPrev.
                    % use only x,y (without prev.).
                    % RE-define xPrev, yPrev:
                    xPrev = xPrev(1);
                    yPrev = yPrev(1);
                    xPrevStart = xPrev(1);
                    % Also, for brown5 in the current trial, we consider almost 
                    % only the values until y(standadized)=1.3 (see: yCurPrev 
                    % and yCorrected) below. Thus, for brown5, we here cut the
                    % trajectory data such that we exclude all data following
                    % the first time when the sum. fly passed y=1.3. Because we
                    % don't consider the full xPrev, this will corresp. to
                    % passing y=1.3 during upward motion:
                    yCurRaw = ( double(y) - minYPolyFly(nF) ) / ...
                        (lenBotFly(nF) * ratioLenWoW); 
                    locExBrown = min([find( yCurRaw > 1.3, 1, 'first' ), ...
                        length(yCurRaw)+1]);
                    % To also keep a copy of the real end location, we can use
                    % y=100 (for yInf TPI values).
                    xBrownLastLoc = x(end);
                    yBrownLastLoc = 100;
                    % Re-define x,y
                    x = x(1:locExBrown-1); 
                    y = y(1:locExBrown-1); 
                    % Re-define last x, using the avove re-def.:
                    xCurEnd = x(end);
                end

                
                % Cut trial data such that: from ~middle of previous arm +
                % entering bottom arm + exiting bottom arm + next arm
                % before cul-de-sac:
                if xPrevStart > xIniBottomEdges(2)
                    locCut = 1 + find( xPrev > xBottomEdge(2), 1, 'last' );
                elseif xPrevStart < xIniBottomEdges(1)
                    locCut = 1 + find( xPrev < xBottomEdge(1), 1, 'last' );
                else
                    locCut = 1e8;
                end


                % Normalize the fly's trimmed trial data (this will result
                % in the bottom arm coordinates starting at y=0 and in y=1
                % corresponding to the upper edge of the bottom arm without
                % the intersection; x is also centered and normalized to
                % have xCenter ~around zero to preserve the original 
                % width/height of the maze:
                xCurPrev = ( double([xPrev(locCut:end); x]) - ...
                    widCenFly(nF) ) / (lenBotFly(nF) * ratioLenWoW);
                yCurPrev = ( double([yPrev(locCut:end); y]) - ...
                    minYPolyFly(nF) ) / ...
                    (lenBotFly(nF) * ratioLenWoW); 
                
                % Add OTHER METHOD for brown5: 
                if strcmp(dataType,'Turndata_ShortABM_brown5')
                    % Store only first loc. of previous trial: 
                    xCurPrev = [(xPrevStart-widCenFly(nF) )/...
                        (lenBotFly(nF)*ratioLenWoW); xCurPrev];
                    yCurPrev = [1.3; yCurPrev];
                end

                % Compute distance from bounding polygon:
                dFromPoly = nan( size(xCurPrev) ); 
                for dPol = 1:length(xCurPrev)
                    [distPol,~,~] = p_poly_dist( ...
                        xCurPrev(dPol), yCurPrev(dPol), ...
                        ( polyFlySav(:,1) - widCenFly(nF) ) / ...
                        (lenBotFly(nF) * ratioLenWoW), ...
                        ( polyFlySav(:,2) -  minYPolyFly(nF) ) / ...
                        (lenBotFly(nF) * ratioLenWoW) );
                    dFromPoly(dPol) = distPol;
                end
                    
                % Further operation on y so that y<0 and y>0 correspods to
                % entering and exiting the bottom arm, resp.
                [yMin, yMinLoc] = min( yCurPrev );
                yZeroed = yCurPrev;
                yCorrected = [-1*yZeroed(1:yMinLoc); ...
                    yZeroed(yMinLoc+1:end)];
                yGood = yCorrected;
                yGoodAbsLoc = yCurPrev;
                xGood = xCurPrev;
                dGood = dFromPoly;
                
                % For relative time:
                % Compute values only for WITHIN the bottom arm (no 
                % intersection): -->
                locCutStart = find( yGood <= 1 & yGood >= -1, 1, 'first' );
                locCutEnd = find( yGood <= 1 & yGood >= -1, 1, 'last' );
                % Verify (below cond.) that all trajectories between start to end 
                % fall within the bottom arm. Otherwise exclude from dataset:
                yGood_try = yGood(locCutStart:locCutEnd);
                check_yCurPrev_try = ( sum( (yGood_try < -1) ...
                    | (yGood_try > 1) ) == 0 ); 
                yCurGoodWithin = yGood_try;
                % Add OTHER METHOD for the brown5 dataset: 
                if strcmp(dataType,'Turndata_ShortABM_brown5')
                    % Exclude check_yCurPrev_try for Brownian5h, such that
                    % T=1 will corresp. to the LAST TIME y=1 was reached
                    % during upward motion:
                    check_yCurPrev_try = (1==1);
                end

                
                % Compute relative time:
                    % zero will correspond to the minimal y location (going 
                    % down the maze),
                    % -100% to the upper most y-loc. within the bottom arm 
                    % (without intersection),
                    % +100% to the upper most y-loc. within the bottom arm 
                    % (without intersection). -->

                % --> Find the minimal location (going down):
                % over trajectories WITHIN:
                allLocsWithin = 1:length(yCurGoodWithin);
                locsNegWithin = allLocsWithin( yCurGoodWithin<0 );
                locMinWithin = locsNegWithin(end);
                % over ALL trajectories: 
                allLocs = 1:length(yGood);
                locsNeg = allLocs( yGood<0 );
                locMin = locsNeg(end);
                % Compute relative time:
                timeVectFromMin = [...
                    flip((-1/(locMinWithin-1))*( 1:(locMin-1) ))'; ...
                    0; ...
                    (1/(length(yCurGoodWithin)-locMinWithin))*...
                    ( 1:(length(yGood)-locMin) )'];
                
                % Store trial processed trajectories and decisions if 
                % cond's are met:
                if ~isempty(yGood) && check_yCurPrev_try 
                    % Trajectories and values for median loc's:
                    meanX_vects.xytyd{t} = [xGood, yGood, ...
                        timeVectFromMin, yGoodAbsLoc, -1*(dGood)];
                    meanX_vects.nFrames(t) = length( xGood );
                    meanX_vects.forMedians.AllInRange(t) = mean( xGood );
                    % Determine current decision:
                    if xCurEnd > xIniBottomEdges(2)
                        meanX_vects.Dec(t) = 1;
                    elseif xCurEnd < xIniBottomEdges(1)
                        meanX_vects.Dec(t) = -1;
                    end
                    % Determine previous decision:
                    if xPrevStart > xIniBottomEdges(2)
                        meanX_vects.PrevDec(t) = -1;
                    elseif xPrevStart < xIniBottomEdges(1)
                        meanX_vects.PrevDec(t) = +1;
                    end 
                end 

            end
        end
    end
    

    % Find loc's of relevant trials:
    locsBotTrialT = ( ~isnan( meanX_vects.Dec ) & ...
        ~isnan( meanX_vects.PrevDec ) & ...
        ~isnan( meanX_vects.forMedians.AllInRange ) );

    if strcmp(dataType,'Turndata_ShortABM_brown5')
       locsBotTrialT = ( ~isnan( meanX_vects.Dec ) & ...
           ~isnan( meanX_vects.forMedians.AllInRange ) ); 
    end

    meanX_vects_ALL.(['f' num2str(f)]).xytyd = meanX_vects.xytyd;
    meanX_vects_ALL.(['f' num2str(f)]).locsBotTrialT = locsBotTrialT;
    meanX_vects_ALL.(['f' num2str(f)]).Dec = meanX_vects.Dec;
    meanX_vects_ALL.(['f' num2str(f)]).PrevDec = meanX_vects.PrevDec;

end     

close all;

% Store processed data:
save(['datafiles/output_drosophila_main_' dataType(10:end) '.mat']);


