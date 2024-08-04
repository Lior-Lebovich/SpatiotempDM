%% MAIN PREPROCESSING STAGE OF HUMAN DATA:
%
% Normalizes trajectories, computes relative time in trial, determines turn 
% decisions, excludes invalid trials.
% Note that the pre-processed human data requires less additional
% corrections post human_selectHumans.m (compared to flies' datasets). 
%
% INPUT FILES [all outputs of human_selectHumans.m]:
%   HumanTurndata_human_visInsp.mat
%   Turndata_humanHuman.mat
%   humanHuman_edgesDef.mat
%
% OUTPUT FILE [required for all analyses in the manuscript]:
%   output_drosophila_main_HumanYy.mat

    

%% Use relevant directory and add path to external fucntions and datafiles:

cdName = 'SpatiotempDM';
currentFolder = cd;
if ~endsWith(currentFolder,cdName)
    cd(cdName);
end
addpath( 'analyses', 'datafiles', 'analyses/customfuns', '-frozen');



%% Load dataset, bins and locations:

% Define dataset:
dataType = 'Turndata_humanHuman'; 
fps = 20;
culDeSacY = 0.3;

% Load and re-define edges:
load([dataType(10:end) '_edgesDef.mat']); 
% Redefine x-edges:
xIniBottomEdges = [-82,88]; 
xBottomEdge = [-205,210];

% Load turn data:
load([dataType '.mat']); 
selectedFlies = selectedHumans;
saved_yMinByFly = yMinGlobVectGood;
Turncell = TurncellGood2;
% TurncellGood2{s,t,k} corresponds to the data of selected human s, trial t
% (trials are ordered, but trial t-1 and t are not necessarily succeseive 
% in the experiment due to the exclusion of, e.g., "change of mind" 
% trials), and k (1)-(7) encodes:
    % (1) x, 
    % (2) y (normalized), 
    % (3) previous turn (-1/1 left/right), 
    % (4) current turn (-1/1 left/right), 
    % (5) arm,
    % (6) normalized trial-time (-1..1),
    % (7) Time (ms) of trial onset.
maxNTrial = size(Turncell,2);
    
% Create vectors for storing nTrials, pRight etc.:
pAll.right = nan( length(selectedFlies), 1 );
pPrevAll.right = nan( length(selectedFlies), 1 );
yMinByFly = nan( length(selectedFlies), 1 );
nTrialVect = nan( length(selectedFlies), 1);
pRepeatGood = nan( length(selectedFlies), 1 );

% Define positions of Y, relative time:
dy = .1;
yEdgeCut = 1.5; 
yEdges = -yEdgeCut:dy:yEdgeCut;
yCenters = .5*diff(yEdges(1:2)) + yEdges(1:end-1);
dt = .075; 
tEdges = -1.5:dt:1.5;
tCenters = ((.5*(tEdges(2)-tEdges(1)))+tEdges(1)):(...
    tEdges(2)-tEdges(1)):tEdges(end);

% Corrected x,y bins for possible later map:
yEdges_new = linspace(yEdges(1),yEdges(end),ceil(length(yEdges)/1)); 
yCenters_new = (yEdges_new(1)+0.5*(yEdges_new(2)-yEdges_new(1))):...
    (yEdges_new(2)-yEdges_new(1)):yEdges_new(end);
tEdges_new = linspace(tEdges(1),tEdges(end),ceil(length(yEdges)/1));
tCenters_new = (tEdges_new(1)+0.5*(tEdges_new(2)-tEdges_new(1))):...
    (tEdges_new(2)-tEdges_new(1)):tEdges_new(end);
numBinsX = length(yCenters_new); 
xEdges = linspace(-160/yBottomEdge(2), 160/yBottomEdge(2), numBinsX+1); 
xCenters = (xEdges(1)+0.5*(xEdges(2)-xEdges(1))):...
    (xEdges(2)-xEdges(1)):xEdges(end);

% For later distance calculations:
yEdgesForBotDist = 1e2*floor(1e-2*min(woInsctnBotArmLenY)) * [-1,1];
allX = [];
allY = [];
for f = 1:size(TurncellGood2,1)
    for t = 1:size(TurncellGood2,2)
        x = TurncellGood2{f,t,1};
        y = TurncellGood2{f,t,2};
        allX = [allX; x];
        allY = [allY; y];
    end
end
xInBot = allX( allY >= yEdgesForBotDist(1) & ...
    allY <= yEdgesForBotDist(2) );
minMaxXBot = [min(xInBot), max(xInBot)];
clear f t allX allY x y xInBot;



%% Run main preprocessing:

for nF = 1:length(selectedFlies) 
    clear meanX_vects; 
    f = selectedFlies(nF);
    % For y-normalization (bottom arm before intesection--> 1)
    flyWoLen = woInsctnBotArmLenY(nF); 
    % define vectors and cells for individual fly:
    meanX_vects.Dec = nan(maxNTrial,1);
    meanX_vects.PrevDec = nan(maxNTrial,1);
    meanX_vects.xytyd = cell(maxNTrial,1);
    meanX_vects.nFrames = nan(maxNTrial,1);
    meanX_vects.forMedians.AllInRange = nan(maxNTrial,1);
    yMinFlyVect = nan(maxNTrial,1);
    median_parts = nan( length(yEdges)-1, 1 );
    
    for t = 1:maxNTrial
        decCur = Turncell{f,t,4};
        decPrev = Turncell{f,t,3};
        xCurPrev = Turncell{f,t,1};
        % Normalize y [x will be centered and normalized later on, see: 
         % xGood]: 
        yCurPrev = Turncell{f,t,2} / flyWoLen;
        
        if ~isempty(xCurPrev) && sum(isnan(xCurPrev)) == 0 
            
            
            % cut: only bottom leg (both by x and y lims):

            % Compute values only for within the bottom arm: -->
            
            % --> where to start and where to end:
            locCutStart = find( yCurPrev <= 1 & ...
                yCurPrev >= -1, 1, 'first' );
            locCutEnd = find( yCurPrev <= 1 & ...
                yCurPrev >= -1, 1, 'last' );
            % Verify that all trajectories between start to end fall within
            % the bottom arm. Otherwise exclude:
            xCurPrev_try = xCurPrev(locCutStart:locCutEnd);
            yCurPrev_try = yCurPrev(locCutStart:locCutEnd);
            check_xCurPrev_try = ( sum( (xCurPrev_try < xBottomEdge(1)) ...
                | (xCurPrev_try > xBottomEdge(2)) ) == 0 ); 
            check_yCurPrev_try = ( sum( (yCurPrev_try < -1) ...
                | (yCurPrev_try > 1) ) == 0 ); 

            % Also verify that there aren't multiple enterances/exits 
            % to cul-de-sac within a given trial:
            locCutStartCUL = find( yCurPrev <= culDeSacY & ...
                yCurPrev >= -culDeSacY, 1, 'first' );
            locCutEndCUL = find( yCurPrev <= culDeSacY & ...
                yCurPrev >= -culDeSacY, 1, 'last' );
            yCurPrev_tryCUL = yCurPrev(locCutStartCUL:locCutEndCUL);
            check_yCurPrev_tryCUL = ( sum( (yCurPrev_tryCUL < -culDeSacY) ...
                | (yCurPrev_tryCUL > culDeSacY) ) == 0 ); 
            

            % Run trial IFF all trajectories that should be are within 
            % bottom arm:
            if check_xCurPrev_try && check_yCurPrev_try && ...
                    check_yCurPrev_tryCUL

                yMinFlyVect(t) = min( abs( double(yCurPrev) ) ); 
                
                % The trajectories within the bottom arm (without the 
                % intersection):
                xCurPrevWithin = xCurPrev_try; 
                yCurPrevWithin = yCurPrev_try;

                % Compute distance from edges:
                dFromPoly = min( abs(xCurPrev - minMaxXBot), [], 2);
                
                % Store x,y,d for the entire trial:
                yGood = yCurPrev;
                yGoodAbsLoc = abs( yGood );
                % Center and normalize X:
                xGood = ( xCurPrev - mean(minMaxXBot) ) / flyWoLen;
                dGood = dFromPoly;

                % Compute relative time:
                    % zero will correspond to the minimal y location (going 
                    % down the maze),
                    % -100% to the upper most y-loc. within the bottom arm 
                    % (without intersection),
                    % +100% to the upper most y-loc. within the bottom arm 
                    % (without intersection). -->
                
                % --> Find the minimal location (going down):
                % over trajectories WITHIN: 
                allLocsWithin = 1:length(yCurPrevWithin);
                locsNegWithin = allLocsWithin( yCurPrevWithin<0 );
                locMinWithin = locsNegWithin(end);
                % over ALL trajectories: 
                allLocs = 1:length(yGood);
                locsNeg = allLocs( yGood<0 );
                locMin = locsNeg(end);
                % Compute relative time:
                timeVectFromMin = [...
                    flip((-1/(locMinWithin-1))*( 1:(locMin-1) ))'; ...
                    0; ...
                    (1/(length(yCurPrevWithin)-locMinWithin))*...
                    ( 1:(length(yGood)-locMin) )'];

                % Store locations:
                meanX_vects.xytyd{t} = [xGood, yGood, ...
                    timeVectFromMin, yGoodAbsLoc, -1*(dGood)];
                meanX_vects.nFrames(t) = length( xGood );
                meanX_vects.forMedians.AllInRange(t) = mean( xGood );
        
                % determine current decision:
                meanX_vects.Dec(t) = Turncell{f,t,4};
                
                % determine previous decision:
                meanX_vects.PrevDec(t) = Turncell{f,t,3};

            end
        end
    end
    yMinByFly(nF) = min(yMinFlyVect);
    
    
    % choose relevant trials:
    locsBotTrialT = ( ~isnan( meanX_vects.Dec ) & ...
        ~isnan( meanX_vects.PrevDec ) & ...
        ~isnan( meanX_vects.forMedians.AllInRange ) );

    if 1==1 

        % compute the medians of x locations in over what's def. trial:
        meanXs_good = meanX_vects.forMedians.AllInRange( locsBotTrialT );
        Dec_good = meanX_vects.Dec( locsBotTrialT );
        pAll.right(nF) = mean( Dec_good == 1);
        prevDec_good = meanX_vects.PrevDec( locsBotTrialT );
        pPrevAll.right(nF) = mean( prevDec_good == 1);
        pRepeatGood(nF) = mean( Dec_good == prevDec_good );
        xytyd_good = meanX_vects.xytyd( locsBotTrialT );
        nFrames_good = meanX_vects.nFrames( locsBotTrialT );
        nGoodTrials = sum(locsBotTrialT); nTrialVect(nF) = nGoodTrials;
        xy_good_trial = cell(nGoodTrials,1);
        xt_good_trial = cell(nGoodTrials,1);
        yd_good_trial = cell(nGoodTrials,1);
        
        for tk = 1:nGoodTrials
            xytyd_good_trial = xytyd_good{tk};
            xy_good_trial{tk} = xytyd_good_trial( : ,1:2 );
            xt_good_trial{tk} = xytyd_good_trial( : ,[1,3] );
            yd_good_trial{tk} = xytyd_good_trial( : ,[2,5] );
        end
        uniqueXY = unique( cell2mat( xy_good_trial ), 'rows' );
        uniqueXT = unique( cell2mat( xt_good_trial ), 'rows' );
        
        meanX_vects_ALL.(['f' num2str(f)]).xytyd = meanX_vects.xytyd;
        meanX_vects_ALL.(['f' num2str(f)]).locsBotTrialT = locsBotTrialT;
        meanX_vects_ALL.(['f' num2str(f)]).Dec = meanX_vects.Dec;

    end
end     

close all;
save('datafiles/output_drosophila_main_HumanYy.mat');

clearvars -except cdName fps dataTypes d;


