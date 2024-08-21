%% INITIAL PREPROCESSING STAGE OF HUMAN DATA:
%
% Defines edges, divides into turns, omits irrelevant trials and 
% participants.
%
% INPUT FILES:
%   rawdata/rawHumanData/Prolific _[SESSIONID_DATE_TIME].csv
%
% OUTPUT FILES [required for human_mainPreprocess.m]:
%   HumanTurndata_human_visInsp.mat
%   Turndata_humanHuman.mat
%   humanHuman_edgesDef.mat
%
% Note that due to the short experimental durartion for the human datasets 
% (10 mins), here we consider all turns (not only from the bottom arm), by 
% rotating the maze for turns from the right and left arm. Eventually 
% (after rotation), all turns will appear as if they were made from the 
% bottom arm.
%
% Copyright (c) Lior Lebovich, 2024
% lebovich.lior@gmail.com



%% Define folder and data type: 

currCdName = cd;
mainCDName = 'SpatiotempDM';

if ~ endsWith(currCdName,mainCDName)
    cd(mainCDName)
end
addpath('analyses', 'datafiles', 'datafiles/rawdata/rawHumanData')

cd('datafiles/rawdata/rawHumanData')

% Define dataType name:
dataType = 'human';



%% Initial run - select humans (full experiment, not bad trajectories):

% List file names w/ individual human data:
filesInDir = dir;
allFilenames = {filesInDir.name};
filenames = allFilenames( endsWith( allFilenames, '.csv' ) );

% Back to main cdName:
cd('../../../')

% Limits: x,y (for plot) and full experimental duration:
xlims = [-1000,1000];
ylims = [-1100,750];
experimentalTime = 600000;


% Visual inspection: plot trajectories:
allX = [];
allY = [];
lastT = [];
notEndExperiment = [];
sizes = nan(size(filenames,2),2);
for i = 1:size(filenames,2)
    if mod(i,20) == 1
        figure;
        d = 20*floor(i/20);
    end
    % Read file of individual human participant:
    C = readtable(filenames{i});
    subplot(4,5,i-d);
    plot( C.XPos, C.YPos, '-' ); axis image; title(i);
    xlim(xlims); ylim(ylims);
    % Store trajectories of all humans:
    allX = [allX; C.XPos];
    allY = [allY; C.YPos];
    lastT = [lastT; max(C.Time_ms_,[],1,'omitnan')];
    if max(C.Time_ms_,[],1,'omitnan') < experimentalTime
        notEndExperiment = [notEndExperiment, i];
    end
    sizes(i,1) = size(max(C.Time_ms_),1);
    sizes(i,2) = size(max(C.Time_ms_),2);
end

% Plot maze dims:
figMaze = figure; 
plot( allY, allX, '.', 'Color', [.5,.5,.5] ); 
axis image;
set(gca,'XColor', 'none','YColor','none')
figStartName = 'datafiles/rawdata/rawHumanData/humanMazwWithCul';
saveas( figMaze, [figStartName '.fig']);

% Exclude humans with bad trajectories (not reaching cul-de-sac) or not 
% ending experiment:
locsBadTraj = [3, 5, 6, 8, 10, 13, 17, 20, 21, 27, 29, 30, 33, 35, 47];
locExclude = unique([notEndExperiment,locsBadTraj]);
selectedHumanFiles = filenames(~ismember(1:numel(filenames),locExclude));
selectedHumanCD = 'datafiles/rawdata/rawHumanData';

% Save selected human files names and cd names:
save(['datafiles/HumanTurndata_' dataType '_visInsp.mat'], ...
    'selectedHumanFiles', 'selectedHumanCD' );



%% Compute boundaries and rotations (for later dividing into turns):


% Load selected humans:
load(['datafiles/HumanTurndata_' dataType '_visInsp.mat'], ...
    'selectedHumanFiles', 'selectedHumanCD' );


% Define maze boundaries, origin, y-bins and time-bins:

dy = 50;

% Define boundaries locations for trial separation:
leftCorUp = [-633, 518];
leftCorDown = [-720, 364];
rightCorUp = [630, 488]; 
rightCorDown = [721, 340]; 
centUp = [13, 143];
centLeft = [-87, -5];
centRight = [87, -5];
bottomCorLeft = [-82, -750];
bottomCorRight = [92, -750];
xEdgesArm = [-82,88]; 
outCheckY = 350;
% Define y-bins:
yEdges = -1225:dy:1225; 

% Compute maze origin (intersection's triangle center of mass):
intrTriangle = [centUp; centRight; centLeft];
intrPoly = polyshape( intrTriangle(:,1), intrTriangle(:,2) );
[xOrigin, yOrigin] = centroid(intrPoly);
originMaze = [xOrigin, yOrigin];

% Define y-bins centers:
yCenters = (yEdges(1)+0.5*(yEdges(2)-yEdges(1))):...
    (yEdges(2)-yEdges(1)):yEdges(end);

% Define time-bins:
dt = 0.05;
tEdges = -1:dt:1;
tCenters = ((.5*(tEdges(2)-tEdges(1)))+tEdges(1)):...
    (tEdges(2)-tEdges(1)):tEdges(end);

% Compute boundaries:
midLeft = 0.5 * ( leftCorUp + centUp );
midRight = 0.5 * ( rightCorUp + centUp );
slopeLeft = ( leftCorUp(2) - leftCorDown(2) ) / ...
    (leftCorUp(1) - leftCorDown(1) );
slopeRight = ( rightCorUp(2) - rightCorDown(2) ) / ...
    (rightCorUp(1) - rightCorDown(1) );
intsptLeft = midLeft(2) - midLeft(1) * slopeLeft;
intsptRight = midRight(2) - midRight(1) * slopeRight;
funLeft = @(x)( (slopeLeft * x) + intsptLeft );
funRight = @(x)( (slopeRight * x) + intsptRight );


% Define maze's arms and corresp. angles for rotations:

arms = {'bottom', 'right', 'left'};
theoreticalRad = [-pi/2, pi/6, -pi/6];
theoreticalRot = [0, -2*pi/3, 2*pi/3];

% Compute corrections (bottom arm not 90 deg):
% Bottom arm:
slopeBottom_left = ( centLeft(2) - bottomCorLeft(2) ) / ...
    ( centLeft(1) - bottomCorLeft(1) );
slopeBottom_right = ( centRight(2) - bottomCorRight(2) ) / ...
    ( centRight(1) - bottomCorRight(1) );
radBottom = .5 * ( atan(slopeBottom_left) + atan(slopeBottom_right) );
% Right arm:
slopeRight_up = ( rightCorUp(2) - centUp(2) ) / ...
    ( rightCorUp(1) - centUp(1) );
slopeRight_down = ( rightCorDown(2) - centRight(2) ) / ...
    ( rightCorDown(1) - centRight(1) );
radRight = .5 * ( atan(slopeRight_up) + atan(slopeRight_down) );
% Left arm
slopeLeft_up = ( leftCorUp(2) - centUp(2) ) / ...
    ( leftCorUp(1) - centUp(1) );
slopeLeft_down = ( leftCorDown(2) - centLeft(2) ) / ...
    ( leftCorDown(1) - centLeft(1) );
radLeft = .5 * ( atan(slopeLeft_up) + atan(slopeLeft_down) );
% Apply corrections:
corrTheta = theoreticalRad - [radBottom, radRight, radLeft];
thetas = theoreticalRot + corrTheta;



%% Encode turns (bottom, right, left) for each participant:
% (first run the above section):

% The way this is carried out is as follows:
% For a given participant and a given arm from out of 3 from which the 
% participant turns, first rotate all coordinates, so that a turn from that
% arm would translate to turning from the bottom arm (bottom arm not 
% rotated).
% Then, change to NaN all coordinated that aren't either within the bottom
% arm and the internal halves of the left and right arms (after rotation).
% Finally, use these NaNs to define separate turns from the current arm.
% Then, repeat the above for the two other arms.


% Data will be divided into turns and storef in Turncell(s,a,t,k) s.t:
% subject s, 
% arm a, 
% trial t, 
% k:
    % (1) x, 
    % (2) y (normalized), 
    % (3) previous turn (-1/1 left/right), 
    % (4) current turn (-1/1 left/right), 
    % (5) arm,
    % (6) normalized trial-time (-1..1),
    % (7) Time (ms) of trial onset.
Turncell = cell( numel(selectedHumanFiles), numel(arms), 100, 7 );

% Also create time1MsCell, a cell storing the experimental time (ms) of the
% first location in each trial. time1MsCell{s,a} is a vector with the
% first-times of human s and arm a.
time1MsCell = cell( numel(selectedHumanFiles), numel(arms) );

% For storing #trials by arm:
nBotTrials = nan(size(selectedHumanFiles,2),1);
nGoodTrials.bottom = zeros(size(selectedHumanFiles,2),1);
nGoodTrials.right = zeros(size(selectedHumanFiles,2),1);
nGoodTrials.left = zeros(size(selectedHumanFiles,2),1);
yMinGlobVect = nan(size(selectedHumanFiles,2),1);

% Run over humans and arms, rotate, save as turns:
for i = 1:size(selectedHumanFiles,2)
    C = readtable(['datafiles/rawdata/rawHumanData/' ...
        selectedHumanFiles{i}]);
    
    % Compute minimal y-loc (for normalization):
    yMinGlob = min( C.YPos );
    yMinGlobVect(i) = yMinGlob;

    for arm = 1:length(thetas)
        theta = thetas(arm);
        armName = arms{arm};
        
        % Rotation for current arm:
        R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
        rotXY = ( R * ([C.XPos, C.YPos]' - originMaze') + originMaze'  )';
        % Define X,Y after rotation (aligned to bottom, rot. around origin):
        X = rotXY(:,1);
        Y = rotXY(:,2);
        % Also read Time (ms):
        timeMs = C.Time_ms_;

        % Omit irrelevant locations (second halves of unoccupied arms):
        fxLeft = funLeft(X);
        fxRight = funRight(X);
        locOmit = (Y > fxLeft) | (Y > fxRight);
        newX = X;
        newY = Y;
        newTimeMs = timeMs;
        newX(locOmit) = NaN;
        newY(locOmit) = NaN;
        newTimeMs(locOmit) = NaN;
        newX = [NaN; newX; NaN];
        newY = [NaN; newY; NaN];
        newTimeMs = [NaN; newTimeMs; NaN];
    
       % Start encoding trials:
    
       % Find transitions of nan--> not nan:
       naX = isnan(newX);
       delta_naX = diff( naX );
       locTrnsT = find( delta_naX == -1 ); % loc of transition from nan to not-nan
       locTrnsT2 = find( delta_naX == 1 ); % loc of transition from not-nan to nan
       nTrials = min([length(locTrnsT), length(locTrnsT2)]);
       nBotTrials(i) = nTrials;

       % Create vector for storign time-ms of initial loc. in each trial:
       time1MsCell{i,arm} = nan(nTrials,1);

       % Use transitions to code turns from current arm:
       for t = 2:nTrials-1

           trLocStart = locTrnsT(t)+1;
           trLocEnd = locTrnsT2(t);

           % Find x,y,time(ms) in a trial:
           xTrial = newX( trLocStart:trLocEnd );
           yTrial = newY( trLocStart:trLocEnd );
           timeMsTrial = newTimeMs( trLocStart:trLocEnd );

           % Transform y s.t min(yTrialZeroed) = 0:
           yTrialZeroed = abs(yMinGlob) + newY( trLocStart:trLocEnd );
           
           % Apply additional y-correction so that going down/up y<0 / y>0:
           [yMin, yMinLoc] = min( abs(yTrialZeroed) ); 
           yCorrected = [-1*yTrialZeroed(1:yMinLoc); ...
               yTrialZeroed(yMinLoc+1:end)];
           
           % Define relative time (-1..1) w.r.t. minimal location:
           timeVectFromMin = [flip((-1/(yMinLoc-1))*( 1:(yMinLoc-1) ))'; ...
               0; (1/(length(yCorrected)-yMinLoc))*...
               ( 1:(length(yCorrected)-yMinLoc) )'];

           % For later exclusion: omit trials of "change of mind":
           % Omit "changes of mind" when the participant turned from the 
           % "bottom" (after rot.) arm, but didn't arrive at the middle of 
           % the arm to which they turned):
           xTrialOut = xTrial( yCorrected > outCheckY );
           isOutBothSides = ( sum( xTrialOut > xEdgesArm(2) ) >= 1 ) & ...
               ( sum( xTrialOut < xEdgesArm(1) ) >= 1 );
           % Omit "changes of mind" when, before turning to the current
           % "bottom" (after rot.) arm, the participant entered to another 
           % arm, but didn't reach it's middle):
           xTrialIn = xTrial( yCorrected < -outCheckY );
           isInBothSides = ( sum( xTrialIn > xEdgesArm(2) ) >= 1 ) & ...
               ( sum( xTrialIn < xEdgesArm(1) ) >= 1 );

           % Save trial - IFF the trial is of the current arm and omit 
           % changes of mind (see above):

           if ( sum(yTrial<-outCheckY) > 1 ) && ( isOutBothSides ~= 1 ) ...
                   && ( isInBothSides ~= 1 ) && ...
                   ( sum( abs(yCorrected) < outCheckY ) > 1 )
               
               % #trials per arm:
               nGoodTrials.(armName)(i) = nGoodTrials.(armName)(i) + 1;

               % Store x,y coordinates:
               Turncell{i,arm,t,1} = xTrial;
               Turncell{i,arm,t,2} = yCorrected; 
               % Store time (ms) of first trial loc.:
               Turncell{i,arm,t,7} = timeMsTrial(1);
               time1MsCell{i,arm}(t) = timeMsTrial(1);
               % Store previous decision:
               if newX(trLocStart) > centUp(1)
                   Turncell{i,arm,t,3} = -1;
               elseif newX(trLocStart) < centUp(1)
                   Turncell{i,arm,t,3} = 1;
               else
                   Turncell{i,arm,t,3} = 0;
               end
               % Store current decision:
               if newX(trLocEnd) > centUp(1)
                   Turncell{i,arm,t,4} = 1;
               elseif newX(trLocEnd) < centUp(1)
                   Turncell{i,arm,t,4} = -1;
               else
                   Turncell{i,arm,t,4} = 0;
               end
               Turncell{i,arm,t,5} = armName;
               Turncell{i,arm,t,6} = timeVectFromMin;
           end

       end

    end

end


% Plot all corrected trajectories of all humans:
%{
figure;
for a =1:3
    armName = arms{a};
    col = zeros(1,3);
    col(a) = 1;
    for h = 1:size(selectedHumanFiles,2)
        for t = 1:100
            plot( Turncell{h,a,t,1}, Turncell{h,a,t,2}, 'Color', col, ...
                'marker', '.', 'lineStyle', 'none' ); hold on;
        end
    end
end
%} 


% Define good participants (by number of overall trials):
nGoodTrials.all3 = nGoodTrials.bottom + nGoodTrials.left + ...
    nGoodTrials.right;
selectedSubs = find( nGoodTrials.all3 >= 28 );
% Also omit (1) participant with minimal yLoc above cul-de-sac edge:
selectedSubs = selectedSubs([1:7,9:end]);


% Merge each selected human's trials (from different arms) and store data
% of selected humans in a TurncellGood s.t:

% TurncellGood{s}{t,k} corresponds to the data of selected human s, trial t
% (trials are ordered, but trial t-1 and t are not necessarily succeseive 
% in the experiment due to the exclusion of, e.g., "change of mind" 
% trials), and k encodes: 
    % (1) x, 
    % (2) y (normalized), 
    % (3) previous turn (-1/1 left/right), 
    % (4) current turn (-1/1 left/right), 
    % (5) arm,
    % (6) normalized trial-time (-1..1),
    % (7) Time (ms) of trial onset.

% TurncellGood2{s,t,k} corresponds to the data of selected human s, trial t
% (trials are ordered, but trial t-1 and t are not necessarily succeseive 
% in the experiment due to the exclusion of, e.g., "change of mind" 
% trials), and k encodes (1)-(7) as above.
    
TurncellGood = cell( length(selectedSubs), 1 );
TurncellGood2 = cell( length(selectedSubs), max(nGoodTrials.all3), 7 );
yMinGlobVectGood = yMinGlobVect(selectedSubs);
for nS = 1:length(selectedSubs)
    s = selectedSubs(nS);

    % Order trials (of all arms) by ms time:
    arm_trAr_timems = [];
    for arm = 1:numel(arms)
        arm_trAr_timems = [arm_trAr_timems; ...
            arm*ones(size(time1MsCell{s,arm})), ...
            (1:length(time1MsCell{s,arm}))', time1MsCell{s,arm}];
    end
    [B,I] = sort( arm_trAr_timems(:,3) );
    notnanI = I( ~isnan( arm_trAr_timems(I,3) ) );

    % Run over ordered trials and write into cells:
    
    TurncellGood_sub = cell( numel(notnanI), 1 );
    for t = 1:numel(notnanI)
        locTrial = notnanI(t);
        for k = 1:size(Turncell,4)
            % Store trial to TurncellGood2{s,t,k}:
            TurncellGood2{nS,t,k} = Turncell{ s, ...
                arm_trAr_timems(locTrial,1), ...
                arm_trAr_timems(locTrial,2), k };
            % Write trial to subject's cell [for TurncellGood{s}{t,k}]:
            TurncellGood_sub{t,k} = Turncell{ s, ...
                arm_trAr_timems(locTrial,1), ...
                arm_trAr_timems(locTrial,2), k };
        end
    end
    % Store human's trials in TurncellGood{s}{t,k}:
    TurncellGood{nS} = TurncellGood_sub;

end

selectedHumans = 1:length(selectedSubs);


% Save data of selected humans:
save( ['datafiles/Turndata_human' upper(dataType(1)) dataType(2:end) '.mat'], ...
    'TurncellGood', 'TurncellGood2', 'yMinGlobVectGood', ... 
    'selectedHumans' );



%% Define global x and yCorrected edges and by-(selcted) subject "bottom" y edges:

% Plot all selected humans' trajectories: 
figure;
for s = 1:size(TurncellGood2,1)
    col = rand(1,3);
    for t = 1:size(TurncellGood2,2)
        plot( TurncellGood2{s,t,1}, TurncellGood2{s,t,2}, '.', ...
            'color', col ); hold on;
    end
end
title(dataType);

% By visual insepction, define the x and y edges of the "bottom" arm 
% (note that for the above, a portion of the participants that don't 
% reach the edge are not considered):
xBottomEdge = [-82,88]; % [-205,210]
yBottomEdge = [-1027.5, 1025.5];

hold on;
plot( xBottomEdge(1)*[1,1], yBottomEdge, 'k', ...
    xBottomEdge(2)*[1,1], yBottomEdge, 'k', ...
    xBottomEdge, yBottomEdge(1)*[1,1], 'k', ...
    xBottomEdge, yBottomEdge(2)*[1,1], 'k', 'lineWidth', 2 );
axis image;

% Also, compute by-subject "bottom" arm yEdges w/o the intersection:
wInsctnBotArmLenY = nan( size(TurncellGood,1), 1 );
woInsctnBotArmLenY = nan( size(TurncellGood,1), 1 );
for s = 1:size(TurncellGood,1)
    selectedSubX = cell2mat( TurncellGood{s}(:,1) );
    selectedSubY = cell2mat( TurncellGood{s}(:,2) );
    wInsctnBotArmLenY(s) = max( abs( selectedSubY( ...
        selectedSubX >= xBottomEdge(1) & ...
        selectedSubX <= xBottomEdge(2) ) ) );
    woInsctnBotArmLenY(s) = yBottomEdge(2) - ...
        abs( min(yMinGlobVectGood) - yMinGlobVectGood(s) );
end

% Save edges definitions:
save( ['datafiles/human' upper(dataType(1)) dataType(2:end) ...
    '_edgesDef.mat'], 'xBottomEdge', 'yBottomEdge', ...
    'wInsctnBotArmLenY', 'woInsctnBotArmLenY' );


