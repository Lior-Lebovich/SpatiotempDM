%% INITIAL PREPROCESSING STAGE OF FLY DATA AND ABM SIMULATION DATA:
%
% Selects flies, defines edges, computes bounding polygons, omits 
% invalid trials anf agents.
%
% INPUT FILE:
%   rawdata/Turndata_[DATASET_NAME].mat
%
% OUTPUT FILE [required for drosophila_mainPreprocess.m]:
%   Turndata_[DATASET_NAME]_visInsp.mat
%   Turndata_[DATASET_NAME]_edgesDef.mat
%   poly360_Turndata_[DATASET_NAME].mat
%   selectedFlies360DBTurndata_[DATASET_NAME].mat
%   saved_yMinByFly360_[DATASET_NAME].mat



%% Read raw data:

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

cdName = 'SpatiotempDM';
currentFolder = cd;
if ~endsWith(currentFolder,cdName)
    cd(cdName);
end

addpath( 'analyses', 'datafiles', 'datafiles/rawdata', ...
    'analyses/customfuns', '-frozen');



%% Initial run - select flies (not missing, not bad trajectories):

% Omit missing flies:
load([dataType '.mat']);
if exist('Turncell_L')
    Turncell = Turncell_L;
end
maxNTrial=size(Turncell,2);
nFlies = size(Turncell,1);
minX_maxX_minY_maxY = nan(nFlies,4);
for f = 1:nFlies
    all_x = [];
    all_y = [];
    col = rand(1,3);
    for t = 1:maxNTrial
        x = double( Turncell{f,t,1} );
        y = double( Turncell{f,t,2} );
        all_x = [all_x; x];
        all_y = [all_y; y];
    end
    if ~isnan(min( all_x ))
        minX_maxX_minY_maxY(f,1) = min( all_x );
        minX_maxX_minY_maxY(f,2) = max( all_x );
        minX_maxX_minY_maxY(f,3) = min( all_y );
        minX_maxX_minY_maxY(f,4) = max( all_y );
    end
end
selectedFliesLocs = find( ~isnan( sum( minX_maxX_minY_maxY, 2 ) ) );
minX = min( minX_maxX_minY_maxY(:,1) );
maxX = max( minX_maxX_minY_maxY(:,2) );
minY = min( minX_maxX_minY_maxY(:,3) );
maxY = max( minX_maxX_minY_maxY(:,4) );


% Visual inspection - plot fly's trajectories:
maxNTrial = size(Turncell,2);
nselectedFlies = length(selectedFliesLocs);
minX_maxX_minY_maxY = nan(nselectedFlies,4);
for f = 1:nselectedFlies 
    if mod(f,20) == 1
        figure;
    end
    k = 20*floor((f-1)/20);
    nF = selectedFliesLocs(f);
    all_x = [];
    all_y = [];
    col = rand(1,3);
    for t = 1:maxNTrial
        x = double( Turncell{nF,t,1} );
        y = double( Turncell{nF,t,2} );
        all_x = [all_x; x];
        all_y = [all_y; y];
    end
    %{
    subplot(4,5,f-k);
    plot( all_x, all_y, '.', 'color', col ); hold on;
    xlim([minX,maxX]); ylim([minY,maxY]);
    title(num2str(nF));
    %}
    minX_maxX_minY_maxY(f,1) = min( all_x );
    minX_maxX_minY_maxY(f,2) = max( all_x );
    minX_maxX_minY_maxY(f,3) = min( all_y );
    minX_maxX_minY_maxY(f,4) = max( all_y );
end
close all;



%% Visual inspection - exclude flies with bad trajectories + apply corrections:

if strcmp(dataType,'Turndata_ShortYy')
    locExclude = [2,24,60,154,159,165,179,183,233,282,283,284,285,301,308];
elseif strcmp(dataType,'Turndata_LongYy')
    locExclude = [3,206,244,245,268,301,317];
elseif strcmp( dataType, 'Turndata_ShortNorpAYy' )
    locExclude = [1, 4, 6:8, 11:12, 15, 18, ...
        21:26, 28, 30:31, 33, 35, 38:40, ...
        41, 43:45, 48:50, 54:57, 62, 66, 73:74, 84, ...
        85, 87, 90, 95:98, 101:103, ...
        107, 108, 111, 113, 115, 116:118, 120:121, ...
        127:128, 130:131, 133, 135, 138, 141, 143:145, ...
        148:154, 156:157, 159, 161, 163, 164, 167, ...
        169:171, 176:183, 185:186, 188, ...
        189, 191:196, 199:201, 203, 210, 219, 241:243, 246:247, ...
        250:251, 253, 256, 259, 264:267, 271, 274, 280:286, ...
        289:301, 303:305, 307, ...
        309:311, 313:315, 317:323, 325:326, ...
        328:335, 337:338, 343:344, 347:348, ...
        349:351, 353:354, 358, 362:363, 366, ...
        369, 371, 374, 377, 379, 380, 381, 384:385, 388, 390, 392, ...
        393:394, 400, 404, 406, 408, 420:421, 424:427, 429, 433:434, ...
        436, 438:439, 443:445, 447, 452, 454:456, 457:458, 461, ...
        467:474, 476, 479:481, 485, 487:488, 491:493, 498:499, 503, ...
        508:515, 518:519, 523:524, 527, 529:530, 534:536, 539:540, ...
        544:545, 555, 557, 573:575];
elseif strcmp( dataType, 'Turndata_ShortDumbYy' )
    locExclude = [2, 6:9, 11:15, 17:18, 20, 24, ...
        29, 35, 50:55, 60:61, 65, ...
        68, 71, 75, 76, 78:81, 89, ...
        93, 95, 97, 102, 104:106, 108, 113:116, ...
        118, 121:124, 126, 128, 129, 132:136, ...
        140, 145:147, 149, 151:152, 153:154, 157, 159, 161, 163, 166, ...
        168, 169:174, 176:178, 180, 185:189, 241:242, 246, 251, 257:258, ...
        262, 268, 271, 289, 291, 295:301, 310:312, 314, 316:317, ...
        322:323, 338, 343, 350, 354:355, 359:360, 361, 365, 372, 385, ...
        389, 390, 394, 396, 398:402, 405:413, 417:418, 420, 433, 435, ...
        442, 449, 451, 454, 457, 462, 468];
elseif strcmp( dataType, 'Turndata_ShortFoxPYy' )
    locExclude = [8:10, 13, 50:51, 54:56, 58:62, 64:70, 75, 77, 80:83, ...
        84, 85, 87, 90, 93, 95, 98:107, 108:111, 113:114, 116:120, ...
        122:124, 126, 128:129, 130, 132, 134:138, 147, 149, 152, 154, ...
        157, 159, 161:163, 166, 168, 170, 172, 174, 176, 178, 180, ...
        182:183, 185:186, 188:189, 191, 193, 194, 196, 200:204, ...
        207:209, 212:214, 216:223, 225, 227:230, 232:233, 235, 237, ...
        240:242, 244, 246:247, 249:252, 254:257, 259:263, 265:268, ...
        269:276, 278, 280, 282:283, 285, 288:289, 292:294, 296, 300, ...
        302:303, 305, 307:310, 314:318, 320:326, 329, 331:332, 333, ...
        335:336, 341, 343:346, 348:349, 351:353, 354, 356:357, 359, ...
        361, 364, 370, 372, 374:376, 377, 379:381, 383:385, 388, 390, ...
        393:397, 399:401, 403:405, 407, 409, 411:414, 416:418, 420, ...
        423:425, 442, 446, 448, 452, 455, 460, 487:492, 494, 498, 500, ... 
        507, 510, 513, 516, 519, 522:524, 527, 530, 537, 539, 548, 551, ...
        555, 556:557, 561:562, 565, 570:571, 579, 582:584, 588:589, ...
        591:593, 596, 601:602, 603:606, 608, 610:612, 614, 616, ...
        620:621, 624, 626, 630, 633, 635:636, 639, 644:646, 648, ...
        651:654, 657, 659, 665:667, 669, 671, 673, 676, 678:680, ...
        683:685, 687:688, 690, 691:692, 695:697, 699:703, 705, 707, ...
        717, 721, 723:725, 727, 730, 733:734, 736, 738, 740, 746:750, ...
        752:769, 770, 773:774, 777:780, 782, 784:785, 788:789, 790:791, ...
        794, 796:798, 801:805, 807:811, 813, 815, 821:822, 825:827, ...
        829, 831:833, 835, 837:838, 845:846, 849:850, 852, 855, 856, ...
        861, 864, 868:871, 875:876, 879, 881, 884:885, 887:890, 891, ...
        893:895, 898, 900:901, 903:906, 912, 918, 920:922, 925:928, ...
        929:931, 933, 937:938, 941, 944:946, 949, 951, 952:953, 957, 960];
elseif strcmp(dataType,'Turndata_ShortNompCYy')
    locExclude = [9, 10, 15, 64, 71, 75:77, 110, 152:153, 155, 158, ...
        164:167, 170, 176:177, 179, 182, 185, 188:189, 191, 209, 212, ...
        214:215, 218, 225, 227, 230:233, 241:242, 244, 246, 248:249, ...
        253:255, 259:260, 265:266, 271, 282, 299, 303:305, 312, 323, ...
        339, 343:345, 348, 350, 353, 355:356, 361:362, 367:368, ...
        373:374, 377:379, 400, 405, 407, 413, 417, 420, 425:426, 430, ...
        436:437, 440, 445, 448:449, 451, 453, 455, 458:459, 463, ...
        465:466, 468];
elseif strcmp(dataType,'Turndata_ShortABM_NoWF')
    locExclude = [27, 33, 39, 59, 101, 187, 203, 214, 241, 246, 253, ...
        265, 306, 348, 375, 393, 396, 424, 426:427, 476, 489, 491];
elseif strcmp(dataType,'Turndata_ShortABM_WallFollowing')
    locExclude = [8,402];
elseif strcmp(dataType,'Turndata_ShortABM_WF_00100')
    locExclude = [7,9,58,85];
elseif strcmp(dataType,'Turndata_ShortABM_WF_00200')
    locExclude = [3, 21, 25, 31, 79, 80];
elseif strcmp(dataType,'Turndata_ShortABM_WF_00500')
    locExclude = [4, 71, 99];
elseif strcmp(dataType,'Turndata_ShortABM_WF_00700')
    locExclude = [3, 26, 39, 67];
elseif strcmp(dataType,'Turndata_ShortABM_WF_00900')
    locExclude = [52, 78, 97];
elseif strcmp(dataType,'Turndata_ShortABM_brown5')
    locExclude = [];
end
selectedFliesLocs = selectedFliesLocs( ~ismember(selectedFliesLocs,...
    locExclude) );
Turncell = Turncell(selectedFliesLocs,:,:); % consider flies after visual inspection


% Apply correction for REAL and SIMULATED data:

if contains(dataType, 'ABM')
    % Simulated data - correcting w.r.t to ALL flies' min(x), min(Y):
    minX = Inf;
    minY = Inf;
    % Compute minimal x loc. over all flies and trials:
    for f = 1:size(Turncell,1)
        for t = 1:size(Turncell,2)
            minX = min( [minX; Turncell{f,t,1}] );
            minY = min( [minY; Turncell{f,t,2}] );
        end
    end
    % Now implement correction:
    TurncellNew = Turncell;
    for f = 1:size(Turncell,1)
        for tt = 1:size(Turncell,2)
            TurncellNew{f,tt,1} = 1 * ( Turncell{f,tt,1} - minX );
            TurncellNew{f,tt,2} = 1 * ( Turncell{f,tt,2} - minY );
        end
    end  
else
    % Real data - apply corrections for mazes (centering X):
    TurncellNew = Turncell;
    for f = 1:size(Turncell,1) 
        % Compute minimal x loc. over all of the fly's trials:
        minX = Inf;
        for t = 1:size(Turncell,2)
            minX = min( [minX; Turncell{f,t,1}] );
        end
        % Now implement correction:
        for tt = 1:size(Turncell,2)
            TurncellNew{f,tt,1} = Turncell{f,tt,1} - minX;
            TurncellNew{f,tt,2} = Turncell{f,tt,2};
        end
    end
end
% Save corrected turn cell:
Turncell = TurncellNew;


% Save visual inspection output (flies that passed initial visual inspection):
save(['datafiles/' dataType '_visInsp.mat'],'Turncell');



%% Further run - compute edges for ranges:
% These edges will be used later in drosophila_mainPreprocess.m:

load([dataType '_visInsp.mat']);

nFlies = size(Turncell,1);
maxNTrial = size(Turncell,2);
ini_x_all = [];
ini_y_all = [];
figAllInitials = figure;
figAllTrajectories = figure;
for f = 1:nFlies  
    nF = f;
    all_x = [];
    all_y = [];
    ini_x = [];
    ini_y = [];
    col = rand(1,3);
    for t = 2:maxNTrial
        x = double( Turncell{nF,t,1} );
        y = double( Turncell{nF,t,2} );
        all_x = [all_x; x];
        all_y = [all_y; y];
        if sum(~isnan(x))>1 && sum(~isnan(y))>1
            ini_x(t) = x(1);
            ini_y(t) = y(1);
        end
    end
    ini_x_all = [ini_x_all, ini_x];
    ini_y_all = [ini_y_all, ini_y];
    %{
    figure(figAllInitials);
    plot( ini_x, ini_y, 'wo', 'markerFaceColor', col ); grid on; hold on;
    figure(figAllTrajectories);
    plot( all_x, all_y, 'b.' ); grid on; hold on;
    %}
end

%{
figure(figAllInitials);
xlabel('x'); ylabel('y');
title('Initial (fly x trial) x,y locations');

figure(figAllTrajectories);
xlabel('x'); ylabel('y');
title('All flies x,y locations');

figure;
histogram(ini_y_all); xlabel('y'); ylabel('#');
title('initial y locations hist.');
figure;
histogram(ini_x_all); xlabel('x'); ylabel('#');
title('initial x locations hist.');
%}
close all;



%% Define edges (will be used later in drosophila_mainPreprocess.m):

if strcmp(dataType,'Turndata_ShortYy')
    yIniBottomEdges = [10,25];
    xBottomEdge = [20,45]; % entire arm
    xInY_minmax = minmax( ini_x_all( ini_y_all>=yIniBottomEdges(1) & ...
        ini_y_all<=yIniBottomEdges(2) ) );
    xIniBottomEdges = [floor(xInY_minmax(1)), ceil(xInY_minmax(2))]; % only initial location
    yBottomEdge = [0,32];% entire arm
elseif strcmp(dataType,'Turndata_ShortNorpAYy') || ...
        strcmp(dataType,'Turndata_ShortDumbYy') || ...
        strcmp(dataType,'Turndata_ShortFoxPYy') || ...
        strcmp(dataType,'Turndata_ShortNompCYy') 
    yIniBottomEdges = [10,25];
    xBottomEdge = [20,45]; % entire arm
    xInY_minmax = minmax( ini_x_all( ini_y_all>=yIniBottomEdges(1) & ...
        ini_y_all<=yIniBottomEdges(2) ) );
    xIniBottomEdges = [floor(xInY_minmax(1)), ceil(xInY_minmax(2))]; % only initial location
    yBottomEdge = [0,32];% entire arm
elseif strcmp(dataType,'Turndata_LongYy')
    yIniBottomEdges = [20,40];
    xBottomEdge = [40,70]; % entire arm
    xInY_minmax = minmax( ini_x_all( ini_y_all>=yIniBottomEdges(1) & ...
        ini_y_all<=yIniBottomEdges(2) ) );
    xIniBottomEdges = [floor(xInY_minmax(1)), ceil(xInY_minmax(2))]; % only initial location
    yBottomEdge = [0,60];% entire arm
elseif strcmp(dataType,'Turndata_ShortABM_NoWF') || ...
        strcmp(dataType,'Turndata_ShortABM_WallFollowing') || ...
        strcmp(dataType,'Turndata_ShortABM_WF_00100') || ...
        strcmp(dataType,'Turndata_ShortABM_WF_00200') || ...
        strcmp(dataType,'Turndata_ShortABM_WF_00500') || ...
        strcmp(dataType,'Turndata_ShortABM_WF_00700') || ...
        strcmp(dataType,'Turndata_ShortABM_WF_00900') || ...
        strcmp(dataType,'Turndata_ShortABM_brown5')
    yIniBottomEdges = [15,20];
    xBottomEdge = [21,35]; % entire arm
    xInY_minmax = minmax( ini_x_all( ini_y_all>=yIniBottomEdges(1) & ...
        ini_y_all<=yIniBottomEdges(2) ) );
    xIniBottomEdges = [floor(xInY_minmax(1)), ceil(xInY_minmax(2))]; % only initial location
    yBottomEdge = [0,30];% entire arm
end


% Save edges defs:
save(['datafiles/' dataType '_edgesDef.mat'], 'yIniBottomEdges', ...
    'xBottomEdge', 'xIniBottomEdges','yBottomEdge');



%% Second analysis - BEHAVIORAL DATA - count bottom trials & compute boundary
% Consider flies that weren't excluded in the initial run (visual inspection).
% For each of these flies, count the number of bottom trials that fall 
% under the edge bounds that were defined above. If the number of trials
% >=minNumTrials (set to 80) then compute a dynamical boundary and omit 
% trials that fall outside it. 
% Finally, save: relevant flies, final bounding polygon, trials in polygon.

if ~contains(dataType, 'ABM')

load([dataType '_visInsp.mat']); 
load([dataType '_edgesDef.mat']);

nFlies = size(Turncell,1);
maxNTrial = size(Turncell,2);
minNumTrials = 80;
fps = 30;
armNames = {'bottom','right','left'};
isOkayFlies = nan( nFlies, 1 );
xEdgesBot = nan( nFlies, 2 );
xEdgesBotCDS = nan( nFlies, 2 );
xyCenter = nan( nFlies, 2 );
polyAllFlies = cell( nFlies, 1 ); %%%
medianX1015AllFlies = nan( nFlies, 1 ); %%%
nTrialsFlies = nan(nFlies,1);
nTrialsFliesBot = nan(nFlies,1);

figure;

for nF = 1:nFlies
    clear meanX_vects; 
    f = nF;
    
    % Define vectors and cells for individual fly:
    AllTurnsX = cell(maxNTrial,1);
    AllTurnsY = cell(maxNTrial,1);
    yMinFlyVect = nan(maxNTrial,1);
    arm = cell(maxNTrial,1);
    isOkayTrial = nan(maxNTrial,1);
    isOkayTrialBot = nan(maxNTrial,1);
    minMax1025Bot = nan(maxNTrial,2);
    minMax05Bot = nan(maxNTrial,2);
    
    % Collect data from fly (all arms):
    for t = 2:maxNTrial
        x = double( Turncell{f,t,1} ); %x = x(~isnan(x));
        y = double( Turncell{f,t,2} ); %y = y(~isnan(y));
        xy = [x, y];
        xyUnique = unique( xy, 'rows' );
        if ~isempty(x) && (sum(isnan(x)) == 0) && ~isempty(y) && ...
                (sum(isnan(y)) == 0)
            isOkayTrial(t) = 1;
            AllTurnsX{t} = xyUnique(:,1);
            AllTurnsY{t} = xyUnique(:,2);
            % determine current arm:
            if y(1) < yIniBottomEdges(2)
                arm{t} = 'bottom';
                isOkayTrialBot(t) = 1;
                for1025 = [min( x(y>=yIniBottomEdges(1) & ...
                    y<=yIniBottomEdges(2)) ), ...
                    max( x(y>=yIniBottomEdges(1) & y<=yIniBottomEdges(2)) )];
                if isempty(for1025)
                    minMax1025Bot(t,:) = nan(1,2);
                else
                    minMax1025Bot(t,:) = for1025;
                end
            elseif x(1) > xIniBottomEdges(2)
                arm{t} = 'right';
            elseif x(1) < xIniBottomEdges(1)
                arm{t} = 'left';
            end
        end
    end
    
    
    % Operations - only for relevant #trials:
    
    nTrialsFlies(nF) = sum(isOkayTrial==1);
    nTrialsFliesBot(nF) = sum(isOkayTrialBot==1);


    if sum(isOkayTrialBot==1) > minNumTrials
        isOkayFlies(f) = 1;
        xAllCell.all = AllTurnsX( isOkayTrial==1 );
        yAllCell.all = AllTurnsY( isOkayTrial==1 );
        xAll.all = cell2mat( AllTurnsX( isOkayTrial==1 ) );
        yAll.all = cell2mat( AllTurnsY( isOkayTrial==1 ) );
        for armNum = 1:length(armNames)
            armName = armNames{armNum};
            xAll.(armName) = cell2mat( AllTurnsX( (isOkayTrial==1) & ...
                strcmp(arm,armName) ) );
            yAll.(armName) = cell2mat( AllTurnsY( (isOkayTrial==1) & ...
                strcmp(arm,armName) ) );
        end
        
        % run over yCleanEdges --> omit points after pdf(x|y)=0 and keep 
        % only good dots --> save boundary --> save polygon (boundary) 
        % --> run over trials again and omit trials with dots outside of 
        % the polygon --> repeat boundary-omitting procedure until there 
        % is no change in the number of points omitted. Then, save the
        % polygon for later matlab codes.
        yCleanEdges = linspace( min(yAll.all), yBottomEdge(2), 10 ); 
        initialPoly = boundary(xAll.all, yAll.all);

        % initialize loop (based on data)
        xCellP = xAllCell.all;
        yCellP = yAllCell.all;
        xP = xAll.all;
        yP = yAll.all;
        isRepeatLoop = true;
        countloops = 0;
        while isRepeatLoop
            omitLocs = logical( zeros( size(xP) ) ); % 1 for locs to omit
            % run over y --> clean x:
            for nY = 1:length(yCleanEdges)-1
                relevantX = xP( (~omitLocs) & (yP>=yCleanEdges(nY)) & ...
                    (yP<=yCleanEdges(nY+1)) );
                [N,EDGES] = histcounts( relevantX );
                CENTERS = .5*(EDGES(2)-EDGES(1)) + EDGES(1:end-1);
                meanX = mean( relevantX, 'omitnan' );
                rightEDGE = min( [EDGES(end), ...
                    EDGES( 1 + find( (CENTERS>meanX) & (N==0), 1 ) )] );
                leftEDGE = max( [EDGES(1), ...
                    EDGES( find( (CENTERS<meanX) & (N==0), 1, 'last' ) )] );
                omitLocsAdd = ( (yP>=yCleanEdges(nY)) & ...
                    (yP<=yCleanEdges(nY+1)) & ...
                    ( (xP > rightEDGE) | (xP < leftEDGE) ) );
                omitLocs = ( omitLocs | omitLocsAdd );
            end
            % run over x|(0<y<20) --> clean y:
            upperY = 20;
            xMinInLowY= min( xP( (yP>=0) & (yP<=upperY) ) );
            xMaxInLowY = max( xP( (yP>=0) & (yP<=upperY) ) );
            xCleanEdges = linspace( xMinInLowY, xMaxInLowY, 2 );
            for nX = 1:length(xCleanEdges)-1
                relevantY = yP( (~omitLocs) & (xP>=xCleanEdges(nX)) & ...
                    (xP<=xCleanEdges(nX+1)) & ...
                    (yP>=0) & (yP<=upperY) );
                [N,EDGES] = histcounts( relevantY );
                CENTERS = .5*(EDGES(2)-EDGES(1)) + EDGES(1:end-1);
                meanY = mean( relevantY, 'omitnan' );
                rightEDGE = min( [EDGES(end), ...
                    EDGES( 1 + find( (CENTERS>meanY) & (N==0), 1 ) )] );
                leftEDGE = max( [EDGES(1), ...
                    EDGES( find( (CENTERS<meanY) & (N==0), 1, 'last' ) )] );
                omitLocsAdd = ( (xP>=xCleanEdges(nX)) & ...
                    (xP<=xCleanEdges(nX+1)) & ...
                    (yP>=0) & (yP<=upperY) & ...
                    ( (yP > rightEDGE) | (yP < leftEDGE) ) );
                omitLocs = ( omitLocs | omitLocsAdd );
            end
            xAllClean = xP(~omitLocs);
            yAllClean = yP(~omitLocs);
            cleanPoly = boundary( xAllClean, yAllClean );

            % Now redo the poly by omitting all irrelevant trials:
            isTrialInPoly = nan( length(xCellP), 1 );
            for tg = 1:length(xCellP)
                isTrialInPoly(tg) = ( min( ...
                    inpolygon( xCellP{tg}, yCellP{tg}, ...
                    xAllClean(cleanPoly), yAllClean(cleanPoly) ) ) == 1 );
            end
            isRepeatLoop = ( (length(xCellP) - sum(isTrialInPoly)) ~= 0 );
            xP = cell2mat( xCellP( isTrialInPoly==1 ) );
            yP = cell2mat( yCellP( isTrialInPoly==1 ) );
            xCellP = xCellP( isTrialInPoly==1 );
            yCellP = yCellP( isTrialInPoly==1 );
            polyLocs = boundary( xP, yP );
            countloops = countloops + 1;
            [nF, countloops]
        end
        % plot bounding polygon for current fly:
        plot( xP(polyLocs), yP(polyLocs), '-', 'lineWidth', 1, 'color', ...
            rand(1,3) );
        hold on;
        % save bounding polygon:
        polyAllFlies{nF} = [xP(polyLocs), yP(polyLocs)];
        medianX1015AllFlies = median( xP( (yP>=10) & (yP<=15) ) );
        yMinByFly = median( yP );
        
    end

end

xlabel('x'); ylabel('y'); 
title('bounding polygons selected flies');


% Opportunity to omit out-of-center flies that passed:
if strcmp(dataType,'Turndata_ShortYy')
    isOkayFlies( [3,4,5,6,14:20,49] ) = NaN;
elseif strcmp(dataType,'Turndata_LongYy')
    isOkayFlies( [5,11,18,22,23,27] ) = NaN;
elseif strcmp(dataType,'Turndata_ShortNorpAYy')
    isOkayFlies( [12,233] ) = NaN;
elseif strcmp(dataType,'Turndata_ShortDumbYy')
    isOkayFlies( [129,136] ) = NaN;
elseif strcmp(dataType,'Turndata_ShortFoxPYy')
    isOkayFlies( [] ) = NaN;
elseif strcmp(dataType,'Turndata_ShortNompCYy')
    isOkayFlies( [4,56] ) = NaN;
end

% Define polygons of all good flies:
polyAllFlies = polyAllFlies( ~isnan(isOkayFlies) );
% Save polygons:
save(['datafiles/poly360_' dataType '.mat'],'polyAllFlies');

% Define good flies:
selectedFlies = find( ~isnan(isOkayFlies) );
% Save good flies:
save( ['datafiles/selectedFlies360DB' dataType '.mat'], 'selectedFlies' );


end



%% Second analysis - ABM SIMULATIONS - count bottom trials & compute boundary
% Consider flies that weren't excluded in the initial run (visual inspection).
% For each of these flies, count the number of bottom trials that fall 
% under the edge bounds that were defined above. If the number of trials
% >=minNumTrials (set to 80, 1 for Brownian) then compute a dynamical 
% boundary and omit trials that fall outside it. 
% Finally, save: relevant flies, trials in polygon.

if contains(dataType, 'ABM')

load([dataType '_visInsp.mat']); 
load([dataType '_edgesDef.mat']);

nFlies = size(Turncell,1);
maxNTrial = size(Turncell,2);
if strcmp(dataType,'Turndata_ShortABM_brown5')
    minNumTrials = 1;
else
    minNumTrials = 80; 
end
armNames = {'bottom','right','left'};
isOkayFlies = nan( nFlies, 1 );
xEdgesBot = nan( nFlies, 2 );
xEdgesBotCDS = nan( nFlies, 2 );
xyCenter = nan( nFlies, 2 );
nTrialsFlies = nan(nFlies,1);
nTrialsFliesBot = nan(nFlies,1);

figure;

xAllFliesCell_all = cell(0,1);
yAllFliesCell_all = cell(0,1);

for nF = 1:nFlies
    clear meanX_vects; 
    f = nF;
    
    % Define vectors and cells for individual fly:
    AllTurnsX = cell(maxNTrial,1);
    AllTurnsY = cell(maxNTrial,1);
    yMinFlyVect = nan(maxNTrial,1);
    arm = cell(maxNTrial,1);
    isOkayTrial = nan(maxNTrial,1);
    isOkayTrialBot = nan(maxNTrial,1);
    minMax1025Bot = nan(maxNTrial,2);
    minMax05Bot = nan(maxNTrial,2);
    
    % Collect data from fly (all arms):
    for t = 3:maxNTrial
        x = double( Turncell{f,t,1} ); 
        y = double( Turncell{f,t,2} ); 
        xy = [x, y];
        xyUnique = unique( xy, 'rows' );
        if ~isempty(x) && (sum(isnan(x)) == 0) && ...
                ~isempty(y) && (sum(isnan(y)) == 0)
            isOkayTrial(t) = 1;
            AllTurnsX{t} = xyUnique(:,1);
            AllTurnsY{t} = xyUnique(:,2);
            % determine current arm:
            if y(1) < yIniBottomEdges(2)
                arm{t} = 'bottom';
                isOkayTrialBot(t) = 1;
                for1025 = [min( x(y>=yIniBottomEdges(1) & ...
                    y<=yIniBottomEdges(2)) ), ...
                    max( x(y>=yIniBottomEdges(1) & ...
                    y<=yIniBottomEdges(2)) )];
                if isempty(for1025)
                    minMax1025Bot(t,:) = nan(1,2);
                else
                    minMax1025Bot(t,:) = for1025;
                end
            elseif x(1) > xIniBottomEdges(2)
                arm{t} = 'right';
            elseif x(1) < xIniBottomEdges(1)
                arm{t} = 'left';
            end
        end
    end
    
    
    % Use fly if #trials >= min.:
    nTrialsFlies(nF) = sum(isOkayTrial==1);
    nTrialsFliesBot(nF) = sum(isOkayTrialBot==1);
    if sum(isOkayTrialBot==1) > minNumTrials
        isOkayFlies(f) = 1;
        xAllFlyCell_all = AllTurnsX( isOkayTrial==1 );
        yAllFlyCell_all = AllTurnsY( isOkayTrial==1 );
        xAllFly_all = cell2mat( AllTurnsX( isOkayTrial==1 ) );
        yAllFly_all = cell2mat( AllTurnsY( isOkayTrial==1 ) );
        % Store fly's data:
        xAllFliesCell_all = [xAllFliesCell_all; xAllFlyCell_all];
        yAllFliesCell_all = [yAllFliesCell_all; yAllFlyCell_all];
    end

end

xAllFlies_all = cell2mat(xAllFliesCell_all);
yAllFlies_all = cell2mat(yAllFliesCell_all);
        

polyLocs = boundary(xAllFlies_all, yAllFlies_all);

% plot bounding polygon for current fly:
plot( xAllFlies_all(polyLocs), yAllFlies_all(polyLocs), '-', ...
    'lineWidth', 1 );
xlabel('x'); ylabel('y'); 
title('1 bounding polygon - all flies');

% save bounding polygon:
polyAllFlies1 = [xAllFlies_all(polyLocs), yAllFlies_all(polyLocs)];

% Save 1 polygon:
save(['datafiles/poly360_' dataType '.mat'],'polyAllFlies1');

% Define good flies:
selectedFlies = find( ~isnan(isOkayFlies) );
% Save good flies:
save( ['datafiles/selectedFlies360DB' dataType '.mat'], 'selectedFlies' );


end



%% Compute the minimal y location of each fly, after omitting the
% minAbsYtrialToConsider condition:

load([dataType '_visInsp.mat']); 
load([dataType '_edgesDef.mat']);
load(['poly360_' dataType '.mat']); % polyAllFlies
load(['selectedFlies360DB' dataType '.mat']); % selectedFlies

% define positions of Y:
minAbsYtrialToConsider = yIniBottomEdges(1);

maxNTrial = size(Turncell,2);
yMinByFly = nan( length(selectedFlies), 1 );

for nF = 1:length(selectedFlies) 
    clear meanX_vects; 
    f = selectedFlies(nF);
    if strcmp(dataType,'Turndata_ShortABM_NoWF') || ...
        strcmp(dataType,'Turndata_ShortABM_WallFollowing') || ...
        strcmp(dataType,'Turndata_ShortABM_WF_00100') || ...
        strcmp(dataType,'Turndata_ShortABM_WF_00200') || ...
        strcmp(dataType,'Turndata_ShortABM_WF_00500') || ...
        strcmp(dataType,'Turndata_ShortABM_WF_00700') || ...
        strcmp(dataType,'Turndata_ShortABM_WF_00900') || ...
        strcmp(dataType,'Turndata_ShortABM_brown5')
        polyFlySav = polyAllFlies1;
    else
        polyFlySav = polyAllFlies{nF};
    end
    
    % define vectors and cells for individual fly:
    AllTurnsX = cell(maxNTrial,1);
    AllTurnsY = cell(maxNTrial,1);
    yMinFlyVect = nan(maxNTrial,1);
    %
    for t = 2:maxNTrial
        [nF,t];
        x = Turncell{f,t,1}; 
        y = Turncell{f,t,2}; 
        xPrevForCheck = Turncell{f,t-1,1};
        if ~isempty(x) && sum(isnan(x)) == 0 && ~isempty(xPrevForCheck) 
            xCurEnd = x(end);
            % save All turns to compute the x-center of the bottom leg.
            AllTurnsX{t} = x;
            AllTurnsY{t} = y;
            xPrev = Turncell{f,t-1,1}; 
            yPrev = Turncell{f,t-1,2};
            trialInPoly = inpolygon( [double(xPrev); double(x)], ...
                [double(yPrev); double(y)], polyFlySav(:,1), polyFlySav(:,2) );
            % cut: only bottom leg and trial(+prev) in polygon:
            if (y(1) < yIniBottomEdges(2)) && (min(double(y)) <= ...
                    minAbsYtrialToConsider) && (sum(trialInPoly==0)==0)
                
                % cut: from last time (in the trial) passing the bottom leg - going up:
                % also consider prev trial
                % find start and end location in y for mutual information:
                xPrevStart = xPrev(1);
                % considering only a portion of locations of prev vect:
                if xPrevStart > xIniBottomEdges(2)
                    locCut = min( find( xPrev < xPrevStart ) );
                elseif xPrevStart < xIniBottomEdges(1)
                    locCut = min( find( xPrev > xPrevStart ) );
                else
                    locCut = 1e8;
                end
                xCurPrev = [xPrev(locCut:end); x]; 
                yCurPrev = [yPrev(locCut:end); y]; 
                
                yMinFlyVect(t) = min( double(yCurPrev) ); 
            end
        end
    end
    yMinByFly(nF) = min(yMinFlyVect);
end

% save min location:
saved_yMinByFly = yMinByFly;
save(['datafiles/saved_yMinByFly360_' dataType(10:end) '.mat'],...
    'saved_yMinByFly');


