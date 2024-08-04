% INDIVIDUAL QUIVER DEMONSTRATION:
%
% Computes and plots a quiver plots of an individual example fly, for 
% upward motion within the cul-de-sac, separately for left and right turns, 
% alongside with trajectories from most visited 2D bin.
% For more details about the code and the selection of datatypes, RoIs 
% and flies to plot read the more MORE INFORMATION below.
% 
% INPUT FILES: 
%   output_drosophila_main_[DATASET_NAME(10:END)].mat
%   saved_exampleFlyLoc.mat
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

% MORE INFORMATION:
%
% Here, a density-based quiver plot is computed and plotted for a 
% specific maze region of interest (namely, the cul-de-sac), based on 
% all bottom trials of an individual fly and conditioned on turn 
% direction. Following, the trajectories passing 
% through the most visited location (2D bin) in each turn-dir. are 
% plotted in a separate figure (depending on 'plotExample' value).
%
% The code uses the function quiverGivenTurn.m to compute the 
% magnitude-weighted avgerage direction (in Radians) and the counts 
% (#frames) in each 2D bin in the RoI, accumulated over all chosen 
% trials (for more details - visit quiverGivenTurn.m).
%
% Finally, for each 2D bin, the quiver plots vectors (converted back to
    % Cartesian coordinates) with:
    % (1) Vector's direction corresp. to the magnitude-weighted angle. 
    % (2) Vector's magnitude corresp. to the relative #frames spent in 
        % that  2D bin.
% The above quiver is plotted 3 times: 
    % (1) For the fly's bottom trials that ended in left turns. The quiver 
        % is plotted on top of a bivariate hist. correp. to the rel. 
        % #frames spent in each bin (in fact redundant given the vectors' 
        % magnitude). The bounding polygon plotted above the hist. corrsp. 
        % to all of the fly's loc, in all trials. 
    % (2) The same, for bottom trials that ended in right turns.
    % (3) The left and right quivers, overlayed (without the backround 
        % bivariate histogram).
%
% Datatypes and RoIs: 
% The above is either plotted for a SINGLE FLY in the short maze 
% (default) or for specific/multiple flies in other dataTypes (requires 
% 1 line code modification, under 'plotExample'). If multiple flies 
% are chosen then only the quiver is plotted (without example
% trajectories).
% The default RoI is the cul-de-sac, but could also use the intersection 
% (requires 1 line code modification, under 'regionName').



%% Use relevant directory and add path to external fucntions:

cdName = 'SpatiotempDM';
currentFolder = cd;
if ~endsWith(currentFolder,cdName)
    cd(cdName);
end
addpath( 'analyses', 'datafiles', 'analyses/customfuns', '-frozen');



%% Define quivers': RoI, datatype (+load required data), flies, bins, turns:

% Determine maze region of interest (can be modified to 'intersection');
regionName = 'culdesac';

% Determine flies and datsets for which quiver plot will be displayed:
plotExample = '1FlyShort';
% To view the quiver of selected 1/multiple example flies in a given
% datasets use (resp.) plotExample = 
    % '1FlyShort' / 'multFliesShort'
    % '1FlyLong'/ 'multFliesLong'
    % '1FlyShortNorpA' / 'multFliesShortNorpA'
    % '1FlyShortDumb' / 'Turndata_ShortDumbYy' 
    % '1FlyShortFoxP' / 'multFliesShortFoxP' 

% Define colors and names denoting left and right decisions (for one quiver 
% of both turns dir's use, instead: decVals=10):
decVals = [-1,1]; % Left and right turns
if decVals==10 % plotting one quiver for both left and right turns
    decCols = {'k'};
    decCols2 = {'k'};
    turnNames = {'both'};
else % plotting separate quivers for left/right turns (default):
    decCols = {'k','k'};
    decCols2 = {'b','r'};
    turnNames = {'LEFT', 'RIGHT'};
end
decCols3 = {'b','r'};

% Define subsampling rate:
subsamp = 1; % no subsampling
 
% Derive data type and example flies' locations:
if strcmp( plotExample, '1FlyShort')
    dataType = 'Turndata_ShortYy';
    load('saved_exampleFlyLoc.mat'); % also load exmple fly's loc. 
    nF_vect = locFly_InSelected; 
elseif strcmp( plotExample, 'multFliesShort')
    dataType = 'Turndata_ShortYy';
    nF_vect = [3, 4, 6, 8:11, 14, 18:19, 27:30, 32, 35, 38, 42, 44, ...
        46, 49, 51, 55];
elseif strcmp( plotExample, '1FlyLong')
    dataType = 'Turndata_LongYy';
    nF_vect = 2;
elseif strcmp( plotExample, 'multFliesLong')
    dataType = 'Turndata_LongYy';
    nF_vect = [2:3, 7, 9, 14:15, 27:29, 31, 34:35, 40, 42:43, 45:46];
elseif strcmp( plotExample, '1FlyShortNorpA')
    dataType = 'Turndata_ShortNorpAYy';
    nF_vect = [33, 157]; 
elseif strcmp( plotExample, 'multFliesShortNorpA')
    dataType = 'Turndata_ShortNorpAYy';
    nF_vect = [1:3, 5:8, 10:13, 17:18, 22:27, 30:34, 36, 40, 44, ...
        46, 49, 52, 54, 57, 59, 60:62, 64, 66, 69:75, 78:79, 82:83, ...
        86, 89:91, 93:94, 96:97, 100:110, 113:117, 119:121, 123:125, ...
        127:131, 134:135, 137:138, 142:147, 150:152, 156:157, ...
        159:160, 162, 167:170, 172, 174:177, 180, 183, 185:186, ...
        189:193, 195, 200, 204:206, 210:211, 216, 218:219, 221:223, ...
        226, 228, 233:236, 239:242, 246]; 
elseif strcmp( plotExample, '1FlyShortDumb')
    dataType = 'Turndata_ShortDumbYy';
    nF_vect = 52; % 
elseif strcmp( plotExample, 'multFliesShortDumb')
    dataType = 'Turndata_ShortDumbYy';
    nF_vect = [1:4, 7:8, 12:13, 15, 17:21, 23:27, 30:31, 33, 35, ...
        37:39, 45:50, 52:53, 57, 59, 62:63, 65:68, 70, 74:75, 80, ...
        83 89]; 
elseif strcmp( plotExample, '1FlyShortFoxP')
    dataType = 'Turndata_ShortFoxPYy';
    nF_vect = [45, 28]; 
elseif strcmp( plotExample, 'multFliesShortFoxP')
    dataType = 'Turndata_ShortFoxPYy';
    nF_vect = [3, 6, 12, 16, 18:19, 23, 28, 31:32, 40, 52, 56, 67, ...
        70, 73, 79, 82, 84, 96:98, 100, 101:103, 105:109, 111:114, ...
        116:119, 121:123, 127, 130, 134, 148:149, 153, 155, 162:164, ...
        167:168, 173:174, 176, 178:179, 181:183, 185, 187, 189, ...
        193:194, 196, 198, 200, 208, 215, 217, 223:224, 228, 232, 236]; 
end


% Load dataset:
matFileName = ['output_drosophila_main_' dataType(10:end) '.mat'];
load(matFileName);


% Define upper edge of cul-de-sac based on datatype:
if strcmp(dataType,'Turndata_LongYy')
        culDeSacEdge = 0.68/3.36;
elseif strcmp(dataType,'Turndata_ShortYy') || ...
        strcmp(dataType,'Turndata_ShortNorpAYy') || ...
        strcmp(dataType,'Turndata_ShortDumbYy') || ...
        strcmp(dataType,'Turndata_ShortFoxPYy') 
    culDeSacEdge = 0.69/2.03;
end

% Define 2D (x,y) grid for bins, based on maze size:

% Define X-edges: based on maze size:
binWidthShort = .0162;
binWidth = ratioSize * binWidthShort;
xEdges = -.67*culDeSacEdge:binWidth:.67*culDeSacEdge; 

% Define Y-edges: based on maze size:
if strcmp(regionName,'culdesac')
    yEdges = 0:binWidth:1.05*culDeSacEdge; 
elseif strcmp(regionName,'intersection')
    yEdges = linspace( 1, 1.05/ratioLenWoW, round(16/ratioSize) );
end

% For edges correction:
dxEdges = diff( xEdges(1:2) );
dyEdges = diff( yEdges(1:2) );

% Create grid for 2D (x,y) bins:
[X, Y] = meshgrid( xEdges, yEdges );



%% Compute and plot individual fly's quivers:

for nF = nF_vect
    f = selectedFlies(nF);

    % Read fly's data:
    xytyd = meanX_vects_ALL.(['f' num2str(f)]).xytyd;
    locsBotTrialT = meanX_vects_ALL.(['f' num2str(f)]).locsBotTrialT;
    xytyd_bottom = xytyd( locsBotTrialT );
    Dec = meanX_vects_ALL.(['f' num2str(f)]).Dec;

    % Compute #turn to each side from bottom arm:
    nTrialsLeftFly = sum( locsBotTrialT & Dec == -1);
    nTrialsRightFly = sum( locsBotTrialT & Dec == 1);

    % Compute qiuver IFF, for bottom arm, #leftTurn>=30 and 
    % #rightTurn>=30:
    if nTrialsLeftFly>=30 && nTrialsRightFly>=30

        % Find fly's med x loc.:
        traj_cell_mat = cell2mat( xytyd_bottom );
        xy_bottom_mat = traj_cell_mat(:,1:2);
        uniqueXY = unique( xy_bottom_mat, 'rows' );
        xGoodMedian_allTrajs = mean( ...
            [min( uniqueXY( abs(uniqueXY(:,2)) < 1, 1 ) ), ...
            max( uniqueXY( abs(uniqueXY(:,2)) < 1, 1 ) )] );

        % Create figure for a single fly's quiver plot:
        figSingleFlyQuiver = figure( 'Units', 'Normalized', 'Position', ...
            [0 .3 1 .7] );
    
        % Create figure for a single fly's quiver-corresp. trajectories:
        if startsWith( plotExample, '1')
            figSingleFlyQuiverCorrespTrajs = figure; 
        end
    
        % Run over selected turn sides and plot corrsp. figures:
        for d = 1:length(decVals)
    
            % Read turn-dir. value and colors:
            decVal = decVals(d);
            decCol = decCols{d};
            decCol2 = decCols2{d};
            turnName = turnNames{d};
    
            % Use only fly's: bottom trials that ended in specific 
            % turn-dir.:
            if decVals==10
                xytydBotDec = xytyd( locsBotTrialT );   
                TurnBotDec = Dec( locsBotTrialT );
            else
                xytydBotDec = xytyd( locsBotTrialT & Dec == decVal);   
                TurnBotDec = Dec( locsBotTrialT & Dec == decVal);   
            end
            
            % Compute counts and magnitude-based avg. dir. in Radians 
            % for each 2D bin:
            [counts, avgDirThetaRad] = quiverGivenTurns( X, Y, xEdges, ...
                yEdges, xytydBotDec, xGoodMedian_allTrajs, subsamp );
                    
            
            % QUIVER PLOT: 
    
            figure(figSingleFlyQuiver);

            % Create quiver plot for trials of sepcific turn dir.:

            subplot(1,1+length(decVals),d);
    
            % Plot bivariate PDF:
            imagesc( (.5*dyEdges) + yEdges, (.5*dxEdges) + xEdges, ...
                counts' / sum(counts(:)) );
            colorbar;
            set(gca, 'YDir', 'Normal');
            hold on;
    
            % Plot bounding polygon (based on all of the fly's traj's):
            pol = polyAllFlies{nF};
            plot( ( pol(:,2) - minYPolyFly(nF) ) / ...
                (lenBotFly(nF) * ratioLenWoW), ...
                ( pol(:,1) - widCenFly(nF) ) / ...
                (lenBotFly(nF) * ratioLenWoW) - xGoodMedian_allTrajs, ...
                'color', 'w', 'lineWidth', 1 ); 
            hold on;
            plot( -( pol(:,2) - minYPolyFly(nF) ) / ...
                (lenBotFly(nF) * ratioLenWoW), ...
                ( pol(:,1) - widCenFly(nF) ) / ...
                (lenBotFly(nF) * ratioLenWoW) - xGoodMedian_allTrajs, ...
                'color', 'w', 'lineWidth', 1 ); 
            hold on;
    
            % Plot the vector field for trials w/ specific turn dir.:
            quiver( (.5*dyEdges) + Y, (.5*dxEdges) + X, ...
                ( counts / sum(counts(:)) ) .* sin( avgDirThetaRad ), ...
                ( counts / sum(counts(:)) ) .* cos( avgDirThetaRad ), ...
                decCol, 'LineWidth', 1 ); 
            hold on;
    
            % Define subplot's axes:
            axis image;
            ylim(xEdges([1,end]));
            xlim(yEdges([1,end]));
            xlabel('y');
            ylabel('x');
            title([turnName ' turns, fly ' num2str(f) ', ' ...
                dataType(10:end-2)]);
            ggg = gca;
            ggg.YDir = 'reverse';
            ggg.Color = 'k';
            ggg.XColor = 'k'; 
            ggg.YColor = 'k'; 
            ggg.GridColor = 'y';
            ggg.GridAlpha = 0.2;
            grid on;
            
            % Create quiver plots for trials of left and left turn. dir's,
            % overlayed:
            
            subplot(1,1+length(decVals),1+length(decVals));
    
            % Plot bounding polygon (based on all pf the fly's traj's):
            plot( ( pol(:,2) - minYPolyFly(nF) ) / ...
                (lenBotFly(nF) * ratioLenWoW), ...
                ( pol(:,1) - widCenFly(nF) ) / ...
                (lenBotFly(nF) * ratioLenWoW) - xGoodMedian_allTrajs, ...
                'color', 'k', 'lineWidth', 1 ); 
            hold on;
            plot( -( pol(:,2) - minYPolyFly(nF) ) / ...
                (lenBotFly(nF) * ratioLenWoW), ...
                ( pol(:,1) - widCenFly(nF) ) / ...
                (lenBotFly(nF) * ratioLenWoW) - xGoodMedian_allTrajs, ...
                'color', 'k', 'lineWidth', 1 ); 
            hold on;
            
            % Plot the vector field for trial w/ specific tirn dir.:
            quiver( (.5*dyEdges) + Y, (.5*dxEdges) + X,...
                ( counts / sum(counts(:)) ) .* sin( avgDirThetaRad ),...
                ( counts / sum(counts(:)) ) .* cos( avgDirThetaRad ), ...
                decCol2, 'LineWidth', 1 ); 
            hold on;
    
            % Define subplot's axes:
            axis image;
            ylim(xEdges([1,end]));
            xlim(yEdges([1,end]));
            xlabel('y');
            ylabel('x');
            title(['fly ' num2str(f) ' (selected ' num2str(nF) '), ' ...
                dataType(10:end-2)]);
            ggg = gca;
            ggg.YDir = 'reverse';
            ggg.GridAlpha = 0.2;
            grid on;
            
            % Store counts vector for plotting example traj. (see below):
            savedCounts.(turnName) = counts;
            
    
            % QUIVER-CORRESP. TRAJECTORIES: 

            % Only for single fly examples - also mark the most visited
            % loc. and plot trajectories pathing it in a separate figure:

            if startsWith( plotExample, '1')

                % Find the loc's of the most visited 2D bin (bottom arm, 
                % current turnning dir.:
                [rowMax,colMax] = find( counts == max(counts,[],'all') ); 
                
                % Find bin edges for most visited 2D bin:
                xMaxEdges = xEdges( colMax:(colMax+1) );
                yMaxEdges = yEdges( rowMax:(rowMax+1) );
        
                % Def. square around most visited bin:
                mostVisitY = [yMaxEdges, flip(yMaxEdges)];
                mostVisitX = [xMaxEdges(1)*[1,1], xMaxEdges(2)*[1,1]];

                % Mark the most visited bin with a square (over quiver 
                % plot):

                figure(figSingleFlyQuiver);
                subplot(1,1+length(decVals),d);
                plot( mostVisitY([1:end,1]), mostVisitX([1:end,1]), ...
                    'm-', 'LineWidth', 1 );
    
                % Find fly's traj's for ALL turn-dir's that pass through 
                % the most visited 2D bin (found above). For these traj's 
                % - plot the traj. starting from when that bin was visited 
                % for the LAST time in the trial:

                figure(figSingleFlyQuiverCorrespTrajs);
                subplot(1,length(decVals),d);
                
                % Run over ALL turns dir's:
                xytydBotBothDec = xytyd( locsBotTrialT );
                TurnBotBothDec = Dec( locsBotTrialT );
                
                for tt = 1:numel(xytydBotBothDec)
                    % Read the x, y coordinates for the current trial:
                    xy_raw = xytydBotBothDec{tt}( :, 1:2 ); 
                    xy = [xy_raw(:,1)-xGoodMedian_allTrajs, xy_raw(:,2)];
                    % Find loc. for last visit of most visited 2D bin 
                    % (found above):
                    locVisLast = find( xy(:,1) >= xMaxEdges(1) & ...
                        xy(:,1) <= xMaxEdges(2) & ...
                        xy(:,2) >= yMaxEdges(1) & ...
                        xy(:,2) <= yMaxEdges(2), 1, 'last' );
                    % If bin visited - Plot traj. from when last visit:
                    if ~isempty(locVisLast)
                        patchline( xy(locVisLast:end,2)', ...
                            xy(locVisLast:end,1)', 'linestyle', ...
                            '-', 'edgecolor', ...
                            decCols3{ .5*TurnBotBothDec(tt) + 1.5}, ...
                            'linewidth', 0.5, 'edgealpha', 0.5 ); 
                        hold on;
                    end
                end
        
                % Also plot bounding polygon (based on all of the fly's 
                % traj's):
                plot( ( pol(:,2) - minYPolyFly(nF) ) / ...
                    (lenBotFly(nF) * ratioLenWoW), ...
                    ( pol(:,1) - widCenFly(nF) ) / ...
                    (lenBotFly(nF) * ratioLenWoW) - ...
                    xGoodMedian_allTrajs, 'color', 'k', 'lineWidth', 1 ); 
                hold on;
            
                % Mark the most visited bin with a square (over traj's 
                % plot):
                plot( mostVisitY([1:end,1]), mostVisitX([1:end,1]), ...
                    'm-', 'LineWidth', 1 );

                % Define subplot's axes:
                axis image;
                xlabel('y');
                ylabel('x');
                title([turnName ' turns from fly ' num2str(f) ...
                    ' most visited bin']);
                ggg = gca;
                ggg.YDir = 'reverse';
                ggg.GridAlpha = 0.2;
                grid on; 
            end

        end

    end

end


