%% EXAMPLE FLY IN SHORT MAZE:
%
% Plots raw and corrected trajectories, pdf(x) and TPI illustration for an 
% example WT fly in the short maze.
%
% INPUT FILE:
%   output_drosophila_main_ShortYy.mat
%
% FUNCTION:
%   TPI.m for computing the TPI.
%
% EXTERNAL FUNCTION:
%   This code uses an external function, patchline.m, Reference: 
%   Brett Shoelson (2023). Patchline 
%   (https://www.mathworks.com/matlabcentral/fileexchange/36953-patchline), 
%   MATLAB Central File Exchange. Retrieved September 18, 2023.
%
% OUTPUT: 
%   saved_exampleFlyLoc.mat: example fly locations.
%
% Figures:
%   Trajectories (raw and corrected):
%       Fig. 1B: exampleFly87__procTrajectories_YX.fig
%       Fig. 1C: exampleFly87__correctedTrajectories_YX.fig
%   Probabability Density Functions of x for example y-bin:
%       Fig. 1D: exampleFly87__yBinEample_xPDF.fig
%       Fig. 1D: exampleFly87__yBinEample_xDeltaPDF.fig
%   Illustration of TPI(Y) computation for example y-bin:
%       Fig. 1E: exampleFly87__yBinEample_pRightGivenMedX.fig
%   Turn-Predictiveness-Index: TPI(Y) and TPI(T) with SEM
%       Fig. 1E: exampleFly87_TPIY.fig
%       Fig. S1A: exampleFly87_TPIT.fig



%% Load short WT dataset:

dataType = 'Turndata_ShortYy';

% Define current directory:
cdName = 'SpatiotempDM';
currentFolder = cd;
if ~endsWith(currentFolder,cdName)
    cd(cdName);
end

% Add path for functions and datafiles:
addpath( 'analyses', 'datafiles', 'analyses/customfuns', '-frozen');

% Load processed WT data:
matFileName = ['output_drosophila_main_' dataType(10:end) '.mat'];
load(matFileName);



%% Select example fly (to modify: change *ONLY* the first line below):

% location of fly in those flies selected for analyses (see: 
% drosophila_selectFlies.m ):
locFly_InSelected = 46; 

% location of fly in those flies selected for analyses (raw data; see: 
% Turndata_ShortYy.mat):
locFly_InAllFlies = selectedFlies(locFly_InSelected); % 87

% define figure(s) start name:
figStartName = ['analyses/figures/exampleFly' num2str(locFly_InAllFlies) '_'];

% Save the fly location [for later example fly figure (quiver plot)]:
save( 'datafiles/saved_exampleFlyLoc.mat', ...
    'locFly_InSelected', 'locFly_InAllFlies' );



%% Read example fly data:

% Read example fly bounding polygon:
polyRawFly = polyAllFlies{locFly_InSelected};

% Also read example fly *PROCESSED* data for bottom trials:

% location of bottom trials:
locsBotTrialT = meanX_vects_ALL.(['f' num2str(locFly_InAllFlies)]...
    ).locsBotTrialT;

% corrected trajectories in bottom trials: 
xytyd = meanX_vects_ALL.(['f' num2str(locFly_InAllFlies)]).xytyd;
xytyd_bottom = xytyd( locsBotTrialT );

% decision vector for bototm trials:
Dec = meanX_vects_ALL.(['f' num2str(locFly_InAllFlies)]).Dec;
Dec_bottom = Dec( locsBotTrialT );



%% Plot example fly *UNPROCESSED* trajectories that lay within example fly 
% bounding polygon (defined in drosophila_selectFlies.m):

figRawXYInPoly = figure;
figRawYXInPoly = figure;
figProcXYInPoly = figure;
figProcYXInPoly = figure;

for t = 2:maxNTrial
    x = Turncell{locFly_InAllFlies,t,1}; 
    y = Turncell{locFly_InAllFlies,t,2};
    xPrev = Turncell{locFly_InAllFlies,t-1,1};  
    yPrev = Turncell{locFly_InAllFlies,t-1,2};
    % Verifing that trial def. trajectories are within bounding polyogn:
    trialInPoly = inpolygon( [double(xPrev); double(x)], ...
        [double(yPrev); double(y)], polyRawFly(:,1), polyRawFly(:,2) );
    xPrevX = double( [xPrev; x] );
    yPrevY = double( [yPrev; y] );
    % Verifing that midlines are too within bounding polyogn (omit "jumps"):
    trialInPoly2 = inpolygon( ...
        .5 * ( xPrevX(1:end-1)+xPrevX(2:end) ), ...
        .5 * ( yPrevY(1:end-1)+yPrevY(2:end) ), ...
        polyFlySav(:,1), polyFlySav(:,2) );
    % Compute processed x,y (y = 1 is the upper bound of the bottom arm,
    % just before the intersection):
    xProcessed = ( double(x) - widCenFly(locFly_InSelected) ) / ...
        (lenBotFly(locFly_InSelected) * ratioLenWoW);
    yProcessed = ( double(y) - minYPolyFly(locFly_InSelected) ) / ...
        (lenBotFly(locFly_InSelected) * ratioLenWoW); 
    if (sum(trialInPoly==0)==0)
        % Plot unproccessed data:
        figure(figRawXYInPoly);
        patchline( double(x)', double( y)', 'linestyle', '-', ...
            'edgecolor', [0,0,0], 'linewidth', 0.5, 'edgealpha', 0.5 ); 
        hold on;
        figure(figRawYXInPoly);
        patchline( double(y)', double( x)', 'linestyle', '-', ...
            'edgecolor', [0,0,0], 'linewidth', 0.5, 'edgealpha', 0.5 ); 
        hold on;
        % Plot proccessed data (y = 1 <-- upper bound of the bottom arm:
        figure(figProcXYInPoly);
        patchline( double(xProcessed)', double( yProcessed)', ...
            'linestyle', '-', 'edgecolor', [0,0,0], 'linewidth', 0.5, ...
            'edgealpha', 0.5 ); 
        hold on;
        figure(figProcYXInPoly);
        patchline( double(yProcessed)', double(xProcessed)', ...
            'linestyle', '-', 'edgecolor', [0,0,0], 'linewidth', 0.5, ...
            'edgealpha', 0.5 ); 
        hold on;
    end
end

% axes + save figure:
figure(figRawXYInPoly); 
axis image; xlabel('x'); ylabel('y'); grid on;
saveas( figRawXYInPoly, [figStartName '_rawTrajectories_XY.fig']);

% axes + save figure:
figure(figRawYXInPoly);
axis image; xlabel('y'); ylabel('x'); grid on;
ggg = gca;
ggg.YDir = 'reverse';
saveas( figRawYXInPoly, [figStartName '_rawTrajectories_YX.fig']);

% axes + save figure:
figure(figProcXYInPoly); 
axis image; xlabel('x'); ylabel('y'); grid on;
saveas( figProcXYInPoly, [figStartName '_procTrajectories_XY.fig']);

% axes + save figure:
figure(figProcYXInPoly);
axis image; xlabel('y'); ylabel('x'); grid on;
ggg = gca;
ggg.YDir = 'reverse';
saveas( figProcYXInPoly, [figStartName '_procTrajectories_YX.fig']);



%% Plot example fly *CORRECTED* trajectories (red/blue: right/left turns):

figCorrectedYX = figure;
colLeftRight = [0,0,1; 1,0,0];

for t = 1:size(xytyd_bottom,1)
    trialData = xytyd_bottom{t};
    trialDecision = Dec_bottom(t);
    xCorrected = trialData(:,1);
    yCorrected = trialData(:,2);
    figure(figCorrectedYX);
    % Down the arm (negative y):
    patchline( yCorrected(yCorrected<0)', xCorrected(yCorrected<0)', ...
        'linestyle','-','edgecolor', ...
        colLeftRight(.5*trialDecision + 1.5, :), 'linewidth', 0.5, ...
        'edgealpha', 0.5 ); 
    hold on;
    % Up the arm (positive y):
    patchline( yCorrected(yCorrected>0)', xCorrected(yCorrected>0)', ...
        'linestyle','-','edgecolor', ...
        colLeftRight(.5*trialDecision + 1.5, :), 'linewidth', 0.5, ...
        'edgealpha', 0.5 ); 
    hold on;
end

% Add corrected bounding polygon:

polyXProcessed = ( polyRawFly(:,1) - widCenFly(locFly_InSelected) ) / ...
    (lenBotFly(locFly_InSelected) * ratioLenWoW);
polyYProcessed = ( polyRawFly(:,2) - minYPolyFly(locFly_InSelected) ) / ...
    (lenBotFly(locFly_InSelected) * ratioLenWoW); 
polyCorrectedFly = [polyXProcessed, polyYProcessed];
% Down the arm (negative y):
plot( -polyCorrectedFly(:,2), polyCorrectedFly(:,1), 'color', 'k', ...
    'lineWidth', 1 );
hold on; 
% Up the arm (positive y):
plot( polyCorrectedFly(:,2), polyCorrectedFly(:,1), 'color', 'k', ...
    'lineWidth', 1 ); 

% axes + save figure:
figure(figCorrectedYX);
axis image; 
xlabel('y'); 
ylabel('x'); 
grid on;
ggg = gca;
ggg.YDir = 'reverse';
saveas( figCorrectedYX, [figStartName '_correctedTrajectories_YX.fig']);

% also plot x,t trajectories:
figCorrectedTX = figure;
colLeftRight = [0,0,1; 1,0,0];

for t = 1:size(xytyd_bottom,1)
    trialData = xytyd_bottom{t};
    trialDecision = Dec_bottom(t);
    xCorrected = trialData(:,1);
    tCorrected = trialData(:,3);
    figure(figCorrectedTX);
    patchline( tCorrected', xCorrected', ...
        'linestyle','-','edgecolor', ...
        colLeftRight(.5*trialDecision + 1.5, :), 'linewidth', 0.5, ...
        'edgealpha', 0.5 ); 
    hold on;
end

% axes + save figure:
figure(figCorrectedTX); 
xlabel('t'); 
ylabel('x'); 
grid on;
ggg = gca;
ggg.YDir = 'reverse';
saveas( figCorrectedTX, [figStartName '_correctedTrajectories_TX.fig']);



%% TPI(Y) and demonstration of its derivation for 1 y-bin:

% Demonstrate delta pdf(x|turn-decision) for example y-bin in example fly:

% Choose sepcific bin for demonstration (can be modified by chaning the 
% line below ONLY):

[~, exBinLoc] = min( abs(yCenters - 0.35) );
exBinEdges = yEdges(exBinLoc:exBinLoc+1);

% Determine x-edges in selected y-bin:
xytyd_bottom_mat = cell2mat( xytyd_bottom );
xInYBin_vect = xytyd_bottom_mat( xytyd_bottom_mat(:,2) >= exBinEdges(1) & ...
        xytyd_bottom_mat(:,2) < exBinEdges(2), 1 );
xInBin_bounds = [1e-1 * floor( 1e1 * min(xInYBin_vect) ), ...
    1e-1 * ceil( max( 1e1 * xInYBin_vect ) )];
xInBin_Edges = linspace( xInBin_bounds(1), xInBin_bounds(2), 18 );
xInBin_Centers = .5 * diff( xInBin_Edges(1:2) ) + xInBin_Edges(1:end-1);

% Compute pdf(x|turn-decision) for x in bottom trials in selected y-bin:
sideNames = {'Left', 'Right'};
sideCols = {'b','r'};
figXPdfByDec = figure;

for dec = unique(Dec_bottom)'
    xytyd_bottom_side_mat = cell2mat( xytyd_bottom( Dec_bottom == dec ) );
    xInYBin_side_vect = xytyd_bottom_side_mat( ...
        xytyd_bottom_side_mat(:,2) >= exBinEdges(1) & ...
        xytyd_bottom_side_mat(:,2) < exBinEdges(2), 1 );
    sideName = sideNames{ .5*dec + 1.5 };
    sideCol = sideCols{ .5*dec + 1.5 };
    
    % Compute pdf(x|turn-decision) for a given side:
    [xPdfInBin_side,~] = histcounts( xInYBin_side_vect, 'binEdges', ...
        xInBin_Edges, 'Normalization', 'pdf' );
    xPdfInBin_sides.(sideName) = xPdfInBin_side;

    % Plot pdf(x|turn-decision) for a given side:
    plot( xInBin_Centers, xPdfInBin_side, sideCol, 'lineWidth', 1 );
    hold on;    
end

% axes + save figure:
figure(figXPdfByDec);
xlim( xInBin_Centers([1,end]) );
xlabel( 'x' );
ylabel('PDF');
title('PDF(x | decision)');
legend( sideNames );
saveas( figXPdfByDec, [figStartName '_yBinEample_xPDF.fig']);

% Plot DELTA pdf(x|turn-decision):
figXDeltaPdf = figure;
figure(figXDeltaPdf);
plot( xInBin_Centers, xPdfInBin_sides.Right - xPdfInBin_sides.Left, ...
    'k', 'lineWidth', 1 );
hold on;
plot( xInBin_bounds, [0, 0], 'k--' );

% axes + save figure:
xlim( xInBin_Centers([1,end]) );
xlabel( 'x' );
ylabel('\DeltaPDF');
saveas( figXDeltaPdf, [figStartName '_yBinEample_xDeltaPDF.fig']);


% Compute TPI(Y) and demonstrate computation for example y-bin in example fly:

% Define (new) y-bin:
yEdges = -1.3:.1:1.3;
yCenters = yEdges(1:end-1) + .5*diff(yEdges(1:2));
infVal = 1.5;

% Compute TPI(Y):
[TPI_y_fly,TPI_y_SEM_fly, pRgR_nLocR_pRgL_nLocL, ~, TPI_y_fly_mpInf, ...
    TPI_y_SEM_fly_mpInf] = TPI( xytyd_bottom, yEdges, Dec_bottom, 'y' );

% Plot computation of TPI(y) for example y-bin:
p_LocR_AND_TurnR = pRgR_nLocR_pRgL_nLocL(1,exBinLoc);
n_LocR = pRgR_nLocR_pRgL_nLocL(2,exBinLoc);
p_LocL_AND_TurnR = pRgR_nLocR_pRgL_nLocL(3,exBinLoc);
n_LocL = pRgR_nLocR_pRgL_nLocL(4,exBinLoc);
pVect = [p_LocL_AND_TurnR, p_LocR_AND_TurnR];
nVect = [n_LocL, n_LocR];
fig_tpiYExampleBin = figure;
figure(fig_tpiYExampleBin);
errorbar( [1,2], pVect, sqrt( pVect .* (1-pVect) ./ nVect ), 'm', ...
    'lineWidth', 1 );

% axes + save figure:
xticks(1:2);
xticklabels( {'Left', 'Right'} );
xlim([.7,2.3]);
ylim([0,1]);
xlabel( 'x median location in y-bin' );
ylabel('P( Right | median(x) )');
saveas( fig_tpiYExampleBin, [figStartName ...
    '_yBinEample_pRightGivenMedX.fig']);


% Plot TPI(y):

figTPIY = figure;

% Add maze regions:
culDeSacY = 0.69/2.03;
regionsBreak = [-1.3, -1, -culDeSacY, culDeSacY, 1, 1.3];
regionPrintNames = {'Pre-arm', 'Down-the-arm', 'Cul-de-sac', ...
    'Up-the-arm', 'Post-arm'};
ax1 = axes();
box(ax1);
for rr = 1:numel(regionsBreak)
    plot(regionsBreak(rr) * [1,1], [-1,1], 'k:' ); 
    hold on;
end

% Add TPI(y):
errorbar( yCenters, TPI_y_fly, TPI_y_SEM_fly, ...
    'CapSize', 0, 'color', 'k', 'linewidth', 1 ); 
hold on;
% Add TPI(+-yInf):
errorbar( infVal*[-1,1], TPI_y_fly_mpInf, TPI_y_SEM_fly_mpInf, ...
    'CapSize', 0, 'color', 'k', 'linewidth', 1, 'Marker', 'o', ...
    'lineStyle', 'none' ); 
hold on;

% Mark example y-bin:
scatter( yCenters(exBinLoc), TPI_y_fly( 1, exBinLoc ), ...
    100, 'filled', 'mo', 'markerFaceAlpha', .5 );

% Axes:
xlim( infVal * [-1,1] ); 
xticks([-infVal, -1:.5:1, infVal]);
xticklabels({['-' 'y' '_{\infty}'], '-1', '-0.5', '0', ...
    '0.5', '1', ['y' '_{\infty}']});
ylim([-1,1]);
yticks(-1:.2:1);  
ggg = gca;
ggg.XGrid = 'off';
ggg.YGrid = 'on';
xlabel('y'); 
ylabel('TPI(y)');
% Add region names:
ax2 = axes('Position', get(ax1,'Position'), 'XAxisLocation', 'top', ...
    'Color', 'none', 'XColor', 'k' );
ax2.YAxis.Visible = 'off';
ax2.XLim = ax1.XLim;
ax2.XTick = regionsBreak(1:end-1) + .5*diff(regionsBreak);
ax2.XTickLabel = regionPrintNames;
ax2.FontSize = 9;

% Save figure:
saveas( figTPIY, [figStartName 'TPIY.fig']);



%% TPI(t) for example fly :

% Define (new) time-bin edges:
tEdges = -1.275:.075:1.275;
tCenters = tEdges(1:end-1) + .5*diff(tEdges(1:2));
infVal = 1.5;

% Compute TPI(t):

% Run tpiT fun. to compute TPI(T):
[TPI_t_fly,TPI_t_SEM_fly, pRgR_nLocR_pRgL_nLocL, ~, TPI_t_fly_mpInf, ...
    TPI_t_SEM_fly_mpInf] = TPI( xytyd_bottom, tEdges, Dec_bottom, 'time' );

% Plot TPI(T):

figTPIT = figure;

% Add maze regions:
regionsBreak = [-1.3, -1, 0, 1, 1.3];
regionPrintNames = {'Pre-arm', 'Down-the-arm & c.d.s<0', ...
    'Up-the-arm & c.d.s>0', 'Post-arm' };
ax1 = axes();
box(ax1);
for rr = 1:numel(regionsBreak)
    plot(regionsBreak(rr) * [1,1], [-1,1], 'k:' ); 
    hold on;
end

% Add TPI(T):
errorbar( tCenters, TPI_t_fly, TPI_t_SEM_fly, 'CapSize', 0, 'color', ...
    'k', 'linewidth', 1 ); 
hold on;
% Add TPI(+-TInf):
errorbar( infVal*[-1,1], TPI_t_fly_mpInf, TPI_t_SEM_fly_mpInf, ...
    'CapSize', 0, 'color', 'k', 'linewidth', 1, 'Marker', 'o', ...
    'lineStyle', 'none' ); 
hold on;

% Axes:
xlim( infVal * [-1,1] ); 
xticks([-infVal, -1:.5:1, infVal]);
xticklabels({['-' 'T' '_{\infty}'], '-1', '-0.5', '0', ...
    '0.5', '1', ['T' '_{\infty}']});
ylim([-1,1]);
yticks(-1:.2:1);  
ggg = gca;
ggg.XGrid = 'off';
ggg.YGrid = 'on';
xlabel('T'); 
ylabel('TPI(T)');
% Add region names:
ax2 = axes('Position', get(ax1,'Position'), 'XAxisLocation', 'top', ...
    'Color', 'none', 'XColor', 'k' );
ax2.YAxis.Visible = 'off';
ax2.XLim = ax1.XLim;
ax2.XTick = regionsBreak(1:end-1) + .5*diff(regionsBreak);
ax2.XTickLabel = regionPrintNames;
ax2.FontSize = 9;

% Save figure:
saveas( figTPIT, [figStartName 'TPIT.fig']);


