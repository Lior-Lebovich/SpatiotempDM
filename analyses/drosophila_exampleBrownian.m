%% EXAMPLE BROWNIAN AGENT IN SHORT MAZE:
%
% Plots raw and corrected trajectories for an example ABM agent.
%
% INPUT FILE:
%   output_drosophila_main_ShortABM_brown5.mat
%
% EXTERNAL FUNCTION:
%   This code uses an external function, patchline.m, Reference: 
%   Brett Shoelson (2023). Patchline 
%   (https://www.mathworks.com/matlabcentral/fileexchange/36953-patchline), 
%   MATLAB Central File Exchange. Retrieved September 18, 2023.
%
% OUTPUT: 
%   saved_exampleBrownianLoc.mat: example agent locations.
%
% Figures:
%   Trajectories (raw untrimmed and and corrected+trimmed):
%       Fig. S14A: exampleBrownian1__procTrajectories_YX.fig
%       Fig. S14B: exampleBrownian1__correctedTrajectories_YX.fig



%% Use relevant dir., add path to fucntions/data, load Brownian dataset:

% Define current directory:
cdName = 'SpatiotempDM';
currentFolder = cd;
if ~endsWith(currentFolder,cdName)
    cd(cdName);
end
addpath( 'analyses', 'datafiles', 'analyses/customfuns', '-frozen');

% Load datatype:
dataType = 'Turndata_ShortABM_brown5';
matFileName = ['output_drosophila_main_' dataType(10:end) '.mat'];
load(matFileName);



%% Select example fly (to modify: change *ONLY* the first line below):

% location of fly in those flies selected for analyses (see: 
% drosophila_selectFlies.m ):
locFly_InSelected = 1; 

% location of fly in those flies selected for analyses (raw data; see: 
% Turndata_ShortYy.mat):
locFly_InAllFlies = selectedFlies(locFly_InSelected); % 87

% define figure(s) start name:
figStartName = ['analyses/figures/exampleBrownian' ...
    num2str(locFly_InAllFlies) '_'];

% Save the fly location [for later example fly figure (quiver plot)]:
save( 'datafiles/saved_exampleBrownianLoc.mat', ...
    'locFly_InSelected', 'locFly_InAllFlies' );


%% Read example fly data:

% Read example fly bounding polygon:
polyRawFly = polyAllFlies{locFly_InSelected};

% Also read example fly *PROCESSED* data for bottom trials:

% location of bottom trials:
locsBotTrialT = meanX_vects_ALL.(['f' num2str(locFly_InAllFlies)]).locsBotTrialT;

% corrected trajectories in bottom trials: 
xytyd = meanX_vects_ALL.(['f' num2str(locFly_InAllFlies)]).xytyd;
xytyd_bottom = xytyd( locsBotTrialT );

% decision vector for bototm trials:
Dec = meanX_vects_ALL.(['f' num2str(locFly_InAllFlies)]).Dec;
Dec_bottom = Dec( locsBotTrialT );



%% Plot example Brownian *PROCESSED, UNCORRECTED* trajectories:
% bounding polygon (defined in drosophila_selectFlies.m):

figProcYXInPoly = figure;

% Plot trajectories of bottom trials:
for t=2:maxNTrial
    x = Turncell{locFly_InAllFlies,t,1}; 
    y = Turncell{locFly_InAllFlies,t,2};
    % Compute processed x,y (y = 1 is the upper bound of the bottom arm,
    % just before the intersection):
    xProcessed = ( double(x) - widCenFly(locFly_InSelected) ) / ...
        (lenBotFly(locFly_InSelected) * ratioLenWoW);
    yProcessed = ( double(y) - minYPolyFly(locFly_InSelected) ) / ...
        (lenBotFly(locFly_InSelected) * ratioLenWoW); 
    if 1==1
        % Plot proccessed data (y = 1 <-- upper bound of the bottom arm):
        figure(figProcYXInPoly);
        patchline( double(yProcessed)', double(xProcessed)', ...
            'linestyle', '-', 'edgecolor', [.5,.5,.5], 'linewidth', 0.5, ...
            'edgealpha', 0.5 ); 
        hold on;
    end
end

% Plot example trial trajectory:
expt = 14;
x = Turncell{locFly_InAllFlies,expt,1};
y = Turncell{locFly_InAllFlies,expt,2};
xProcessed = ( double(x) - widCenFly(locFly_InSelected) ) / ...
    (lenBotFly(locFly_InSelected) * ratioLenWoW);
yProcessed = ( double(y) - minYPolyFly(locFly_InSelected) ) / ...
    (lenBotFly(locFly_InSelected) * ratioLenWoW); 
c = [turbo(length(xProcessed)+1)];
patch('xdata', [yProcessed', nan]', ...   
    'ydata', [xProcessed', nan]', ...
    'facevertexcdata', c, ...
    'edgecolor', 'flat', 'EdgeAlpha', .3, 'lineWidth', 1 );
hold on; 
% Add colorbar
xColBar = linspace(0.1,.9,length(xProcessed))';
yColBar = .8 * ones(size(xColBar));
patch('xdata', [xColBar', nan]', ...   
    'ydata', [yColBar', nan]', ...
    'facevertexcdata', c, ...
    'edgecolor', 'flat', 'EdgeAlpha', .75, 'lineWidth', 7 );
hold on; 
text( xColBar(1), yColBar(1)+.1, 'Raw trial duration \rightarrow' );


% axes + save figure:
figure(figProcYXInPoly);
axis image; xlabel('y'); ylabel('x'); grid on;
ggg = gca;
ggg.YDir = 'reverse';
saveas( figProcYXInPoly, [figStartName '_procTrajectories_YX.fig']);



%% Plot example fly *CORRECTED AND TRIMMED* trajectories:

figTrimmedCorrectedYX = figure;
colLeftRight = [0,0,1; 1,0,0];

% Plot the bounding poly of all trajectories / actual maze bounds:
xyAllSim = cell2mat( squeeze( Turncell(locFly_InAllFlies,:,1:2) ) );
boundaryLocs = boundary( xyAllSim(:,1), xyAllSim(:,2), 0.9 );
xProcessedBound = ( double(xyAllSim(boundaryLocs,1)) - ...
    widCenFly(locFly_InSelected) ) / ...
    (lenBotFly(locFly_InSelected) * ratioLenWoW);
yProcessedBound = ( double(xyAllSim(boundaryLocs,2)) - ...
    minYPolyFly(locFly_InSelected) ) / ...
    (lenBotFly(locFly_InSelected) * ratioLenWoW); 
patch( yProcessedBound, xProcessedBound, [.5,.5,.5], 'faceColor', ...
    'none', 'edgeColor', 'k', 'lineWidth', 1 ); hold on;
patch( -yProcessedBound, xProcessedBound, [.5,.5,.5], 'faceColor', ...
    'none', 'edgeColor', 'k', 'lineWidth', 1 ); hold on;


% Plot trajectories of bottom trials:
for t = 1:size(xytyd_bottom,1)
    trialData = xytyd_bottom{t};
    trialDecision = Dec_bottom(t);
    % Note that the 1st coordinate in 'trialData' is used only for -yInf or 
    % -tInf TPI derivations and is thus not plotted below:
    xCorrectedTrimmed = trialData(2:end,1);
    yCorrectedTrimmed = trialData(2:end,2);
    figure(figTrimmedCorrectedYX);
    % Down the arm (negative y):
    patchline( yCorrectedTrimmed(yCorrectedTrimmed<0)', ...
        xCorrectedTrimmed(yCorrectedTrimmed<0)', ...
        'linestyle','-','edgecolor', ...
        [.5,.5,.5], 'linewidth', 0.5, ...
        'edgealpha', 0.5 ); 
    hold on;
    % Up the arm (positive y):
    patchline( yCorrectedTrimmed(yCorrectedTrimmed>0)', ...
        xCorrectedTrimmed(yCorrectedTrimmed>0)', ...
        'linestyle','-','edgecolor', ...
        [.5,.5,.5], 'linewidth', 0.5, ...
        'edgealpha', 0.5 ); 
    hold on;
end

% Plot example trial trajectory (corrected+trimmed):
trialData = xytyd{expt};
xCorrectedTrimmed = trialData(2:end,1);
yCorrectedTrimmed = trialData(2:end,2);
tCorrectedTrimmed = trialData(2:end,3);
c = [turbo(length(xCorrectedTrimmed)+1)];
patch('xdata', [yCorrectedTrimmed', nan]', ...   
    'ydata', [xCorrectedTrimmed', nan]', ...
    'facevertexcdata', c, ...
    'edgecolor', 'flat', 'EdgeAlpha', .3, 'lineWidth', 1 );
hold on; 
% Add colorbar
xColBar = linspace(-0.8,.8,length(xCorrectedTrimmed))';
yColBar = .8 * ones(size(xColBar));
patch('xdata', [xColBar', nan]', ...   
    'ydata', [yColBar', nan]', ...
    'facevertexcdata', c, ...
    'edgecolor', 'flat', 'EdgeAlpha', .75, 'lineWidth', 7 );
hold on; 
text( xColBar(1), yColBar(1)+.1, 'Trial duration (trimmed) \rightarrow' );

% axes + save figure:
figure(figTrimmedCorrectedYX);
axis image; 
xlabel('y'); 
ylabel('x'); 
grid on;
ggg = gca;
ggg.YDir = 'reverse';
saveas( figTrimmedCorrectedYX, [figStartName '_correctedTrajectories_YX.fig']);


