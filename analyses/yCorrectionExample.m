%% ILLUSTRATION OF TRAJECTORIES NORMALIZATION: 
%
% Plots the noramlized bounding polygons in the short WT vs long WT 
% datasets as well as the further correction required so that the short and
% long mazes size will reflect absolute proportions.
%
% INPUT FILES: 
%   Turndata_[DATASET_NAME]_visInsp.mat
%   Turndata_[DATASET_NAME]_edgesDef.mat
%   poly360_Turndata_[DATASET_NAME].mat
%   selectedFlies360DBTurndata_[DATASET_NAME].mat
%   saved_yMinByFly360_[DATASET_NAME].mat
%
% EXTERNAL FUNCTION:
%   This code uses an external function, patchline.m, Reference: 
%   Brett Shoelson (2023). Patchline 
%   (https://www.mathworks.com/matlabcentral/fileexchange/36953-patchline), 
%   MATLAB Central File Exchange. Retrieved September 18, 2023.
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



%% Plots noramlized bounding polygons with corrected y s.t end of arm 
% (before intersection) = 1.

dataTypes = {'Turndata_ShortYy', 'Turndata_LongYy'};
cols = {'k', 'b', 'r'};
ratioLenWoWs = [.8024, .8705, .8024];

figure;
for d = 1:numel(dataTypes)
    dataType = dataTypes{d};
    ratioLenWoW = ratioLenWoWs(d);
    col = cols{d};
    load([dataType '_edgesDef.mat']); 
    load([dataType '_visInsp']); 
    load(['poly360_' dataType '.mat']);
    load(['selectedFlies360DB' dataType '.mat']); 
    load(['saved_yMinByFly360_' dataType(10:end) '.mat']);
    % For later y-normalization (bottom arm BEFORE intesection--> 1) and
    % x-centering:
    lenBotFly = nan(size(selectedFlies));
    widCenFly = nan(size(selectedFlies));
    for nF = 1:length(selectedFlies)
        f = selectedFlies(nF);
        % Compute fly's bottom arm length WITH intersection:
        pol = polyAllFlies{nF}; % The fly's raw traj's bounding polygon.
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
        yUpLimBot_fly = min( polY(polyLocs) );
        % Compute fly's bottom arm length WITHOUT intersection:
        lenBotFly(nF) = yUpLimBot_fly - min(polY);
        widCenFly(nF) = .5 * max(polX);
        % Plot:
        patchline( ( (polY-min(polY)) / ( lenBotFly(nF) * ratioLenWoW)) ...
            / 1, (polX - widCenFly(nF) ) / (lenBotFly(nF) * ...
            ratioLenWoW), 'edgecolor', col, 'linewidth', .5, ...
            'edgealpha',0.2);
        hold on;
    end
end
axis image;



%% Plot the further correction required s.t the short and long mazes size 
% will reflect absolute proportions.

dataTypes = {'Turndata_ShortYy', 'Turndata_LongYy'};
cols = {'k', 'b', 'r'};
ratioLenWoWs = [.8024, .8705, .8024];
ratioSizes = [1, 2.53/3.86];

figure;
for d = 1:numel(dataTypes)
    dataType = dataTypes{d};
    ratioLenWoW = ratioLenWoWs(d);
    ratioSize = ratioSizes(d);
    col = cols{d};
    load([dataType '_edgesDef.mat']); 
    load([dataType '_visInsp']); 
    load(['poly360_' dataType '.mat']);
    load(['selectedFlies360DB' dataType '.mat']); 
    load(['saved_yMinByFly360_' dataType(10:end) '.mat']);
    % For later y-normalization (bottom arm BEFORE intesection--> 1) and
    % x-centering:
    lenBotFly = nan(size(selectedFlies));
    widCenFly = nan(size(selectedFlies));
    for nF = 1:length(selectedFlies)
        f = selectedFlies(nF);
        % Compute fly's bottom arm length WITH intersection:
        pol = polyAllFlies{nF}; % The fly's raw traj's bounding polygon.
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
        yUpLimBot_fly = min( polY(polyLocs) );
        % Compute fly's bottom arm length WITHOUT intersection:
        lenBotFly(nF) = yUpLimBot_fly - min(polY);
        widCenFly(nF) = .5 * max(polX);
        % Plot:
        patchline( ( (polY-min(polY)) / ( lenBotFly(nF) * ratioLenWoW) ...
            ) / ratioSize, ((polX - widCenFly(nF) ) / (lenBotFly(nF) * ...
            ratioLenWoW) ) / ratioSize, 'edgecolor', col, 'linewidth', ...
            .5, 'edgealpha',0.2);
        hold on;
    end
end
axis image;


