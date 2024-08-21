%% AVERAGE BOUNDING POLYGON ACROSS FLIES IN DATASET:
%
% Computes normalized and x-centered avg. bounding polygon across flies in 
% dataset.
%
% INPUT FILE: 
%   output_drosophila_main_[DATASET_NAME(10:END)].mat
%
% OUTPUT FILE:
%   flies_boundPolProc_[DATASET_NAME(10:END)].mat: avg. bounding polygon
%       across flies in dataset.
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



%% Compute norm. and x-centered bounding polygon, per fly and dataset avg.:

dataTypes = {'Turndata_ShortYy', 'Turndata_LongYy', ...
    'Turndata_ShortNorpAYy', 'Turndata_ShortDumbYy', ...
    'Turndata_ShortFoxPYy', 'Turndata_ShortNompCYy'};

% Define 2D bin threshold (see below):
percentThresh = .1;

for d = 1:numel(dataTypes) 

    dataType = dataTypes{d};
    
    % Load processed data:
    matFileName = ['output_drosophila_main_' dataType(10:end) '.mat'];
    load(matFileName);

    % Create cell for storing normalized and x-centered bounding polygonS
    % of all flies in dataset: 
    polyNormAllFlies = cell(numel(selectedFlies),1);

    % Store x-(horizontal) centering:
    xForCentering = nan( numel(selectedFlies), 1 );

    % Read fly data, compute x-centered and normalized bounding polygon, 
    % combine with those of all flies:
    for nF = 1:length(selectedFlies)
        f = selectedFlies(nF);
        flyName = ['f' num2str(f)];

        % Read individual animal's processed data for bottom trials: 
        % Location of bottom trials:
        locsBotTrialT = meanX_vects_ALL.(flyName).locsBotTrialT;
        % Corrected trajectories in bottom trials: 
        xytyd = meanX_vects_ALL.(flyName).xytyd;
        xytyd_bottom = xytyd( locsBotTrialT );

        % Compute fly's med x loc. (over bottom trials):
        traj_cell_mat = cell2mat( xytyd_bottom );
        xy_bottom_mat = traj_cell_mat(:,1:2);
        uniqueXY = unique( xy_bottom_mat, 'rows' );
        xGoodMedian_allTrajs = mean( ...
            [min( uniqueXY( abs(uniqueXY(:,2)) < 1, 1 ) ), ...
            max( uniqueXY( abs(uniqueXY(:,2)) < 1, 1 ) )] );
        % Store xMed required for x-centering:
        xForCentering(nF) = xGoodMedian_allTrajs;

        % Read fly's not normalized, not centered bounding polygon:
        pol = polyAllFlies{nF};
        
        % Compute x-centered and normalized bounding polygon:
        xPol = ( pol(:,1) - widCenFly(nF) ) / ...
            (lenBotFly(nF) * ratioLenWoW) - xForCentering(nF);
        yPol = ( pol(:,2) - minYPolyFly(nF) ) / ...
            (lenBotFly(nF) * ratioLenWoW);
        % Store:
        polyNormAllFlies{nF} = [xPol, yPol];
    end

    % Read norm. and x-centered bounding polygons of all flies:
    combinedXY = cell2mat(polyNormAllFlies);
    combinedX = combinedXY(:,1);
    combinedY = combinedXY(:,2);
    
    % Create a 2D grid representing the entire maze area
    xGrid = min(combinedX):0.001:max(combinedX);
    yGrid = min(combinedY):0.001:max(combinedY);
    [X, Y] = meshgrid(xGrid, yGrid);
    
    % Initialize a counter for each grid cell
    cellCounter = zeros(size(X));
    
    % Count the number of flies with norm&centered x,y trajectories (all 
    % trials) that visit each bin in the X,Y grid:
    for nF = 1:numel(selectedFlies)
        % Read individual animal's normed and x-centered bounding polygon:
        polyNormAllFlies_fly = polyNormAllFlies{nF};
        xNormCenPol = polyNormAllFlies_fly(:,1);
        yNormCenPol = polyNormAllFlies_fly(:,2);
        % Find all 2D bins within the fly's polygon:
        inPolygon = inpolygon(X, Y, xNormCenPol, yNormCenPol);
        cellCounter = cellCounter + inPolygon;
    end

    % Find 2D bins that are visited by at least %percentThresh of the
    % flies:
    avgX = X( cellCounter > (numel(selectedFlies) * percentThresh) );
    avgY = Y( cellCounter > (numel(selectedFlies) * percentThresh) );
    
    % Plot these bins:
    figure;
    plot(avgX, avgY, 'k.');
    xlabel('X');
    ylabel('Y');
    title([dataType(10:end-2) ' - avg. Bounding Polygon']);
    axis image;
    grid on;
    
    % Compute the boundary of the above bins:
    polyLocs = boundary( avgX, avgY );
    avgPolyX = avgX(polyLocs);
    avgPolyY = avgY(polyLocs);
    
    % Plot the average bounding polygon
    hold on;
    plot(avgPolyX, avgPolyY, 'b', 'LineWidth', 2);
    
    % Save the normalized and x-centered bounding polygon of indiv. flies
    % and the avg. over flies dataset:
    polyNormAvg = [avgPolyX, avgPolyY];
    
    save( ['datafiles/flies_boundPolProc_' dataType(10:end-2) '.mat'], ...
        'polyNormAvg', 'polyNormAllFlies' );

    % Clear unnecesary vars:
    clearvars -except cdName dataTypes percentThresh d;

end


