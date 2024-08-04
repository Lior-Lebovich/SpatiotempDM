function [avgAvgSpeedByWBin, avgSpeedAllByWBin, avgAvgDistByWBin, ...
    avgAvgRelTimeByWBin] = ySpeed( traj_cell, edgesVect, binType, fps )


% COMPUTES AVG. KINEMATIC PARAMETERS, OVER Y-BINS OR RELATIVE-TIME-BINS, 
% FOR AN INDIVIDUAL ANIMAL.
%
% INPUTS:
%   traj_cell: Tx1 cell array, where  T = # trials. traj_cell{t} is FxK
%       matrix, where F is the number of frames in that trial and columns  
%       1:3 are the (x,y,time) trajectories (ordered), resp.. 
%       Note that time here denotes realtive (rather than absolute) time.
%   edgesVect: 1x(N+1) vector of ordered bin edges, where N = # bins.
%   binType: 'y' for y-bin speed,
%           'time' for time-bin speed.
%   fps: frames per second.
%
% OUTPUTS:
%   avgAvgSpeedByWBin: 1xN vector of avg. speed values (a.u), defined by
%       edgesVect. avgAvgSpeedByWBin(k) is the avg. over avg. speed in each
%       trial and w-bin defined by edgesVect(k:k+1).
%   avgSpeedAllByWBin: 1xN vector of avg. speed values (a.u), defined by
%       edgesVect. avgSpeedAllByWBin(k) is the avg. speed over all 
%       trajectories (all trials) for the the w-bin defined by
%       edgesVect(k:k+1).
%   avgAvgDistByWBin: as in avgAvgSpeedByWBin, for distance.
%   avgAvgRelTimeByWBin: as in avgAvgSpeedByWBin, for relative time.


% Define the colum from which data should be read (and TPI computed). In
% what follows, this dimension is termed w.
if strcmp( binType, 'y' )
    binCol = 2;
elseif strcmp( binType, 'time' )
    binCol = 3;
end


% Compute median x:

% Load unique pairs of x,w trajectories of all relevant trials: 
nGoodTrials = size(traj_cell,1);
traj_cell_mat = cell2mat( traj_cell );
if strcmp( binType, 'y' ) || strcmp( binType, 'time' )
    xw_bottom_mat = traj_cell_mat(:,[1,binCol]);
end
uniqueXW = unique( xw_bottom_mat, 'rows' );
% Median over all traj's before intersection:
xGoodMedian_allTrajs = mean( ...
    [min( uniqueXW( abs(uniqueXW(:,2)) < 1, 1 ) ), ...
    max( uniqueXW( abs(uniqueXW(:,2)) < 1, 1 ) )] ); 


% Compute, for each trial and w-bin, the average velocity [average: 
% (1) over trial positions in w-bin of CURRENT trial, then over trials:
avgSpeedByWBin_trialMat = nan(nGoodTrials,numel(edgesVect)-1);
avgDistByWBin_trialMat = nan(nGoodTrials,numel(edgesVect)-1);
avgRelTimeByWBin_trialMat = nan(nGoodTrials,numel(edgesVect)-1);
% (2) all positions in bin (all trials)]:
speedW_trialsCell = cell(nGoodTrials,2);

for tt = 1:nGoodTrials
    
    % load x,y,w data of current trial:
    xytyd_trial = traj_cell{tt};
    x = xytyd_trial(:,1);
    y = abs( xytyd_trial(:,2) );
    w = xytyd_trial(:,binCol);
    
    % Center x so that xCentered = x - median(x within bottom arm):
    xCentered = x - xGoodMedian_allTrajs; 

    % Compute speed:
        
    % Calculate distances between consecutive points
    distances = sqrt(diff(xCentered).^2 + diff(y).^2);
    
    % Find consecutive identical points
    identical_points = (distances == 0);
    
    % Exclude consecutive identical points from distance calculation
    distances(identical_points) = [];
    
    % Calculate speed (distance traveled per unit time)
    % Assuming time intervals are uniform
    speedY = (1/fps) * distances ./ mean(diff(find(~identical_points)));

    
    % Associate with corresponding w-bins:
    wOfSpeedRaw = .5*w(1:end-1) + .5*w(2:end);
    wOfSpeed = wOfSpeedRaw(~identical_points);
    wOfSpeed2 = wOfSpeedRaw(wOfSpeedRaw>=edgesVect(1) & wOfSpeedRaw<=(end));

    
    % (2) Store speed and w's for current trial:
    speedW_trialsCell{tt} = [speedY,wOfSpeed,distances];

    % (1) Compute, for each w-bin, the average speed in CURRENT trial:
    for dwPos = 1:(length(edgesVect)-1)
        locsp = (wOfSpeed >= edgesVect(dwPos)) & ...
            (wOfSpeed < edgesVect(dwPos+1)); 
        % Store speed is positions in w-bin exist:
        avgSpeedCurrTrialW = mean( speedY(locsp) );
        if ~isnan(avgSpeedCurrTrialW)
            avgSpeedByWBin_trialMat(tt,dwPos) = avgSpeedCurrTrialW;
        end
        avgDistCurrTrialW = sum( distances(locsp) );
        if ~isnan(avgDistCurrTrialW)
            avgDistByWBin_trialMat(tt,dwPos) = avgDistCurrTrialW;
        end
        % Relative time:
        locsp2 = (wOfSpeed2 >= edgesVect(dwPos)) & ...
            (wOfSpeed2 < edgesVect(dwPos+1)); 
        avgRelTimeByWBin_trialMat(tt,dwPos) = sum(locsp2) / ...
            length(wOfSpeed2);
    end

end


% (1) Compute avg. speed by averaging over positions in w-bin of current 
    % trial, then over trials:
% Speed:
avgAvgSpeedByWBin = mean( avgSpeedByWBin_trialMat, 'omitnan' );
% Distance:
avgAvgDistByWBin = mean( avgDistByWBin_trialMat, 'omitnan' );
% Relative time:
avgAvgRelTimeByWBin = mean( avgRelTimeByWBin_trialMat, 'omitnan' );


% (2) all positions in bin (all trials):
speedW_allMat = cell2mat( speedW_trialsCell );
speedAll = speedW_allMat(:,1);
wAll = speedW_allMat(:,2);
avgSpeedAllByWBin = nan(length(edgesVect)-1,1);
for dwPos = 1:(length(edgesVect)-1)
    locsp = (wAll >= edgesVect(dwPos)) & ...
        (wAll < edgesVect(dwPos+1)); 
    % Compute and store avg. speed over all positions in bins (all trials):
    speedAllInWBin = speedAll(locsp); 
    avgSpeedAllByWBin(dwPos) = mean( speedAllInWBin );
end



end