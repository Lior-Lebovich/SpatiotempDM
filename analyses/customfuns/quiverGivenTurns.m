function [counts, avgDirThetaRad] = quiverGivenTurns( X, Y, xEdges, ...
    yEdges, xytydBotDec, xGoodMedian, subsamp )

% COMPUTES A DENSITY-BASED QUIVER OF A SINGLE FLY, FOR A GIVEN TURN 
% DIRECTION AND REGION OF INTEREST (RoI).
%
% INPUTS:
%   X, Y, xEdges, yEdges: 2D grid for bins of the region of interest for 
%       which the quiver will be computed such that [X, Y] = 
%       meshgrid( xEdges, yEdges ).
%   xytydBotDec: Tx1 cell array, where  T = # trials. xytydBotDec{t} is FxK
%       matrix, where F is the number of frames in that trial and columns  
%       1:2 are the (x,y) trajectories (ordered), resp.. 
%   xGoodMedian: 1x1 double required for centering x values s.t the entire
%       maze will be horizontally aligned (x=0 corresp. to maze's center). 
%   subsamp: for subsampling rate. If =1 then no subsampling.
%
% OUTPUTS:
%   counts: size(X) matrix of #frames spent in each 2D bin.
%   avgDirThetaRad: size(X) matrix of magnitude-weighted avgerage direction
%       in each 2D bin.

% ADDITIONAL INFORMATION:
% The region of interest (typically, the cul-de-sac) was binned into 2D
    % bins and the grid is (represented by X ans Y, inputs).
% The current function runs over specific trials (turns) of a specific 
    % fly from the bottom arm (all trials in xytydBotDec, input) and over 
    % datapoints within that trial.
% For each n=2:N datapoint in the trial (N=#frames in trial; n=current 
    % frame) that falls within the region of interest:
    % (1) The velocity vector is computed with the previous datapoint.
    % (2) The dir. of the vel. vector's is coverted to angle in radians and
        % appended to the vector of all angles of the 2D bin corresp. to
        % the x,y loc. of the n-1 frame.
    % (3) The magnitude of the vel. vector's is computed and then stored 
        % similary.
% The above results in 2 cell arrays, of angles and magnitudes (dirThetaRad 
    % and magnitudes, resp.) where the {yBin,xBin} loc. in each cell is a 
    % vector storing all angles or magnitudes in that bin, for all of the 
    % trials in xytydBotDec. 
% The code then runs over all 2D bins in the RoI. 
% For each bin:
    % (1) The bin's angles are read in complex plane representation.
    % (2) The corresp. bin's magnitudes are read. 
    % (3) The bin's magnitude-weighted avgerage is computed and 
        % representation is coverted back to (and stored as) average 
        % angle in radians (avgDirThetaRad, output). 
% The number of frames spent in each 2D bin (counts, output) is also stored. 
% Finally, 
% drosophila_quiverIndivs.m will use avgDirThetaRad and counts to
    % plot, for the 2D bins in the RoI, vectors (converted back to 
    % Cartesian coordinates) with: 
    % (1) Vector's direction corresp. to the magnitude-weighted angle. 
    % (2) Vector's magnitude corresp. to the relative #frames spent in that  
        % 2D bin.
    % See: drosophila_quiverIndivs.m.
% drosophila_quiverAvgs.m will weight the avgDirThetaRad and counts of all
    % flies to compute the AVERAGE quiver over all flies with:
    % (1) Vector's direction corresp. to the frequncy-weighted angle. 
    % (2) Vector's magnitude corresp. to the relative avg. frequency spent 
        % in that 2D bin.
    % See: drosophila_quiverAvgs.m.
%
% Copyright (c) Lior Lebovich, 2024
% lebovich.lior@gmail.com


% Initialize variables for accumulating the vectors and counts:
counts = zeros(size(X));
dirThetaRad = cell( size(X) );
magnitudes = cell( size(X) );
avgDirThetaRad = nan(size(X));


% Loop over trials:

for t = 1:numel(xytydBotDec)

    % Read the x, y coordinates for the current trial:
    xy_raw = xytydBotDec{t}( :, 1:2 );

    % Downsample:
    xyNotCentered = xy_raw( 1:subsamp:end, : ); 

    % Center x coordinates s.t x=0 corresp.to center (horizontal):
    xy = [xyNotCentered(:,1)-xGoodMedian, xyNotCentered(:,2)];

    % Find the bin indices for each coordinate:

    xBin = discretize( xy(:, 1), xEdges, 'IncludedEdge', 'right');
    yBin = discretize( xy(:, 2), yEdges, 'IncludedEdge', 'right');


    % Accumulate the velocity vectors and counts (counts denote 
    % #frames spent in a specific bin) in each trial by running 
    % over datapoints (x,y pairs) in traj. of current trial:

    for n = 2:size(xy, 1)
        % Consider only datapoints that are within def. bins:
        if ~isnan(xBin(n)) && ~isnan(yBin(n)) && ...
                ~isnan(xBin(n-1)) && ~isnan(yBin(n-1))
            % Compute counts (#frames spent in a specific bin):
            counts(yBin(n-1), xBin(n-1)) = counts( yBin(n-1), ...
                xBin(n-1) ) + 1;
            % Compute dx, dy:
            dx = xy(n, 1) - xy(n-1, 1);
            dy = xy(n, 2) - xy(n-1, 2);
            % Store dir. of velocity vect. as angle in radians:
            dirThetaRad{yBin(n-1),xBin(n-1)} = ...
                [dirThetaRad{yBin(n-1), xBin(n-1)}; ...
                atan2(dy, dx)];
            % Store velocity vect. magnitude:
            magnitudes{yBin(n-1),xBin(n-1)} = ...
                [magnitudes{yBin(n-1), xBin(n-1)}; ...
                sqrt(dx^2 + dy^2)];
        end
    end
end


% Compute average angle in each bin given entire data:

for i = 1:size(avgDirThetaRad,1)
    for j = 1:size(avgDirThetaRad,2)

        % Read angles and magnitudes of current bin:
        anglesRadInBin = dirThetaRad{i,j};
        magsInBin_weights = magnitudes{i,j};
        
        % Use complex representation for angles of current bin:
        complexAngInBin = exp(1i * anglesRadInBin);

        % Compute (magnitude-based) weighted avg. of the angle in 
        % complex representation:
        weightedMeanAngleComplex = sum( complexAngInBin .* ...
            magsInBin_weights ) / sum(magsInBin_weights);

        % Convert avg. to angle in radians:
        weightedMeanAngleRad = angle(weightedMeanAngleComplex);
        
        % Store bin's avg. angle in radians:
        avgDirThetaRad(i,j) = weightedMeanAngleRad;

    end
end


end