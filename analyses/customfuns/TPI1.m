function [TPI_fly,TPI_SEM_fly, pRgR_nLocR_pRgL_nLocL, edgesVect, ...
    TPI_fly_mpInf, TPI_SEM_fly_mpInf] = ...
    TPI1( traj_cell, edgesVect, DecisionVect, binType, sub100Y )

% COMPUTES THE TPI - TURN PREDICTIVENESS INDEX - OVER Y-BINS OR TIME-BINS -
% FOR AN INDIVIDUAL ANIMAL.
%
% Note that TPI centers x-locations around zero whereas TPI1 is used for 
% centered data.
%
% INPUTS:
%   traj_cell: Tx1 cell array, where  T = # trials. traj_cell{t} is FxK
%       matrix, where F is the number of frames in that trial and columns  
%       1:3 are the (x,y,time) trajectories (ordered), resp.. 
%       Note that time here denotes realtive (rather than absolute) time.
%   edgesVect: 1x(N+1) vector of ordered bin edges, where N = # bins.
%   DecisionVect: 1xT vector with Right/Left turn decisions (corresp. 
%       +1/-1).
%   binType: 'y' for y-bin TPI,
%           'time' for time-bin TPI,
%           'percentArm' for y-bin in %[arm traversed]. Not that %arm is
%           w.r.t the subject's traj's, i.e., not reaching down edge would
%           translate to 100% that span a shorter range of y-values.
%   sub100Y: [relevant only for binType='percentArm'] upper edge of 
%       subject's bottom arm (excluding intersection).
%
% OUTPUTS:
%   TPI_fly: 1xN vector of TPI values for each bin (def. by edgesVect).
%   TPI_SEM_fly: 1xN vector of corresp. SEM for TPI values.
%   pRgR_nLocR_pRgL_nLocL: saved by-bin probabilities and N's (for 
%       illustration of how TPI is computed).
%   TPI_fly_mpInf: 1x2 vector of TPI values conditioned on the first/last
%       xlocations in a trial.
%   TPI_SEM_fly_mpInf: 1x2 vector of corresp. SEM for TPI values.
%
% Copyright (c) Lior Lebovich, 2024
% lebovich.lior@gmail.com


if ~exist('sub100Y','var')
    sub100Y = []; 
end

nGoodTrials = size(traj_cell,1);

% Define the column from which data should be read (and TPI computed). 
% In what follows, this dimension is termed w.
if strcmp( binType, 'y' ) || strcmp( binType, 'percentArm' )
    binCol = 2;
elseif strcmp( binType, 'time' )
    binCol = 3;
end


% Compute x-center param:
xGoodMedian_allTrajs = 0; % no centering required


% Compute, for each trial and w-bin, a binary variable, 
% binaryXLoc_w_byTrial, denoting whehter mean(x in w-bin and trial) is 
% greater/samller than (+1/-1 resp.) median(x in w-bin and ALL trials):

binaryXLoc_w_byTrial = nan(nGoodTrials,length(edgesVect)-1);
binaryXLoc_w_byTrial_minusInf = nan(nGoodTrials,1);
binaryXLoc_w_byTrial_plusInf = nan(nGoodTrials,1);
for tt = 1:nGoodTrials
    
    % load x,w data of current trial:
    xytyd_trial = traj_cell{tt};
    x = xytyd_trial(:,1);
    w = xytyd_trial(:,binCol);

    % Center x around 0:
    xCentered = x - xGoodMedian_allTrajs; 

    % Use the centered X to compute, for each w-bin, binary mean(x), 
    % denoting whether mean(x in given trial and w-bin) is larger(+1) or 
    % smaller(-1) than zero (the horizontal midline):
    for dwPosTrial = 1:(length(edgesVect)-1)
        locsp = (w >= edgesVect(dwPosTrial)) & ...
            (w < edgesVect(dwPosTrial+1)); 
        % Store binary loc. w.r.t xCentered IFF bin was visited in that
        % trial:
        if sum(locsp)>=1
            binaryXLoc_w_byTrial(tt,dwPosTrial) = ...
                ( mean( xCentered(locsp)) > 0 ) - ...
                ( mean( xCentered(locsp) ) < 0 );
        end
    end

    % Also compute binaryXLoc_w_byTrial for first and last loc. (-/+Inf):
    binaryXLoc_w_byTrial_minusInf(tt) = ( mean( xCentered(1)) > 0 ) - ...
                ( mean( xCentered(1) ) < 0 );
    binaryXLoc_w_byTrial_plusInf(tt) = ( mean( xCentered(end)) > 0 ) - ...
                ( mean( xCentered(end) ) < 0 );

end


% Compute TPI(w) by running over w-bins:

TPI_fly = nan(1,length(edgesVect)-1);
TPI_SEM_fly = nan(1,length(edgesVect)-1);
pRgR_nLocR_pRgL_nLocL = nan(4,length(edgesVect)-1);


for dwPosTrial = 1:length(edgesVect)-1
    
    % Compute P( decision = Right | xsInBin = Right ):
    n_LocR = sum( binaryXLoc_w_byTrial(:,dwPosTrial) == 1, ...
        'omitnan' ); % denominator
    n_LocR_AND_TurnR = sum( (binaryXLoc_w_byTrial(:,dwPosTrial) == 1) & ...
        (DecisionVect == 1), 'omitnan' ); % numerator
    pRight_givenLocR = n_LocR_AND_TurnR ./ n_LocR;
    
    % Compute P( decision = Right | xsInBin = Left ):
    n_LocL = sum( binaryXLoc_w_byTrial(:,dwPosTrial) == -1, ...
        'omitnan' ); % denominator
    n_LocL_AND_TurnR = sum( (binaryXLoc_w_byTrial(:,dwPosTrial) == -1) & ...
        (DecisionVect == 1), 'omitnan' ); % numerator
    pRight_givenLocL = n_LocL_AND_TurnR ./ n_LocL;

    % TPI = P( dec. = R | xsInBin = R ) - P( dec. = R | xsInBin = L ):
    TPI_fly(1,dwPosTrial) = pRight_givenLocR - ...
        pRight_givenLocL;

    % Compute SEM for corresp. TPI:
    TPI_SEM_fly(1,dwPosTrial) = sqrt( ...
        ( pRight_givenLocR * (1-pRight_givenLocR) / n_LocR ) + ...
        ( pRight_givenLocL * (1-pRight_givenLocL) / n_LocL ) );

    % Save probabilities and N's (for illustration):
    pRgR_nLocR_pRgL_nLocL(:,dwPosTrial) = [pRight_givenLocR; n_LocR; ...
        pRight_givenLocL; n_LocL];

end

% Also compute TPI(w) for first loc. (-Inf):
n_LocR = sum( binaryXLoc_w_byTrial_minusInf == 1, 'omitnan' ); % denominator
n_LocR_AND_TurnR = sum( (binaryXLoc_w_byTrial_minusInf == 1) & ...
    (DecisionVect == 1), 'omitnan' ); % numerator
pRight_givenLocR = n_LocR_AND_TurnR ./ n_LocR;
n_LocL = sum( binaryXLoc_w_byTrial_minusInf == -1, 'omitnan' ); % denominator
n_LocL_AND_TurnR = sum( (binaryXLoc_w_byTrial_minusInf == -1) & ...
    (DecisionVect == 1), 'omitnan' ); % numerator
pRight_givenLocL = n_LocL_AND_TurnR ./ n_LocL;
TPI_fly_minusInf = pRight_givenLocR - pRight_givenLocL;
TPI_SEM_fly_minusInf = sqrt( ...
        ( pRight_givenLocR * (1-pRight_givenLocR) / n_LocR ) + ...
        ( pRight_givenLocL * (1-pRight_givenLocL) / n_LocL ) );

% Also compute TPI(w) for last loc. (+Inf):
n_LocR = sum( binaryXLoc_w_byTrial_plusInf == 1, 'omitnan' ); % denominator
n_LocR_AND_TurnR = sum( (binaryXLoc_w_byTrial_plusInf == 1) & ...
    (DecisionVect == 1), 'omitnan' ); % numerator
pRight_givenLocR = n_LocR_AND_TurnR ./ n_LocR;
n_LocL = sum( binaryXLoc_w_byTrial_plusInf == -1, 'omitnan' ); % denominator
n_LocL_AND_TurnR = sum( (binaryXLoc_w_byTrial_plusInf == -1) & ...
    (DecisionVect == 1), 'omitnan' ); % numerator
pRight_givenLocL = n_LocL_AND_TurnR ./ n_LocL;
TPI_fly_plusInf = pRight_givenLocR - pRight_givenLocL;
TPI_SEM_fly_plusInf = sqrt( ...
        ( pRight_givenLocR * (1-pRight_givenLocR) / n_LocR ) + ...
        ( pRight_givenLocL * (1-pRight_givenLocL) / n_LocL ) );

TPI_fly_mpInf = [TPI_fly_minusInf, TPI_fly_plusInf];
TPI_SEM_fly_mpInf = [TPI_SEM_fly_minusInf, TPI_SEM_fly_plusInf];

end