%% LONG MAZE CROSSINGS AND SIMULATIONS:
%
% Considers the post cul-de-sac midline-crossings made by flies in the long 
% maze and either simulates the expected #crossings in shorter mazes or 
% plots #crossings for example flies in the long maze.
%
% INPUT FILE: 
%   flies_pEvenTrans.mat: w/ outputs of drosophila_probEvenTransition.m.
%
% OUTPUT FILE:
%   simTransByItiDistY.mat: simulated #crossings by drawing from 
%       inter-crossings-deltaY-distances.
%
% FIGURES:
%   Observed #crossing distributions of WT flies in short/long maze vs
%       expected by drawing from long maze inter-crossing-deltaY-distances:
%       Fig. S5B-C: flies_bootTransByLongTrasIntervals.fig
%   Expected #crossing in shorter mazes, considering 1st crossings in
%       long maze (vs observed short/long):
%       Fig. S5D: flies_bootTransByLongCONSECUTIVETrasIntervals.fig
%   Observed #crossings distributions for 5x3 individual flies in long
%       maze:
%       Fig. S5E: fliesLong_indivTransDist.fig
%
% MORE INFORMATION:
%
% The first and second section use the #crossings made in the lond maze to 
%   simulate the #crossings expected in the:
%       First section: short and long mazes.
%       Second section: shorter mazes with different arm length.
%   For these sim's, we use: the #crossings, inter-crossing-y-distances 
%       and y|crossings in the long maze (computed in 
%       drosophila_probEvenTransition.m).
%   Sampling [more info. in the sections]: 
%       First section: computes expected #crossings by sampling 
%           *inter-crossing-y-distances* from the dist. of EACH FLY until
%           the desired INTERVAL is exceeded.
%       Second section: computes expected #crossings by sampling *trials* 
%           from ALL TRIALS (of all flies) in the long maze. The expected 
%           #crossings is given by the first #crossings in the sampled
%           trials before the INTERVAL is exceeded.
%   INTERVAL def. for crossing and resulting simulations (ran in both 1st 
%   and 2nd sections): 
%       Note that in the long maze, the crossings of interest occur in a 
%       2d+k interval between the upper edge of the cul-de-sac and the 
%       upper edge of the intersection, where d is the length of the arm 
%       in a short maze (2d in the long) and k is the length of the 
%       intersection (similar for all mazes). The corresponding INTERVAL 
%       in any maze with increased/decreased arm length is given by
%       r*d+k, where r is the arm-length ratio of any-maze to short-maze. 
%       In each section (1 or 2) we either:
%           (i) 'armAndInter': use the crossing of the long maze 
%               occurring in the 2d+k interval and consider the expected 
%               #crossings that will occur in a r*d+k interval.
%           (ii) 'armOnly': first use the crossing of the long maze 
%               occurring in the 2d interval (between cul-de-sac upper 
%               edge to BOTTOM edge of the intersection) and consider the 
%               expected #crossings that will occur in r*d interval. 
%               Then, add (to this expected number) the real number of 
%               crossings occuring in the k interval (within the 
%               intersection).



%% Use relevant directory and add path to external fucntions and datafiles:

cdName = 'SpatiotempDM';
currentFolder = cd;
if ~endsWith(currentFolder,cdName)
    cd(cdName);
end
addpath( 'analyses', 'datafiles', 'analyses/customfuns', '-frozen');



%% Simulate #crossings by drawing from inter-crossings-deltaY-distances
% (of the long maze):
 
% Here, we simulate #crossings (the number of crossings after the 
% cul-de-sac until the end of the trial) in the short and long mazes by 
% sampling from the between-crossings deltaY's distribution of each fly in 
% the long maze.

% Possible simulations:

% Each simulation runs over all flies in the long maze. 
% For each fly f in a given sim., the expected #crossings is computed n(f) 
% times, where n(f) is the number of trials made by fly f. 
% For each trial, there are two possible simulations (computed separately):
% (1) 'armOnly': simulate #crossings WITHIN THE ARM - by randomly drawing 
    % values from all inter-crossings-deltaY-distances made by that fly 
    % WITHIN THE ARM (deltaY_armOnly_IntTransInt_flies) and stop once the 
    % accumulated distances reach the arm length (2d for the long sim's, or
    % half the length, d for the short sim's). The # of steps - 1 in the 
    % sim is the expected #crossings in the arm. Then, we add to this sim. 
    % value the real number of #crossings occuring IN THE INTERSECTION for
    % that trial.
% (2) 'armAndInter': simulate #crossings AFTER THE CUL-DE-SAC 
    % (arm+intersection) - by randomly drawing values from all 
    % inter-crossings-deltaY-distances made by that fly AFTER THE
    % CUL-DE-SAC (deltaY_IntTransInt_flies, arm+intersection) and stop once 
    % the accumulated distances reach the length between the cul-de-sac and 
    % the intersection (2d+k for the long sim's, or half of the arm length 
    % plus the full length of the intersection, d+k, for the short sim's). 
    % The # of steps - 1 is the sim #crossings in the arm+intersection 
    % (i.e., after the cul-de-sac). Then we add zero this sim. value.

% Note that for both options above, the sampling weight is given by 
    % 1/(#distances in trial) which was computed in 
    % drosophila_probEvenTransition.m (wDeltaY_IntTransInt_flies AND
    % wDeltaY_armOnly_IntTransInt_flies). This guarantees that, within a 
    % given fly, each of its trials contributes equally. 
% Also, note that it is extremely rare to observe a transition after the 
% end of the intersection.

% Load data (outputs of drosophila_probEvenTransition.m):
load( 'flies_pEvenTrans.mat' );
dataType = 'Turndata_LongYy';

% Define #flies to simulate (each fly f is sampled once, with 
% nTrials[fly=f] in each simulation):
numFliesLong = numel( numTrans.Long.yUpArm );

% Define cul-de-sac loc:
if strcmp(dataType,'Turndata_LongYy')
    culDeSacY = 0.68/3.36;
end

% Define maze regions for long:
longArmLength = 1-culDeSacY;
longIntersectionLength = 1.158-1;

% Define simulation types:
mazeTypes = {'short', 'long'};
regionTypes = {'armAndInter', 'armOnly'};
regionPrintNames = {'arm & inter', 'arm'};
mazeTypeCols.short = 'k';
mazeTypeCols.long = 'b';
regionTypeCols.armAndInter = "#D95319";
regionTypeCols.armOnly = "#D95319";

% Define #sims and bin edges for Prob.|#crossings (which will be stored 
% for each simulation).
nSim = 1e3;
binEdges = [-.5:29.5, Inf];
binCenters = 0:30;
binCenters4Plt = reshape( [-.5+binCenters; .5+binCenters], [], 1 )';
% Create the random number stream for reproducibility.
sss = RandStream('mlfg6331_64');

% Load the real #crossings and compute the difference: 
real_numTrans_arm_flies = numTrans.(dataType(10:end-2)).yUpArm;
real_numTrans_intAndArm_flies = numTrans.(dataType(10:end-2)).yAfterCul;
% Iterate over each element in the cell arrays
real_numTrans_int_flies = cell(size(real_numTrans_arm_flies));
real_numTrans_zero_flies = cell(size(real_numTrans_arm_flies));
for k = 1:numel(real_numTrans_int_flies)
    % Compute the difference between the corresponding elements
    real_numTrans_int_flies{k} = real_numTrans_intAndArm_flies{k} - ...
        real_numTrans_arm_flies{k};
    real_numTrans_zero_flies{k} = zeros( size( ...
        real_numTrans_intAndArm_flies{k} ) );
end


% Run/ load this section's simulations:

if exist('simTransByItiDistY.mat')
    % Load simulation output:
    load('simTransByItiDistY.mat','sim_probTrans_allFlies');
else
    % Produce + save simulation output:
    for m = 1:numel(mazeTypes)
        mazeType = mazeTypes{m};
    
        % Define maze type 
        if strcmp(mazeType,'short')
            simArmLength = 0.5 * longArmLength;
        elseif strcmp(mazeType,'long')
            simArmLength = longArmLength;
        end
    
        for r = 1:numel(regionTypes)
            regionType = regionTypes{r};
    
            % Define simulation type
            if strcmp(regionType,'armOnly') % sim. arm & add real #crossings inter
                simTotalLength = simArmLength;
                real_yDistITI = deltaY_armOnly_IntTransInt_flies;
                real_yDistWeights = wDeltaY_armOnly_IntTransInt_flies;
                real_numTrans_add = real_numTrans_int_flies;
    
            elseif strcmp(regionType,'armAndInter') % sim. inter+arm & add zero
                simTotalLength = longIntersectionLength + simArmLength;
                real_yDistITI = deltaY_IntTransInt_flies;
                real_yDistWeights = wDeltaY_IntTransInt_flies;
                real_numTrans_add = real_numTrans_zero_flies;
            end
    
            % Create array to store Prob.|#crossings for each sim
            sim_probTrans_allFlies.(mazeType).(regionType) = nan( ...
                nSim, length(binCenters) );
    
            for s = 1:nSim
                
                % Create empty cell for storing the sim. #crossings of each fly:
                sim_numTrans_flies = cell(numel(selectedFlies),1);
    
                for nF = 1:numFliesLong
    
                    nTrials_fly = numel(real_numTrans_add{nF});
    
                    % Create column vector to store fly's sim. #crossings in 
                    % each trial: 
                    sim_numTrans_fly = nan(nTrials_fly,1); 
    
                    % Load real yDistITI, weights and relevant #crossings to 
                    % add to the sim: 
                    real_numTrans_add_fly = real_numTrans_add{nF};
                    real_yDistITI_fly = real_yDistITI{nF};
                    real_yDistWeights_fly_raw = real_yDistWeights{nF};
                    real_yDistWeights_fly = real_yDistWeights_fly_raw / ...
                        sum(real_yDistWeights_fly_raw);
    
                    % Run over trials and simulate #crossings by yDist. (then add):
    
                    for t = 1:nTrials_fly
    
                        numTrans_simTrial = -1;
                        yDist_simTrial = 0;
                        
                        % Sample distances. We use 1000 to avoid using
                        % "datasample" multiple times in a sim trial:
                        sampledYDistances = datasample( sss, ...
                            real_yDistITI_fly, 1000, 'Weights', ...
                            real_yDistWeights_fly);
    
                        % Run the simulation until yDidt_simTrial >= simTotalLength
                        while yDist_simTrial < simTotalLength
                            % Randomly draw yDistITI value to add to yDidt_simTrial
                            yDist_simTrial = yDist_simTrial + ...
                                sampledYDistances( 2 + numTrans_simTrial );
                            % Add 1 to numTrans_simTrial
                            numTrans_simTrial = numTrans_simTrial + 1;
                        end
                        sim_numTrans_fly(t) = real_numTrans_add_fly(t) + ...
                            numTrans_simTrial;
                    end
    
                    % Store fly's #crossings in trials of current sim.
                    sim_numTrans_flies{nF} = sim_numTrans_fly;
    
                end
    
                % Compute the Prob.|#crossings for current sim:
                sim_numTrans_flies_vect = cell2mat(sim_numTrans_flies);
                P = histcounts( sim_numTrans_flies_vect, 'normalization', ...
                    'probability', 'BinEdges', binEdges );
                sim_probTrans_allFlies.(mazeType).(regionType)(s,:) = P;
    
            end
    
        end
        
    end
    
    save('datafiles/simTransByItiDistY.mat','sim_probTrans_allFlies');
end


% Plot the real and simulated #crossings:

fig_simTransProbByLongItiDistY = figure;

for m = 1:numel(mazeTypes)
    mazeType = mazeTypes{m};
    mazeTypeCol = mazeTypeCols.(mazeType);

    for r = 1:numel(regionTypes)
        regionType = regionTypes{r};
        regionTypeCol = regionTypeCols.(regionType);

        clear plt_leg str_leg;

        subplot(numel(mazeTypes),numel(regionTypes),2*(m-1)+r);
        
        % Plot the real distributions of #crossings after cul-de-sac:
        realDataType = ['Turndata_' upper(mazeType(1)) mazeType(2:end) 'Yy'];
        real_numTransAfterCul_flies = numTrans.(realDataType(10:end-2) ...
            ).yAfterCul;
        real_numTransAfterCul_mat = cell2mat(real_numTransAfterCul_flies);
        P = histcounts( real_numTransAfterCul_mat, 'normalization', ...
            'probability', 'BinEdges', binEdges );
        P4plt = reshape( repmat(P, [2,1]), [], 1 )';
        plt_leg(1) = fill( [-.5, binCenters4Plt, 30.5], [0, P4plt, 0], ...
            'r', 'FaceColor', mazeTypeCol, 'FaceAlpha', 0.3, ...
            'EdgeColor', 'none' );
        str_leg{1} = [mazeType ' maze - data'];
        hold on;
        %avgTransReal = sum(P .*  0:30);
        
        % Plot the simulated distributions of #crossings after cul-de-sac:

        sim_pTrans_mat = sim_probTrans_allFlies.(mazeType).(regionType);
        sim_pTrans_mean = mean(sim_pTrans_mat);
        sim_pTrans_std = std(sim_pTrans_mat);
        sim_pTrans_mean4plt = reshape( repmat(sim_pTrans_mean, [2,1]), [], 1 );
        sim_pTrans_std4plt = reshape( repmat(sim_pTrans_std, [2,1]), [], 1 );

        plt_leg(2) = plot( binCenters4Plt, sim_pTrans_mean4plt, ...
            'Color', regionTypeCol, 'LineWidth', 1 );
        str_leg{2} = ['Sim by long (' regionPrintNames{r} ')'];
        hold on;

        fillPlot = fill( [binCenters4Plt, flip(binCenters4Plt)]', ...
            [sim_pTrans_mean4plt + sim_pTrans_std4plt; ...
            flip(sim_pTrans_mean4plt - sim_pTrans_std4plt)], 'r', ...
            'FaceColor', regionTypeCol, 'FaceAlpha', 0.3, 'EdgeColor', ...
            'none' ); 

        xlim([binCenters(1)-.5, binCenters(end)+.5]);
        legend(plt_leg, str_leg);
        title([mazeType ' maze']);
        xlabel('#midline-crossings');
        ylabel('Probability');
        ylim([0,0.27])

    end

end

% Save the figure:

figure(fig_simTransProbByLongItiDistY);
figStartName = 'analyses/figures/flies_bootTransByLongTrasIntervals';
saveas( fig_simTransProbByLongItiDistY, [figStartName '.fig']);



%% Consider #(1st crossings) in long maze to compute expected #crossing in shorter mazes 

% Here, we use the #midline-crossing in the long maze (during upward 
% motion, between the upper edge of the cul-de-sac and the upper edge of
% the intersection) to bootstrap trials and consider the #midline-crossing 
% expected in shorter mazes.
% Specifically, we consider all trials made by all flies in the long maze.
% In each simulation, we boostrap n trials (the number of trials made by 
% all flies in the long maze, with replacements). 
% For each sampled trial, the computed number of crossings expected from a
% shorter maze is given by the number of crossings in the long maze that
% occured before the INTERVAL is exceeded (see above for INTERVAL def).
% In the 'armAndInter' sim. type: 
    % we use the number of crossings in the long maze occuring on the 2d+k 
    % interval (between top end of the cul-de-sac to top end of the 
    % intersection) and consider for the shorter maze only the nunber of 
    % crossings that occur before r*d+k is exceeded.
% In the 'armOnly' sim. type: 
    % we use the number of crossing in the long maze occuring on the 2d 
    % interval (between top end of the cul-de-sac to BOTTOM end of the 
    % intersection) and consider for the shorter maze only the number of 
    % crossings that occur before r*d is exceeded. We add to this derived 
    % number of crossings the number of crossings in the long maze that 
    % occured within intersection (the k interval, similar across all 
    % mazes). 

% Load data (outputs of drosophila_probEvenTransition.m):
load( 'flies_pEvenTrans.mat' );
dataType = 'Turndata_LongYy'; 

% Define #flies (only for computing the trimmed #crossings of trials in 
% each fly. The Bootstrap procedure will sample trimmed #crossings from all 
% trials of all trials, irrespective of flies' identities):
numFliesLong = numel(selectedFlies);

% Define cul-de-sac loc:
if strcmp(dataType,'Turndata_LongYy')
    culDeSacY = 0.68/3.36;
end

% Define maze regions for long:
longArmLength = 1-culDeSacY;
longIntersectionLength = 1.158-1;

% Define trim types:
mazeTypes = {'long', 'med190', 'med167', 'med150', 'med133', ...
    'med110', 'short', 'med090', 'med067', 'med050', 'med033', 'med010'};
simArmLength_vect = nan(size(mazeTypes));
simArmRatioShort_vect = nan(size(mazeTypes));
regionTypes = {'armAndInter', 'armOnly'};
regionPrintNames = {'arm & inter', 'arm'};
mazeTypeCols.short = 'k';
mazeTypeCols.long = 'b';
mazeTypeCols.med = 'r';
regionTypeCols.armAndInter = "#D95319";
regionTypeCols.armOnly = "#D95319";

% Load the real #crossings and compute the difference: 
real_numTrans_arm_flies = numTrans.(dataType(10:end-2)).yUpArm;
real_numTrans_intAndArm_flies = numTrans.(dataType(10:end-2)).yAfterCul;

% Create figures:
fig_simNumTrans_cutYInterTransDist.armAndInter = figure;
fig_simNumTrans_cutYInterTransDist.armOnly = figure;


% Compute #crossings from Long that fit in [r * Short_arm_length], where r 
% is given by mazeTypes:
for m = 1:numel(mazeTypes)
    mazeType = mazeTypes{m};

    % Define maze type 
    if strcmp(mazeType,'short')
        simArmLength = 0.5 * longArmLength;
        simArmRatioShort_vect(m) = 1;
    elseif strcmp(mazeType,'long')
        simArmLength = longArmLength;
        simArmRatioShort_vect(m) = 2;
    elseif strcmp(mazeType(1:3),'med')
        simArmLength = ( str2double( mazeTypes{m}(4:end) )/ 100 ) * ...
            0.5 * longArmLength; 
        simArmRatioShort_vect(m) = str2double( mazeTypes{m}(4:end) )/ 100;
    end
    simArmLength_vect(m) = simArmLength;

    for r = 1:numel(regionTypes)
        regionType = regionTypes{r};

        % Create array (in struct.) to store trim #crossings:
        trim_numTrans.(mazeType).(regionType) = cell( numFliesLong, 1 );

        for nF = 1:numFliesLong

            % Read real #crossings of current fly in long maze:
            real_numTrans_arm_fly = real_numTrans_arm_flies{nF};
            real_numTrans_intAndArm_fly = real_numTrans_intAndArm_flies{nF};
            real_numTrans_int_fly = real_numTrans_intAndArm_fly - ...
                real_numTrans_arm_fly;
            nTrials_fly = numel(real_numTrans_arm_fly);

            % Create vector to store trimmed #crossings:
            trim_numTrans_fly = nan(nTrials_fly,1);

            % Define trim type
            if strcmp(regionType,'armOnly') % trim arm & add real #crossings inter
                simTotalLength = simArmLength;
                real_yDistITIcc_fly = deltaY_armOnly_IntTransInt_flies_cc{nF};
                real_numTrans_add_fly = real_numTrans_int_fly;
            elseif strcmp(regionType,'armAndInter') % trim inter+arm & add zero
                simTotalLength = longIntersectionLength + simArmLength;
                real_yDistITIcc_fly = deltaY_IntTransInt_flies_cc{nF};
                real_numTrans_add_fly = zeros(size(real_numTrans_int_fly));
            end

            % Run over trials and compute #crossings that fit in the trimmed 
            % long maze: 

            for t = 1:nTrials_fly

                % Compute cumulative inter-crossings-y-distances:
                cums_yDistITI_trial = cumsum( real_yDistITIcc_fly{t} );
                locTrans = find( cums_yDistITI_trial >= simTotalLength, ...
                    1, 'first' );
                % Store trimmed #crossings:
                if isempty(locTrans) && strcmp(mazeType,'long') % sim. long maze
                    trim_numTrans_fly(t) = real_numTrans_add_fly(t);
                else % sim. maze with arm < 2*arm_short
                    trim_numTrans_fly(t) = locTrans - 1 + ...
                        real_numTrans_add_fly(t);
                end

            end

            % Store fly's #crossings in trials of current sim.
            trim_numTrans.(mazeType).(regionType){nF} = trim_numTrans_fly;

        end

    end
    
end


% Now, bootstrap of the trimmed #crossings and plot the distributions:

% Define #sims and bin edges for Prob.|#crossings.
nSim = 1e4;
nTrials2Sim = length( cell2mat( numTrans.Short.yUpArm ) );
binEdges = [-.5:29.5, Inf];
binCenters = 0:30;
binCenters4Plt = reshape( [-.5+binCenters; .5+binCenters], [], 1 )';
% Create the random number stream for reproducibility.
sss = RandStream('mlfg6331_64');

% Create arrays to also store the sim sample pEven:
sim_pEvenSample.armAndInter = nan(nSim,numel(mazeTypes));
sim_pEvenSample.armOnly = nan(nSim,numel(mazeTypes));

for m = 1:numel(mazeTypes)
    mazeType = mazeTypes{m};

    for r = 1:numel(regionTypes)
        regionType = regionTypes{r};
        regionTypeCol = regionTypeCols.(regionType);

        % Load trimmed #crossings as vect:
        trim_numTrans_vect = cell2mat( trim_numTrans.(mazeType).( ...
            regionType) );

        % Create array to store Prob.|#crossings for each sim
        sim_probTrans_allFlies2.(mazeType).(regionType) = nan(nSim, ...
            length(binCenters));

        for s = 1:nSim
            % Bootstrap trimmed #crossings
            sampledNumTrans = datasample( sss, trim_numTrans_vect, ...
                nTrials2Sim );
            % Compute Prob.|#crossings
            P = histcounts( sampledNumTrans, 'normalization', ...
                'probability', 'BinEdges', binEdges );
            % Store Prob.|#crossings
            sim_probTrans_allFlies2.(mazeType).(regionType)(s,:) = P;
            % Compute and store pEven in sampled trimmed #crossings
            sim_pEvenSample.(regionType)(s,m) = mean( ...
                rem(sampledNumTrans,2)==0 );
        end


        % Plot #crossings distribution:

        clear plt_leg str_leg;

        figure(fig_simNumTrans_cutYInterTransDist.(regionType));
        subplot(3,5,m);
        
        % Plot the real distributions of #crossings after cul-de-sac:
        if strcmp(mazeType,'short') || strcmp(mazeType,'long')
            realDataType = ['Turndata_' upper(mazeType(1)) ...
                mazeType(2:end) 'Yy'];
            real_numTransAfterCul_flies = numTrans.(realDataType( ...
                10:end-2) ).yAfterCul;
            real_numTransAfterCul_mat = cell2mat( ...
                real_numTransAfterCul_flies);
            P = histcounts( real_numTransAfterCul_mat, 'normalization', ...
                'probability', 'BinEdges', binEdges );
            P4plt = reshape( repmat(P, [2,1]), [], 1 )';
            plt_leg(1) = fill( [-.5, binCenters4Plt, 30.5], [0, P4plt, 0], ...
                'r', 'FaceColor', mazeTypeCols.(mazeType), ...
                'FaceAlpha', 0.3, 'EdgeColor', 'none' );
            str_leg{1} = [mazeType ' maze - data'];
            hold on;
        end
        
        % Plot the simulated distributions of #crossings after cul-de-sac:

        sim_pTrans_mat = sim_probTrans_allFlies2.(mazeType).(regionType);
        sim_pTrans_mean = mean(sim_pTrans_mat);
        sim_pTrans_std = std(sim_pTrans_mat);
        sim_pTrans_mean4plt = reshape( repmat(sim_pTrans_mean, [2,1]), [], 1 );
        sim_pTrans_std4plt = reshape( repmat(sim_pTrans_std, [2,1]), [], 1 );

        plt_leg(2) = plot( binCenters4Plt, sim_pTrans_mean4plt, ...
            'Color', regionTypeCol, 'LineWidth', 1 );
        str_leg{2} = ['Sim by long (' regionPrintNames{r} ')'];
        hold on;

        fillPlot = fill( [binCenters4Plt, flip(binCenters4Plt)]', ...
            [sim_pTrans_mean4plt + sim_pTrans_std4plt; ...
            flip(sim_pTrans_mean4plt - sim_pTrans_std4plt)], 'r', ...
            'FaceColor', regionTypeCol, 'FaceAlpha', 0.3, 'EdgeColor', ...
            'none' ); 

        xlim([binCenters(1)-.5, binCenters(end)+.5]);
        if strcmp(mazeType,'long') || strcmp(mazeType,'short')
            legend(plt_leg, str_leg);
            title([mazeType ' maze']);
        else
            title([num2str( str2double( mazeTypes{m}(4:end) )/ 100 ) ...
                ' x short arm maze']);
        end
        xlabel('#midline-crossings');
        ylabel('Probability');
        

        % Also plot pEven in trimmed #crossings samples:
        figure(fig_simNumTrans_cutYInterTransDist.(regionType));
        subplot(3,5,13:15);
        sim_pEvenSample_regionMaze = sim_pEvenSample.(regionType)(:,m);
        errorbar( simArmRatioShort_vect(m), mean(sim_pEvenSample_regionMaze), ...
            std(sim_pEvenSample_regionMaze), 'Color', regionTypeCol, ...
            'LineWidth', 2, 'LineStyle', 'none', 'CapSize', 0 );
        hold on;

    end

end


% Add real data pEven to pEven in trimmed #crossings samples:

for r = 1:numel(regionTypes)
    regionType = regionTypes{r};
    figure(fig_simNumTrans_cutYInterTransDist.(regionType));
    subplot(3,5,13:15);
    plot( simArmRatioShort_vect, mean( sim_pEvenSample.(regionType) ), ...
        'Color', regionTypeCol, 'LineWidth', 2 );
    hold on;
    pl(1) = yline( mean( rem( cell2mat(numTrans.Short.yUpArm), 2 ) == 0 ), ...
        'k:', 'LineWidth', 2 );
    hold on;
    pl(2) = yline( mean( rem( cell2mat(numTrans.Long.yUpArm), 2 ) == 0 ), ...
        'b:', 'LineWidth', 2 );
    title(regionType);
    xticks(flip(simArmRatioShort_vect))
    xlabel('sim maze arm = Short arm x ');
    ylabel('pEven sim. sample');
    legend(pl,{'pEven Short', 'pEven Long'});
end


% Save the figure (for armAndInter): 
figure(fig_simNumTrans_cutYInterTransDist.armAndInter);
figStartName = ...
    'analyses/figures/flies_bootTransByLongCONSECUTIVETrasIntervals';
saveas( fig_simNumTrans_cutYInterTransDist.armAndInter, [figStartName '.fig']);



%% Plot #crossings distributions in long maze for 5x3 example flies:
% Example flies: 5 with largest/lowest/intermediate pEven:


load( 'flies_pEvenTrans.mat' );
dataType = 'Turndata_LongYy';

% Define cul-de-sac and intersection upper edge loc's.
if strcmp(dataType,'Turndata_LongYy')
    culDeSacY = 0.68/3.36;
    intersectionEnd = 1.158;
end

% Define #crossings bin edges:
binEdges = [-.5:29.5, Inf];

% Read post cul-de-sac #crossings, pEven ( P(#crossing is even) ) and 
% y-values of crossing:
regName = 'yAfterCul';
numTranCell = numTrans.(dataType(10:end-2)).(regName);
pEvenVect = pEven_with0_strct.yAfterCul.(dataType(10:end-2));
nFlies = length(pEvenVect);
yOfTrans = yVals_trans_flies;

% Select 5x3 flies, based on pEven values:
qLarge_pEven = quantile( pEvenVect, 1-(1/11) );
qSmall_pEven = quantile( pEvenVect, 1/11 );
qMidD_pEven = quantile( pEvenVect, .5 - (.5/11) );
qMidU_pEven = quantile( pEvenVect, .5 + (.5/10) );
exmplFliesNums_pEven.large = find( pEvenVect >= qLarge_pEven );
exmplFliesNums_pEven.small = find( pEvenVect <= qSmall_pEven );
exmplFliesNums_pEven.mid = find( pEvenVect <= qMidU_pEven & ...
    pEvenVect >= qMidD_pEven );
pSizeNames = {'large', 'mid', 'small'};
pSizeCols.large = [.05 .23 .6];
pSizeCols.mid = [.2 .6 .9];
pSizeCols.small = [.4 .86 .86];

fig_numTrans = figure; % For plotting #crossing dist's of indiv. flies
fig_yTransDist = figure; % For plotting y|crossing dist's of indiv. flies

for p = 1:numel(pSizeNames)
    pSizeName = pSizeNames{p};
    pSizeCol = pSizeCols.(pSizeName);
    fliesNums = exmplFliesNums_pEven.(pSizeName);
    for f = 1:numel(fliesNums)
        flyNum = fliesNums(f);
        numTransFly = numTranCell{flyNum};
        yOfTransFly = yOfTrans{flyNum};

        % Plot #crossing dist's of indiv. flies:
        figure(fig_numTrans);
        subplot(3,5,5*(p-1)+f);
        histogram( numTransFly, 'normalization', 'count', ...
            'BinEdges', binEdges, 'FaceColor', pSizeCol );
        xlabel(['#crossings | ' regName(2:end)]);
        title(['pEven = ' num2str( 1e-3 * round( 1e3 * ...
            pEvenVect(flyNum) ) )]);

        % Plot y|crossing dist's of indiv. flies
        figure(fig_yTransDist);
        subplot(3,5,5*(p-1)+f);
        histogram( yOfTransFly, 'normalization', 'count', ...
            'BinEdges', -.0125:.05:1.4125, 'FaceColor', pSizeCol );
        hold on;
        xline(intersectionEnd, 'k--'); 
        hold on;
        xline(1, 'k--'); 
        hold on;
        xline(culDeSacY, 'k--'); 
        xlabel(['y | crossing ' regName(2:end)]);
        title(['pEven = ' num2str( 1e-3 * round( 1e3 * ...
            pEvenVect(flyNum) ) )]);

    end
end


% Save the #crossings figure:
figure(fig_numTrans);
figStartName = ['analyses/figures/flies' dataType(10:end-2) ...
    '_indivTransDist'];
saveas( fig_numTrans, [figStartName '.fig']);


