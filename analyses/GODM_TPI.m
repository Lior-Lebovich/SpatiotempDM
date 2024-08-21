%% EXTERNAL DATASET: SPATIAL TPI OF FLIES  IN CIRCULAR ARENA:
%
% Computes and plots the trajectories and spatial TPI (Turn Predictiveness 
% Index) of flies in the study of Sridhar and colleagues. Reference: 
% Sridhar, V. H., Li, L., Gorbonos, D., Nagy, M., Schell, B. R., Sorochkin,
% T., Gov, N. S., & Couzin, I. D. (2021). The geometry of decision-making 
% in individuals and collectives. Proceedings of the National Academy of 
% Sciences of the United States of America, 118(50), e2102157118. 
% https://doi.org/10.1073/pnas.2102157118
%
% TPIs are computed across all trials of all flies, separately for two
% conditions (with 2 targets each).
%
% INPUT FILE:
%   rawdata/dataBifurcation.csv: 2 targets dataset.
%
% FUNCTION:
%   TPI1.m: computes the TPI without x-centering (not required).
%
% FIGURE:
%   Trajectories of experimental conditions with 2 targets and corresp. 
%       spatial TPIs of target angle ~= 180:
%       Fig. S12: GODM_trajectoriesAndTpis_byTarAngles.fig
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



%% Load data and defines relevant experimental conditions:

% Load GODM data - rotated x,y for 2 targets: 
T = readtable( 'datafiles/rawdata/dataBifurcation.csv' );

% Select only stimuli experimental sets: 
uniNStimuli = 1:3;
T = T( ismember( T.nStimuli, uniNStimuli), :  );

% Define rotated_post0_x conditions: 
rangePost0X = [0, 4, 4.2, Inf];
post0CondNames = {'strange', 'small', 'large'};

% Compute y-edges (raw: rotated x): 
ylims = [1e-1*floor(1e1*min(T.rotated_x)), 1e-1*ceil(1e1*max(T.rotated_x))];
dy = .3;
yEdges = 0:dy:4.5;
yCenters = .5 * diff(yEdges(1:2)) + yEdges(1:end-1);
infVal = 5;

% Define left/right colors:
colLR = {'b','r'};

% Define target colors:
colTar = {"#EDB120", "#7E2F8E", "#77AC30"};



%% Runs over cond's and flies, plots traj's, computes and plots spatial TPI:

% Reads and stores data of 1<=nStim<=3 trials in two rotated_post0_x 
    % conditions: 4<rotated_post0_x<4.2 and rotated_post0_x>=4.2:
% Data is stored in structures for each conditions. 
% Computes TPI, separately for the 2 targets conditions that aren't 180 deg 
    % and trials that don't pass through rotated-x < 0 (only trials that 
    % appear in black). To plot over all trials change the value of 
    % TPI_trial_type from 'trialsWithRotX>=0' to 'allTrials'.
% Note that due to the small #trials per fly, the spatial TPI is computed
    % across trials made by ALL flies (rather than separately for each 
    % fly).

TPI_trial_type = 'trialsWithRotX>=0';

targetsLocs = [];
targetsAngles = [];
targetsRe = [];

fig_trajTPI = figure;

for r = 1:length(rangePost0X)-1

    post0CondName = post0CondNames{r};

    trajSubplot = subplot(2,numel(uniNStimuli),r);
    
    % Read rotated_post0_x limits:
    pst0x_d = rangePost0X(r);
    pst0x_u = rangePost0X(r+1);
    
    % Read data for current rotated_post0_x condition:
    T_post0Cond = T( (T.rotated_post0_x>=pst0x_d) & ...
        (T.rotated_post0_x<pst0x_u), : );

    % Compute #trials in rotated_post0_x condition:
    [~,IA,~] = unique( [T_post0Cond.uuid, T_post0Cond.nStimuli, ...
        T_post0Cond.event], 'rows' );
    nTrialsPost0 = numel( IA );
    
    % Create empty cell and nan vector for storing x,y and decision(+-1),
    % resp. (note: x,y <--> y,x for TPI.m):
    traj_cell = cell(nTrialsPost0,1); 
    dec_vect = nan(nTrialsPost0,1);
    goodTrialsLocs = nan(nTrialsPost0,1);

    % Initialize trial counter:
    c = 0;

    % Compute unique flies' ID:
    uniFlies = unique(T_post0Cond.uuid);

    % Store end locations:
    trEndLocs = [];

    for f = 1:length(uniFlies)
        flyID = uniFlies(f);
    
        % Read fly data in postCond: 
        TPost_f = T_post0Cond( T_post0Cond.uuid==flyID, : );
        
        for s = 1:length(uniNStimuli)
            stim = uniNStimuli(s);
    
            % Read fly data in postCond x stim: 
            TPost_f_s = TPost_f( TPost_f.nStimuli == stim, : );
            
            uniEvent_f_s = unique(TPost_f_s.event);
    
            for e = 1:length(uniEvent_f_s)
                event = uniEvent_f_s(e);

                % Read fly data in postCond x stim x event: 
                TPost_f_s_e = TPost_f_s( TPost_f_s.event == event, :);
                
                % Read x.y in trial
                x = TPost_f_s_e.rotated_x;
                y = TPost_f_s_e.rotated_y;
                
                % Compute corrected x, y:
                xCutLoc = 1 + max( [0, find( x<0, 1, 'last' )] );
                xCorrected = x; 
                yCorrected = y; 

                % ADD TO COUNTER:
                c = c + 1;

                % Store x,y <--> y,x in cell:
                traj_cell{c} = [yCorrected, xCorrected];

                % Compute and store decision:
                dec_vect(c) = sign( TPost_f_s_e.rotated_y(end) );

                % Compute relative time:
                trEndLocs = [trEndLocs; [TPost_f_s_e.rotated_x(end), ...
                    TPost_f_s_e.rotated_y(end)]];

                % Plot in black IFF not going below x<0:

                if sum(xCorrected<0)>=1
                    colTraj = [.5,.5,.5];
                else
                    colTraj = [0,0,0];
                    goodTrialsLocs(c) = 1;
                end
                    
                % Plot path in post0x cond. figure:
                plot( xCorrected, yCorrected, 'Color', colTraj );
                hold on;
                plot( TPost_f_s_e.rotated_x(end), ...
                    TPost_f_s_e.rotated_y(end), 'r.', ...
                    TPost_f_s_e.rotated_x(1), ...
                    TPost_f_s_e.rotated_y(1), 'b.' );
                hold on;

            end
        end
    end 
    
    % Add estimates for targets locations to the traj. plot:
    for tar = 1:2
        % Read end locations for which y>0 or y<0:
        tarYsign = 2*tar-3;
        locsTarg = ( sign(trEndLocs(:,2))==tarYsign );
        tarX = trEndLocs(locsTarg,1);
        tarY = trEndLocs(locsTarg,2);
        % Fit circle:
        [xc,yc,Re,a] = circfit(tarX,tarY);

        % Compute circle:
        th=0:0.01:2*pi; 
        xe = Re*cos(th)+xc; 
        ye = Re*sin(th)+yc;
        pgon = polyshape( xe, ye );
        % Store target's center locations:
        targetsLocs = [targetsLocs; xc,yc];
        targetsRe = [targetsRe; Re];
        % Plot and compute angle between ini. loc. (0,0) and targets' 
        % center:
        plot( pgon, 'EdgeColor', 'none', 'FaceColor', colTar{r}, ...
            'FaceAlpha', .5 );
        hold on;
        % Compute angle between target and ini. loc.:
        if tar==2
            v_1 = [targetsLocs(1+2*(r-1),1),targetsLocs(1+2*(r-1),2),0];
            v_2 = [targetsLocs(2+2*(r-1),1),targetsLocs(2+2*(r-1),2),0];
            Theta = atan2(norm(cross(v_1, v_2)), dot(v_1, v_2));
            estAngle = rad2deg(Theta);
            targetsAngles = [targetsAngles; estAngle];
        end      
        
    end
        
    % Store trajectories and decisions in structure:
    traj_strct.(post0CondName) = traj_cell;
    dec_strct.(post0CondName) = dec_vect;
    trEndLocs_strct.(post0CondName) = trEndLocs;

    % Axes for paths figure:
    xlabel('y (rotated-x)'); 
    ylabel('x (rotated-y)'); 
    grid on;
    title([num2str(pst0x_d) ' <= rotated-post0-x < ' num2str(pst0x_u) ...
        ' (n_{trials}= ' num2str(numel(dec_vect)) ')']);
    ggg = gca;
    axis equal;
    
    % Compute spatial TPI over all trials with y>=0 trials (or all trials):
    if strcmp(TPI_trial_type,'trialsWithRotX>=0')
        % Compute TPI(y) only with trajectories of trials with rot-x>=0:
        dec_vect_good = dec_vect(goodTrialsLocs == 1);
        traj_cell_good = traj_cell(goodTrialsLocs == 1);
        [TPI_flies,TPI_SEM_flies,~,~,TPI_flies_inf,TPI_SEM_flies_inf] = ...
            TPI1( traj_cell_good, yEdges, dec_vect_good, 'y' );
    elseif strcmp(TPI_trial_type,'allTrials')
        % Compute TPI(y) with trajectories of all trials:
        dec_vect_good = dec_vect;
        traj_cell_good = traj_cell;
        [TPI_flies,TPI_SEM_flies,~,~,TPI_flies_inf,TPI_SEM_flies_inf] = ...
            TPI1( traj_cell, yEdges, dec_vect, 'y' );
    end
    
    % Plot TPI(y):
    if r ~=1
        subplot(2,numel(uniNStimuli),r+numel(uniNStimuli));
        errorbar( yCenters, TPI_flies, TPI_SEM_flies, 'k', 'capSize', 0, ...
            'LineWidth', 1 );
        hold on;
        % Also plot TPI(inf):
        errorbar( [-infVal,infVal], TPI_flies_inf, ...
            TPI_SEM_flies_inf, 'ko', 'capSize', 0, 'LineWidth', 1, ...
            'LineStyle', 'none' );
        % Axes for TPI(y):
        xlabel('y (rotated-x)');
        ylabel('TPI(y)');
        ylim([-1,1])
        grid on;
        %axis square;
        title([num2str(pst0x_d) ' <= rotated-post0-x < ' ...
            num2str(pst0x_u) ' (n_{trials}= ' num2str(numel( ...
            dec_vect_good)) ')']);
        % Use the same x-lims as in the figure above:
        xlim( trajSubplot.XLim);
    end
    

end
   

% Also plot all targets:

subplot(2,numel(uniNStimuli),1+numel(uniNStimuli));
for r = 1:length(rangePost0X)-1
    estAngle = targetsAngles(r);
    for tar = 1:2
        % Compute circle:
        targetLocs = targetsLocs(tar+2*(r-1),:);
        xc = targetLocs(1); 
        yc = targetLocs(2); 
        Re = targetsRe(tar+2*(r-1));
        th=0:0.01:2*pi; 
        xe = Re*cos(th)+xc; 
        ye = Re*sin(th)+yc;
        pgon = polyshape( xe, ye );
        % Compute and plot arena edges: 
        if r == 1 && tar == 1
            mRadiusTarget = mean( targetsRe );
            mDistToTarCenter = mean( sqrt( targetsLocs(:,1).^2 + ...
                targetsLocs(:,2).^2 ) );
            mRadiusArena = mRadiusTarget + mDistToTarCenter;
            xe_arena = mRadiusArena*cos(th); 
            ye_arena = mRadiusArena*sin(th);
            pgon_arena = polyshape( ...
                {xe_arena, mRadiusArena * [-1,1,1,-1]}, ...
                {ye_arena, mRadiusArena * [-1,-1,1,1]} );
            plot( pgon_arena, 'EdgeColor', 'none', 'FaceColor', ...
                [.5,.5,.5] );
            hold on;
        end
        % Plot circle:
        if tar == 2
            forPlot(r) = plot( pgon, 'EdgeColor', 'none', 'FaceColor', ...
                colTar{r}, 'FaceAlpha', .5 );
            forLeg{r} = ['angle ~ ' num2str(round(estAngle))];
        else
            plot( pgon, 'EdgeColor', 'none', 'FaceColor', ...
                colTar{r}, 'FaceAlpha', .5 );
        end
        hold on;
        % Plot circle-to-ini. loc.:
        plot( [0, xc], [0,yc], 'Color', colTar{r}, 'lineWidth', 1 );
        hold on;
    end
    
end

% Axes def's:
legend( forPlot, forLeg);
axis equal; grid on;
xlabel('y (rotated-x)');  
ylabel('x (rotated-y)'); 
xlim(mRadiusArena*[-1,1]);
ylim(mRadiusArena*[-1,1]);
xticks(-5:5);
yticks(-5:5);

% Plot arena edges:
for r = 1:3
    trajSubplot = subplot(2,numel(uniNStimuli),r);
    plot( pgon_arena, 'EdgeColor', 'none', 'FaceColor', ...
        [.5,.5,.5] );
    hold on;
    xlim(mRadiusArena*[-1,1]);
    ylim(mRadiusArena*[-1,1]);
    xticks(-5:5);
    yticks(-5:5);
end

% Add ini. loc:
for s = 5:6
    subplot(2,numel(uniNStimuli),s);
    hold on; 
    plot( [0,0], [-1,1], 'b--' );
    xlim([0,Inf]);
    yticks(-1:.2:1);

    allTicks = -infVal:infVal;
    xticks(allTicks);
    midTicks = allTicks(1:end-1);
    midLabels = arrayfun(@num2str, midTicks, 'UniformOutput', false);
    lastLabel = 'y_{\infty}';
    xticklabels([midLabels,lastLabel]);
end
for s = 2:3
    subplot(2,numel(uniNStimuli),s);
    plot( [0,0], mRadiusArena*[-1,1], 'b--' );
end

% Save figure:
saveas( fig_trajTPI, ...
    'analyses/figures/GODM_trajectoriesAndTpis_byTarAngles.fig' );


