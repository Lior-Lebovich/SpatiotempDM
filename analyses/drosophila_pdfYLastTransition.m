% Y PROBABILITY DENSITY FUNCTION OF LAST MIDLINE-CROSSINGS: 
%
% Plots PDF(y|last midline-crossing) for all fly datasets.
%
% INPUTS: 
%   output_drosophila_main_[DATASET_NAME(10:END)].mat
%
% EXTERNAL FUNCTION:
%   This code uses an external function, patchline.m, Reference: 
%   Brett Shoelson (2023). Patchline 
%   (https://www.mathworks.com/matlabcentral/fileexchange/36953-patchline), 
%   MATLAB Central File Exchange. Retrieved September 18, 2023.
%
% OUTPUTS:
%   flies_PdfLastTrans.mat: last midline-crossing Probability Density
%       Functions.
%
% Figures:
%   Avg. PDF(y|last crossing) overlaid on individual curves in dataset:
%       Fig. S1C: flies_pdfLastTrans_Short.fig
%   Datasets comparispon of PDF(y|last crossing) over ALL trials in
%       dataset. Note that this is not the average PDF, but computes the
%       PDF over all trials in the dataset:
%       Figs. 1F, 4A, S3C (modify code for S3B), S13A: 
%       fliesAll_[COMPNAME]Comp_pdfLastTrans.fig
%   Datasets comparispon of PDF(T|last crossing) as above:
%       Fig. 1G: fliesAll_[COMPNAME]Comp_pdfLastTransT.fig
%   Region-aligned short/long comparispon of PDF(y|last crossing) as 
%   above:
%       Fig. 3B: fliesAll_LongAlignComp_pdfLastTrans.fig
%
% Note in the current midline-crossings are termed "x transitions".
% Note that for the long dataset the PDF(y|last x crossing) uses y 
% normalized values (|y|=1 is upper edge of bottom arm, before 
% intersection) rather than region-aligned. For region-alignment - see the 
% last section below. 
% For comp. in absolute size --> un-comments the 7 lines following: 
% "% For abs. long/short comparison:".
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



%% Define datasets for which the PDFs will be computed:

dataTypes_stable = {'Turndata_ShortYy', 'Turndata_LongYy', ...
    'Turndata_ShortNorpAYy', 'Turndata_ShortDumbYy', ...
    'Turndata_ShortFoxPYy', 'Turndata_ShortNompCYy', ...
    'Turndata_ShortABM_NoWF', ...
    'Turndata_ShortABM_WF_00100', 'Turndata_ShortABM_WF_00200', ...
    'Turndata_ShortABM_WF_00500', 'Turndata_ShortABM_WF_00700', ...
    'Turndata_ShortABM_WF_00900', ...
    'Turndata_ShortABM_WallFollowing', 'Turndata_ShortABM_brown5'};
dataCols_stable = {'k', 'b', 'r', 'm', 'g', '#77AC30', ... % fly data
    "#A2142F", ... % no WF 
    "#D95319", "#EDB120", "#77AC30", "#4DBEEE", "#0072BD", ... % WF parm.
    "#7E2F8E", "#7E2F8E" }; % WF, brownian

dataTypes = dataTypes_stable; 
dataCols = dataCols_stable;

kSubSamp = 1; % subsampling frames (1= no subsamling)

regionPrintNames = {'Pre-arm', 'Down-the-arm', 'Cul-de-sac', ...
    'Up-the-arm', 'Post-arm'};

infVal_global = 1.5;
fps = 30;



%% Compute and plot PDF for each fly and dataset:
% The dataset avg. PDF, PDF over all trials in dataset and corresp. 
% yCenters will be stored in structures and used for datasets comparisons
% (following section).


for d = 1:length(dataTypes)

    d_stable = d;
    dataType = dataTypes{d};

    % Load processed data:
    matFileName = ['output_drosophila_main_' dataType(10:end) '.mat'];
    load(matFileName);

    % Keep dataType and Cols despite val's in matFileName.m 
    d = d_stable;
    dataTypes = dataTypes_stable; 
    dataCols = dataCols_stable;
    dataType = dataTypes{d};


    % Add bins to yEdges:
    yEdges = [yEdges, yEdges(end) + (1:3)*diff(yEdges(1:2))];
    yCenters = yEdges(1:end-1) + .5*diff(yEdges(1:2));

    if strcmp(dataType,'Turndata_LongYy')
        culDeSacY = 0.68/3.36; 
        yEdges_final = yEdges;
        yCenters_final = yCenters;
        infVal = 1.5; 
        regStartLoc = 1;
        % For abs. long/short comparison:
        %{
        rel_yEdges = yEdges;
        delta_rel_yEdges = diff(yEdges(1:2));
        abs_yEdges_pos = delta_rel_yEdges:delta_rel_yEdges:2.2;
        abs_yEdges = [-1*flip(abs_yEdges_pos), 0, abs_yEdges_pos];
        yEdges_final = abs_yEdges;
        yEdges = abs_yEdges * ratioSize;
        yCenters_final = yEdges_final(1:end-1) + .5 * diff(yEdges_final(1:2));
        %}
    elseif strcmp(dataType,'Turndata_ShortYy') || ...
            strcmp(dataType,'Turndata_ShortNorpAYy') || ...
            strcmp(dataType,'Turndata_ShortDumbYy') || ...
            strcmp(dataType,'Turndata_ShortFoxPYy') || ...
            strcmp(dataType,'Turndata_ShortNompCYy') 
        culDeSacY = 0.69/2.03;
        yEdges_final = yEdges;
        yCenters_final = yCenters;
        infVal = 1.5; 
        regStartLoc = 1;
    elseif strcmp(dataType,'Turndata_ShortABM_NoWF') || ...
            strcmp(dataType,'Turndata_ShortABM_WallFollowing') || ...
            strcmp(dataType,'Turndata_ShortABM_Brownian') || ...
            strcmp(dataType,'Turndata_ShortABM_WF_00100') || ...
            strcmp(dataType,'Turndata_ShortABM_WF_00200') || ...
            strcmp(dataType,'Turndata_ShortABM_WF_00500') || ...
            strcmp(dataType,'Turndata_ShortABM_WF_00700') || ...
            strcmp(dataType,'Turndata_ShortABM_WF_00900') || ...
            strcmp(dataType,'Turndata_ShortABM_brown5')
        culDeSacY = 7.5/30;
        locNegCul = find(yEdges<=-(.2 + culDeSacY),1,'last');
        yEdges_final = yEdges(locNegCul:end);
        yEdges = yEdges_final;
        yCenters_final = yCenters(locNegCul:end);
        infVal = .2 + culDeSacY;
        regStartLoc = 3;
    end

    % Define region breaks (for plotting):
    regionsBreak = [-1.3, -1, -culDeSacY, culDeSacY, 1, 1.3];

    % Create arrays for storing pdf's and Y,T value of last tranitions:
    lastYOfXSwitch_cell = cell(numel(selectedFlies),1);
    lastTOfXSwitch_cell = cell(numel(selectedFlies),1);
    pdfYLastTrans = nan( numel(selectedFlies), numel(yEdges)-1 );

    % Create figure for all flies in dataset:
    figDataPdfYTrans = figure;
    ax1 = axes();
    box(ax1);
    
    % Compute and plot indiv's pdf(y|last x transition):

    for nF = 1:length(selectedFlies) 

        % Read fly's data:
        f = selectedFlies(nF);
        fly_xytyd = meanX_vects_ALL.(['f' num2str(f)]).xytyd;
        fly_locsBotTrialT = meanX_vects_ALL.(['f' num2str(f)]).locsBotTrialT;
        fly_Dec = meanX_vects_ALL.(['f' num2str(f)]).Dec;
        

        % Read fly's botttom trials:
        fly_xytyd_good = fly_xytyd( fly_locsBotTrialT );
        fly_dec_good = fly_Dec( fly_locsBotTrialT );

        % Create cells for storing x,y coordinates:
        flyTrlsX = cell(size(fly_xytyd_good,1),1);
        flyTrlsY = cell(size(fly_xytyd_good,1),1);
        flyTrlsT = cell(size(fly_xytyd_good,1),1);
        flyTrlsID = cell(size(fly_xytyd_good,1),1);
        flyTrlsDec = cell(size(fly_xytyd_good,1),1);

        flyTrialsConsidered = zeros( numel(fly_xytyd_good), 1);

        for t = 1:numel(fly_xytyd_good)
    
            % read fly's x,y in trial (downsampling):
            fly_xytyd_good_trial = fly_xytyd_good{t};
            fly_xGood = fly_xytyd_good_trial(1:kSubSamp:end,1);
            fly_yGood = fly_xytyd_good_trial(1:kSubSamp:end,2);
            fly_tGood = fly_xytyd_good_trial(1:kSubSamp:end,3);
            fly_decGood = fly_dec_good(t);
    
            % Only use trials with full data:
            if min(fly_yGood)<=-1 && max(fly_yGood)>=1
                flyTrialsConsidered(t) = 1; 
                flyTrlsX{t} = fly_xGood;
                flyTrlsY{t} = fly_yGood;
                flyTrlsT{t} = fly_tGood;
                flyTrlsID{t} = t * ones(size(fly_yGood));
                flyTrlsDec{t} = fly_decGood * ones(size(fly_yGood));
            end
        end


        % Centering x (x=0 corresp. horizontal midline):

        % First, find the fly's min. and max. x values for each y-bin:
        flyMinXInBin = nan(length(yEdges)-1,1);
        flyMaxXInBin = nan(length(yEdges)-1,1);
        flyCentXInBin = nan(length(yEdges)-1,1);
        flyDiffXInBin = nan(length(yEdges)-1,1);
        trials_xytIdDecMat = [cell2mat( flyTrlsX(flyTrialsConsidered==1) ), ...
            cell2mat( flyTrlsY(flyTrialsConsidered==1) ), ...
            cell2mat( flyTrlsT(flyTrialsConsidered==1) ), ...
            cell2mat( flyTrlsID(flyTrialsConsidered==1) ), ...
            cell2mat( flyTrlsDec(flyTrialsConsidered==1) )];
        % Then, compute midline:
        yyy = trials_xytIdDecMat(:,2);
        xinArm = trials_xytIdDecMat( abs(yyy)<=1 & abs(yyy)>culDeSacY,1);
        minmax_xinArm = [min(xinArm), max(xinArm)];
        forCentX = mean(minmax_xinArm);
        % Store centered x:
        trials_xytIdDecXNormMat = [trials_xytIdDecMat, ...
            trials_xytIdDecMat(:,1) - forCentX];
        
        
        % Now, use the normX * decision for each trial (this will be used
        % as an estimate for the predictability: for meaduring the last 
        % "transition" y-location in each trial: 
        % Also compute the rate of transitions:

        uniTrialNums = unique(trials_xytIdDecXNormMat(:,4));
        lastYOfXSwitch = nan(size(uniTrialNums));
        lastTOfXSwitch = nan(size(uniTrialNums));
        numTransArm = nan(size(uniTrialNums));
        % Create matrices for storing #trans. and #frames for trans. rate:
        flyNumTrans_yBin = nan(length(uniTrialNums), length(yEdges)-1);
        flyNumFrames_yBin = nan(length(uniTrialNums), length(yEdges)-1);
        flyNumTrans_tBin = nan(length(uniTrialNums), length(tEdges)-1);
        flyNumFrames_tBin = nan(length(uniTrialNums), length(tEdges)-1);

        for tt = 1:length(uniTrialNums)
            trialNum = uniTrialNums(tt);
            locTrial = ( trials_xytIdDecXNormMat(:,4) == trialNum );
            decsTr = trials_xytIdDecXNormMat(locTrial,5);
            decTr = decsTr(1);
            yTr = trials_xytIdDecXNormMat(locTrial,2);
            tTr = trials_xytIdDecXNormMat(locTrial,3);
            xNormTr = trials_xytIdDecXNormMat(locTrial,end);
            xNormTrMultDec = xNormTr * decTr;

            % Compute last y in which x-transition:
            locLastNegX = find( xNormTrMultDec < 0, 1, 'last' );
            lastYOfXSwitch(tt) =  mean( yTr( ...
                max([1,locLastNegX]) + (0:1) ) );
            lastTOfXSwitch(tt) =  mean( tTr( ...
                max([1,locLastNegX]) + (0:1) ) );

            % Also compute the rate of transitions by bin 
            % and #transitions WITHIN the bottom arm:
            isTransition = ( [xNormTrMultDec(1:end-1) .* ...
                xNormTrMultDec(2:end); 0] < 0 );
            isTransitionBinary = (isTransition == 1);
            armLocs = (abs(yTr)<=1);
            numTransArm(tt) = sum(isTransitionBinary);
            % Rate over y-bins:
            for dyPos = 1:(length(yEdges)-1)
                locsYp = (yTr >= yEdges(dyPos)) & (yTr < yEdges(dyPos+1)); 
                flyNumTrans_yBin(tt,dyPos) = sum(isTransitionBinary(locsYp));
                flyNumFrames_yBin(tt,dyPos) =  sum(locsYp);
            end
            % Rate over relative time-bins:
            for dtPos = 1:(length(tEdges)-1)
                locsTp = (tTr >= tEdges(dtPos)) & (tTr < tEdges(dtPos+1)); 
                flyNumTrans_tBin(tt,dtPos) = sum(isTransitionBinary(locsTp));
                flyNumFrames_tBin(tt,dtPos) =  sum(locsTp);
            end

        end

        % Store fly's Y's of last x-switch in cell array:
        lastYOfXSwitch_cell{nF} = lastYOfXSwitch;
        lastTOfXSwitch_cell{nF} = lastTOfXSwitch;

        % Compute + plot pdf(y| last x switch in trial):
        figure(figDataPdfYTrans);
        N = histcounts( lastYOfXSwitch, 'Normalization', ...
            'pdf', 'binEdges', yEdges );
        plot( yCenters_final, N );
        hold on;

        % Store pdf(y| last x switch in trial):
        pdfYLastTrans(nF,:) = N;

    end

    % Plot average pdf(y| last x switch in trial):
    figure(figDataPdfYTrans);
    if strcmp(dataType,'Turndata_ShortABM_brown5')
        avgCol = [.5, .3, .1];
    else
        avgCol = dataCols{d};
    end
    errorbar( yCenters_final, mean(pdfYLastTrans), ...
        std(pdfYLastTrans)/sqrt(numel(selectedFlies)), 'lineWidth', 1, ...
        'Color', avgCol );
    xlabel('y'); 
    ylabel('pdf(y | last crossing)');
    currentYLim = get(ax1,'YLim');
    xlim([-infVal, infVal_global]);
    xticks(-1:.5:1);
    ylim([0,Inf]);
    
    % Plot region-breaks:
    for rr = 1:numel(regionsBreak)
        plot(regionsBreak(rr) * [1,1], [0,currentYLim(2)], 'k:' ); 
        hold on;
    end

    % Add region names:
    ax2 = axes('Position', get(ax1,'Position'), ...
        'XAxisLocation', 'top', 'Color', 'none', 'XColor', 'k' );
    ax2.YAxis.Visible = 'off';
    ax2.XLim = ax1.XLim;
    ax2.XTick = regionsBreak(regStartLoc:end-1) + ...
        .5*diff(regionsBreak(regStartLoc:end));
    ax2.XTickLabel = regionPrintNames(regStartLoc:end);
    ax2.FontSize = 9;

    % Save indiv' figure:
    figStartName = ['analyses/figures/flies_pdfLastTrans_' ...
        dataType(10:end-2)];
    saveas( figDataPdfYTrans, [figStartName '.fig']);

    % Compute and store the avg. and sem PDF over all flies in dataset:
    savedPDFs.(dataType).meanPDF = mean(pdfYLastTrans);
    savedPDFs.(dataType).semPDF = std(pdfYLastTrans) / ...
        sqrt(numel(selectedFlies));

    % Compute and store pdf(y| last x switch in trial) over TRIALS OF ALL 
    % FLIES in dataset:
    lastYOfXSwitch_all = cell2mat(lastYOfXSwitch_cell);
    N = histcounts( lastYOfXSwitch_all, 'Normalization', 'pdf', ...
        'binEdges', yEdges );
    savedPDFs.(dataType).allTrialsPDF = N;
    savedPDFs.(dataType).yCenters = yCenters_final;
    savedPDFs.(dataType).regions = regionsBreak(regStartLoc:end);
    savedPDFs.(dataType).regionsNames = regionPrintNames(regStartLoc:end);
    savedPDFs.(dataType).xLims = [-infVal, infVal_global];

    % Compute and store pdf(T| last x switch in trial) over TRIALS OF ALL 
    % FLIES in dataset:
    lastTOfXSwitch_all = cell2mat(lastTOfXSwitch_cell);
    N = histcounts( lastTOfXSwitch_all, 'Normalization', 'pdf', ...
        'binEdges', tEdges );
    savedPDFs.(dataType).allTrialsPDF_T = N;
    savedPDFs.(dataType).tCenters = tCenters;

end


% Save dataset PDFs:
save( 'datafiles/flies_PdfLastTrans.mat', 'savedPDFs' );



%% Plot datasets comparison of PDF(y|last x transition):

load( 'flies_PdfLastTrans.mat', 'savedPDFs' );

% Determines which dataset PDF to plot. Change below pdfGlobalType value to 
% 'meanPDF' to plot the avg. PDF over flies +-sem.
pdfGlobalType = 'allTrialsPDF';


% Note that short/long comp. plots long in norm. rather than 
% region-aligned. For region-alignment - see next section. For comp. in 
% abs. size the section after.

% Define dataset comparisons:
compNames = {'Short', 'Models', 'LongNorm'};
compLoc.Short = [1, 3:6]; % short fly comparison
compLoc.Models = 7:14; % model comparison
compLoc.LongNorm = 1:2; % short/long Norm. comparison

for c = 1:numel(compNames)
    compName = compNames{c};
    compLocs = compLoc.(compName);
    dataTypes_c = dataTypes_stable(compLocs);
    dataCols_c = dataCols_stable(compLocs);

    % compute maximal Y value: 
    maxY = 0;
    for dd = 1:numel(dataTypes_c)
        maxY = max([maxY, max(savedPDFs.(dataTypes_c{dd}).(pdfGlobalType))]);
    end

    % Create dataset PDF figure:
    figAvgPdfYTrans = figure;
    ax1 = axes();
    box(ax1);          

    % Initialize legends:
    clear legs legNames;
    
    % Plot PDF(y|last transition) for each dataType:
    
    for d = 1:numel(dataTypes_c)
        dataType = dataTypes_c{d};
        dataCol = dataCols_c{d};
        
        % Plot PDF(y| last x switch in trial):

        figure(figAvgPdfYTrans);
        yCenters_final = savedPDFs.(dataType).yCenters;

        if strcmp(dataType,'Turndata_ShortABM_brown5')
            avgCol = [.5, .3, .1];
        else
            avgCol = dataCol;
        end

        if strcmp(pdfGlobalType,'allTrialsPDF')
            allTrialsPDF = savedPDFs.(dataType).allTrialsPDF;
            legs(d) = plot( yCenters_final, allTrialsPDF, 'lineWidth', ...
                1, 'Color', avgCol, 'Marker', '.' );
        elseif strcmp(pdfGlobalType,'meanPDF')
            meanPDF = savedPDFs.(dataType).meanPDF;
            semPDF = savedPDFs.(dataType).semPDF;
            legs(d) = errorbar( yCenters_final, meanPDF, semPDF, ...
                'lineWidth', 1, 'Color', avgCol );
        end
        hold on;
        legNames{d} = dataType(10:end-2);


        % Load region breaks:
        regionsBreak = savedPDFs.(dataType).regions;
        regionPrintNames = savedPDFs.(dataType).regionsNames;
        xLims = savedPDFs.(dataType).xLims;

        % Plot region-breaks:
        if d==1
            for rr = 1:numel(regionsBreak)
                plot(regionsBreak(rr) * [1,1], [0,maxY], 'k:' ); 
                hold on;
            end
        end
        % Also plot Long region breaks if Long/short comparsion:
        if strcmp(dataType,'Turndata_LongYy') && ...
                strcmp(compName,'LongNorm')
            for rr = 1:numel(regionsBreak)
                plot(regionsBreak(rr) * [1,1], [0,maxY], 'b:' ); 
                hold on;
            end
        end
        
    end
    % axes: 
    xlim(xLims);
    xticks(-1:.5:1);
    xlabel( 'y' ); 
    ylabel('pdf(y | last crossing)');
    legend( legs, legNames, 'Location', 'NorthWest' );
    ylim([0,Inf]);
    % Add region names: y-regions are only relevant for SHORT fly 
    % datasets; Time-regions are relevant for all datasets:
    if d == numel(dataTypes_c)
        ax2 = axes('Position', get(ax1,'Position'), 'XAxisLocation', ...
            'top', 'Color', 'none', 'XColor', 'k' );
        ax2.YAxis.Visible = 'off';
        ax2.XLim = ax1.XLim;
        ax2.XTick = regionsBreak(1:end-1) + .5*diff(regionsBreak);
        ax2.XTickLabel = regionPrintNames;
        ax2.FontSize = 9;
    end
    % axes + save figure: 
    figStartName = ['analyses/figures/fliesAll_' compName ...
        'Comp_pdfLastTrans'];
    saveas( figAvgPdfYTrans, [figStartName '.fig']);
      
end



%% Plot datasets comparison of PDF(T|last x transition):

load( 'flies_PdfLastTrans.mat', 'savedPDFs' );

% Determines which dataset PDF to plot. Change below pdfGlobalType value to 
% 'meanPDF' to plot the avg. PDF over flies +-sem.
pdfGlobalType = 'allTrialsPDF_T';


% Note that short/long comp. plots long in norm. rather than 
% region-aligned. For region-alignment - see next section. For comp. in 
% abs. size the section after.

% Define dataset comparisons:
compNames = {'Short', 'Models', 'LongNorm'};
compLoc.Short = [1, 3:6]; % short fly comparison
compLoc.Models = 7:14; % model comparison
compLoc.LongNorm = 1:2; % short/long Norm. comparison

for c = 1:numel(compNames)
    compName = compNames{c};
    compLocs = compLoc.(compName);
    dataTypes_c = dataTypes_stable(compLocs);
    dataCols_c = dataCols_stable(compLocs);

    % compute maximal T value: 
    maxY = 0;
    for dd = 1:numel(dataTypes_c)
        maxY = max([maxY, max(savedPDFs.(dataTypes_c{dd}).(pdfGlobalType))]);
    end

    % Create dataset PDF figure:
    figAvgPdfYTrans = figure;
    ax1 = axes();
    box(ax1);          

    % Initialize legends:
    clear legs legNames;
    
    % Plot PDF(T|last transition) for each dataType:
    
    for d = 1:numel(dataTypes_c)
        dataType = dataTypes_c{d};
        dataCol = dataCols_c{d};
        
        % Plot PDF(T| last x switch in trial):

        figure(figAvgPdfYTrans);
        tCenters_final = savedPDFs.(dataType).tCenters;

        if strcmp(dataType,'Turndata_ShortABM_brown5')
            avgCol = [.5, .3, .1];
        else
            avgCol = dataCol;
        end

        allTrialsPDF = savedPDFs.(dataType).(pdfGlobalType);
        legs(d) = plot( tCenters_final, allTrialsPDF, 'lineWidth', ...
            1, 'Color', avgCol, 'Marker', '.' );

        hold on;
        legNames{d} = dataType(10:end-2);


        % Load region breaks:
        regionsBreak = [-1.3, -1, 0, 1, 1.3];
        regionPrintNames = {'Pre-arm', 'Down-the-arm & c.d.s<0', ...
            'C.d.s>0 up-the-arm' ,'post-arm'};
        xLims = savedPDFs.(dataType).xLims;

        % Plot region-breaks:
        if d==1
            for rr = 1:numel(regionsBreak)
                plot(regionsBreak(rr) * [1,1], [0,maxY], 'k:' ); 
                hold on;
            end
        end
        
    end
    % axes: 
    xlim(xLims);
    xticks(-1:.5:1);
    xlabel( 'T' ); 
    ylabel('pdf(T | last crossing)');
    legend( legs, legNames, 'Location', 'NorthWest'  );
    ylim([0,Inf]);
    % Add region names: y-regions are only relevant for SHORT fly 
    % datasets; Time-regions are relevant for all datasets:
    if d == numel(dataTypes_c)
        ax2 = axes('Position', get(ax1,'Position'), 'XAxisLocation', ...
            'top', 'Color', 'none', 'XColor', 'k' );
        ax2.YAxis.Visible = 'off';
        ax2.XLim = ax1.XLim;
        ax2.XTick = regionsBreak(1:end-1) + .5*diff(regionsBreak);
        ax2.XTickLabel = regionPrintNames;
        ax2.FontSize = 9;
    end
    % axes + save figure: 
    figStartName = ['analyses/figures/fliesAll_' compName ...
        'Comp_pdfLastTransT'];
    saveas( figAvgPdfYTrans, [figStartName '.fig']);
      
   % clear legs legNames;
end



%% Short vs Long PDF(y|last x transition) REGION-ALIGNED

% Determines which dataset PDF to plot. Change below pdfGlobalType value to 
% 'meanPDF' to plot the avg. PDF over flies +-sem.
pdfGlobalType = 'allTrialsPDF';

% Define datasets for which the PDFs will be computed:
dataTypes_stable = {'Turndata_ShortYy', 'Turndata_LongYy'};
dataCols_stable = {'k', 'b'}; 
dataTypes = dataTypes_stable; 
dataCols = dataCols_stable;
kSubSamp = 1; % subsampling frames (1= no subsamling)

% Read ratio size in long maze and cul-de-sac loc. of short:
ratioForLong = load( 'output_drosophila_main_LongYy.mat', 'ratioSize' );
ratioLong = ratioForLong.ratioSize;
culDeSacY = 0.69/2.03;
% Define short sub-edges:
edgesy.Short.beforeBottom = -1.5:.1:-1;
edgesy.Short.downArm = linspace(-1,-culDeSacY,7);
edgesy.Short.culdesac = linspace(-culDeSacY,culDeSacY,7);
edgesy.Short.upArm = linspace(culDeSacY,1,7);
edgesy.Short.afterBottom = 1:.1:1.5;
% Define long sub-edges:
edgesy.Long.beforeBottom = -(1+(0.5*ratioLong)):(.1*ratioLong):-1;
edgesy.Long.downArm = linspace(-1,-culDeSacY*ratioLong,15);
edgesy.Long.culdesac = ratioLong * linspace(-culDeSacY,culDeSacY,7);
edgesy.Long.upArm = linspace(culDeSacY*ratioLong,1,15);
edgesy.Long.afterBottom = 1:(.1*ratioLong):(1+(0.5*ratioLong));

% Define short and long full edges and centers:
regionNames = {'beforeBottom', 'downArm', 'culdesac', 'upArm', ...
    'afterBottom'};
regionsBreak = [-1.3, -1, -culDeSacY, culDeSacY, 1, 1.3];


for d = 1:numel(dataTypes)
    dataType = dataTypes{d};
    dataType0 = dataType(10:end-2);
    edgesy_full = [];
    centersy_full = [];
    centersFinaly_full = [];
    if strcmp(dataType,'Turndata_LongYy')
        yEdgesCrrtn = 1 / ratioLong;
    else
        yEdgesCrrtn = 1;
    end
    for r = 1:numel(regionNames)
        regionName = regionNames{r};
        currentEdge = edgesy.(dataType0).(regionName);
        currentCenter = currentEdge(1:end-1) + ...
            (.5 * diff(currentEdge(1:2)));
        if r == 1
            addEdges = currentEdge;
        else
            addEdges = currentEdge(2:end);
        end
        edgesy_full = [edgesy_full, addEdges];
        centersy_full = [centersy_full, currentCenter];
        if strcmp(dataType,'Turndata_LongYy') && ( ...
                strcmp(regionName,'beforeBottom') || ...
                strcmp(regionName,'afterBottom') )
            theShort = edgesy.Short.(regionName);
            yCentersFinal = theShort(1:end-1) + .5 * diff(theShort(1:2));
        elseif strcmp(dataType,'Turndata_LongYy') && ( ...
                strcmp(regionName,'downArm') || ...
                strcmp(regionName,'upArm') )
            theShort = edgesy.Short.(regionName);
            theShort2 = linspace( theShort(1), theShort(end), numel(currentEdge) );
            yCentersFinal = theShort2(1:end-1) + .5 * diff(theShort2(1:2));
        else
            yCentersFinal = yEdgesCrrtn * ( currentEdge(1:end-1) + ...
                .5* diff(currentEdge(1:2)) );
        end
        centersFinaly_full = [centersFinaly_full, yCentersFinal];
    end
    edgesy.(dataType0).All = edgesy_full;
    centersy.(dataType0).All = centersy_full;
    centersFinaly.(dataType0).All = centersFinaly_full;
end


% Compute PDF(y|last x transition) for each dataset:

figAvgPdfYTrans = figure;
ax1 = axes();
box(ax1);
hold on;

for d = 1:numel(dataTypes)

    d_stable = d;
    dataType = dataTypes{d};
    dataType0 = dataType(10:end-2);
    % Load dataset:
    matFileName = ['output_drosophila_main_' dataType(10:end) '.mat'];
    load(matFileName);
    % Keep dataType and Cols despite val's in matFileName.m.
    d = d_stable;
    dataTypes = dataTypes_stable; 
    dataCols = dataCols_stable;
    dataType = dataTypes{d};

    % Read cul-de-sac in data (for mid loc. est.):
    if strcmp(dataType,'Turndata_LongYy')
        culDeSacY = 0.68/3.36; 
    elseif strcmp(dataType,'Turndata_ShortYy')
        culDeSacY = 0.69/2.03;
    end

    % Read y-bins:
    yEdges = edgesy.(dataType0).All;
    yCenters = centersy.(dataType0).All;
    yFinalCenters = centersFinaly.(dataType0).All;

    lastYOfXSwitch_cell = cell(numel(selectedFlies),1);

    pdfYLastTrans = nan( numel(selectedFlies), numel(yEdges)-1 );

    for nF = 1:length(selectedFlies) 

        % Read fly's data:
        f = selectedFlies(nF);
        fly_xytyd = meanX_vects_ALL.(['f' num2str(f)]).xytyd;
        fly_locsBotTrialT = meanX_vects_ALL.(['f' num2str(f)]).locsBotTrialT;
        fly_Dec = meanX_vects_ALL.(['f' num2str(f)]).Dec;
        

        % Read fly's botttom trials:
        fly_xytyd_good = fly_xytyd( fly_locsBotTrialT );
        fly_dec_good = fly_Dec( fly_locsBotTrialT );

        % Create cells for storing x,y coordinates:
        flyTrlsX = cell(size(fly_xytyd_good,1),1);
        flyTrlsY = cell(size(fly_xytyd_good,1),1);
        flyTrlsID = cell(size(fly_xytyd_good,1),1);
        flyTrlsDec = cell(size(fly_xytyd_good,1),1);

        flyTrialsConsidered = zeros( numel(fly_xytyd_good), 1); %%%

        for t = 1:numel(fly_xytyd_good)
    
            % read fly's x,y in trial (downsampling):
            fly_xytyd_good_trial = fly_xytyd_good{t};
            fly_xGood = fly_xytyd_good_trial(1:kSubSamp:end,1);
            fly_yGood = fly_xytyd_good_trial(1:kSubSamp:end,2);
            fly_decGood = fly_dec_good(t);
    
            % Only use trials with full data:
            if min(fly_yGood)<=-1 && max(fly_yGood)>=1
                flyTrialsConsidered(t) = 1; 
                flyTrlsX{t} = fly_xGood;
                flyTrlsY{t} = fly_yGood;
                flyTrlsID{t} = t * ones(size(fly_yGood));
                flyTrlsDec{t} = fly_decGood * ones(size(fly_yGood));
            end
        end


        % Centering x (x=0 corresp. horizontal midline):

        % First, find the fly's min. and max. x values for each y-bin:
        flyMinXInBin = nan(length(yEdges)-1,1);
        flyMaxXInBin = nan(length(yEdges)-1,1);
        flyCentXInBin = nan(length(yEdges)-1,1);
        flyDiffXInBin = nan(length(yEdges)-1,1);
        trials_xytIdDecMat = [cell2mat( flyTrlsX(flyTrialsConsidered==1) ), ...
            cell2mat( flyTrlsY(flyTrialsConsidered==1) ), ...
            cell2mat( flyTrlsID(flyTrialsConsidered==1) ), ...
            cell2mat( flyTrlsDec(flyTrialsConsidered==1) )];
        % Then, compute midline:
        yyy = trials_xytIdDecMat(:,2);
        xinArm = trials_xytIdDecMat( abs(yyy)<=1 & abs(yyy)>culDeSacY,1);
        minmax_xinArm = [min(xinArm), max(xinArm)];
        forCentX = mean(minmax_xinArm);
        % Store centered x:
        trials_xytIdDecXNormMat = [trials_xytIdDecMat, ...
            trials_xytIdDecMat(:,1) - forCentX];
        

        % Now, use the normX * decision for each trial (this will be used
        % as an estimate for the predictability: for meaduring the last 
        % "transition" y-location in each trial: 

        uniTrialNums = unique(trials_xytIdDecXNormMat(:,3));
        lastYOfXSwitch = nan(size(uniTrialNums));
        for tt = 1:length(uniTrialNums)
            trialNum = uniTrialNums(tt);
            locTrial = ( trials_xytIdDecXNormMat(:,3) == trialNum );
            decsTr = trials_xytIdDecXNormMat(locTrial,4);
            decTr = decsTr(1);
            yTr = trials_xytIdDecXNormMat(locTrial,2);
            xNormTr = trials_xytIdDecXNormMat(locTrial,end);
            xNormTrMultDec = xNormTr * decTr;
            % Compute last y in which x-transition:
            locLastNegX = find( xNormTrMultDec < 0, 1, 'last' );
            lastYOfXSwitch(tt) =  mean( yTr( ...
                max([1,locLastNegX]) + (0:1) ) );
        end

        % Store fly's Y's of last x-switch in cell array:
        lastYOfXSwitch_cell{nF} = lastYOfXSwitch;

        % Compute + pdf(y| last x switch in trial):
        N = histcounts( lastYOfXSwitch, 'Normalization', ...
            'pdf', 'binEdges', yEdges );

        % Store pdf(y| last x switch in trial):
        pdfYLastTrans(nF,:) = N;

    end


    % Plot average pdf(y| last x switch in trial):
    figure(figAvgPdfYTrans);
    if strcmp(pdfGlobalType,'allTrialsPDF')
        lastYOfXSwitch_all = cell2mat(lastYOfXSwitch_cell);
        N = histcounts( lastYOfXSwitch_all, 'Normalization', 'pdf', ...
            'binEdges', yEdges );
        legs(d) = plot( yFinalCenters, N, 'lineWidth', 1, 'Color', ...
            dataCols{d}, 'Marker', '.' );
    elseif strcmp(pdfGlobalType,'meanPDF')
        legs(d) = errorbar( yFinalCenters, mean(pdfYLastTrans), ...
            std(pdfYLastTrans)/sqrt(numel(selectedFlies)), 'lineWidth', 1, ...
            'Color', dataCols{d} );
    end
    hold on;
    legNames{d} = dataType(10:end-2);

    if d == numel(dataTypes)
        currentYLim = get(ax1,'YLim');
        % Plot region-breaks:
        for rr = 1:numel(regionsBreak)
            plot(regionsBreak(rr) * [1,1], [0,currentYLim(2)], 'k:' ); 
            hold on;
        end
        xlabel('y'); 
        ylabel('pdf(y | last crossing)');
        xlim([-1.5, 1.5]);
        xticks(-1:.5:1);
        ylim([0,Inf]);
        legend( legs, legNames, 'Location', 'NorthWest' );

         % Add region names:
        ax2 = axes('Position', get(ax1,'Position'), ...
            'XAxisLocation', 'top', 'Color', 'none', 'XColor', 'k' );
        ax2.YAxis.Visible = 'off';
        ax2.XLim = ax1.XLim;
        ax2.XTick = regionsBreak(1:end-1) + .5*diff(regionsBreak);
        ax2.XTickLabel = regionPrintNames;
        ax2.FontSize = 9;
    end

end

% Save figure: 
figStartName = 'analyses/figures/fliesAll_LongAlignComp_pdfLastTrans';
saveas( figAvgPdfYTrans, [figStartName '.fig']);


