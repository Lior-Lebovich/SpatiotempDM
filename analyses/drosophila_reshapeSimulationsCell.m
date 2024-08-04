%% RESTRUCTURING OF ABM SIMULATIONS DATASETSS:
%
% Reshapes each model simulation dataset (cell array) to real flies data 
% structures.
%
% INPUT FILE:
%   rawdata/rawABMData/dataABM_[DATATYPE_NAME].mat.
%   'data' is a cell such that data{fl,fr,p} stores the values of 
%   parameter p, in frame fr, of (simulated) fly fl.
%   Parameters (p) for every frame are: 
%   1: X coordinate; 
%   2: Y coordinate; 
%   3: current Turn number; 
%   4: current Turn direction (0: left, 1: right); 
%   5: current Turn arm start (0: bottom arm, 1: left arm, 2: right arm); 
%   6: current angular velocity; 
%   7: current fly heading decision (change in angular velocity).
% 
% OUTPUT FILE:
%   Turndata_Short[DATATYPE_NAME(5:end)].mat
%   'Turncell_L' is a cell such that Turncell_L{fl,tr,p} stores 
%   the values of parameter p, in trial (turn) tr, for (simulated) fly 
%   fl.
%   Parameters (p) for every trial are:
%   1: x-coordinateS;
%   2: y-coordinateS; 
%   3: start frame #;
%   4: current Turn direction (0: left, 1: right). 



%% Define working directory:

cdName = 'SpatiotempDM';
currentFolder = cd;

if ~endsWith(currentFolder,cdName)
    cd(cdName);
end

addpath( 'analyses', 'datafiles', 'datafiles/rawdata/rawABMData', ...
    'analyses/customfuns', '-frozen');



%% For each datatype: read data and store each fly*trial*param in Turncell_L:


% Run over datatypes:

dataTypes = {'dataABM_NoWF', 'dataABM_WallFollowing', ...
    'dataABM_WF_00100', 'dataABM_WF_00200', 'dataABM_WF_00500', ...
    'dataABM_WF_00700', 'dataABM_WF_00900', 'dataABM_brown5'};

for d = 1:length(dataTypes) 

    dataType = dataTypes{d};
    

    % Read raw data and create an empty cell array in which data will be 
    % stored:
    
    % Load raw data:
    load([dataType '.mat']);
    TurncellRaw = data;
    
    % Compute desired dimensions of output cell array:
    nFlies = size(TurncellRaw,1); % # simulated flies.
    % Compute maximal number of trials:
    flyNumTrials = nan(nFlies,1);
    for f = 1:nFlies
        flyTrials = TurncellRaw(f,:,3);
        flyNumTrials(f) = 1 + flyTrials(end);
    end
    maxNTrial = max(flyNumTrials);
    flyNotMissing = zeros(nFlies,1);
    
    % Create the output (empty) cell array:
    Turncell_L = cell(nFlies, maxNTrial, 4);
    
    
    % Run over flies and store the data of each fly*trial*param in 
    % Turncell_L:
    
    for f = 1:nFlies
        
        % Read entire fly's data:
        xAll = .1 * TurncellRaw(f,:,1)';
        yAll = .1 * TurncellRaw(f,:,2)';
        frameNumAll = (1:length(xAll))';
        turnNumberAll = 1 + TurncellRaw(f,:,3)';
        turnDecisionAll = TurncellRaw(f,:,4)';
        
        % Include only flies without missing datapoints:
        if sum( isnan(xAll) ) == 0
            flyNotMissing(f) = 1;

    
            % Run over trials (turns) and store each parameter p values in 
                % trial t in: Turncell_L(f,t,p):
            for t = 1:flyNumTrials(f)-1
                locsCurrTrial = (turnNumberAll == t);
                locsNextTrial = (turnNumberAll == t+1);
                % Store x,y coordinates:
                Turncell_L{f,t,1} = xAll(locsCurrTrial);
                Turncell_L{f,t,2} = yAll(locsCurrTrial);
                % Store trial's first frame #:
                framesInTrial = frameNumAll(locsCurrTrial);
                Turncell_L{f,t,3} = framesInTrial(1);
                % Store trial's turn direction:
                decisionsInTrial = turnDecisionAll(locsNextTrial);
                Turncell_L{f,t,4} = decisionsInTrial(1);
                % verify that turn decision is unique within a given trial:
                if decisionsInTrial(1) ~= 1 && decisionsInTrial(1) ~= 0
                    disp(['error - fly' num2str(f) ', trial ' ...
                        num2str(t) '- no decision val.']);
                end
            end


        end
    end
    
    
    % Omit sim flies with missing datapoints:
    Turncell_L = Turncell_L( flyNotMissing==1, :, : ); 

    % Save cell array of given datatype:
    save( ['datafiles/rawdata/Turndata_Short' dataType(5:end) '.mat'], ...
        'Turncell_L' );
    
    % Clear workspace variables:
    clearvars -except cdName dataTypes dataType d;

end

