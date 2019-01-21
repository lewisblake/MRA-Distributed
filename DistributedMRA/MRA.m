function [ sumLogLikelihood, predictions ] = MRA( theta, outputData, knots, ...
    NUM_LEVELS_M, NUM_PARTITIONS_J, nRegions, indexMatrix, isPredicting, nLevelsInSerial, nWorkersUsed, verbose, varargin)
%% MRA.m is the main Multi-resolution approximation function
%
% Input: theta, data, knots, MAX_LEVEL_M, NUM_PARTITIONS_J, nRegions, varargin
%
% Output: sumLogLikelihood, predictions
%%
% Check number of optional input arguments
numVarArgsIn = length(varargin);
if numVarArgsIn > 2
    error('myfuns:create_prior:TooManyInputs', ...
        'requires at most 2 optional inputs');
end
% set defaults for optional inputs
optionalArguments = {0 NaN};
% overwrite the ones specified in varargin.
optionalArguments(1 : numVarArgsIn) = varargin;
[ varEps, predictionLocations ] = optionalArguments{:};
% Calculate key quantities
%totalRegions = sum(nRegions);
cumulativeRegions = cumsum(nRegions);
nLevelToBeginInParallel = nLevelsInSerial + 1;
nTilesAssignedToEachWorker = nRegions(nLevelToBeginInParallel)/nWorkersUsed; % Set number of tiles each worker is responsible for
firstIndexOfLevelToBeginInParallel = NUM_PARTITIONS_J^(nLevelToBeginInParallel-1); % Find a cleaner way to find correct index!!
indexsOfAllChildrenOnThisBranch = find_branch_children(firstIndexOfLevelToBeginInParallel,nRegions, NUM_PARTITIONS_J, NUM_LEVELS_M);
nRegionsInBranch = length(indexsOfAllChildrenOnThisBranch);
nRegionsAboveBranch = length(find_ancestry(firstIndexOfLevelToBeginInParallel, nRegions, NUM_PARTITIONS_J));
%nTotalRegionsInBranch = nRegionsAboveBranch + nRegionsInBranch;
nRegionsAtFinestLevelForEachWorker = (nRegions(NUM_LEVELS_M)/nWorkersUsed);
% 1/7 LB: Calculate quantities needed for nTotalRegionsAssignedToEachWorker
maxLevelOnASingleRow = sum(nRegions <= nWorkersUsed); % How many times indices from a level are assigned to a worker
counter = 1:(NUM_LEVELS_M - maxLevelOnASingleRow);
nTotalRegionsAssignedToEachWorker = maxLevelOnASingleRow + sum(NUM_PARTITIONS_J.^counter);
%nTotalRegionsAssignedToEachWorker = nRegionsAboveBranch + nTilesAssignedToEachWorker*nRegionsInBranch; % Should Work? Need verification
%nRowsWithRepeatedEntriesInIndexMatrix = sum(nRegions < nTilesAssignedToEachWorker);
nTotalRegionsInSerial = cumulativeRegions(nLevelsInSerial);


%% Pre-allocate space for codistributed arrays
if verbose
    disp('In MRA.m: Pre-allocating space for objects ...')
end
spmd(nWorkersUsed)
    codistributionScheme = codistributor1d(2); % Distribute across the second dimension
    % Create codistributed arrays. LB: CHECK SIZES OF THESE AGAIN ONCE COMPLETE!
    RpriorChol = cell(nTotalRegionsAssignedToEachWorker,  nWorkersUsed, codistributionScheme);
    KcB = cell(nTotalRegionsAssignedToEachWorker,  nWorkersUsed, codistributionScheme);
    AtildePrevious = cell(nTotalRegionsAssignedToEachWorker,  nWorkersUsed, codistributionScheme);
    wtildePrevious = cell(nTotalRegionsAssignedToEachWorker,  nWorkersUsed, codistributionScheme);
    logLikelihood = nan(nTotalRegionsAssignedToEachWorker, nWorkersUsed, codistributionScheme);
    if isPredicting % If predicting, pre-allocate space for necessary quantities
        posteriorPredictionMean =  cell(nRegionsAtFinestLevelForEachWorker, nWorkersUsed, codistributionScheme); % Only need to store values at finest resolution
        posteriorPredictionVariance = cell(nRegionsAtFinestLevelForEachWorker, nWorkersUsed, codistributionScheme); % Only need to store values at finest resolution
        Btilde = cell(nRegionsAtFinestLevelForEachWorker, nWorkersUsed, codistributionScheme); % Only need to store values at finest resolution
        predictions = cell(nRegionsAtFinestLevelForEachWorker, nWorkersUsed, codistributionScheme); % Only need to store values at finest resolution
        RposteriorChol = cell(nTotalRegionsAssignedToEachWorker-nRegionsAtFinestLevelForEachWorker,  nWorkersUsed, codistributionScheme); % LB: Doesn't need to include entries for finest resolution
        KcholA = cell(nTotalRegionsAssignedToEachWorker-nRegionsAtFinestLevelForEachWorker,  nWorkersUsed, codistributionScheme); % LB: Doesn't need to include finest resolution
        Kcholw = cell(nTotalRegionsAssignedToEachWorker-nRegionsAtFinestLevelForEachWorker,  nWorkersUsed, codistributionScheme); % LB: Doesn't need to include entries for finest resolution
    else % Workaround to pass correct object to create_prior
        % LB: Check this. Not sure the dimensions needed here.
        predictionLocations = num2cell(nan(nTotalRegionsAssignedToEachWorker, nWorkersUsed, codistributionScheme));
    end  
end
%% Pre-allocate space for cell arrays used in serial computations
% For creating the prior, except for knots, which are made by
% process_knots_for_serial()
RpriorCholSerial = cell(nTotalRegionsInSerial,1);
KcBSerial = cell(nTotalRegionsInSerial, 1);
% For calculating the posterior...

%% Create the prior distribution
lastIndexOfSerialLevel = nRegions(nLevelsInSerial+1)-1;
[lastRowInSerial, lastColInSerial] = find(indexMatrix(:,:) == lastIndexOfSerialLevel,1);% Assuming here that the last region to compute in serial is located in the last indexMatrix column
knotsSubset = gather(knots(1:lastRowInSerial,:));
indexMatrixSubset = indexMatrix(1:lastRowInSerial,:);
% Process knots for serial computation
[ knotsSerial ] = process_knots_for_serial( knotsSubset, indexMatrixSubset, nLevelsInSerial, nRegions);
disp('Creating the prior in serial ...')
for iLevel = 1:nLevelsInSerial
   for jRegion  = (cumulativeRegions(iLevel) - nRegions(iLevel) + 1): cumulativeRegions(iLevel) 
    % Find the ancestry for this jRegion
    indexAncestry = find_ancestry(jRegion, nRegions, NUM_PARTITIONS_J);
    % Create the prior with create_prior() function
    [thisRpriorChol, thisKcholBchol, ~, ~, ~] = create_prior(theta, ...
            NUM_LEVELS_M, knotsSerial([indexAncestry;jRegion],1),RpriorCholSerial(indexAncestry,1), ...
            KcBSerial(indexAncestry,1), [], varEps, []);
    % Fina quantities that determine where to place thisRpriorChol,
    % thisKcholBchol into their codistributed arrays
    [firstRowContainingThisRegion,firstColContainingThisRegion] = find(indexMatrixSubset(:,:)==jRegion, 1 );   
    nTimesEachIndexIsRepeatedThisRow = ceil(nWorkersUsed/nRegions(iLevel)); % 1/4 LB: Added ceiling function around to ensure positive integer double    
    % Place objects from create_prior() in their serial containers
    RpriorCholSerial{jRegion} = thisRpriorChol;
    KcBSerial{jRegion} = thisKcholBchol;
    % Place objects from create_prior in their codistributed array
    % containers for execution later on
    RpriorChol(firstRowContainingThisRegion,firstColContainingThisRegion:firstColContainingThisRegion + nTimesEachIndexIsRepeatedThisRow - 1) = {thisRpriorChol};
    KcB(firstRowContainingThisRegion, firstColContainingThisRegion:firstColContainingThisRegion + nTimesEachIndexIsRepeatedThisRow - 1) = {thisKcholBchol};    
   end
end
knotsSubset = []; knotsSerial = []; % To save memory?
disp('Serial section of the prior complete.');
disp('Creating the prior in parallel ...');
% Finish calculating prior in parallel
spmd(nWorkersUsed) %nWorkersUsed
    % In parallel, loop over indexMatrix rows
    mCounterIndex = 1;
    for iIndexMatrixRow = lastRowInSerial+1:nTotalRegionsAssignedToEachWorker
        index = indexMatrix(iIndexMatrixRow, labindex);
        indexAncestry = find_ancestry( index, nRegions, NUM_PARTITIONS_J );
        % Fill vector with location of indexAncestry in indexMatrix - an area
        %for optimization
        indexAncestryInIndexMatrix = nan(length(indexAncestry),1);
        for k = 1:length(indexAncestry)
           indexAncestryInIndexMatrix(k) = find(indexMatrix(1:iIndexMatrixRow,labindex) == indexAncestry(k)); 
        end 
        
        % Get local part of the objects needed to create the prior
        knotsLocalPart = getLocalPart(knots([indexAncestryInIndexMatrix;iIndexMatrixRow], labindex));        
        RpriorCholLocalPart = getLocalPart(RpriorChol(indexAncestryInIndexMatrix, labindex));       
        KcBLocalPart = getLocalPart(KcB(indexAncestryInIndexMatrix, labindex));
        
        if index < nRegions(NUM_LEVELS_M) % If not dealing with an index on the finest resolution...
            % don't send create_prior any data or predictionLocations
            [thisRpriorChol, thisKcholBchol, ~, ~, ~] = create_prior(theta, ...
                NUM_LEVELS_M, knotsLocalPart, RpriorCholLocalPart,...
                KcBLocalPart, [], varEps, []);
        else
            dataLocalPart = getLocalPart(outputData(mCounterIndex, labindex));
            predictionLocationsLocalPart = getLocalPart(predictionLocations(mCounterIndex, labindex));
            % Create the prior
            [thisRpriorChol, thisKcholBchol, thisAtj, thiswtj, thisRetLikPred] = create_prior(theta, ...
                NUM_LEVELS_M, knotsLocalPart, RpriorCholLocalPart,...
                KcBLocalPart, dataLocalPart{:}, varEps, predictionLocationsLocalPart{:});

            AtildePrevious(iIndexMatrixRow, labindex) = {thisAtj};
            wtildePrevious(iIndexMatrixRow, labindex) = {thiswtj};
            if isPredicting  % If predicting
                posteriorPredictionMean(mCounterIndex, labindex) = {thisRetLikPred{1}};
                posteriorPredictionVariance(mCounterIndex, labindex) = {thisRetLikPred{2}};
                Btilde(mCounterIndex, labindex) = {thisRetLikPred{3}};
            else % If not predicting
                logLikelihood(iIndexMatrixRow,labindex) = thisRetLikPred;
            end     
            mCounterIndex = mCounterIndex + 1;
        end
        RpriorChol(iIndexMatrixRow, labindex) = {thisRpriorChol};
        KcB(iIndexMatrixRow, labindex) = {thisKcholBchol};
    end 
end
disp('Parallel section of the prior complete.');

%% Create the posterior distribution
disp('Calculating parallel section of the posterior ...');
lastIndexOfSecondFinestLevel = nRegions(NUM_LEVELS_M)-1;
lastRowBeforeFinestLevel = find(indexMatrix(:,end)==lastIndexOfSecondFinestLevel);
spmd(nWorkersUsed)
%     % Pre-allocate posterior quantities. LB: don't need??
%     AtildeCurrent = cell(nTotalRegionsAssignedToEachWorker,  nWorkersUsed, codistributionScheme);
%     wtildeCurrent = cell(nTotalRegionsAssignedToEachWorker,  nWorkersUsed, codistributionScheme);
       
    for iIndexMatrixRow = lastRowBeforeFinestLevel:-1:lastRowInSerial+1
        index = indexMatrix(iIndexMatrixRow, labindex);
        [indexChildren] = find_children(index, nRegions, NUM_PARTITIONS_J);
        
        indexChildrenInIndexMatrix = nan(length(indexChildren),1);
        for k = 1:length(indexChildren)
        indexChildrenInIndexMatrix(k) = find(indexMatrix(:, labindex) == indexChildren(k));
        end
        
        RpriorCholj = getLocalPart(RpriorChol(iIndexMatrixRow, labindex));
        wtildeChildren = getLocalPart(wtildePrevious(indexChildrenInIndexMatrix, labindex));
        AtildeChildren = getLocalPart(AtildePrevious(indexChildrenInIndexMatrix, labindex));
        
        % Calculate posterior_inference()
        [ wtildeCurrentj, AtildeCurrentj, logLikelihoodj, ...
            RposteriorCholj, Kcholwj, KcholAj ] = posterior_inference(RpriorCholj{:}, ...
            wtildeChildren, AtildeChildren);
        
        % LB: changed from wtileCurrent to wtildePrevious. Seems to be
        % working however need to check back later
        wtildePrevious(iIndexMatrixRow,labindex) = {wtildeCurrentj};
        AtildePrevious(iIndexMatrixRow, labindex) = {AtildeCurrentj};
        
        if isPredicting
            RposteriorChol(iIndexMatrixRow, labindex) = {RposteriorCholj};
            Kcholw(iIndexMatrixRow, labindex) = {Kcholwj};
            KcholA(iIndexMatrixRow, labindex) = {KcholAj};
        else
            logLikelihood(iIndexMatrixRow,labindex) = logLikelihoodj;
        end
    end
end
disp('Parallel section of the posterior complete.');
disp('Calculating the serial section of the posterior ...');
% Find the last index of the children of the level at which we begin
% computing in serial. These children are regions that are within the
% parallel portion of computations
lastIndexOfSerialChildren = nRegions(nLevelToBeginInParallel+1)-1;
% Get only the row at which the last seruial children entry exists
[lastRowSerialChildren, ~] = find(indexMatrix(:,:) == lastIndexOfSerialChildren);
% Select the subset of the indexMatrix, wtildePrevious, and AtildePrevious which contains all serial levels and
% their children
indexMatrixSubset = indexMatrix(1:lastRowSerialChildren,:); % Redefine from above to include children ot nLevelsInSerial
wtildeChildrenSubset = gather(wtildePrevious(1:lastRowSerialChildren,:));
AtildeChildrenSubset = gather(AtildePrevious(1:lastRowSerialChildren,:));
% Process AtildePrevious, wtildePrevious for serial computation
[ AtildePreviousSerial, wtildePreviousSerial ] = process_A_and_w_for_serial(AtildeChildrenSubset, ...
    wtildeChildrenSubset, indexMatrixSubset, nLevelsInSerial, nRegions, NUM_PARTITIONS_J); 
% Loop through remaining regions in serial
for iLevel = nLevelsInSerial:-1:1 % 1/7 LB: changed going from nLevelsInSerial to lastRowInSerial. Actually should be nLevelsInSerial since loop over jRegion
    for jRegion = (cumulativeRegions(iLevel) - nRegions(iLevel) + 1) : cumulativeRegions(iLevel)
        [ indexChildren ] = find_children( jRegion, nRegions, NUM_PARTITIONS_J );
        % Calculate posterior quantities
        RpriorCholj = RpriorCholSerial{jRegion};
        wtildeChildren = wtildePreviousSerial(indexChildren);
        AtildeChildren = AtildePreviousSerial(indexChildren);
        
        % Calculate posterior_inference()
        [ wtildeCurrentj, AtildeCurrentj, logLikelihoodj, ...
            RposteriorCholj, Kcholwj, KcholAj ] = posterior_inference( RpriorCholj, ...
            wtildeChildren, AtildeChildren );
        
        wtildePreviousSerial{jRegion} = wtildeCurrentj;
        AtildePreviousSerial{jRegion} = AtildeCurrentj;
        
        % 1/7 LB : commeneted out row two below this one and added below
        [firstRowContainingThisRegion,firstColContainingThisRegion] = find(indexMatrixSubset(:,:)==jRegion, 1 );  
        %firstColContainingThisRegion = find(indexMatrixSubset(iLevel,:)==jRegion, 1 );   
        nTimesEachIndexIsRepeatedThisRow = ceil(nWorkersUsed/nRegions(iLevel)); % 1/7 LB: added ceil() function around calculation
    
        % This is where some trouble may lie....
        % Check that iLevel is the row we want..... % 1/7 Lb changed from
        % iLevel to firstRowContainingThisRegion
        if isPredicting
            RposteriorChol(firstRowContainingThisRegion, firstColContainingThisRegion:firstColContainingThisRegion + nTimesEachIndexIsRepeatedThisRow - 1) = {RposteriorCholj};
            Kcholw(firstRowContainingThisRegion, firstColContainingThisRegion:firstColContainingThisRegion + nTimesEachIndexIsRepeatedThisRow - 1) = {Kcholwj};
            KcholA(firstRowContainingThisRegion, firstColContainingThisRegion:firstColContainingThisRegion + nTimesEachIndexIsRepeatedThisRow - 1) = {KcholAj};
        else
            logLikelihood(firstRowContainingThisRegion,firstColContainingThisRegion:firstColContainingThisRegion + nTimesEachIndexIsRepeatedThisRow - 1) = logLikelihoodj;
        end 
        
    end
end
disp('Serial section of the posterior complete.');
% Calculate the sum loglikelihood
sumLogLikelihood = sum(sum(gather(logLikelihood)));


%% Spatial prediction
if isPredicting
    disp('Beginning the spatial prediction ...')
    spmd(nWorkersUsed)
        for iIndexMatrixRow = lastRowBeforeFinestLevel+1:nTotalRegionsAssignedToEachWorker
            mCounterIndex = iIndexMatrixRow - lastRowBeforeFinestLevel;
            if NUM_LEVELS_M > 0
                % Set up appropriate indicies
                index = indexMatrix(iIndexMatrixRow, labindex);
                indexAncestry = find_ancestry(index, nRegions,NUM_PARTITIONS_J);
                % Fill vector with location of indexAncestry in indexMatrix - an area
                %for optimization
                indexAncestryInIndexMatrix = nan(length(indexAncestry),1);
                for k = 1:length(indexAncestry)
                    indexAncestryInIndexMatrix(k) = find(indexMatrix(1:iIndexMatrixRow,labindex) == indexAncestry(k));
                end
                % Collect the appropriate inputs to make predictions at
                % this region
                thisPosteriorPredictionMean = getLocalPart(posteriorPredictionMean(mCounterIndex, labindex));
                thisPosteriorPredictionVariance = getLocalPart(posteriorPredictionVariance(mCounterIndex, labindex));
                thisBtilde = getLocalPart(Btilde(mCounterIndex, labindex));
                thisRposteriorChol = getLocalPart(RposteriorChol(indexAncestryInIndexMatrix, labindex));
                thisKcholA = getLocalPart(KcholA(indexAncestryInIndexMatrix, labindex));
                thisKcholw = getLocalPart(Kcholw(indexAncestryInIndexMatrix, labindex));
                % Use spatial_prediction() to make predictions. Place
                % output into a cell by wrapping function in {}.
                predictions(mCounterIndex,labindex) = {spatial_prediction(thisPosteriorPredictionMean{:}, ...
                    thisPosteriorPredictionVariance{:}, thisBtilde{:}, thisRposteriorChol, ...
                    thisKcholA, thisKcholw)};
            else
                predictions(mCounterIndex, labindex) = {[posteriorPredictionMean(mCounterIndex,labindex), posteriorVariance(mCounterIndex,labindex)]};
            end
        end        
    end
    disp('Spatial prediction complete.');
end
disp('MRA.m execution complete.');
end