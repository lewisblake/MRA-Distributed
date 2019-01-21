function [ knots, partitions, nRegions, outputData, predictionLocations, indexMatrix, nWorkersUsed ] = build_structure_in_parallel( NUM_LEVELS_M, ...
    NUM_PARTITIONS_J, NUM_KNOTS_r, domainBoundaries, offsetPercentage, nWorkersAssigned, nLevelsInSerial, varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% Check inputs and display progress check

% Check number of optional input arguments does not exceed two
numVarArgs = length(varargin); % LB: outputs 2
if numVarArgs > 2
    error('myfuns:build_structure:TooManyInputs', ...
        'requires at most 2 optional inputs');
end
% Display progress check
disp('Begining to build hierarchical grid structure in parallel...');

% Optional argument that can be passed is data
optArgs(1:numVarArgs) = varargin;
[ data, predictionVector ] = optArgs{:};
%% Set finest knot level
if numVarArgs == 2
    finestKnotLevel = NUM_LEVELS_M-1;
    indexEndFinestKnotLevel = NUM_PARTITIONS_J^(finestKnotLevel)-1;
else
    finestKnotLevel = NUM_LEVELS_M;
    indexEndFinestKnotLevel = NUM_PARTITIONS_J^(finestKnotLevel)-1;
end

%% Calculate quantities of interest
mLevels = 0:NUM_LEVELS_M-1; % Create a vector of levels
nRegions = NUM_PARTITIONS_J.^mLevels; % Vector of regions (partitions) at each level
totalRegions = sum(nRegions); % Calculate total number of regions
cummulativeRegions = cumsum(nRegions);


%% Calculate number of knots in each direction
if isinteger(sqrt(NUM_KNOTS_r))  % Assign knots.
    nKnotsX0 = sqrt(NUM_KNOTS_r); nKnotsX = sqrt(NUM_KNOTS_r); % Number of knots in x-direction
    nKnotsY0 = sqrt(NUM_KNOTS_r); nKnotsY = sqrt(NUM_KNOTS_r); % Number of knots in y-direction
else
    nKnotsX0 = ceil(sqrt(NUM_KNOTS_r)); nKnotsX = ceil(sqrt(NUM_KNOTS_r));
    nKnotsY0 = NUM_KNOTS_r/nKnotsX0; nKnotsY = NUM_KNOTS_r/nKnotsX;
end

%% NEW STRATEGY

% From nWorkersAssigned, calculate nWorkersUsed
nWorkersUsed = nWorkersAssigned; % Hardcoded: need to generalize
nLevelToBeginInParallel = nLevelsInSerial + 1; % Level to begin computing in parallel is the one directly after we finish computing in serial

% Make assumption that tiles are equally distributed across workers
nTilesAssignedToEachWorker = nRegions(nLevelToBeginInParallel)/nWorkersUsed; % Set number of tiles each worker is responsible for
firstIndexOfLevelToBeginInParallel = NUM_PARTITIONS_J^(nLevelsInSerial); % Find a cleaner way to find correct index!!
% Use the first index from level at which we begin computing in parallel to
% find the length of the branch "below" that region. The same number will
% work for all regions at the level at which we begin computing in parallel
indexsOfAllChildrenOnThisBranch = find_branch_children(firstIndexOfLevelToBeginInParallel,nRegions, NUM_PARTITIONS_J, NUM_LEVELS_M);
nRegionsInBranch = length(indexsOfAllChildrenOnThisBranch);
nRegionsAboveBranch = length(find_ancestry(firstIndexOfLevelToBeginInParallel, nRegions, NUM_PARTITIONS_J));
nTotalRegionsInBranch = nRegionsAboveBranch + nRegionsInBranch;
%nRegionsAtFinestLevel = nRegions(NUM_LEVELS_M);
nRegionsAtFinestLevelForEachWorker = (nRegions(NUM_LEVELS_M)/nWorkersUsed);
% Calculate quantities needed for nTotalRegionsAssignedToEachWorker
maxLevelOnASingleRow = sum(nRegions <= nWorkersUsed); % How many times indices from a level are assigned to a worker
counter = 1:(NUM_LEVELS_M - maxLevelOnASingleRow); % Count number of times indicies from a level are going to be repeated in each indexMatrix column
nTotalRegionsAssignedToEachWorker = maxLevelOnASingleRow + sum(NUM_PARTITIONS_J.^counter);
%nTotalRegionsAssignedToEachWorker = nRegionsAboveBranch + nTilesAssignedToEachWorker*nRegionsInBranch; % Should Work? Need verification



%% Create matrix to store continuous index for all regions
[indexMatrix] = create_indexMatrix( NUM_LEVELS_M, NUM_PARTITIONS_J, nRegions, nWorkersUsed, nLevelsInSerial, nTotalRegionsAssignedToEachWorker);
% Find the index within the indexMatrix corresponding to the finest level at which the knots are not set to the data.
indexOfFinestKnotLevelWithinIndexMatrix = find(indexMatrix(:,end)==indexEndFinestKnotLevel);
nRowsWithRepeatedEntriesInIndexMatrix = sum(nRegions < nWorkersUsed);
%% Pre-allocate memory for codistributed arrays
spmd(nWorkersUsed)
    codistributionScheme = codistributor1d(2); % Distribute across the second dimension
    % Create codistributed cell arrays
    knots = cell(nTotalRegionsAssignedToEachWorker,  nWorkersUsed, codistributionScheme);
    outputData = cell(nRegionsAtFinestLevelForEachWorker, nWorkersUsed, codistributionScheme);
    partitions = cell(indexOfFinestKnotLevelWithinIndexMatrix+1,  nWorkersUsed, codistributionScheme); % 12/3: LB changed from nTotalRegionsAssignedToEachWorker. 12/5: last entry in each column is empty
    predictionLocations = cell(nRegionsAtFinestLevelForEachWorker, nWorkersUsed, codistributionScheme);
end
%% Construct zeroth level
xMin0 = domainBoundaries(1); xMax0 = domainBoundaries(2);
yMin0 = domainBoundaries(3); yMax0 = domainBoundaries(4);
% Edge buffer added to xMax0 and yMax0 to
% include all observations on the boundary at the zeroth level
xMax0 = xMax0 + (offsetPercentage/2)*(xMax0 - xMin0);
yMax0 = yMax0 + (offsetPercentage/2)*(yMax0 - yMin0);
% Create the knots at the coarsest resolution
[ knotsX, knotsY ] = create_knots(xMin0, xMax0, nKnotsX0, yMin0, yMax0, nKnotsY0, offsetPercentage);
% Each branch at the coarsest resolution will have the same knots
knots(1,1:nWorkersUsed) = {[knotsX(:), knotsY(:)]}; % Knots at coarsest resolution are part of the parental hierarchy for all branches
% Create the paritition at the coarsest resolution
[ xMin, xMax, yMin, yMax ] = create_partition(xMin0, xMax0, yMin0, yMax0, NUM_PARTITIONS_J);
% Each branch at the coarest resolution will be contained within same
% region
partitions(1,1:nWorkersUsed) = {[ xMin, xMax, yMin, yMax ]};
%% Create portion of partitions that has repeated entries across columns
%vectorOfIndiciesForNonUniqueIndexMatrixRows = 1:NUM_PARTITIONS_J:nWorkersUsed;
% LB: Be sure to check this condition for other test cases. 12/1:  Doesn't work.
% Lb: Changing from sum(nRegions < nTilesAssignedToEachWorker) to sum(nRegions < nWorkersUsed)
% 12/1 NEW ATTEMPT
for iRow = 2:nRowsWithRepeatedEntriesInIndexMatrix
    nTimesEachIndexIsRepeatedThisRow = nWorkersUsed/nRegions(iRow);
    jBeginningColumnEntries = 1:nTimesEachIndexIsRepeatedThisRow:nWorkersUsed;
    parentsOfThisRow = unique(indexMatrix(iRow-1,:));
    columnToEnterMatrix = reshape(jBeginningColumnEntries, [], length(parentsOfThisRow));
    for jj = 1:length(parentsOfThisRow)
        jParent = parentsOfThisRow(jj);
        % WLOG, Get first column entry (not specific) of the iRow above
        % containing this jParent
        firstColContainingThisParent = find(indexMatrix(iRow-1,:)==jParent, 1 );
        % Get partitions for this jParent
        thisParentsPartitions = gather(partitions(iRow-1,firstColContainingThisParent));
        xMin = thisParentsPartitions{:}(:,1);
        xMax = thisParentsPartitions{:}(:,2);
        yMin = thisParentsPartitions{:}(:,3);
        yMax = thisParentsPartitions{:}(:,4);       
        for lPartition = 1:NUM_PARTITIONS_J
            % Get correct column to enter
            columnToEnter = columnToEnterMatrix(lPartition, jj);           
            % Create knots & partitions for...WHAT????
            [knotsX,knotsY] = create_knots(xMin(lPartition), xMax(lPartition), nKnotsX, yMin(lPartition), yMax(lPartition), nKnotsY, offsetPercentage);
            [ xMinTemp, xMaxTemp, yMinTemp, yMaxTemp ] = create_partition(xMin(lPartition), xMax(lPartition), yMin(lPartition), yMax(lPartition), NUM_PARTITIONS_J);
            % Place knots partitions in array entries corresponding to
            % this index
            knots(iRow, columnToEnter:columnToEnter+nTimesEachIndexIsRepeatedThisRow-1) = {[knotsX(:),knotsY(:)]};
            partitions(iRow, columnToEnter:columnToEnter+nTimesEachIndexIsRepeatedThisRow-1) = {[ xMinTemp, xMaxTemp, yMinTemp, yMaxTemp ]};            
        end
    end
end

%% Create partitions and knots up until finestKnotLevel
spmd(nWorkersUsed)
    for iIndex = nRowsWithRepeatedEntriesInIndexMatrix+1:indexOfFinestKnotLevelWithinIndexMatrix
        indexCurrent = indexMatrix(iIndex, labindex); % Assign (continuous) indexCurrent from indexMatrix.
        [~, ~, indexParent] = find_parent(indexCurrent, nRegions, NUM_PARTITIONS_J); % Find (continuous) parent of indexCurrent
        indexOfParentWithinIndexMatrix = find(indexMatrix(:, labindex)==indexParent); % Find the index of indexParent within the indexMatrix
        % Get entire vectors of partitions from the parent partition
        parentPartitionsLocalPart = getLocalPart(partitions(indexOfParentWithinIndexMatrix, labindex));
        xMin = parentPartitionsLocalPart{:}(:,1);
        xMax = parentPartitionsLocalPart{:}(:,2);
        yMin = parentPartitionsLocalPart{:}(:,3);
        yMax = parentPartitionsLocalPart{:}(:,4);
        % Assign correct xMin,...,yMax for THIS indexCurrent
        thesePartitionBoundaries = find(indexCurrent == find_children(indexParent, nRegions, NUM_PARTITIONS_J));
        thisXMin = xMin(thesePartitionBoundaries); thisXMax = xMax(thesePartitionBoundaries);
        thisYMin = yMin(thesePartitionBoundaries); thisYMax = yMax(thesePartitionBoundaries);
        % With the correct thisXMin,...,thisYMax, create the knots and
        % partitions for this indexCurrent and store them in their
        % corrresponding locations within indexMatrix
        [knotsX,knotsY] = create_knots(thisXMin, thisXMax, nKnotsX, thisYMin, thisYMax, nKnotsY, offsetPercentage);
        [ xMinTemp, xMaxTemp, yMinTemp, yMaxTemp ] = create_partition(thisXMin, thisXMax, thisYMin, thisYMax, NUM_PARTITIONS_J);
        indexOfIndexCurrentWithinIndexMatrix = find(indexMatrix(:, labindex)==indexCurrent);
        knots(indexOfIndexCurrentWithinIndexMatrix, labindex) = {[knotsX(:),knotsY(:)]};
        partitions(indexOfIndexCurrentWithinIndexMatrix, labindex) = {[ xMinTemp, xMaxTemp, yMinTemp, yMaxTemp ]};
    end
end

%% Special construct to find knots for finest resolution region
if numVarArgs == 2 % If data is sent to build_structure_in_parallel
    spmd(nWorkersUsed)
        for iIndex = indexOfFinestKnotLevelWithinIndexMatrix + 1 : length(indexMatrix)
            mCounterIndex = iIndex - indexOfFinestKnotLevelWithinIndexMatrix; % Create a counter index to get into outputData and predictionLocations starting at index 1 to maintain minimal size
            indexCurrent = indexMatrix(iIndex, labindex); % Find the indexCurrent from indexMatrix
            [~, ~, indexParent] = find_parent(indexCurrent, nRegions, NUM_PARTITIONS_J); % Find unique parent of indexCurrent
            indexOfParentWithinIndexMatrix = find(indexMatrix(:, labindex)==indexParent); % Find the index within the indexMatrix of indexParent
            % On this worker, get the local part of the partitions
            % corresponding to indexParent as determined by
            % indexParent's index within the indexMatrix
            parentPartitionsLocalPart = getLocalPart(partitions(indexOfParentWithinIndexMatrix, labindex));
            xMin = parentPartitionsLocalPart{:}(:,1);
            xMax = parentPartitionsLocalPart{:}(:,2);
            yMin = parentPartitionsLocalPart{:}(:,3);
            yMax = parentPartitionsLocalPart{:}(:,4);
            % Assign correct xMin,...,yMax for THIS indexCurrent
            thesePartitionBoundaries = find(indexCurrent == find_children(indexParent, nRegions, NUM_PARTITIONS_J));
            thisXMin = xMin(thesePartitionBoundaries); thisXMax = xMax(thesePartitionBoundaries);
            thisYMin = yMin(thesePartitionBoundaries); thisYMax = yMax(thesePartitionBoundaries);
            % Collect the data needed for this region
            ind = find(data(:,1) >= thisXMin & data(:,1) < thisXMax & data(:,2) >= thisYMin & data(:,2) < thisYMax); % Find the indicies for the data within the data matrix that are within thisXMin,...,thisYMax
            knotsX=data(ind,1); knotsY=data(ind,2); % Set the knots to these data
            indexOfIndexCurrentWithinIndexMatrix = find(indexMatrix(:, labindex)==indexCurrent); % Find the index of indexCurrent within the indexMatrix
            knots(indexOfIndexCurrentWithinIndexMatrix, labindex) = {[knotsX(:),knotsY(:)]}; % Place knotsX(:) and knotsY(:) in their corresponding location within the knots codistributed cell
            outputData(mCounterIndex,labindex) = {data(ind,3)}; % Place data within thisXMin,...,thisYMax within outputData
            data(ind,:) = []; % Eliminate the data that has already been assigned to a region, speeds up subsequent searching
            % Vinay's addition
            if ~isnan(predictionVector) % If predicting
                predictionIndex = find(predictionVector(:,1) >= thisXMin & predictionVector(:,1) < thisXMax & predictionVector(:,2) >= thisYMin & predictionVector(:,2) < thisYMax); % Find the predictionVector locations within this region
                predictionLocations(mCounterIndex, labindex) = {predictionVector(predictionIndex,:)}; % Assign predictionVector locations within this region to corresponding entry of predictionLocations codistributed cell
            else
                predictionLocations = NaN;
            end          
        end
    end
    disp('Building the hierarchical structure complete.')
end


