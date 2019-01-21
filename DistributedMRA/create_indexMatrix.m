function [indexMatrix] = create_indexMatrix(NUM_LEVELS_M, NUM_PARTITIONS_J, nRegions, nWorkersUsed, nLevelsInSerial,nTotalRegionsAssignedToEachWorker)
%% Create the indexMatrix
%  Input: NUM_LEVELS_M, NUM_PARTITIONS_J, nRegions, nWorkersUsed, nTotalRegionsAssignedToEachWorker
%  nLevelsInSerial, 
% 
%  Output: indexMatrix
%

nLevelToBeginParallel = nLevelsInSerial +1;
% Calculate necessary quantities
cummulativeRegions = cumsum(nRegions);
vectorOfRegionsAtNLevelsInSerial = nRegions(nLevelsInSerial):cummulativeRegions(nLevelsInSerial); % Create a vector of alll the regions where we will compute in serial
%matrixOfRegionsAtNLevelsInSerial = reshape(vectorOfRegionsAtNLevelsInSerial, [], nWorkersUsed); % Reshape that vector to extract the correct indexes in the following for-loop
vecOfRegionsAtFirstParallelLevel = nRegions(nLevelsInSerial+1):cummulativeRegions(nLevelsInSerial+1);
matrixOfRegionsAtFirstParallelLevel = reshape(vecOfRegionsAtFirstParallelLevel, [], nWorkersUsed);
% Calculate quantities needed for nTotalRegionsAssignedToEachWorker
maxLevelOnASingleRow = sum(nRegions <= nWorkersUsed); % How many times indices from a level are assigned to a worker
counter = 1:(NUM_LEVELS_M - maxLevelOnASingleRow);
nTotalRegionsAssignedToEachWorker = maxLevelOnASingleRow + sum(NUM_PARTITIONS_J.^counter);
% Pre-allocated memory for indexMatrix
indexMatrix = nan(nTotalRegionsAssignedToEachWorker, nWorkersUsed);


%spmd(nWorkersUsed) % Construct the indexMatrix in parallel. Works!
for iWorker = 1:nWorkersUsed  % Keep for testing 
    %tempVecOfRegionsAtNLevelsSerial = matrixOfRegionsAtNLevelsInSerial(:,iWorker); % Collect columns of matrix
    % Double check this later on
    % LB 12/1: Changed everything to be based on nLevelToBeginInParallel
    if nWorkersUsed == length(vecOfRegionsAtFirstParallelLevel)%length(vectorOfRegionsAtNLevelsInSerial) 
        indexMatrix(nLevelToBeginParallel, iWorker) = vecOfRegionsAtFirstParallelLevel(iWorker); % Regions above will not take vertical space in indexMatrix
        indexMatrix(1:nLevelsInSerial,iWorker) = find_ancestry(indexMatrix(nLevelToBeginParallel,  iWorker), nRegions, NUM_PARTITIONS_J);
        indexMatrix(nLevelToBeginParallel:end, iWorker) = find_branch_children(indexMatrix(nLevelToBeginParallel, iWorker), nRegions, NUM_PARTITIONS_J, NUM_LEVELS_M);
        
    %elseif nWorkersUsed > length(vecOfRegionsAtFirstParallelLevel) % LB: Even allow this case??    
    elseif nWorkersUsed < length(vecOfRegionsAtFirstParallelLevel)%length(vectorOfRegionsAtNLevelsInSerial)
        % LB: This finally seems to be working now, however I would not be
        % surprised if for some cases it fails. Right now it works for
        % nLevelsInSerial = 4, nWorkerUsed=4
         %Correct. Could be same as nRowsWithRepeatedEntries
        % Assign the regions allocated to this worker in a vector
        thisWorkerParallelAssignment = matrixOfRegionsAtFirstParallelLevel(:, iWorker);
        % Find all parents of the first parallel level. Store in vector
        tempVecOfParentsOfFirstParallelLevel=[];
        for i=1:length(thisWorkerParallelAssignment)           
            thisRegion  = thisWorkerParallelAssignment(i);
            [~, ~, thisRegionParents]= find_parent(thisRegion, nRegions, NUM_PARTITIONS_J);
            tempVecOfParentsOfFirstParallelLevel = [tempVecOfParentsOfFirstParallelLevel thisRegionParents];
        end
        tempVecOfParentsOfFirstParallelLevel = unique(tempVecOfParentsOfFirstParallelLevel);
        % Take the transpose of vector of all parents
        tempVecOfParentsOfFirstParallelLevel = tempVecOfParentsOfFirstParallelLevel';
        % Loop through the parent regions of the first level at which to compute in parallel.
        % Find all decendants of the parallel regions parents -  these are
        % unique to each column. Moreover, find the ancestors of these regions and store them both in vectors with repeats       
        tempVecOfParentsDescendents=[];
        tempVecOfParentsAncestors=[];
        for j = 1: length(tempVecOfParentsOfFirstParallelLevel)
            % Parent descendents
            theseParentsDescendents = find_branch_children(tempVecOfParentsOfFirstParallelLevel(j), nRegions, NUM_PARTITIONS_J, NUM_LEVELS_M);
            theseParentsDescendents = theseParentsDescendents(2:end); % Cut off first entry as it contains the index which was given as input
            tempVecOfParentsDescendents = [tempVecOfParentsDescendents theseParentsDescendents]; 
            % Parent ancestors
            theseParentsAncestors = find_ancestry(tempVecOfParentsOfFirstParallelLevel(j), nRegions, NUM_PARTITIONS_J);
            tempVecOfParentsAncestors = [tempVecOfParentsAncestors theseParentsAncestors];
        end
        
        % Select only the unique entries from both vectors
        tempVecOfParentsDescendents = unique(tempVecOfParentsDescendents);
        tempVecOfParentsAncestors = unique(tempVecOfParentsAncestors);
        
        % Stack them all in a vector and store it into the indexMatrix
        indexMatrix(:, iWorker) = [tempVecOfParentsAncestors; tempVecOfParentsOfFirstParallelLevel; tempVecOfParentsDescendents];
      
    end    
end

end

