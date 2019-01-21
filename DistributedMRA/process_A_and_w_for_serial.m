function [ AtildePreviousSerial, wtildePreviousSerial ] = process_A_and_w_for_serial(AtildeChildrenSubset, wtildeChildrenSubset, indexMatrixSubset, nLevelsInSerial, nRegions, NUM_PARTITIONS_J  )
%% PROCESS_A_AND_W_FOR_SERIAL
% This function takes the AtildeChildrenSubset and wtildeChildren subset
% corresponding to all regions to be calculated in serial and the children
% of nLevelsInSerial

% Calculate needed quantities
cumulativeRegions = cumsum(nRegions);
% Pre-allocate space for serial cell arrays
AtildePreviousSerial = cell(cumulativeRegions(nLevelsInSerial+1),1);
wtildePreviousSerial = cell(cumulativeRegions(nLevelsInSerial+1),1);

for iLevel = nLevelsInSerial:-1:1
    for jRegion = (cumulativeRegions(iLevel) - nRegions(iLevel) + 1): cumulativeRegions(iLevel)
        [indexChildren] = find_children(jRegion, nRegions, NUM_PARTITIONS_J);
        for kChild = 1:length(indexChildren)
            % Find the column and row containing the first occurance of this indexChild in the indexMatrix
            [childrenRow, childrenCol] = find(indexMatrixSubset(:,:) == indexChildren(kChild),1);
            % Find the correct row to enter
            rowToEnterForThisRegion = indexMatrixSubset(childrenRow, childrenCol);
            % Make cell arrays for serial computations
            wtildePreviousSerial(rowToEnterForThisRegion) = wtildeChildrenSubset(childrenRow, childrenCol);
            AtildePreviousSerial(rowToEnterForThisRegion) = AtildeChildrenSubset(childrenRow, childrenCol);
        end
    end
end

end

