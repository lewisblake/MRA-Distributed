function [ knotsSerial ] = process_knots_for_serial( knotsSubset, indexMatrixSubset, nLevelsInSerial, nRegions)
%% PROCESS_KNOTS_FOR_SERIAL
% This function take the subsets of knots
% codistributed arrays and processes then into cell arrays on the client
% for serial computation. The benefit of this is that it reduces the
% communication overhead caused by repeatedly calling the gather() function
% in the serial portion of the code. Then we can place the contents of
% these cell arrays into their associated portions of the codistributed
% arrays for parallel computation.
%%
% Calculate needed quantities
cummulativeRegions = cumsum(nRegions);
nTotalRegionsInSerial = cummulativeRegions(nLevelsInSerial);
% Pre-allocate space for serial cell arrays

knotsSerial = cell(nTotalRegionsInSerial,1);

for iLevel = 1:nLevelsInSerial
    for jRegion  = (cummulativeRegions(iLevel) - nRegions(iLevel) + 1): cummulativeRegions(iLevel)
        [firstRowContainingThisRegion, firstColContainingThisRegion] = find(indexMatrixSubset(:,:)==jRegion, 1 );
        knotsSerial(jRegion,1) = gather(knotsSubset(firstRowContainingThisRegion,firstColContainingThisRegion));
    end
end
end

