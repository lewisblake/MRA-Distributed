%% Playground for writing a new serial prior:

knotsSubset = gather(knots(1:nLevelsInSerial,:));
indexMatrixSerialSubset = indexMatrix(1:nLevelsInSerial,:);

[ knotsSerial ] = process_knots_for_serial( knotsSubset, indexMatrix(1:nLevelsInSerial,:), nLevelsInSerial, nRegions);
% Serial portion of the Prior
for iLevel = 1:nLevelsInSerial
   for jRegion  = (cumulativeRegions(iLevel) - nRegions(iLevel) + 1): cumulativeRegions(iLevel) 
    indexAncestry = find_ancestry(jRegion, nRegions, NUM_PARTITIONS_J);
 
    [thisRpriorChol, thisKcholBchol, ~, ~, ~] = create_prior(theta, ...
            NUM_LEVELS_M, knotsSerial([indexAncestry;jRegion],1),RpriorCholSerial(indexAncestry,1), ...
            KcBSerial(indexAncestry,1), [], varEps, []);
    
    firstColContainingThisRegion = find(indexMatrixSerialSubset(iLevel,:)==jRegion, 1 );   
    nTimesEachIndexIsRepeatedThisRow = nWorkersUsed/nRegions(iLevel);
    
    % Place objects from create_prior() in their serial containers
    RpriorCholSerial{jRegion} = thisRpriorChol;
    KcBSerial{jRegion} = thisKcholBchol;
    % Place objects from create_prior in their codistributed array
    % containers for execution later on
    RpriorChol(iLevel,firstColContainingThisRegion:firstColContainingThisRegion + nTimesEachIndexIsRepeatedThisRow - 1) = {thisRpriorChol};
    KcB(iLevel, firstColContainingThisRegion:firstColContainingThisRegion + nTimesEachIndexIsRepeatedThisRow - 1) = {thisKcholBchol};
    
   end
end

% Serial portion of the posterior

% Need to allocate space for a wtildePreviousSerial and
% AtildePreviousSerial since we use the children 1

% Level at which to begin parallel computations
nLevelToBeginInParallel = nLevelsInSerial + 1; % Should already have in MRA.m
% Find the last index of the children of the level at which we begin
% computing in serial. These children are regions that are within the
% parallel portion of computations
lastIndexOfSerialChildren = nRegions(nLevelToBeginInParallel+1)-1;
% Get only the row at which the last seruial children entry exists
[lastRowSerialChildren, ~] = find(indexMatrix(:,:) == lastIndexOfSerialChildren);
% Select the subset of the indexMatrix, wtildePrevious, and AtildePrevious which contains all serial levels and
% their children
indexMatrixSubset = indexMatrix(1:lastRowSerialChildren,:);
wtildeChildrenSubset = gather(wtildePrevious(1:lastRowSerialChildren,:));
AtildeChildrenSubset = gather(AtildePrevious(1:lastRowSerialChildren,:));


for iLevel = nLevelsInSerial:-1:1
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
        
        
        firstColContainingThisRegion = find(indexMatrixSubset(iLevel,:)==jRegion, 1 );   
        nTimesEachIndexIsRepeatedThisRow = nWorkersUsed/nRegions(iLevel);
    
        % This is where some trouble may lie....
        % Check that iLevel is the row we want.....
        if isPredicting
            RposteriorChol(iLevel, firstColContainingThisRegion:firstColContainingThisRegion + nTimesEachIndexIsRepeatedThisRow - 1) = {RposteriorCholj};
            Kcholw(iLevel, firstColContainingThisRegion:firstColContainingThisRegion + nTimesEachIndexIsRepeatedThisRow - 1) = {Kcholwj};
            KcholA(iLevel, firstColContainingThisRegion:firstColContainingThisRegion + nTimesEachIndexIsRepeatedThisRow - 1) = {KcholAj};
        else
            logLikelihood(iLevel,firstColContainingThisRegion:firstColContainingThisRegion + nTimesEachIndexIsRepeatedThisRow - 1) = logLikelihoodj;
        end 
        
    end
end




