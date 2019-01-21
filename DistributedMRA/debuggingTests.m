spmd 
    parentPartitionsLocalPart = getLocalPart(partitions(20, labindex));
    parentPartitionsLocalPart{:}(:,:)
%     knotsLocalPart = getLocalPart(knots(2, labindex));
%     knotsLocalPart{:}(:,:)
%     outputDataLocalPart = getLocalPart(outputData(100, 1, labindex));
%     outputDataLocalPart{:}(:,:)
end


spmd(nWorkersUsed)
thisKcholw
end

testMat = nan(4,4);
testMat(1:2,1) = [0;0]