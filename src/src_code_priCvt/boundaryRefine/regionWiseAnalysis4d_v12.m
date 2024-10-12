function [newIdMap_refined,thresholdMap] = regionWiseAnalysis4d_v12(seedMap,eig3d,vid,option)
%% 20230616, Wei Zheng
%% refine boundary
newIdMap_boundary=boundaryRefineModule_v3(seedMap,eig3d, vid,option);
%% remove false positive
newIdMap_remove=removeFalsePositiveModule_v2(newIdMap_boundary,vid,option);
%% refine boundary again
idMap_current = region_sanity_check(seedMap.*(newIdMap_remove>0),1);
[newIdMap_refined,thresholdMap]=boundaryRefineModule_v3(idMap_current,eig3d, vid,option);
end