function [featsW] = FeatureSelectionOptimization(source, samples, Kfs, mcCost, FEAT_SEL_METHOD)
%Optimizes the set of features to be used according to different criteria.

switch FEAT_SEL_METHOD
  case 0
    featsW = ones(1,size(source.fv_sp,1));
  case 1
  %SGA Feature Selection
    SGAFeatsSubset = speedyGA(source, samples, Kfs, mcCost);
    featsW = SGAFeatsSubset;
  case 2
  %MOEA Feature Space Optimization
    [Accs, MOEAFeatsWeigths] = crossValidationAccuracy([source.fv_sp' source.sp_clusters'], ...
      Kfs, mcCost);
    [~, midx] = max(Accs);
    MOEAFeatsWeigths = MOEAFeatsWeigths(:,midx)';
    featsW = MOEAFeatsWeigths;
  otherwise
    disp('Invalid FEAT_SEL_METHOD');
end

end

