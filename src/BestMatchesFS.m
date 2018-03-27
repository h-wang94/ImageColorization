function [best_matches, ties] = BestMatchesFS(targetFV, featWeights, PoolFV)
%Find the best matches in feature space between current target and pool of
%source samples.
    kTol = 1e-10;

    [~, pool_len] = size(PoolFV);
    enlarged_target = repmat(targetFV, 1, pool_len);
    square_diff = (enlarged_target - PoolFV).^2;
  
    weighted_sums = featWeights'*square_diff;
    [sorted_sums, best_matches] = sort(weighted_sums);
    
    ties = sum( (weighted_sums - sorted_sums(1)).^2 < kTol);
    best_matches = best_matches(1:max([ties 10]));
end