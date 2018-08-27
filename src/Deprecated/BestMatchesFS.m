function [bestMatches, bestDists, ties] = BestMatchesFS(targetFV, featWeights, PoolFVs, nCandidates)
%Find the best matches in feature space between current target and pool of
%source samples.
    kTol = 1e-10;

    [~, pool_len] = size(PoolFVs);
    enlarged_target = repmat(targetFV, 1, pool_len);
    square_diff = (enlarged_target - PoolFVs).^2;
  
    weighted_sums = featWeights'*square_diff;
    [sorted_sums, bestMatches] = sort(weighted_sums);

    %Return
    ties = sum( (weighted_sums - sorted_sums(1)).^2 < kTol);
    bestDists = sorted_sums(1:max([ties nCandidates-1]));
    bestMatches = bestMatches(1:max([ties nCandidates-1]));
end