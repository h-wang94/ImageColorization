function [best_matches, ties] = BestMatchesFS(targetFV, featWeights, PoolFV)
%Find the best matches in feature space between target and pool elements.
    [~, pool_len] = size(PoolFV);
    enlarged_target = repmat(targetFV, 1, pool_len);
    square_diff = (enlarged_target - PoolFV).^2;
  
    %Index of smaller distance
    weighted_sum = featWeights'*square_diff;
%     [~, best_match] = min(weighted_sum);
    sorted_sum = sort(weighted_sum);
    best_matches = sorted_sum(1);

    %% Debug 
    ties = sum(weighted_sum == weighted_sum(best_matches));

end