function [best_match, ties] = compute_best_match(targetFV, featWeights, PoolFV)
%Find the vector in source that best represents the current target.
    [~, pool_len] = size(PoolFV);
    enlarged_target = repmat(targetFV, 1, pool_len);
    square_diff = (enlarged_target - PoolFV).^2;
  
    %Index of smaller distance
    weighted_sum = featWeights'*square_diff;
    [~, best_match] = min(weighted_sum);

    %% Debug 
    ties = sum(weighted_sum == weighted_sum(best_match));

end