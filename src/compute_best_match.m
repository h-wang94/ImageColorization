function best_match = compute_best_match(tgt_value, src_values_pool)
%Find the vector in source that best represents the current target.
    [feat_len, pool_size] = size(src_values_pool);
    enlarged_target = repmat(tgt_value, 1, pool_size);
    square_diff = (enlarged_target - src_values_pool).^2;

	weights = ones(1, feat_len)/feat_len;
    if (length(weights) ~= feat_len)
        display('Corrigir tamanho do vetor de pesos');
    end
    
    weighted_sum = weights*square_diff;
    %Best match is the index of smaller distance
    [~, best_match] = min(weighted_sum);
end