function best_match = compute_best_match(target_value, source)
   %Index of best match feature
   
   %Euclidean distance
   [y, ~] = size(source);
   enlarged_target = repmat(target_value, y, 1);
   square_diff = (enlarged_target - source).^2;
   
   weight_1 = 0.3;
   weight_2 = 1-weight_1;
   weighted_sum = weight_1 * square_diff(:, 1) + weight_2 * square_diff(:, 2);
   
   %Best match is the index of smaller distance
   [~, best_match] = min(weighted_sum);
end