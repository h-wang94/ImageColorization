function best_match = compute_best_match(target_value, source)
   [y, x] = size(source);
   enlarged_target = repmat(target_value, y, 1);
   [~, best_match] = min(sum((enlarged_target - source).^2, 2));
end