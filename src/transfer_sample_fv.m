function transferred = transfer_sample_fv(jittered, target)
    target_fv = target.fv;
    transferred = zeros(size(target.image));
    transferred(:, :, 1) = target.lab;
    [c_y, c_x] = size(target.lab);
    for j = 1:c_x
        for i = 1:c_y
            idx = (i-1) + (j-1)*c_y + 1;
            best_match = compute_best_match(target_fv(idx, :), jittered.fv);
            transferred(i, j, 2:3) = jittered.ab(best_match, :);
        end
    end
    transferred(:,:,1) = transferred(:,:,1) * 100;
end