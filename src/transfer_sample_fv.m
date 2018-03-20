function transferred = transfer_sample_fv(samples, target)
    transferred = zeros(size(target.image));
    transferred(:, :, 1) = target.luminance;
    
    [c_y, c_x] = size(target.luminance);
    for j = 1:c_x
        for i = 1:c_y
            idx = (i-1) + (j-1)*c_y + 1;
            best_idx = compute_best_match(target.fv(:, idx), target.fv_w, samples.fv);
            transferred(i, j, 2:3) = samples.ab(:, best_idx);
        end
    end
    
    transferred(:,:,1) = transferred(:,:,1)*100;
end