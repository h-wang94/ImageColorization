function transferred = transfer_fv(csource, gsource, target)
    target_fv = target.fv;
    gsource_fv = gsource.fv;
    transferred = zeros(size(target.image));
    transferred(:, :, 1) = target.lab;
    [c_y, c_x, ~] = size(csource.image);
    num_fv = size(target_fv, 2)-1;
    for j = 1:c_x
        for i = 1:c_y
            idx = (i-1) + (j-1)*c_y + 1;
            best_match = compute_best_match(target_fv(idx, :), gsource_fv);
            new_j = ceil(best_match/c_y);
            new_i = best_match - (new_j-1)*c_y;
            transferred(i, j, 2:3) = csource.lab(new_i, new_j, 2:3);
        end
    end
    transferred(:,:,1) = transferred(:,:,1) * 100;
end