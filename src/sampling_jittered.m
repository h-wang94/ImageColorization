function jittered = sampling_jittered(image_luminance, image_lab, num_samples)
    % assumes num_samples is square number. num_samples = n^2 where n is an
    % integer
    fv = {};
    jittered = {};
    n = sqrt(num_samples);
    [y, x] = size(image_luminance);
    y_step = floor(y / n);
    x_step = floor(x / n);
    sample_idx = randi([1, y_step*x_step], 1, num_samples*2);
    fv.sds = zeros([1,num_samples]);
    fv.luminance = zeros([1, num_samples]);
    fv.ab = zeros([num_samples, 2]);
    jittered.points = zeros([num_samples, 2]);
    idx = 1;
    amt_to_pad = 2;
    padded = padarray(image_luminance, [amt_to_pad, amt_to_pad]);
    for j = 1:x_step:x-x_step-1
        for i = 1:y_step:y-y_step-1
            value = sample_idx(idx);
            % determines index inside block
            new_j = ceil(value / y_step);
            new_i = value - (new_j-1)*y_step;

            actual_i = i+new_i-1;
            actual_j = j+new_j-1;           
            actual_pixel = image_luminance(actual_i, actual_j);
            ab = image_lab(actual_i, actual_j, 2:3);
            block = padded(actual_i:actual_i+amt_to_pad*2, actual_j:actual_j+amt_to_pad*2);
            fv.sds(idx) = std(block(:));
            fv.luminance(idx) = actual_pixel(:, :);
            fv.ab(idx, :) = ab;
            jittered.points(idx, :) = [actual_i, actual_j];
            idx = idx + 1;
            
        end
    end 
    jittered.fv = [fv.luminance(:), fv.sds(:)];
    jittered.ab = fv.ab;
end