function sample_points = sampling_jittered(image, num_samples)
    % assumes num_samples is square number. num_samples = n^2 where n is an
    % integer
    
    n = sqrt(num_samples);
    [y, x, z] = size(image);
    y_step = floor(y / n);
    x_step = floor(x / n);
    sample_idx = randi([1, y_step*x_step], 1, num_samples);
end