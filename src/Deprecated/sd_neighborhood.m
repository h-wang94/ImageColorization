function sds = sd_neighborhood(image, neighborhood_size)
    if nargin < 2
        neighborhood_size = 5;
    end
    amt_to_pad = (neighborhood_size - 1)/2;
    [y, x] = size(image);
    sds = zeros([y, x]);
    padded = padarray(image, [amt_to_pad, amt_to_pad]);
    
    for i = amt_to_pad+1:amt_to_pad+y
        for j = amt_to_pad+1:amt_to_pad+x
            region = padded(i-amt_to_pad:i+amt_to_pad, j-amt_to_pad:j+amt_to_pad);
            sd = std(region(:));
            sds(i-floor(neighborhood_size/2), ...
				j-floor(neighborhood_size/2)) = sd;
        end
    end
end