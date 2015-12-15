function transferred = image_colorization_jitter_sampling(target, csource, GRAPH)
    jittered = sampling_jittered(csource.luminance, csource.lab, 256);
    if GRAPH
        imshow(csource.image);
        hold on;
        scatter(jittered.points(:,2), jittered.points(:,1), 'r');
    end
    target.sds = sd_neighborhood(target.luminance, 5);
    target.fv = compute_fv(target);
    transferred = transfer_sample_fv(jittered, target);
end