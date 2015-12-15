function transferred = image_colorization_brute_force(target, gsource, csource)
    %% get std in neighborhood
    csource.sds = sd_neighborhood(csource.luminance, 5);
    gsource.sds = sd_neighborhood(gsource.luminance, 5);
    target.sds = sd_neighborhood(target.luminance, 5);
    %%
    csource.fv = compute_fv(csource);
    gsource.fv = compute_fv(gsource);
    target.fv = compute_fv(target);

    transferred = transfer_fv(csource, gsource, target);
end