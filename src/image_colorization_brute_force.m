function transferred = image_colorization_brute_force(target, gsource, csource)

%     %> For each pixel, compute std in neighboorhood
%     csource.sds = sd_neighborhood(csource.luminance, 5);
%     gsource.sds = sd_neighborhood(gsource.luminance, 5);
%     target.sds = sd_neighborhood(target.luminance, 5);
%     
%     %> For each pixel, generate pair of (luminance,std)
%     csource.fv = compute_fv(csource);
%     gsource.fv = compute_fv(gsource);
%     target.fv = compute_fv(target);

    %TEST:
    csource.fv = compute_fv_mod(csource.luminance);
    gsource.fv = compute_fv_mod(gsource.luminance);
    target.fv = compute_fv_mod(target.luminance);
    
    %> Transfer the chrominance from source to target.
    transferred = transfer_fv(csource, gsource, target);
end