function transferred = image_colorization_pyramid(target, gsource, csource)
    %%
    csource.gpyramid = create_gauss_pyramid(csource.image);
    gsource.gpyramid = create_gauss_pyramid(gsource.image);
    target.gpyramid = create_gauss_pyramid(target.image);


end