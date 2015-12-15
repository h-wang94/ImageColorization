function pyramids = create_gauss_pyramid(image)
    if ndims(image) == 3
        [y, x, ~] = size(image);
    elseif ndims(image) == 2
        [y, x] = size(image);
    end
    pyramids = {};
    level = 1;
    while y > 100 & x > 100
        pyramids(level).pyramid = impyramid(image, 'reduce');
        image = pyramids(level).pyramid;
        level = level + 1;
        [y, x] = size(image);
    end
end