function [gauss_stack, lap_stack] = pyramids_rgb(img, level)
    N = ndims(img);
    gauss_stack = [];
    lap_stack = [];
    for i = 1:N
        [g_stack, l_stack] = pyramids(img(:,:,i), level);
        gauss_stack = cat(4, gauss_stack, g_stack);
        lap_stack = cat(4, lap_stack, l_stack);
    end
    gauss_stack = permute(gauss_stack, [1, 2, 4, 3]);
    lap_stack = permute(lap_stack, [1, 2, 4, 3]);
end