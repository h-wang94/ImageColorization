function [src_image, tgt_image] = LoadImages(src_name, tgt_name, images_folder) 
%% load source and target images
    if nargin < 2
        images_folder = '../images/';
    end
    
    tgt_full = strcat(images_folder, tgt_name);
    src_full = strcat(images_folder, src_name);
    
    tgt_image = imread(tgt_full);
    if (size(tgt_image,3) > 1)
        tgt_image = rgb2gray(tgt_image);
    end
    tgt_image = im2double(tgt_image);
    src_image = im2double(imread(src_full));
end

