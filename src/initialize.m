%%
clc; clear all;
%%
tic;
GRAPH = false;
SAVE = false;
OPTION = 0; % 0 for brute-force, 1, for jitter-sampling

%% load images
img_name = '8.jpg';
imgs = load_images(img_name);
target = {};
gsource = {};
csource = {};
target.image = imgs.target_image;
gsource.image = imgs.gsource_image;
csource.image = imgs.csource_image;
tic;
%% to LAB color space
csource.lab = rgb2lab(csource.image);
if ndims(gsource.image) == 3
    gsource.lab = rgb2lab(gsource.image);
else
    gsource.lab = gsource.image;
end
if ndims(target.image) == 3
    target.lab = rgb2gray(target.image); 
else
    target.lab = target.image;
end

%% map luminance to target luminance
csource.luminance = luminance_remap(csource.lab, target.lab);
gsource.luminance = luminance_remap(gsource.lab, target.lab);
target.luminance = target.lab;
% pixel values are luminance

%% 
if OPTION == 0
    transferred = image_colorization_brute_force(target, gsource, csource);
    new_name = strcat('bf_', img_name);
elseif OPTION == 1
    transferred = image_colorization_jitter_sampling(target, csource, GRAPH);
    new_name = strcat('jitter_', img_name);
end
%%
new_image = lab2rgb(transferred);

toc;
%%
imshow(new_image);
%% save images
success = save_image(new_image, new_name, SAVE);
toc;