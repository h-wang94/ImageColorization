%%
clc; clear all;
%%
GRAPH = true;
SAVE = false;
%% load images
img_name = '1.jpg';
imgs = load_images(img_name);
target = {};
gsource = {};
csource = {};
target.image = imgs.target_image;
gsource.image = imgs.gsource_image;
csource.image = imgs.csource_image;


%% to LAB color space
csource.lab = rgb2lab(csource.image);
gsource.lab = rgb2lab(gsource.image);
target.lab = rgb2gray(target.image); % pixel values are luminance

%% map luminance to target luminance
csource.luminance = luminance_remap(csource.lab, target.lab);
gsource.luminance = luminance_remap(gsource.lab, target.lab);
target.luminance = target.lab;
%% get std in neighborhood
csource.sds = sd_neighborhood(csource.luminance, 5);
gsource.sds = sd_neighborhood(gsource.luminance, 5);
target.sds = sd_neighborhood(target.luminance, 5);
%%
csource.fv = compute_fv(csource);
gsource.fv = compute_fv(gsource);
target.fv = compute_fv(target);
%% 
transferred = transfer_fv(csource, gsource, target);
transferred(:,:,1) = transferred(:,:,1) * 100;
%%
new_image = lab2rgb(transferred);
%%
imshow(new_image);
%% save images
success = save_image(ctarget_image, img_name, SAVE);
