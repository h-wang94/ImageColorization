%%
clc; clear all;
%%
GRAPH = true;
SAVE = false;
%% load images
img_name = '1.jpg';
imgs = load_images(img_name);
target_image = imgs.target_image;
gsource_image = imgs.gsource_image;
csource_image = imgs.csource_image;

%% to LAB color space
lab_csource = rgb2lab(csource_image);
lab_gsource = rgb2lab(gsource_image);
lab_target = rgb2gray(target_image); % pixel values are luminance
%% map luminance to target luminance
lab_csource = luminance_remap(lab_csource, lab_target);
lab_gsource = luminance_remap(lab_gsource, lab_target);

%% sample pixels in colored images
%% save images
success = save_image(ctarget_image, img_name, SAVE);
