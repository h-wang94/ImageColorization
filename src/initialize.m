%%
% clc
clear all; close all;

%%
tic;
GRAPH = true;
SAVE = false;
SAMPLE_METHOD = 1; % 0 for brute-force, 1, for jitter-sampling

%TODO:INPUT FROM FILE

%% Input data (source and target images)
src_name = 'dog2.jpg'
tgt_name = 'dog1.jpg'

imgs = LoadImages(src_name, tgt_name, '../data/');

target = {};
source = {};
target.image = imgs.target_image;
source.image = imgs.source_image;

tic;
%% Color space conversion
source.lab = rgb2lab(source.image);
target.lab = target.image;

%% Map luminance to target luminance
%TODO: verificar necessidade e algoritmo.
source.luminance = luminance_remap(source.lab, target.lab);
target.luminance = target.lab;

%% Main colorization

%TODO: separar feature extraction da otimização.
if SAMPLE_METHOD == 0
    tgt_color = image_colorization_brute_force(target, gsource, source);
    new_name = strcat('bf_', img_name);
elseif SAMPLE_METHOD == 1
    tgt_color = image_colorization_jitter_sampling(target, source, GRAPH);
    new_name = strcat('jitter_', img_name);
end

%% Color space reconversion
tgt_color = lab2rgb(tgt_color);

toc;
%% Show results
figure;
imshow(source.image)
title('Source image');

figure;
imshow(tgt_color);
title('Colorized result');

%% save images
success = save_image(tgt_color, new_name, SAVE);
toc;